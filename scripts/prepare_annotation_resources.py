#!/usr/bin/env python3
"""Prepare only public/local annotation resources; no patient variants are queried."""
import argparse
import ctypes
import concurrent.futures
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import time
import threading
import urllib.parse
import urllib.request
import uuid


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def bsd_sum(path):
    value = 0
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            for byte in block:
                value = ((value >> 1) | ((value & 1) << 15))
                value = (value + byte) & 65535
    return value


def command(argv, timeout=None):
    return subprocess.run(argv, check=True, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, timeout=timeout).stdout


def ensure_vep_image(cfg, pull=False):
    image = cfg['vep']['image']
    need(re.fullmatch(r'ensemblorg/ensembl-vep@sha256:[0-9a-f]{64}', image), 'A pinned official VEP image is required.')
    docker = ['docker']
    try:
        command(docker + ['info'], timeout=30)
    except (subprocess.CalledProcessError, FileNotFoundError):
        docker = ['sudo', '-n', 'docker']
        command(docker + ['info'], timeout=30)
    if pull:
        print('START pinned VEP image pull/check before resource download', flush=True)
        command(docker + ['pull', image], timeout=900)
    version = command(docker + ['run', '--rm', '--network', 'none', image, 'vep', '--help'], timeout=120)
    need(re.search(r'\b115(?:\.2)?\b', version), 'Installed VEP does not match the fixed cache release.')
    print('PASS pinned VEP CLI version before large download', flush=True)
    return version


def publisher_check(path, item):
    checks = {'actual_bytes': path.stat().st_size}
    if item.get('bytes') is not None:
        need(path.stat().st_size == item['bytes'], 'Resource byte count mismatch: ' + item['name'])
    if item.get('sha256'):
        need(sha256(path) == item['sha256'], 'Resource SHA256 mismatch: ' + item['name'])
    if item.get('md5'):
        h = hashlib.md5()
        with path.open('rb') as stream:
            for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
                h.update(chunk)
        need(h.hexdigest() == item['md5'], 'Publisher MD5 mismatch: ' + item['name'])
    if item.get('bsd_sum') is not None:
        # GNU sum is the publisher's BSD rotating checksum, not a cryptographic hash.
        got = command(['sum', str(path)], timeout=600).split()
        need(int(got[0]) == item['bsd_sum'] and int(got[1]) == item['bsd_blocks'],
             'Publisher BSD sum/block count mismatch: ' + item['name'])
        checks['publisher_checksum'] = dict(algorithm='BSD sum', expected_sum=item['bsd_sum'],
            actual_sum=int(got[0]), expected_blocks=item['bsd_blocks'], actual_blocks=int(got[1]))
    if item.get('required_text'):
        prefix = path.read_text(encoding='utf-8-sig')[:10000]
        need(all(text in prefix for text in item['required_text']), 'Snapshot header/date changed: ' + item['name'])
    return checks


def download(item, target):
    url = item['url']
    host = urllib.parse.urlsplit(url).hostname
    need(url.startswith('https://') and host in {
        'ftp.ensembl.org', 'ftp.ncbi.nlm.nih.gov', 'search.clinicalgenome.org',
        'ftp.clinicalgenome.org', 'storage.googleapis.com', 'github.com'}, 'Unapproved public resource host.')
    cap = item.get('bytes') or item.get('max_bytes', 5 * 1024 * 1024)
    print('START ' + item['name'] + ' version=' + item['version'], flush=True)
    h = hashlib.sha256()
    count = 0
    previous = 0
    started = time.monotonic()
    with urllib.request.urlopen(url, timeout=120) as response, target.open('xb') as stream:
        need(urllib.parse.urlsplit(response.url).scheme == 'https', 'Refusing non-TLS redirect.')
        for chunk in iter(lambda: response.read(256 * 1024), b''):
            count += len(chunk)
            need(count <= cap, 'Resource exceeds explicit download limit.')
            stream.write(chunk)
            h.update(chunk)
            need(time.monotonic() - started < item.get('max_download_seconds', 1800), 'Resource download time limit reached; partial preserved.')
            if count - previous >= 64 * 1024 ** 2:
                print(item['name'] + ': ' + str(count) + ' bytes', flush=True)
                previous = count
    publisher_check(target, item)
    print('PASS ' + item['name'] + ' bytes=' + str(count), flush=True)
    return dict(name=item['name'], version=item['version'], source=url,
                path=target.name, bytes=count, sha256=h.hexdigest(),
                publisher_validation=('MD5' if item.get('md5') else
                    'BSDsum+bytes; noncryptographic' if item.get('bsd_sum') is not None else
                    'SHA256' if item.get('sha256') else 'explicit dated snapshot; local SHA256'),
                retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(), status='available')


def prefetch_cache(item, destination, workers=16, reserve_gib=150):
    import fcntl
    destination.parent.mkdir(parents=True, exist_ok=True)
    lockpath = destination.with_name(destination.name + '.agent-download.lock')
    need(not lockpath.is_symlink(), 'Refusing a download-lock symlink.')
    with lockpath.open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('WAITING for the existing direct-server download process; no second transfer started.', flush=True)
            fcntl.flock(lock, fcntl.LOCK_EX)
        return _prefetch_cache_locked(item, destination, workers, reserve_gib)


def _prefetch_cache_locked(item, destination, workers=16, reserve_gib=150):
    """Bounded HTTPS ranges into ONE owned partial, with resumable per-range hashes.

    This optional route is useful when a workstation's path to Ensembl is faster
    than the VM's. It never contacts a patient-specific endpoint.
    """
    need(1 <= workers <= 32 and destination.name == item['destination'], 'Invalid prefetch destination/workers.')
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        checks = publisher_check(destination, item)
        final_sha = sha256(destination)
        validation = dict(status='CACHE_PREFETCH_PASS', bytes=item['bytes'], sha256=final_sha, **checks)
        destination.with_name(destination.name + '.agent-validation.json').write_text(json.dumps(validation, indent=2) + '\n')
        print(json.dumps(dict(status='CACHE_ALREADY_PRESENT', bytes=destination.stat().st_size, sha256=final_sha)))
        return final_sha
    partial = destination.with_name(destination.name + '.agent-partial')
    state_path = destination.with_name(destination.name + '.agent-parts.json')
    request = urllib.request.Request(item['url'], method='HEAD')
    with urllib.request.urlopen(request, timeout=30) as response:
        identity = dict(url=item['url'], bytes=item['bytes'], etag=response.headers.get('ETag'),
                        last_modified=response.headers.get('Last-Modified'), chunk_bytes=64 * 1024 ** 2)
        need(int(response.headers['Content-Length']) == item['bytes'], 'Prefetch source size changed.')
    need(identity['etag'] or identity['last_modified'], 'Range transfers require a stable ETag or Last-Modified identity.')
    if partial.exists() or state_path.exists():
        need(partial.is_file() and state_path.is_file() and not partial.is_symlink() and not state_path.is_symlink(), 'Unrecognized partial; refusing overwrite.')
        state = json.loads(state_path.read_text())
        need(state['identity'] == identity and partial.stat().st_size == item['bytes'], 'Partial belongs to a different public object.')
        fd = os.open(partial, os.O_RDWR)
    else:
        need(shutil.disk_usage(destination.parent).free >= item['bytes'] + reserve_gib * 1024 ** 3,
             'Insufficient workstation space for cache and explicit reserve.')
        fd = os.open(partial, os.O_RDWR | os.O_CREAT | os.O_EXCL, 0o600)
        os.ftruncate(fd, item['bytes'])
        state = dict(identity=identity, completed={})
        state_path.write_text(json.dumps(state))
    lock = threading.Lock()
    deadline = time.monotonic() + item.get('max_download_seconds', 21600)
    chunk_size = identity['chunk_bytes']
    chunks = (item['bytes'] + chunk_size - 1) // chunk_size
    def save_state():
        temporary = state_path.with_name(state_path.name + '.tmp')
        temporary.write_text(json.dumps(state))
        os.replace(temporary, state_path)
    def fetch(index):
        start, end = index * chunk_size, min((index + 1) * chunk_size, item['bytes']) - 1
        previous = state['completed'].get(str(index))
        if previous:
            h = hashlib.sha256()
            for position in range(start, end + 1, 1024 ** 2):
                h.update(os.pread(fd, min(1024 ** 2, end + 1 - position), position))
            need(h.hexdigest() == previous, 'A completed owned prefetch range changed.')
            return
        headers = {'Range': f'bytes={start}-{end}'}
        if identity['etag']:
            headers['If-Match'] = identity['etag']
        for attempt in range(3):
            try:
                with urllib.request.urlopen(urllib.request.Request(item['url'], headers=headers), timeout=60) as response:
                    need(response.status == 206 and response.headers.get('Content-Range') == f'bytes {start}-{end}/{item["bytes"]}',
                         'Server did not honor the exact requested public range.')
                    need(response.headers.get('ETag') == identity['etag'] and response.headers.get('Last-Modified') == identity['last_modified'],
                         'Object identity changed during range transfer.')
                    if response.headers.get('Content-Length'):
                        need(int(response.headers['Content-Length']) == end - start + 1, 'Range Content-Length mismatch.')
                    position, h = start, hashlib.sha256()
                    for block in iter(lambda: response.read(1024 ** 2), b''):
                        need(time.monotonic() < deadline, 'Direct download time limit reached; owned partial is resumable.')
                        need(position + len(block) <= end + 1, 'Range response exceeded its bounds.')
                        need(os.pwrite(fd, block, position) == len(block), 'Short local range write.')
                        h.update(block)
                        position += len(block)
                    need(position == end + 1, 'Short public range response.')
                with lock:
                    state['completed'][str(index)] = h.hexdigest()
                    save_state()
                    complete_bytes = sum(min(chunk_size, item['bytes'] - int(i) * chunk_size) for i in state['completed'])
                    print('CACHE_PROGRESS ' + str(complete_bytes) + '/' + str(item['bytes']), flush=True)
                    need(shutil.disk_usage(destination.parent).free >= reserve_gib * 1024 ** 3,
                         'Workstation free space fell below the configured reserve.')
                return
            except (OSError, ValueError):
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(fetch, range(chunks)))
        os.fsync(fd)
    finally:
        os.close(fd)
    checks = publisher_check(partial, item)
    final_sha = sha256(partial)
    os.link(partial, destination)  # atomic exclusive publication, same inode/no second 25GB copy
    partial.unlink()
    destination.with_name(destination.name + '.sha256').write_text(final_sha + '  ' + destination.name + '\n')
    validation = dict(status='CACHE_PREFETCH_PASS', bytes=item['bytes'], sha256=final_sha, **checks)
    destination.with_name(destination.name + '.agent-validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(dict(status='CACHE_PREFETCH_PASS', bytes=item['bytes'], sha256=final_sha)), flush=True)
    return final_sha


def extract_cache(archive, directory, *, limit, reserve):
    print('START cache tar header validation', flush=True)
    total, names = 0, set()
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            path = Path(member.name)
            need(not path.is_absolute() and '..' not in path.parts and path.parts
                 and path.parts[0] == 'homo_sapiens', 'Unsafe cache archive path.')
            need(member.name not in names and (member.isfile() or member.isdir()),
                 'Duplicate path or link/special file in cache archive.')
            names.add(member.name)
            total += member.size
            need(total <= limit, 'Expanded cache exceeds configured size limit.')
    need(total + reserve <= shutil.disk_usage(directory.parent).free,
         'Insufficient space for expanded cache plus analysis reserve.')
    print('PASS tar headers; expanded bytes=' + str(total), flush=True)
    directory.mkdir()
    count = 0
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            destination = directory / member.name
            if member.isdir():
                destination.mkdir(parents=True, exist_ok=True)
            else:
                destination.parent.mkdir(parents=True, exist_ok=True)
                with tar.extractfile(member) as source, destination.open('xb') as out:
                    shutil.copyfileobj(source, out, 8 * 1024 * 1024)
                os.chmod(destination, 0o644)
                count += 1
    print('PASS extracted cache files=' + str(count), flush=True)
    return total


def recipe_sources(config_path):
    cfg = json.loads(config_path.read_text())
    result = []
    if cfg.get('base_config'):
        base = (config_path.parent / cfg['base_config']).resolve(strict=True)
        need(base.parent == config_path.parent.resolve(), 'Base recipe must be beside final recipe.')
        result += recipe_sources(base)
    return result + [item for item in cfg['resources'] if item.get('enabled', False)]


def source_cache_directory(root, item):
    identity = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()
    return root / '.public-source-cache' / identity


def prepare_small_sources(config_path, root):
    """Stage independently useful public databases before waiting for the large cache."""
    root.mkdir(parents=True, exist_ok=True)
    need(shutil.disk_usage(root).free >= 150 * 1024 ** 3, 'Small-resource staging requires 150 GiB free reserve.')
    for item in recipe_sources(config_path):
        if item['name'] == 'vep_cache':
            continue
        need(Path(item['destination']).name == item['destination'], 'Unsafe public source basename.')
        target = source_cache_directory(root, item)
        target.parent.mkdir(exist_ok=True)
        if target.exists():
            receipt = json.loads((target / 'resource.json').read_text())
            file = target / item['destination']
            need(file.is_file() and sha256(file) == receipt['sha256'] and file.stat().st_size == receipt['bytes'], 'An existing public source changed.')
            continue
        stage = target.parent / ('.partial-source-' + uuid.uuid4().hex)
        stage.mkdir(mode=0o700)
        receipt = download(item, stage / item['destination'])
        (stage / 'resource.json').write_text(json.dumps(receipt, indent=2) + '\n')
        promote(stage, target)
        print('SMALL_RESOURCE_READY ' + item['name'] + ' ' + str(target / item['destination']), flush=True)


def tree_manifest(directory):
    entries = []
    for path in sorted(directory.rglob('*')):
        need(not path.is_symlink(), 'Symlink in prepared resource tree.')
        if path.is_file():
            entries.append(dict(path=path.relative_to(directory).as_posix(),
                                bytes=path.stat().st_size, sha256=sha256(path)))
    return entries


def promote(source, target):
    libc = ctypes.CDLL(None, use_errno=True)
    need(hasattr(libc, 'renameat2'), 'Atomic no-replace promotion requires Linux renameat2.')
    result = libc.renameat2(-100, os.fsencode(source), -100, os.fsencode(target), 1)
    if result != 0:
        raise OSError(ctypes.get_errno(), 'No-replace resource promotion failed.')


def cache_metadata_records(info, stage):
    """Require actual official cache metadata; MANE release is not always reported."""
    text = info.read_text(errors='replace')
    values = {}
    for line in text.splitlines():
        if '\t' not in line:
            continue
        key, value = line.split('\t', 1)
        key = key.lower()
        need(key not in values, 'Duplicate cache metadata key: ' + key)
        values[key] = value.strip()
    required = {'species': 'homo_sapiens', 'assembly': 'GRCh38', 'source_assembly': 'GRCh38.p14',
        'source_gencode': 'GENCODE 49', 'source_sift': '6.2.1', 'source_polyphen': '2.2.3',
        'source_gnomade': 'v4.1', 'source_gnomadg': 'v4.1', 'source_dbsnp': '156'}
    for key, expected in required.items():
        need(values.get(key) == expected, 'Cache metadata version unexpected for ' + key + '; inspect preserved info.txt.')
    identity = dict(path=info.relative_to(stage).as_posix(), bytes=info.stat().st_size, sha256=sha256(info),
        documented_source='https://sep2025.archive.ensembl.org/info/docs/tools/vep/script/vep_cache.html')
    result = {}
    for name, keys in {'gnomad': ['source_gnomade', 'source_gnomadg'], 'sift': ['source_sift'], 'polyphen2': ['source_polyphen']}.items():
        result[name] = dict(identity, status='bundled_in_verified_cache',
            version='; '.join(key + '=' + values[key] for key in keys),
            scope='versioned_cache_allele_AF_only_no_AN_coverage_constraint_claim' if name == 'gnomad' else 'cache_component')
    mane = values.get('source_mane')
    result['mane_cache'] = dict(identity,
        status='bundled_in_verified_cache' if mane else 'source_version_unreported_in_cache_metadata',
        version=mane or None, scope='cache_MANE_source_metadata_only; runtime transcript labels require public VEP controls',
        external_mane_summary_is_distinct_resource=True,
        limitation='Missing source_MANE does not establish MANE1.4, absence of all MANE labels, or runtime mapping readiness.')
    return result


def prepare(config_path, root, execute=False, pull_image=False, prefetched_cache=None, prefetched_cache_sha256=None):
    config_bytes = config_path.read_bytes()
    cfg = json.loads(config_bytes)
    need(cfg.get('schema_version') == 1 and cfg['assembly'] == 'GRCh38', 'Unsupported resource configuration.')
    image = cfg['vep']['image']
    need(re.fullmatch(r'ensemblorg/ensembl-vep@sha256:[0-9a-f]{64}', image), 'An immutable official VEP image is required.')
    need(str(cfg['vep']['cache_version']) == '115' and cfg['vep']['version'] == '115.2', 'Expected matched VEP/cache 115.')
    label = cfg['bundle_id']
    need(re.fullmatch(r'[A-Za-z0-9_.-]+', label), 'Unsafe resource bundle ID.')
    sources = [item for item in cfg['resources'] if item.get('enabled', False)]
    need(len({item['name'] for item in sources}) == len(sources), 'Duplicate resource names.')
    for item in sources:
        need(Path(item['destination']).name == item['destination'], 'Source destinations must be basenames.')
        need(item.get('version'), 'Every resource requires an explicit version/snapshot date.')
    lock_sha = hashlib.sha256(config_bytes).hexdigest()
    target = root / label
    base_target, base_manifest = None, None
    if cfg.get('base_config'):
        base_config = (config_path.parent / cfg['base_config']).resolve(strict=True)
        need(base_config.parent == config_path.parent.resolve(), 'Base recipe must be beside the final recipe.')
        # Recursion occurs BEFORE acquiring the final preparation lock.
        base_target = prepare(base_config, root, execute, pull_image, prefetched_cache, prefetched_cache_sha256)
        if execute:
            base_manifest = json.loads((base_target / 'resource-manifest.json').read_text())
    if not execute:
        print(json.dumps(dict(status='DRY_RUN', target=str(target), image=image,
            resources=[{k: i[k] for k in ('name', 'version', 'destination')} for i in sources],
            known_download_bytes=sum(i.get('bytes', 0) for i in sources), analysis_started=False)))
        return target
    root.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (root / '.annotation-prep.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if target.exists():
            receipt = json.loads((target / 'resource-manifest.json').read_text())
            need(receipt['config_sha256'] == lock_sha, 'Existing resource bundle belongs to a different lock.')
            for entry in receipt['artifacts']:
                path = target / entry['path']
                need(path.is_file() and path.stat().st_size == entry['bytes'] and sha256(path) == entry['sha256'],
                     'Existing resource artifact changed; refusing overwrite.')
            print('PASS existing resource bundle verified: ' + str(target), flush=True)
            return target
        need(shutil.disk_usage(root).free >= max(150, cfg.get('min_free_gib', 150)) * 1024 ** 3,
             'At least 150 GiB free is required before resource preparation.')
        docker = ['docker']
        try:
            command(docker + ['info'], timeout=30)
        except (subprocess.CalledProcessError, FileNotFoundError):
            docker = ['sudo', '-n', 'docker']
            command(docker + ['info'], timeout=30)
        if pull_image:
            print('START pinned VEP image pull', flush=True)
            command(docker + ['pull', image], timeout=900)
        version = command(docker + ['run', '--rm', '--network', 'none', image, 'vep', '--help'], timeout=120)
        need(re.search(r'\b115(?:\.2)?\b', version), 'VEP tool version does not match cache 115.')
        stage = root / ('.partial-annotation-' + uuid.uuid4().hex)
        stage.mkdir(mode=0o700)
        (stage / 'resources.lock.json').write_bytes(config_bytes)
        (stage / 'vep-version.txt').write_text(version)
        downloaded = []
        if base_manifest:
            need(base_manifest['assembly'] == cfg['assembly'] and base_manifest['vep'] == cfg['vep'],
                 'Base recipe image/assembly/cache differs from final recipe.')
            print('START hardlink reuse of verified immutable base artifacts', flush=True)
            for name, item in base_manifest['resources'].items():
                if item.get('status') != 'available':
                    continue
                original = base_target / item['path']
                need(original.parent == base_target and not original.is_symlink(), 'Unexpected base resource layout.')
                os.link(original, stage / original.name)
                downloaded.append(dict(item, reused_from_base_manifest_sha256=sha256(base_target / 'resource-manifest.json')))
            shutil.copytree(base_target / 'cache', stage / 'cache', copy_function=os.link, symlinks=False)
            (stage / 'base-resource-manifest.json').write_bytes((base_target / 'resource-manifest.json').read_bytes())
            print('PASS validated immutable base hardlinks reused', flush=True)
        for item in sources:
            if item['name'] == 'vep_cache' and prefetched_cache is not None:
                need(re.fullmatch('[0-9a-f]{64}', prefetched_cache_sha256 or ''), 'Prefetched cache needs its expected transfer SHA256.')
                need(sha256(prefetched_cache) == prefetched_cache_sha256, 'Prefetched cache transfer SHA256 mismatch.')
                publisher_check(prefetched_cache, item)
                os.link(prefetched_cache, stage / item['destination'])
                downloaded.append(dict(name=item['name'], version=item['version'], source=item['url'],
                    path=item['destination'], bytes=prefetched_cache.stat().st_size, sha256=prefetched_cache_sha256,
                    publisher_validation='BSDsum+bytes; noncryptographic; transfer SHA256 verified',
                    status='available', retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
            elif (source_cache_directory(root, item) / 'resource.json').is_file():
                cached = source_cache_directory(root, item)
                record = json.loads((cached / 'resource.json').read_text())
                file = cached / item['destination']
                need(file.stat().st_size == record['bytes'] and sha256(file) == record['sha256'], 'A precached public database changed.')
                os.link(file, stage / item['destination'])
                downloaded.append(record)
            else:
                downloaded.append(download(item, stage / item['destination']))
        if base_manifest:
            expanded = base_manifest['cache_expanded_bytes']
        else:
            cache_item = next(i for i in sources if i['name'] == 'vep_cache')
            expanded = extract_cache(stage / cache_item['destination'], stage / 'cache',
                limit=cfg.get('max_cache_expanded_gib', 128) * 1024 ** 3,
                reserve=max(150, cfg.get('analysis_reserve_gib', 150)) * 1024 ** 3)
        info = stage / 'cache/homo_sapiens/115_GRCh38/info.txt'
        need(info.is_file(), 'The expected matching cache info.txt is absent.')
        # Preserve exact source annotation versions; do not relabel bundled AF as a newer gnomAD.
        cache_info = info.read_text(errors='replace')
        available = {i['name']: i for i in downloaded}
        for item in cfg.get('optional_adapters', []):
            record = dict(item)
            record.setdefault('status', 'unavailable_manual_or_licensed')
            available.setdefault(item['name'], record)
        available.update(cache_metadata_records(info, stage))
        manifest = dict(schema_version=1, status='PASS', assembly='GRCh38',
            config_sha256=lock_sha, vep=cfg['vep'], cache_info=cache_info,
            cache_expanded_bytes=expanded, resources=available,
            artifacts=tree_manifest(stage),
            completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
            patient_variants_sent_to_remote_services=False)
        (stage / 'resource-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        promote(stage, target)
        print('PASS resource manifest: ' + str(target / 'resource-manifest.json'), flush=True)
        return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--root', type=Path, default=Path('/srv/marina/resources/annotation'))
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--pull-image', action='store_true', help='Pull only the configured immutable official VEP image.')
    parser.add_argument('--prefetched-cache', type=Path, help='Optional already transferred public cache; publisher bytes/checksum are still verified.')
    parser.add_argument('--prefetched-cache-sha256', help='Required expected end-to-end transfer SHA256 for a prefetched cache.')
    parser.add_argument('--cache-only', type=Path, help='Download the locked public cache only to this basename, without Docker or preparation.')
    parser.add_argument('--download-workers', type=int, default=16)
    parser.add_argument('--reserve-gib', '--local-reserve-gib', dest='reserve_gib', type=int, default=150)
    parser.add_argument('--small-only', action='store_true', help='Stage verified public ClinVar/ClinGen/HPO/HGNC/MANE without waiting for the VEP cache.')
    parser.add_argument('--direct-cache-dir', type=Path, help='Direct SERVER download/resume directory; default sibling resources/.downloads.')
    args = parser.parse_args()
    os.umask(0o077)
    if args.cache_only:
        cfg = json.loads(args.config.read_text())
        if cfg.get('base_config'):
            cfg = json.loads((args.config.parent / cfg['base_config']).read_text())
        item = next(i for i in cfg['resources'] if i['name'] == 'vep_cache')
        if args.execute:
            prefetch_cache(item, args.cache_only, args.download_workers, args.reserve_gib)
        else:
            print(json.dumps(dict(status='DRY_RUN', public_cache=item['url'], destination=str(args.cache_only), bytes=item['bytes'])))
    elif args.small_only:
        if args.execute:
            prepare_small_sources(args.config, args.root)
        else:
            print(json.dumps(dict(status='DRY_RUN', small_resources=[item['name'] for item in recipe_sources(args.config) if item['name'] != 'vep_cache'])))
    else:
        if args.execute:
            ensure_vep_image(json.loads(args.config.read_text()), args.pull_image)
            prepare_small_sources(args.config, args.root)
            if args.prefetched_cache is None:
                item = next(i for i in recipe_sources(args.config) if i['name'] == 'vep_cache')
                directory = args.direct_cache_dir or args.root.parent / '.downloads'
                args.prefetched_cache = directory / item['destination']
                args.prefetched_cache_sha256 = prefetch_cache(item, args.prefetched_cache,
                    args.download_workers, max(150, args.reserve_gib))
        prepare(args.config, args.root, args.execute, args.pull_image, args.prefetched_cache, args.prefetched_cache_sha256)


if __name__ == '__main__':
    main()
