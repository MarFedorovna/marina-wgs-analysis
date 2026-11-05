#!/usr/bin/env python3
"""Install canonical/isoform AlphaMissense public scores directly on the server."""
import argparse
import base64
import csv
import ctypes
import datetime as dt
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
import uuid


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 ** 2), b''):
            h.update(chunk)
    return h.hexdigest()


def command(argv, timeout=180):
    return subprocess.run(argv, check=True, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, timeout=timeout).stdout


def promote(source, target):
    libc = ctypes.CDLL(None, use_errno=True)
    need(hasattr(libc, 'renameat2'), 'Linux atomic no-replace rename is required.')
    if libc.renameat2(-100, os.fsencode(source), -100, os.fsencode(target), 1):
        raise OSError(ctypes.get_errno(), 'Refusing resource promotion over an existing target.')


def download(item, path):
    url = ('https://storage.googleapis.com/storage/v1/b/dm_alphamissense/o/'
           + urllib.parse.quote(item['object'], safe='') + '?alt=media&generation=' + item['generation'])
    h, md5, count, previous = hashlib.sha256(), hashlib.md5(), 0, 0
    started = time.monotonic()
    print('START AlphaMissense ' + item['name'] + ' bytes=' + str(item['bytes']), flush=True)
    with urllib.request.urlopen(url, timeout=60) as response, path.open('xb') as stream:
        for chunk in iter(lambda: response.read(1024 ** 2), b''):
            count += len(chunk)
            need(count <= item['bytes'] and time.monotonic() - started < 1800, 'AlphaMissense size/time guard reached; partial preserved.')
            h.update(chunk)
            md5.update(chunk)
            stream.write(chunk)
            if count - previous >= 128 * 1024 ** 2:
                print('ALPHAMISSENSE_PROGRESS ' + item['name'] + ' ' + str(count), flush=True)
                previous = count
    expected_md5 = item.get('md5') or base64.b64decode(item['md5_base64']).hex()
    need(count == item['bytes'] and md5.hexdigest() == expected_md5, 'AlphaMissense publisher size/MD5 mismatch.')
    print('PASS AlphaMissense ' + item['name'], flush=True)
    return dict(name=item['name'], path=path.name, bytes=count, sha256=h.hexdigest(), md5=md5.hexdigest(),
                source=url, generation=item['generation'], status='available')


def schema(path, canonical):
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.startswith('#CHROM'):
                header = line.lstrip('#').rstrip('\r\n').split('\t')
                first = next(csv.DictReader(stream, fieldnames=header, delimiter='\t'))
                required = {'CHROM', 'POS', 'REF', 'ALT', 'genome', 'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class'}
                need(required.issubset(header) and (not canonical or 'uniprot_id' in header), 'AlphaMissense schema differs from the lock.')
                need(first['genome'] == 'hg38' and 0 <= float(first['am_pathogenicity']) <= 1
                     and re.fullmatch(r'ENST\d+\.\d+', first['transcript_id']), 'AlphaMissense build/transcript/score validation failed.')
                return dict(columns=header, public_first_record=first)
    raise ValueError('AlphaMissense #CHROM schema header not found.')


def prepare(lock_path, root, execute=False):
    lock_bytes = lock_path.read_bytes()
    cfg = json.loads(lock_bytes)
    need(cfg.get('schema_version') == 1 and cfg['assembly'] == 'GRCh38'
         and re.fullmatch(r'[A-Za-z0-9_.-]+', cfg['bundle_id']), 'Unsupported AlphaMissense source recipe.')
    target = root / 'resources/annotation-supplements' / cfg['bundle_id']
    if not execute:
        print(json.dumps(dict(status='DRY_RUN', target=str(target), bytes=sum(x['bytes'] for x in cfg['resources']), patient_queries=False)))
        return
    os.umask(0o077)
    target.parent.mkdir(parents=True, exist_ok=True)
    with (target.parent / '.alphamissense.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if target.exists():
            manifest = json.loads((target / 'resource-manifest.json').read_text())
            need(manifest['lock_sha256'] == hashlib.sha256(lock_bytes).hexdigest(), 'Existing AlphaMissense recipe differs; refusing overwrite.')
            for item in manifest['artifacts']:
                file = target / item['path']
                need(file.is_file() and file.stat().st_size == item['bytes'] and sha256(file) == item['sha256'], 'Existing AlphaMissense artifact changed.')
            print('PASS existing AlphaMissense resources verified', flush=True)
            return
        need(shutil.disk_usage(target.parent).free >= sum(x['bytes'] for x in cfg['resources']) + cfg['min_free_gib'] * 1024 ** 3,
             'AlphaMissense preparation would violate the 150 GiB reserve.')
        tabix = shutil.which('tabix')
        if not tabix:
            command(['sudo', '-n', 'apt-get', 'update'], 300)
            command(['sudo', '-n', 'apt-get', 'install', '-y', '--no-install-recommends', 'tabix'], 300)
            tabix = shutil.which('tabix')
        need(tabix, 'Official Ubuntu tabix is required.')
        version = command([tabix, '--version'])
        need(re.search(r'\b' + re.escape(cfg['tabix_version']) + r'\b', version), 'Review a new recipe for a different tabix version.')
        stage = target.parent / ('.partial-alphamissense-' + uuid.uuid4().hex)
        stage.mkdir(mode=0o700)
        (stage / 'sources.lock.json').write_bytes(lock_bytes)
        resources, schemas = {}, {}
        for item in cfg['resources']:
            need(Path(item['object']).name == item['object'], 'Unsafe AlphaMissense source object.')
            path = stage / item['object']
            record = download(item, path)
            record.update(version=cfg['dataset_version'], assembly='GRCh38', transcript_annotation=cfg['transcript_annotation'])
            if item['name'] in ('canonical', 'isoforms'):
                schemas[item['name']] = schema(path, item['name'] == 'canonical')
                print('START tabix AlphaMissense ' + item['name'], flush=True)
                command([tabix, '-s', '1', '-b', '2', '-e', '2', str(path)], 1800)
                record['index_path'] = path.name + '.tbi'
                record['contigs'] = command([tabix, '-l', str(path)]).splitlines()
                need('chr1' in record['contigs'], 'AlphaMissense index lacks primary chromosomes.')
                probe = schemas[item['name']]['public_first_record']
                result = command([tabix, str(path), probe['CHROM'] + ':' + probe['POS'] + '-' + probe['POS']])
                need(probe['transcript_id'] in result and probe['protein_variant'] in result, 'Public AlphaMissense local index probe failed.')
                print('PASS tabix/public lookup ' + item['name'], flush=True)
            resources[item['name']] = record
        artifacts = [dict(path=p.relative_to(stage).as_posix(), bytes=p.stat().st_size, sha256=sha256(p))
                     for p in sorted(stage.iterdir()) if p.is_file()]
        manifest = dict(schema_version=1, status='PASS', assembly='GRCh38', dataset_version=cfg['dataset_version'],
            source_reference=cfg['source_reference'], transcript_annotation=cfg['transcript_annotation'],
            license=cfg['license'], license_source=cfg['license_source'], citation=cfg['citation'],
            lock_sha256=hashlib.sha256(lock_bytes).hexdigest(), resources=resources, artifacts=artifacts, schemas=schemas,
            tabix=dict(binary=tabix, version=version, sha256=sha256(Path(tabix))),
            installer_sha256=sha256(Path(__file__)), completed_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
            prediction_model_labels_are_not_ACMG_classes=True, patient_queries=False)
        (stage / 'resource-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        promote(stage, target)
        print('ALPHAMISSENSE_READY ' + str(target / 'resource-manifest.json'), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lock', type=Path, default=Path(__file__).resolve().parents[1] / 'config/alphamissense-resources.lock.json')
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    prepare(args.lock, args.root, args.execute)
