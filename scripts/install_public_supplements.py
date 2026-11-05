#!/usr/bin/env python3
"""Bounded public-only supplements; default status. Does not invoke annotation."""
import argparse
import concurrent.futures
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET

CODE = Path(__file__).resolve().parents[1]


def utc():
    return dt.datetime.now(dt.timezone.utc)


def load(path):
    return json.loads(path.read_text()) if path.exists() else {}


def save(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    if path.is_symlink() or tmp.is_symlink():
        raise ValueError('Refusing symlink')
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def download(item, folder, deadline):
    receipt_path = folder / (item['id'] + '-receipt.json')
    target = folder / item['file']; partial = target.with_suffix(target.suffix + '.partial')
    identity_path = target.with_suffix(target.suffix + '.identity.json')
    status_path = folder / (item['id'] + '-status.json')
    receipt = load(receipt_path)
    if receipt.get('status') == 'PASS':
        assert target.is_file() and target.stat().st_size == receipt['bytes'] and digest(target) == receipt['sha256']
        assert receipt['source']['url'] == item['url']
        return receipt
    assert not target.exists(), 'Unreceipted target will not be overwritten'
    for attempt in range(1, 5):
        try:
            assert utc() < deadline, 'Original budget/runtime deadline reached'
            request = urllib.request.Request(item['url'], method='HEAD', headers={'Accept-Encoding': 'identity'})
            try:
                response = urllib.request.urlopen(request, timeout=30)
            except urllib.error.HTTPError as error:
                # Documented CSpec JSON GET works while its HEAD currently returns 400.
                if error.code not in (400, 405) or not item['id'].startswith('clingen_cspec_'):
                    raise
                response = urllib.request.urlopen(urllib.request.Request(item['url'], headers={'Accept-Encoding':'identity'}), timeout=30)
            with response:
                identity = dict(url=response.url, etag=response.headers.get('ETag'),
                    last_modified=response.headers.get('Last-Modified'), bytes=response.headers.get('Content-Length'))
            assert identity['url'] == item['url'], 'Unreviewed redirect'
            old = load(identity_path)
            assert not old or old == identity, 'Public source changed; review a new version'
            save(identity_path, identity)
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {'Accept-Encoding': 'identity'}
            if offset:
                headers['Range'] = 'bytes=%d-' % offset
                if identity['last_modified']:
                    headers['If-Range'] = identity['last_modified']
            save(status_path, dict(status='DOWNLOADING', attempt=attempt, received_bytes=offset,
                source_url=item['url'], updated_utc=utc().isoformat(), pid=os.getpid()))
            with urllib.request.urlopen(urllib.request.Request(item['url'], headers=headers), timeout=30) as response:
                assert response.url == item['url'], 'Unreviewed redirect'
                if offset:
                    assert response.status == 206 and response.headers.get('Content-Range', '').startswith('bytes %d-' % offset), 'Server did not honor safe resume; partial preserved'
                with partial.open('ab' if offset else 'xb') as out:
                    last = time.monotonic()
                    while True:
                        assert utc() < deadline, 'Original budget/runtime deadline reached'
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        offset += len(block)
                        assert offset <= item['max_bytes'], 'Resource exceeds reviewed byte ceiling'
                        out.write(block)
                        if time.monotonic() - last >= 5:
                            save(status_path, dict(status='DOWNLOADING', attempt=attempt, received_bytes=offset,
                                source_url=item['url'], updated_utc=utc().isoformat(), pid=os.getpid()))
                            last = time.monotonic()
            expected = item.get('expected_bytes') or (int(identity['bytes']) if identity['bytes'] else None)
            assert expected is None or offset == expected, 'Publisher length differs'
            if item['format'] == 'xml':
                doc = ET.parse(partial)
                assert doc.getroot().tag == 'JDBOR', 'Unexpected Orphadata schema'
                validation = dict(root=doc.getroot().tag, root_attributes=doc.getroot().attrib)
            else:
                doc = load(partial)
                if item['format'] == 'json_object':
                    assert isinstance(doc.get('data'), dict) and doc['data'].get('entCount'), 'Unexpected CSpec metadata'
                    validation = dict(registry_entities=doc['data']['entCount'])
                else:
                    assert isinstance(doc.get('data'), list) and doc['data'], 'Unexpected CSpec schema'
                    validation = dict(entities=len(doc['data']), pagination={k:v for k,v in doc.items() if k != 'data'})
            receipt = dict(status='PASS', scope='PUBLIC_SOURCE_DATA_ONLY_NOT_PATIENT_ANNOTATION', source=item,
                transport_identity=identity, bytes=offset, sha256=digest(partial), retrieved_utc=utc().isoformat(),
                validation=validation, publisher_sha256_provided=False)
            partial.replace(target); save(receipt_path, receipt)
            save(status_path, dict(status='PASS', bytes=offset, sha256=receipt['sha256'], updated_utc=utc().isoformat()))
            return receipt
        except urllib.error.HTTPError as error:
            if error.code not in (408, 429) and error.code < 500:
                save(status_path, dict(status='FAIL_ACCESS_OR_REQUEST_REVIEW', http_status=error.code,
                    error=str(error)[:200], updated_utc=utc().isoformat()))
                raise
            save(status_path, dict(status='RETRY_TRANSIENT', attempt=attempt, error=str(error)[:200], updated_utc=utc().isoformat()))
            if attempt < 4:
                time.sleep(min(15 * attempt, max(0, (deadline - utc()).total_seconds())))
        except (OSError, ValueError) as error:
            save(status_path, dict(status='RETRY_TRANSIENT', attempt=attempt, error=str(error)[:200], updated_utc=utc().isoformat()))
            if attempt < 4:
                time.sleep(min(15 * attempt, max(0, (deadline - utc()).total_seconds())))
        except BaseException as error:
            save(status_path, dict(status='FAIL_REVIEW_REQUIRED', error=str(error)[:200], updated_utc=utc().isoformat()))
            raise
    save(status_path, dict(status='FAIL_BOUNDED_RETRIES_EXHAUSTED', updated_utc=utc().isoformat()))
    raise RuntimeError('Bounded resource retries exhausted: ' + item['id'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=CODE / 'config/public-supplements-v1.json')
    parser.add_argument('--root', type=Path, default=Path('/srv/marina/resources/public-supplements-v1'))
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--only', help='One exact resource ID from the reviewed config')
    args = parser.parse_args(); os.umask(0o077)
    cfg = load(args.config)
    assert args.root.is_absolute() and args.root.name == cfg['bundle_id'] and not args.root.is_symlink(), 'Use the dedicated absolute bundle directory'
    assert all(Path(i['file']).name == i['file'] and not (args.root/i['file']).is_symlink() for i in cfg['resources'])
    assert not args.only or args.only in {i['id'] for i in cfg['resources']}
    assert sum(i['max_bytes'] for i in cfg['resources']) <= cfg['total_byte_ceiling'] <= 1000000000
    assert all(i['url'].startswith(('https://www.orphadata.com/data/xml/', 'https://cspec.genome.network/cspec/SequenceVariantInterpretation')) for i in cfg['resources'])
    if not args.execute:
        print(json.dumps({i['id']:load(args.root / (i['id'] + '-status.json')) or {'status':'MISSING'} for i in cfg['resources']}, indent=2)); return
    deadline = min(dt.datetime.fromisoformat(cfg['deadline_utc']), utc() + dt.timedelta(minutes=30))
    assert utc() < deadline
    args.root.mkdir(mode=0o700, parents=True, exist_ok=True)
    with (args.root / '.installer.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('ALREADY_RUNNING'); return
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            jobs = [pool.submit(download, i, args.root, deadline) for i in cfg['resources'] if not args.only or i['id'] == args.only]
            receipts = [job.result() for job in jobs]
        receipts = [load(args.root / (i['id'] + '-receipt.json')) for i in cfg['resources']]
        assert all(i.get('status') == 'PASS' for i in receipts), 'Not all reviewed sources are ready'
        save(args.root / 'source-manifest.json', dict(status='PASS_PUBLIC_SOURCES_ONLY',
            resources=receipts, no_patient_analysis=True, annotation_binding_unchanged=True, updated_utc=utc().isoformat()))
        print('PASS_PUBLIC_SOURCES_ONLY', flush=True)


if __name__ == '__main__':
    main()
