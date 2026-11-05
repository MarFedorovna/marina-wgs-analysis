#!/usr/bin/env bash
# One-command entry point. Status never starts a workload.
set -euo pipefail
export MARINA_CODE_ROOT
MARINA_CODE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec python3 - "$@" <<'PY'
import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

os.umask(0o077)
ROOT = Path('/srv/marina')
PRIVATE = ROOT / 'private'
CODE = Path(os.environ['MARINA_CODE_ROOT'])
RUNNER = CODE / 'scripts/run_parabricks.py'
TRANSFER = PRIVATE / 'pilot-inputs-verified.json'
READY = PRIVATE / 'pilot-ready.json'
CALLSET = 'chr1_diploid_technical_pilot'

def need(ok, message):
    if not ok:
        raise ValueError(message)

def public_error(error):
    # Filesystem errors can contain original patient FASTQ names.
    return str(error) if isinstance(error, ValueError) else 'A required private input is unavailable or malformed; inspect the private configuration on the server.'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def readiness():
    need(TRANSFER.is_file(), 'Waiting for the first FASTQ pair to pass SHA256 verification.')
    transfer = load(TRANSFER)
    need(transfer.get('SchemaVersion') == 1 and transfer.get('VerifiedFiles') == 2
         and transfer.get('Status') == 'PASS'
         and re.fullmatch(r'[0-9a-f]{64}', transfer.get('PlanId', ''))
         and transfer.get('VerifiedUtc'), 'Transfer verification receipt is incomplete.')
    need(READY.is_file(), 'Waiting for reviewed reference and private chr1 pilot configuration.')
    ready = load(READY)
    need(ready.get('status') == 'PASS' and ready.get('scope') == CALLSET,
         'Pilot readiness receipt has the wrong status or scope.')
    paths = {
        'run': PRIVATE / 'pilot-run.json', 'samples': PRIVATE / 'pilot-samples.tsv',
        'calls': PRIVATE / 'pilot-calls.tsv', 'runner': RUNNER,
    }
    for key, path in paths.items():
        need(path.is_file() and digest(path) == ready.get(key + '_sha256'),
             'Prepared ' + key + ' differs from the reviewed pilot receipt.')
    reference_receipt = Path(ready['reference_receipt'])
    need(reference_receipt.is_relative_to(ROOT / 'references') and reference_receipt.is_file()
         and digest(reference_receipt) == ready.get('reference_receipt_sha256'),
         'Reference verification receipt changed or is missing.')
    ref = load(reference_receipt)
    need(ref.get('schema_version') == 1 and ref.get('index_mode') == 'prebuilt'
         and ref.get('index_smoke_test') == 'PASS'
         and ref.get('validation', {}).get('sequence_md5_all_contigs') == 'PASS',
         'Reference preparation has not passed validation.')
    cfg = load(paths['run'])
    need(cfg.get('bqsr_enabled') is False and cfg.get('known_sites') == [],
         'This prepared technical pilot requires BQSR disabled.')
    need(Path(cfg['samples_manifest']) == paths['samples'] and Path(cfg['calls_manifest']) == paths['calls'],
         'Unexpected pilot manifests.')
    need(cfg['reference']['sha256'] == ready.get('fasta_sha256')
         and cfg['image_digest'] == ready.get('image_digest'), 'Reference or image identity changed.')
    samples, calls = rows(paths['samples']), rows(paths['calls'])
    need(len(samples) == len(calls) == 1, 'Expected one prepared sample and one call-plan row.')
    ids = {s['sample_id'] for s in samples}
    need(len(ids) == 1 and {c['sample_id'] for c in calls} == ids and transfer.get('SampleId') in ids,
         'Sample/call plan matching failed.')
    need(all(c['callset'] == CALLSET and c['ploidy'] == '2' for c in calls)
         and len({c['intervals'] for c in calls}) == 1, 'Only the prepared diploid chr1 callset is allowed.')
    reference_root = Path(cfg['reference_root']).resolve(strict=True)
    interval = (reference_root / calls[0]['intervals']).resolve(strict=True)
    need(interval.is_relative_to(reference_root) and digest(interval) == calls[0]['intervals_sha256']
         and all(c['intervals_sha256'] == calls[0]['intervals_sha256'] for c in calls),
         'Interval identity mismatch.')
    data = [line.split('\t') for line in interval.read_text().splitlines() if line and not line.startswith('@')]
    fasta = (reference_root / cfg['reference']['fasta']).resolve(strict=True)
    fai = [line.split('\t') for line in Path(str(fasta) + '.fai').read_text().splitlines()]
    chr1 = [r for r in fai if r[0] == 'chr1']
    need(len(chr1) == 1 and len(data) == 1 and data[0][:4] == ['chr1', '1', chr1[0][1], '+'],
         'Intervals must contain exactly the full chr1 autosome and no other contig.')
    transferred = {f['RelativePath']: f for f in transfer.get('Files', [])}
    need(len(transferred) == 2, 'Expected exactly two verified pilot FASTQ files.')
    for sample in samples:
        for mate in ('r1', 'r2'):
            path = (Path(cfg['reads_root']) / sample[mate]).resolve(strict=True)
            need(path.is_relative_to(Path(cfg['reads_root']).resolve()) and path.is_file() and path.stat().st_size > 0,
                 'A prepared FASTQ input is absent.')
            key = 'inputs/raw/' + sample[mate]
            need(key in transferred and path.stat().st_size == transferred[key]['SizeBytes']
                 and digest(path) == transferred[key]['Sha256'], 'Pilot FASTQ differs from the verified input receipt.')
    return paths['run'], samples[0]['sample_id']

parser = argparse.ArgumentParser(description='Status, dry-run plan, or manually start ONE sample technical chr1 pilot.')
parser.add_argument('action', nargs='?', choices=('status', 'plan', 'pilot'), default='status')
parser.add_argument('--wait', action='store_true', help='For pilot only: wait up to 60 minutes for verified inputs/config.')
args = parser.parse_args()
if args.wait and args.action != 'pilot':
    parser.error('--wait is only supported with pilot')
deadline = time.monotonic() + (3600 if args.wait else 0)
try:
    while True:
        try:
            config, first = readiness()
            break
        except (ValueError, OSError, KeyError, TypeError) as error:
            message = public_error(error)
            if args.action == 'status':
                print('NOT READY: ' + message)
                raise SystemExit(0)
            if not args.wait or time.monotonic() >= deadline:
                raise ValueError(message) from None
            print('WAITING: ' + message, flush=True)
            time.sleep(min(30, max(0, deadline - time.monotonic())))
    print('READY: verified first input pair, reference, and one-sample chr1 technical pilot plan.')
    print('Scope: first sample only; full input-pair alignment, then diploid chr1 gVCF. No final whole-WGS interpretation.')
    if args.action == 'status':
        state = PRIVATE / 'pilot-state.json'
        if state.is_file():
            print('Last launcher state: ' + load(state).get('status', 'unknown'))
        raise SystemExit(0)
    lock = (PRIVATE / '.pilot-launch.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise ValueError('Another launcher is already running; use status.') from None
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    log = PRIVATE / (args.action + '-' + stamp + '.log')
    mode = '--execute' if args.action == 'pilot' else '--dry-run'
    command = [sys.executable, str(RUNNER), '--config', str(config), '--stage', 'all', '--sample', first, mode]
    state = PRIVATE / 'pilot-state.json'
    def save_state(status):
        temporary = state.with_suffix('.json.tmp')
        temporary.write_text(json.dumps({'status': status, 'action': args.action,
            'scope': CALLSET, 'sample_count': 1, 'log': str(log),
            'updated_utc': dt.datetime.now(dt.timezone.utc).isoformat()}, indent=2) + '\n')
        temporary.replace(state)
    save_state('RUNNING' if args.action == 'pilot' else 'PLANNING')
    print('Private log: ' + str(log), flush=True)
    with log.open('x') as stream:
        result = subprocess.run(command, stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT)
    save_state('PASS' if result.returncode == 0 else 'FAIL')
    need(result.returncode == 0, 'Runner failed; inspect the private log on the server.')
    print('PASS: ' + ('one-sample technical pilot completed or verified by resume.' if args.action == 'pilot'
                      else 'dry-run plan validated; no analysis started.'))
except (ValueError, OSError, KeyError, TypeError) as error:
    print('ERROR: ' + public_error(error), file=sys.stderr)
    raise SystemExit(1)
PY
