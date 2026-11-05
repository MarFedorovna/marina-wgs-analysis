#!/usr/bin/env bash
# Standalone sidecar only. Never starts/restarts an annotation workflow.
set -euo pipefail
umask 077
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 - "$SCRIPT_DIR" "$@" <<'PY'
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

script_dir = Path(sys.argv[1])
parser = argparse.ArgumentParser(description='Require completed stage 6 AND 7 before separate consequence export.')
parser.add_argument('--workflow-root', required=True, help='Directory containing stages/06_annotate and stages/07_report')
parser.add_argument('--output', required=True, help='New private sidecar directory outside workflow stages')
parser.add_argument('--wait-seconds', type=int, default=0, help='Bounded wait; default is fail immediately when not ready')
parser.add_argument('--synthetic', action='store_true', help='Label synthetic control exports explicitly')
args = parser.parse_args(sys.argv[2:])
if args.wait_seconds < 0 or args.wait_seconds > 43200:
    raise SystemExit('wait-seconds must be 0..43200')
root = Path(args.workflow_root).resolve()
output = Path(args.output).resolve()
stages = root / 'stages'
if output == root or root in output.parents:
    raise SystemExit('Output must be outside the workflow root.')
if output.exists():
    raise SystemExit('Output already exists; refusing overwrite.')
source = stages / '06_annotate/transcript_samples.tsv.gz'
receipt_paths = [stages / name / 'workflow-receipt.json' for name in ('06_annotate', '07_report')]
deadline = time.monotonic() + args.wait_seconds
while True:
    try:
        receipts = [json.loads(path.read_text()) for path in receipt_paths]
    except FileNotFoundError:
        receipts = []
    if len(receipts) == 2 and all(r.get('status') == 'PASS' for r in receipts) and source.is_file():
        break
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise SystemExit('NOT_READY: require stage 6 and 7 PASS receipts and source master. No export started.')
    time.sleep(min(15, remaining))
expected = receipts[0].get('artifacts', {}).get('transcript_samples.tsv.gz')
digest = hashlib.sha256()
with source.open('rb') as stream:
    for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
        digest.update(block)
expected_hash = expected.get('sha256') if isinstance(expected, dict) else None
if not expected_hash or digest.hexdigest() != expected_hash:
    raise SystemExit('Source master does not match completed stage 6 artifact SHA256.')
command = [sys.executable, str(script_dir / 'export_consequence_tables.py'),
           '--source', str(source), '--receipt', str(receipt_paths[0]), '--output', str(output)]
if args.synthetic:
    command.append('--synthetic')
subprocess.run(command, check=True)
provenance = {'stage6_receipt': receipts[0], 'stage7_receipt': receipts[1],
              'note': 'Sidecar export; existing workflow outputs unchanged.'}
(output / 'continuation-receipts.json').write_text(json.dumps(provenance, indent=2) + '\n')
extra = output / 'continuation-receipts.json'
manifest_path = output / 'manifest.json'
manifest = json.loads(manifest_path.read_text())
manifest['artifacts'][extra.name] = {'bytes': extra.stat().st_size,
                                    'sha256': hashlib.sha256(extra.read_bytes()).hexdigest()}
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
PY
