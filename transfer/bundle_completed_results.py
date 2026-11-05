#!/usr/bin/env python3
"""Build one server-side small-results TAR from completed receipts and frozen metadata.

Run on the server. Excludes BAM, raw reads, individual gVCF, public databases,
uncompressed/intermediate VCF and credentials. No patient-file hashes are computed.
The expected manifest should reach the receiver through pinned SSH separately.
"""
import argparse
import contextlib
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import sys
import tarfile
import uuid
from prepare_native_plan import HELPER
from pull_verified_results import REMOTE_SNAPSHOT


def run_helper(code, config):
    previous = sys.argv
    try:
        sys.argv = ['helper', json.dumps(config)]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            exec(compile(code, '<receipt-helper>', 'exec'), {})
        return json.loads(output.getvalue())
    finally:
        sys.argv = previous


def selected(relative):
    name = PurePosixPath(relative).name.lower()
    parts = PurePosixPath(relative).parts
    if len(parts)>1 and parts[1].startswith('mitochondria'):
        return parts[1]=='mitochondria-recovery-v2' and len(parts)==6 and parts[-2]=='calls' and name in {
            'mitochondria.filtered.vcf.gz','mitochondria.filtered.vcf.gz.tbi',
            'mitochondria.pass.vcf.gz','mitochondria.pass.vcf.gz.tbi',
            'mitochondria.raw.vcf.gz','mitochondria.raw.vcf.gz.tbi',
            'research-report.json','numt-assessment.json','contamination.json',
            'heteroplasmy.tsv','circular-coverage.tsv'}
    if name.endswith('.bam') or '.g.vcf' in name:
        return False
    if '/checkpoints/' in relative and any(s in name for s in ('.vcf', '.bcf', '.tbi', '.csi')):
        return False
    if '.vcf' in name or name.endswith(('.bcf', '.tbi', '.csi')):
        return name.startswith(('cohort.filtered.vcf.gz', 'cohort.quality.vcf.gz'))
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--remote-root', required=True)
    ap.add_argument('--workspace-root', required=True)
    ap.add_argument('--jupyter-root', required=True)
    ap.add_argument('--code-root', required=True)
    ap.add_argument('--execute', action='store_true')
    ap.add_argument('--max-mib', type=int, default=200)
    args = ap.parse_args()
    if not args.execute:
        print('DRY RUN: no snapshots or archive written')
        return
    os.umask(0o077)
    cfg = {'root': args.remote_root, 'workspace': args.workspace_root,
           'jupyter': args.jupyter_root, 'code': args.code_root}
    native = run_helper(HELPER, cfg)
    metadata = run_helper(REMOTE_SNAPSHOT, cfg)['Manifest']
    entries = {e['RelativePath']: e for e in native['Files'] if selected(e['RelativePath'])}
    # Mutable private configs/logs are frozen by REMOTE_SNAPSHOT; raw sources are never archived.
    for entry in metadata['Files']:
        rel = entry['RelativePath']
        if rel.startswith('metadata/') or rel.endswith(('receipt.json', 'workflow-receipt.json')):
            entries[rel] = dict(entry, ServerPath=entry['Source'])
    total = sum(e['SizeBytes'] for e in entries.values())
    if total > args.max_mib * 2**20:
        raise RuntimeError('Selected bundle exceeds bounded size; review selection before transfer')
    name = 'completed-results-' + dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8] + '.tar'
    target = Path(args.workspace_root) / 'marina-verified-downloads' / name
    temporary = target.with_suffix('.tar.partial')
    plan = {'SchemaVersion': 1, 'Verification': 'size_only', 'PatientFileHashesComputed': False,
            'GeneratedUtc': dt.datetime.now(dt.timezone.utc).isoformat(), 'Files': []}
    with temporary.open('xb') as stream, tarfile.open(fileobj=stream, mode='w') as archive:
        for index, (rel, entry) in enumerate(sorted(entries.items())):
            p = Path(entry['ServerPath'])
            if p.is_symlink() or not p.is_file() or p.stat().st_size != entry['SizeBytes']:
                raise RuntimeError('Source changed or unsafe')
            before = p.stat()
            member = f'files/{index:05d}'
            info = tarfile.TarInfo(member)
            info.size, info.mode, info.mtime = before.st_size, 0o600, int(before.st_mtime)
            with p.open('rb') as source:
                archive.addfile(info, source)
            after = p.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise RuntimeError('Source changed during bundle creation')
            plan['Files'].append({'RelativePath': rel, 'SizeBytes': before.st_size,
                                  'Sha256': entry['Sha256'], 'ArchiveMember': member})
    os.link(temporary, target)  # atomic no-replace; never overwrite an earlier download.
    temporary.unlink()
    plan.update({'DownloadFileName': name, 'SizeBytes': target.stat().st_size,
                 'JupyterRelativePath': target.relative_to(Path(args.jupyter_root)).as_posix(),
                 'PayloadBytes': total, 'Count': len(entries)})
    manifest = target.with_suffix('.manifest.json')
    with manifest.open('x') as stream:
        json.dump(plan, stream, indent=2)
    print(json.dumps(plan))


if __name__ == '__main__':
    main()
