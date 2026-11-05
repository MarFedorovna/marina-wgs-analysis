#!/usr/bin/env python3
"""Prepare a separate full-reference alignment queue after a completed one-sample pilot."""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tempfile


def need(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def read_samples(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        return reader.fieldnames, list(reader)


def write_exclusive(path, content):
    need(not path.is_symlink(), 'Refusing a generated-file symlink.')
    if path.exists():
        need(path.read_bytes() == content, 'Existing queue configuration differs; refusing overwrite.')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix='.remaining-', dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--execute', action='store_true', help='Write private configuration only; never launch analysis.')
    args = parser.parse_args()
    os.umask(0o077)
    root = args.root.resolve(strict=True)
    private = root / 'private'
    pilot = read_json(private / 'pilot-run.json')
    fields, canonical = read_samples(private / 'samples.tsv')
    _, selected = read_samples(Path(pilot['samples_manifest']))
    need(len(selected) == 1 and len(canonical) > 1, 'Expected a one-sample pilot and a larger canonical cohort.')
    first = selected[0]
    need(sum(row == first for row in canonical) == 1, 'Pilot metadata differs from the canonical manifest.')
    remaining = [row for row in canonical if row['sample_id'] != first['sample_id']]
    need(len({row['sample_id'] for row in canonical}) == len(canonical), 'Duplicate cohort sample identifiers.')
    completed = Path(pilot['output_root']) / first['sample_id'] / 'align'
    receipt = read_json(completed / 'receipt.json')
    need(receipt.get('sample') == first['sample_id'] and receipt.get('config') == pilot,
         'Completed pilot alignment provenance differs.')
    for name, expected in receipt['outputs'].items():
        path = (completed / name).resolve(strict=True)
        need(path.is_relative_to(completed.resolve()) and digest(path) == expected,
             'An existing pilot alignment artifact changed.')
    transfer = read_json(private / 'transfer-verified.json')
    need(transfer.get('Failed') == 0 and transfer.get('VerifiedFiles', 0) > 0,
         'The complete input transfer has not passed verification.')
    plan_id = transfer.get('PlanId', '')
    need(re.fullmatch(r'[0-9a-f]{64}', plan_id), 'Invalid private transfer plan identity.')
    expected_path = root / '.transfer' / (plan_id + '.json')
    need(digest(expected_path) == plan_id, 'Expected-source transfer manifest changed.')
    expected = {row['RelativePath']: row for row in read_json(expected_path)['Files']}
    reads = Path(pilot['reads_root']).resolve(strict=True)
    for row in canonical:
        for mate in ('r1', 'r2'):
            path = (reads / row[mate]).resolve(strict=True)
            need(path.is_relative_to(reads), 'An input escapes reads_root.')
            source = expected[path.relative_to(root).as_posix()]
            need(path.stat().st_size == source['SizeBytes'] and digest(path) == source['Sha256'],
                 'A read file differs from its expected source.')
    manifest = private / 'remaining-align-samples.tsv'
    config = private / 'remaining-align.json'
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(remaining)
    cfg = dict(pilot)
    cfg.update(samples_manifest=str(manifest), output_root=str(root / 'outputs/remaining-align'),
               scratch_root=str(root / 'scratch/remaining-align'), max_run_hours=4)
    cfg['processing_assumptions'] = [text for text in pilot['processing_assumptions']
                                   if not text.startswith(('Limited technical pilot:', 'ONE first canonical sample only:'))
                                   and 'final all-sample queue' not in text
                                   and 'final four-sample queue' not in text]
    cfg['processing_assumptions'].append(
        'Alignment-only queue for the remaining canonical samples against the complete fixed reference; '
        'the completed first alignment remains in its original verified pilot directory. '
        'No new variant calls or chromosome/ploidy policy are selected by this configuration.')
    cfg['calls_manifest'] = str(private / 'NOT_CONFIGURED_FOR_ALIGNMENT_ONLY.tsv')
    artifacts = {manifest: stream.getvalue().encode(),
                 config: (json.dumps(cfg, indent=2) + '\n').encode()}
    for path, content in artifacts.items():
        need(not path.is_symlink() and (not path.exists() or path.read_bytes() == content),
             'An existing private queue configuration differs; refusing overwrite.')
    if args.execute:
        for path, content in artifacts.items():
            write_exclusive(path, content)
    print(json.dumps({'status': 'PREPARED' if args.execute else 'DRY_RUN',
                      'cohort_samples': len(canonical), 'existing_alignments_preserved': 1,
                      'new_alignments': len(remaining), 'scope': 'full-reference alignment only',
                      'analysis_started': False}))


if __name__ == '__main__':
    main()
