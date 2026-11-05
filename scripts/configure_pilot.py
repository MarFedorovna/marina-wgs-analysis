#!/usr/bin/env python3
"""Prepare private ONE-sample chr1 technical-pilot configuration; never run analysis."""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import re
import sys
import tempfile

CALLSET = 'chr1_diploid_technical_pilot'
TAG = 'nvcr.io/nvidia/clara/clara-parabricks:4.8.0-1'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def encoded_json(value):
    return (json.dumps(value, indent=2) + '\n').encode()


def tsv(rows, fields):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def prepare(args):
    root = args.root.resolve(strict=True)
    refs = root / 'references'
    private = (args.private_dir or root / 'private').resolve()
    receipt_path = args.reference_receipt.resolve(strict=True)
    require(receipt_path.is_relative_to(refs.resolve(strict=True)), 'Reference receipt must be below root/references.')
    receipt = json.loads(receipt_path.read_text())
    require(receipt.get('schema_version') == 1 and receipt.get('index_mode') == 'prebuilt'
            and receipt.get('index_smoke_test') == 'PASS'
            and receipt.get('validation', {}).get('sequence_md5_all_contigs') == 'PASS',
            'A validated prebuilt reference receipt is required.')
    image = receipt.get('image', '')
    require(re.fullmatch(r'nvcr\.io/nvidia/clara/clara-parabricks@sha256:[0-9a-f]{64}', image),
            'Reference receipt must identify the pinned immutable Parabricks image.')
    artifacts = receipt['artifacts']
    candidates = [a for a in artifacts if Path(a['path']).suffix in ('.fa', '.fasta', '.fna')]
    require(len(candidates) == 1, 'Expected one uncompressed FASTA in the reference receipt.')
    artifact = candidates[0]
    fasta = (receipt_path.parent / artifact['path']).resolve(strict=True)
    require(fasta.is_relative_to(refs.resolve()) and fasta.stat().st_size == artifact['bytes']
            and sha256(fasta) == artifact['sha256'], 'FASTA does not match its verified reference receipt.')
    expected = {a['path']: a for a in artifacts}
    paths = [Path(str(fasta) + s) for s in ('.fai', '.amb', '.ann', '.bwt', '.pac', '.sa', '.alt')]
    paths += [fasta.with_suffix('.dict')]
    for path in paths:
        rel = path.relative_to(receipt_path.parent).as_posix()
        require(path.is_file() and rel in expected and path.stat().st_size == expected[rel]['bytes']
                and sha256(path) == expected[rel]['sha256'], 'A required reference index does not match its receipt.')
    with args.samples.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames
        samples = list(reader)
    columns = ['sample_id', 'r1', 'r2', 'read_group_id', 'library', 'library_origin', 'platform', 'platform_unit']
    require(fields == columns and samples and all(set(s) == set(columns) and all(v is not None for v in s.values()) for s in samples),
            'Canonical samples must use the documented eight-column TSV contract.')
    require(len({s['sample_id'] for s in samples}) == len(samples), 'Duplicate canonical sample IDs.')
    selected = samples[0]
    for key in ('sample_id', 'read_group_id', 'library'):
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,99}', selected[key]), 'Invalid selected sample metadata.')
    require(selected['library_origin'] in ('confirmed', 'processing_assumption'), 'Set library_origin explicitly.')
    for mate in ('r1', 'r2'):
        relative = Path(selected[mate])
        require(not relative.is_absolute() and '..' not in relative.parts and selected[mate], 'FASTQ paths must be relative to reads_root.')
    require(selected['r1'] != selected['r2'], 'R1/R2 must be distinct.')
    # FASTQs/transfer marker may still be in transit. The launcher verifies them before analysis.
    fai = [line.split('\t') for line in Path(str(fasta) + '.fai').read_text().splitlines()]
    chr1 = [line for line in fai if line[0] == 'chr1']
    require(len(chr1) == 1 and int(chr1[0][1]) > 0, 'Reference must contain exactly one chr1 autosome.')
    require(re.fullmatch(r'[A-Za-z0-9_.-]+\.interval_list', args.interval_name), 'Invalid interval filename.')
    interval_path = refs / 'analysis-plans' / args.interval_name
    headers = [line for line in fasta.with_suffix('.dict').read_text().splitlines() if line.startswith('@')]
    interval = ('\n'.join(headers) + '\nchr1\t1\t' + chr1[0][1] + '\t+\tchr1_technical_pilot\n').encode()
    interval_sha = hashlib.sha256(interval).hexdigest()
    sample_path, call_path, run_path = (private / name for name in ('pilot-samples.tsv', 'pilot-calls.tsv', 'pilot-run.json'))
    call = dict(sample_id=selected['sample_id'], callset=CALLSET,
                intervals=interval_path.relative_to(refs).as_posix(), intervals_sha256=interval_sha, ploidy='2')
    code = Path(__file__).resolve().parent.parent
    runner = code / 'scripts/run_parabricks.py'
    cfg = json.loads((code / 'config/run.example.json').read_text(encoding='utf-8-sig'))
    cfg.pop('_instructions', None)
    cfg.update(image=TAG, image_digest=image, reads_root=str(root / 'inputs/raw'), reference_root=str(refs),
               output_root=str(root / 'outputs/chr1-technical-pilot'), scratch_root=str(root / 'scratch/chr1-technical-pilot'),
               samples_manifest=str(sample_path), calls_manifest=str(call_path), bqsr_enabled=False, known_sites=[])
    cfg['reference'] = dict(fasta=fasta.relative_to(refs).as_posix(), sha256=artifact['sha256'],
        assembly_label=receipt['assembly_label'],
        source=json.dumps([{'bucket': s['bucket'], 'object': s['object'], 'generation': s['generation']}
                           for s in receipt['sources']], sort_keys=True))
    cfg['processing_assumptions'] = [
        'ONE first canonical sample only: full FASTQ-pair alignment, then chr1 autosome-only diploid gVCF; not final whole-WGS interpretation.',
        'This is a new uniform reference verified by the reference receipt, not a claim to reproduce the earlier laboratory bundle.',
        'A library marked processing_assumption is technical grouping of one available FASTQ pair; original wet-lab library provenance remains unknown.',
        'Empty PL remains unknown and is omitted. Empty physical PU remains unknown; runner uses a technical analysis_unit_<read_group_id> label, not a confirmed flowcell.',
        'BQSR is disabled and no known-sites are supplied. Duplicates are marked rather than deleted.',
        'No sex inference, sex-chromosome/MT ploidy policy, CNV policy, clinical interpretation, or final four-sample queue is defined.'
    ]
    outputs = {interval_path: interval, sample_path: tsv([selected], columns),
               call_path: tsv([call], list(call)), run_path: encoded_json(cfg)}
    if args.transfer_manifest:
        transfer_bytes = args.transfer_manifest.read_bytes()
        transfer = json.loads(transfer_bytes)
        require(transfer.get('SchemaVersion') == 1 and isinstance(transfer.get('Files'), list),
                'Expected a SchemaVersion 1 server upload manifest with Files.')
        files = {item['RelativePath']: item for item in transfer['Files']}
        require(len(files) == len(transfer['Files']), 'Duplicate relative paths in the upload manifest.')
        verified = []
        for mate in ('r1', 'r2'):
            relative = 'inputs/raw/' + selected[mate]
            require(relative in files, 'Selected FASTQ is absent from the expected-source upload manifest.')
            expected_file = files[relative]
            path = (root / relative).resolve(strict=True)
            require(path.is_relative_to((root / 'inputs/raw').resolve(strict=True))
                    and path.stat().st_size == expected_file['SizeBytes']
                    and re.fullmatch(r'[0-9a-f]{64}', expected_file.get('Sha256', ''))
                    and sha256(path) == expected_file['Sha256'],
                    'Selected FASTQ does not match the expected source size/SHA256.')
            verified.append({k: expected_file[k] for k in ('RelativePath', 'SizeBytes', 'Sha256')})
        marker = dict(SchemaVersion=1, Status='PASS', VerifiedFiles=2, SampleId=selected['sample_id'],
                      PlanId=hashlib.sha256(transfer_bytes).hexdigest(),
                      VerifiedUtc=dt.datetime.now(dt.timezone.utc).isoformat(), Files=verified)
        marker_path = private / 'pilot-inputs-verified.json'
        if marker_path.exists():
            existing = json.loads(marker_path.read_text())
            require(all(existing.get(k) == v for k, v in marker.items() if k not in ('VerifiedUtc', 'Files'))
                    and sorted(existing.get('Files', []), key=lambda x: x['RelativePath'])
                    == sorted(verified, key=lambda x: x['RelativePath']) and existing.get('VerifiedUtc'),
                    'Existing pilot input marker belongs to a different expected-source plan; refusing overwrite.')
            marker = existing  # preserve its verified timestamp and bytes
        outputs[marker_path] = encoded_json(marker) if not marker_path.exists() else marker_path.read_bytes()
    ready = dict(status='PASS', scope=CALLSET, sample_count=1, canonical_sample_count=len(samples),
                 reference_receipt=str(receipt_path), reference_receipt_sha256=sha256(receipt_path),
                 fasta_sha256=artifact['sha256'], image_digest=image, runner_sha256=sha256(runner))
    for name, path in [('run', run_path), ('samples', sample_path), ('calls', call_path)]:
        ready[name + '_sha256'] = hashlib.sha256(outputs[path]).hexdigest()
    outputs[private / 'pilot-ready.json'] = encoded_json(ready)
    # Check every collision before writing anything. Keep an existing changed plan for review.
    for path, content in outputs.items():
        require(not path.is_symlink(), 'Refusing a symlink in a generated configuration destination.')
        require(not path.exists() or path.read_bytes() == content, 'An existing generated artifact differs; use a separate --private-dir/--interval-name for review.')
    if args.execute:
        for path, content in outputs.items():
            if path.exists():
                continue
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix='.configure-', dir=path.parent)
            temporary = Path(name)
            try:
                with os.fdopen(fd, 'wb') as stream:
                    stream.write(content)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.link(temporary, path)  # atomic, exclusive: never replace a concurrently created file
            finally:
                temporary.unlink(missing_ok=True)
    print(json.dumps(dict(status='PREPARED' if args.execute else 'DRY_RUN', scope=CALLSET,
                         selected_samples=1, canonical_samples=len(samples), generated_artifacts=len(outputs),
                         chr1_bases=int(chr1[0][1]), analysis_started=False)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--samples', type=Path, required=True, help='Canonical private eight-column samples TSV; never modified.')
    parser.add_argument('--reference-receipt', type=Path, required=True)
    parser.add_argument('--transfer-manifest', type=Path,
                        help='Expected-source SchemaVersion 1 upload manifest; verify first pair and create the launcher input marker.')
    parser.add_argument('--private-dir', type=Path, help='Optional separate destination for private review/validation.')
    parser.add_argument('--interval-name', default='chr1-diploid-technical-pilot.interval_list')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--dry-run', action='store_true', help='Default: validate and plan without any writes.')
    mode.add_argument('--execute', action='store_true', help='Create private configuration only; never start analysis.')
    args = parser.parse_args()
    os.umask(0o077)
    try:
        prepare(args)
    except (ValueError, OSError, KeyError, TypeError) as error:
        message = str(error) if isinstance(error, ValueError) else 'A required input is unavailable or malformed; inspect private inputs on the server.'
        print('ERROR: ' + message, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
