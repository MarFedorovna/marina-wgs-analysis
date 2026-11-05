#!/usr/bin/env python3
"""Review-only germline X/Y plan. There is deliberately no patient-calling executor."""
import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import re

CODE = Path(__file__).resolve().parents[1]
COLUMNS = ['sample_id', 'karyotype', 'x_nonpar_ploidy', 'y_nonpar_ploidy', 'x_par_ploidy', 'evidence_source']


def need(ok, message):
    if not ok:
        raise ValueError(message)


def load(path):
    return json.loads(Path(path).read_text())


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 ** 2), b''):
            result.update(block)
    return result.hexdigest()


def tsv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def complement(length, ranges):
    result, start = [], 1
    for left, right in sorted(ranges):
        need(1 <= left <= right <= length and left >= start, 'Overlapping/out-of-bounds PAR intervals.')
        if start < left:
            result.append([start, left - 1])
        start = right + 1
    if start <= length:
        result.append([start, length])
    return result


def partitions(lock):
    par, lengths = lock['par_one_based_closed'], lock['primary_contig_lengths']
    return {'X_NONPAR': [('chrX', a, b) for a, b in complement(lengths['chrX'], par['chrX'])],
            'Y_NONPAR': [('chrY', a, b) for a, b in complement(lengths['chrY'], par['chrY'])],
            'X_PAR': [('chrX', a, b) for a, b in par['chrX']]}


def count_range(fasta, row, start, end):
    length, offset, bases, width = map(int, row[1:5])
    need(1 <= start <= end <= length and bases > 0 and width >= bases, 'Invalid FAI region.')
    position = lambda n: offset + ((n - 1) // bases) * width + ((n - 1) % bases)
    with fasta.open('rb') as stream:
        left, right = position(start), position(end)
        stream.seek(left)
        sequence = stream.read(right - left + 1).replace(b'\n', b'').replace(b'\r', b'').upper()
    need(len(sequence) == end - start + 1, 'FASTA/FAI lengths differ.')
    return {'bases': len(sequence), 'N': sequence.count(b'N'),
            'ACGT': sum(sequence.count(bytes([letter])) for letter in b'ACGT')}


def reference_preflight(cfg, lock):
    fasta = Path(cfg['reference']['fasta']).resolve(strict=True)
    need(digest(fasta) == cfg['reference']['sha256'], 'Analysis FASTA SHA mismatch.')
    fai = Path(str(fasta) + '.fai')
    rows = {row[0]: row for row in csv.reader(fai.open(), delimiter='\t')}
    need(all(chrom in rows and int(rows[chrom][1]) == length for chrom, length in lock['primary_contig_lengths'].items()),
         'Expected GRCh38 chrX/chrY lengths are absent; do not use GRCh37 or another assembly.')
    y = [count_range(fasta, rows['chrY'], start, end) for start, end in lock['par_one_based_closed']['chrY']]
    return {'status': 'PASS_REFERENCE_IDENTITY', 'fasta_sha256': cfg['reference']['sha256'], 'fai_sha256': digest(fai),
        'Y_PAR_counts': y, 'Y_PAR_all_hard_masked_N': all(item['N'] == item['bases'] for item in y),
        'mask_at_original_alignment_not_proven_by_current_FASTA_alone': True}


def sample_plan(row, blockers):
    reasons = list(blockers)
    if row.get('karyotype') not in ('XX', 'XY'):
        reasons.append('unknown_or_unsupported_karyotype; technical_XY_ratios_are_not_confirmation')
    if row.get('evidence_source', 'unknown') in ('', 'unknown'):
        reasons.append('copy_number_evidence_not_documented')
    values = [row.get(name, 'unknown') for name in COLUMNS[2:5]]
    if any(value not in ('0', '1', '2') for value in values):
        reasons.append('undefined_ploidy; no_diploid_default')
    expected = {'XX': ['2', '0', '2'], 'XY': ['1', '1', '2']}
    if row.get('karyotype') in expected and values != expected[row['karyotype']]:
        reasons.append('explicit_ploidy_inconsistent_with_supported_model; specialized_review_required')
    result = []
    for name, value in zip(('X_NONPAR', 'Y_NONPAR', 'X_PAR'), values):
        if name == 'Y_NONPAR' and value == '0':
            result.append({'partition': name, 'declared_ploidy': 0, 'status': 'SKIP_DECLARED_ZERO_COPY_NOT_REFERENCE',
                'genotype': 'not_assessed_no_GT; never_0/0_or_homref_gVCF', 'review_required': bool(reasons), 'command_argv': None})
        else:
            result.append({'partition': name, 'declared_ploidy': int(value) if value in ('1', '2') else None,
                'status': 'BLOCKED_BEFORE_CALLING' if reasons else 'READY_FOR_COMMAND_REVIEW_ONLY',
                'blockers': sorted(set(reasons)), 'command_argv': None})
    return result


def make_plan(cfg, rows, contract=None, preflight=None):
    lock = load(CODE / 'config/sex-chromosomes.resources.lock.json')
    samples = [row['sample_id'] for row in tsv(cfg['samples_manifest'])]
    need(samples and len(set(samples)) == len(samples) and all(re.fullmatch(r'[A-Za-z0-9_.-]{1,100}', name) for name in samples), 'Invalid sample manifest.')
    need(len(rows) == len(samples) and {row['sample_id'] for row in rows} == set(samples), 'Ploidy table must contain each cohort sample exactly once.')
    alignments = {row['sample_id']: row for row in tsv(cfg['alignment_manifest'])}
    need(set(alignments) == set(samples), 'Alignment manifest/cohort mismatch.')
    blockers = []
    if not preflight or preflight.get('status') != 'PASS_REFERENCE_IDENTITY':
        blockers.append('reference_identity_and_Y_PAR_mask_not_preflight_verified')
    elif not preflight['Y_PAR_all_hard_masked_N']:
        blockers.append('Y_PAR_not_hard_masked; reviewed_remapping_reference_alignment_contract_required')
    contract = contract or {}
    if not (contract.get('status') == 'PASS' and contract.get('reviewer') not in (None, '', 'SET_REVIEWER')
            and contract.get('reference_sha256') == cfg['reference']['sha256']
            and contract.get('official_par_gff_sha256') == lock['retrieved_gff_sha256']
            and contract.get('par_policy') == 'Y_PAR_HARD_MASKED_AT_ALIGNMENT'
            and contract.get('mask_present_at_original_alignment') is True
            and contract.get('par_report_coordinates') == 'X_only'):
        blockers.append('PAR_reference_alignment_contract_missing_or_unverified; X_only_diploid_is_not_sufficient')
    image = cfg['tools']['gatk_image']
    need(re.fullmatch(r'[^\s@]+@sha256:[0-9a-f]{64}', image), 'Reuse an immutable reviewed GATK image.')
    planned = []
    for row in rows:
        sample = row['sample_id']
        entry = alignments[sample]
        receipt_path = Path(entry['receipt'])
        sample_blockers = list(blockers)
        saved = load(receipt_path)
        need(saved.get('sample') == sample and saved.get('config', {}).get('reference', {}).get('sha256') == cfg['reference']['sha256'],
             'Original alignment receipt/reference mismatch.')
        if contract.get('sample_alignment_receipt_sha256', {}).get(sample) != digest(receipt_path):
            sample_blockers.append('sample_original_alignment_receipt_not_bound_to_PAR_contract')
        parts = sample_plan(row, sample_blockers)
        for item in parts:
            # Command is a JSON review artifact only. This program never invokes Docker/GATK.
            item['command_template'] = ['docker', 'run', '--rm', '--network', 'none', '--cpus', '4', '--memory', '16g',
                '--mount', 'type=bind,src=SET_VERIFIED_REFERENCE_DIR,dst=/ref,readonly',
                '--mount', 'type=bind,src=SET_VERIFIED_BAM_DIR,dst=/bam,readonly',
                '--mount', 'type=bind,src=SET_NEW_XY_OUTPUT_DIR,dst=/out', image,
                'gatk', '--java-options', '-Xmx12g', 'HaplotypeCaller', '-R', '/ref/SET_FASTA', '-I', '/bam/SET_VERIFIED_BAM',
                '-L', '/out/' + item['partition'] + '.intervals', '--interval-padding', '0',
                '--sample-ploidy', str(item['declared_ploidy']) if item['declared_ploidy'] else 'PLOIDY_REQUIRED_OR_SKIP',
                '--native-pair-hmm-threads', '4', '-ERC', 'GVCF', '-O', '/out/SET_SAMPLE.' + item['partition'] + '.g.vcf.gz']
            if item['status'].startswith('SKIP'):
                item['command_template'] = None
        planned.append({'sample_id': sample, 'karyotype_declared': row['karyotype'], 'evidence_source': row['evidence_source'], 'partitions': parts})
    return {'schema_version': 1, 'status': 'REVIEW_ONLY_NO_CALLING_EXECUTOR',
        'sample_count': len(samples), 'any_calling_started': False, 'sex_inferred': False,
        'scope': 'germline_nuclear_XY; not somatic; no MT or autosomal changes',
        'reference': cfg['reference'], 'reference_preflight': preflight or {'status': 'NOT_REQUESTED'},
        'PAR_source': lock, 'partitions_one_based_closed': partitions(lock), 'samples': planned,
        'global_blockers': blockers, 'ploidy_table_sha256': hashlib.sha256(json.dumps(rows, sort_keys=True).encode()).hexdigest(),
        'tools': {'gatk_image': image, 'GATK_baseline': '4.3.0.0', 'cpus': 4, 'memory_gib': 16, 'gpu': False},
        'budget': dict(cfg.get('budget', {}), note='Existing absolute VM allowance must be rechecked by a future executor; no new allowance is created.'),
        'limitations': ['No BAM realignment/masking is performed.', 'No genotype is manufactured for zero-copy or unknown regions.',
            'Joint genotyping/merging mixed-ploidy X/Y and clinical validation are not implemented by this planner.',
            'X-transposed regions, segmental homology, CNVs/aneuploidy/mosaicism and low coverage require separate review.']}


def write_plan(target, plan, rows):
    need(not target.exists(), 'Refusing to overwrite an existing private X/Y plan.')
    target.mkdir(parents=True, mode=0o700)
    (target / 'plan.json').write_text(json.dumps(plan, indent=2, ensure_ascii=False) + '\n')
    with (target / 'ploidy-review.tsv').open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    for name, ranges in plan['partitions_one_based_closed'].items():
        # Picard interval_list needs full sequence headers; use unambiguous GATK text intervals instead.
        (target / (name + '.intervals')).write_text(''.join(f'{chrom}:{start}-{end}\n' for chrom, start, end in ranges))
        (target / (name + '.bed')).write_text(''.join(f'{chrom}\t{start-1}\t{end}\n' for chrom, start, end in ranges))
    (target / 'NOT_EXECUTED.txt').write_text('Review artifacts only. No calls, realignment or cohort merge executed.\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--ploidy-manifest', type=Path)
    parser.add_argument('--alignment-contract', type=Path)
    parser.add_argument('--preflight', action='store_true', help='Read-only FASTA SHA/FAI and Y-PAR hard-mask check; no caller.')
    parser.add_argument('--write-plan', type=Path, help='Write a new private review directory; never execute commands.')
    args = parser.parse_args()
    os.umask(0o077)
    cfg = load(args.config)
    lock = load(CODE / 'config/sex-chromosomes.resources.lock.json')
    rows = tsv(args.ploidy_manifest) if args.ploidy_manifest else [dict.fromkeys(COLUMNS, 'unknown') | {'sample_id': row['sample_id']} for row in tsv(cfg['samples_manifest'])]
    plan = make_plan(cfg, rows, load(args.alignment_contract) if args.alignment_contract else None,
        reference_preflight(cfg, lock) if args.preflight else None)
    if args.write_plan:
        write_plan(args.write_plan, plan, rows)
    print(json.dumps({'status': plan['status'], 'sample_count': plan['sample_count'], 'calling_started': False,
        'blocked_partitions': sum(item['status'] == 'BLOCKED_BEFORE_CALLING' for sample in plan['samples'] for item in sample['partitions']),
        'global_blockers': plan['global_blockers'], 'output': str(args.write_plan) if args.write_plan else None,
        'reference_preflight': plan['reference_preflight']}, ensure_ascii=False, indent=2))
