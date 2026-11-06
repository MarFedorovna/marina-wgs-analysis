"""Patient-independent GENCODE49 exonic XY controls and common X SNP panel."""
import argparse
from collections import Counter, defaultdict
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import random
import signal
import time
import pysam

PAR = {'chrX': [(10000, 2781479), (155701382, 156030895)],
       'chrY': [(10000, 2781479), (56887902, 57217415)]}
CHROMS = ['chr' + str(n) for n in range(1, 23)] + ['chrX', 'chrY']


def pin(path):
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def merge(regions):
    result = []
    for a, b in sorted(regions):
        if result and a <= result[-1][1]:
            result[-1] = (result[-1][0], max(b, result[-1][1]))
        else:
            result.append((a, b))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['gtf', 'cache', 'reference-directory', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    def expired(*unused):
        raise TimeoutError('Public target panel preparation max300s')
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(300)
    began = time.monotonic()
    args.output.mkdir(mode=0o700, parents=True, exist_ok=False)
    counts = Counter()
    intervals = defaultdict(list)
    headers = []
    with gzip.open(args.gtf, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                headers.append(line.rstrip())
                continue
            row = line.rstrip().split('\t')
            if row[0] not in CHROMS or row[2] != 'exon' or 'gene_type "protein_coding";' not in row[8]:
                continue
            a, b = int(row[3]) - 1, int(row[4])
            if any(a < end and b > start for start, end in PAR.get(row[0], [])):
                counts['PAR_overlapping_exon_rows_excluded'] += 1
                continue
            intervals[row[0]].append((a, b))
    if not any('version 49 (Ensembl 115)' in line and 'GRCh38' in line for line in headers):
        raise ValueError('Exact public GENCODE49 GRCh38 header required')
    intervals = {c: merge(intervals[c]) for c in CHROMS}
    reference = args.reference_directory / 'Homo_sapiens_assembly38.fasta'
    receipt = args.reference_directory / 'reference-receipt.json'
    source = next(r for r in json.loads(receipt.read_text())['sources'] if r['destination'] == reference.name)
    if source['sha256'] != '93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564' or source['bytes'] != reference.stat().st_size:
        raise ValueError('Original exact reference source identity required')
    candidates = []
    with pysam.FastaFile(str(reference)) as fasta:
        for chrom in CHROMS:
            for start, end in intervals[chrom]:
                for a in range(start, end, 500):
                    b = min(a + 500, end)
                    if b - a < 50:
                        counts['segments_shorter_than50bp_excluded'] += 1
                        continue
                    seq = fasta.fetch(chrom, a, b).upper()
                    if set(seq) - set('ACGT'):
                        counts['segments_with_N_or_IUPAC_excluded'] += 1
                        continue
                    gc = (seq.count('G') + seq.count('C')) / len(seq)
                    candidates.append([chrom, a, b, gc])
    sex = [r for r in candidates if r[0] in PAR]
    buckets = defaultdict(list)
    for row in candidates:
        if row[0] not in PAR:
            buckets[((row[2] - row[1] - 1) // 50, min(19, int(row[3] * 20)))].append(row)
    rng = random.Random(20261106)
    for values in buckets.values():
        rng.shuffle(values)
    controls = []
    for row in sex:
        key = ((row[2] - row[1] - 1) // 50, min(19, int(row[3] * 20)))
        if buckets[key]:
            controls.append(buckets[key].pop())
        else:
            counts['sex_windows_without_exact_length50bp_GC5percent_stratum_control'] += 1
    targets = sorted(sex + controls, key=lambda r: (CHROMS.index(r[0]), r[1]))
    if not any(r[0] == 'chrY' for r in targets):
        raise ValueError('No public nonPAR ACGT Y exonic windows')
    target_file = args.output / 'targets.tsv'
    with target_file.open('x', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['chrom', 'start0', 'end0', 'group', 'gc_fraction'])
        for chrom, a, b, gc in targets:
            group = 'X_NONPAR_EXON' if chrom == 'chrX' else 'Y_NONPAR_EXON' if chrom == 'chrY' else 'AUTOSOME_MATCHED'
            writer.writerow([chrom, a, b, group, format(gc, '.12g')])
    info = args.cache / 'info.txt'
    meta = dict(line.split('\t', 1) for line in info.read_text().splitlines())
    if meta['assembly'] != 'GRCh38' or meta['source_1000genomes'] != 'phase3':
        raise ValueError('Pinned GRCh38 phase3 frequency source required')
    fields = meta['variation_cols'].split(',')
    ix = {k: fields.index(k) for k in ['start', 'end', 'failed', 'somatic', 'allele_string', 'AF', 'variation_name']}
    variant_file = args.cache / 'X' / 'all_vars.gz'
    csi = Path(str(variant_file) + '.csi')
    retained = []
    previous = -1000
    with pysam.TabixFile(str(variant_file), index=str(csi)) as variants, pysam.FastaFile(str(reference)) as fasta:
        for a, b in intervals['chrX']:
            seq = fasta.fetch('chrX', a, b).upper()
            for line in variants.fetch('X', a, b):
                row = line.split('\t')
                if len(row) != len(fields):
                    raise ValueError('Public cache row schema changed')
                if row[ix['failed']] != '.' or row[ix['somatic']] != '.':
                    continue
                pos = int(row[ix['start']])
                if not a <= pos - 1 < b or pos - previous < 1000 or row[ix['end']] not in ['.', str(pos)]:
                    continue
                alleles = row[ix['allele_string']].split('/')
                if len(alleles) != 2 or len(set(alleles)) != 2 or any(len(v) != 1 or v not in 'ACGT' for v in alleles):
                    continue
                ref = seq[pos - 1 - a]
                if ref not in alleles:
                    continue
                alt = next(v for v in alleles if v != ref)
                raw = row[ix['AF']]
                if raw.count(':') != 1 or raw.split(':')[0] != alt:
                    continue
                try:
                    af = float(raw.split(':')[1])
                except ValueError:
                    continue
                if not 0.1 <= af <= 0.9:
                    continue
                retained.append(['chrX', pos, ref, alt, raw.split(':')[1], row[ix['variation_name']], format(2 * af * (1 - af), '.12g')])
                previous = pos
                if len(retained) >= 20000:
                    break
            if len(retained) >= 20000:
                break
    site_file = args.output / 'x_nonpar_common_snps.tsv'
    with site_file.open('x', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['chrom', 'pos', 'ref', 'alt', 'af', 'source_id', 'expected_he'])
        writer.writerows(retained)
    if not retained:
        raise ValueError('No exonic public SNPs; do not claim panel')
    manifest = dict(schema_version=1, status='PASS_PUBLIC_TARGET_AWARE_XY_PANEL', assembly='GRCh38',
                    sites=pin(site_file), row_count=len(retained), target_bed=pin(target_file),
                    target_coordinate_system='0-based half-open TSV with header',
                    target_fields=['chrom', 'start0', 'end0', 'group', 'gc_fraction'],
                    target_group_counts=dict(Counter('X_NONPAR_EXON' if r[0] == 'chrX' else 'Y_NONPAR_EXON' if r[0] == 'chrY' else 'AUTOSOME_MATCHED' for r in targets)),
                    target_group_bases={group: sum(r[2] - r[1] for r in targets if ('X_NONPAR_EXON' if r[0] == 'chrX' else 'Y_NONPAR_EXON' if r[0] == 'chrY' else 'AUTOSOME_MATCHED') == group) for group in ['X_NONPAR_EXON', 'Y_NONPAR_EXON', 'AUTOSOME_MATCHED']},
                    exon_scope='Union of all GENCODE49 gene_type protein_coding exons, includes UTR; PAR-overlapping exon rows removed; no patient capture/genotype ascertainment',
                    controls='Public autosomal protein-coding exonic windows sampled without replacement within same50bp length/5percent GC stratum using fixed seed20261106; not bait-matched or proven uniquely mappable',
                    window_width_range=[50, 500], min_spacing_bp=1000, max_sites=20000,
                    allele_frequency_source='1000GENOMES_PHASE3_GLOBAL_AF', frequency_range_inclusive=[0.1, 0.9],
                    expected_he_formula='2*p*(1-p)', named_cache_af_zero_based_column=ix['AF'],
                    population_scope='Global public AF proxy; ancestry and capture ascertainment can bias F; residual LD not modeled',
                    patient_data_used=False, no_patient_selection=True, all_sites_exact_fasta_ref_verified=True,
                    source_gtf=pin(args.gtf), source_gtf_headers=headers, source_gtf_url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz',
                    reference_source=source, reference_receipt=pin(receipt), fai=pin(Path(str(reference) + '.fai')),
                    cache_info=pin(info), cache_source_manifest=pin(args.cache.parents[2] / 'resource-manifest.json'), cache_variant_index=pin(csi),
                    par_source='https://www.ncbi.nlm.nih.gov/grc/human', excluded_par_bed0=PAR,
                    counts=dict(counts), producer=pin(Path(__file__)), elapsed_seconds=round(time.monotonic() - began, 3))
    with (args.output / 'manifest.json').open('x') as handle:
        json.dump(manifest, handle, indent=2)
        handle.write('\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
