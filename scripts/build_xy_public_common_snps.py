"""Public-only sparse X-nonPAR panel from installed VEP cache; no BAM input."""
import argparse
from collections import Counter
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
import time

import pysam

PAR = [(10000, 2781479), (155701382, 156030895)]


def pin(path):
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--reference-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    hard_end = dt.datetime.fromisoformat('2026-11-06T04:10:00+00:00')

    def guard():
        if time.monotonic() - started > 180 or dt.datetime.now(dt.timezone.utc) >= hard_end:
            raise TimeoutError('Public derivative bounded180s/original approved XY task deadline')

    args.output.mkdir(mode=0o700, parents=True, exist_ok=False)
    info = args.cache / 'info.txt'
    metadata = dict(line.split('\t', 1) for line in info.read_text().splitlines())
    if metadata['assembly'] != 'GRCh38' or metadata['source_1000genomes'] != 'phase3':
        raise ValueError('Exact named public frequency source required')
    fields = metadata['variation_cols'].split(',')
    ix = {name: fields.index(name) for name in ['chr', 'variation_name', 'failed', 'somatic',
                                               'start', 'end', 'allele_string', 'AF']}
    reference = args.reference_directory / 'Homo_sapiens_assembly38.fasta'
    receipt = args.reference_directory / 'reference-receipt.json'
    ref_source = next(r for r in json.loads(receipt.read_text())['sources']
                      if r['destination'] == reference.name)
    if (ref_source['sha256'] != '93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564'
            or ref_source['bytes'] != reference.stat().st_size):
        raise ValueError('Original pinned reference source identity changed')
    variant_file = args.cache / 'X' / 'all_vars.gz'
    csi = Path(str(variant_file) + '.csi')
    counts = Counter()
    retained = []
    previous = -10000
    # Fixed binning is independent of all patient genotypes, coverage and phenotype.
    # Central10kb of31kb bins yields >20kb spacing between different bins.
    windows = [(start + 10000, min(start + 20000, PAR[1][0]))
               for start in range(PAR[0][1], PAR[1][0] - 20000, 31000)]
    with pysam.TabixFile(str(variant_file), index=str(csi)) as variants, \
            pysam.FastaFile(str(reference)) as fasta:
        if fasta.get_reference_length('chrX') != 156040895:
            raise ValueError('GRCh38 chrX length mismatch')
        for left, right in windows:
            guard()
            sequence = fasta.fetch('chrX', left, right).upper()
            counts['windows_examined'] += 1
            for line in variants.fetch('X', left, right):
                counts['public_rows_examined'] += 1
                row = line.split('\t')
                if len(row) != len(fields):
                    raise ValueError('Cache row does not match declared variation_cols')
                if row[ix['failed']] != '.' or row[ix['somatic']] != '.':
                    continue
                alleles = row[ix['allele_string']].split('/')
                if len(alleles) != 2 or len(set(alleles)) != 2 or any(a not in 'ACGT' or len(a) != 1 for a in alleles):
                    continue
                pos = int(row[ix['start']])
                if not left <= pos - 1 < right or pos - previous < 10000:
                    continue
                if row[ix['end']] not in ['.', str(pos)]:
                    continue
                ref = sequence[pos - 1 - left]
                if ref not in alleles:
                    continue
                alt = next(a for a in alleles if a != ref)
                # VEP stores AF as allele:value. Never confuse MAF with ALT AF.
                value = row[ix['AF']]
                if value == '.' or value.count(':') != 1:
                    continue
                allele, raw_frequency = value.split(':')
                if allele != alt:
                    continue
                try:
                    af = float(raw_frequency)
                except ValueError:
                    continue
                if not 0.1 <= af <= 0.9:
                    continue
                retained.append(['chrX', str(pos), ref, alt, raw_frequency,
                                 row[ix['variation_name']], format(2 * af * (1 - af), '.12g')])
                previous = pos
                counts['windows_with_retained_snp'] += 1
                break
    guard()
    target = args.output / 'x_nonpar_common_snps.tsv'
    with target.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['chrom', 'pos', 'ref', 'alt', 'af', 'source_id', 'expected_he'])
        writer.writerows(retained)
    if not retained:
        raise ValueError('No public common exact-reference SNPs; do not claim a panel')
    manifest = dict(schema_version=1, status='PASS_PUBLIC_X_NONPAR_COMMON_SNP_PANEL',
                    patient_data_used=False, no_patient_selection=True,
                    output=pin(target), row_count=len(retained), counts=dict(counts),
                    frequency_source='1000 Genomes phase3 global ALT AF as represented in VEP115 cache',
                    population_scope='Global public population proxy; ancestry mismatch may bias expected heterozygosity',
                    frequency_range_inclusive=[0.1, 0.9], expected_he_formula='2*p*(1-p)',
                    genomic_scope='chrX nonPAR; fixed central10kb per31kb bin, first eligible SNP, at mostone/window',
                    min_spacing_bp=10000, linkage_scope='Sparse spacing reduces local linkage; LD not measured; sites are not guaranteed independent',
                    capture_scope='No capture/target/BAM selection; sparse off-target sites may be uncallable; callable denominator must be reported per sample',
                    reference_scope='Exact FASTA REF; original public reference SHA from prior source receipt, no whole reference rehash',
                    reference_source=ref_source, reference_receipt=pin(receipt),
                    fai=pin(Path(str(reference) + '.fai')), cache_info=pin(info),
                    cache_source_manifest=pin(args.cache.parents[2] / 'resource-manifest.json'),
                    cache_variant=dict(path=str(variant_file), bytes=variant_file.stat().st_size,
                                       whole_file_rehashed=False), cache_variant_index=pin(csi),
                    par_source='https://www.ncbi.nlm.nih.gov/grc/human',
                    excluded_par_bed0=PAR, producer=pin(Path(__file__)),
                    elapsed_seconds=round(time.monotonic() - started, 3),
                    created_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    with (args.output / 'manifest.json').open('x', encoding='utf-8') as handle:
        json.dump(manifest, handle, indent=2)
        handle.write('\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()
