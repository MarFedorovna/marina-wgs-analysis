"""Independent all-site validation and explicit public-panel contract; no BAM input."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import pysam


def pin(path):
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.source_manifest.read_text())
    if (source['status'] != 'PASS_PUBLIC_X_NONPAR_COMMON_SNP_PANEL'
            or source['patient_data_used'] is not False or source['no_patient_selection'] is not True
            or source['frequency_source'] != '1000 Genomes phase3 global ALT AF as represented in VEP115 cache'):
        raise ValueError('Exact public source contract required')
    tsv = Path(source['output']['path'])
    if pin(tsv) != source['output']:
        raise ValueError('Public site artifact changed')
    metadata = dict(line.split('\t', 1) for line in Path(source['cache_info']['path']).read_text().splitlines())
    if (pin(Path(source['cache_info']['path'])) != source['cache_info']
            or metadata['assembly'] != 'GRCh38' or metadata['source_1000genomes'] != 'phase3'):
        raise ValueError('Pinned public frequency metadata changed')
    ref = Path(source['reference_receipt']['path']).parent / source['reference_source']['destination']
    previous = None
    count = 0
    minimum_spacing = None
    with pysam.FastaFile(str(ref)) as fasta, tsv.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if reader.fieldnames != ['chrom', 'pos', 'ref', 'alt', 'af', 'source_id', 'expected_he']:
            raise ValueError('Unexpected public TSV fields')
        for row in reader:
            pos = int(row['pos'])
            p = float(row['af'])
            if (row['chrom'] != 'chrX' or not 2781479 <= pos - 1 < 155701382
                    or row['ref'] not in 'ACGT' or len(row['ref']) != 1
                    or row['alt'] not in 'ACGT' or len(row['alt']) != 1
                    or row['ref'] == row['alt'] or not 0.1 <= p <= 0.9
                    or abs(float(row['expected_he']) - 2 * p * (1 - p)) > 1e-10
                    or fasta.fetch('chrX', pos - 1, pos).upper() != row['ref']):
                raise ValueError('Public SNP allele/AF/PAR/FASTA invariant failed')
            if previous is not None:
                spacing = pos - previous
                if spacing < 10000:
                    raise ValueError('Public panel spacing/order invariant failed')
                minimum_spacing = spacing if minimum_spacing is None else min(minimum_spacing, spacing)
            previous = pos
            count += 1
    if count != source['row_count'] or not 100 <= count <= 5000:
        raise ValueError('Public panel denominator invariant failed')
    contract = dict(schema_version=1, status='PASS_PUBLIC_X_NONPAR_COMMON_SNP_PANEL',
                    assembly='GRCh38', allele_frequency_source='1000GENOMES_PHASE3_GLOBAL_AF',
                    sites=source['output'], site_count=count, row_count=count,
                    patient_data_used=False, no_patient_selection=True,
                    all_sites_exact_fasta_ref_verified=True, all_sites_nonpar_verified=True,
                    all_sites_biallelic_common_af_verified=True,
                    actual_min_spacing_bp=minimum_spacing,
                    source_manifest=pin(args.source_manifest), producer=pin(Path(__file__)),
                    population_scope=source['population_scope'], linkage_scope=source['linkage_scope'],
                    capture_scope=source['capture_scope'], cache_info=source['cache_info'],
                    cache_source_manifest=source['cache_source_manifest'],
                    cache_variant_index=source['cache_variant_index'],
                    reference_receipt=source['reference_receipt'], reference_source=source['reference_source'],
                    expected_he_formula='2*p*(1-p)', frequency_range_inclusive=[0.1, 0.9],
                    named_cache_af_zero_based_column=metadata['variation_cols'].split(',').index('AF'))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(contract, handle, indent=2)
        handle.write('\n')
    print(json.dumps(contract, indent=2))


if __name__ == '__main__':
    main()
