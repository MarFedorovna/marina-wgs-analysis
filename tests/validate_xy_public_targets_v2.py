"""Independent public target/SNP/FASTA checks and exact downstream contract; no BAM."""
import argparse
import bisect
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import pysam


def pin(path):
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.manifest.read_text())
    assert source['status'] == 'PASS_PUBLIC_TARGET_AWARE_XY_PANEL'
    assert source['patient_data_used'] is False and source['no_patient_selection'] is True
    assert source['assembly'] == 'GRCh38'
    assert source['allele_frequency_source'] == '1000GENOMES_PHASE3_GLOBAL_AF'
    for key in ['source_gtf', 'sites', 'target_bed', 'reference_receipt', 'fai', 'cache_info', 'cache_variant_index']:
        assert pin(Path(source[key]['path'])) == source[key], key
    exons = defaultdict(list)
    with gzip.open(source['source_gtf']['path'], 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            r = line.split('\t')
            if r[2] == 'exon' and 'gene_type "protein_coding";' in r[8]:
                exons[r[0]].append((int(r[3]) - 1, int(r[4])))
    union = {}
    for chrom, regions in exons.items():
        merged = []
        for start, end in sorted(regions):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        union[chrom] = merged
    starts = {c: [a for a, b in regions] for c, regions in union.items()}

    def exonic(chrom, a, b):
        i = bisect.bisect_right(starts.get(chrom, []), a) - 1
        return i >= 0 and union[chrom][i][0] <= a < b <= union[chrom][i][1]

    def nonpar(chrom, a, b):
        return not any(a < end and b > start for start, end in source['excluded_par_bed0'].get(chrom, []))

    ref = Path(source['reference_receipt']['path']).parent / source['reference_source']['destination']
    previous = {}
    target_counts = Counter()
    strata = defaultdict(Counter)
    with pysam.FastaFile(str(ref)) as fasta:
        with Path(source['target_bed']['path']).open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            assert reader.fieldnames == ['chrom', 'start0', 'end0', 'group', 'gc_fraction']
            for r in reader:
                c, a, b, group = r['chrom'], int(r['start0']), int(r['end0']), r['group']
                assert group == ('X_NONPAR_EXON' if c == 'chrX' else 'Y_NONPAR_EXON' if c == 'chrY' else 'AUTOSOME_MATCHED')
                assert 50 <= b - a <= 500 and a >= previous.get(c, -1)
                assert exonic(c, a, b) and nonpar(c, a, b)
                sequence = fasta.fetch(c, a, b).upper()
                assert len(sequence) == b - a and not set(sequence) - set('ACGT')
                gc = (sequence.count('C') + sequence.count('G')) / len(sequence)
                assert abs(gc - float(r['gc_fraction'])) < 1e-9
                previous[c] = b
                target_counts[group] += 1
                strata['auto' if group == 'AUTOSOME_MATCHED' else 'sex'][((b-a-1)//50, min(19, int(gc*20)))] += 1
        assert all(n <= strata['sex'][key] for key, n in strata['auto'].items())
        count = 0
        prior = None
        with Path(source['sites']['path']).open() as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            assert reader.fieldnames == ['chrom', 'pos', 'ref', 'alt', 'af', 'source_id', 'expected_he']
            for r in reader:
                c, pos, p = r['chrom'], int(r['pos']), float(r['af'])
                assert c == 'chrX' and nonpar(c, pos-1, pos) and exonic(c, pos-1, pos)
                assert all(len(r[k]) == 1 and r[k] in 'ACGT' for k in ['ref', 'alt']) and r['ref'] != r['alt']
                assert .1 <= p <= .9 and abs(float(r['expected_he']) - 2*p*(1-p)) < 1e-10
                assert fasta.fetch(c, pos-1, pos).upper() == r['ref']
                assert prior is None or pos - prior >= 1000
                prior = pos
                count += 1
    assert count == source['row_count'] and 100 <= count <= 20000
    assert dict(target_counts) == source['target_group_counts']
    contract = dict(source)
    contract.update(targets=source['target_bed'], source_manifest=pin(args.manifest),
                    independent_validation_producer=pin(Path(__file__)),
                    all_targets_disjoint_exonic_nonpar_ACGT_GC_verified=True,
                    all_sites_exonic_nonpar_refalt_AF_spacing_verified=True,
                    autosome_control_length_GC_strata_verified=True)
    with args.output.open('x') as handle:
        json.dump(contract, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status=contract['status'], sites=count, targets=dict(target_counts),
                          contract=pin(args.output)), indent=2))


if __name__ == '__main__':
    main()
