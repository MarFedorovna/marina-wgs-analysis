"""Read-only small public-reference mask/metadata proof; no BAM or patient input."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


PAR = {
    'chrX': [(10000, 2781479), (155701382, 156030895)],
    'chrY': [(10000, 2781479), (56887902, 57217415)],
}
REFERENCE_SHA = '93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564'


def pin(path):
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def fetch(handle, record, start, end):
    length, offset, line_bases, line_bytes = record
    if not 0 <= start < end <= length or end - start > 3_000_000:
        raise ValueError('Public-reference interval out of bounds')
    left = offset + start // line_bases * line_bytes + start % line_bases
    right = offset + (end - 1) // line_bases * line_bytes + (end - 1) % line_bases + 1
    handle.seek(left)
    sequence = handle.read(right - left).replace(b'\n', b'').replace(b'\r', b'').upper()
    if len(sequence) != end - start:
        raise ValueError('FAI sequence read length mismatch')
    return sequence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-directory', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    fasta = args.reference_directory / 'Homo_sapiens_assembly38.fasta'
    fai = Path(str(fasta) + '.fai')
    receipt = args.reference_directory / 'reference-receipt.json'
    sources = json.loads(receipt.read_text())['sources']
    source = next(r for r in sources if r['destination'] == fasta.name)
    if source['sha256'] != REFERENCE_SHA or source['bytes'] != fasta.stat().st_size:
        raise ValueError('Pinned original reference source identity required')
    index = {r[0]: tuple(map(int, r[1:])) for line in fai.read_text().splitlines()
             if (r := line.split('\t'))}
    if index['chrX'][0] != 156040895 or index['chrY'][0] != 57227415:
        raise ValueError('Unexpected GRCh38 chromosome lengths')
    evidence = []
    with fasta.open('rb') as handle:
        for chrom, regions in PAR.items():
            for i, (start, end) in enumerate(regions, 1):
                sequence = fetch(handle, index[chrom], start, end)
                evidence.append(dict(chrom=chrom, label='PAR' + str(i), start0=start, end0=end,
                                     bases=len(sequence), base_counts=dict(Counter(map(chr, sequence))),
                                     all_N=all(base == 78 for base in sequence),
                                     sequence_sha256=hashlib.sha256(sequence).hexdigest()))
    info = args.cache / 'info.txt'
    metadata = dict(line.split('\t', 1) for line in info.read_text().splitlines())
    report = dict(status='PASS_PUBLIC_XY_REFERENCE_METADATA_PROBE', patient_data_used=False,
                  source_url='https://www.ncbi.nlm.nih.gov/grc/human',
                  masking_documentation='https://gatk.broadinstitute.org/hc/en-us/articles/360041155232-Reference-Genome-Components',
                  coordinate_system='0-based half-open; converted from GRC 1-based closed PAR',
                  original_reference_sha256_from_verified_receipt=REFERENCE_SHA,
                  reference_whole_file_rehashed=False, reference_receipt=pin(receipt), fai=pin(fai),
                  par_reference_probes=evidence, cache_info=pin(info),
                  cache_assembly=metadata.get('assembly'), variation_columns=metadata['variation_cols'].split(','),
                  frequency_sources={k: v for k, v in metadata.items() if k in
                                     ['source_1000genomes', 'source_gnomADe', 'source_gnomADg', 'source_dbSNP']},
                  local_x_variant_cache_bytes=(args.cache / 'X' / 'all_vars.gz').stat().st_size,
                  local_x_variant_cache_index=pin(args.cache / 'X' / 'all_vars.gz.csi'),
                  interpretation='Reference-only proof; masked or uncallable bases cannot establish biological absence; no sex/karyotype inference')
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2)
        handle.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
