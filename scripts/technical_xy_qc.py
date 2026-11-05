#!/usr/bin/env python3
"""Read-only filtered chromosome-density QC; never assigns sex or changes ploidy."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

PAR_SOURCE = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/pseudoautosomal_region/par_align.gff'
PAR_SHA256 = 'cc1b9b6dac63774e8e22ec87cfb71edfea59483534520094a1758a4a509f7311'
PAR = {'chrX': [(10001, 2781479), (155701383, 156030895)],
       'chrY': [(10001, 2781479), (56887903, 57217415)]}


def require(ok, text):
    if not ok:
        raise ValueError(text)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def complement(length, ranges):
    result, start = [], 1
    for left, right in ranges:
        if start < left:
            result.append((start, left-1))
        start = right+1
    if start <= length:
        result.append((start, length))
    return result


def reference_counts(reference, chrom, ranges):
    counts = {'ACGT': 0, 'N': 0, 'Other': 0}
    for start, stop in ranges:
        process = subprocess.Popen(['samtools', 'faidx', str(reference), f'{chrom}:{start}-{stop}'], stdout=subprocess.PIPE)
        process.stdout.readline()  # Exclude the FASTA header from nucleotide counts.
        while True:
            data = process.stdout.read(8*1024*1024)
            if not data:
                break
            data = data.upper().replace(b'\n', b'').replace(b'\r', b'')
            acgt = sum(data.count(base) for base in (b'A', b'C', b'G', b'T'))
            n = data.count(b'N')
            counts['ACGT'] += acgt
            counts['N'] += n
            counts['Other'] += len(data)-acgt-n
        require(process.wait() == 0, 'Reference faidx failed')
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--output-root', required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    ref, par_file = Path(cfg['reference']), Path(cfg['par_source'])
    require(sha(par_file) == PAR_SHA256, 'Pinned primary PAR artifact changed')
    fai = Path(str(ref) + '.fai')
    lengths = {p[0]: int(p[1]) for p in (line.split('\t') for line in fai.read_text().splitlines())}
    require(lengths.get('chrX') == 156040895 and lengths.get('chrY') == 57227415, 'GRCh38 primary X/Y lengths required')
    ranges = {'chr'+str(n): [(1, lengths['chr'+str(n)])] for n in range(1, 23)}
    ranges.update({chrom: complement(lengths[chrom], PAR[chrom]) for chrom in PAR})
    for sample in cfg['samples']:
        require(re.fullmatch('[A-Za-z0-9_.-]+', sample['sample_id']), 'Invalid sample code')
        bam, receipt = Path(sample['bam']), Path(sample['receipt'])
        require(bam.is_file() and Path(str(bam)+'.bai').is_file() and receipt.parent == bam.parent, 'Completed BAM/index/receipt required')
        saved = json.loads(receipt.read_text())
        require(saved.get('sample') == sample['sample_id'] and saved.get('completed_utc') and bam.name in saved.get('outputs', {}), 'BAM receipt mismatch')
        sample['expected_bam_sha256_from_receipt'] = saved['outputs'][bam.name]
    if not args.execute:
        print(f'DRY RUN: {len(cfg["samples"])} completed BAMs; MQ20/-F0xF04; primary-source nonPAR; no sex/ploidy change')
        return
    os.umask(0o077)
    os.nice(15)
    cpus = sorted(os.sched_getaffinity(0))[-2:]
    os.sched_setaffinity(0, cpus)
    extra_threads = 1 if len(cpus) > 1 else 0
    out = Path(args.output_root)
    out.mkdir(parents=True, exist_ok=False, mode=0o700)
    started = dt.datetime.now(dt.timezone.utc).isoformat()
    version = subprocess.check_output(['samtools', '--version'], text=True, errors='replace').splitlines()[0]
    report = {'SchemaVersion': 1, 'StartedUtc': started, 'Status': 'RUNNING', 'SexAssigned': False, 'PloidyChanged': False,
              'PARSourceURL': PAR_SOURCE, 'PARSourceSha256': PAR_SHA256, 'PAROneBasedClosed': PAR,
              'ReferenceFaiSha256': sha(fai), 'ReferenceExpectedSha256': cfg.get('reference_sha256'),
              'SamtoolsVersion': version, 'CodeSha256': sha(__file__), 'CPUCountLimit': len(cpus),
              'Method': {'MinimumMappingQuality': 20, 'ExcludedFlags': '0xF04', 'ExtraSamtoolsThreads': extra_threads,
                         'Count': 'primary nonduplicate QC-passing mapped read-segments; not fragments or depth',
                         'BaseQualityFilter': None, 'PairedOverlapRemoved': False, 'RegionsAreOneBasedClosed': True},
              'ReferenceDenominators': {}, 'Samples': [], 'Limitations': [
                  'Strongly nonuniform/target-biased coverage prevents confident XX/XY assignment from ratios alone',
                  'NonPAR is not a unique-mappability mask; MQ filtering reduces but does not eliminate homologous/repetitive mapping',
                  'Read-segment density is not depth, breadth, callable fraction or a karyotype',
                  'No sex label, phenotype or calling ploidy was inferred or changed']}
    for chrom in ranges:
        report['ReferenceDenominators'][chrom] = {'FullSpan': sum(b-a+1 for a,b in ranges[chrom]),
                                                 **reference_counts(ref, chrom, ranges[chrom])}
    for sample in cfg['samples']:
        row = {'SampleId': sample['sample_id'], 'BAM': sample['bam'],
               'ExpectedBamSha256FromCompletedReceipt': sample['expected_bam_sha256_from_receipt'], 'Counts': {}}
        report['Samples'].append(row)
        for chrom in ['chrX', 'chrY'] + ['chr'+str(n) for n in range(1,23)]:
            regions = [f'{chrom}:{a}-{b}' for a,b in ranges[chrom]]
            command = ['samtools', 'view', '-@', str(extra_threads), '-c', '-M', '-q', '20', '-F', '0xF04', sample['bam'], *regions]
            began = time.monotonic()
            count = int(subprocess.check_output(command, text=True).strip())
            row['Counts'][chrom] = {'MappedReadSegments': count, 'Seconds': round(time.monotonic()-began, 3)}
            save(out/'progress.json', report)
            print(f'Count complete {sample["sample_id"]} {chrom}', flush=True)
        auto_count = sum(row['Counts']['chr'+str(n)]['MappedReadSegments'] for n in range(1,23))
        auto_span = sum(report['ReferenceDenominators']['chr'+str(n)]['FullSpan'] for n in range(1,23))
        auto_acgt = sum(report['ReferenceDenominators']['chr'+str(n)]['ACGT'] for n in range(1,23))
        row['Autosomes'] = {'MappedReadSegments': auto_count, 'FullSpan': auto_span, 'ACGTSpan': auto_acgt}
        densities = [row['Counts']['chr'+str(n)]['MappedReadSegments']/report['ReferenceDenominators']['chr'+str(n)]['FullSpan'] for n in range(1,23)]
        row['AutosomalDensityMaxMin'] = max(densities)/min(densities) if min(densities)>0 else None
        row['NonPARComparisons'] = {}
        for chrom in ('chrX','chrY'):
            count = row['Counts'][chrom]['MappedReadSegments'];denom=report['ReferenceDenominators'][chrom]
            row['NonPARComparisons'][chrom] = {
                'MappedReadSegments': count, 'FullSpan': denom['FullSpan'], 'ACGTSpan': denom['ACGT'],
                'DensityRatioToAutosomesFullSpan': (count/denom['FullSpan'])/(auto_count/auto_span) if auto_count else None,
                'DensityRatioToAutosomesACGTSpan': (count/denom['ACGT'])/(auto_count/auto_acgt) if auto_count and denom['ACGT'] else None}
    report.update(Status='PASS', CompletedUtc=dt.datetime.now(dt.timezone.utc).isoformat())
    save(out/'technical-xy-profile.json', report)
    with (out/'technical-xy-profile.tsv').open('w') as stream:
        stream.write('sample_id\tautosomal_mapped_segments\tX_nonPAR_mapped\tY_nonPAR_mapped\tX_A_density_ratio_fullspan\tY_A_density_ratio_fullspan\tX_A_density_ratio_ACGT\tY_A_density_ratio_ACGT\tautosomal_density_max_min\tsex_assigned\n')
        for row in report['Samples']:
            x,y=row['NonPARComparisons']['chrX'],row['NonPARComparisons']['chrY']
            stream.write('\t'.join(map(str,[row['SampleId'],row['Autosomes']['MappedReadSegments'],x['MappedReadSegments'],y['MappedReadSegments'],x['DensityRatioToAutosomesFullSpan'],y['DensityRatioToAutosomesFullSpan'],x['DensityRatioToAutosomesACGTSpan'],y['DensityRatioToAutosomesACGTSpan'],row['AutosomalDensityMaxMin'],False]))+'\n')
    artifacts = {name: {'bytes': (out/name).stat().st_size, 'sha256': sha(out/name)} for name in ('technical-xy-profile.json','technical-xy-profile.tsv')}
    identity = hashlib.sha256(json.dumps({'config':cfg,'code':report['CodeSha256']},sort_keys=True).encode()).hexdigest()
    save(out/'workflow-receipt.json', {'schema_version':1,'status':'PASS','stage':'technical_xy_qc_only','identity':identity,
         'started_utc':started,'completed_utc':report['CompletedUtc'],'artifacts':artifacts,'metadata':{'SexAssigned':False,'PloidyChanged':False}})
    print('PASS technical XY QC; no sex/ploidy assignment', flush=True)


if __name__ == '__main__':
    main()
