#!/usr/bin/env python3
"""Read-only autosomal protein-coding CDS union coverage with pinned mosdepth.

No clinical classification, sex inference, canonical-transcript claim or base-
quality threshold. One complete BAM pass per sample; default mate-overlap handling.
"""
import argparse
import collections
import csv
import datetime as dt
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import subprocess
import sys
import tempfile
import concurrent.futures
import time
import urllib.request

AUTOSOMES = {f'chr{i}' for i in range(1, 23)}
LIMITATIONS = [
    'All-transcript protein-coding CDS union, not a canonical clinical exon panel.',
    'Only chr1-22 assessed; chrX/chrY and mitochondrial/noncoding targets are NA, not zero.',
    'MQ>=20; excluded UNMAP/SECONDARY/QCFAIL/DUP/SUPPLEMENTARY (0xF04).',
    'No BQ filter: low base-quality aligned bases may contribute; coverage is not genotype confidence.',
    'Default mosdepth overlapping-mate correction; no fast-mode or fragment-mode.',
    'Per-gene denominators include reference N/gaps inside CDS; not a callable-genome fraction.',
    'ClinGen HCM includes all exact MONDO:0005045 curations, including disputed/limited evidence; no diagnosis.',
]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*2**20), b''):
            h.update(block)
    return h.hexdigest()


def download_gtf(gtf, source):
    """Two bounded HTTPS ranges; public-resource publisher size/MD5 final check."""
    total=source['bytes'];parts=gtf.parent/(gtf.name+'.parts')
    parts.mkdir(exist_ok=True)
    ranges=[(start,min(start+2**21,total)-1) for start in range(0,total,2**21)]
    def one(bounds):
        start,end=bounds;p=parts/f'{start:012d}-{end:012d}'
        if p.exists() and not p.is_symlink() and p.stat().st_size==end-start+1:return p
        require(not p.is_symlink(),'Unsafe public resource range path')
        tmp=p.with_suffix('.partial')
        for attempt in range(4):
            try:
                req=urllib.request.Request(source['url'],headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'})
                with urllib.request.urlopen(req,timeout=60) as response:
                    require(response.status==206 and response.headers.get('Content-Range')==f'bytes {start}-{end}/{total}','Server did not honor exact bounded HTTP range')
                    with tmp.open('wb') as out:
                        remaining=end-start+1
                        while remaining:
                            block=response.read(min(2**18,remaining))
                            require(block,'Truncated GTF HTTP range');out.write(block);remaining-=len(block)
                        require(not response.read(1),'Oversized GTF HTTP range')
                require(tmp.stat().st_size==end-start+1,'GTF range size mismatch')
                os.replace(tmp,p);return p
            except Exception:
                if attempt==3:raise
                time.sleep(30*(attempt+1))
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        paths=list(pool.map(one,ranges))
    partial=gtf.with_suffix(gtf.suffix+'.partial')
    with partial.open('wb') as out:
        for p in paths:
            with p.open('rb') as stream:
                for block in iter(lambda:stream.read(2**20),b''):out.write(block)
    require(partial.stat().st_size==total and digest(partial,'md5')==source['publisher_md5'],'Publisher GTF identity mismatch')
    os.link(partial,gtf);partial.unlink()


def save(path, value):
    temporary = path.with_name(path.name+f'.{os.getpid()}.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    os.replace(temporary, path)


def union(intervals):
    result = []
    for start, end in sorted(intervals):
        require(0 <= start < end, 'Invalid BED interval')
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def model(gtf, lengths):
    genes, intervals = {}, collections.defaultdict(list)
    with gzip.open(gtf, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            require(len(f) == 9, 'Invalid GTF columns')
            chrom, feature = f[0], f[2]
            if chrom not in lengths or chrom not in AUTOSOMES | {'chrX', 'chrY', 'chrM'}:
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)";', f[8]))
            gid = attrs.get('gene_id')
            require(gid, 'Missing gene_id')
            key = gid+'|'+chrom  # retain chromosome and full version/PAR identity
            if feature == 'gene':
                genes[key] = {'GeneId': gid, 'GeneSymbol': attrs.get('gene_name', ''),
                              'HGNC': attrs.get('hgnc_id', ''), 'Chromosome': chrom,
                              'GeneType': attrs.get('gene_type', '')}
            if feature == 'CDS' and attrs.get('gene_type') == 'protein_coding':
                start, end = int(f[3])-1, int(f[4])
                require(0 <= start < end <= lengths[chrom], 'GTF outside pinned BAM reference')
                intervals[key].append((start, end))
    require(genes and intervals, 'Empty GENCODE model')
    return genes, {key: union(values) for key, values in intervals.items()}


def hcm_curations(path, lock):
    require(digest(path) == lock['csv_sha256'], 'ClinGen snapshot SHA mismatch')
    with Path(path).open(newline='') as stream:
        lines = stream.readlines()
    start = next((i for i, line in enumerate(lines) if '"GENE SYMBOL","GENE ID (HGNC)"' in line), None)
    require(start is not None, 'ClinGen CSV header absent')
    records = [r for r in csv.DictReader(lines[start:]) if r.get('DISEASE ID (MONDO)') == lock['disease_id']]
    require(len(records) == lock['expected_curation_rows'], 'Unexpected exact HCM curation count')
    return records


def write_bed(path, genes, intervals, lengths):
    rows, global_ranges = [], collections.defaultdict(list)
    for key, ranges in intervals.items():
        chrom = genes[key]['Chromosome']
        if chrom not in AUTOSOMES:
            continue
        for start, end in ranges:
            rows.append((chrom, start, end, key))
            global_ranges[chrom].append((start, end))
    for chrom, ranges in global_ranges.items():
        rows.extend((chrom, start, end, '__GLOBAL__') for start, end in union(ranges))
    order = {c: i for i, c in enumerate(lengths)}
    rows.sort(key=lambda r: (order[r[0]], r[1], r[2], r[3]))
    require(rows, 'No autosomal CDS intervals')
    path.write_text(''.join('\t'.join(map(str, row))+'\n' for row in rows))
    return rows


def aggregate(prefix, expected_rows):
    expected = {(c, str(a), str(b), n) for c, a, b, n in expected_rows}
    result = collections.defaultdict(lambda: {'CDSBases': 0, 'DepthSumApprox': 0.0, 'BasesGE1': 0, 'BasesGE10': 0, 'BasesGE20': 0})
    means, seen = {}, set()
    with gzip.open(str(prefix)+'.regions.bed.gz', 'rt') as stream:
        for line in stream:
            f = line.rstrip().split('\t'); key = tuple(f[:4])
            require(key in expected and key not in means, 'Unexpected/duplicate region mean')
            means[key] = float(f[4])
    with gzip.open(str(prefix)+'.thresholds.bed.gz', 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t'); key = tuple(f[:4])
            require(key in expected and key not in seen and len(f) == 7, 'Unexpected/duplicate threshold interval')
            seen.add(key)
            length = int(f[2])-int(f[1]); counts = list(map(int, f[4:]))
            require(0 <= counts[2] <= counts[1] <= counts[0] <= length, 'Invalid threshold counts')
            value = result[f[3]]
            value['CDSBases'] += length
            value['DepthSumApprox'] += means[key]*length
            for name, count in zip(('BasesGE1', 'BasesGE10', 'BasesGE20'), counts):
                value[name] += count
    require(seen == expected == set(means), 'Incomplete mosdepth regions/thresholds')
    for value in result.values():
        total = value['CDSBases']
        value['MeanDepth'] = value.pop('DepthSumApprox')/total
        value['ZeroBases'] = total-value['BasesGE1']
        for field in ('BasesGE1', 'BasesGE10', 'BasesGE20', 'ZeroBases'):
            value['Pct'+field] = 100*value[field]/total
    return dict(result)


def write_tsv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def gene_rows(genes, intervals, stats):
    rows = []
    empty = {name: 'NA' for name in ('CDSBases', 'MeanDepth', 'BasesGE1', 'BasesGE10', 'BasesGE20', 'ZeroBases', 'PctBasesGE1', 'PctBasesGE10', 'PctBasesGE20', 'PctZeroBases')}
    for key, gene in sorted(genes.items()):
        if gene['GeneType'] != 'protein_coding':
            continue
        status = 'ASSESSED' if key in stats else ('NOT_ASSESSED_CHROMOSOME' if gene['Chromosome'] not in AUTOSOMES else 'NO_CDS_MODEL')
        rows.append({**gene, 'Status': status, **(stats[key] if key in stats else empty)})
    return rows


def subset_hcm(curations, rows):
    hgnc = collections.defaultdict(list)
    for row in rows:
        if row['HGNC']:
            hgnc[row['HGNC']].append(row)
    result = []
    for curation in curations:
        matches = hgnc.get(curation['GENE ID (HGNC)'], [])
        if matches:
            for match in matches:
                result.append({**curation, **match, 'MatchMethod': 'exact_HGNC_id'})
        else:
            result.append({**curation, 'Status': 'NA_NO_PROTEIN_CODING_MODEL_MATCH', 'MatchMethod': 'unmatched_no_symbol_guess', 'CDSBases': 'NA', 'MeanDepth': 'NA'})
    return result


def tool_fixture(mosdepth, output):
    output.mkdir(parents=True, exist_ok=False)
    sam = output/'fixture.sam'; bam = output/'fixture.bam'; bed = output/'fixture.bed'
    seq, qual = 'A'*100, '!'*100  # explicitly BQ0; counted because mosdepth has no BQ filter
    rows = [f'pair\t99\tchr1\t1\t60\t100M\t=\t51\t150\t{seq}\t{qual}',
            f'pair\t147\tchr1\t51\t60\t100M\t=\t1\t-150\t{seq}\t{qual}']
    for index, (flag, mq) in enumerate(((1024,60),(2048,60),(256,60),(512,60),(0,19))):
        rows.append(f'excluded{index}\t{flag}\tchr1\t1\t{mq}\t100M\t*\t0\t0\t{seq}\t{qual}')
    sam.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:200\n'+'\n'.join(rows)+'\n')
    subprocess.run(['samtools','sort','-o',str(bam),str(sam)],check=True)
    subprocess.run(['samtools','index',str(bam)],check=True)
    bed.write_text('chr1\t0\t100\tgeneA\nchr1\t0\t200\tfixture\nchr1\t50\t150\tgeneB\n')
    prefix = output/'fixture'
    subprocess.run([str(mosdepth),'-n','-t','1','-Q','20','-F','3844','--by',str(bed),'--thresholds','1,10,20',str(prefix),str(bam)],env={**os.environ,'MOSDEPTH_PRECISION':'8'},check=True)
    all_metrics = aggregate(prefix,[('chr1',0,100,'geneA'),('chr1',0,200,'fixture'),('chr1',50,150,'geneB')])
    metrics = all_metrics['fixture']
    require(abs(metrics['MeanDepth']-.75)<1e-8 and metrics['BasesGE1']==150 and metrics['BasesGE10']==0, 'Actual mosdepth mate-overlap/filter/BQ0 fixture failed')
    require(all_metrics['geneA']['MeanDepth']==1 and all_metrics['geneB']['MeanDepth']==1, 'Overlapping BED regions incorrectly aggregated')
    save(output/'validation-summary.json',{'status':'PASS','metrics':metrics,'proves':['default_mate_overlap_removed','BQ0_included_no_BQ_filter','MQ19_DUP_SUPP_SECONDARY_QCFAIL_excluded']})
    return output/'validation-summary.json'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',required=True)
    ap.add_argument('--execute',action='store_true')
    args=ap.parse_args()
    cfg=json.loads(Path(args.config).read_text())
    if not args.execute:
        print('DRY RUN: no resource download, mosdepth run or result files written');return
    os.umask(0o077)
    cpus=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,cpus[:4]);os.nice(10)
    resource.setrlimit(resource.RLIMIT_AS,(16*2**30,16*2**30))
    output=Path(cfg['output_root']);output.mkdir(parents=True,exist_ok=True)
    lockfile=(output/'.run.lock').open('a');fcntl.flock(lockfile,fcntl.LOCK_EX|fcntl.LOCK_NB)
    status=Path(cfg['status_file'])
    sources=json.loads(Path(cfg['sources_lock']).read_text()); mosdepth=Path(cfg['mosdepth'])
    require(digest(mosdepth)==sources['mosdepth']['sha256_observed'],'mosdepth executable mismatch')
    gtf=Path(cfg['gtf']);gtf.parent.mkdir(parents=True,exist_ok=True)
    save(status,{'status':'RUNNING','step':'prepare_resources','updated_utc':dt.datetime.now(dt.timezone.utc).isoformat()})
    if not gtf.exists():
        download_gtf(gtf,sources['gencode_gtf'])
    require(gtf.stat().st_size==sources['gencode_gtf']['bytes'] and digest(gtf,'md5')==sources['gencode_gtf']['publisher_md5'],'GTF identity mismatch')
    lengths={f[0]:int(f[1]) for line in Path(cfg['reference']+'.fai').read_text().splitlines() if (f:=line.split('\t'))}
    genes,intervals=model(gtf,lengths);bed=output/'protein-coding-cds-union.bed';bed_rows=write_bed(bed,genes,intervals,lengths)
    curations=hcm_curations(cfg['clingen_csv'],sources['clingen_hcm'])
    fixture=tool_fixture(mosdepth,output/('tool-fixture-'+dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S')))
    source={'gtf_sha256':digest(gtf),'gtf_publisher_md5':sources['gencode_gtf']['publisher_md5'],'sources_lock_sha256':digest(cfg['sources_lock']),'script_sha256':digest(__file__),'bed_sha256':digest(bed),'fixture_sha256':digest(fixture),'reference_fasta_sha256_from_manifest':cfg['reference_sha256'],'reference_fai_sha256':digest(cfg['reference']+'.fai'),'limitations':LIMITATIONS}
    for index,sample in enumerate(cfg['samples']):
        sid=sample['sample_id'];require(re.fullmatch(r'[A-Za-z0-9_.-]+',sid),'Unsafe sample ID')
        dest=output/sid
        receipt=dest/'workflow-receipt.json'
        if receipt.exists():
            existing=json.loads(receipt.read_text());require(existing.get('status')=='PASS' and existing.get('sources')==source,'Existing receipt belongs to another model/run');continue
        dest.mkdir(exist_ok=True)
        bam=Path(sample['bam']);rp=Path(sample['receipt']);r=json.loads(rp.read_text())
        require(r.get('completed_utc') and r.get('fingerprint') and r.get('sample')==sid,'Source alignment receipt not completed or wrong sample')
        bamsha=r.get('outputs',{}).get(bam.name);require(isinstance(bamsha,str) and re.fullmatch(r'[0-9a-f]{64}',bamsha),'No source BAM receipt identity')
        subprocess.run(['samtools','quickcheck','-v',str(bam)],check=True)
        before=bam.stat();prefix=dest/'coverage'
        header=subprocess.check_output(['samtools','view','-H',str(bam)],text=True)
        bamlengths={d['SN']:int(d['LN']) for line in header.splitlines() if line.startswith('@SQ\t') for d in [dict(p.split(':',1) for p in line.split('\t')[1:])]}
        require(all(bamlengths.get(c)==lengths[c] for c in AUTOSOMES),'BAM reference lengths mismatch')
        save(status,{'status':'RUNNING','step':'mosdepth','sample_index':index+1,'sample_count':len(cfg['samples']),'updated_utc':dt.datetime.now(dt.timezone.utc).isoformat()})
        command=[str(mosdepth),'-n','-t','3','-Q','20','-F','3844','--by',str(bed),'--thresholds','1,10,20',str(prefix),str(bam)]
        with (dest/'mosdepth.log').open('w') as log:
            subprocess.run(command,env={**os.environ,'MOSDEPTH_PRECISION':'8'},stdout=log,stderr=subprocess.STDOUT,check=True)
        require((bam.stat().st_size,bam.stat().st_mtime_ns)==(before.st_size,before.st_mtime_ns),'Source BAM changed during read-only QC')
        stats=aggregate(prefix,bed_rows);rows=gene_rows(genes,intervals,stats)
        write_tsv(dest/'gene-cds-coverage.tsv',rows)
        write_tsv(dest/'clingen-hcm-cds-coverage.tsv',subset_hcm(curations,rows))
        save(dest/'coverage-summary.json',{'status':'PASS','global_unique_autosomal_cds':stats['__GLOBAL__'],'assessed_genes':sum(x['Status']=='ASSESSED' for x in rows),'clingen_hcm_curation_rows':len(curations),'limitations':LIMITATIONS})
        artifacts={p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in dest.iterdir() if p.is_file() and p.name!='workflow-receipt.json'}
        save(receipt,{'schema_version':1,'status':'PASS','identity':hashlib.sha256(json.dumps([source,bamsha,sid],sort_keys=True).encode()).hexdigest(),'completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'sources':source,'source_bam':{'path':str(bam),'size_bytes':before.st_size,'sha256_from_completed_receipt':bamsha,'receipt_sha256':digest(rp),'rehashed':False},'command':command,'artifacts':artifacts})
        print(f'PASS gene coverage {index+1}/{len(cfg["samples"])}',flush=True)
    save(status,{'status':'PASS','completed_samples':len(cfg['samples']),'updated_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'sources':source})


if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        # Retain partial outputs for diagnosis, never emit a PASS receipt on failure.
        if '--config' in sys.argv:
            try:
                config=json.loads(Path(sys.argv[sys.argv.index('--config')+1]).read_text())
                save(Path(config['status_file']),{'status':'FAIL','error':type(exc).__name__+': '+str(exc),'updated_utc':dt.datetime.now(dt.timezone.utc).isoformat()})
            except Exception:
                pass
        raise
