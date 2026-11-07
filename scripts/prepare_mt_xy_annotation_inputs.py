#!/usr/bin/env python3
"""Prepare existing local MT/XY calls and annotation sidecars for further annotation.

Streams source files; never recalls variants, changes GT/ploidy, or uploads data.
Nuclear QC uses existing site PASS + alternate genotype FT PASS, DP>=10,GQ>=20.
MT uses literal site PASS and retains Mutect2 AF/AD/DP without nuclear GQ gates.
"""
import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


def read_vcf(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else path.open(encoding='utf-8')


def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return -1


def alt_gt(gt):
    a = re.split(r'[/|]', gt)
    return all(x.isdigit() for x in a) and any(int(x) > 0 for x in a)


def prepare_vcf(source, target, branch):
    selected, counts = set(), Counter()
    target.parent.mkdir(parents=True, exist_ok=True)
    with read_vcf(source) as src, target.open('w', encoding='utf-8', newline='') as dst:
        for line in src:
            if line.startswith('#'):
                dst.write(line)
                if line.startswith('#CHROM'):
                    samples = line.rstrip().split('\t')[9:]
                continue
            cells = line.rstrip('\n').split('\t')
            counts['input_records'] += 1
            if cells[6] != 'PASS':
                continue
            fmt = cells[8].split(':')
            good = []
            for sample, value in zip(samples, cells[9:]):
                g = dict(zip(fmt, value.split(':')))
                ok = branch == 'MT' or (alt_gt(g.get('GT', '.')) and g.get('FT') == 'PASS'
                     and numeric(g.get('DP')) >= 10 and numeric(g.get('GQ')) >= 20)
                if ok:
                    good.append(sample)
            if not good:
                continue
            dst.write(line)
            counts['selected_records'] += 1
            counts[cells[0]] += 1
            for sample in good:
                counts['sample:' + sample] += 1
                for alt in cells[4].split(','):
                    selected.add((sample, cells[0], cells[1], cells[3], alt))
    return selected, dict(counts)


def meaningful(value):
    return value not in ('', '.', '[]', '{}', 'null', None)


def export_annotations(source, target, selected, branch):
    counts = Counter()
    populated = defaultdict(set)
    fields = ['SIFT','PolyPhen','FATHMM','CADD_PHRED','REVEL','SpliceAI','AlphaMissense_score',
              'clinvar_significance','clingen_evidence_json','af_max_observed']
    extra = ['gnomAD_MT_match_status','gnomAD_MT_AN','gnomAD_MT_AC_hom','gnomAD_MT_AC_het',
             'gnomAD_MT_AF_hom','gnomAD_MT_AF_het','Helix_match_status','Helix_AN',
             'Helix_counts_hom','Helix_counts_het','Helix_AF_hom','Helix_AF_het',
             'ClinVar_CLNSIG','ClinVar_CLNREVSTAT','ClinVar_CLNDN','ClinVar_ID'] if branch == 'MT' else []
    with gzip.open(source, 'rt', encoding='utf-8', newline='') as src, gzip.open(target, 'wt', encoding='utf-8', newline='') as dst:
        reader = csv.DictReader(src, delimiter='\t')
        writer = csv.DictWriter(dst, fieldnames=reader.fieldnames + extra, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for row in reader:
            counts['source_transcript_sample_rows'] += 1
            key = (row['sample'],row['chrom'],row['pos'],row['ref'],row.get('matched_alt',row['alt']))
            if key not in selected:
                continue
            # The cohort VCF retains failed sample genotypes: exclude them from this TSV.
            if branch == 'XY' and not (row['FT'] == 'PASS' and alt_gt(row['GT'])
                    and numeric(row['DP']) >= 10 and numeric(row['GQ']) >= 20):
                continue
            if branch == 'MT':
                ctx = json.loads(row['frequency_context_json'])
                csq = json.loads(row['csq_json'])
                for prefix, tag in [('gnomAD_MT','gnomAD_MT_reference_json'),('Helix','Helix_reference_json')]:
                    raw = ctx.get(tag)
                    data = json.loads(raw) if raw and raw not in ('null','{}') else {}
                    for field in extra:
                        if field.startswith(prefix + '_'):
                            suffix = field[len(prefix)+1:]
                            row[field] = ctx.get(field,'') if suffix == 'match_status' else data.get(suffix,'')
                for field in extra:
                    if field.startswith('ClinVar_'):
                        row[field] = csq.get('ClinVar' if field == 'ClinVar_ID' else field,'')
            writer.writerow(row)
            counts['selected_transcript_sample_rows'] += 1
            counts['sample:' + row['sample']] += 1
            for field in fields + extra:
                if meaningful(row.get(field)):
                    populated[field].add(key)
    return {'counts':dict(counts), 'populated_sample_alleles':{k:len(v) for k,v in populated.items()}}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    out=args.output
    out.mkdir(parents=True,exist_ok=False)
    report={'status':'PREPARED_LOCAL_ONLY','branches':{},'uploaded':False}
    reports=args.results/'Reports'
    xy_source=reports/'xy-germline-v1/cohort.xy.soft-filtered.vcf.gz'
    keys,counts=prepare_vcf(xy_source,out/'XY.quality.vcf','XY')
    report['branches']['XY']={'source':str(xy_source),'vcf':counts,
        'annotation':export_annotations(reports/'xy-supplementary-functional-views-v1/masters/cohort.master.tsv.gz',out/'XY.quality.annotations.tsv.gz',keys,'XY')}
    # Small independent chromosome inputs retain original cohort sample columns.
    for chrom in ['chrX','chrY']:
        with (out/'XY.quality.vcf').open(encoding='utf-8') as src, (out/(chrom+'.quality.vcf')).open('w',encoding='utf-8',newline='') as dst:
            for line in src:
                if line.startswith('#') or line.startswith(chrom+'\t'): dst.write(line)
    mt_keys=set()
    mt_counts={}
    samples=json.loads((reports/'mitochondrial-consequence-annotation-v1/manifest.json').read_text(encoding='utf-8'))['samples']
    for item in samples:
        sample=item['sample']
        source=args.results/'VCF/mitochondria-recovery-v2/checkpoints'/sample/'calls/mitochondria.filtered.vcf.gz'
        keys,counts=prepare_vcf(source,out/(sample+'.MT.site_PASS.vcf'),'MT')
        mt_keys.update(keys)
        mt_counts[sample]=counts
    report['branches']['MT']={'vcf':mt_counts,
        'annotation':export_annotations(reports/'mitochondrial-functional-views-v1/masters/cohort.master.tsv.gz',out/'MT.site_PASS.annotations.tsv.gz',mt_keys,'MT')}
    (out/'preparation-audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
