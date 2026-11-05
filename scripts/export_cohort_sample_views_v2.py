#!/usr/bin/env python3
"""Reviewed joined-master extension of the unchanged cohort/sample views v1 engine."""
import argparse
import csv
import ctypes
import gzip
import hashlib
import html
import importlib.util
import itertools
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import uuid


def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def pin(path,guard):
    path=Path(path)
    if path.is_symlink() or not path.is_file():raise ValueError('Missing or symlink joined-master provenance input')
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):guard.check();h.update(block)
    return dict(path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest())


def check_pin(expected,guard):
    actual=pin(expected['path'],guard)
    if any(expected.get(k)!=actual[k] for k in ('path','bytes','sha256')):raise ValueError('Joined-master source/receipt pin differs')
    return actual


class DigestText:
    def __init__(self):self.digest=hashlib.sha256()
    def write(self,value):self.digest.update(value.encode('utf-8'));return len(value)


def verify_joined(manifest_path,workflow,guard):
    csv.field_size_limit(64*1024**2)
    manifest_path=Path(manifest_path).resolve();workflow=Path(workflow).resolve()
    expected_folder=workflow.parent/'combined-annotation-master-v1'
    if manifest_path!=expected_folder/'manifest.json':raise ValueError('Exact approved combined-master output root required')
    manifest_pin=pin(manifest_path,guard);value=load(manifest_path)
    required={'schema_version','status','base_master','enriched_master','original_fields','appended_fields','row_count','invariant_digest',
              'original_rows_fields_order_preserved','ranking_eligibility_unchanged','stage_receipts','upstream_receipts','created_utc'}
    if not required<=set(value) or value.get('schema_version')!=1 or value.get('status')!='PASS_COMBINED_ANNOTATION_MASTER':raise ValueError('Completed reviewed joined-master schema required')
    if value['original_rows_fields_order_preserved'] is not True or value['ranking_eligibility_unchanged'] is not True:raise ValueError('Joined-master original/ranking invariants required')
    base=workflow/'stages/07_report/master_evidence.tsv.gz'
    if value['base_master'].get('path')!=str(base):raise ValueError('Joined base must be exact original stage7 master')
    stage7=load(base.parent/'workflow-receipt.json');expected=stage7['artifacts'][base.name]
    if stage7.get('status')!='PASS' or not stage7.get('identity') or any(value['base_master'].get(k)!=expected[k] for k in ('bytes','sha256')):raise ValueError('Joined base pin differs from stage7 receipt')
    base_pin=check_pin(value['base_master'],guard)
    enriched=expected_folder/'master_evidence.enriched.tsv.gz'
    if value['enriched_master'].get('path')!=str(enriched):raise ValueError('Unexpected enriched-master artifact path')
    enriched_pin=check_pin(value['enriched_master'],guard)
    stage_paths={str(workflow/'stages'/s/'workflow-receipt.json') for s in ('06_annotate','07_report')}
    if {x.get('path') for x in value['stage_receipts']}!=stage_paths or len(value['stage_receipts'])!=2:raise ValueError('Exact stage6/7 receipt pins required')
    upstream_paths={str(workflow.parent/'tissue-protein-context-v1/manifest.json'),str(workflow.parent/'population-reference-context-v1/workflow-receipt.json')}
    if {x.get('path') for x in value['upstream_receipts']}!=upstream_paths or len(value['upstream_receipts'])!=2:raise ValueError('Exact tissue/population completion receipt pins required')
    pins=[manifest_pin,base_pin,enriched_pin]
    for item in value['stage_receipts']+value['upstream_receipts']:pins.append(check_pin(item,guard))
    tissue=load(workflow.parent/'tissue-protein-context-v1/manifest.json')
    population=load(workflow.parent/'population-reference-context-v1/workflow-receipt.json')
    if tissue.get('status')!='PASS_SEPARATE_RESEARCH_REFERENCE_CONTEXT_EXPORT' or tissue.get('original_pipeline_changed') is not False or tissue.get('network_requests')!=0:raise ValueError('Tissue context is not a completed unchanged-primary sidecar')
    if population.get('status')!='PASS' or population.get('stage')!='population-reference-context-v1' or population.get('metadata',{}).get('status')!='PASS_LOCAL_POPULATION_MT_CONTEXT_EXPORT':raise ValueError('Population context is not a completed reviewed sidecar')
    receipt_path=expected_folder/'workflow-receipt.json';receipt=load(receipt_path);pins.append(pin(receipt_path,guard))
    if receipt.get('status')!='PASS' or receipt.get('stage')!='combined-annotation-master-v1' or not re.fullmatch('[0-9a-f]{64}',receipt.get('identity','')):raise ValueError('Joined-master own completion receipt required')
    for name,item in [('manifest.json',manifest_pin),(enriched.name,enriched_pin)]:
        expected=receipt.get('artifacts',{}).get(name,{})
        if any(expected.get(k)!=item[k] for k in ('bytes','sha256')):raise ValueError('Joined workflow receipt artifact digest differs')
    original=value['original_fields'];appended=value['appended_fields']
    if not original or not appended or len(set(original+appended))!=len(original+appended) or len(appended)>256 or any(not re.fullmatch('allref_[A-Za-z0-9_]{1,120}',f) for f in appended):raise ValueError('Only unique namespaced appended context fields accepted')
    digest_text=DigestText();digest_writer=csv.writer(digest_text,delimiter='\t',lineterminator='\n')
    with gzip.open(base,'rt',encoding='utf-8',newline='') as base_stream,gzip.open(enriched,'rt',encoding='utf-8',newline='') as joined_stream:
        left=csv.reader(base_stream,delimiter='\t');right=csv.reader(joined_stream,delimiter='\t')
        if next(left)!=original or next(right)!=original+appended:raise ValueError('Original column identity/order or append-only header differs')
        digest_writer.writerow(original);rows=0
        for rows,(before,after) in enumerate(itertools.zip_longest(left,right),1):
            if rows%10000==1:guard.check()
            if before is None or after is None or len(before)!=len(original) or len(after)!=len(original)+len(appended) or after[:len(original)]!=before:
                raise ValueError('Joined master changed an original row/field/order or row count')
            digest_writer.writerow(before)
    if type(value['row_count']) is not int or rows!=value['row_count'] or digest_text.digest.hexdigest()!=value['invariant_digest']:raise ValueError('Joined-master original-row invariant digest/count differs')
    return value,pins


def promote(source,target):
    libc=ctypes.CDLL(None,use_errno=True)
    if not hasattr(libc,'renameat2'):raise ValueError('Atomic no-replace Linux promotion required')
    if libc.renameat2(-100,os.fsencode(source),-100,os.fsencode(target),1)!=0:
        error=ctypes.get_errno();raise OSError(error,os.strerror(error),str(target))


def copy_applied_evidence(joined,staging,guard):
    """Portable copies of declared applied result tables, never source databases."""
    folder=Path(joined['enriched_master']['path']).parent
    receipt=load(folder/'workflow-receipt.json');population_path=folder.parent/'population-reference-context-v1/workflow-receipt.json'
    population=load(population_path);expected_upstream=pin(population_path,guard);copied={}
    def copy_one(relative,expected):
        rel=PurePosixPath(relative)
        if rel.is_absolute() or '..' in rel.parts or rel.as_posix()!=relative or '\\' in relative:raise ValueError('Unsafe applied evidence path')
        source=folder.joinpath(*rel.parts);target=staging.joinpath(*rel.parts)
        if any(p.is_symlink() for p in [source,*source.parents]):raise ValueError('Symlink applied evidence source')
        if any(receipt.get('artifacts',{}).get(relative,{}).get(k)!=expected.get(k) for k in ('bytes','sha256')):raise ValueError('Applied evidence differs from joined completion receipt')
        check_pin(dict(expected,path=str(source)),guard);target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
        with source.open('rb') as inp,target.open('xb') as out:
            for block in iter(lambda:inp.read(8*1024*1024),b''):
                guard.check()
                if shutil.disk_usage(staging).free<64*1024**3:raise ValueError('Applied evidence copy reached sidecar disk reserve')
                out.write(block)
        checked=check_pin(dict(expected,path=str(target)),guard);copied[relative]={k:checked[k] for k in ('bytes','sha256')}
    for entry in joined.get('auxiliary_artifacts',[]):
        relative=entry['relative_path'];source_relative=entry.get('source_relative_path','')
        if relative!='applied-evidence/population/'+source_relative or source_relative=='cohort/nuclear-sample-rows.tsv.gz':raise ValueError('Only applied population result evidence can accompany views')
        if entry.get('path')!=str(folder/relative) or entry.get('source_upstream_receipt')!=expected_upstream:raise ValueError('Applied evidence producer provenance differs')
        source_expected=population.get('artifacts',{}).get(source_relative,{})
        if any(entry.get(k)!=source_expected.get(k) for k in ('bytes','sha256')):raise ValueError('Applied evidence identity differs from population receipt')
        role='MT_APPENDIX' if 'MT-' in source_relative or '/MT-QC/' in source_relative else 'APPLIED_NUCLEAR_REFERENCE_LONG_TABLE'
        if entry.get('role')!=role:raise ValueError('Unknown/mismatched applied evidence role')
        # Reject corpus/database paths even if a manifest claims they are applied.
        cohort={'nuclear-population-long','nuclear-exome-AN-v41','ClinVar-validated-allele-map','ClinVar-SCV-long','ClinVar-citation-links','gene-coverage-context-long','MT-population-reference-long'}
        if not (source_relative in {'cohort/'+name+'.tsv.gz' for name in cohort} or re.fullmatch(r'samples/[A-Za-z0-9_.-]+/(MT-population-reference\.tsv\.gz|MT-QC/[A-Za-z0-9_.-]+\.txt)',source_relative)):
            raise ValueError('Unrecognized applied result-table path')
        copy_one(relative,{k:entry[k] for k in ('bytes','sha256')})
    if copied:
        for relative in ('auxiliary-index.json','applied-evidence/MT-appendix-index.json'):
            copy_one(relative,receipt['artifacts'][relative])
        page=staging/'index.html'
        with page.open('a',encoding='utf-8') as stream:
            stream.write('<details><summary>Applied reference evidence and separate MT appendix</summary><ul>')
            for relative in sorted(copied):stream.write('<li><a href="'+html.escape(relative,quote=True)+'">'+html.escape(relative)+'</a></li>')
            stream.write('</ul></details>\n')
        copied['index.html']={k:pin(page,guard)[k] for k in ('bytes','sha256')}
    return copied


def export(workflow,output,ranking_module,ranking_sha256,base_exporter,base_exporter_sha256,joined_manifest,*,deadline=None,synthetic=False):
    if hasattr(os,'sched_getaffinity'):os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
    try:
        import resource
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=4*1024**3
        resource.setrlimit(resource.RLIMIT_AS,(min(cap,soft) if soft>0 else cap,hard))
    except ImportError:pass
    base_exporter=Path(base_exporter).resolve()
    if hashlib.sha256(base_exporter.read_bytes()).hexdigest()!=base_exporter_sha256:raise ValueError('Frozen views-v1 engine SHA differs')
    spec=importlib.util.spec_from_file_location('_reviewed_views_v1_private_instance',base_exporter)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    guard=module.Guard(deadline or module.ABSOLUTE_DEADLINE);guard.check()
    joined,pins=verify_joined(joined_manifest,workflow,guard)
    original_inputs=module.checked_inputs
    def joined_inputs(root,current_guard):
        sources,source_pins=original_inputs(root,current_guard)
        sources['master']=Path(joined['enriched_master']['path'])
        source_pins.update({p['path']:{k:p[k] for k in ('bytes','sha256')} for p in pins})
        return sources,source_pins
    # Extend only this freshly loaded module instance. The frozen v1 file and all
    # active controllers stay unchanged. Ranking/eligibility code is reused verbatim.
    module.checked_inputs=joined_inputs
    module.PROJECTION=list(module.PROJECTION)+joined['appended_fields']
    module.SUMMARY=list(module.SUMMARY)+joined['appended_fields']
    output=Path(output).resolve();workflow=Path(workflow).resolve()
    if output.exists() or output==workflow or workflow in output.parents:raise ValueError('New canonical sidecar output outside primary workflow required')
    staging=output.parent/('.'+output.name+'.joined-partial-'+uuid.uuid4().hex)
    result=module.export(workflow,staging,ranking_module,ranking_sha256,deadline=guard.deadline,synthetic=synthetic)
    applied=copy_applied_evidence(joined,staging,guard)
    result['artifacts'].update(applied)
    result.update(report_input_kind='COMBINED_ANNOTATION_MASTER_V1',
        joined_master=dict(manifest=pin(joined_manifest,guard),enriched_master=joined['enriched_master'],base_master=joined['base_master'],
                           original_fields=joined['original_fields'],appended_fields=joined['appended_fields'],row_count=joined['row_count'],invariant_digest=joined['invariant_digest'],
                           original_rows_fields_order_verified=True,rank_eligibility_columns_unchanged_verified=True),
        views_v1_engine=dict(path=str(base_exporter),sha256=base_exporter_sha256),views_v2_exporter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        copied_applied_evidence=sorted(k for k in applied if k!='index.html'))
    for item in pins:check_pin(item,guard)
    (staging/'manifest.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    promote(staging,output)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('workflow-root','output','ranking-module','base-exporter','joined-manifest'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ranking-sha256',required=True);p.add_argument('--base-exporter-sha256',required=True);a=p.parse_args()
    print(json.dumps(export(a.workflow_root,a.output,a.ranking_module,a.ranking_sha256,a.base_exporter,a.base_exporter_sha256,a.joined_manifest),indent=2))
