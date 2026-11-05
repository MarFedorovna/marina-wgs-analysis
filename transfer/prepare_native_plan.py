#!/usr/bin/env python3
"""Prepare private native-download links for receipt-completed outputs over pinned SSH.

Creates only same-filesystem hardlinks inside an existing authenticated Jupyter
root. Does not change authentication, open ports, copy raw/public databases or
read browser cookies/tokens. Keep the generated plan private and outside Git.
"""
import argparse
import base64
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shlex
import subprocess
import time
from urllib.parse import quote


HELPER = r'''
import base64,datetime,hashlib,html,json,os,pathlib,re,sys
from urllib.parse import quote
os.umask(0o077)
cfg=json.loads(sys.argv[1]);root=pathlib.Path(cfg['root']).resolve(strict=True)
workspace=pathlib.Path(cfg['workspace']).resolve(strict=True);jroot=pathlib.Path(cfg['jupyter']).resolve(strict=True)
workspace.relative_to(jroot)
links=workspace/'marina-verified-downloads';links.mkdir(mode=0o700,exist_ok=True)
if links.is_symlink():raise RuntimeError('Native link directory must not be a symlink')
view=workspace/'RESULTS';view.mkdir(mode=0o700,exist_ok=True)
if view.is_symlink():raise RuntimeError('RESULTS must not be a symlink')
for group in ('BAM','gVCF','VCF','QC','Reports'):(view/group).mkdir(mode=0o700,exist_ok=True)
files={};receipts=[];names={};view_rows={}
def human_path(p):
 parts=p.relative_to(root).parts;name=p.name
 if name.endswith(('.bam','.bam.bai','.bai')) and len(parts)>=4 and parts[-2]=='align':
  return pathlib.Path('BAM',parts[-3],name),'full-reference alignment'
 if '.g.vcf.gz' in name:
  if 'gvcf' in parts:
   i=parts.index('gvcf');scope=parts[i+1];sample=parts[i+2]
   return (pathlib.Path('gVCF',sample,name) if scope=='autosomes' else pathlib.Path('gVCF',scope,sample,name)),scope
  return pathlib.Path('Pilot',*parts[1:]),'technical pilot; separate scope'
 if '.vcf' in name or name.endswith(('.bcf','.tbi','.csi')):
  if 'stages' in parts:return pathlib.Path('VCF',name),'completed cohort stage'
  return pathlib.Path('VCF','checkpoints',*parts[1:]),'partition/checkpoint'
 if name.endswith(('.html','.pdf','.xlsx','.docx','.csv','.tsv')):
  return pathlib.Path('Reports',*parts[1:]),'completed report/table'
 return pathlib.Path('QC',*parts[1:]),'completed QC/log/provenance'
def expose(p,entry):
 rel,scope=human_path(p);target=view/rel
 if target.is_symlink():raise RuntimeError('Human view cannot contain symlinks')
 if target.exists() and not os.path.samefile(p,target):
  rel=pathlib.Path('Versions',entry['Sha256'][:16],rel);target=view/rel
 target.parent.mkdir(mode=0o700,parents=True,exist_ok=True)
 if target.exists():
  if not os.path.samefile(p,target):raise RuntimeError('Refusing different human-view file')
 else:
  try:os.link(p,target)
  except FileExistsError:
   if not os.path.samefile(p,target):raise
 view_rows[entry['RelativePath']]={'ViewPath':rel.as_posix(),'Scope':scope,'Source':entry['RelativePath'],'SizeBytes':entry['SizeBytes'],'ExpectedServerSha256':entry['Sha256'],'Receipt':entry['SourceReceipt'],'DownloadFileName':entry['DownloadFileName']}
for rp in sorted((root/'outputs').rglob('*receipt.json')):
 if rp.name not in ('receipt.json','workflow-receipt.json'):continue
 if any(p.startswith('.') or 'partial' in p.lower() for p in rp.relative_to(root).parts):continue
 if rp.is_symlink() or root not in rp.resolve().parents:raise RuntimeError('Unsafe stage receipt')
 raw=rp.read_bytes();d=json.loads(raw)
 if rp.name=='workflow-receipt.json':
  if d.get('status')!='PASS' or not d.get('identity'):continue
  artifacts=d.get('artifacts',{})
 else:
  if not d.get('completed_utc') or not d.get('fingerprint'):continue
  artifacts=d.get('outputs',{})
 if not artifacts:continue
 receipts.append({'RelativePath':rp.relative_to(root).as_posix(),'SizeBytes':len(raw),'Sha256':hashlib.sha256(raw).hexdigest(),'ContentBase64':base64.b64encode(raw).decode()})
 for name,value in artifacts.items():
  rel=pathlib.PurePosixPath(name)
  if not rel.parts or rel.is_absolute() or '..' in rel.parts or '\\' in name:raise RuntimeError('Unsafe artifact path')
  p=rp.parent/name
  if p.is_symlink() or root not in p.resolve(strict=True).parents or not p.is_file():raise RuntimeError('Unsafe completed artifact')
  # Needed biological outputs, indices, QC, tables and reports; large public
  # resource caches and FASTQ/FASTA are deliberately absent.
  if not p.name.endswith(('.bam','.bai','.vcf','.vcf.gz','.bcf','.tbi','.csi','.tsv','.csv','.txt','.html','.pdf','.xlsx','.docx','.md','.json','.log')):continue
  sha=value['sha256'] if isinstance(value,dict) else value
  if not re.fullmatch('[0-9a-f]{64}',sha):raise RuntimeError('Invalid receipt digest')
  size=p.stat().st_size
  if isinstance(value,dict) and size!=value.get('bytes',value.get('size_bytes',size)):raise RuntimeError('Completed artifact changed size')
  source=p.relative_to(root).as_posix()
  if source in files and files[source]['Sha256']!=sha:raise RuntimeError('Conflicting receipts')
  filename='__'.join(p.relative_to(root).parts[1:])
  if len(filename.encode())>220:filename=hashlib.sha256(source.encode()).hexdigest()[:16]+'__'+p.name
  if filename in names and names[filename]!=source:raise RuntimeError('Native download filename collision')
  names[filename]=source;target=links/filename
  if target.is_symlink():raise RuntimeError('Native download target must be a hardlink, not symlink')
  if target.exists():
   if not os.path.samefile(p,target):raise RuntimeError('Refusing differing native download link')
  else:
   try:os.link(p,target)
   except FileExistsError:
    if not os.path.samefile(p,target):raise
  jrel=target.relative_to(jroot).as_posix()
  files[source]={'RelativePath':source,'ServerPath':str(p),'SizeBytes':size,'Sha256':sha,'SourceReceipt':rp.relative_to(root).as_posix(),'DownloadFileName':filename,'JupyterRelativePath':jrel}
  expose(p,files[source])
excluded_gvcfs=sum('.g.vcf' in pathlib.PurePosixPath(x['RelativePath']).name.lower() for x in files.values())
def local_selected(entry):
 name=pathlib.PurePosixPath(entry['RelativePath']).name.lower()
 parts=pathlib.PurePosixPath(entry['RelativePath']).parts
 if len(parts)>1 and parts[1].startswith('mitochondria'):
  if parts[1]!='mitochondria-recovery-v2' or len(parts)!=6 or parts[-2]!='calls':return False
  return name in {'mitochondria.filtered.vcf.gz','mitochondria.filtered.vcf.gz.tbi','mitochondria.pass.vcf.gz','mitochondria.pass.vcf.gz.tbi','mitochondria.raw.vcf.gz','mitochondria.raw.vcf.gz.tbi','research-report.json','numt-assessment.json','contamination.json','heteroplasmy.tsv','circular-coverage.tsv'}
 if '.g.vcf' in name:return False
 if '/checkpoints/' in entry['RelativePath'] and ('.vcf' in name or name.endswith(('.bcf','.tbi','.csi'))):return False
 return True
selected=[x for x in files.values() if local_selected(x)]
ordered=sorted(selected,key=lambda x:(not x['RelativePath'].endswith('.bam'),-x['SizeBytes']))
rows=sorted(view_rows.values(),key=lambda x:x['ViewPath'])
bams=sum(x['ViewPath'].startswith('BAM/') and x['ViewPath'].endswith('.bam') for x in rows)
gvcfs=sum(x['ViewPath'].startswith('gVCF/') and len(pathlib.PurePosixPath(x['ViewPath']).parts)==3 and x['ViewPath'].endswith('.g.vcf.gz') for x in rows)
cohort=any('/stages/04_joint/' in x['Source'] and x['ViewPath']=='VCF/cohort.raw.vcf.gz' for x in rows)
filtered=any('/stages/05_filter/' in x['Source'] and x['ViewPath']=='VCF/cohort.quality.vcf.gz' for x in rows)
status='Общий VCF готов' if cohort else 'Общего VCF пока нет: стадия объединения ещё не завершена'
body=['<!doctype html><html lang="ru"><meta charset="utf-8"><meta http-equiv="refresh" content="60"><title>Готовые результаты</title><style>body{font:16px system-ui;max-width:1300px;margin:32px auto;padding:0 20px}td,th{padding:8px;border-bottom:1px solid #ddd;text-align:left}table{border-collapse:collapse;width:100%}small{color:#555}</style>','<h1>Готовые результаты на сервере</h1>',f'<p>BAM: {bams}; gVCF образцов: {gvcfs} (текущий scope — chr1–22/autosomes). '+status+'.</p>','<p>Показаны файлы только из завершённых стадий с квитанциями. BAM — выравнивания; gVCF — пообразцовые вызовы; VCF — общий набор. Технические результаты ещё требуют аннотации и интерпретации.</p>','<p><a href="BAM/">BAM</a> · <a href="gVCF/">gVCF</a> · <a href="VCF/">VCF</a> · <a href="QC/">QC</a> · <a href="Reports/">Reports</a></p>','<table><tr><th>Файл</th><th>GiB</th><th>Область</th></tr>']
for row in rows:body.append('<tr><td><a href="../marina-verified-downloads/'+quote(row['DownloadFileName'],safe='')+'?download=1">'+html.escape(row['ViewPath'])+'</a></td><td>'+f'{row["SizeBytes"]/2**30:.4f}'+'</td><td>'+html.escape(row['Scope'])+'</td></tr>')
body.append('</table></html>')
temporary=view/('INDEX.html.'+str(os.getpid())+'.tmp');temporary.write_text('\n'.join(body),encoding='utf-8');os.replace(temporary,view/'INDEX.html')
view_manifest={'GeneratedUtc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'CompletedBAMs':bams,'AutosomalGVCFs':gvcfs,'CohortVCFReady':cohort,'FilteredCohortVCFReady':filtered,'Files':rows}
temporary=view/('manifest.json.'+str(os.getpid())+'.tmp');temporary.write_text(json.dumps(view_manifest,ensure_ascii=False,indent=2),encoding='utf-8');os.replace(temporary,view/'manifest.json')
print(json.dumps({'SchemaVersion':1,'GeneratedUtc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'ReceiptPolicy':'Completed PASS/fingerprint receipts; source size checked; local backup verification mode is recorded separately','LocalSelection':'BAM/index, cohort VCF/index, QC, reports; individual gVCFs stay on server','ExcludedIndividualGVCFArtifacts':excluded_gvcfs,'Files':ordered,'Receipts':receipts,'ServerView':{'Path':str(view),'CompletedBAMs':bams,'AutosomalGVCFs':gvcfs,'CohortVCFReady':cohort,'FilteredCohortVCFReady':filtered}}))
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True)
    parser.add_argument('--ssh-config', required=True)
    parser.add_argument('--known-hosts', required=True)
    parser.add_argument('--remote-root', required=True)
    parser.add_argument('--workspace-root', required=True)
    parser.add_argument('--jupyter-root', required=True)
    parser.add_argument('--output', required=True, help='Private ignored expected-identity plan')
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--watch', action='store_true')
    parser.add_argument('--interval', type=int, default=60)
    parser.add_argument('--max-hours', type=float, default=8)
    parser.add_argument('--deadline-utc', help='Absolute stop time for local plan refresh')
    parser.add_argument('--pause-flag', help='Persistent local user-pause marker to respect')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.host):
        raise ValueError('Confirmed SSH alias required')
    if args.interval < 60 or not 0 < args.max_hours <= 8:
        raise ValueError('Interval >=60 seconds and max-hours <=8 required')
    for value in (args.remote_root, args.workspace_root, args.jupyter_root):
        p = PurePosixPath(value)
        if not p.is_absolute() or '..' in p.parts or len(p.parts) < 3 or any(ord(c) < 32 for c in value):
            raise ValueError('Dedicated absolute remote directories required')
    if not args.execute:
        print('DRY RUN: no SSH, hardlinks or plan files created')
        return
    config, known = Path(args.ssh_config).expanduser().resolve(strict=True), Path(args.known_hosts).expanduser().resolve(strict=True)
    ssh = ['ssh', '-T', '-F', str(config), '-o', 'BatchMode=yes', '-o', 'StrictHostKeyChecking=yes',
           '-o', 'UserKnownHostsFile=' + str(known), '-o', 'ConnectTimeout=30',
           '-o', 'ServerAliveInterval=30', '-o', 'ServerAliveCountMax=4']
    command = shlex.join(['python3', '-', json.dumps({'root': args.remote_root, 'workspace': args.workspace_root,
                                                    'jupyter': args.jupyter_root})])
    deadline=time.monotonic()+args.max_hours*3600
    wall_deadline=time.time()+args.max_hours*3600
    if args.deadline_utc:
        import datetime as dt
        absolute=dt.datetime.fromisoformat(args.deadline_utc.replace('Z','+00:00'))
        if absolute.tzinfo is None:raise ValueError('Deadline must contain a UTC offset')
        deadline=min(deadline,time.monotonic()+max(0,absolute.timestamp()-time.time()))
        wall_deadline=min(wall_deadline,absolute.timestamp())
    while True:
        if time.monotonic()>=deadline or time.time()>=wall_deadline or (args.pause_flag and Path(args.pause_flag).exists()):
            print('PAUSED local plan refresh; server outputs unchanged',flush=True);break
        try:
            result = subprocess.run(ssh + [args.host, command], input=HELPER, text=True, capture_output=True,
                                    timeout=min(120,max(1,deadline-time.monotonic())))
        except subprocess.TimeoutExpired:
            if time.monotonic()>=deadline:break
            raise
        if result.returncode:
            raise RuntimeError(result.stderr)
        plan = json.loads(result.stdout)
        for receipt in plan['Receipts']:
            data = base64.b64decode(receipt['ContentBase64'], validate=True)
            if len(data) != receipt['SizeBytes'] or hashlib.sha256(data).hexdigest() != receipt['Sha256']:
                raise ValueError('Receipt metadata integrity check failed')
        for entry in plan['Files']:
            rel = quote(entry['JupyterRelativePath'], safe='/')
            entry['JupyterDownloadPath'] = '/files/' + rel + '?download=1'
            entry['JupyterUIPath'] = '/lab/tree/' + rel
        output = Path(args.output).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        temp = output.with_name(output.name + '.tmp')
        temp.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        temp.replace(output)
        print(f'PASS native plan: {len(plan["Files"])} completed artifacts; {len(plan["Receipts"])} receipts; view {plan["ServerView"]}',flush=True)
        if not args.watch or time.monotonic()+args.interval>=deadline:break
        time.sleep(min(args.interval,max(0,deadline-time.monotonic()),max(0,wall_deadline-time.time())))


if __name__ == '__main__':
    main()
