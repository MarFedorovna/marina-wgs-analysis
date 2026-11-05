#!/usr/bin/env python3
"""Isolated CPU runtime and wholly synthetic SpliceAI control; no patient inputs."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import urllib.request
import zipfile

PINS=['tensorflow-cpu==2.15.1','keras==2.15.0','numpy==1.26.4','h5py==3.10.0','pandas==2.2.3','pyfaidx==0.9.0.3','pysam==0.23.3','setuptools==75.8.2']
DEADLINE=dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
ROOT=Path('/srv/marina/resources/spliceai-research-runtime-v1')
def save(path,data):
 temp=path.with_name(path.name+'.tmp'); temp.write_text(json.dumps(data,indent=2)+'\n'); os.replace(temp,path)
def stamp():return dt.datetime.now(dt.timezone.utc).isoformat()
def command(args,limit=900):
 seconds=min(limit,int((DEADLINE-dt.datetime.now(dt.timezone.utc)).total_seconds()))
 if seconds<=0:raise TimeoutError('Original budget deadline')
 subprocess.run(args,check=True,timeout=seconds)
def main():
 import fcntl
 from download_external_reference_source import safe,sha256
 os.umask(0o077);safe(ROOT);ROOT.mkdir(exist_ok=True,parents=True)
 with (ROOT/'.runtime.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  status=ROOT/'status.json'
  def state(label,**extra):save(status,dict(status=label,pid=os.getpid(),updated_utc=stamp(),patient_inputs=False,**extra))
  try:
   purpose=json.loads((ROOT/'research-purpose-approval.json').read_text())
   if purpose.get('purpose')!='noncommercial_research' or purpose.get('source')!='explicit_user_confirmation':raise ValueError('Purpose confirmation required')
   if sys.version_info[:2]!=(3,10):raise ValueError('Pinned runtime requires Python 3.10')
   env=ROOT/'venv';wheelhouse=ROOT/'wheelhouse';wheelhouse.mkdir(exist_ok=True)
   if not env.exists():command([sys.executable,'-m','venv',str(env)])
   python=str(env/'bin/python');wheel=Path('/srv/marina/resources/spliceai-source-only-v1/spliceai-1.3.1-py2.py3-none-any.whl')
   if sha256(wheel)!='63c633b8de6803d4ffb613e6ef62f1a896fec9fc13245c0f000ff7adea279776':raise ValueError('Source wheel hash changed')
   state('DOWNLOADING_PINNED_CPU_WHEELS',top_level_pins=PINS)
   command([python,'-m','pip','download','--only-binary=:all:','--index-url','https://pypi.org/simple','--dest',str(wheelhouse),*PINS,str(wheel)],1200)
   inventory=[];requirements=[]
   for path in sorted(wheelhouse.glob('*.whl')):
    with zipfile.ZipFile(path) as archive:
     text=archive.read(next(x for x in archive.namelist() if x.endswith('.dist-info/METADATA'))).decode()
    name=next(x[6:] for x in text.splitlines() if x.startswith('Name: '));version=next(x[9:] for x in text.splitlines() if x.startswith('Version: '));digest=sha256(path)
    with urllib.request.urlopen(f'https://pypi.org/pypi/{name}/{version}/json',timeout=30) as r:metadata=json.load(r)
    official=next(x for x in metadata['urls'] if x['filename']==path.name)
    if official['digests']['sha256']!=digest or official['size']!=path.stat().st_size:raise ValueError('Publisher wheel hash mismatch: '+name)
    inventory.append(dict(name=name,version=version,filename=path.name,bytes=path.stat().st_size,sha256=digest,url=official['url'],publisher_sha256_verified=True))
    requirements.append(f'{name}=={version} --hash=sha256:{digest}')
   (ROOT/'requirements.lock.txt').write_text('\n'.join(requirements)+'\n');save(ROOT/'wheel-manifest.json',inventory)
   state('INSTALLING_VERIFIED_OFFLINE_CPU_RUNTIME',wheels=len(inventory))
   command([python,'-m','pip','install','--no-index','--find-links',str(wheelhouse),'--require-hashes','-r',str(ROOT/'requirements.lock.txt')],900)
   command([python,'-m','pip','check']);command([str(env/'bin/spliceai'),'--help'])
   state('RUNNING_SYNTHETIC_PUBLIC_CONTROL')
   command([python,str(Path(__file__).resolve()),'--synthetic-control'],300)
   control=json.loads((ROOT/'synthetic-control-receipt.json').read_text())
   if control['status']!='PASS_SYNTHETIC_RUNTIME_CONTROL':raise ValueError('Control failed')
   state('READY_LOCAL_RESEARCH_RUNTIME_PUBLIC_SYNTHETIC_CONTROL_PASS',clinical_pipeline_enabled=False,default_annotation='GENCODE_v24_bundled_NOT_GENCODE49',control_receipt=str(ROOT/'synthetic-control-receipt.json'),wheel_manifest_sha256=sha256(ROOT/'wheel-manifest.json'))
  except BaseException as exc:
   state('FAIL_RUNTIME_PREPARATION',error=str(exc)[:300]);raise
def synthetic():
 os.environ.update(CUDA_VISIBLE_DEVICES='-1',TF_NUM_INTEROP_THREADS='1',TF_NUM_INTRAOP_THREADS='2',OMP_NUM_THREADS='2',TF_CPP_MIN_LOG_LEVEL='2')
 import numpy as np
 from spliceai.utils import Annotator,get_delta_scores,one_hot_encode
 from types import SimpleNamespace
 from pyfaidx import Fasta
 # Deterministic artificial chromosome, no person or real-genome allele.
 sequence=('ACGT'*5000)
 reference=ROOT/'synthetic.fa'; reference.write_text('>SYNTHETIC\n'+sequence+'\n')
 Fasta(str(reference)).close()
 annotation=ROOT/'synthetic-annotation.tsv'
 annotation.write_text('#NAME\tCHROM\tSTRAND\tTX_START\tTX_END\tEXON_START\tEXON_END\nSYNTHETIC_CONTROL\tSYNTHETIC\t+\t100\t19900\t100,9900,\t10100,19900,\n')
 ann=Annotator(str(reference),str(annotation));pos=10001;ref=sequence[pos-1];alt='C' if ref!='C' else 'A'
 aligned=ann.ref_fasta['SYNTHETIC'][pos-1:pos].seq==ref
 if not aligned:raise ValueError('Synthetic reference mismatch')
 window=sequence[pos-5051:pos+5050];x=one_hot_encode(window)[None,:]
 probabilities=np.mean([model.predict(x,verbose=0) for model in ann.models],axis=0)
 if probabilities.shape!=(1,101,3) or not np.isfinite(probabilities).all() or np.any(probabilities<0) or np.any(probabilities>1) or not np.allclose(probabilities.sum(axis=-1),1,atol=1e-5):raise ValueError('Invalid model probabilities')
 record=SimpleNamespace(chrom='SYNTHETIC',pos=pos,ref=ref,alts=(alt,))
 result=get_delta_scores(record,ann,50,0)
 if len(result)!=1 or len(result[0].split('|'))!=10:raise ValueError('Missing delta scores')
 bad=get_delta_scores(SimpleNamespace(chrom='SYNTHETIC',pos=pos,ref='N',alts=(alt,)),ann,50,0)
 if bad:raise ValueError('Ref mismatch was not rejected')
 save(ROOT/'synthetic-control-receipt.json',dict(status='PASS_SYNTHETIC_RUNTIME_CONTROL',test_kind='wholly_synthetic_nonpatient_sequence',reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest(),reference_alignment=aligned,reference_mismatch_rejected=True,models_loaded=len(ann.models),probabilities_center=probabilities[0,50].tolist(),probability_output_shape=list(probabilities.shape),ds_fields=['DS_AG','DS_AL','DS_DG','DS_DL'],dp_fields=['DP_AG','DP_AL','DP_DG','DP_DL'],annotation=result,clinical_validation=False,patient_inputs=False,completed_utc=stamp()))
if __name__=='__main__':
 if '--synthetic-control' in sys.argv:synthetic()
 else:main()
