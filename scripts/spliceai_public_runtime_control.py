#!/usr/bin/env python3
"""Official-model CPU control on a fixed public MYBPC3 allele; no patient input."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

def main(root,output):
 os.environ.update(CUDA_VISIBLE_DEVICES='-1',TF_NUM_INTEROP_THREADS='1',TF_NUM_INTRAOP_THREADS='2',OMP_NUM_THREADS='2',TF_CPP_MIN_LOG_LEVEL='2')
 import numpy as np
 from spliceai.utils import Annotator,get_delta_scores,one_hot_encode
 import tensorflow as tf
 reference=root/'references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/Homo_sapiens_assembly38.fasta';ann=Annotator(str(reference),'grch38');record=SimpleNamespace(chrom='chr11',pos=47332146,ref='T',alts=('G',));aligned=ann.ref_fasta['chr11'][47332145:47332146].seq=='T'
 if not aligned:raise ValueError('Fixed public MYBPC3 control REF mismatch')
 begin=time.monotonic();standard=get_delta_scores(record,ann,500,0);standard_time=time.monotonic()-begin
 original=ann.models
 class DirectModel:
  def __init__(self,model):self.model=model
  def predict(self,x):return self.model(x,training=False).numpy()
 sequence=ann.ref_fasta['chr11'][47332146-5501:47332146+5500].seq;x=one_hot_encode(sequence)[None,:];predict=original[0].predict(x,verbose=0);direct=original[0](x,training=False).numpy();error=float(np.max(np.abs(predict-direct)))
 if not np.allclose(predict,direct,atol=1e-7,rtol=1e-5):raise ValueError('Direct dispatch differs from official Keras predict')
 ann.models=[DirectModel(model) for model in original];timings=[];results=[]
 for _ in range(3):
  begin=time.monotonic();result=get_delta_scores(record,ann,500,0);timings.append(time.monotonic()-begin);results.append(result)
 if not standard or any(r!=standard for r in results):raise ValueError('Official DS/DP control differs with direct dispatch')
 result=dict(status='PASS_PUBLIC_SPLICEAI_CPU_DISPATCH_EQUIVALENCE',patient_inputs=False,network_requests=0,public_control=dict(chrom='chr11',pos=47332146,ref='T',alt='G',source='fixed public MYBPC3 MyVariant reference control'),ref_aligned=True,annotation=standard,model_count=5,distance=500,mask=0,standard_seconds=standard_time,direct_seconds=timings,median_direct_seconds=float(np.median(timings)),max_probability_absolute_error=error,tensorflow_version=tf.__version__,numpy_version=np.__version__,gpu_devices=[d.name for d in tf.config.list_physical_devices('GPU')],annotation_scope='GENCODE24_gene_union_NOT_VEP_GENCODE49_transcript_mapping',completed_utc=dt.datetime.now(dt.timezone.utc).isoformat())
 with output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path('/srv/marina'));p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.root,a.output)
