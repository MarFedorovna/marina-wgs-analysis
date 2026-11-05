"""Public-panel and public-reference helpers, with no patient data inputs."""
import bisect
import hashlib
import json
from pathlib import Path

PUBLIC_PANEL=Path('/srv/marina/resources/public-hcm-cds-v2')
PUBLIC_FASTA=Path('/srv/marina/references/broad-hg38-v0-altaware-g1575676516681666-prebuilt/Homo_sapiens_assembly38.fasta')
def panel():
 manifest=json.loads((PUBLIC_PANEL/'panel-manifest.json').read_text())
 if manifest.get('no_patient_selection') is not True or manifest.get('no_patient_variants_read') is not True or manifest.get('assembly')!='GRCh38' or manifest.get('gencode_release')!=49:raise ValueError('Fixed public GENCODE49 panel required')
 if manifest['source_clingen_sha256']!='5cc793e4c476776d7d60f326d8250b65198c65e576cfd771df28d320d02fb10e':raise ValueError('Public ClinGen provenance changed')
 bed=PUBLIC_PANEL/'hcm-cds-union-cadd.bed'
 if hashlib.sha256(bed.read_bytes()).hexdigest()!=manifest['bed_sha256']:raise ValueError('Public BED provenance mismatch')
 intervals=[('chr'+v[0].removeprefix('chr'),int(v[1]),int(v[2])) for v in (line.split('\t') for line in bed.read_text().splitlines())]
 if len(manifest['genes'])!=39 or any(start<0 or end<=start for chrom,start,end in intervals):raise ValueError('Expected fixed public39-gene interval panel')
 return manifest,intervals
class IntervalLookup:
 def __init__(self,intervals):
  self.ranges={};self.starts={}
  for chrom,start,end in intervals:self.ranges.setdefault(chrom,[]).append((start,end))
  for chrom,values in self.ranges.items():values.sort();self.starts[chrom]=[a for a,b in values]
 def contains(self,chrom,pos1):
  values=self.ranges.get(chrom,[]);index=bisect.bisect_right(self.starts.get(chrom,[]),pos1-1)-1
  return index>=0 and pos1-1<values[index][1]
class PublicFasta:
 def __init__(self,path=PUBLIC_FASTA):
  self.stream=path.open('rb');self.rows={}
  for line in path.with_name(path.name+'.fai').read_text().splitlines():
   v=line.split('\t');self.rows[v[0]]=tuple(map(int,v[1:5]))
 def sequence(self,chrom,pos1,length=1):
  candidates=[chrom,chrom.removeprefix('chr')]+(['MT','chrM','M'] if chrom in ['chrM','MT','M'] else [])
  key=next((x for x in candidates if x in self.rows),None)
  if key is None:raise ValueError('Public reference contig missing')
  total,offset,bases,width=self.rows[key]
  if not 1<=pos1<=pos1+length-1<=total:raise ValueError('Public reference bounds')
  result=bytearray();cursor=pos1-1
  while len(result)<length:
   amount=min(length-len(result),bases-cursor%bases);self.stream.seek(offset+cursor//bases*width+cursor%bases);result.extend(self.stream.read(amount));cursor+=amount
  return result.decode().upper()
 def close(self):self.stream.close()
def info_fields(text):return dict((part.split('=',1) if '=' in part else (part,True)) for part in text.split(';') if part)
def population_info_projection(text):
 """Keep every population/quality field; omit only redundant VEP consequence payload."""
 return ';'.join(part for part in text.split(';') if part.split('=',1)[0].lower() not in {'vep','csq'})
def population_fields(row):
 if len(row)<8:raise ValueError('Truncated public VCF record')
 chrom,pos,identifier,ref,alt,qual,filters,raw=row[:8];fields=info_fields(raw)
 alts=alt.split(',')
 for name in ['AC','AN','AF','nhomalt']:fields.setdefault(name,'.')
 def numbers(value,cast):return [None]*len(alts) if value=='.' else [None if x=='.' else cast(x) for x in value.split(',')]
 ac=numbers(fields['AC'],int);hom=numbers(fields['nhomalt'],int);af=numbers(fields['AF'],float);an=None if fields['AN']=='.' else int(fields['AN'])
 if not len(alts)==len(ac)==len(hom)==len(af) or (an is not None and an<0) or any(x is not None and (x<0 or (an is not None and x>an)) for x in ac) or any(h is not None and (h<0 or (a is not None and h*2>a)) for h,a in zip(hom,ac)) or any(f is not None and not 0<=f<=1 for f in af):raise ValueError('Invalid gnomAD allele/denominator arrays')
 return dict(chrom=chrom,pos=int(pos),ref=ref,alt=alt,filter=filters,AC=fields['AC'],AN=an,AF=fields['AF'],nhomalt=fields['nhomalt'],info=raw)
