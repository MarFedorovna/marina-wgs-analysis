#!/usr/bin/env python3
"""REVEL archive only; downloading does not enable scores or accept a new license."""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.error
import urllib.request
import zipfile

CODE=Path(__file__).resolve().parents[1]


def load(path):return json.loads(path.read_text()) if path.is_file() else {}
def now():return dt.datetime.now(dt.timezone.utc)


def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    if path.is_symlink() or temporary.is_symlink():raise ValueError('Refusing symlink')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def run(cfg,root,execute=False):
    assert cfg['download_authorized_source_only'] is True and cfg['enable_variant_scores'] is False
    assert cfg['url']=='https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip'
    assert cfg['expected_bytes']<=cfg['max_bytes']<=700000000
    assert root.is_absolute() and root.name==cfg['bundle_id'] and not root.is_symlink()
    status=root/'status.json'; receipt=root/'resource-receipt.json'
    if not execute:return load(status) or dict(status='MISSING_SOURCE_ONLY_NOT_ENABLED')
    os.umask(0o077);root.mkdir(mode=0o700,parents=True,exist_ok=True)
    with (root/'.download.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return dict(status='ALREADY_RUNNING_SOURCE_ONLY')
        target=root/cfg['file'];partial=root/(cfg['file']+'.partial');identity_path=root/'source-identity.json'
        assert not target.is_symlink() and not partial.is_symlink()
        prior=load(receipt)
        if prior.get('status')=='DOWNLOADED_NOT_ENABLED_LICENSE_AND_ADAPTER_REVIEW_REQUIRED':
            assert target.is_file() and target.stat().st_size==cfg['expected_bytes']
            return prior
        assert not target.exists(),'Unreceipted target cannot be overwritten'
        deadline=min(dt.datetime.fromisoformat(cfg['deadline_utc']),now()+dt.timedelta(minutes=30))
        for attempt in range(1,5):
            try:
                assert now()<deadline,'Original deadline/runtime cap reached'
                with urllib.request.urlopen(urllib.request.Request(cfg['url'],method='HEAD'),timeout=30) as r:
                    identity=dict(url=r.url,bytes=int(r.headers['Content-Length']),last_modified=r.headers.get('Last-Modified'),etag=r.headers.get('ETag'))
                assert identity['url']==cfg['url'] and identity['bytes']==cfg['expected_bytes']
                old=load(identity_path);assert not old or old==identity,'Publisher identity changed; new review required'
                save(identity_path,identity)
                offset=partial.stat().st_size if partial.exists() else 0
                sha=hashlib.sha256();md5=hashlib.md5()
                if offset:
                    with partial.open('rb') as existing:
                        for block in iter(lambda:existing.read(1024*1024),b''):sha.update(block);md5.update(block)
                headers={}
                if offset:headers.update(Range='bytes=%d-'%offset,**{'If-Range':identity['last_modified'] or identity['etag'] or ''})
                save(status,dict(status='DOWNLOADING_SOURCE_ONLY_NOT_ENABLED',pid=os.getpid(),attempt=attempt,received_bytes=offset,expected_bytes=cfg['expected_bytes'],updated_utc=now().isoformat()))
                with urllib.request.urlopen(urllib.request.Request(cfg['url'],headers=headers),timeout=30) as r:
                    assert r.url==cfg['url'],'Unreviewed redirect'
                    if offset:assert r.status==206 and r.headers.get('Content-Range','').startswith('bytes %d-'%offset),'Safe resume unavailable; partial retained'
                    with partial.open('ab' if offset else 'xb') as out:
                        last=time.monotonic()
                        while True:
                            assert now()<deadline,'Original deadline/runtime cap reached'
                            block=r.read(1024*1024)
                            if not block:break
                            offset+=len(block);assert offset<=cfg['expected_bytes']
                            out.write(block);sha.update(block);md5.update(block)
                            if time.monotonic()-last>=5:
                                save(status,dict(status='DOWNLOADING_SOURCE_ONLY_NOT_ENABLED',pid=os.getpid(),attempt=attempt,received_bytes=offset,expected_bytes=cfg['expected_bytes'],updated_utc=now().isoformat()));last=time.monotonic()
                assert offset==cfg['expected_bytes'] and md5.hexdigest()==cfg['publisher_md5'],'Publisher checksum/size mismatch'
                assert zipfile.is_zipfile(partial),'Unexpected archive format'
                with zipfile.ZipFile(partial) as archive:
                    members=[dict(name=i.filename,compressed_bytes=i.compress_size,uncompressed_bytes=i.file_size) for i in archive.infolist()]
                    assert members and all(not Path(i['name']).is_absolute() and '..' not in Path(i['name']).parts for i in members)
                result=dict(status='DOWNLOADED_NOT_ENABLED_LICENSE_AND_ADAPTER_REVIEW_REQUIRED',source=cfg,bytes=offset,sha256=sha.hexdigest(),publisher_md5=md5.hexdigest(),retrieved_utc=now().isoformat(),transport_identity=identity,archive_members=members,extracted=False,patient_scores_enabled=False)
                partial.replace(target);save(receipt,result);save(status,{k:v for k,v in result.items() if k not in ['source','archive_members','transport_identity']});return result
            except urllib.error.HTTPError as e:
                if e.code not in (408,429) and e.code<500:
                    save(status,dict(status='BLOCKED_ACCESS_OR_REQUEST_REVIEW',http_status=e.code,updated_utc=now().isoformat()));raise
                save(status,dict(status='TRANSIENT_NETWORK_RETRY_SOURCE_ONLY',attempt=attempt,error=str(e),updated_utc=now().isoformat()))
            except (OSError,ValueError) as e:
                save(status,dict(status='TRANSIENT_NETWORK_RETRY_SOURCE_ONLY',attempt=attempt,error=str(e)[:200],updated_utc=now().isoformat()))
            except BaseException as e:
                save(status,dict(status='FAIL_REVIEW_REQUIRED',error=str(e)[:200],updated_utc=now().isoformat()));raise
            if attempt<4:time.sleep(min(15*attempt,max(0,(deadline-now()).total_seconds())))
        save(status,dict(status='FAIL_BOUNDED_RETRIES_EXHAUSTED',updated_utc=now().isoformat()))
        raise RuntimeError('REVEL bounded retries exhausted')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,default=CODE/'config/revel-download-only-v1.json');p.add_argument('--root',type=Path,default=Path('/srv/marina/resources/revel-source-only-v1'));p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(run(load(a.config),a.root,a.execute),indent=2))
