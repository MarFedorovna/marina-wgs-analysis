#!/usr/bin/env python3
"""Explicit verified reuse of an owned extracted public cache stage; no auto retry or queue changes."""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import tarfile
import time
import prepare_annotation_resources as prep


def now():return dt.datetime.now(dt.timezone.utc)


def write_status(directory,status,**extra):
    directory.mkdir(mode=0o700,parents=True,exist_ok=True);path=directory/'status.json';tmp=directory/'status.tmp'
    assert not path.is_symlink() and not tmp.is_symlink()
    value=dict(status=status,pid=os.getpid(),updated_utc=now().isoformat(),deadline_utc='2026-11-06T03:46:48.046025+00:00',**extra)
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path);print(json.dumps(value),flush=True)


def stream_sha(stream):
    value=hashlib.sha256()
    for block in iter(lambda:stream.read(8*1024*1024),b''):value.update(block)
    return value.hexdigest()


def verify_extracted(archive,cache,limit,progress=None):
    """Compare every extracted regular file's exact bytes to the checksum-verified tar."""
    entries=[];names=set();total=0
    prep.need(cache.is_dir() and not cache.is_symlink(),'Owned extracted cache directory absent/symlink')
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar:
            relative=Path(member.name)
            prep.need(not relative.is_absolute() and '..' not in relative.parts and relative.parts and relative.parts[0]=='homo_sapiens','Unsafe archive path')
            normalized=relative.as_posix();prep.need(normalized not in names and (member.isfile() or member.isdir()),'Duplicate/link/special archive member')
            names.add(normalized);path=cache/relative
            prep.need(not path.is_symlink() and path.resolve(strict=True).is_relative_to(cache.resolve()),'Extracted path symlink/escape')
            if member.isdir():prep.need(path.is_dir(),'Extracted directory missing');continue
            prep.need(path.is_file() and path.stat().st_size==member.size,'Extracted file missing/truncated')
            total+=member.size;prep.need(total<=limit,'Expanded cache exceeds unchanged configured limit')
            with tar.extractfile(member) as source:expected=stream_sha(source)
            actual=prep.sha256(path);prep.need(actual==expected,'Extracted file differs from verified archive: '+normalized)
            entries.append(dict(path='cache/'+normalized,bytes=member.size,sha256=actual))
            if progress and len(entries)%250==0:progress(len(entries),total)
    actual=set()
    for path in cache.rglob('*'):
        prep.need(not path.is_symlink(),'Unexpected extracted symlink')
        if path.is_file():actual.add(path.relative_to(cache).as_posix())
        else:prep.need(path.is_dir(),'Unexpected extracted special file')
    expected={e['path'].removeprefix('cache/') for e in entries}
    prep.need(actual==expected and entries,'Extracted tree has missing or unexpected files')
    return dict(status='PASS_EXACT_ARCHIVE_TO_EXTRACTED_FILES',files=len(entries),expanded_bytes=total,
        tree_entries=sorted(entries,key=lambda i:i['path']),
        tree_entries_sha256=hashlib.sha256(json.dumps(sorted(entries,key=lambda i:i['path']),sort_keys=True).encode()).hexdigest())


def recover(config,root,stage,archive_sha,status_dir,execute=False,max_minutes=30):
    if not execute:
        return dict(status='REVIEW_ONLY_NOT_EXECUTED',config=str(config),root=str(root),resume_extracted_stage=str(stage),expected_archive_sha256=archive_sha,no_download_or_extraction_requested=True)
    os.umask(0o077);prep.need(0<max_minutes<=30 and re.fullmatch('[0-9a-f]{64}',archive_sha or ''),'Bounded runtime and expected archive SHA required')
    deadline=min(dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00'),now()+dt.timedelta(minutes=max_minutes))
    prep.need(now()<deadline,'Original budget allowance exhausted')
    def timeout(*_):raise TimeoutError('Original VM deadline or explicit recovery runtime bound reached')
    config=config.resolve(strict=True);root=root.resolve(strict=True);stage=stage.absolute()
    final_cfg=json.loads(config.read_text());base_config=(config.parent/final_cfg['base_config']).resolve(strict=True) if final_cfg.get('base_config') else config
    cfg=json.loads(base_config.read_text());target=root/cfg['bundle_id']
    prep.need(stage.parent==root and re.fullmatch(r'\.partial-annotation-[0-9a-f]{32}',stage.name) and not stage.is_symlink(),'Use the explicitly named owned stage inside annotation root')
    prep.need(stage.is_dir() and not target.exists(),'Owned stage absent or immutable base target already exists; inspect state instead of overwriting')
    prep.need((stage/'resources.lock.json').read_bytes()==base_config.read_bytes(),'Resume stage recipe differs byte-for-byte from the exact base config')
    sources=[i for i in cfg['resources'] if i.get('enabled',False)];cache_item=next(i for i in sources if i['name']=='vep_cache');archive=stage/cache_item['destination']
    expected_names={'resources.lock.json','vep-version.txt','cache',*[i['destination'] for i in sources]}
    prep.need({p.name for p in stage.iterdir()}==expected_names,'Unexpected owned stage files; review before reuse')
    status_dir.mkdir(mode=0o700,parents=True,exist_ok=True)
    signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,(deadline-now()).total_seconds())
    try:
        with (root/'.annotation-prep.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            write_status(status_dir,'VERIFYING_PINNED_TOOL_RECIPE_AND_ARCHIVE',stage=str(stage),base_recipe_sha256=prep.sha256(base_config),code_sha256=prep.sha256(Path(__file__)))
            prep.need(shutil.disk_usage(root).free>=max(150,cfg.get('min_free_gib',150))*1024**3,'Analysis free-space reserve below unchanged 150GiB minimum')
            fresh_help=prep.ensure_vep_image(cfg,pull=False)
            prep.need(archive.stat().st_size==cache_item['bytes'] and prep.sha256(archive)==archive_sha,'Completed archive identity differs')
            publisher=prep.publisher_check(archive,cache_item)
            write_status(status_dir,'VERIFYING_ALL_EXTRACTED_FILES_AGAINST_ARCHIVE',archive_sha256=archive_sha,archive_bytes=archive.stat().st_size)
            checked=verify_extracted(archive,stage/'cache',cfg.get('max_cache_expanded_gib',128)*1024**3,
                lambda count,total:write_status(status_dir,'VERIFYING_ALL_EXTRACTED_FILES_AGAINST_ARCHIVE',verified_files=count,verified_expanded_bytes=total))
            downloaded=[]
            for item in sources:
                file=stage/item['destination']
                if item['name']=='vep_cache':
                    downloaded.append(dict(name=item['name'],version=item['version'],source=item['url'],path=file.name,
                        bytes=file.stat().st_size,sha256=archive_sha,status='available',publisher_validation='BSDsum+bytes; archive SHA256 and exact extracted bytes revalidated',retrieved_utc=now().isoformat()))
                else:
                    cached=prep.source_cache_directory(root,item)/'resource.json';record=json.loads(cached.read_text())
                    prep.need(record['name']==item['name'] and record['version']==item['version'] and record.get('source')==item['url'],'Small source receipt identity differs')
                    prep.need(file.is_file() and not file.is_symlink() and file.stat().st_size==record['bytes'] and prep.sha256(file)==record['sha256'],'Extracted-stage small source differs from original receipt')
                    prep.publisher_check(file,item);downloaded.append(record)
            available={i['name']:i for i in downloaded}
            for item in cfg.get('optional_adapters',[]):available.setdefault(item['name'],dict(item))
            info=stage/'cache/homo_sapiens/115_GRCh38/info.txt';available.update(prep.cache_metadata_records(info,stage))
            validation=dict(status=checked['status'],files=checked['files'],expanded_bytes=checked['expanded_bytes'],tree_entries_sha256=checked['tree_entries_sha256'],archive_sha256=archive_sha,
                publisher_validation=publisher,fresh_tool_help_sha256=hashlib.sha256(fresh_help.encode()).hexdigest(),
                source_recipe_sha256=prep.sha256(base_config),original_stage=str(stage),validated_utc=now().isoformat(),no_archive_download=True,no_reextraction=True)
            (stage/'extracted-resume-validation.json').write_text(json.dumps(validation,indent=2)+'\n')
            write_status(status_dir,'HASHING_COMPLETE_PUBLIC_RESOURCE_ARTIFACT_TREE',verified_files=checked['files'],expanded_bytes=checked['expanded_bytes'])
            manifest=dict(schema_version=1,status='PASS',assembly='GRCh38',config_sha256=prep.sha256(base_config),vep=cfg['vep'],cache_info=info.read_text(),cache_expanded_bytes=checked['expanded_bytes'],resources=available,
                artifacts=prep.tree_manifest(stage),completed_utc=now().isoformat(),patient_variants_sent_to_remote_services=False,
                extracted_resume_validation=validation,metadata_correction='Missing cache source_MANE is unreported; external MANE summary remains a distinct source')
            (stage/'resource-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
            prep.need(now()<deadline,'Recovery bound reached before atomic base promotion')
            prep.promote(stage,target)
        # Original final preparation code re-verifies immutable base artifacts, then
        # reuses hardlinks and the already-ready small source receipts. No re-extraction.
        write_status(status_dir,'PREPARING_FINAL_BUNDLE_FROM_VERIFIED_BASE',base_manifest=str(target/'resource-manifest.json'))
        final=prep.prepare(config,root,execute=True,pull_image=False,prefetched_cache=target/cache_item['destination'],prefetched_cache_sha256=archive_sha)
        write_status(status_dir,'PASS_REAL_RESOURCE_MANIFEST_READY_AWAITING_PUBLIC_VEP_CONTROLS',final_manifest=str(final/'resource-manifest.json'),final_manifest_sha256=prep.sha256(final/'resource-manifest.json'),
            final_recipe_sha256=prep.sha256(config),base_recipe_sha256=prep.sha256(base_config),cache_MANE_source_version_reported=False,external_MANE_summary_separate=True,
            no_patient_analysis_started_by_recovery=True,existing_controller_not_changed=True)
        return json.loads((status_dir/'status.json').read_text())
    except BaseException as error:
        write_status(status_dir,'FAILED_RECOVERY_NO_AUTOMATIC_RETRY',error_type=type(error).__name__,error=str(error)[:300],preserved_stage_or_base=str(stage if stage.exists() else target));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',type=Path,required=True);p.add_argument('--root',type=Path,default=Path('/srv/marina/resources/annotation'));p.add_argument('--resume-extracted-stage',type=Path,required=True);p.add_argument('--archive-sha256',required=True);p.add_argument('--status-dir',type=Path,default=Path('/srv/marina/private/annotation-resource-recovery-v1'));p.add_argument('--max-minutes',type=float,default=30);p.add_argument('--execute',action='store_true');a=p.parse_args();print(json.dumps(recover(a.config,a.root,a.resume_extracted_stage,a.archive_sha256,a.status_dir,a.execute,a.max_minutes),indent=2))
