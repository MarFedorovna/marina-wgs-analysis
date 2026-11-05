#!/usr/bin/env python3
"""Native Windows SFTP with an existing managed key and exact pinned host key.

The private JSON configuration holds endpoints/paths; no credentials are copied.
Downloads are size-only, resumable in .part files and atomically finalized.
Only the public SSH host key is hashed for its trust fingerprint.
"""
import argparse
import base64
import datetime as dt
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import time
import paramiko


def require(ok, message):
    if not ok:
        raise ValueError(message)


def save(path, value):
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')
    os.replace(temporary,path)


def connect(config):
    client=paramiko.SSHClient()
    client.load_system_host_keys()
    client.load_host_keys(config['known_hosts'])
    keys=paramiko.HostKeys(config['known_hosts'])
    pinned=keys.lookup(config['pin_alias'])
    require(pinned is not None,'Existing pinned host alias required')
    wanted=[]
    for key in pinned.values():
        fingerprint='SHA256:'+base64.b64encode(hashlib.sha256(key.asbytes()).digest()).decode().rstrip('=')
        if fingerprint==config['expected_fingerprint']:
            wanted.append(key)
    require(len(wanted)==1,'Pinned public host-key fingerprint differs')
    key=wanted[0]
    address=config['host'] if config['port']==22 else f'[{config["host"]}]:{config["port"]}'
    # Reuse the SAME verified public key for the configured network endpoint,
    # in memory only. Never AutoAddPolicy and never persist new trust.
    client.get_host_keys().add(address,key.get_name(),key)
    client.set_missing_host_key_policy(paramiko.RejectPolicy())
    private=paramiko.PKey.from_path(config['identity_file'])
    client.connect(config['host'],port=config['port'],username=config['user'],pkey=private,
                   allow_agent=False,look_for_keys=False,timeout=15,banner_timeout=15,
                   auth_timeout=15,channel_timeout=60)
    require(client.get_transport().get_remote_server_key()==key,'Actual server key differs from pinned key')
    client.get_transport().set_keepalive(30)
    sftp=paramiko.SFTPClient.from_transport(client.get_transport(),window_size=8*2**20)
    sftp.get_channel().settimeout(60)
    return client,sftp


def blocks(remote,start,size):
    end=start+size
    requests=[(offset,min(32768,end-offset)) for offset in range(start,end,32768)]
    yield from remote.readv(requests,max_concurrent_prefetch_requests=32)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',required=True,help='Private endpoint/known-host/key paths and receipt-backed files')
    ap.add_argument('--probe-only',action='store_true')
    ap.add_argument('--probe-mib',type=int,default=4)
    ap.add_argument('--execute',action='store_true')
    args=ap.parse_args()
    cfg=json.loads(Path(args.config).read_text(encoding='utf-8-sig'))
    require(1<=args.probe_mib<=32 and 0<cfg.get('max_hours',8)<=8,'Bounded probe and runtime required')
    if not args.execute:
        print('DRY RUN: no network or key read');return
    destination=Path(cfg['download_directory']).resolve(strict=True)
    require(not any((p/'.git').exists() for p in (destination,*destination.parents)),'Destination cannot be Git')
    os.umask(0o077)
    status_path=Path(cfg['status_file'])
    deadline=time.monotonic()+cfg.get('max_hours',8)*3600
    for entry in cfg['files']:
        require(Path(entry['name']).name==entry['name'] and '\\' not in entry['name'] and ':' not in entry['name'],'Unsafe final basename')
        remote=PurePosixPath(entry['remote_path'])
        require(remote.is_absolute() and '..' not in remote.parts and str(remote).startswith('/srv/marina/outputs/'),'Only known completed project results allowed')
        require(isinstance(entry['size_bytes'],int) and entry['size_bytes']>0,'Expected receipt size required')
    client,sftp=connect(cfg)
    try:
        first=cfg['files'][0]
        actual=sftp.stat(first['remote_path'])
        require(actual.st_size==first['size_bytes'],'Remote result size differs from receipt')
        size=min(args.probe_mib*2**20,actual.st_size);start=time.monotonic();read=0
        with sftp.open(first['remote_path'],'rb') as remote:
            for block in blocks(remote,0,size):read+=len(block)
        require(read==size,'Incomplete bounded probe')
        probe={'Status':'PROBE_PASS','Bytes':read,'ElapsedSeconds':time.monotonic()-start,'PinnedHostKey':True,'PatientHashesComputed':False}
        probe['MiBPerSecond']=read/2**20/probe['ElapsedSeconds']
        save(status_path,probe);print(json.dumps(probe),flush=True)
        if args.probe_only:return
        for number,entry in enumerate(cfg['files'],1):
            final=destination/entry['name'];partial=final.with_name(final.name+'.part')
            provenance=partial.with_name(partial.name+'.source.json')
            if final.exists():
                require(not final.is_symlink() and final.stat().st_size==entry['size_bytes'],'Differing completed download; refusing overwrite')
                continue
            before=sftp.stat(entry['remote_path'])
            require(before.st_size==entry['size_bytes'],'Remote file size differs from expected')
            identity={'RemotePath':entry['remote_path'],'ExpectedBytes':entry['size_bytes'],'RemoteMtime':before.st_mtime,'PinnedHostKey':cfg['expected_fingerprint'],'Verification':'size_only'}
            if partial.exists():
                require(not partial.is_symlink() and provenance.is_file() and json.loads(provenance.read_text())==identity,'Unrelated partial retained; refusing to reuse')
            else:
                require(not provenance.exists(),'Stale source metadata without partial; review before new download')
                partial.touch(exist_ok=False);save(provenance,identity)
            offset=partial.stat().st_size
            require(offset<=entry['size_bytes'],'Partial exceeds expected size')
            require(shutil.disk_usage(destination).free>=entry['size_bytes']-offset+15*2**30,'15GiB reserve required')
            started=time.monotonic();initial=offset;last=started
            with sftp.open(entry['remote_path'],'rb') as remote,partial.open('ab') as local:
                while offset<entry['size_bytes']:
                    if cfg.get('pause_file') and Path(cfg['pause_file']).exists():
                        print('PAUSED; partial and manifest preserved',flush=True);return
                    require(time.monotonic()<deadline,'Bounded download deadline reached; partial preserved')
                    amount=min(8*2**20,entry['size_bytes']-offset);received=0
                    for block in blocks(remote,offset,amount):
                        local.write(block);received+=len(block)
                    require(received==amount,'Incomplete SFTP block; partial retained')
                    local.flush();offset+=received
                    now=time.monotonic()
                    if now-last>=5:
                        status={'Status':'DOWNLOADING','FileIndex':number,'FileCount':len(cfg['files']),'CompletedBytes':offset,'ExpectedBytes':entry['size_bytes'],'MiBPerSecond':(offset-initial)/2**20/max(.001,now-started),'PatientHashesComputed':False,'UpdatedUtc':dt.datetime.now(dt.timezone.utc).isoformat()}
                        save(status_path,status);print(json.dumps(status),flush=True);last=now
                os.fsync(local.fileno())
            after=sftp.stat(entry['remote_path'])
            require((before.st_size,before.st_mtime)==(after.st_size,after.st_mtime) and partial.stat().st_size==entry['size_bytes'],'Source changed or download incomplete')
            os.rename(partial,final)  # Windows rename refuses an existing target.
            complete={'Status':'COMPLETED_SIZE_ONLY','Name':entry['name'],'SizeBytes':entry['size_bytes'],'PatientHashesComputed':False,'CompletedUtc':dt.datetime.now(dt.timezone.utc).isoformat()}
            save(final.with_name(final.name+'.completed.json'),complete);save(status_path,complete);print(json.dumps(complete),flush=True)
        print('PASS all requested native SFTP downloads completed; SIZE ONLY',flush=True)
    finally:
        sftp.close();client.close()


if __name__=='__main__':
    main()
