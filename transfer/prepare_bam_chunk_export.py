#!/usr/bin/env python3
"""Prepare bounded, resumable size-only BAM chunks in an existing protected file root.

No server/listener/authentication is created, and no BAM content hash is computed.
The private manifest must reach the PC separately through pinned managed SSH.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import time

CHUNK_BYTES = 64 * 2**20
RESERVE_BYTES = 15 * 2**30
DEADLINE = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def basename(name):
    require(isinstance(name, str) and name and Path(name).name == name and
            not any(c in name for c in '/\\:') and not any(ord(c) < 32 for c in name),
            'Unsafe basename')
    return name


def atomic_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def identity(path):
    stat = path.stat()
    return dict(size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                device=stat.st_dev, inode=stat.st_ino)


def make_manifest(cfg):
    source = Path(cfg['source_bam'])
    receipt_path = Path(cfg['source_receipt'])
    require(source.is_absolute() and not source.is_symlink() and source.is_file(), 'Regular absolute source BAM required')
    approved = Path(cfg['approved_source_root']).resolve(strict=True)
    resolved = source.resolve(strict=True)
    require(resolved.is_relative_to(approved) and resolved.suffix == '.bam', 'Source outside approved completed BAM root')
    require(not receipt_path.is_symlink() and receipt_path.is_file(), 'Completed source receipt required')
    receipt = json.loads(receipt_path.read_text())
    require(receipt.get('completed_utc') and receipt.get('fingerprint'), 'Source receipt is incomplete')
    old_hash = receipt.get('outputs', {}).get(source.name)
    require(isinstance(old_hash, str) and re.fullmatch(r'[0-9a-f]{64}', old_hash), 'Existing historical BAM checksum missing')
    original = identity(source)
    require(original['size_bytes'] == cfg['expected_bytes'] and original['size_bytes'] > 0, 'Receipt-backed expected size differs')
    name = basename(cfg['download_filename'])
    jupyter = Path(cfg['jupyter_root']).resolve(strict=True)
    export_base = Path(cfg['existing_protected_export_root']).resolve(strict=True)
    require(export_base.is_relative_to(jupyter) and not export_base.is_symlink(), 'Existing protected export root outside Jupyter root')
    metadata = dict(source_path=str(resolved), source_identity=original,
                    historical_source_sha256=old_hash, download_filename=name,
                    chunk_bytes=CHUNK_BYTES)
    # Hash of metadata only; never opens or hashes BAM sequence content.
    export_id = hashlib.sha256(json.dumps(metadata, sort_keys=True).encode()).hexdigest()[:24]
    folder = export_base / ('bam-chunks-' + export_id)
    chunks = []
    for index, offset in enumerate(range(0, original['size_bytes'], CHUNK_BYTES)):
        chunk_name = f'bam-{export_id}.chunk{index:05d}.bin'
        chunks.append(dict(index=index, offset=offset,
                           size_bytes=min(CHUNK_BYTES, original['size_bytes']-offset),
                           filename=chunk_name,
                           jupyter_relative_path=(folder / chunk_name).relative_to(jupyter).as_posix()))
    manifest = dict(schema_version=1, export_id=export_id,
                    verification='size_only', content_hashes_computed=False,
                    historical_source_sha256_preserved_not_recomputed=old_hash,
                    source_path=str(resolved), source_identity=original,
                    original_completed_utc=receipt['completed_utc'],
                    original_receipt_fingerprint=receipt['fingerprint'],
                    download_filename=name, size_bytes=original['size_bytes'],
                    chunk_bytes=CHUNK_BYTES, chunk_count=len(chunks), chunks=chunks,
                    protected_export_directory=str(folder),
                    jupyter_root=str(jupyter), authentication='existing Jupyter only; no new endpoint/token/cookie',
                    absolute_deadline_utc=DEADLINE.isoformat())
    return source, folder, manifest


def export(cfg, max_chunks=2, execute=False, max_runtime_seconds=1800):
    require(1 <= max_chunks <= 128, 'Bounded max-chunks1..128 required')
    require(1 <= max_runtime_seconds <= 1800, 'Maximum runtime must be1..1800seconds')
    source, folder, manifest = make_manifest(cfg)
    if not execute:
        return dict(status='DRY_RUN', export_id=manifest['export_id'],
                    source_bytes=manifest['size_bytes'], chunk_count=manifest['chunk_count'],
                    next_chunks=manifest['chunks'][:max_chunks], content_hashes_computed=False)
    require(dt.datetime.now(dt.timezone.utc) < DEADLINE, 'Original absolute VM deadline reached')
    os.umask(0o077)
    folder.mkdir(exist_ok=True)
    require(not folder.is_symlink(), 'Unsafe export directory')
    manifest_path = folder / 'private-chunk-manifest.json'
    if manifest_path.exists():
        require(json.loads(manifest_path.read_text()) == manifest, 'Existing export metadata differs')
    else:
        with manifest_path.open('x') as stream:
            json.dump(manifest, stream, indent=2)
    selected_missing=[c for c in manifest['chunks'] if not (folder/c['filename']).exists()][:max_chunks]
    require(shutil.disk_usage(folder).free >= sum(c['size_bytes'] for c in selected_missing)+RESERVE_BYTES,
            'Server requires planned new chunk bytes plus15GiB reserve')
    started=time.monotonic()
    made, ready = [], []
    for chunk in manifest['chunks']:
        final = folder / chunk['filename']
        if final.exists():
            require(not final.is_symlink() and final.stat().st_size == chunk['size_bytes'], 'Existing completed chunk has wrong size')
            ready.append(chunk['index'])
            continue
        if len(made) >= max_chunks:
            continue
        require(dt.datetime.now(dt.timezone.utc) < DEADLINE, 'Original absolute deadline reached')
        require(time.monotonic()-started < max_runtime_seconds, 'Bounded export runtime reached')
        require(identity(source) == manifest['source_identity'], 'Source BAM changed')
        partial = final.with_name(final.name + '.export.part')
        partial_state = partial.with_name(partial.name + '.json')
        state = dict(export_id=manifest['export_id'], chunk=chunk)
        if partial.exists():
            require(not partial.is_symlink() and partial_state.is_file() and
                    json.loads(partial_state.read_text()) == state, 'Unrelated export partial retained')
        else:
            require(not partial_state.exists(), 'Unrelated source state without partial')
            partial.touch(exist_ok=False)
            atomic_json(partial_state, state)
        offset = partial.stat().st_size
        require(offset <= chunk['size_bytes'], 'Oversized partial')
        require(shutil.disk_usage(folder).free >= chunk['size_bytes']-offset+RESERVE_BYTES, 'Server15GiB reserve required')
        with source.open('rb') as src, partial.open('ab') as out:
            src.seek(chunk['offset']+offset)
            while offset < chunk['size_bytes']:
                require(time.monotonic()-started < max_runtime_seconds and dt.datetime.now(dt.timezone.utc)<DEADLINE,
                        'Export deadline reached; resumable partial retained')
                payload = src.read(min(4*2**20, chunk['size_bytes']-offset))
                require(bool(payload), 'Unexpected source EOF')
                out.write(payload)
                offset += len(payload)
            out.flush()
            os.fsync(out.fileno())
        require(identity(source) == manifest['source_identity'], 'Source changed while exporting')
        require(partial.stat().st_size == chunk['size_bytes'], 'Chunk size differs')
        os.link(partial, final)  # Atomic no-replace, same filesystem.
        partial.unlink()
        made.append(chunk['index'])
        ready.append(chunk['index'])
    status = dict(status='BOUNDED_CHUNKS_READY', export_id=manifest['export_id'],
                  created_chunk_indices=made, ready_chunk_indices=ready,
                  manifest_path=str(manifest_path),
                  source_bytes=manifest['size_bytes'], chunk_count=manifest['chunk_count'],
                  ready_chunks=[c for c in manifest['chunks'] if c['index'] in ready],
                  content_hashes_computed=False, source_stat_unchanged=identity(source)==manifest['source_identity'],
                  all_chunks_ready=len(ready)==len(manifest['chunks']))
    status['elapsed_seconds']=round(time.monotonic()-started,3)
    atomic_json(folder / 'private-chunk-export-status.json', status)
    return status


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True,help='Private receipt-backed source and existing Jupyter export root')
    parser.add_argument('--max-chunks',type=int,default=2)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--max-runtime-seconds',type=int,default=1800)
    args=parser.parse_args()
    cfg=json.loads(args.config.read_text(encoding='utf-8-sig'))
    print(json.dumps(export(cfg,args.max_chunks,args.execute,args.max_runtime_seconds)))


if __name__=='__main__':main()
