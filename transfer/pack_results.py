#!/usr/bin/env python3
"""Create one ZIP64 of completed friendly result folders, without file SHA hashing.

Dry-run by default. Execute only after the cohort VCF has been accepted. The
archive gets a standalone INDEX and a size/verification manifest. Content stores,
staging, raw reads, individual gVCFs, public databases and credentials are excluded.
"""
import argparse
import datetime as dt
import html
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import zipfile
from urllib.parse import quote

VIEWS = ('BAM', 'VCF', 'QC', 'Reports', 'Settings', 'Code', 'Logs')
EXCLUDED_DIRS = {'objects', 'size-only', 'snapshots', '.staging', 'gvcf', 'raw', 'raw_data',
                 'resources', 'references', 'cache', 'caches', '.git', '.ssh', '.brev', '.aws', '.azure'}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def local_path(value):
    if os.name != 'nt' and re.match(r'^[A-Za-z]:[\\/]', value):
        value = subprocess.check_output(['wslpath', '-u', value], text=True).strip()
    return Path(value).expanduser().resolve()


def collect(root):
    history = {}
    state_path = root / 'verified-objects.json'
    state = json.loads(state_path.read_text(encoding='utf-8-sig')) if state_path.is_file() else {}
    p = root / 'verified-files.json'
    if p.is_file():
        for entry in json.loads(p.read_text(encoding='utf-8-sig')).values():
            identity = state.get(entry['Sha256'], {})
            history[entry['LocalPath']] = {'VerificationAtImport': 'sha256', 'HistoricalSHA256': entry['Sha256'],
                                         'ExpectedSizeBytes': entry['SizeBytes'], 'ExpectedMtimeNs': identity.get('mtime_ns')}
    p = root / 'size-only-files.json'
    if p.is_file():
        for entry in json.loads(p.read_text(encoding='utf-8-sig')).values():
            history[entry['LocalPath']] = {'VerificationAtImport': 'size_only', 'HistoricalSHA256': None,
                                         'ExpectedServerSha256': entry['ExpectedServerSha256'],
                                         'ExpectedSizeBytes': entry['SizeBytes'], 'ExpectedMtimeNs': entry['MtimeNs']}
    rows = []
    for group in VIEWS:
        folder = root / group
        if not folder.is_dir():
            continue
        require(not folder.is_symlink(), 'Friendly view cannot be a symlink')
        for current, directories, names in os.walk(folder, followlinks=False):
            directories[:] = [d for d in directories if d.lower() not in EXCLUDED_DIRS and not d.startswith('.')
                              and not (Path(current) / d).is_symlink()]
            for name in names:
                path = Path(current) / name
                lower = name.lower()
                if path.is_symlink() or not path.is_file():
                    continue
                if lower.startswith(('.', '~$')) or '.g.vcf' in lower or lower.endswith(('.part', '.crdownload', '.tmp', '.lock',
                    '.fastq', '.fastq.gz', '.fq', '.fq.gz', '.fasta', '.fasta.gz', '.fa', '.fa.gz',
                    '.db', '.sqlite', '.sqlite3', '.pem', '.key', '.p12', '.pfx')):
                    continue
                if any(word in lower for word in ('credential', 'private-key', 'secret', 'token')):
                    continue
                require(root in path.resolve().parents, 'File escaped results root')
                relative = path.relative_to(root).as_posix()
                stat = path.stat()
                if relative not in history:
                    if group in ('BAM', 'VCF'):
                        continue  # Never package a preallocated/unimported biological output.
                    require(time.time() - stat.st_mtime >= 60, 'Local report/config is still changing')
                else:
                    expected = history[relative]
                    require(stat.st_size == expected['ExpectedSizeBytes'] and stat.st_mtime_ns == expected['ExpectedMtimeNs'],
                            'Imported file changed size/mtime; refusing stale verification label')
                row = {'RelativePath': relative, 'SizeBytes': stat.st_size, 'MtimeNs': stat.st_mtime_ns,
                       'SHA256CalculatedByPacker': False,
                       **history.get(relative, {'VerificationAtImport': 'local_size_only', 'HistoricalSHA256': None})}
                rows.append(row)
    return sorted(rows, key=lambda x: x['RelativePath'])


def archive_index(rows):
    body = ['<!doctype html><html lang="ru"><meta charset="utf-8"><title>Результаты анализа</title>',
            '<style>body{font:16px system-ui;max-width:1200px;margin:36px auto;padding:0 20px}td,th{border-bottom:1px solid #ddd;padding:9px;text-align:left}table{width:100%;border-collapse:collapse}</style>',
            '<h1>Результаты анализа</h1><p>Файлы сгруппированы по папкам. Состав и метод проверки записаны в ARCHIVE_MANIFEST.json.</p>']
    for group in VIEWS:
        items = [r for r in rows if r['RelativePath'].startswith(group + '/')]
        body.append('<h2>' + group + '</h2><table><tr><th>Файл</th><th>GiB</th><th>Проверка при получении</th></tr>')
        for row in items:
            label = 'SHA256 ранее проверен' if row['VerificationAtImport'] == 'sha256' else 'Только размер'
            body.append('<tr><td><a href="' + quote(row['RelativePath'], safe='/') + '">' +
                        html.escape(row['RelativePath']) + '</a></td><td>' + f'{row["SizeBytes"]/2**30:.4f}' +
                        '</td><td>' + label + '</td></tr>')
        body.append('</table>')
    body.append('</html>')
    return '\n'.join(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', required=True)
    parser.add_argument('--output', required=True, help='New ZIP outside the friendly input views')
    parser.add_argument('--folder-name', default='analysis_results')
    parser.add_argument('--expected-bams', type=int)
    parser.add_argument('--cohort-accepted', action='store_true')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    root, output = local_path(args.input_root), local_path(args.output)
    require(root.is_dir(), 'Results root missing')
    require(output.suffix.lower() == '.zip', 'Output must be a new ZIP file')
    require(re.fullmatch(r'[A-Za-z0-9_.-]{1,100}', args.folder_name) and args.folder_name not in ('.', '..'), 'Safe ZIP root folder name required')
    require(not any((p / '.git').exists() for p in (output.parent, *output.parents)), 'Patient archive cannot be inside Git')
    for group in VIEWS:
        view = (root / group).resolve()
        require(output != view and view not in output.parents, 'Archive cannot be inside an input view')
    rows = collect(root)
    total = sum(r['SizeBytes'] for r in rows)
    bams = [r for r in rows if r['RelativePath'].startswith('BAM/') and r['RelativePath'].endswith('.bam')]
    require(rows, 'No completed friendly files to pack')
    estimated = total + 64 * 2**20 + len(rows) * 2048
    print(f'Plan: {len(rows)} files, {len(bams)} BAMs, {total/2**30:.3f} GiB; archive upper estimate {estimated/2**30:.3f} GiB; no file SHA hashing', flush=True)
    if not args.execute:
        print('DRY RUN: no ZIP created. Execute only after the accepted cohort VCF and all BAM/index pairs are present.', flush=True)
        return
    require(args.cohort_accepted, 'Cohort VCF must be accepted before building the one final archive')
    if args.expected_bams is not None:
        require(len(bams) == args.expected_bams, 'Expected BAM count is not complete')
    paths = {r['RelativePath'] for r in rows}
    for bam in bams:
        name = bam['RelativePath']
        require(name + '.bai' in paths or name[:-4] + '.bai' in paths, 'BAM index missing')
    cohorts = [r for r in rows if r['RelativePath'].startswith('VCF/') and Path(r['RelativePath']).name.startswith('cohort')
               and '.vcf' in r['RelativePath'] and not r['RelativePath'].endswith(('.tbi', '.csi'))]
    require(cohorts, 'Combined cohort VCF missing')
    for cohort in cohorts:
        name = cohort['RelativePath']
        if name.endswith('.vcf.gz'):
            require(name + '.tbi' in paths or name + '.csi' in paths, 'Combined VCF index missing')
    output.parent.mkdir(parents=True, exist_ok=True)
    require(not output.exists(), 'Refusing to overwrite an existing archive')
    temporary = output.with_name(output.name + '.partial')
    require(not temporary.exists(), 'Existing partial archive requires review; refusing overwrite')
    require(shutil.disk_usage(output.parent).free >= estimated + 15 * 2**30,
            'Insufficient space for archive upper estimate plus 15 GiB reserve')
    manifest = {'SchemaVersion': 1, 'CreatedUtc': dt.datetime.now(dt.timezone.utc).isoformat(),
                'FileCount': len(rows), 'TotalSourceBytes': total, 'CohortAccepted': True,
                'SHA256CalculatedByPacker': False, 'FullDatasetSHA256Pass': False,
                'Policy': 'Sizes/stability checked while packaging; historical import verification retained; ZIP64 with standard ZIP CRC',
                'Files': rows}
    with zipfile.ZipFile(temporary, 'x', compression=zipfile.ZIP_DEFLATED, allowZip64=True,
                         compresslevel=6, strict_timestamps=False) as archive:
        for group in VIEWS:
            archive.writestr(args.folder_name + '/' + group + '/', b'')
        archive.writestr(args.folder_name + '/INDEX.html', archive_index(rows))
        archive.writestr(args.folder_name + '/ARCHIVE_MANIFEST.json', json.dumps(manifest, ensure_ascii=False, indent=2))
        remaining = total
        for number, row in enumerate(rows, 1):
            require(shutil.disk_usage(output.parent).free >= remaining + 15 * 2**30,
                    'Disk reserve changed while packing; completed source files remain intact')
            source = root / row['RelativePath']
            before = source.stat()
            require((before.st_size, before.st_mtime_ns) == (row['SizeBytes'], row['MtimeNs']), 'Source changed before packaging')
            stored = source.name.lower().endswith(('.bam', '.cram', '.vcf.gz', '.bcf', '.gz', '.zip'))
            archive.write(source, args.folder_name + '/' + row['RelativePath'],
                          compress_type=zipfile.ZIP_STORED if stored else zipfile.ZIP_DEFLATED)
            after = source.stat()
            require((after.st_size, after.st_mtime_ns) == (before.st_size, before.st_mtime_ns), 'Source changed during packaging')
            remaining -= row['SizeBytes']
            print(f'Packed {number}/{len(rows)}', flush=True)
    os.link(temporary, output)  # Atomic no-replace publication; same output directory.
    temporary.unlink()
    print(f'PASS ZIP64 complete: {output.stat().st_size} bytes; source files preserved; no SHA hashing', flush=True)


if __name__ == '__main__':
    main()
