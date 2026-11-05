#!/usr/bin/env python3
"""Index completed local files, explicitly separating SHA256 and size-only checks."""
import argparse
import datetime as dt
import html
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
from urllib.parse import quote


def atomic_json(path, value):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    os.replace(temp, path)


def category_path(relative):
    parts = PurePosixPath(relative).parts
    if not parts or '..' in parts or PurePosixPath(relative).is_absolute():
        raise ValueError('Unsafe verified relative path')
    name = parts[-1]
    if name.endswith(('.bam', '.bam.bai', '.bai')):
        sample = parts[-3] if len(parts) >= 3 and parts[-2] == 'align' else parts[-2]
        return 'BAM', Path('BAM', sample, name)
    if '.vcf' in name or name.endswith(('.bcf', '.tbi', '.csi')):
        if 'stages' in parts and name.startswith('cohort.'):
            return 'VCF', Path('VCF', name)
        return 'VCF', Path('VCF', *parts[1:] if parts[0] == 'outputs' else parts)
    if parts[:2] == ('metadata', 'code'):
        return 'Code', Path('Code', *parts[2:])
    if parts[:2] == ('metadata', 'private'):
        return 'Settings', Path('Settings', *parts[2:])
    if parts[:2] == ('metadata', 'logs') or name.endswith('.log'):
        return 'Logs', Path('Logs', *parts)
    if name.endswith(('.html', '.pdf', '.docx', '.xlsx', '.tsv', '.csv', '.tsv.gz', '.csv.gz', '.png', '.jpg', '.jpeg', '.svg')) or 'report' in name.lower():
        return 'Reports', Path('Reports', *parts)
    return 'QC and provenance', Path('QC', *parts)


def publish_report_alias(destination, entry, stored):
    """Keep offline report-relative links usable without duplicating file bytes."""
    relative = PurePosixPath(entry['RelativePath'])
    parts = relative.parts
    if relative.is_absolute() or '..' in parts:
        raise ValueError('Unsafe report-relative path')
    bundles = {'consequence-tables-v1', 'tissue-protein-context-v1',
               'population-reference-context-v1', 'cohort-sample-views-v1',
               'combined-annotation-master-v1'}
    allowed = len(parts) > 2 and parts[0] == 'outputs' and parts[1] in bundles
    allowed = allowed or (len(parts) > 4 and parts[:3] == ('outputs', 'workflow-v1', 'stages')
                          and parts[3] in ('06_annotate', '07_report'))
    if not allowed:
        return None
    alias = destination / 'Reports' / Path(*parts)
    if alias.is_symlink() or not alias.resolve().is_relative_to(destination.resolve()):
        raise ValueError('Report bundle alias escapes local result storage')
    alias.parent.mkdir(parents=True, exist_ok=True)
    if alias.exists():
        if not os.path.samefile(alias, stored):
            raise ValueError('Different report bundle file exists; refusing overwrite')
    else:
        os.link(stored, alias)
    return alias.relative_to(destination).as_posix()


def publish_verified(destination, entry, mapping):
    """Caller invokes only after its successful checksum check against receipt."""
    sha, rel = entry['Sha256'], entry['RelativePath']
    if not re.fullmatch('[0-9a-f]{64}', sha):
        raise ValueError('Invalid object hash')
    obj = destination / 'objects' / sha
    if not obj.is_file() or obj.is_symlink() or obj.stat().st_size != entry['SizeBytes']:
        raise ValueError('Verified object missing')
    group, friendly = category_path(rel)
    target = destination / friendly
    if target.exists() and not os.path.samefile(target, obj):
        friendly = Path('Versions', sha, friendly)
        target = destination / friendly
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if not os.path.samefile(target, obj):
            raise ValueError('Differing existing friendly file; refusing overwrite')
    else:
        os.link(obj, target)
    bundle_alias = publish_report_alias(destination, entry, obj)
    mapping[sha + ':' + rel] = {'Category': group, 'RelativePath': rel, 'LocalPath': friendly.as_posix(),
                                'Sha256': sha, 'SizeBytes': entry['SizeBytes'], 'ReportBundleAlias': bundle_alias}
    atomic_json(destination / 'verified-files.json', mapping)


def size_key(entry):
    # Hashes metadata only, never file bytes. It is deliberately outside objects/.
    identity = json.dumps([entry['RelativePath'], entry['SizeBytes'], entry['Sha256']], separators=(',', ':'))
    return 'source-' + hashlib.sha256(identity.encode()).hexdigest()[:32]


def size_record_valid(destination, entry):
    path = destination / 'size-only-files.json'
    if not path.is_file():
        return False
    records = json.loads(path.read_text(encoding='utf-8-sig'))
    record = records.get(size_key(entry))
    if not record or record['ExpectedServerSha256'] != entry['Sha256']:
        return False
    stored, friendly = destination / record['StoragePath'], destination / record['LocalPath']
    return stored.is_file() and not stored.is_symlink() and friendly.is_file() and not friendly.is_symlink() and \
        stored.stat().st_size == entry['SizeBytes'] and stored.stat().st_mtime_ns == record['MtimeNs'] and os.path.samefile(stored, friendly)


def publish_size_only(destination, entry, stored, completion, records=None, persist=True):
    key = size_key(entry)
    group, friendly = category_path(entry['RelativePath'])
    target = destination / friendly
    if target.is_symlink():
        raise ValueError('Friendly file cannot be a symlink')
    if target.exists() and not os.path.samefile(target, stored):
        friendly = Path('Versions', 'size-only', key, friendly)
        target = destination / friendly
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        if not os.path.samefile(target, stored):
            raise ValueError('Differing friendly file exists; refusing overwrite')
    else:
        os.link(stored, target)
    bundle_alias = publish_report_alias(destination, entry, stored)
    path = destination / 'size-only-files.json'
    if records is None:
        records = json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else {}
    record = {'Verification': 'size_only', 'SourceKey': key, 'Category': group,
              'RelativePath': entry['RelativePath'], 'LocalPath': friendly.as_posix(),
              'StoragePath': stored.relative_to(destination).as_posix(), 'SizeBytes': entry['SizeBytes'],
              'ActualSizeBytes': stored.stat().st_size, 'MtimeNs': stored.stat().st_mtime_ns,
              'ExpectedServerSha256': entry['Sha256'], 'ActualSha256': None, 'CompletionEvidence': completion,
              'ReportBundleAlias': bundle_alias}
    records[key] = record
    if persist:
        atomic_json(path, records)
    atomic_json(stored.parent / 'expected-metadata.json', record)
    return record


def write_index(destination):
    map_path, state_path = destination / 'verified-files.json', destination / 'verified-objects.json'
    mapping = json.loads(map_path.read_text(encoding='utf-8')) if map_path.exists() else {}
    state = json.loads(state_path.read_text(encoding='utf-8')) if state_path.exists() else {}
    groups, total = {}, 0
    for entry in mapping.values():
        obj = destination / 'objects' / entry['Sha256']
        target = destination / entry['LocalPath']
        if not obj.is_file() or not target.is_file() or obj.is_symlink() or target.is_symlink():
            continue
        identity = {'size': obj.stat().st_size, 'mtime_ns': obj.stat().st_mtime_ns}
        if state.get(entry['Sha256']) != identity or not os.path.samefile(obj, target):
            continue  # Never advertise partial or changed files as verified.
        groups.setdefault(entry['Category'], []).append(dict(entry, Verification='sha256', IntegrityLabel='SHA256 проверен'))
        total += entry['SizeBytes']
    size_path = destination / 'size-only-files.json'
    if size_path.is_file():
        for entry in json.loads(size_path.read_text(encoding='utf-8-sig')).values():
            stored, target = destination / entry['StoragePath'], destination / entry['LocalPath']
            if not stored.is_file() or not target.is_file() or stored.is_symlink() or target.is_symlink():
                continue
            if stored.stat().st_size != entry['SizeBytes'] or stored.stat().st_mtime_ns != entry['MtimeNs'] or not os.path.samefile(stored, target):
                continue
            groups.setdefault(entry['Category'], []).append(dict(entry, IntegrityLabel='Только размер; SHA256 не проверен'))
    paused = (destination / 'Settings' / 'LOCAL_DOWNLOADS_PAUSED.json').is_file()
    activity = 'Локальная загрузка на паузе. Возобновление — только по команде пользователя.' if paused else 'Остальные файлы могут ещё ожидать получения.'
    body = ['<!doctype html><html lang="ru"><meta charset="utf-8"><title>Локальная копия результатов</title>',
            '<style>body{font:16px system-ui;max-width:1200px;margin:36px auto;padding:0 20px}table{border-collapse:collapse;width:100%}td,th{border-bottom:1px solid #ddd;padding:10px;text-align:left}small{color:#555}code{font-size:12px}</style>',
            '<h1>Локальная копия результатов</h1>',
            '<p>' + activity + '</p>',
            '<p>Показаны завершённые локальные файлы; способ проверки указан в каждой строке. '
            'BAM — выравнивания; VCF/gVCF — варианты; QC — метрики качества. '
            'Техническое завершение анализа не является клиническим заключением.</p>',
            '<p>Снапшоты и их manifest/receipt находятся в папке <a href="snapshots/">snapshots</a>. '
            'Одинаковые файлы связаны hard links и занимают место один раз. Сохраняйте весь каталог резервной копии.</p>',
            '<small>Обновлено: ' + html.escape(dt.datetime.now(dt.timezone.utc).isoformat()) + '</small>']
    for group, entries in sorted(groups.items()):
        body.append('<h2>' + html.escape(group) + '</h2><table><tr><th>Файл</th><th>GiB</th><th>Проверка</th><th>SHA256 файла</th></tr>')
        for entry in sorted(entries, key=lambda x: (x['LocalPath'], x.get('Sha256', ''))):
            body.append('<tr><td><a href="' + quote(entry['LocalPath'], safe='/') + '">' +
                        html.escape(entry['LocalPath']) + '</a></td><td>' +
                        f'{entry["SizeBytes"] / 2**30:.4f}' + '</td><td>' + html.escape(entry['IntegrityLabel']) +
                        '</td><td><code>' + (entry['Sha256'] if entry['Verification'] == 'sha256' else 'не проверялся') + '</code></td></tr>')
        body.append('</table>')
    body.append('</html>')
    temp = destination / 'INDEX.html.tmp'
    temp.write_text('\n'.join(body), encoding='utf-8')
    os.replace(temp, destination / 'INDEX.html')
    return sum(len(v) for v in groups.values())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', required=True)
    args = parser.parse_args()
    value = args.destination
    if re.match(r'^[A-Za-z]:[\\/]', value):
        value = subprocess.check_output(['wslpath', '-u', value], text=True).strip()
    print('Indexed completed files:', write_index(Path(value).resolve(strict=True)))
