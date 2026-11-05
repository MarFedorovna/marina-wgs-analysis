#!/usr/bin/env python3
"""Write a timestamped public-resource status snapshot without patient inputs."""
import argparse
import datetime as dt
import html
import json
from pathlib import Path


def load(path, default=None):
    return json.loads(path.read_text()) if path.is_file() else (default or {})


def snapshot(root, output):
    now = dt.datetime.now(dt.timezone.utc).isoformat()
    state = load(root / 'private/databases-supervisor-status.json')
    queue = load(root / 'private/annotation-queue-status.json')
    partial = root / 'resources/.downloads/homo_sapiens_vep_115_GRCh38.tar.gz.agent-partial'
    parts = load(partial.with_name(partial.name.replace('.agent-partial', '.agent-parts.json')))
    identity = parts.get('identity', {})
    total, chunk = identity.get('bytes', 25272957721), identity.get('chunk_bytes', 67108864)
    verified = sum(min(chunk, total - int(index) * chunk) for index in parts.get('completed', {}))
    allocated = partial.stat().st_blocks * 512 if partial.is_file() else 0
    rows = []
    labels = {'clinvar': 'ClinVar GRCh38', 'clinvar_index': 'Индекс ClinVar',
        'clingen_validity': 'ClinGen: связь ген–заболевание', 'clingen_dosage': 'ClinGen: dosage sensitivity',
        'hpo_gene': 'HPO: гены и фенотипы', 'hpo_ontology': 'HPO: онтология',
        'hgnc': 'HGNC: официальные символы генов', 'mane': 'MANE: метаданные транскриптов'}
    for receipt_path in sorted((root / 'resources/annotation/.public-source-cache').glob('*/resource.json')):
        item = load(receipt_path)
        source = receipt_path.parent / item['path']
        exists = source.is_file() and source.stat().st_size == item['bytes']
        rows.append([labels.get(item['name'], item['name']), item.get('version'),
            'Готов по проверенной квитанции' if item.get('status') == 'available' and exists else 'Требует проверки',
            f"{item['bytes']:,} байт; SHA256 записан при подготовке"])
    alpha = load(root / 'resources/annotation-supplements/alphamissense-2023-v2-grch38-v3/resource-manifest.json')
    alpha_control = load(root / 'private/alphamissense-control-v1/public-controls-receipt.json')
    alpha_status = ('Готов; отдельные 20 публичных контролей PASS' if alpha_control.get('status') == 'PASS'
                    else 'Готов по квитанции ресурса; отдельный контроль здесь не подтверждён')
    rows.append(['AlphaMissense canonical + isoforms', alpha.get('dataset_version', 'не установлен'),
        alpha_status if alpha.get('status') == 'PASS' else 'Не готов',
        'GRCh38.p13 / GENCODE v32; точное сопоставление transcript.version и белковой замены'])
    esc = lambda value: html.escape(str(value))
    table = '<table><tr><th>Ресурс</th><th>Версия</th><th>Состояние</th><th>Объём / ограничение</th></tr>'
    table += ''.join('<tr>' + ''.join('<td>' + esc(value) + '</td>' for value in row) + '</tr>' for row in rows) + '</table>'
    missing = ['FATHMM — публичный legacy SQL/HMM backend пока не подключён',
        'REVEL, Varity, MetaLR, MaxEntScan — отдельные оценки пока не загружены и не вычислены',
        'SpliceAI — не установлен; нужны разрешённый локальный ресурс и проверка лицензии',
        'CADD — не установлен; название четвёртого предиктора в видео не подтверждено',
        'HGMD Professional, OMIM, dbNSFP — локальный разрешённый доступ/лицензии и adapters не предоставлены',
        'GeneReviews, LOVD, Franklin, VarSome — ручная проверка; пациентские API-запросы отключены']
    body = f'''<!doctype html><html lang="ru"><meta charset="utf-8"><title>Базы и аннотация: состояние</title>
<style>body{{font:16px system-ui;max-width:1200px;margin:32px auto;padding:0 20px;color:#17242d}}h1,h2{{color:#16475b}}table{{border-collapse:collapse;width:100%;font-size:14px}}td,th{{border:1px solid #ccd6dc;padding:10px;text-align:left;vertical-align:top}}th{{background:#eaf3f7}}.box{{background:#eef6fa;padding:18px;border-left:5px solid #2882a1}}small{{color:#52616c}}li{{margin:9px 0}}code{{background:#edf0f2;padding:2px 5px}}</style>
<h1>Базы данных и аннотация: фактическое состояние</h1>
<p><b>Снимок на {esc(now)} UTC.</b> Это статический снимок, без автоматического обновления.</p>
<div class="box"><b>VEP 115.2 как программа установлен.</b> Большой GRCh38-кэш VEP 115 ещё подготавливается.
Проверено блоками: <b>{verified:,} / {total:,} байт ({verified/total*100:.1f}%)</b>.
Выделено в частичном файле: {allocated:,} байт; незавершённые блоки ещё не имеют подтверждённого SHA.</div>
<p>Супервизор: <code>{esc(state.get('status', 'UNKNOWN'))}</code>; PID загрузчика: {esc(state.get('current_pid'))}.
Политика: {esc(state.get('retry_policy', 'fixed attempt limit'))}. Исходный бюджетный дедлайн: {esc(state.get('deadline_utc', 'unknown'))}.</p>
<h2>Уже подготовлено на сервере</h2>{table}
<p><small>«Готов» здесь относится к квитанции подготовки ресурса и наличию файла, а не к завершённой аннотации наших вариантов. Этот снимок не пересчитывает SHA всех гигабайтов повторно.</small></p>
<h2>Что даст ожидаемый кэш VEP</h2>
<p>Рецепт проверяет SIFT 6.2.1, PolyPhen-2 2.2.3, gnomAD 4.1 AF и MANE 1.4 по реальному cache info.txt.
Эти компоненты <b>ещё не считаются готовыми</b>, пока весь архив не получен, не проверен и не установлен.
Кэш частот не заменяет полный gnomAD с AN, покрытием, constraint и всеми клиническими доказательствами.</p>
<h2>Почему загрузка занимает время</h2>
<p>Архив — 25 272 957 721 байт (около 23.54 GiB). Прямая сеть сервера до официального Ensembl даёт переменную скорость и отказы <code>Connection refused</code>.
Предыдущие две попытки сохранили 8 GiB проверенных блоков, затем старый лимит двух попыток остановил процесс.
Теперь временные сетевые ошибки автоматически повторяются с паузой 60–900 секунд до прежнего абсолютного бюджетного дедлайна.
При несовпадении данных, версии или контрольной суммы повтор запрещён. Срок окончания не гарантируется; маршрут через компьютер не используется.</p>
<h2>Дальнейшая очередь</h2><p>Состояние: <code>{esc(queue.get('status', 'UNKNOWN'))}</code>.
После полного PASS ресурсов: двадцать публичных контрольных вариантов через VEP, затем этап 6 аннотации и этап 7 таблиц/HTML/Word.
Это исследовательский отчёт, без автоматического диагноза или клинического класса ACMG.</p>
<h2>Отдельно запланировано / пока не установлено</h2><ul>{''.join('<li>'+esc(x)+'</li>' for x in missing)}</ul>
<p>Упоминание инструмента в лекции не означает его установку. Недоступные оценки остаются missing и не превращаются в «доброкачественно» или нулевую частоту.</p>
<p>Логи сервера: <code>/srv/marina/logs/annotation-databases.log</code>, <code>/srv/marina/logs/annotation-supervisor.log</code>.</p></html>'''
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(body, encoding='utf-8')
    print(json.dumps({'status': 'SNAPSHOT_WRITTEN', 'path': str(output), 'utc': now,
        'verified_bytes': verified, 'total_bytes': total, 'public_small_receipts': len(rows) - 1}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    snapshot(args.root, args.output)
