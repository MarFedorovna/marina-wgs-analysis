#!/usr/bin/env python3
"""Write a local Russian START.html from result metadata only; keep INDEX.html."""
import argparse
from collections import defaultdict
import datetime as dt
import html
import json
import os
from pathlib import Path
import re
import uuid
from urllib.parse import quote


def read_json(path):
    try:
        if path.is_file() and not path.is_symlink() and path.stat().st_size<=8*1024**2:
            return json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError,ValueError):pass
    return {}


def size(value):
    for unit in ('Б','KiB','MiB','GiB','TiB'):
        if value<1024 or unit=='TiB':return f'{value:.1f} {unit}'
        value/=1024


def write_home(results,transfer_status=None):
    root=Path(results).resolve(strict=True)
    if not root.is_dir() or root.is_symlink() or (root/'.git').exists():raise ValueError('Select an existing local results directory outside Git')
    target=root/'START.html'
    if target.is_symlink():raise ValueError('Refusing a symlink START.html')
    e=html.escape
    def href(path):return quote(os.path.relpath(path,root).replace('\\','/'),safe='/')
    def final(path,expected=None):
        return path.is_file() and not path.is_symlink() and '.part' not in path.name and (expected is None or path.stat().st_size==expected.get('bytes',-1))
    def link(path,label=None,expected=None):
        label=label or path.name
        if final(path,expected):return f'<a href="{href(path)}">{e(label)}</a> <small>{size(path.stat().st_size)}</small>'
        return f'<span class="pending">{e(label)} — ожидается</span>'
    def card(title,body):return '<section class="card"><h2>'+e(title)+'</h2>'+body+'</section>'
    reports=root/'Reports';views=reports/'cohort-sample-views-v1';manifest=read_json(views/'manifest.json');artifacts=manifest.get('artifacts',{})
    ready=sum(final(views/name,value) for name,value in artifacts.items() if isinstance(value,dict))
    status=read_json(Path(transfer_status)) if transfer_status else {}
    progress=''
    if isinstance(status.get('committed_selected_bytes'),int) and isinstance(status.get('selected_bytes'),int) and status['selected_bytes']>0:
        got,total=status['committed_selected_bytes'],status['selected_bytes'];pct=min(100,100*got/total)
        progress=f'<p><strong>{pct:.1f}% подтверждено</strong> · {size(got)} из {size(total)}</p><progress value="{got}" max="{total}"></progress><p class="note">По карте завершённых блоков, включая проверенные локальные ссылки. Полный размер выделенного .part не является прогрессом.</p>'
        if status.get('status')=='PASS_ALL_CLOSED_DATASET_FILES_SIZE_ONLY':progress+='<p class="ok">Весь закрытый набор сохранён — SIZE ONLY.</p>'
        else:progress+='<p>Перенос ещё не завершён. Статус: '+e(str(status.get('status','неизвестен')))+'</p>'
    cards=[]
    main=root/'VCF/cohort.fully-annotated.vcf.gz';main_note=read_json(main.with_name(main.name+'.completed.json'))
    checked=final(main) and main_note.get('Status')=='COMPLETED_SIZE_ONLY' and main_note.get('LocalMtimeNs')==main.stat().st_mtime_ns and main_note.get('SourceIdentity',{}).get('ExpectedBytes')==main.stat().st_size
    cards.append(card('1 · Основной аннотированный VCF','<p>'+link(main,'Общий VCF')+'<br>'+link(Path(str(main)+'.tbi'),'Индекс TBI')+'</p><p>'+('Локальная копия завершена, проверка SIZE ONLY.' if checked else 'Статус локальной копии сверяйте с квитанцией.')+'</p><p class="note">Здесь полные транскриптные строки, применённые источники и SpliceAI. Это ядерный набор; MT и X/Y — отдельные выдачи.</p>'))
    rows=[]
    for i,bam in enumerate(sorted((root/'BAM').rglob('*.bam')),1):
        if final(bam):rows.append('<li>'+link(bam,f'BAM {i}')+' · '+link(Path(str(bam)+'.bai'),'BAI')+f' · <a href="{href(bam.parent)}/">Папка</a></li>')
    cards.append(card('2 · BAM и индексы','<ul>'+''.join(rows)+'</ul>' if rows else '<p class="pending">BAM пока не найдены.</p>'))
    summary=f'<p>Готовые артефакты на диске: <strong>{ready} / {len(artifacts)}</strong>. .part не учитываются.</p>' if artifacts else '<p>Манифест таблиц пока не сохранён.</p>'
    summary+=progress+'<p>'+link(views/'index.html','Страница полной выдачи')+'</p><ul>'
    for name,label in [('masters/cohort.master.tsv.gz','Общий master — все транскрипты'),('tables/cohort.top20.tsv','Общий top-20'),('tables/shared_ge2.top20.tsv','Общие ≥2: top-20'),('tables/all4.top20.tsv','Общие всем четырём: top-20'),('tables/allele-membership-matrix.tsv.gz','Матрица точных аллелей / доля 25–100%'),('tables/spliceai/cohort.top20.tsv','SpliceAI top-20 — рейтинг предсказаний'),('tables/ACMG-manual-top20-worksheet.tsv.gz','ACMG: лист ручной проверки')]:
        summary+='<li>'+link(views/name,label,artifacts.get(name))+'</li>'
    summary+='</ul><p class="note">Top-20 — исследовательский рейтинг. Соответствие фенотипу, механизм и ACMG проверяются отдельно; это не автоматический диагноз. Общность означает тот же точный аллель, а не только общий ген.</p>'
    cards.append(card('3 · Общий набор и приоритизация',summary))
    sample_rows=[]
    labels=manifest.get('sample_file_labels',{})
    for i,(sample,label) in enumerate(labels.items() if isinstance(labels,dict) else [],1):
        if not isinstance(label,str) or not re.fullmatch(r'[A-Za-z0-9_.-]+',label):continue
        body=f'<h3>Образец {i} · {e(str(sample))}</h3><ul>'
        for name,title in [(f'vcf/{label}.annotated.vcf.gz','VCF только этого образца'),(f'vcf/{label}.annotated.vcf.gz.tbi','TBI'),(f'masters/{label}.master.tsv.gz','Полный master'),(f'tables/{label}.top20.tsv','Top-20'),(f'tables/{label}.functional.tsv.gz','Функциональные последствия')]:body+='<li>'+link(views/name,title,artifacts.get(name))+'</li>'
        sample_rows.append(body+'</ul>')
    cards.append(card('4 · Четыре индивидуальных набора',''.join(sample_rows) or '<p class="pending">Метки образцов появятся после манифеста.</p>'))
    categories=defaultdict(list)
    for name,expected in artifacts.items():
        parts=Path(name).parts
        if len(parts)==3 and parts[0]=='tables' and name.endswith('.tsv.gz') and parts[1]!='spliceai':categories[parts[1]].append((name,expected))
    body='<p>'+link(views/'tables/functional-category-index.tsv','Указатель всех функциональных категорий')+'</p>'
    for scope,items in sorted(categories.items()):
        done=sum(final(views/name,expected) for name,expected in items);scope_label={'cohort':'Весь набор','shared_ge2':'Общие ≥2','all4':'Все четыре'}.get(scope,scope)
        body+=f'<details><summary>{e(scope_label)} · {done}/{len(items)} файлов</summary><ul>'
        for name,expected in sorted(items):
            label=re.sub(r'-[0-9a-f]{8}\.tsv\.gz$','',Path(name).name);body+='<li>'+link(views/name,label,expected)+'</li>'
        body+='</ul></details>'
    body+='<p class="note">Категории могут пересекаться. Stop gained/lost, frameshift и splice — последствия, а не доказанный механизм LoF/GoF или клиническая тяжесть.</p>'
    cards.append(card('5 · Функциональные типы',body))
    extra=[]
    for folder in sorted(reports.iterdir() if reports.is_dir() else []):
        if not folder.is_dir() or folder.is_symlink() or folder==views:continue
        name=folder.name.lower()
        if not any(token in name for token in ('xy','mitochond','completed-mt')):continue
        receipt=read_json(folder/'LOCAL-COPY-RECEIPT.json');stage=read_json(folder/'manifest.json');wf=read_json(folder/'workflow-receipt.json')
        label='X/Y calling' if name=='xy-germline-v1' else 'MT: последствия' if 'mitochondrial-consequence' in name else 'MT: функциональные типы' if 'mitochond' in name and 'functional' in name else 'X/Y: аннотация' if 'xy' in name and ('annotation' in name or 'consequence' in name) else 'MT и X/Y: готовые QC и частотный контекст' if 'completed-mt' in name else folder.name
        local=receipt.get('status','Локальная папка; состав проверяйте по квитанции')
        extra.append(f'<h3>{e(label)}</h3><p><a href="{href(folder)}/">Открыть папку</a> · '+e(str(local))+'</p>')
        indexes=sorted(folder.glob('**/index.html'))
        if indexes:extra.append('<p>'+link(indexes[0],'Страница набора')+'</p>')
        vcf=[x for x in folder.rglob('*.vcf.gz') if final(x) and ('fully' in x.name or x.name=='cohort.xy.soft-filtered.vcf.gz')]
        if vcf:extra.append('<ul>'+''.join('<li>'+link(x,x.name)+'</li>' for x in vcf)+'</ul>')
    if (root/'QC').is_dir():extra.append('<p><a href="'+href(root/'QC')+'/">Все технические QC</a></p>')
    extra.append('<p class="note">Calling находит варианты; аннотация добавляет последствия и источники. X/Y-модель обработки по QC не является клиническим определением пола или кариотипа. MT сохраняет свою политику гетероплазмии и QC.</p>')
    cards.append(card('6 · X/Y, MT и технический QC',''.join(extra) or '<p class="pending">Дополнительные наборы пока не сохранены.</p>'))
    archives=sorted((root.parent/'code').glob('*.zip'),key=lambda p:p.stat().st_mtime,reverse=True)
    tail='<p>'+link(archives[0],'Последний архив общего кода')+'</p>' if archives else '<p>Архив кода пока не найден рядом с results.</p>'
    tail+='<p>'+link(root/'INDEX.html','Подробный прежний INDEX.html')+'</p><p class="note">Исходные файлы, старый INDEX и результаты не изменялись. На этой странице только ссылки и метаданные, без повторной аннотации, хеширования больших файлов или сетевых запросов.</p>'
    cards.append(card('7 · Код и подробный инвентарь',tail))
    phenotype=reports/'phenotype-hcm-rerank-v1'
    if read_json(phenotype/'manifest.json').get('status')=='PASS_FULL_TRANSCRIPT_HCM_RESEARCH_RERANK':
        workbook=phenotype/'HCM-review-with-read-support.xlsx'
        if not workbook.is_file():workbook=phenotype/'HCM-review.xlsx'
        body='<p>'+link(workbook,'Excel: общий и индивидуальные HCM top-20, MT и доступная проверка чтений')+'</p>'
        support=reports/'hcm-read-support-review-v1'
        if read_json(support/'status.json').get('status')=='PASS_LOCAL_INDEXED_SNV_READ_SUPPORT_REVIEW':
            body+='<p>'+link(support/'cohort.read-support.tsv','Поддержка чтениями: все четыре образца')+' · '+link(support/'targets.igv.bed','Позиции для ручного IGV-просмотра')+'</p>'
        body+='<p>'+link(phenotype/'index.html','Страница фенотипного отчёта')+' · '+link(phenotype/'functional-category-index.tsv','Ранжированные функциональные категории кандидатов')+'</p>'
        body+='<p><a href="'+href(phenotype/'tables')+'/">Полные очереди кандидатов и HCM top-20</a></p>'
        body+='<p class="note">ГКМП-контекст HPO/ClinGen — очередь для разбора, а не доказанный причинный вариант. Пустой общий HCM top означает отсутствие совпадений в выбранной очереди. ACMG и наследование требуют отдельной проверки.</p>'
        cards.append(card('8 · Фенотипное ранжирование ГКМП',body))
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    document='''<!doctype html><html lang="ru"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Результаты анализа</title><style>body{font:16px/1.5 system-ui,sans-serif;background:#f0f4f8;color:#172b3a;margin:0}header,main,footer{max-width:1120px;margin:auto;padding:24px}header{padding-bottom:4px}h1{font-size:30px;margin:0}h2{font-size:21px;margin-top:0}h3{font-size:17px}main{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:18px}.card{background:white;border:1px solid #d9e3eb;border-radius:14px;padding:22px;overflow-wrap:anywhere}a{color:#056998}li{margin:7px 0}small,.note{color:#566c7e;font-size:13px}.pending{color:#91622d}.ok{color:#197043}progress{width:100%;height:15px}details{border-top:1px solid #e5edf2;padding:10px 0}summary{cursor:pointer;font-weight:600}ul{padding-left:20px}</style><header><h1>Результаты анализа</h1><p>Общий аннотированный VCF → таблицы и четыре индивидуальных набора. X/Y и MT показаны отдельно.</p><p class="note">Обновлено '''+stamp+''' · SIZE ONLY означает проверку размера и неизменности источника, а не новый локальный SHA.</p></header><main>'''+''.join(cards)+'</main><footer>Для свежего прогресса заново сформируйте START.html после обновления файлов.</footer></html>'
    temp=target.with_name('.START-'+uuid.uuid4().hex+'.tmp')
    with temp.open('x',encoding='utf-8') as f:f.write(document);f.flush();os.fsync(f.fileno())
    os.replace(temp,target)
    return dict(status='WROTE_START_HTML_ONLY',html_bytes=target.stat().st_size,manifest_artifacts=len(artifacts),completed_artifacts=ready,linked_sample_sets=len(sample_rows),functional_scopes=len(categories),legacy_index_untouched=True,patient_payloads_read=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True);p.add_argument('--transfer-status',type=Path);a=p.parse_args();print(json.dumps(write_home(a.results,a.transfer_status)))
