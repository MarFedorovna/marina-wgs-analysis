#!/usr/bin/env python3
"""Export completed research tables to Word, without inventing clinical conclusions."""
import argparse
import csv
import json
from pathlib import Path

from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.image.image import Image


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def table(document, headings, rows):
    if len(rows) > 10:
        for start in range(0, len(rows), 10):
            document.add_page_break()
            table(document, headings, rows[start:start + 10])
        return
    output = document.add_table(rows=1, cols=len(headings))
    output.autofit = True
    for cell, value in zip(output.rows[0].cells, headings):
        cell.text = str(value)
        shade = OxmlElement('w:shd'); shade.set(qn('w:fill'), '25476A'); cell._tc.get_or_add_tcPr().append(shade)
        for run in cell.paragraphs[0].runs:
            run.bold = True; run.font.color.rgb = RGBColor(255, 255, 255)
    repeat = OxmlElement('w:tblHeader'); output.rows[0]._tr.get_or_add_trPr().append(repeat)
    for index, values in enumerate(rows):
        cells = output.add_row().cells
        for cell, value in zip(cells, values):
            cell.text = str(value if value is not None else 'Не указано')
            if index % 2 == 0:
                shade = OxmlElement('w:shd'); shade.set(qn('w:fill'), 'F0F3F7'); cell._tc.get_or_add_tcPr().append(shade)
    for row in output.rows:
        for cell in row.cells:
            props = cell._tc.get_or_add_tcPr()
            margins = OxmlElement('w:tcMar')
            for side in ('top', 'left', 'bottom', 'right'):
                edge = OxmlElement('w:' + side); edge.set(qn('w:w'), '85'); edge.set(qn('w:type'), 'dxa'); margins.append(edge)
            props.append(margins)
            borders = OxmlElement('w:tcBorders')
            for side in ('top', 'left', 'bottom', 'right'):
                edge = OxmlElement('w:' + side); edge.set(qn('w:val'), 'single'); edge.set(qn('w:sz'), '4'); edge.set(qn('w:color'), 'D9D9D9'); borders.append(edge)
            props.append(borders)
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_after = Pt(3)
                for run in paragraph.runs:
                    run.font.size = Pt(9)
    document.add_paragraph()


def export(source, output, synthetic=False):
    source, output = Path(source).resolve(strict=True), Path(output).resolve()
    if output.exists():
        raise ValueError('Refusing to overwrite an existing Word report.')
    methods = load(source / 'methods_and_limitations.json')
    counts = load(source / 'consequence_counts.json')
    document = Document()
    section = document.sections[0]
    section.page_width, section.page_height = Inches(8.5), Inches(11)
    section.top_margin = section.bottom_margin = Inches(1)
    section.left_margin = section.right_margin = Inches(0.8)
    for name in ('Normal', 'Title', 'Heading 1', 'Heading 2'):
        document.styles[name].font.name = 'Arial'
        document.styles[name].font.color.rgb = RGBColor(0, 0, 0)
    document.styles['Normal'].font.size = Pt(11)
    document.styles['Normal'].paragraph_format.space_after = Pt(7)
    document.add_heading('Проверка шаблона отчёта' if synthetic else 'Исследовательский отчёт по вариантам', 0)
    if synthetic:
        document.add_paragraph('СИНТЕТИЧЕСКИЙ ПРИМЕР. Эти значения не являются результатами анализа пациентов.')
    document.add_paragraph('Автоматически собранные результаты для экспертного рассмотрения. Клиническая классификация и заключение требуют отдельной оценки доказательств.')
    document.add_heading('Область исследования', 1)
    scope = methods.get('analysis_scope', {})
    if isinstance(scope, dict):
        for key, label in [('default', 'Исследованные интервалы'), ('sex_chromosomes', 'X/Y'),
                           ('not_assessed', 'Не оценивалось'), ('coverage', 'Покрытие')]:
            if key in scope:
                value = ', '.join(map(str, scope[key])) if isinstance(scope[key], list) else str(scope[key])
                document.add_paragraph(label + ': ' + value)
    else:
        document.add_paragraph(str(scope))
    document.add_paragraph('Клинические сведения, родословная, фенотип и семейная сегрегация включаются только при наличии проверенных исходных данных. Отсутствующие сведения не восстанавливаются автоматически.')
    document.add_heading('Методы', 1)
    document.add_paragraph('Готовые BAM используются для поиска малых вариантов в заданных интервалах. После совместного генотипирования выполняются нормализация, маркировка качества сайта и генотипа, офлайн-аннотация и формирование очередей рассмотрения. Точные версии, команды и контрольные суммы сохраняются в квитанциях этапов и манифесте ресурсов.')
    document.add_paragraph('Правило приоритизации: ' + str(methods.get('ranking', {}).get('rule', 'не указано')) + '. Порядок критериев: ' + ', '.join(methods.get('ranking', {}).get('order', [])) + '.')
    document.add_heading('Сводка вариантов', 1)
    headings = ['Образец', 'Исходные', 'После QC', 'Кандидаты']
    summary = []
    for sample, levels in sorted(counts.get('samples', {}).items()):
        summary.append([sample, *[levels.get(level, {}).get('unique_normalized_alleles', 0) for level in ('raw', 'qc', 'candidate')]])
    table(document, headings, summary or [['Нет строк', 0, 0, 0]])
    document.add_paragraph('Числа относятся к уникальным нормализованным аллелям. Несколько транскриптов одной аллели не увеличивают этот счётчик. Полные определения сохранены в consequence_counts.json.')
    flags = ['stop_gained', 'stop_lost', 'start_lost', 'frameshift', 'splice_donor_acceptor', 'plof_candidate']
    document.add_heading('Функциональные последствия', 2)
    rows = [[flag, *[counts.get('cohort', {}).get(level, {}).get('nonexclusive_flags', {}).get(flag, 0)
                    for level in ('raw', 'qc', 'candidate')]] for flag in flags]
    table(document, ['Флаг', 'Исходные', 'После QC', 'Кандидаты'], rows)
    document.add_paragraph('Флаги пересекаются между собой; их сумма не равна общему числу вариантов. pLoF-кандидат не означает доказанный LoF, PVS1 или патогенность; stop-loss не означает автоматически GoF.')
    for image in sorted(source.glob('figures/cohort__*.png')):
        if image.is_symlink():
            continue
        dimensions = Image.from_file(str(image))
        width = min(6.6, 7.5 * dimensions.px_width / dimensions.px_height)
        document.add_picture(str(image), width=Inches(width))
        document.add_paragraph(image.name)
    document.add_heading('7.2 Приоритетные варианты', 1)
    top = source / 'cohort_top20.tsv'
    records = []
    if top.is_file():
        with top.open(encoding='utf-8', newline='') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                records.append([row.get('rank'), row.get('gene'), row.get('variant_id'), row.get('Consequence') or row.get('consequence', ''), row.get('carrier_samples', '')])
    if records:
        table(document, ['Ранг', 'Ген', 'Вариант', 'Последствие', 'Носители'], records[:20])
    else:
        document.add_paragraph('В очереди нет подходящих вариантов. Это не является отрицательным клиническим результатом.')
    document.add_paragraph('Топ-10 и топ-20 по каждому образцу, полные приложения, отдельные записи для срочной ручной проверки и лист критериев ACMG находятся в соседних таблицах. Очереди не дополняются искусственно до заданного числа строк.')
    document.add_heading('Доступность источников', 1)
    document.add_paragraph('Фактически загруженные версии и статусы всех дополнительных источников приведены в resource_availability.json. Упоминание базы в курсе не означает её установку или наличие лицензии. Недоступные оценки не превращаются в доброкачественные предсказания или нулевую частоту.')
    resources = methods.get('resources', {})
    installed = resources.get('resources', {})
    if installed:
        table(document, ['Источник', 'Версия', 'Статус'], [
            [name, item.get('version', ''), item.get('status', '')] for name, item in sorted(installed.items())])
    audited = resources.get('audited_source_coverage', {}).get('sources', [])
    if audited:
        document.add_heading('Покрытие требований по базам', 2)
        groups = {}
        for item in audited:
            groups.setdefault(item['status'], []).append(item['resource_id'])
        for status, names in sorted(groups.items()):
            document.add_paragraph(status + ': ' + ', '.join(names) + '.')
        document.add_paragraph('Точные версии, ссылки, лицензии и ограничения каждого источника сохранены в resource_availability.json; все источники из реестра перечислены выше по фактическому статусу.')
    document.add_paragraph('Фенотип: ' + str(methods.get('phenotype_status', 'не указан')) + '. Путевые анализы: ' + str(methods.get('pathway_status', 'не выполнялись')) + '.')
    document.add_heading('Ограничения и дальнейшая проверка', 1)
    for item in methods.get('limitations', []):
        document.add_paragraph(item, style='List Bullet')
    document.add_heading('Экспертное заключение', 1)
    document.add_paragraph('Не сформировано автоматически. Необходимо проверить покрытие целевых областей, качество выбранных генотипов, соответствие фенотипу и наследованию, первичные утверждения баз и литературу; применимость критериев ACMG/ClinGen документируется отдельно.')
    document.add_heading('Файлы приложения', 1)
    for path in sorted(source.iterdir()):
        if path.is_file() and not path.name.startswith('.'):
            document.add_paragraph(path.name, style='List Bullet')
    footer = section.footer.paragraphs[0]
    footer.text = 'Исследовательская обработка | Не является автоматическим диагнозом'
    footer.style = document.styles['Normal']; footer.runs[0].font.size = Pt(8)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('xb') as stream:
        document.save(stream)
    print(json.dumps({'status': 'PASS', 'output': str(output), 'synthetic': synthetic, 'automatic_diagnosis': False}, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--synthetic', action='store_true')
    args = parser.parse_args()
    export(args.source, args.output, args.synthetic)


if __name__ == '__main__':
    main()
