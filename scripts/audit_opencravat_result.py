#!/usr/bin/env python3
"""Read-only coverage audit of an OpenCRAVAT result SQLite database.

Column prefixes identify modules. A SQL NULL never proves benign status,
out-of-scope applicability, or a failed annotator; these need logs/source context.
https://docs.opencravat.org/en/latest/Annotator-Tutorial.html
"""
import argparse
import csv
import json
import sqlite3
from pathlib import Path

EXPECTED = 'alphamissense cadd fathmm revel sift polyphen2 spliceai gnomad4 clingen clinvar clinvar_acmg dbsnp hpo ncbigene omim loftool phi prec uniprot pubmed gtex'.split()


def q(name):
    return '"' + name.replace('"', '""') + '"'


def filled(col):
    # Zero is a valid numeric score, not missing.
    return f'({q(col)} IS NOT NULL AND TRIM(CAST({q(col)} AS TEXT)) NOT IN (\'\',\'.\',\'NA\',\'N/A\',\'null\',\'None\',\'[]\',\'{{}}\'))'


def audit(db, output):
    output.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db.resolve().as_uri() + '?mode=ro', uri=True)
    conn.execute('PRAGMA query_only=ON')
    conn.execute('PRAGMA cache_size=-8192')
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if 'variant' not in tables:
        raise ValueError('Expected OpenCRAVAT variant table is absent')
    metadata = {}
    for table in ('info', 'variant_annotator', 'gene_annotator'):
        if table in tables:
            cur = conn.execute('SELECT * FROM ' + q(table))
            metadata[table] = [dict(zip([c[0] for c in cur.description], r)) for r in cur]
    summaries, columns = [], []
    found = set()
    for table in ('variant', 'gene', 'sample', 'mapping'):
        if table not in tables:
            continue
        total = conn.execute('SELECT COUNT(*) FROM ' + q(table)).fetchone()[0]
        names = [r[1] for r in conn.execute('PRAGMA table_info(' + q(table) + ')')]
        module_columns = {m: [c for c in names if c.startswith(m+'__')] for m in EXPECTED}
        evidence_columns = {m: [c for c in cs if c.split('__',1)[1] not in ('bp4_benign','pp3_pathogenic','link','gene_link')] for m,cs in module_columns.items() if cs}
        aggregate = ['SUM(CASE WHEN '+filled(c)+' THEN 1 ELSE 0 END)' for c in names]
        aggregate += ['SUM(CASE WHEN '+(' OR '.join(filled(c) for c in cs) if cs else '0')+' THEN 1 ELSE 0 END)' for cs in evidence_columns.values()]
        # Aggregate once for all columns: bounded memory, read-only input.
        counts = conn.execute('SELECT ' + ','.join(aggregate) + ' FROM ' + q(table)).fetchone()
        module_counts = dict(zip(evidence_columns, counts[len(names):]))
        for c, n in zip(names, counts):
            columns.append({'table': table, 'module': c.split('__')[0], 'column': c, 'total_rows': total, 'filled_rows': n or 0, 'empty_rows': total-(n or 0)})
        for module in EXPECTED:
            cols = module_columns[module]
            if not cols:
                continue
            found.add(module)
            # Default ACMG flags and generated hyperlinks alone do not prove a
            # source score/assertion was matched. Count evidence columns instead.
            n = module_counts[module] or 0
            summaries.append({'module': module, 'table': table, 'total_rows': total, 'rows_with_any_data': n, 'rows_without_data': total-n, 'columns': len(cols), 'status': 'PRESENT_WITH_DATA' if n else 'PRESENT_NO_POPULATED_VALUES', 'empty_reason': 'UNRESOLVED_NO_HIT_OR_OUT_OF_SCOPE_OR_ERROR_NOT_INFERRED_FROM_NULL'})
    for module in EXPECTED:
        if module not in found:
            summaries.append({'module': module, 'table': '', 'total_rows': '', 'rows_with_any_data': '', 'rows_without_data': '', 'columns': 0, 'status': 'MODULE_COLUMNS_ABSENT', 'empty_reason': 'Not found in result; compare completed-job module metadata and logs'})
    for filename, rows in [('module-coverage.tsv', summaries), ('column-coverage.tsv', columns)]:
        with (output/filename).open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    sample_membership = {'status': 'SAMPLE_TABLE_ABSENT'}
    if 'sample' in tables:
        sample_columns = [r[1] for r in conn.execute('PRAGMA table_info(sample)')]
        id_col = next((c for c in sample_columns if c.endswith('__sample_id')), None)
        uid_col = next((c for c in sample_columns if c.endswith('__uid')), None)
        if id_col and uid_col:
            grouped = conn.execute('SELECT '+q(id_col)+',COUNT(*),COUNT(DISTINCT '+q(uid_col)+') FROM sample GROUP BY '+q(id_col)).fetchall()
            sample_membership = {'status':'SAMPLE_UID_MAPPING_PRESENT','sample_id_column':id_col,'variant_uid_column':uid_col,'samples':[{'sample_id':r[0],'mapping_rows':r[1],'distinct_variant_uids':r[2]} for r in grouped]}
            variant_columns = [r[1] for r in conn.execute('PRAGMA table_info(variant)')]
            variant_uid = next((c for c in variant_columns if c.endswith('__uid')), None)
            if variant_uid:
                sample_membership['orphan_mapping_rows'] = conn.execute('SELECT COUNT(*) FROM sample s LEFT JOIN variant v ON s.'+q(uid_col)+'=v.'+q(variant_uid)+' WHERE v.'+q(variant_uid)+' IS NULL').fetchone()[0]
        else:
            sample_membership = {'status':'UNRECOGNIZED_SAMPLE_SCHEMA','columns':sample_columns}
    report = {'source': str(db.resolve()), 'source_bytes': db.stat().st_size, 'database_open_mode': 'READ_ONLY', 'expected_modules': EXPECTED, 'metadata': metadata, 'module_coverage': summaries, 'null_semantics': 'NULL is not a benign prediction and does not by itself distinguish no hit, not applicable, module failure or missing data.', 'source': str(db.resolve()), 'official_schema_reference': 'https://docs.opencravat.org/en/latest/Annotator-Tutorial.html'}
    report['sample_membership'] = sample_membership
    (output/'annotation-audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    conn.close()
    print(json.dumps({'database':str(db), 'modules_with_columns':len(found),'modules_absent':[m for m in EXPECTED if m not in found], 'audit':str(output/'annotation-audit.json')},ensure_ascii=False))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('database',type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();audit(a.database,a.output)
