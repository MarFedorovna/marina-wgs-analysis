#!/usr/bin/env python3
"""Standalone, local-only consequence tables. Does not modify workflow stages."""
import argparse
import csv
import gzip
import hashlib
import html
import itertools
import json
from pathlib import Path
import re
import sqlite3
import sys
from contextlib import contextmanager

CORE = ('stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant',
        'splice_donor_variant', 'splice_acceptor_variant', 'splice_region_variant',
        'missense_variant', 'synonymous_variant', 'inframe_insertion', 'inframe_deletion')
PLOF = {'stop_gained', 'frameshift_variant', 'splice_donor_variant', 'splice_acceptor_variant'}
SPECIAL = 'pLoF_candidates'
REQUIRED = {'sample', 'variant_id', 'chrom', 'pos', 'ref', 'alt', 'GT', 'consequence'}
SUMMARY = ['export_allele_id', 'variant_id', 'chrom', 'pos', 'ref', 'alt',
           'observed_alt_samples_json', 'unknown_gt_samples_json', 'reference_gt_samples_json',
           'partial_alt_samples_json', 'invalid_gt_samples_json', 'genotypes_json',
           'quality_levels_json', 'source_FT_json', 'source_site_filter_json',
           'source_sample_state_json', 'source_candidate_eligible_json',
           'source_mechanism_status_json', 'mechanism_review_statuses_json',
           'source_transcript_rows', 'evidence_file', 'evidence_sample_id',
           'clinical_interpretation_status']


def js(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def gt_flags(gt):
    """Observed ALT is independent of missingness, quality, and source sample_state."""
    tokens = re.split(r'[/|]', gt or '.')
    valid = all(t == '.' or re.fullmatch(r'\d+', t) for t in tokens)
    alt = any(t.isdigit() and int(t) > 0 for t in tokens)
    missing = not valid or any(t == '.' for t in tokens)
    return alt, missing, valid and not missing and not alt, alt and missing, not valid


def category_slug(term):
    if term in CORE or term == SPECIAL:
        return term
    slug = re.sub(r'[^A-Za-z0-9_-]+', '_', term).strip('_')[:65] or 'term'
    return slug + '-' + hashlib.sha256(term.encode()).hexdigest()[:12]


def safe_label(value):
    """Only derived display labels: exact clinical evidence stays in raw gzip TSV."""
    return "'" + value if value.lstrip().startswith(('=', '+', '-', '@')) else value


@contextmanager
def tsv_writer(path, fields):
    with gzip.open(path, 'wt', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(fields)
        yield writer


def input_stream(path):
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if str(path).endswith('.gz') else Path(path).open(encoding='utf-8', newline='')


def summarize(rows, allele, scope_sid, samples):
    selected = rows if scope_sid == 0 else [r for r in rows if r['sid'] == scope_sid]
    flags = {key: set() for key in ('alt', 'unknown', 'reference', 'partial', 'invalid')}
    genotypes, quality = {}, set()
    for row in selected:
        sample = samples[row['sid']]
        genotypes.setdefault(sample, set()).add(row['gt'])
        for name, flag in zip(flags, gt_flags(row['gt'])):
            if flag:
                flags[name].add(sample)
        if row['ft'] in ('', '.') or row['filter'] in ('', '.'):
            quality.add('UNKNOWN_QUALITY')
        if row['ft'] not in ('', '.', 'PASS') or row['filter'] not in ('', '.', 'PASS'):
            quality.add('FILTERED')
        if row['ft'] == 'PASS' and row['filter'] == 'PASS':
            quality.add('PASS')
    values = lambda column: js(sorted({r[column] for r in selected}))
    review = sorted({r['review'] for r in selected if r['review']}) or ['UNASSESSED']
    return [*allele, *(js(sorted(flags[k])) for k in flags),
            js({s: sorted(v) for s, v in genotypes.items()}), js(sorted(quality)),
            values('ft'), values('filter'), values('state'), values('candidate'),
            values('mechanism'), js(review), sum(r['n'] for r in selected),
            '../../evidence/transcript_effects.tsv.gz', scope_sid or '*',
            'RESEARCH_GROUPING_NOT_DIAGNOSIS']


def export(source, output, *, synthetic=False, receipt=None):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError('Output must be a new directory; refusing overwrite.')
    receipt_data = None
    if receipt:
        receipt_data = json.loads(Path(receipt).read_text(encoding='utf-8'))
        if receipt_data.get('status') != 'PASS':
            raise ValueError('Completed workflow receipt must have status PASS.')
    source_stat = source.stat()
    source_hash = sha256(source)
    output.mkdir(parents=True)
    (output / 'evidence').mkdir()
    db_path = output / '.working.sqlite'
    db = sqlite3.connect(db_path)
    db.row_factory = sqlite3.Row
    db.executescript('''PRAGMA journal_mode=DELETE; PRAGMA temp_store=FILE;
        PRAGMA cache_size=-32768;
        CREATE TABLE alleles (aid INTEGER PRIMARY KEY, variant_id TEXT UNIQUE,
            chrom TEXT, pos TEXT, ref TEXT, alt TEXT);
        CREATE TABLE samples (sid INTEGER PRIMARY KEY, sample TEXT UNIQUE);
        CREATE TABLE categories (term TEXT PRIMARY KEY);
        CREATE TABLE membership (aid INTEGER, sid INTEGER, term TEXT,
            PRIMARY KEY(term,sid,aid)) WITHOUT ROWID;
        CREATE TABLE obs (aid INTEGER, sid INTEGER, gt TEXT, ft TEXT, filter TEXT,
            state TEXT, candidate TEXT, mechanism TEXT, review TEXT, n INTEGER,
            PRIMARY KEY(aid,sid,gt,ft,filter,state,candidate,mechanism,review)) WITHOUT ROWID;
        CREATE TABLE summaries (aid INTEGER, sid INTEGER, payload TEXT,
            PRIMARY KEY(sid,aid)) WITHOUT ROWID;''')
    db.executemany('INSERT INTO categories VALUES (?)', [(x,) for x in (*CORE, SPECIAL)])
    source_rows = 0
    try:
        with input_stream(source) as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            fields = reader.fieldnames or []
            if not REQUIRED.issubset(fields) or len(fields) != len(set(fields)):
                raise ValueError('Missing required fields or duplicate source column names.')
            prefix = ['export_row_id', 'export_allele_id', 'export_sample_id']
            if set(prefix) & set(fields):
                raise ValueError('Reserved export column already present in source.')
            with tsv_writer(output / 'evidence/transcript_effects.tsv.gz', prefix + fields) as raw:
                for number, row in enumerate(reader, 1):
                    if None in row or any(value is None for value in row.values()):
                        raise ValueError('Malformed TSV row %d' % number)
                    if not row['sample'] or not row['variant_id']:
                        raise ValueError('Empty sample or variant_id at row %d' % number)
                    allele = tuple(row[k] for k in ('variant_id', 'chrom', 'pos', 'ref', 'alt'))
                    db.execute('INSERT OR IGNORE INTO alleles(variant_id,chrom,pos,ref,alt) VALUES (?,?,?,?,?)', allele)
                    found = db.execute('SELECT * FROM alleles WHERE variant_id=?', (row['variant_id'],)).fetchone()
                    if tuple(found)[1:] != allele:
                        raise ValueError('variant_id refers to inconsistent genomic alleles.')
                    aid = found['aid']
                    db.execute('INSERT OR IGNORE INTO samples(sample) VALUES (?)', (row['sample'],))
                    sid = db.execute('SELECT sid FROM samples WHERE sample=?', (row['sample'],)).fetchone()[0]
                    raw.writerow([number, aid, sid] + [row[f] for f in fields])
                    terms = {x.strip() for x in re.split(r'[&,;|]', row['consequence']) if x.strip()}
                    terms = terms or {'missing_consequence'}
                    if terms & PLOF or row.get('plof_candidate', '').lower() in ('1', 'true', 'yes'):
                        terms.add(SPECIAL)
                    observed_alt = gt_flags(row['GT'])[0]
                    for term in terms:
                        db.execute('INSERT OR IGNORE INTO categories VALUES (?)', (term,))
                        db.execute('INSERT OR IGNORE INTO membership VALUES (?,?,?)', (aid, 0, term))
                        if observed_alt:
                            db.execute('INSERT OR IGNORE INTO membership VALUES (?,?,?)', (aid, sid, term))
                    values = [aid, sid, row['GT'], row.get('FT', ''), row.get('site_filter', ''),
                              row.get('sample_state', ''), row.get('candidate_eligible', ''),
                              row.get('mechanism_status', ''), row.get('mechanism_review_status', ''), 1]
                    db.execute('''INSERT INTO obs VALUES (?,?,?,?,?,?,?,?,?,?)
                        ON CONFLICT DO UPDATE SET n=n+1''', values)
                    source_rows = number
                    if number % 25000 == 0:
                        db.commit()
            db.commit()
        if source.stat().st_size != source_stat.st_size or sha256(source) != source_hash:
            raise ValueError('Source changed during export; output is incomplete.')
        samples = dict(db.execute('SELECT sid,sample FROM samples ORDER BY sid').fetchall())
        cursor = db.execute('SELECT * FROM obs ORDER BY aid,sid')
        for aid, group in itertools.groupby(cursor, key=lambda r: r['aid']):
            rows = list(group)  # Distinct call states for one allele, not transcript rows.
            allele = tuple(db.execute('SELECT * FROM alleles WHERE aid=?', (aid,)).fetchone())
            for sid in [0] + sorted({r['sid'] for r in rows if gt_flags(r['gt'])[0]}):
                payload = summarize(rows, allele, sid, samples)
                db.execute('INSERT INTO summaries VALUES (?,?,?)', (aid, sid, js(payload)))
        db.commit()
        index = []
        terms = [r[0] for r in db.execute('SELECT term FROM categories ORDER BY term')]
        scopes = [(0, 'cohort', 'cohort')] + [(sid, 'sample', name) for sid, name in sorted(samples.items())]
        for sid, kind, label in scopes:
            scope_dir = 'cohort' if sid == 0 else 'sample_%04d' % sid
            folder = output / 'tables' / scope_dir
            folder.mkdir(parents=True)
            for term in terms:
                path = folder / (category_slug(term) + '.tsv.gz')
                count = 0
                with tsv_writer(path, SUMMARY) as writer:
                    query = '''SELECT s.payload FROM membership m JOIN summaries s
                               ON s.aid=m.aid AND s.sid=m.sid WHERE m.term=? AND m.sid=? ORDER BY m.aid'''
                    for result in db.execute(query, (term, sid)):
                        writer.writerow(json.loads(result[0]))
                        count += 1
                index.append({'scope_kind': kind, 'scope_id': label, 'export_sample_id': sid,
                              'category': term, 'unique_alleles': count,
                              'path': path.relative_to(output).as_posix()})
        with (output / 'index.tsv').open('w', encoding='utf-8', newline='') as stream:
            writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
            writer.writerow(['scope_kind', 'scope_id', 'export_sample_id', 'category', 'unique_alleles', 'path'])
            for item in index:
                writer.writerow([safe_label(str(v)) if k in ('scope_id', 'category') else v for k, v in item.items()])
        note = ('SYNTHETIC FIXTURE — no patient results. ' if synthetic else '') + (
            'Counts are unique alleles per category and scope. Categories overlap across transcripts; '
            'do not sum them. Sample tables require an observed ALT token; partial ALT calls remain flagged. '
            'Cohort tables retain all source alleles, including unknown and reference calls. '
            'pLoF_candidates is an overlapping research group, not proven LoF, PVS1 or pathogenicity. '
            'Stop loss does not imply gain of function. Missing samples are not inferred as reference. '
            'Quality and candidate values are retained separately; filtered calls are not classified benign. '
            'Genotypes and raw evidence are text; use a text import with columns set to text in Excel.')
        body = ['<!doctype html><meta charset="utf-8"><title>Consequence tables</title>',
                '<h1>Consequence tables</h1><p>' + html.escape(note) + '</p>',
                '<p><a href="manifest.json">Manifest</a> · <a href="index.tsv">TSV index</a> · '
                '<a href="evidence/transcript_effects.tsv.gz">Complete transcript evidence</a></p>',
                '<table><tr><th>Scope</th><th>Category</th><th>Unique alleles</th><th>Table</th></tr>']
        for item in index:
            body.append('<tr><td>%s: %s</td><td>%s</td><td>%d</td><td><a href="%s">TSV.gz</a></td></tr>' %
                        (item['scope_kind'], html.escape(item['scope_id']), html.escape(item['category']),
                         item['unique_alleles'], html.escape(item['path'], quote=True)))
        (output / 'index.html').write_text('\n'.join(body) + '</table>\n', encoding='utf-8')
        db.close()
        db_path.unlink()
        artifacts = {p.relative_to(output).as_posix(): {'sha256': sha256(p), 'bytes': p.stat().st_size}
                     for p in sorted(output.rglob('*')) if p.is_file()}
        manifest = {'schema_version': 1, 'status': 'PASS', 'synthetic': synthetic,
                    'exporter_sha256': sha256(Path(__file__)),
                    'source': {'path': str(source), 'sha256': source_hash, 'bytes': source_stat.st_size,
                               'transcript_rows': source_rows, 'fields': fields},
                    'workflow_receipt': receipt_data, 'counts': index, 'artifacts': artifacts,
                    'evidence_join': ['export_allele_id', 'export_sample_id (or * for cohort)'],
                    'clinical_interpretation_status': 'RESEARCH_GROUPING_NOT_DIAGNOSIS',
                    'source_mechanism_review': 'Copied only if supplied; otherwise UNASSESSED',
                    'formula_handling': 'Display labels in index.tsv escaped; raw evidence is exact text; import TSV as text.',
                    'ordering': 'Allele and sample IDs follow first source occurrence; terms sorted lexically.'}
        (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        return manifest
    except BaseException:
        db.close()
        (output / 'INCOMPLETE.txt').write_text('No PASS manifest: do not treat this export as complete.\n', encoding='utf-8')
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--receipt', help='Optional completed source workflow receipt; must be PASS')
    parser.add_argument('--synthetic', action='store_true', help='Mark every artifact index as synthetic fixture results')
    args = parser.parse_args()
    result = export(args.source, args.output, synthetic=args.synthetic, receipt=args.receipt)
    print(js({'status': result['status'], 'source_rows': result['source']['transcript_rows'], 'output': args.output}))


if __name__ == '__main__':
    main()
