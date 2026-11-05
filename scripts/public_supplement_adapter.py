#!/usr/bin/env python3
"""Offline, exact-ID public gene/disease evidence. No patient API or ACMG inference."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import xml.etree.ElementTree as ET


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require_id(value, prefix, digits=None):
    pattern = re.escape(prefix) + ':' + (r'\d{%d}' % digits if digits else r'\d+')
    if not isinstance(value, str) or not re.fullmatch(pattern, value):
        raise ValueError('Exact canonical ID required: ' + prefix)
    return value


def text(node, path):
    return (node.findtext(path) or '').strip()


def disorders(path):
    for _, node in ET.iterparse(path, events=('end',)):
        if node.tag == 'Disorder':
            yield node
            node.clear()


def gene_rows(path):
    for d in disorders(path):
        code = 'ORPHA:' + text(d, 'OrphaCode')
        for a in d.findall('./DisorderGeneAssociationList/DisorderGeneAssociation'):
            g = a.find('Gene')
            if g is None:
                continue
            ids = [text(e, 'Reference') for e in g.findall('./ExternalReferenceList/ExternalReference') if text(e, 'Source') == 'HGNC']
            for value in ids:
                if not value.isdigit():
                    continue
                yield ('HGNC:' + value, text(g, 'Symbol'), code, text(d, 'Name'),
                    text(a, 'DisorderGeneAssociationType/Name'), text(a, 'DisorderGeneAssociationStatus/Name'),
                    text(a, 'SourceOfValidation'))


def xref_rows(path):
    for d in disorders(path):
        for e in d.findall('./ExternalReferenceList/ExternalReference'):
            source, value = text(e, 'Source'), text(e, 'Reference')
            if not source or not value:
                continue
            target = source + ':' + value
            yield ('ORPHA:' + text(d, 'OrphaCode'), text(d, 'Name'), target,
                   text(e, 'DisorderMappingRelation/Name'), text(e, 'DisorderMappingValidationStatus/Name'))


def hpo_rows(path):
    for d in disorders(path):
        for a in d.findall('./HPODisorderAssociationList/HPODisorderAssociation'):
            term = text(a, 'HPO/HPOId')
            require_id(term, 'HP', 7)
            yield ('ORPHA:' + text(d, 'OrphaCode'), term, text(a, 'HPO/HPOTerm'),
                   text(a, 'HPOFrequency/Name'), text(a, 'DiagnosticCriteria'))


def released_cspec(document, mondo='MONDO:0005045'):
    require_id(mondo, 'MONDO', 7)
    seen = set()
    for entity in document['data']:
        content, linked = entity.get('entContent', {}), entity.get('ld', {})
        if not any(s.get('current') is True and s.get('name') == 'Released' for s in content.get('states', [])):
            continue
        diseases = {d.get('entId') for d in linked.get('Disease', [])}
        if mondo not in diseases:
            continue
        for gene in linked.get('Gene', []):
            hgnc = gene.get('entContent', {}).get('HGNC', {})
            hgnc_id = require_id(hgnc.get('hgnc_id'), 'HGNC')
            key = (entity['entId'], content.get('version'), hgnc_id, mondo)
            if key in seen:
                raise ValueError('Duplicate released CSpec gene/version assertion')
            seen.add(key)
            yield dict(spec_id=entity['entId'], version=content.get('version'), title=content.get('title'),
                hgnc_id=hgnc_id, symbol=hgnc.get('symbol'), mondo_id=mondo,
                approved_on=content.get('approvedOn'), state='Released',
                released_under_revision=bool(content.get('releasedUnderRevision')),
                source_url=entity.get('ldhIri'), rules=linked.get('CriteriaCode', []),
                rule_sets=linked.get('RuleSet', []), notes=content.get('releaseNotes'),
                references=content.get('references', []), doi=content.get('doi'),
                automatic_acmg_application=False)


def build(source, target):
    if target.exists():
        raise ValueError('Refusing to overwrite an existing candidate; choose a new version')
    manifest = json.loads((source / 'source-manifest.json').read_text())
    if manifest.get('status') != 'PASS_PUBLIC_SOURCES_ONLY':
        raise ValueError('Public source manifest is not PASS')
    receipts = {i['source']['id']: i for i in manifest['resources']}
    required = ['orphadata_genes', 'orphadata_disease_ids', 'orphadata_phenotypes',
                'clingen_cspec_public', 'clingen_cspec_page2', 'clingen_cspec_page3', 'clingen_cspec_registry_metadata']
    for key in required:
        item = receipts[key]; path = source / item['source']['file']
        if item.get('status') != 'PASS' or path.stat().st_size != item['bytes'] or sha(path) != item['sha256']:
            raise ValueError('Source integrity differs: ' + key)
    target.mkdir(mode=0o700, parents=True)
    dbpath = target / 'public-gene-evidence.sqlite'
    conn = sqlite3.connect(dbpath)
    conn.executescript('''
      CREATE TABLE gene_disease(hgnc_id TEXT,symbol TEXT,orpha_id TEXT,disease TEXT,association_type TEXT,assessment TEXT,validation TEXT,PRIMARY KEY(hgnc_id,orpha_id,association_type,assessment,validation));
      CREATE TABLE disease_xref(orpha_id TEXT,disease TEXT,target_id TEXT,relation TEXT,validation TEXT,PRIMARY KEY(orpha_id,target_id,relation,validation));
      CREATE TABLE disease_hpo(orpha_id TEXT,hpo_id TEXT,label TEXT,frequency TEXT,diagnostic_criteria TEXT,PRIMARY KEY(orpha_id,hpo_id,frequency,diagnostic_criteria));
      CREATE TABLE cspec(hgnc_id TEXT,mondo_id TEXT,spec_id TEXT,version TEXT,payload TEXT,PRIMARY KEY(hgnc_id,mondo_id,spec_id,version));
      CREATE TABLE provenance(resource_id TEXT PRIMARY KEY,payload TEXT);
      CREATE INDEX gene_disease_hgnc ON gene_disease(hgnc_id);
      CREATE INDEX disease_xref_target ON disease_xref(target_id);
    ''')
    stats = {}
    for key, table, rows, count in [
        ('orphadata_genes','gene_disease',gene_rows(source/'en_product6.xml'),7),
        ('orphadata_disease_ids','disease_xref',xref_rows(source/'en_product1.xml'),5),
        ('orphadata_phenotypes','disease_hpo',hpo_rows(source/'en_product4.xml'),5)]:
        attempted = 0
        for row in rows:
            conn.execute('INSERT OR IGNORE INTO %s VALUES(%s)' % (table, ','.join('?' * count)), row); attempted += 1
        stored = conn.execute('SELECT COUNT(*) FROM ' + table).fetchone()[0]
        stats[key] = dict(source_rows=attempted, unique_rows=stored, exact_duplicate_rows=attempted-stored)
    entities = []
    for name in ('cspec-public.json','cspec-public-page2.json','cspec-public-page3.json'):
        entities.extend(json.loads((source/name).read_text())['data'])
    expected = json.loads((source/'cspec-registry-metadata.json').read_text())['data']['entCount']
    ids = [i['entId'] for i in entities]
    if len(ids) != len(set(ids)) or len(ids) != expected:
        raise ValueError('CSpec pagination incomplete, duplicated, or registry changed; new snapshot review required')
    specs = list(released_cspec({'data':entities}))
    for item in specs:
        conn.execute('INSERT INTO cspec VALUES(?,?,?,?,?)', (item['hgnc_id'],item['mondo_id'],item['spec_id'],item['version'],json.dumps(item)))
        (target/(item['spec_id']+'-'+item['version']+'.json')).write_text(json.dumps(item,indent=2)+'\n')
    stats['clingen_cspec_public'] = dict(registry_entities=len(ids), registry_expected=expected,
        released_hcm_gene_specifications=len(specs), genes=[s['symbol'] for s in specs])
    for key, value in receipts.items():
        conn.execute('INSERT INTO provenance VALUES(?,?)',(key,json.dumps(value)))
    conn.commit(); conn.close()
    files = [dict(path=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(target.iterdir()) if p.is_file()]
    result = dict(status='CANDIDATE_PUBLIC_ADAPTER_NOT_BOUND', automatic_acmg=False,
        source_manifest_sha256=sha(source/'source-manifest.json'), statistics=stats, files=files,
        source_versions={k:v['source']['release'] for k,v in receipts.items()},
        annotation_binding_updated=False, no_patient_variants_processed=True,
        limitations=['Disease/gene context is not variant pathogenicity or phenotype fit.',
            'Exact validated MONDO mappings only for disease joins; other mappings retained as source assertions.',
            'CSpec applicability and ACMG strength need expert review; no criteria are auto-applied.',
            'Orphanet OMIM cross-reference IDs are not downloaded OMIM content.'])
    (target/'candidate-manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def gene_evidence(database, hgnc_id, mondo_id):
    require_id(hgnc_id, 'HGNC'); require_id(mondo_id, 'MONDO', 7)
    conn = sqlite3.connect(database.resolve().as_uri() + '?mode=ro', uri=True)
    conn.row_factory = sqlite3.Row
    rows = conn.execute('''SELECT g.*,x.target_id,x.relation,x.validation AS mapping_validation
      FROM gene_disease g JOIN disease_xref x USING(orpha_id)
      WHERE g.hgnc_id=? AND x.target_id=? AND x.relation LIKE 'E (%' AND x.validation='Validated' ''', (hgnc_id,mondo_id)).fetchall()
    cspec = [json.loads(row[0]) for row in conn.execute('SELECT payload FROM cspec WHERE hgnc_id=? AND mondo_id=?',(hgnc_id,mondo_id))]
    provenance = {k:json.loads(v) for k,v in conn.execute('SELECT resource_id,payload FROM provenance')}
    conn.close()
    return dict(status='PUBLIC_GENE_CONTEXT' if rows or cspec else 'MISSING_EXACT_GENE_DISEASE_CONTEXT',
        hgnc_id=hgnc_id,mondo_id=mondo_id,orphadata_assertions=[dict(row) for row in rows],cspec=cspec,
        provenance=provenance, pathogenicity_classification=None, automatic_acmg=False)


def annotate_variant_context(database, variant):
    # Offline future integration adapter. Exact IDs are required; no symbol/position guessing.
    result = gene_evidence(database, variant['hgnc_id'], variant['mondo_id'])
    result['variant_assertion'] = 'NOT_PROVIDED_BY_GENE_CONTEXT_DATABASE'
    result['variant_key'] = variant.get('variant_key')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=Path('/srv/marina/resources/public-supplements-v1'))
    parser.add_argument('--build-candidate',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(build(args.source,args.build_candidate),ensure_ascii=False,indent=2))
