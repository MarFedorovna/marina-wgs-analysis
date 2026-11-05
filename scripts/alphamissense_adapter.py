#!/usr/bin/env python3
"""Strict offline AlphaMissense identity matching; model labels are not ACMG classes."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

FIELDS = ['AlphaMissense_score', 'AlphaMissense_model_label', 'AlphaMissense_match_status',
          'AlphaMissense_source_transcript', 'AlphaMissense_source_protein_variant',
          'AlphaMissense_source_uniprot', 'AlphaMissense_source_dataset',
          'AlphaMissense_version', 'AlphaMissense_research_rank', 'AlphaMissense_allele_predictions_json']


def need(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 ** 2), b''):
            h.update(block)
    return h.hexdigest()


def key(chrom, pos, ref, alt):
    chrom = chrom[3:] if chrom.startswith('chr') else chrom
    return ':'.join(('MT' if chrom == 'M' else chrom, str(int(pos)), ref, alt))


def path_in(root, relative):
    path = root / relative
    need(not Path(relative).is_absolute() and '..' not in Path(relative).parts
         and not path.is_symlink() and path.resolve(strict=True).is_relative_to(root), 'Unsafe AlphaMissense resource path.')
    return path


def load_manifest(path):
    path = Path(path)
    manifest = json.loads(path.read_text())
    need(manifest.get('schema_version') == 1 and manifest.get('status') == 'PASS'
         and manifest.get('assembly') == 'GRCh38' and {'canonical', 'isoforms'}.issubset(manifest.get('resources', {})),
         'Both verified GRCh38 AlphaMissense canonical and isoform resources are required.')
    root = path.parent.resolve()
    for artifact in manifest['artifacts']:
        file = path_in(root, artifact['path'])
        need(file.stat().st_size == artifact['bytes'] and digest(file) == artifact['sha256'], 'AlphaMissense resource/index hash changed.')
    return manifest


def manifest_location(cfg):
    return Path(cfg['annotation'].get('alphamissense_manifest') or
        Path(cfg.get('data_root', '/srv/marina')) / 'resources/annotation-supplements/alphamissense-2023-v2-grch38-v3/resource-manifest.json')


def prepare_matches(cfg, workdir, runtime, input_vcf):
    config = cfg['annotation']
    manifest_path = manifest_location(cfg)
    manifest = load_manifest(manifest_path)
    root = manifest_path.parent.resolve()
    files, streams, previous = {}, {}, None
    try:
        for name in ('canonical', 'isoforms'):
            region = workdir / ('alphamissense-' + name + '.bed')
            files[name] = dict(regions=region, matches=workdir / ('alphamissense-' + name + '.matches.tsv'))
            streams[name] = region.open('x')
        with gzip.open(input_vcf, 'rt') as stream:
            for line in stream:
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                if len(ref) != 1 or len(alt) != 1 or not re.fullmatch('[ACGT]{2}', ref + alt):
                    continue
                normalized = chrom[3:] if chrom.startswith('chr') else chrom
                source_chrom = 'chrM' if normalized in ('M', 'MT') else 'chr' + normalized
                coordinate = source_chrom, pos
                if coordinate == previous:
                    continue
                previous = coordinate
                for name, output in streams.items():
                    if source_chrom in manifest['resources'][name]['contigs']:
                        output.write(f'{source_chrom}\t{pos-1}\t{pos}\n')
    finally:
        for stream in streams.values():
            stream.close()
    tabix = manifest['tabix']['binary']
    need(digest(Path(tabix)) == manifest['tabix']['sha256'], 'Prepared AlphaMissense tabix binary changed.')
    for name, paths in files.items():
        if paths['regions'].stat().st_size:
            runtime.run([tabix, '-h', '-R', str(paths['regions']), str(path_in(root, manifest['resources'][name]['path']))],
                logname='alphamissense-local-query-' + name, stdout=paths['matches'], timeout=config.get('alphamissense_query_timeout', 3600))
        else:
            paths['matches'].write_text('#' + '\t'.join(manifest['schemas'][name]['columns']) + '\n')
    return manifest, files, manifest_path


class Lookup:
    def __init__(self, connection, files, manifest):
        self.db, self.manifest = connection, manifest
        connection.execute('CREATE TABLE alpha (variant TEXT, transcript TEXT, protein_variant TEXT, uniprot TEXT, score REAL, label TEXT, dataset TEXT)')
        self.source_rows = 0
        for name, paths in files.items():
            with paths['matches'].open() as stream:
                header = None
                for line in stream:
                    if line.startswith('#CHROM'):
                        header = line.lstrip('#').rstrip('\r\n').split('\t')
                        break
                need(header, 'AlphaMissense matched-file schema missing.')
                for row in csv.DictReader(stream, fieldnames=header, delimiter='\t'):
                    score = float(row['am_pathogenicity'])
                    need(row['genome'] == 'hg38' and 0 <= score <= 1
                         and row['am_class'] in ('likely_benign', 'ambiguous', 'likely_pathogenic'), 'Invalid AlphaMissense matched data.')
                    connection.execute('INSERT INTO alpha VALUES (?,?,?,?,?,?,?)',
                        (key(row['CHROM'], row['POS'], row['REF'], row['ALT']), row['transcript_id'], row['protein_variant'],
                         row.get('uniprot_id', ''), score, row['am_class'], name))
                    self.source_rows += 1
        connection.execute('CREATE INDEX alpha_exact_variant ON alpha (variant)')
        connection.commit()

    def match(self, variant, csq, threshold=.564):
        rows = [dict(zip(('transcript', 'protein_variant', 'uniprot', 'score', 'model_label', 'dataset'), row))
            for row in self.db.execute('SELECT transcript,protein_variant,uniprot,score,label,dataset FROM alpha WHERE variant=?', (variant,))]
        result = dict.fromkeys(FIELDS, '')
        result.update(AlphaMissense_version=self.manifest['dataset_version'], AlphaMissense_research_rank=0,
            AlphaMissense_allele_predictions_json=json.dumps(rows, separators=(',', ':'), sort_keys=True))
        if 'missense_variant' not in csq.get('Consequence', '').split('&'):
            result['AlphaMissense_match_status'] = 'not_applicable_to_this_nonmissense_transcript; allele_context_preserved'
            return result
        if not rows:
            result['AlphaMissense_match_status'] = 'not_covered_in_source; not_benign'
            return result
        transcript = csq.get('Feature', '')
        exact_tx = [row for row in rows if row['transcript'] == transcript]
        if not re.fullmatch(r'ENST\d+\.\d+', transcript):
            result['AlphaMissense_match_status'] = 'current_transcript_version_unavailable; no_score_assignment'
            return result
        if not exact_tx:
            same_stable = any(row['transcript'].split('.')[0] == transcript.split('.')[0] for row in rows)
            result['AlphaMissense_match_status'] = 'transcript_version_mismatch; allele_context_only' if same_stable else 'no_exact_transcript_match; allele_context_only'
            return result
        aa = csq.get('Amino_acids', '').split('/')
        position = csq.get('Protein_position', '').split('/')[0]
        if len(aa) != 2 or not all(re.fullmatch('[ARNDCEQGHILKMFPSTWYV]', value) for value in aa) or not position.isdigit():
            result['AlphaMissense_match_status'] = 'current_protein_substitution_unavailable; no_score_assignment'
            return result
        protein_variant = aa[0] + position + aa[1]
        exact = [row for row in exact_tx if row['protein_variant'] == protein_variant]
        if not exact:
            result['AlphaMissense_match_status'] = 'protein_substitution_mismatch; no_score_assignment'
            return result
        accessions = set()
        for field in ('SWISSPROT', 'TREMBL', 'UNIPROT_ISOFORM'):
            accessions.update(filter(None, re.split('[&,]', csq.get(field, ''))))
        compatible = [row for row in exact if not row['uniprot'] or row['uniprot'] in accessions]
        if not compatible:
            result['AlphaMissense_match_status'] = 'uniprot_accession_mismatch_or_mapping_unavailable; no_score_assignment'
            return result
        unique = {(row['score'], row['model_label']) for row in compatible}
        if len(unique) != 1:
            result['AlphaMissense_match_status'] = 'conflicting_exact_source_predictions; no_score_assignment'
            return result
        score, label = next(iter(unique))
        result.update(AlphaMissense_score=score, AlphaMissense_model_label=label,
            AlphaMissense_match_status='exact_allele_transcript_version_protein_change; model_not_ACMG',
            AlphaMissense_source_transcript=transcript, AlphaMissense_source_protein_variant=protein_variant,
            AlphaMissense_source_uniprot='&'.join(sorted({row['uniprot'] for row in compatible if row['uniprot']})),
            AlphaMissense_source_dataset='&'.join(sorted({row['dataset'] for row in compatible})),
            AlphaMissense_research_rank=int(score > threshold))
        return result
