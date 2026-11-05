#!/usr/bin/env python3
"""Local P09 predictor-alignment audit and P11 phenotype/pedigree completeness checks."""
import argparse
import json
import math
from pathlib import Path
import re


def allele(value):
    chrom = str(value['chrom']).removeprefix('chr')
    position = int(value['pos'])
    ref, alt = str(value['ref']).upper(), str(value['alt']).upper()
    if position < 1 or not re.fullmatch('[ACGT]+', ref) or not re.fullmatch('[ACGT]+', alt) or ref == alt:
        raise ValueError('Explicit normalized genomic allele required')
    return chrom, position, ref, alt


def revel_alignment(source, current):
    result = dict(status='UNALIGNED_RAW_REFERENCE_ONLY', canonical_transcript_score=None,
                  source='REVEL', version=source.get('REVEL_version'), automatic_PP3_BP4=False,
                  independent_clinical_vote=False, source_transcript_versions_invented=False)
    if source.get('assembly') != 'GRCh38' or current.get('assembly') != 'GRCh38' or allele(source) != allele(current):
        result['status'] = 'ALLELE_OR_ASSEMBLY_MISMATCH'; return result
    transcript = current.get('transcript', '')
    if not re.fullmatch(r'ENST\d+\.\d+', transcript):
        result['status'] = 'CURRENT_VERSIONED_TRANSCRIPT_REQUIRED'; return result
    source_transcripts = source.get('REVEL_source_transcript_ids', '').split(';')
    if transcript not in source_transcripts and transcript.split('.')[0] not in source_transcripts:
        result['status'] = 'SOURCE_TRANSCRIPT_NOT_MATCHED'; return result
    if current.get('amino_acids', '').split('/') != [source.get('REVEL_source_aa_ref'), source.get('REVEL_source_aa_alt')]:
        result['status'] = 'AMINO_ACID_MISMATCH'; return result
    score = float(source['REVEL_reference_score'])
    if not math.isfinite(score) or not 0 <= score <= 1:
        raise ValueError('REVEL source score outside unit interval')
    result.update(status='ALLELE_AND_AA_MATCH_SOURCE_TRANSCRIPT_VERSION_UNREPORTED_RESEARCH_ONLY', raw_reference_score=score,
                  score_scope='Source allele/AA context; not a validated scalar for every current VEP transcript')
    if transcript in source_transcripts:
        # Version identity improves the audit, but does not itself calibrate PP3/BP4.
        result['status'] = 'EXACT_ALLELE_AA_TRANSCRIPT_VERSION_MATCH_RESEARCH_ONLY'
    return result


def dbnsfp_alignment(raw, predictor, current):
    """Audit source vectors without assuming their ordering defines a valid zip."""
    result = dict(predictor=predictor, status='RAW_SOURCE_ONLY', canonical_transcript_score=None,
                  automatic_PP3_BP4=False, independent_clinical_vote=False)
    if raw.get('reference_genome') != 'GRCh38' or current.get('assembly') != 'GRCh38' or allele(raw) != allele(current):
        result['status'] = 'ALLELE_OR_ASSEMBLY_MISMATCH'; return result
    data = raw.get('dbnsfp', {})
    if predictor not in data:
        result['status'] = 'MISSING_SOURCE_NOT_ZERO_OR_BENIGN'; return result
    transcripts = data.get('ensembl', {}).get('transcriptid', [])
    positions = data.get('aa', {}).get('pos', [])
    scores = data[predictor].get('score', [])
    arrays = [v if isinstance(v, list) else [v] for v in (transcripts, positions, scores)]
    result['vector_lengths'] = dict(zip(('transcripts', 'aa_positions', 'scores'), map(len, arrays)))
    if len({len(v) for v in arrays}) != 1:
        result['status'] = 'NONALIGNED_CARDINALITY_DO_NOT_ZIP'; return result
    # Equal lengths alone do not establish the published source array contract.
    # No scalar is emitted until ordering/protein/transcript version provenance
    # has been independently documented for this predictor's source schema.
    result.update(status='EQUAL_LENGTH_ARRAYS_PARALLEL_CONTRACT_AND_VERSIONS_UNVERIFIED_DO_NOT_ZIP',
                  raw_source_preserved=True, legacy_FATHMM_XF_MKL_are_separate=True)
    return result


def validate_phenotype_pedigree(document):
    errors = []
    if document.get('no_external_upload') is not True:
        errors.append('no_external_upload must be true')
    subjects = document.get('subjects', [])
    ids = [s.get('subject_id') for s in subjects]
    if not ids or len(set(ids)) != len(ids) or any(not isinstance(i, str) or not i for i in ids):
        errors.append('Unique local subject IDs required')
    samples = [sample for subject in subjects for sample in subject.get('sample_ids', [])]
    if len(samples) != len(set(samples)):
        errors.append('A sample cannot silently identify multiple subjects')
    for subject in subjects:
        sid = str(subject.get('subject_id'))
        for feature in subject.get('phenotypes', []):
            if not re.fullmatch(r'HP:\d{7}', str(feature.get('hpo_id', ''))):
                errors.append(sid + ': HPO syntax invalid; ontology membership needs local HPO release check')
            if feature.get('status') not in ('present', 'absent', 'unknown'):
                errors.append(sid + ': phenotype status must be present/absent/unknown')
            if feature.get('status') in ('present', 'absent') and not feature.get('source'):
                errors.append(sid + ': observed or excluded phenotype requires source/date/examination context')
        for field in ('reported_sex', 'karyotype', 'age_at_last_review'):
            value = subject.get(field, {})
            if value.get('value') not in (None, '', 'unknown') and not value.get('source'):
                errors.append(sid + ': ' + field + ' requires source; technical XY is not clinical sex')
        for examination in subject.get('cardiac_evaluations', []):
            if not all(examination.get(k) for k in ('date', 'method', 'source')):
                errors.append(sid + ': cardiac evaluation needs date/method/source')
            for measurement in examination.get('measurements', []):
                if not measurement.get('name') or measurement.get('value') is None or not measurement.get('unit'):
                    errors.append(sid + ': cardiac measurement needs name/value/unit')
    for relationship in document.get('relationships', []):
        if relationship.get('subject_a') not in ids or relationship.get('subject_b') not in ids or relationship.get('subject_a') == relationship.get('subject_b'):
            errors.append('Relationship refers to missing or identical subjects')
        if relationship.get('status') not in ('unknown', 'reported', 'confirmed'):
            errors.append('Relationship status must be unknown/reported/confirmed')
        if relationship.get('status') in ('reported', 'confirmed') and not relationship.get('source'):
            errors.append('Relationship assertion needs source; four samples imply no pedigree')
    for evidence in document.get('segregation_records', []):
        if evidence.get('subject_id') not in ids:
            errors.append('Segregation refers to missing subject')
        if evidence.get('de_novo_status') not in ('unknown', 'suspected', 'confirmed', 'not_de_novo'):
            errors.append('Explicit de novo status required')
        if evidence.get('de_novo_status') == 'confirmed' and not (
                evidence.get('parentage_confirmed') is True and evidence.get('both_parents_tested') is True
                and len(set(evidence.get('parent_sample_ids', []))) == 2 and evidence.get('source')):
            errors.append('Confirmed de novo requires both tested parents, parentage and source evidence')
        if evidence.get('phase') not in ('unknown', 'cis', 'trans'):
            errors.append('Explicit phase state required')
        if evidence.get('phase') in ('cis', 'trans') and not evidence.get('phase_evidence_source'):
            errors.append('cis/trans requires phase evidence; sample overlap is insufficient')
    return dict(status='PASS_INPUT_STRUCTURE_ONLY_NO_CLINICAL_INFERENCE' if not errors else 'INCOMPLETE_OR_INCONSISTENT',
                errors=errors, subjects=len(subjects), hpo_ontology_membership='NOT_CHECKED_BY_SYNTAX_VALIDATOR',
                inferred_relationships=0, clinical_classification=None, network_requests=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', required=True, choices=['phenotype-pedigree', 'revel-alignment', 'dbnsfp-alignment'])
    parser.add_argument('--input', type=Path, required=True)
    args = parser.parse_args(); data = json.loads(args.input.read_text(encoding='utf-8-sig'))
    if args.mode == 'phenotype-pedigree': result = validate_phenotype_pedigree(data)
    elif args.mode == 'revel-alignment': result = revel_alignment(data['source'], data['current'])
    else: result = dbnsfp_alignment(data['source'], data['predictor'], data['current'])
    print(json.dumps(result, indent=2))
