#!/usr/bin/env python3
"""Local sequence-derived stop geometry for an explicitly supplied transcript model.

Not an NMD/PVS1 classifier. Does not fetch sequence or infer a frameshift stop
from the variant position. All coordinates are 1-based on spliced transcript RNA
represented as DNA, already oriented 5' to 3'.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re

STOPS = {'TAA', 'TAG', 'TGA'}


def first_stop(sequence, start):
    return next((position + 3 for position in range(start - 1, len(sequence) - 2, 3)
                 if sequence[position:position + 3] in STOPS), None)


def assess(model):
    result = dict(status='UNASSESSED_REQUIRES_VALIDATED_TRANSCRIPT_INPUTS', nmd_prediction='NOT_ASSIGNED',
                  automatic_PVS1=False, automatic_LoF_GoF=False, per_sample_effect='NOT_INFERRED_FROM_GEOMETRY')
    required = {'transcript', 'hgnc_id', 'source_url_or_local_reference', 'reference_sha256',
                'spliced_cdna', 'cds_start', 'cds_end', 'exon_ends', 'edit', 'consequence_kind',
                'complete_cds', 'reference_genome', 'orientation', 'biological_relevance_review', 'chromosome', 'genetic_code'}
    if not required <= model.keys():
        result['missing_inputs'] = sorted(required - model.keys()); return result
    if not re.fullmatch(r'(ENST\d+|NM_\d+)\.\d+', model['transcript']):
        raise ValueError('Exact versioned transcript required')
    if not re.fullmatch(r'HGNC:\d+', model['hgnc_id']) or model['reference_genome'] != 'GRCh38' or model['orientation'] != 'transcript_5prime_to_3prime':
        raise ValueError('Exact nuclear HGNC/GRCh38 transcript orientation required')
    if model['complete_cds'] is not True:
        result['status'] = 'UNASSESSED_PARTIAL_CDS'; return result
    if str(model['chromosome']).removeprefix('chr') not in [str(i) for i in range(1, 23)] + ['X', 'Y'] or model['genetic_code'] != 'standard_nuclear':
        result['status'] = 'UNASSESSED_NON_NUCLEAR_GENETIC_CODE'; return result
    if model['consequence_kind'] not in ('stop_gained', 'frameshift'):
        result['status'] = 'UNASSESSED_UNSUPPORTED_SPLICE_START_LOST_STOP_LOST_MT_OR_COMPLEX'; return result
    sequence = model['spliced_cdna'].upper()
    if len(sequence) > 1000000:
        raise ValueError('Shortlist transcript model exceeds bounded 1Mb method input')
    if not re.fullmatch('[ACGT]+', sequence) or hashlib.sha256(sequence.encode()).hexdigest() != model['reference_sha256']:
        raise ValueError('Transcript sequence/hash invalid')
    start, end = model['cds_start'], model['cds_end']
    ends = model['exon_ends']
    if not isinstance(start, int) or not isinstance(end, int) or not 1 <= start < end <= len(sequence) or (end-start+1) % 3:
        raise ValueError('Complete CDS bounds/frame invalid')
    if not ends or ends != sorted(set(ends)) or ends[-1] != len(sequence) or any(not isinstance(e, int) or e <= 0 for e in ends):
        raise ValueError('All spliced exon ends including final transcript base required')
    if sequence[start-1:start+2] != 'ATG' or first_stop(sequence, start) != end:
        raise ValueError('Reference CDS start or annotated terminal stop differs')
    edit = model['edit']; position, ref, alt = edit['cdna_start'], edit['ref'].upper(), edit['alt'].upper()
    if not ref or not re.fullmatch('[ACGT]+', ref) or not re.fullmatch('[ACGT]*', alt):
        raise ValueError('Explicit cDNA reference replacement required')
    last = position + len(ref) - 1
    if not start+3 <= position <= last <= end or sequence[position-1:last] != ref:
        raise ValueError('cDNA REF mismatch or unsupported noncoding/start edit')
    if any(position <= boundary < last for boundary in ends[:-1]):
        raise ValueError('Edit crosses exon junction; splice outcome requires a separate model')
    delta = len(alt) - len(ref)
    if model['consequence_kind'] == 'frameshift' and delta % 3 == 0:
        raise ValueError('Frameshift label differs from explicit edit length')
    mutant = sequence[:position-1] + alt + sequence[last:]
    stop = first_stop(mutant, start)
    boundaries = [boundary + (delta if boundary >= last else 0) for boundary in ends[:-1]]
    result.update(status='SEQUENCE_DERIVED_STOP_GEOMETRY_REQUIRES_GENE_SPECIFIC_REVIEW', transcript=model['transcript'],
                  hgnc_id=model['hgnc_id'], reference_sha256=model['reference_sha256'],
                  mutant_sequence_sha256=hashlib.sha256(mutant.encode()).hexdigest(), reference_stop_end_cdna=end,
                  mutant_stop_end_cdna=stop, mutant_final_junction_cdna=boundaries[-1] if boundaries else None,
                  coordinate_convention='1-based stop last nucleotide; junction is last nucleotide of preceding exon; edited transcript coordinates',
                  biological_relevance_review=model['biological_relevance_review'])
    if stop is None:
        result['geometry'] = 'NO_IN_FRAME_STOP_IN_SUPPLIED_TRANSCRIPT'; return result
    # Compare protein stop codon index, rather than unadjusted variant position.
    result['new_stop_before_reference_protein_stop'] = (stop-start)//3 < (end-start)//3
    if not boundaries:
        result['geometry'] = 'SINGLE_EXON_NO_EXON_JUNCTION_HEURISTIC'; return result
    distance = boundaries[-1] - stop
    result.update(distance_stop_end_to_final_junction_nt=distance,
                  geometry='STOP_IN_FINAL_EXON_OR_AT_FINAL_JUNCTION' if distance <= 0 else 'STOP_UPSTREAM_OF_FINAL_JUNCTION',
                  more_than_50nt_upstream=distance > 50, more_than_55nt_upstream=distance > 55,
                  threshold_scope='Descriptive distances only; biological transcript relevance, applicable CSpec, RNA evidence and exceptions still require review')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--model', type=Path, required=True)
    args = parser.parse_args(); print(json.dumps(assess(json.loads(args.model.read_text(encoding='utf-8-sig'))), indent=2))
