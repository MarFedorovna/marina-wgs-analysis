#!/usr/bin/env python3
"""Synthetic-only FAA filtered-row and circular-boundary regression."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path

import recover_mitochondria as recovery
from calling_stages import require
from workflow_common import json_write


def validate(fixture, resource_path, destination):
    cfg = json.loads((fixture / 'config.json').read_text())
    require(cfg.get('synthetic_only') is True, 'Requires an explicit synthetic fixture.')
    manifest = recovery.read_tsv(cfg['alignment_manifest'])
    require(len(manifest) == 1, 'Expected one synthetic sample.')
    original = json.loads(Path(manifest[0]['receipt']).read_text())
    require(original.get('synthetic_fixture_not_patient_data') is True, 'Not a synthetic source receipt.')
    sid = manifest[0]['sample_id']
    source_root = Path(cfg['output_root']) / 'checkpoints' / sid
    source = source_root / 'calls/contamination.filtered.vcf'
    bam = source_root / 'alignment-normal/marked.bam'
    require(source.is_file() and bam.is_file(), 'Complete synthetic MT fixture required.')
    destination.mkdir(mode=0o700)
    cfg = dict(cfg, output_root=str(destination / 'output'))
    runtime = recovery.RecoveryRuntime(cfg, True, destination, 'mitochondria')
    resources = recovery.mt.validate_resources(cfg, resource_path, runtime)
    headers, rows = recovery.read_vcf(source)
    require({int(row[1]) for row in rows} == {100, 5000}, 'Unexpected synthetic candidate fixture.')
    for row in rows:
        if int(row[1]) == 5000:
            row[6] = 'weak_evidence'
    # Policy-only end-boundary row: explicitly fabricated, not a Mutect2 call.
    end = list(rows[0])
    end[1] = '16560'
    rows.append(end)
    test_input, output = destination / 'synthetic-input.vcf', destination / 'synthetic-output.vcf'
    test_input.write_text(''.join(headers) + ''.join('\t'.join(r) + '\n' for r in rows))
    recovery.fixed_artifact_check(runtime, resources, destination,
        ['FilterAlignmentArtifacts', '-V', test_input, '-I', bam, '-O', output])
    _, result = recovery.read_vcf(output)
    by_position = {int(r[1]): r for r in result}
    require(len(result) == 3 and set(by_position) == {100, 5000, 16560}, 'Candidate rows were lost.')
    for row in rows:
        got = by_position[int(row[1])]
        require(got[8:] == row[8:], 'Synthetic FORMAT/genotypes changed.')
    require(by_position[5000][6] == 'weak_evidence', 'Original soft filter lost or false alignment artifact.')
    require('NUMT_ASSESSMENT=ASSESSED_NORMAL_INTERIOR' in by_position[5000][7], 'Filtered interior row not assessed.')
    require(all('NUMTArtifactCheckUnavailable' in by_position[p][6] for p in (100, 16560)), 'Boundary falsely passed.')
    report = {'status': 'PASS', 'synthetic_only': True, 'patient_reads_used': False,
              'already_filtered_interior_row_assessed': True, 'original_filters_retained': True,
              'three_rows_and_FORMAT_preserved': True, 'boundary_100_and_16560_explicitly_unavailable': True,
              'end_boundary_16560_is_fabricated_policy_fixture': True, 'clinical_LoD_validated': False}
    json_write(destination / 'validation-summary.json', report)
    print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', required=True, type=Path)
    parser.add_argument('--resources', required=True, type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    require(args.execute, 'Explicit --execute required for synthetic regression.')
    os.umask(0o077)
    output = args.output or args.fixture / ('recovery-regression-' + dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    validate(args.fixture.resolve(), args.resources.resolve(), output.resolve())


if __name__ == '__main__':
    main()
