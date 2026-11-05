#!/usr/bin/env python3
"""Read-only technical MT output audit; no clinical variant interpretation."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import subprocess

from calling_stages import IDENTIFIER, read_tsv, require, sha256
from workflow_common import json_write


def query(argv):
    return subprocess.check_output(argv, text=True, timeout=120).splitlines()


def audit(output_root, samples_manifest, report, bcftools='bcftools'):
    ids = [r['sample_id'] for r in read_tsv(samples_manifest)]
    require(ids and len(set(ids)) == len(ids) and all(IDENTIFIER.fullmatch(s) for s in ids), 'Invalid sample manifest.')
    results = []
    for ordinal, sid in enumerate(ids, 1):
        directory = output_root / 'checkpoints' / sid / 'calls'
        receipt_path = directory / 'workflow-receipt.json'
        receipt = json.loads(receipt_path.read_text())
        require(receipt['status'] == 'PASS', 'A sample has no completed technical checkpoint.')
        for name, expected in receipt['artifacts'].items():
            path = (directory / name).resolve(strict=True)
            require(path.is_relative_to(directory.resolve()) and path.stat().st_size == expected['bytes']
                    and sha256(path) == expected['sha256'], 'Output artifact SHA256/size changed.')
        files = {label: directory / ('mitochondria.' + label + '.vcf.gz') for label in ('raw', 'filtered', 'pass')}
        record_counts = {}
        for label, path in files.items():
            require(query([bcftools, 'query', '-l', str(path)]) == [sid], 'VCF sample header mismatch.')
            count = int(query([bcftools, 'index', '-n', str(path)])[0])
            indexed = query([bcftools, 'view', '-H', '-r', 'chrM', str(path)])
            require(len(indexed) == count and all(r.startswith('chrM\t') for r in indexed), 'MT VCF index/contig mismatch.')
            record_counts[label] = count
        fmt = '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD\t%AF]\n'
        data = {}
        for label in ('raw', 'filtered', 'pass'):
            rows = [line.split('\t') for line in query([bcftools, 'query', '-f', fmt, str(files[label])])]
            mapping = {tuple(row[:4]): row[4:] for row in rows}
            require(len(mapping) == len(rows) == record_counts[label], 'Duplicate/lost MT allele row.')
            data[label] = mapping
        require(data['raw'] == data['filtered'], 'Raw and filtered allele/GT/DP/AD/AF changed.')
        require(all(data['filtered'].get(k) == v for k, v in data['pass'].items()), 'PASS subset changed genotypes.')
        states = query([bcftools, 'query', '-f', '%FILTER\t%INFO/NUMT_ASSESSMENT\n', str(files['filtered'])])
        counts = {'assessed': 0, 'unavailable': 0, 'alignment_artifact': 0, 'rows': len(states)}
        for line in states:
            flag, state = line.split('\t')
            flags = flag.split(';')
            require(state in ('ASSESSED_NORMAL_INTERIOR', 'UNAVAILABLE_CIRCULAR_BOUNDARY'), 'Missing NUMT assessment status.')
            if state == 'UNAVAILABLE_CIRCULAR_BOUNDARY':
                require('NUMTArtifactCheckUnavailable' in flags, 'Boundary incorrectly accepted as NUMT-assessed.')
                counts['unavailable'] += 1
            else:
                counts['assessed'] += 1
            counts['alignment_artifact'] += 'alignment' in flags
        recorded = json.loads((directory / 'numt-assessment.json').read_text())
        require(recorded['counts'] == counts, 'NUMT summary differs from VCF.')
        pass_states = query([bcftools, 'query', '-f', '%FILTER\t%INFO/NUMT_ASSESSMENT\n', str(files['pass'])])
        require(all(line == 'PASS\tASSESSED_NORMAL_INTERIOR' for line in pass_states), 'Unassessed/filtered row in PASS subset.')
        research = json.loads((directory / 'research-report.json').read_text())
        results.append({'sample_ordinal': ordinal, 'calls_directory': str(directory),
                        'receipt_sha256': sha256(receipt_path), 'artifact_hashes_verified': len(receipt['artifacts']),
                        'record_counts': record_counts, 'NUMT_assessment': counts,
                        'coverage': research['coverage'], 'contamination_status': research['contamination']['status'],
                        'raw_filtered_GT_DP_AD_AF_preserved': True, 'indices_and_single_sample_header': 'PASS'})
    value = {'status': 'PASS_RESEARCH_WITH_LIMITATIONS', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
             'samples': len(results), 'samples_manifest_sha256': sha256(samples_manifest),
             'audit_code_sha256': sha256(Path(__file__)), 'results': results,
             'clinical_LoD_validated': False, 'record_counts_are_not_pathogenic_variant_counts': True}
    require(not report.exists(), 'Audit report already exists; choose a new report path.')
    json_write(report, value)
    print(json.dumps({'status': value['status'], 'samples': len(results), 'artifact_hashes': 'PASS',
                      'GT_DP_AD_AF_preservation': 'PASS', 'indices': 'PASS', 'clinical_LoD_validated': False}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--samples-manifest', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--bcftools', default='bcftools')
    args = parser.parse_args()
    os.umask(0o077)
    audit(args.output_root.resolve(), args.samples_manifest.resolve(), args.report.resolve(), args.bcftools)


if __name__ == '__main__':
    main()
