#!/usr/bin/env python3
"""Versioned research MT recovery; preserve original MT checkpoints unchanged.

GATK4.3 FAA pads beyond chrM end and compares dictionary numeric indices.
Use the full reference dictionary for FAA and explicitly withhold boundary
NUMT assessment. This is an incomplete research assessment, not clinical QC.
"""
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path

import mitochondria as mt
from calling_stages import read_tsv, require, sha256
from workflow_common import json_write

ORIGINAL_GATK = mt.gatk
POLICY = {
    'version': 2,
    'already_filtered_variants': 'Assess with --dont-skip-filtered-variants true; preserve original filters',
    'normal_interior': 'chrM:576-16024',
    'boundary_policy': 'Retain shifted/lifted candidates, mark NUMTArtifactCheckUnavailable',
    'dictionary_policy': 'Textual SAM conversion to full reference SQ dictionary; identical chrM sequence/coordinates',
    'gatk_source': 'https://github.com/broadinstitute/gatk/blob/4.3.0.0/src/main/java/org/broadinstitute/hellbender/tools/walkers/realignmentfilter/FilterAlignmentArtifacts.java',
    'known_issue': 'Unclipped lastReadEnd+50; full image reference IDs require matching -R dictionary order',
    'clinical_LoD_validated': False,
}


def read_vcf(path):
    headers, rows = [], []
    for line in Path(path).read_text().splitlines(True):
        if line.startswith('#'):
            headers.append(line)
        else:
            rows.append(line.rstrip('\n').split('\t'))
    return headers, rows


def key(row):
    return tuple(row[i] for i in (0, 1, 3, 4))


def flags(value):
    return [] if value in ('.', 'PASS') else value.split(';')


def merge_assessment(original, assessed, output):
    """Preserve every original allele/genotype and original soft FILTER."""
    headers, rows = read_vcf(original)
    ah, ar = read_vcf(assessed)
    amap = {key(r): r for r in ar}
    require(len(amap) == len(ar), 'FAA returned duplicate allele keys.')
    require(set(amap).issubset({key(r) for r in rows}), 'FAA returned unknown candidates.')
    extra = [h for h in ah if h.startswith(('##INFO=', '##FILTER='))
             and h.split(',')[0] not in {x.split(',')[0] for x in headers}]
    extra += ['##INFO=<ID=NUMT_ASSESSMENT,Number=1,Type=String,Description="Research whole-reference alignment artifact assessment status">\n',
              '##FILTER=<ID=NUMTArtifactCheckUnavailable,Description="Whole-reference NUMT artifact check unavailable at circular boundary; retained for research review">\n']
    counts = {'assessed': 0, 'unavailable': 0, 'alignment_artifact': 0, 'rows': len(rows)}
    with Path(output).open('x') as out:
        for h in headers:
            if h.startswith('#CHROM'):
                out.writelines(extra)
            out.write(h)
        for source in rows:
            row = list(source)
            result = amap.get(key(row))
            state = 'ASSESSED_NORMAL_INTERIOR' if result else 'UNAVAILABLE_CIRCULAR_BOUNDARY'
            combined_flags = flags(row[6])
            if result:
                require(result[8:] == row[8:], 'FAA changed FORMAT/genotypes; refusing recovery.')
                combined_flags += flags(result[6])
                counts['assessed'] += 1
                counts['alignment_artifact'] += 'alignment' in flags(result[6])
            else:
                combined_flags += ['NUMTArtifactCheckUnavailable']
                counts['unavailable'] += 1
            row[6] = ';'.join(dict.fromkeys(combined_flags)) or 'PASS'
            row[7] = ('' if row[7] == '.' else row[7] + ';') + 'NUMT_ASSESSMENT=' + state
            out.write('\t'.join(row) + '\n')
    return counts


def fixed_artifact_check(runtime, resources, directory, args, readonly=(), heap=6, log='gatk'):
    if args[0] != 'FilterAlignmentArtifacts':
        return ORIGINAL_GATK(runtime, resources, directory, args, readonly, heap, log)
    directory = Path(directory)
    source = Path(args[args.index('-V') + 1])
    bam = Path(args[args.index('-I') + 1])
    output = Path(args[args.index('-O') + 1])
    whole = Path(runtime.cfg['reference']['fasta'])
    dictionary = whole.with_suffix('.dict')
    require(dictionary.is_file(), 'Missing full reference dictionary.')
    sam = directory / 'numt-source.sam'
    runtime.run(['samtools', 'view', '-h', str(bam)], 'numt-source-sam', stdout=sam)
    remapped = directory / 'numt-full-dictionary.sam'
    sq = [x for x in dictionary.read_text().splitlines(True) if x.startswith('@SQ')]
    require(any('SN:chrM\t' in x and 'LN:16569' in x for x in sq), 'Full chrM dictionary differs.')
    with sam.open() as inp, remapped.open('x') as out:
        wrote = False
        for line in inp:
            if line.startswith('@SQ'):
                if not wrote:
                    out.writelines(sq)
                    wrote = True
            else:
                out.write(line)
    require(wrote, 'Source BAM has no sequence dictionary.')
    fullbam = directory / 'numt-full-dictionary.bam'
    runtime.run(['samtools', 'view', '-b', '-o', str(fullbam), str(remapped)], 'numt-remap-numeric-tids')
    runtime.run(['samtools', 'index', str(fullbam)], 'numt-remapped-index')
    runtime.run(['samtools', 'quickcheck', '-v', str(fullbam)], 'numt-remapped-check')
    before = runtime.run(['samtools', 'view', str(bam)], 'numt-before-records', capture=True)
    after = runtime.run(['samtools', 'view', str(fullbam)], 'numt-after-records', capture=True)
    require(before == after, 'Full dictionary conversion changed alignment records.')
    headers, rows = read_vcf(source)
    safe = [r for r in rows if 576 <= int(r[1]) <= 16024 and int(r[1]) + len(r[3]) - 1 <= 16024]
    interior = directory / 'numt-interior.vcf'
    interior.write_text(''.join(headers) + ''.join('\t'.join(r) + '\n' for r in safe))
    assessed = directory / 'numt-interior.assessed.vcf'
    if safe:
        ORIGINAL_GATK(runtime, resources, directory,
                      ['FilterAlignmentArtifacts', '-R', whole, '-V', interior, '-I', fullbam,
                       '--dont-skip-filtered-variants', 'true',
                       '--bwa-mem-index-image', resources['whole_genome_image'], '-O', assessed],
                      readonly=[whole.parent, Path(resources['whole_genome_image']).parent], heap=heap, log=log)
        _, got = read_vcf(assessed)
        require({key(r) for r in got} == {key(r) for r in safe}, 'FAA omitted interior candidate(s).')
    else:
        assessed.write_text(''.join(headers))
    counts = merge_assessment(source, assessed, output)
    json_write(directory / 'numt-assessment.json', dict(POLICY, counts=counts,
               source_bam_sha256=sha256(bam), full_dictionary_bam_sha256=sha256(fullbam),
               whole_reference_sha256=runtime.cfg['reference']['sha256']))


class RecoveryRuntime(mt.MitoRuntime):
    def code_identity(self):
        return {name: self.file_hash(self.code_root / 'scripts' / name)
                for name in ('recover_mitochondria.py', 'mitochondria.py', 'workflow_common.py')}

    def child(self, directory, stage=None):
        child = RecoveryRuntime(self.cfg, self.execute, directory, 'mitochondria', self.started)
        child._hashes, child._containers, child._docker = self._hashes, self._containers, self._docker
        return child

    def checkpoint(self, keyname, inputs, parameters, producer):
        def wrapped(directory):
            result = producer(directory)
            p = Path(directory) / 'research-report.json'
            report = json.loads(p.read_text())
            report['NUMT_recovery'] = json.loads((Path(directory) / 'numt-assessment.json').read_text())
            report['whole_genome_NUMT_artifact_filter'] = 'partial: normal interior only; boundaries explicitly unavailable'
            report['limits'] = [x for x in report['limits'] if not x.startswith('Boundary artifact')]
            report['limits'].append('Circular boundary candidates lack validated NUMT artifact assessment; excluded from PASS subset by explicit soft flag')
            json_write(p, report)
            result['metadata'] = report
            result['outputs'] += ['numt-assessment.json', 'numt-full-dictionary.bam', 'numt-full-dictionary.bam.bai']
            return result
        return super().checkpoint(keyname, inputs, dict(parameters, recovery_policy=POLICY), wrapped)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-config', type=Path, required=True)
    parser.add_argument('--resources', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    os.umask(0o077)
    require(args.source_root.resolve() != args.output_root.resolve(), 'Recovery must use a new output namespace.')
    original = json.loads(args.workflow_config.read_text())
    cfg = dict(original, output_root=str(args.source_root.resolve()),
               mt_resource_manifest_sha256=sha256(args.resources), source_workflow_config_sha256=sha256(args.workflow_config),
               scope={'label': 'Research mtDNA small-variant/heteroplasmy workflow', 'reference': 'rCRS chrM16569',
                      'circular_control_region': '8000-base shifted alignment/calling and liftback',
                      'not_assessed': ['clinical LoD', 'large mtDNA structural deletions', 'tissue-specific clinical interpretation']})
    if original.get('synthetic_only'):
        cfg = dict(original, output_root=str(args.source_root.resolve()))
    source_runtime = mt.MitoRuntime(cfg, args.execute, stage='mitochondria')
    resources = mt.validate_resources(cfg, args.resources, source_runtime)
    rows = mt.validate_alignments(cfg, source_runtime)
    if not args.execute:
        print(json.dumps({'status': 'PLAN', 'samples': len(rows), 'new_namespace': True, 'policy': POLICY}))
        return
    state = Path(cfg['data_root']) / 'private/mitochondria-recovery-status.json'
    recovery_cfg = dict(cfg, output_root=str(args.output_root.resolve()))
    runtime = RecoveryRuntime(recovery_cfg, True, stage='mitochondria')
    mt.gatk = fixed_artifact_check
    done = []
    with (Path(cfg['data_root']) / 'private/.mitochondria.lock').open('a') as guard:
        fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            for i, sample in enumerate(rows):
                json_write(state, {'status': 'RUNNING', 'sample_index': i + 1, 'samples': len(rows),
                                  'completed': len(done), 'pid': os.getpid(), 'policy': POLICY})
                selected = mt.selection(source_runtime, cfg, resources, sample)
                normal = mt.align(source_runtime, resources, sample, selected, 'normal')
                shifted = mt.align(source_runtime, resources, sample, selected, 'shifted')
                done.append(str(mt.final_calls(runtime, recovery_cfg, resources, sample, normal, shifted)))
            json_write(state, {'status': 'PASS_RESEARCH_WITH_LIMITATIONS', 'samples': len(done), 'completed': done,
                              'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'policy': POLICY})
            print(json.dumps({'status': 'PASS_RESEARCH_WITH_LIMITATIONS', 'samples': len(done)}))
        except BaseException as error:
            json_write(state, {'status': 'FAIL', 'completed': len(done), 'error_type': type(error).__name__,
                              'updated_utc': dt.datetime.now(dt.timezone.utc).isoformat()})
            raise


if __name__ == '__main__':
    main()
