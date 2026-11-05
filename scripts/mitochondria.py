#!/usr/bin/env python3
"""Research mtDNA Mutect2 workflow derived from the pinned GATK4.3 mitochondrial WDL.

Normal + 8000-base shifted alignment/calling, liftback, contamination and nuclear
alignment-artifact checks. Original whole-genome BAMs are never modified.
"""
import argparse
import concurrent.futures
import csv
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

from calling_stages import IDENTIFIER, mounts, read_tsv, require, sha256
from workflow_common import Runtime, json_write, safe_key

NORMAL_INTERVAL = 'chrM:576-16024'
SHIFTED_INTERVAL = 'chrM:8025-9144'


class MitoRuntime(Runtime):
    def code_identity(self):
        return {name: self.file_hash(self.code_root / 'scripts' / name)
                for name in ('mitochondria.py', 'workflow_common.py')}

    def child(self, directory, stage=None):
        child = MitoRuntime(self.cfg, self.execute, directory, 'mitochondria', self.started)
        child._hashes, child._containers, child._docker = self._hashes, self._containers, self._docker
        return child

    def checkpoint(self, key, inputs, parameters, producer):
        target = self.output_root / 'checkpoints' / safe_key(key)
        return self._checkpoint_at(target, key, inputs, parameters, producer)


def docker_at(runtime, image, arguments, mapping, directory, cpus=4, memory=16):
    command = runtime.docker(image, arguments, mounts=mapping, gpu=False, cpus=cpus, memory_gib=memory)
    position = command.index(image)
    return command[:position] + ['--workdir', str(directory)] + command[position:]


def gatk(runtime, resources, directory, arguments, readonly=(), heap=6, log='gatk'):
    mapping = mounts(*((Path(p), True) for p in readonly), (Path(directory), False))
    command = ['gatk', '--java-options', f'-Xmx{heap}g -XX:ParallelGCThreads=2', *map(str, arguments)]
    return runtime.run(docker_at(runtime, resources['gatk_image'], command, mapping, directory), log)


def validate_resources(cfg, path, runtime):
    resources = json.loads(Path(path).read_text())
    require(resources.get('status') == 'PASS' and resources['whole_reference_sha256'] == cfg['reference']['sha256'],
            'Mitochondrial resources are not PASS for this exact whole-genome reference.')
    require(resources['normal_sequence_md5'] == 'c68f52674c9fb33aef52dcf399755519' and resources['shift_bases'] == 8000,
            'Expected rCRS16569 and an exact8000-base rotation.')
    for name, item in resources['files'].items():
        p = Path(name).resolve(strict=True)
        require(p.stat().st_size == item['bytes'] and runtime.file_hash(p) == item['sha256'], 'Mitochondrial resource hash changed.')
    return resources


def validate_alignments(cfg, runtime):
    samples = read_tsv(cfg['samples_manifest'])
    ids = [r['sample_id'] for r in samples]
    require(ids and len(set(ids)) == len(ids) and all(IDENTIFIER.fullmatch(s) for s in ids), 'Invalid source cohort.')
    rows = read_tsv(cfg['alignment_manifest'])
    require(len(rows) == len(ids) and {r['sample_id'] for r in rows} == set(ids), 'Alignment cohort mismatch.')
    ordered = [next(r for r in rows if r['sample_id'] == s) for s in ids]
    for row in ordered:
        bam, receipt_path = Path(row['bam']).resolve(strict=True), Path(row['receipt']).resolve(strict=True)
        require(bam.suffix == '.bam' and bam.parent == receipt_path.parent, 'Only the original receipt-associated BAM is accepted.')
        receipt = json.loads(receipt_path.read_text())
        require(receipt['sample'] == row['sample_id']
                and receipt['config']['reference']['sha256'] == cfg['reference']['sha256'], 'Original BAM provenance mismatch.')
        for p in (bam, Path(str(bam) + '.bai')):
            require(runtime.file_hash(p) == receipt['outputs'][p.name], 'Original BAM/index differs from its source receipt.')
        row.update(bam=str(bam), receipt=str(receipt_path), bam_sha256=runtime.file_hash(bam))
    return ordered


def parse_contamination(path):
    with Path(path).open() as stream:
        lines = [line.replace('"', '') for line in stream if line.strip()]
    reader = csv.DictReader(lines, delimiter='\t')
    required = ('SampleID', 'Contamination', 'HgMajor', 'HgMinor', 'MeanHetLevelMajor', 'MeanHetLevelMinor')
    require(reader.fieldnames and all(k in reader.fieldnames for k in required), 'Unexpected haplochecker output schema.')
    rows = list(reader)
    require(len(rows) == 1, 'Haplochecker must report exactly one isolated input sample.')
    row = rows[0]
    status = row['Contamination'].strip().upper()
    estimate = None
    if status == 'NO':
        estimate = 0.0
    elif status == 'YES':
        try:
            major, minor = float(row['MeanHetLevelMajor']), float(row['MeanHetLevelMinor'])
            if all(math.isfinite(v) and 0 <= v <= 1 for v in (major, minor)):
                estimate = minor if major == 0 else 1 - major
        except ValueError:
            pass
    return {'status': status, 'estimate': estimate, 'major_haplogroup': row['HgMajor'],
            'minor_haplogroup': row['HgMinor'], 'raw_report': row,
            'method': 'Pinned Broad GATK4.3 WDL haplochecker image; no nuclear VerifyBamID supplied',
            'clinical_contamination_detection_limit_validated': False}


def add_unknown_contamination(source, destination):
    with Path(source).open() as inp, Path(destination).open('x') as out:
        for line in inp:
            if line.startswith('#CHROM'):
                out.write('##FILTER=<ID=ContaminationEstimateUnavailable,Description="Mitochondrial contamination could not be estimated; all candidates retained for research review">\n')
            if line.startswith('#'):
                out.write(line)
            else:
                row = line.rstrip('\n').split('\t')
                flags = [] if row[6] in ('.', 'PASS') else row[6].split(';')
                row[6] = ';'.join(dict.fromkeys(flags + ['ContaminationEstimateUnavailable']))
                out.write('\t'.join(row) + '\n')


def selection(runtime, cfg, resources, sample):
    sid, bam = sample['sample_id'], Path(sample['bam'])
    inputs = [bam, Path(str(bam) + '.bai'), sample['receipt'], resources['selection_intervals']]
    parameters = {'source_bam_sha256': sample['bam_sha256'], 'source_receipt_sha256': sha256(sample['receipt']),
                  'selection': 'Official chrM+NUMT intervals seed query names, then recover ALL primary mates from the full BAM',
                  'primary_exclude_flags': '0x900', 'drop_duplicate_reads': False,
                  'resource_manifest_sha256': cfg['mt_resource_manifest_sha256']}
    def produce(directory):
        directory = Path(directory)
        bed = directory / 'seed-regions.bed'
        intervals = []
        for line in Path(resources['selection_intervals']).read_text().splitlines():
            if line and not line.startswith('@'):
                c, start, end = line.split('\t')[:3]
                intervals.append((c, int(start) - 1, int(end)))
        require(any(c == 'chrM' and a == 0 and b == 16569 for c, a, b in intervals), 'Official seed intervals lack full chrM.')
        bed.write_text(''.join(f'{c}\t{a}\t{b}\n' for c, a, b in intervals))
        seed = directory / 'seed.sam'
        runtime.run(['samtools', 'view', '-M', '-F', '2304', '-L', str(bed), str(bam)], 'seed-selection', stdout=seed)
        names, nuclear_seeds, mt_seeds = set(), 0, 0
        with seed.open() as stream:
            for line in stream:
                row = line.split('\t')
                names.add(row[0])
                if row[2] == 'chrM': mt_seeds += 1
                else: nuclear_seeds += 1
        require(names, 'No chrM/NUMT seed reads; mtDNA is not assessed for this sample.')
        namefile = directory / 'read-names.txt'
        namefile.write_text('\n'.join(sorted(names)) + '\n')
        selected = directory / 'selected.bam'
        runtime.run(['samtools', 'view', '-b', '-F', '2304', '-N', str(namefile), '-o', str(selected), str(bam)], 'recover-all-mates')
        records = directory / 'selected.sam'
        runtime.run(['samtools', 'view', str(selected)], 'pair-check', stdout=records)
        counts, hardclips = {}, 0
        with records.open() as stream:
            for line in stream:
                row = line.split('\t'); flag = int(row[1])
                counts.setdefault(row[0], []).append(flag)
                hardclips += 'H' in row[5]
        require(set(counts) == names and all(len(v) == 2 and sum(bool(f & 64) for f in v) == 1
                and sum(bool(f & 128) for f in v) == 1 for v in counts.values()),
                'Selected query names do not have complete primary R1/R2 pairs; refusing silent mate loss.')
        runtime.run(['samtools', 'quickcheck', '-v', str(selected)], 'selected-bam-check')
        unmapped = directory / 'unmapped.bam'
        gatk(runtime, resources, directory, ['RevertSam', '-I', selected, '-O', unmapped,
             '--OUTPUT_BY_READGROUP', 'false', '--VALIDATION_STRINGENCY', 'LENIENT',
             '--ATTRIBUTE_TO_CLEAR', 'FT', '--ATTRIBUTE_TO_CLEAR', 'CO', '--SORT_ORDER', 'queryname',
             '--RESTORE_ORIGINAL_QUALITIES', 'false'], readonly=[directory], heap=2, log='revert')
        metadata = dict(parameters, selected_pairs=len(names), selected_primary_records=sum(map(len, counts.values())),
                        mt_seed_records=mt_seeds, nuclear_NUMT_seed_records=nuclear_seeds,
                        selected_hardclipped_records=hardclips, original_source_bam_modified=False)
        json_write(directory / 'selection-qc.json', metadata)
        return {'outputs': ['selected.bam', 'unmapped.bam', 'selection-qc.json'], 'metadata': metadata}
    return Path(runtime.checkpoint(sid + '/selection', inputs, parameters, produce))


def align(runtime, resources, sample, selected, branch):
    reference = Path(resources['normal_fasta' if branch == 'normal' else 'shifted_fasta'])
    inputs = [selected / 'unmapped.bam', reference, *[Path(str(reference) + s) for s in ('.amb', '.ann', '.bwt', '.pac', '.sa')]]
    params = {'branch': branch, 'reference_sha256': sha256(reference), 'bwa_threads': 2,
              'bwa': '0.7.17', 'duplicate_reads': 'marked, not removed', 'optical_duplicate_detection': 'disabled; original instrument provenance unknown'}
    def produce(directory):
        directory = Path(directory)
        fastq, aligned, merged, marked = (directory / n for n in ('interleaved.fastq', 'aligned.sam', 'merged.bam', 'duplicate-marked.bam'))
        gatk(runtime, resources, directory, ['SamToFastq', '-I', selected / 'unmapped.bam',
             '--FASTQ', fastq, '--INTERLEAVE', 'true', '--NON_PF', 'true'], readonly=[selected], heap=2, log='sam-to-fastq')
        runtime.run(['bwa', 'mem', '-K', '100000000', '-p', '-v', '3', '-t', '2', '-Y', str(reference), str(fastq)], 'bwa', stdout=aligned)
        gatk(runtime, resources, directory, ['MergeBamAlignment', '--ALIGNED_BAM', aligned,
             '--UNMAPPED_BAM', selected / 'unmapped.bam', '-O', merged, '-R', reference,
             '--VALIDATION_STRINGENCY', 'SILENT', '--EXPECTED_ORIENTATIONS', 'FR',
             '--ATTRIBUTES_TO_RETAIN', 'X0', '--ATTRIBUTES_TO_REMOVE', 'NM', '--ATTRIBUTES_TO_REMOVE', 'MD',
             '--PAIRED_RUN', 'true', '--SORT_ORDER', 'queryname', '--CLIP_ADAPTERS', 'false',
             '--ADD_MATE_CIGAR', 'true', '--MAX_INSERTIONS_OR_DELETIONS', '-1',
             '--PRIMARY_ALIGNMENT_STRATEGY', 'MostDistant', '--UNMAPPED_READ_STRATEGY', 'COPY_TO_TAG',
             '--ALIGNER_PROPER_PAIR_FLAGS', 'true', '--UNMAP_CONTAMINANT_READS', 'true',
             '--ADD_PG_TAG_TO_READS', 'false'], readonly=[selected, reference.parent], heap=4, log='merge-alignments')
        gatk(runtime, resources, directory, ['MarkDuplicates', '-I', merged, '-O', marked,
             '--METRICS_FILE', directory / 'duplicate-metrics.txt', '--ASSUME_SORT_ORDER', 'queryname',
             '--READ_NAME_REGEX', 'null', '--VALIDATION_STRINGENCY', 'SILENT', '--ADD_PG_TAG_TO_READS', 'false'],
             readonly=[], heap=4, log='mark-duplicates')
        target = directory / 'marked.bam'
        runtime.run(['samtools', 'sort', '-@', '2', '-o', str(target), str(marked)], 'sort')
        runtime.run(['samtools', 'index', '-@', '2', str(target)], 'index')
        runtime.run(['samtools', 'quickcheck', '-v', str(target)], 'bam-check')
        runtime.run(['samtools', 'depth', '-aa', '-s', '-q', '10', '-Q', '20', '-G', '3844', '-r', 'chrM', str(target)],
                    'depth', stdout=directory / 'depth.tsv')
        return {'outputs': ['marked.bam', 'marked.bam.bai', 'duplicate-metrics.txt', 'depth.tsv'], 'metadata': params}
    return Path(runtime.checkpoint(sample['sample_id'] + '/alignment-' + branch, inputs, params, produce))


def merge_depth(normal, shifted, destination):
    def read(path):
        return {int(row[1]): int(row[2]) for row in csv.reader(Path(path).open(), delimiter='\t')}
    a, b = read(normal), read(shifted)
    require(set(a) == set(b) == set(range(1, 16570)), 'Expected all16569 mitochondrial coverage positions.')
    depths = []
    with Path(destination).open('x') as out:
        out.write('position\tdepth\treference_branch\n')
        for position in range(1, 16570):
            use_normal = 576 <= position <= 16024
            shifted_position = (position - 8001) % 16569 + 1
            depth = a[position] if use_normal else b[shifted_position]
            depths.append(depth)
            out.write(f'{position}\t{depth}\t' + ('normal' if use_normal else 'shifted-liftback') + '\n')
    return {'positions': 16569, 'mean_depth': sum(depths) / 16569,
            'breadth_pct': {str(t): 100 * sum(d >= t for d in depths) / 16569 for t in (1, 10, 20, 50, 100)},
            'zero_pct': 100 * depths.count(0) / 16569, 'min_base_quality': 10, 'min_mapping_quality': 20,
            'mate_overlap_counted_once': True, 'excluded_flags': '0xF04', 'clinical_LoD_or_callability_validated': False}


def final_calls(runtime, cfg, resources, sample, normal, shifted):
    normal_ref, shifted_ref = Path(resources['normal_fasta']), Path(resources['shifted_fasta'])
    inputs = [normal / 'marked.bam', normal / 'marked.bam.bai', shifted / 'marked.bam', shifted / 'marked.bam.bai',
              normal_ref, shifted_ref, resources['shift_back_chain'], resources['blacklist_bed'], resources['whole_genome_image']]
    params = {'caller': 'GATK4.3 Mutect2 mitochondria-mode', 'normal_interval': NORMAL_INTERVAL,
              'shifted_interval': SHIFTED_INTERVAL, 'max_reads_per_alignment_start': 75, 'min_allele_fraction_filter': 0,
              'whole_genome_NUMT_artifact_filter': True, 'resource_manifest_sha256': cfg['mt_resource_manifest_sha256']}
    def produce(directory):
        directory = Path(directory)
        refs = [normal_ref.parent, Path(resources['whole_genome_image']).parent]
        for name, aligned, ref, interval in [('normal', normal, normal_ref, NORMAL_INTERVAL), ('shifted', shifted, shifted_ref, SHIFTED_INTERVAL)]:
            gatk(runtime, resources, directory, ['Mutect2', '-R', ref, '-I', aligned / 'marked.bam',
                 '-O', directory / (name + '.raw.vcf'), '-L', interval,
                 '--read-filter', 'MateOnSameContigOrNoMappedMateReadFilter', '--read-filter', 'MateUnmappedAndUnmappedReadFilter',
                 '--annotation', 'StrandBiasBySample', '--mitochondria-mode', '--max-reads-per-alignment-start', '75',
                 '--max-mnp-distance', '0', '--native-pair-hmm-threads', '2'], readonly=refs + [aligned], heap=8, log='mutect-' + name)
        lifted, rejected, combined = (directory / n for n in ('shifted.lifted.vcf', 'shifted.rejected.vcf', 'mitochondria.raw.vcf'))
        gatk(runtime, resources, directory, ['LiftoverVcf', '-I', directory / 'shifted.raw.vcf', '-O', lifted,
             '-R', normal_ref, '--CHAIN', resources['shift_back_chain'], '--REJECT', rejected], readonly=refs, heap=2, log='liftback')
        rejected_count = sum(not line.startswith('#') for line in rejected.open())
        require(rejected_count == 0, 'Mitochondrial liftback rejected candidates; preserved for explicit boundary review.')
        gatk(runtime, resources, directory, ['MergeVcfs', '-I', directory / 'normal.raw.vcf', '-I', lifted, '-O', combined],
             readonly=refs, heap=2, log='combine-regions')
        stats = directory / 'combined.stats'
        gatk(runtime, resources, directory, ['MergeMutectStats', '--stats', directory / 'normal.raw.vcf.stats',
             '--stats', directory / 'shifted.raw.vcf.stats', '-O', stats], heap=2, log='merge-stats')
        def filtering(source, target, estimate=None, label='initial'):
            arguments = ['FilterMutectCalls', '-V', source, '-R', normal_ref, '-O', directory / (label + '.pre-mask.vcf'),
                         '--stats', stats, '--mitochondria-mode', '--max-alt-allele-count', '4', '--min-allele-fraction', '0']
            if estimate is not None:
                arguments += ['--contamination-estimate', str(estimate)]
            gatk(runtime, resources, directory, arguments, readonly=refs, heap=4, log=label + '-filter')
            gatk(runtime, resources, directory, ['VariantFiltration', '-V', directory / (label + '.pre-mask.vcf'), '-O', target,
                 '--apply-allele-specific-filters', '--mask', resources['blacklist_bed'], '--mask-name', 'blacklisted_site'],
                 readonly=refs, heap=2, log=label + '-blacklist')
        initial = directory / 'initial.filtered.vcf'
        filtering(combined, initial)
        haplo_input, haplo_output = directory / 'haplo-input', directory / 'haplo-output'
        haplo_input.mkdir(); haplo_output.mkdir()
        split = directory / 'initial.split.vcf'
        gatk(runtime, resources, directory, ['LeftAlignAndTrimVariants', '-R', normal_ref, '-V', initial, '-O', split,
             '--split-multi-allelics', '--dont-trim-alleles', '--keep-original-ac'], readonly=refs, heap=2, log='split-haplo-input')
        gatk(runtime, resources, directory, ['SelectVariants', '-V', split, '-O', haplo_input / 'pass.vcf', '--exclude-filtered'],
             heap=2, log='select-haplo-input')
        pass_count = sum(not line.startswith('#') for line in (haplo_input / 'pass.vcf').open())
        if pass_count:
            command = docker_at(runtime, resources['haplochecker_image'], ['java', '-Xmx2g', '-jar',
                        '/usr/mtdnaserver/haplocheckCLI.jar', str(haplo_input)],
                        mounts((haplo_input, True), (haplo_output, False)), haplo_output, cpus=2, memory=4)
            runtime.run(command, 'haplochecker', timeout=600)
            contamination = parse_contamination(haplo_output / 'output')
        else:
            contamination = {'status': 'NO_PASS_VARIANTS', 'estimate': None, 'clinical_contamination_detection_limit_validated': False}
        json_write(directory / 'contamination.json', contamination)
        second = directory / 'contamination.filtered.vcf'
        filtering(initial, second, contamination['estimate'], 'contamination')
        artifact = directory / 'whole-genome-artifacts.vcf'
        gatk(runtime, resources, directory, ['FilterAlignmentArtifacts', '-R', normal_ref, '-V', second,
             '-I', normal / 'marked.bam', '--bwa-mem-index-image', resources['whole_genome_image'], '-O', artifact],
             readonly=refs + [normal], heap=6, log='whole-genome-NUMT-artifacts')
        final = artifact
        if contamination['estimate'] is None:
            final = directory / 'mitochondria.contamination-unknown.vcf'
            add_unknown_contamination(artifact, final)
        tool = cfg['tools'].get('bcftools', 'bcftools')
        for source, basename in [(combined, 'mitochondria.raw.vcf.gz'), (final, 'mitochondria.filtered.vcf.gz')]:
            target = directory / basename
            runtime.run([tool, 'view', '-Oz', '-o', str(target), str(source)], 'compress-' + basename)
            runtime.run([tool, 'index', '--force', '--tbi', str(target)], 'index-' + basename)
            actual = runtime.run([tool, 'query', '-l', str(target)], 'sample-check', capture=True).splitlines()
            require(actual == [sample['sample_id']], 'Mitochondrial output sample mismatch.')
        runtime.run([tool, 'view', '-f', 'PASS', '-Oz', '-o', str(directory / 'mitochondria.pass.vcf.gz'),
                     str(directory / 'mitochondria.filtered.vcf.gz')], 'pass-subset')
        runtime.run([tool, 'index', '--force', '--tbi', str(directory / 'mitochondria.pass.vcf.gz')], 'pass-index')
        table = runtime.run([tool, 'query', '-f', '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%DP\t%AD\t%AF]\n',
                            str(directory / 'mitochondria.filtered.vcf.gz')], 'heteroplasmy-table', capture=True)
        (directory / 'heteroplasmy.tsv').write_text('contig\tposition\tref\talt\tfilter\tFORMAT_DP\tFORMAT_AD\tFORMAT_AF\n' + table)
        coverage = merge_depth(normal / 'depth.tsv', shifted / 'depth.tsv', directory / 'circular-coverage.tsv')
        report = dict(params, contamination=contamination, coverage=coverage, liftover_rejections=rejected_count,
                      scope='Research mitochondrial small variants; no validated clinical sensitivity/LoD',
                      limits=['No mtDNA structural-deletion assay or tissue-specific inference',
                              'Reference NUMT alignment checking does not prove exclusion of every non-reference polymorphic NUMT',
                              'Boundary artifact checking uses the normal-reference BAM; raw shifted calls are retained',
                              'No personalized self-reference round is implemented; pinned GATK4.3 WDL normal/shift design'])
        json_write(directory / 'research-report.json', report)
        return {'outputs': ['mitochondria.raw.vcf.gz', 'mitochondria.raw.vcf.gz.tbi', 'mitochondria.filtered.vcf.gz',
                            'mitochondria.filtered.vcf.gz.tbi', 'mitochondria.pass.vcf.gz', 'mitochondria.pass.vcf.gz.tbi',
                            'heteroplasmy.tsv', 'research-report.json', 'contamination.json', 'circular-coverage.tsv'], 'metadata': report}
    return Path(runtime.checkpoint(sample['sample_id'] + '/calls', inputs, params, produce))


def run(cfg, resource_path, execute, workers=1):
    require(1 <= workers <= 4, 'Use1..4 mitochondrial sample workers.')
    runtime = MitoRuntime(cfg, execute, stage='mitochondria')
    resources = validate_resources(cfg, resource_path, runtime)
    rows = validate_alignments(cfg, runtime)
    runtime.remaining_seconds(); runtime.check_space()
    if not execute:
        print(json.dumps({'status': 'PLAN', 'samples': len(rows), 'workers': workers, 'CPU_per_worker': 4,
                          'GiB_per_worker': 16, 'caller': 'Mutect2 mitochondria-mode', 'source_BAMs_modified': False}))
        return
    runtime.output_root.mkdir(parents=True, exist_ok=True)
    state_path = Path(cfg['data_root']) / 'private/mitochondria-status.json'
    json_write(state_path, {'status': 'RUNNING', 'samples': len(rows), 'workers': workers, 'pid': os.getpid(),
                           'started_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'output_root': str(runtime.output_root)})
    def one(sample):
        local = runtime.child(runtime.output_root)
        selected = selection(local, cfg, resources, sample)
        normal = align(local, resources, sample, selected, 'normal')
        shifted = align(local, resources, sample, selected, 'shifted')
        return str(final_calls(local, cfg, resources, sample, normal, shifted))
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            completed = list(pool.map(one, rows))
        json_write(state_path, {'status': 'PASS', 'samples': len(rows), 'completed': completed,
                               'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                               'clinical_classification_complete': False})
        print(json.dumps({'status': 'PASS', 'samples': len(completed), 'clinical_LoD_validated': False}))
    except BaseException as error:
        json_write(state_path, {'status': 'FAIL', 'error_type': type(error).__name__, 'updated_utc': dt.datetime.now(dt.timezone.utc).isoformat()})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-config', type=Path, required=True)
    parser.add_argument('--resources', type=Path, required=True)
    parser.add_argument('--label', default='v1')
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args(); os.umask(0o077)
    require(IDENTIFIER.fullmatch(args.label), 'Unsafe mitochondrial output label.')
    original = json.loads(args.workflow_config.read_text())
    cfg = dict(original)
    cfg.update(output_root=str(Path(original['data_root']) / 'outputs' / ('mitochondria-' + args.label)),
               mt_resource_manifest_sha256=sha256(args.resources), source_workflow_config_sha256=sha256(args.workflow_config),
               scope={'label': 'Research mtDNA small-variant/heteroplasmy workflow', 'reference': 'rCRS chrM16569',
                      'circular_control_region': '8000-base shifted alignment/calling and liftback',
                      'not_assessed': ['clinical LoD', 'large mtDNA structural deletions', 'tissue-specific clinical interpretation']})
    if args.execute:
        private = Path(cfg['data_root']) / 'private'
        private.mkdir(parents=True, exist_ok=True)
        with (private / '.mitochondria.lock').open('a') as guard:
            fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
            run(cfg, args.resources.resolve(), True, args.workers)
    else:
        run(cfg, args.resources.resolve(), False, args.workers)


if __name__ == '__main__':
    main()
