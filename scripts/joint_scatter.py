#!/usr/bin/env python3
"""Joint genotype completed full-autosome gVCFs by chromosome; CPU only."""
import argparse
import concurrent.futures
import fcntl
import json
import os
from pathlib import Path
import re
import subprocess

from calling_stages import checked_vcf, context, mounts, read_tsv, require, sha256
from workflow_common import Runtime


def source_plan(cfg, runtime):
    reference, samples, callsets, scope = context(cfg)
    require(len(callsets) == 1 and callsets[0].get('scope', 'autosomes') == 'autosomes',
            'Scatter currently supports one explicit diploid autosome callset.')
    callset = callsets[0]
    ids = [s['sample_id'] for s in samples]
    require(callset.get('samples', ids) == ids and all(
        str(callset.get('per_sample_ploidy', {}).get(s, callset.get('ploidy', 2))) == '2' for s in ids),
        'Scatter requires the entire cohort and compatible diploid source gVCFs.')
    from calling_stages import interval_rows
    lengths = {r[0]: int(r[1]) for r in (line.split('\t') for line in Path(str(reference) + '.fai').read_text().splitlines())}
    spans = interval_rows(callset['intervals'], lengths)
    chromosomes = [c for c in lengths if any(c == chrom for chrom, start, end in spans)]
    require(len(spans) == len(chromosomes) and all((c, 1, lengths[c]) in spans for c in chromosomes),
            'Chromosome scatter requires exactly one full interval per declared autosome.')
    require(all(re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', c) for c in chromosomes), 'Only autosomes may be scattered.')
    call_root = runtime.stage_path('call')
    call_receipt = call_root / 'workflow-receipt.json'
    saved = json.loads(call_receipt.read_text())
    require(saved.get('status') == 'PASS', 'Stage3 must complete successfully before CPU scatter.')
    for relative, item in saved['artifacts'].items():
        p = (call_root / relative).resolve(strict=True)
        require(p.is_relative_to(call_root.resolve()) and p.stat().st_size == item['bytes']
                and runtime.file_hash(p) == item['sha256'], 'Completed stage3 artifact changed.')
    manifest = call_root / 'gvcfs.tsv'
    rows = read_tsv(manifest)
    require(len(rows) == len(ids) and {r['sample_id'] for r in rows} == set(ids), 'Incomplete source gVCF cohort.')
    rows = [next(r for r in rows if r['sample_id'] == sid) for sid in ids]
    inputs = [reference, Path(str(reference) + '.fai'), reference.with_suffix('.dict'), manifest, call_receipt, callset['intervals']]
    sources = []
    for row in rows:
        require(row['callset'] == callset['name'] and row['ploidy'] == '2'
                and row['intervals_sha256'] == callset['intervals_sha256']
                and row['reference_sha256'] == cfg['reference']['sha256'], 'Source gVCF full-interval identity mismatch.')
        gvcf = Path(row['gvcf']).resolve(strict=True)
        index = Path(str(gvcf) + '.tbi')
        receipt = gvcf.parent / 'workflow-receipt.json'
        original = json.loads(receipt.read_text())
        require(original.get('status') == 'PASS', 'Source gVCF checkpoint has no PASS receipt.')
        for p in (gvcf, index):
            item = original['artifacts'][p.name]
            require(p.stat().st_size == item['bytes'] and runtime.file_hash(p) == item['sha256'], 'Source gVCF/index hash changed.')
        actual = runtime.run([cfg['tools'].get('bcftools', 'bcftools'), 'query', '-l', str(gvcf)],
                             'source-sample-check', capture=True).splitlines()
        if not runtime.dry_run:
            require(actual == [row['sample_id']], 'Source sample identity mismatch; no renaming is permitted.')
        inputs += [gvcf, index, receipt]
        sources.append({'sample_id': row['sample_id'], 'gvcf': str(gvcf), 'sha256': runtime.file_hash(gvcf),
                        'source_receipt_sha256': runtime.file_hash(receipt),
                        'source_full_intervals_sha256': row['intervals_sha256']})
    return reference, ids, rows, chromosomes, scope, inputs, sources


def scatter(cfg, workdir, runtime, workers=4, cpus=6, memory=24, heap=16):
    require(1 <= workers <= 4 and 1 <= cpus and workers * cpus <= 24
            and 4 <= heap < memory and workers * memory <= 96, 'Scatter exceeds the 24CPU/96GiB allocation.')
    reference, ids, rows, chromosomes, scope, inputs, sources = source_plan(cfg, runtime)
    # Populate the shared hash cache before threads start; avoid re-reading each full gVCF per chromosome.
    runtime.input_identity(inputs)
    code_hash = sha256(Path(__file__).resolve())
    parameters = {'scatter_code_sha256': code_hash, 'workers': workers, 'cpus_per_job': cpus,
                  'memory_gib_per_job': memory, 'java_heap_gib_per_job': heap,
                  'reference_sha256': cfg['reference']['sha256'], 'source_full_gvcfs': sources,
                  'gatk_image': cfg['tools']['gatk_image'], 'full_source_manifest_sha256': sha256(runtime.stage_path('call') / 'gvcfs.tsv'),
                  'source_intervals_unchanged': True, 'subset_method': 'GATK -L chromosome only',
                  'reference_order': chromosomes, 'no_missing_to_reference': True}
    def job(chrom):
        local = runtime.child(workdir, 'joint')
        job_parameters = dict(parameters, chromosome=chrom)
        def produce(directory):
            directory = Path(directory)
            (directory / 'tmp').mkdir()
            sample_map = directory / 'sample-map.tsv'
            sample_map.write_text(''.join(r['sample_id'] + '\t' + r['gvcf'] + '\n' for r in rows))
            mapping = mounts((reference.parent, True), *((Path(r['gvcf']).parent, True) for r in rows), (directory, False))
            image = cfg['tools']['gatk_image']
            version = local.run(local.docker(image, ['gatk', '--version'], mounts=mapping,
                gpu=False, cpus=1, memory_gib=4), 'gatk-version', capture=True)
            require(re.search(r'(?<![0-9])4\.3\.0\.0(?![0-9])', version), 'GATK must be version 4.3.0.0.')
            database = directory / 'genomicsdb'
            target = directory / 'chromosome.raw.vcf.gz'
            java = f'-Xmx{heap}g -XX:ParallelGCThreads={min(cpus, 4)}'
            commands = [
                ['gatk', '--java-options', java, 'GenomicsDBImport', '--genomicsdb-workspace-path', str(database),
                 '--sample-name-map', str(sample_map), '--reader-threads', str(min(cpus, 4)), '--batch-size', '50',
                 '-L', chrom, '--tmp-dir', str(directory / 'tmp')],
                ['gatk', '--java-options', java, 'GenotypeGVCFs', '-R', str(reference),
                 '-V', 'gendb://' + str(database), '-O', str(target), '-L', chrom, '--sample-ploidy', '2',
                 '--standard-min-confidence-threshold-for-calling', '0', '--tmp-dir', str(directory / 'tmp')],
            ]
            for i, command in enumerate(commands):
                local.run(local.docker(image, command, mounts=mapping, gpu=False, cpus=cpus, memory_gib=memory), 'gatk-' + str(i))
            tool = cfg['tools'].get('bcftools', 'bcftools')
            ordered = directory / 'ordered.raw.vcf.gz'
            local.run([tool, 'view', '-s', ','.join(ids), '-Oz', '-o', str(ordered), str(target)], 'order-samples')
            checked_vcf(local, cfg, ordered, ids, 'chromosome')
            stats = local.run([tool, 'index', '--stats', str(ordered)], 'chromosome-contigs', capture=True)
            require(all(line.split('\t')[0] == chrom for line in stats.splitlines()), 'A scatter result escaped its chromosome subset.')
            return {'outputs': ['ordered.raw.vcf.gz', 'ordered.raw.vcf.gz.tbi'], 'metadata': job_parameters}
        return Path(local.checkpoint('scatter/' + chrom, inputs, job_parameters, produce)) / 'ordered.raw.vcf.gz'
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        # executor.map returns in input/reference order even when jobs finish in another order.
        files = list(pool.map(job, chromosomes))
    cohort = Path(workdir) / 'cohort.raw.vcf.gz'
    runtime.run([cfg['tools'].get('bcftools', 'bcftools'), 'concat', '-Oz', '-o', str(cohort), *map(str, files)], 'concat-reference-order')
    checked_vcf(runtime, cfg, cohort, ids, 'cohort')
    (Path(workdir) / 'scatter-plan.json').write_text(json.dumps(parameters, indent=2) + '\n')
    (Path(workdir) / 'scope.json').write_text(json.dumps(scope, indent=2) + '\n')
    return {'outputs': ['cohort.raw.vcf.gz', 'cohort.raw.vcf.gz.tbi', 'scatter-plan.json', 'scope.json'],
            'metadata': parameters}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--cpus-per-job', type=int, default=6)
    parser.add_argument('--memory-gib', type=int, default=24)
    parser.add_argument('--heap-gib', type=int, default=16)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    os.umask(0o077)
    cfg = json.loads(args.config.read_text())
    runtime = Runtime(cfg, args.execute, stage='joint')
    if not args.execute:
        plan = source_plan(cfg, runtime)
        print(json.dumps({'status': 'PLAN', 'chromosomes': plan[3], 'workers': args.workers,
                          'total_CPUs': args.workers * args.cpus_per_job, 'total_memory_GiB': args.workers * args.memory_gib}))
        return
    busy = subprocess.check_output(['nvidia-smi', '--query-compute-apps=pid', '--format=csv,noheader'], text=True).strip()
    require(not busy, 'GPU jobs are active; CPU scatter must wait for stage3 completion.')
    lock_path = Path(cfg['data_root']) / 'private/.downstream-workflow.lock'
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        reference, ids, rows, chromosomes, scope, inputs, sources = source_plan(cfg, runtime)
        params = {'implementation': 'chromosome CPU scatter', 'code_sha256': sha256(Path(__file__).resolve()),
                  'workers': args.workers, 'cpus': args.cpus_per_job, 'memory': args.memory_gib, 'heap': args.heap_gib}
        producer = lambda directory: scatter(cfg, directory, runtime.child(directory, 'joint'), args.workers,
                                              args.cpus_per_job, args.memory_gib, args.heap_gib)
        runtime.run_stage('joint', inputs, params, producer)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, subprocess.SubprocessError) as error:
        print('SCATTER FAILED: ' + str(error), flush=True)
        raise SystemExit(1)
