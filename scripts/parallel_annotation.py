#!/usr/bin/env python3
"""Bounded chromosome VEP execution around an unchanged, sealed annotation runtime.

Only `run --execute` hands over an active queue. `prove --execute` uses verified
public controls with synthetic genotypes; prepare and default run are read-only.
"""
import argparse
import concurrent.futures
import contextlib
import copy
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import sys
import time
import uuid


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp-' + uuid.uuid4().hex)
    with temporary.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
    os.replace(temporary, path)


def opened(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()


def guard_deadline(deadline):
    """Include pure Python table/report work in the original absolute VM bound."""
    seconds = (deadline - dt.datetime.now(dt.timezone.utc)).total_seconds()
    need(seconds > 0, 'Original VM allowance exhausted.')
    def expired(signum, frame):
        raise KeyboardInterrupt('Original VM allowance deadline reached')
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)


def header_key(line):
    """Normalize only recorded time and input/output command provenance."""
    if line.startswith('##bcftools_viewCommand='):
        return '##bcftools_viewCommand=EXECUTION_PATH_AND_TIME\n'
    if line.startswith('##VEP='):
        return re.sub(r' time="[^"]+"', ' time="EXECUTION_TIME"', line)
    if line.startswith('##VEP-command-line='):
        return re.sub(r'(--(?:input_file|output_file)) \S+', r'\1 EXECUTION_PATH', line)
    return line


def partition(source, directory):
    """Preserve input chromosome blocks, order, headers, alleles and sample fields."""
    directory.mkdir()
    headers, shards, current, handle, seen = [], [], None, None, set()
    try:
        with opened(source) as stream:
            for line in stream:
                if line.startswith('#'):
                    need(not shards, 'Unexpected header after a data record.')
                    headers.append(line)
                    continue
                fields = line.rstrip('\r\n').split('\t')
                chrom = fields[0]
                need(len(fields) >= 10 and re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2])', chrom),
                     'Parallel primary annotation supports autosomes 1-22 only.')
                if current != chrom:
                    need(chrom not in seen, 'Repeated chromosome block; cannot preserve sorted input.')
                    if handle:
                        handle.close()
                    seen.add(chrom)
                    current = chrom
                    destination = directory / (chrom + '.vcf.gz')
                    handle = gzip.open(destination, 'xt')
                    handle.writelines(headers)
                    shards.append(dict(chromosome=chrom, input=str(destination), records=0))
                handle.write(line)
                shards[-1]['records'] += 1
    finally:
        if handle:
            handle.close()
    need(shards and any(h.startswith('#CHROM') for h in headers), 'Empty/malformed sequence partition.')
    return shards


def merge_shards(source, shards, destination):
    """Validate complete records and scientific headers before ordered concatenation."""
    reference_headers, total = None, 0
    with contextlib.ExitStack() as stack:
        original = stack.enter_context(opened(source))
        source_records = (line for line in original if not line.startswith('#'))
        out = stack.enter_context(Path(destination).open('x'))
        for shard in shards:
            headers, count = [], 0
            with opened(shard['output']) as stream:
                for line in stream:
                    if line.startswith('#'):
                        headers.append(line)
                        continue
                    if count == 0:
                        canonical = [header_key(h) for h in headers]
                        if reference_headers is None:
                            reference_headers = canonical
                            out.writelines(headers)
                        else:
                            need(canonical == reference_headers, 'Shard scientific headers differ.')
                    raw = next(source_records, None)
                    need(raw is not None, 'VEP added records.')
                    before, after = raw.rstrip('\r\n').split('\t'), line.rstrip('\r\n').split('\t')
                    need(before[:7] == after[:7] and before[8:] == after[8:],
                         'VEP changed input order, CHROM/POS/ID/REF/ALT/QUAL/FILTER, FORMAT or sample values.')
                    need(after[0] == shard['chromosome'], 'Shard chromosome mismatch.')
                    old_info = set(before[7].split(';')) - {'.'}
                    need(old_info.issubset(set(after[7].split(';'))), 'VEP discarded original INFO fields.')
                    out.write(line)
                    count += 1
                    total += 1
            need(count == shard['records'], 'VEP omitted an input record.')
        need(next(source_records, None) is None, 'VEP omitted trailing input records.')
    return total


def frozen(binding_path):
    binding = load(binding_path)
    package = Path(binding['code_package']['path'])
    need(sha(package) == binding['code_package']['sha256'], 'Frozen package manifest changed.')
    for item in load(package)['files']:
        path = (package.parent / item['path']).resolve(strict=True)
        need(path.is_relative_to(package.parent.resolve()) and not path.is_symlink()
             and path.stat().st_size == item['bytes'] and sha(path) == item['sha256'], 'Frozen helper changed.')
    sys.path.insert(0, str(package.parent / 'scripts'))
    import annotation_queue as queue
    import annotation_stages as stages
    import workflow_common as common
    import workflow
    queue.verify_bindings(binding)
    return binding, queue, stages, common, workflow


def runtime_class(common):
    class ParallelRuntime(common.Runtime):
        workers = 5
        memory_gib = 16

        def code_identity(self):
            return dict(super().code_identity(), parallel_annotation=sha(Path(__file__)))

        def run(self, argv, logname='command', **kwargs):
            if logname != 'offline-vep':
                return super().run(argv, logname, **kwargs)
            argv = list(map(str, argv))
            need('--network=none' in argv and '--pull=never' in argv and argv[argv.index('--fork') + 1] == '4',
                 'Expected exact frozen offline VEP command with four forks.')
            shards = partition(self.workdir / 'sequence_input.vcf.gz', self.workdir / 'chromosome-scatter')
            # Initialize Docker access before threads start.
            self.docker_prefix()
            def job(shard):
                command = argv[:]
                command[command.index('--name') + 1] = 'marina-wf-' + uuid.uuid4().hex[:12]
                command[command.index('--memory') + 1] = str(self.memory_gib) + 'g'
                command[command.index('--cpus') + 1] = '4'
                relative = 'chromosome-scatter/' + shard['chromosome']
                command[command.index('--input_file') + 1] = '/input/' + relative + '.vcf.gz'
                command[command.index('--output_file') + 1] = '/output/' + relative + '.annotated.vcf.gz'
                shard['output'] = str(self.workdir / (relative + '.annotated.vcf.gz'))
                local = common.Runtime(self.cfg, True, self.workdir, 'annotate', self.started)
                local._docker = self._docker
                local._containers = self._containers
                local.run(command, 'chromosome-scatter/' + shard['chromosome'] + '-vep', **kwargs)
                return shard
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                futures = [executor.submit(job, shard) for shard in shards]
                try:
                    for future in concurrent.futures.as_completed(futures):
                        future.result()
                except BaseException:
                    for future in futures:
                        future.cancel()
                    for name in list(self._containers):
                        import subprocess
                        subprocess.run(self.docker_prefix() + ['stop', '--time', '10', name],
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=20)
                    raise
            merged = self.workdir / 'chromosome-scatter/ordered.annotated.vcf'
            total = merge_shards(self.workdir / 'sequence_input.vcf.gz', shards, merged)
            super().run([self.cfg['tools'].get('bcftools', 'bcftools'), 'view', '--no-version', '-Oz', '-o',
                         str(self.workdir / 'sequence.annotated.vcf.gz'), str(merged)], 'compress-scatter')
            save(self.workdir / 'parallel-annotation-producer.json', dict(schema_version=1,
                 producer='parallel_annotation.py around unchanged frozen annotation_stages.annotate',
                 producer_sha256=sha(Path(__file__)), frozen_helpers=super().code_identity(),
                 workers=self.workers, cpus_per_worker=4, memory_gib_per_worker=self.memory_gib,
                 network='none', records=total, shards=shards,
                 scientific_fields_unchanged=True, nonsequence_preservation='frozen merge_vep_preserving_all'))
            return ''
    return ParallelRuntime


def compare_vcfs(left, right):
    count = 0
    import itertools
    with opened(left) as first, opened(right) as second:
        for a, b in itertools.zip_longest(first, second):
            need(a is not None and b is not None and header_key(a) == header_key(b),
                 'Serial/scatter VCF scientific content differs.')
            count += int(not a.startswith('#'))
    return count


def prove(args):
    binding, queue, stages, common, _ = frozen(args.bindings)
    public = load(args.public_controls / 'public-control-receipt.json')
    need(public['status'] == 'PASS' and public['patient_data_used'] is False and public['public_variants'] == 20,
         'A verified public twenty-control receipt is required.')
    source = args.public_controls / 'public-controls.normalized.vcf'
    need(sha(source) == public['artifacts'][source.name]['sha256'], 'Public synthetic input changed.')
    with source.open() as stream:
        sample_line = next(line for line in stream if line.startswith('#CHROM'))
    need(all(s.startswith('PUBLIC_') for s in sample_line.rstrip().split('\t')[9:]), 'Only public synthetic samples allowed.')
    if not args.execute:
        return dict(status='PUBLIC_PROOF_PLAN_ONLY', patient_data_used=False)
    guard_deadline(queue.budget_deadline(binding['budget']['vm_created_utc'],
                   binding['budget']['hourly_usd'], binding['budget']['planned_vm_usd']))
    args.output_root.mkdir(mode=0o700, parents=True, exist_ok=False)
    cfg = load(binding['annotation_config']['path'])
    cfg['annotation'].update(input_vcf=str(source), forks=4, memory_gib=16, make_plots=False, export_docx=False)
    results = []
    for mode in ('serial', 'scatter'):
        directory = args.output_root / mode
        directory.mkdir()
        cfg['output_root'] = str(directory)
        cls = common.Runtime if mode == 'serial' else runtime_class(common)
        runtime = cls(cfg, True, directory, 'annotate')
        if mode == 'scatter':
            runtime.workers = 2
        started = time.monotonic()
        stages.annotate(cfg, directory, runtime)
        results.append(dict(mode=mode, seconds=time.monotonic() - started))
    serial, scattered = args.output_root / 'serial', args.output_root / 'scatter'
    count = compare_vcfs(serial / 'annotated.vcf.gz', scattered / 'annotated.vcf.gz')
    need(count == 20, 'Expected twenty public alleles.')
    scatter = load(scattered / 'parallel-annotation-producer.json')
    need(len(scatter['shards']) >= 2, 'Public scatter proof requires at least two contigs.')
    for name in ('transcript_samples.tsv.gz', 'intersections.tsv.gz'):
        with gzip.open(serial / name, 'rb') as a, gzip.open(scattered / name, 'rb') as b:
            need(a.read() == b.read(), 'Frozen downstream tables differ: ' + name)
    for name in ('annotation_counts.json', 'consequence_counts.json', 'resource_availability.json'):
        need(load(serial / name) == load(scattered / name), 'Frozen downstream metadata differs: ' + name)
    # Verify bgzip/TBI usability, including random access for every scattered chromosome.
    runtime = common.Runtime(cfg, True, args.output_root, 'annotate')
    bcftools = cfg['tools'].get('bcftools', 'bcftools')
    for directory in (serial, scattered):
        for shard in scatter['shards']:
            runtime.run([bcftools, 'view', '-H', '-r', shard['chromosome'], str(directory / 'annotated.vcf.gz')],
                        directory.name + '-index-' + shard['chromosome'], capture=True)
    receipt = dict(schema_version=1, status='PASS', patient_data_used=False, synthetic_genotypes=True,
        public_alleles=count, contigs=len(scatter['shards']), producer_sha256=sha(Path(__file__)),
        binding_sha256=sha(args.bindings), source_sha256=sha(source),
        resource_manifest_sha256=sha(cfg['annotation']['resource_manifest']),
        frozen_package=binding['code_package'], timings=results,
        exact_VCF_records_all_INFO_CSQ_GT_FORMAT_FILTER_and_order=True,
        exact_transcript_sample_and_intersection_tables=True, bgzip_TBI_random_access=True,
        header_normalization=['VEP execution time only', 'VEP input/output filenames only', 'bcftools view command paths/time only'],
        artifacts=common.artifact_manifest(args.output_root))
    save(args.output_root / 'public-scatter-proof.json', receipt)
    return receipt


def process_identity(pid):
    path = Path('/proc') / str(pid)
    try:
        stat = (path / 'stat').read_text().rsplit(')', 1)[1].split()
        if stat[0] == 'Z':
            return None
        return dict(pid=pid, start_ticks=stat[19], argv=[a.decode() for a in (path / 'cmdline').read_bytes().split(b'\0') if a])
    except FileNotFoundError:
        return None


def verify_proof(path, binding_path):
    proof = load(path)
    need(proof['status'] == 'PASS' and proof['patient_data_used'] is False
         and proof['producer_sha256'] == sha(Path(__file__)) and proof['binding_sha256'] == sha(binding_path),
         'Public proof must match this exact producer and binding.')
    for name, item in proof['artifacts'].items():
        p = (path.parent / name).resolve(strict=True)
        need(p.is_relative_to(path.parent.resolve()) and p.stat().st_size == item['bytes'] and sha(p) == item['sha256'],
             'Public proof artifact changed.')
    return proof


def prepare(args):
    binding, queue, _, _, _ = frozen(args.bindings)
    verify_proof(args.proof, args.bindings)
    need(not args.output.exists(), 'Choose a new immutable candidate rule.')
    identities = [process_identity(pid) for pid in (args.controller_pid, args.queue_pid)]
    need(all(identities), 'Expected live controller and queue for explicit user handoff.')
    need('continue_annotation_reference_v2.py' in {Path(a).name for a in identities[0]['argv']}
         and 'annotation_queue.py' in {Path(a).name for a in identities[1]['argv']}
         and str(args.bindings) in identities[1]['argv'], 'Wrong controller/queue process identity.')
    deadline = queue.budget_deadline(binding['budget']['vm_created_utc'], binding['budget']['hourly_usd'], binding['budget']['planned_vm_usd'])
    rule = dict(schema_version=1, status='CANDIDATE_USER_LAUNCH_ONLY',
        binding_path=str(args.bindings.resolve()), binding_sha256=sha(args.bindings),
        producer_path=str(Path(__file__).resolve()), producer_sha256=sha(Path(__file__)),
        proof_path=str(args.proof.resolve()), proof_sha256=sha(args.proof),
        old_processes=identities, output_root=binding['output_root'],
        workers=5, cpus_per_worker=4, memory_gib_per_worker=16, deadline_utc=deadline.isoformat(),
        output_policy='same primary output; unchanged completed stages; explicit new producer receipt',
        no_old15_modification=True, no_patient_API_queries=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    save(args.output, rule)
    return dict(status='PREPARED_NOT_STARTED', rule=str(args.output), rule_sha256=sha(args.output))


def run(args):
    rule = load(args.rule)
    need(args.reviewed_rule_sha256 and sha(args.rule) == args.reviewed_rule_sha256, 'Exact reviewed rule SHA256 required.')
    need(rule['producer_sha256'] == sha(Path(__file__)) and rule['producer_path'] == str(Path(__file__).resolve()), 'Producer changed.')
    need(sha(rule['binding_path']) == rule['binding_sha256'] and sha(rule['proof_path']) == rule['proof_sha256'], 'Reviewed input changed.')
    binding, queue, stages, common, workflow = frozen(rule['binding_path'])
    verify_proof(Path(rule['proof_path']), rule['binding_path'])
    cfg = load(binding['annotation_config']['path'])
    deadline = queue.budget_deadline(binding['budget']['vm_created_utc'], binding['budget']['hourly_usd'], binding['budget']['planned_vm_usd'])
    need(deadline.isoformat() == rule['deadline_utc'], 'Unchanged budget deadline required.')
    root = Path(binding['output_root'])
    # Do this before any handoff: never disturb an already complete annotation/report.
    if queue.verify_stage(root / 'stages/07_report'):
        return dict(status='ALREADY_COMPLETE_VERIFIED_NO_HANDOFF', report=str(root / 'stages/07_report/report.html'))
    if queue.verify_stage(root / 'stages/06_annotate'):
        return dict(status='ANNOTATION_ALREADY_COMPLETE_NO_HANDOFF_KEEP_EXISTING_REPORT_QUEUE')
    if not args.execute:
        return dict(status='PLAN_ONLY_SERIAL_CONTINUES', rule=rule,
                    note='Execute only after root code/proof review. Partial serial VEP is preserved but not adopted as completed.')
    guard_deadline(deadline)
    import fcntl
    with (args.rule.parent / '.parallel-user-launch.lock').open('a') as launch_lock:
        fcntl.flock(launch_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        need(dt.datetime.now(dt.timezone.utc) < deadline, 'Original VM allowance exhausted.')
        # Exact PID + start-time + argv prevents PID-reuse or unrelated-process kills.
        # The controller has no restart policy; stop it first, then its owned queue.
        for expected in rule['old_processes']:
            actual = process_identity(expected['pid'])
            need(actual is None or actual == expected, 'Old process identity changed; prepare a fresh handoff.')
        # Recheck completion immediately before taking control.
        if queue.verify_stage(root / 'stages/06_annotate'):
            return dict(status='ANNOTATION_COMPLETED_DURING_REVIEW_NO_HANDOFF')
        for expected in rule['old_processes']:
            if process_identity(expected['pid']) == expected:
                os.kill(expected['pid'], signal.SIGTERM)
        wait_end = time.monotonic() + 75
        while any(process_identity(p['pid']) is not None for p in rule['old_processes']):
            need(time.monotonic() < wait_end, 'Old queue did not exit cleanly; refusing concurrent producers.')
            time.sleep(1)
        with (Path(binding['data_root']) / 'private/.annotation-continuation.lock').open('a') as annotation_lock:
            fcntl.flock(annotation_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            need(not any('workflow.py' in {Path(a).name for a in (process_identity(int(p.name)) or {}).get('argv', [])}
                         and any(a in (process_identity(int(p.name)) or {}).get('argv', []) for a in ('annotate', 'report'))
                         for p in Path('/proc').glob('[0-9]*')), 'Another annotation/report producer remains active.')
            # SIGTERM/INT become exceptions so the unchanged runtime records a failed
            # partial and stops only its exact named containers.
            def cancelled(signum, frame):
                raise KeyboardInterrupt('User cancellation')
            signal.signal(signal.SIGTERM, cancelled)
            signal.signal(signal.SIGINT, cancelled)
            cls = runtime_class(common)
            runtime = cls(cfg, True, stage='annotate')
            inputs = workflow.dependencies('annotate', cfg, runtime)
            parameters = dict(stage='annotate', producer='parallel_annotation.py',
                producer_sha256=sha(Path(__file__)), rule_sha256=sha(args.rule),
                public_proof_sha256=sha(rule['proof_path']), workers=5, cpus_per_worker=4, memory_gib_per_worker=16,
                frozen_package=binding['code_package'])
            def produce(directory):
                child = cls(cfg, True, directory, 'annotate', runtime.started)
                return stages.annotate(cfg, directory, child)
            # A final serial promotion may race with the user handoff. Preserve it
            # and continue reporting rather than checking it against a scatter identity.
            if not queue.verify_stage(root / 'stages/06_annotate'):
                runtime.run_stage('annotate', inputs + [args.rule, Path(rule['proof_path']), Path(__file__)], parameters, produce)
            need(queue.verify_stage(root / 'stages/06_annotate'), 'Parallel annotation did not verify.')
            if queue.verify_stage(root / 'stages/07_report'):
                return dict(status='REPORT_COMPLETED_DURING_HANDOFF_VERIFIED_NO_REPLACEMENT',
                            report=str(root / 'stages/07_report/report.html'))
            report = cls(cfg, True, stage='report', started=runtime.started)
            report_inputs = workflow.dependencies('report', cfg, report)
            report.run_stage('report', report_inputs, dict(stage='report', continuation=parameters),
                lambda directory: stages.report(cfg, directory, report.child(directory, 'report')))
            need(queue.verify_stage(root / 'stages/07_report'), 'Report did not verify.')
            save(args.rule.parent / 'parallel-run-status.json', dict(status='RESEARCH_REPORT_COMPLETE',
                rule_sha256=sha(args.rule), output_root=str(root), clinical_classification_complete=False))
            return dict(status='RESEARCH_REPORT_COMPLETE', report=str(root / 'stages/07_report/report.html'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    proof = sub.add_parser('prove')
    proof.add_argument('--bindings', type=Path, required=True)
    proof.add_argument('--public-controls', type=Path, required=True)
    proof.add_argument('--output-root', type=Path, required=True)
    proof.add_argument('--execute', action='store_true')
    setup = sub.add_parser('prepare')
    setup.add_argument('--bindings', type=Path, required=True)
    setup.add_argument('--proof', type=Path, required=True)
    setup.add_argument('--output', type=Path, required=True)
    setup.add_argument('--controller-pid', type=int, required=True)
    setup.add_argument('--queue-pid', type=int, required=True)
    launch = sub.add_parser('run')
    launch.add_argument('--rule', type=Path, required=True)
    launch.add_argument('--reviewed-rule-sha256', required=True)
    launch.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    os.umask(0o077)
    print(json.dumps({'prove': prove, 'prepare': prepare, 'run': run}[args.action](args), indent=2))


if __name__ == '__main__':
    main()
