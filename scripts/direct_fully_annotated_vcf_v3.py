#!/usr/bin/env python3
"""Reviewed direct VCF producer; no legacy 06/07 receipt or process actuation."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import multiprocessing
import os
import re
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import uuid

CORE_SHA = '0c4532fdfe316441cf1eb762297a63a725450a7a84ba082db1fda92407b42fce'
END = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
C = CONTEXT = CODEC = RULE = CORE_SINK = CORE_INIT = None
PARENT_DRIVER_SHA = '18b25475689a2ebe2d1a18b191e99f8e44562a0d94fd8a4797abc7415e7cbeb5'


def need(ok, message):
    if not ok: raise ValueError(message)


def load(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''): digest.update(block)
    return digest.hexdigest()


def pin(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def verify(item):
    path = Path(item['path'])
    need(not path.is_symlink() and path.is_file() and path.stat().st_size == item['bytes'] and sha(path) == item['sha256'], 'A reviewed direct producer input changed.')
    return path


def module(item, name):
    path = verify(item); sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec); sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2); stream.write('\n')


def process_identity(pid):
    try:
        path = Path('/proc') / str(pid); parts = (path / 'stat').read_text().rsplit(')', 1)[1].split()
        if parts[0] in ('Z', 'X'): return None
        return dict(pid=pid, start_ticks=parts[19], argv=[a.decode() for a in (path / 'cmdline').read_bytes().split(b'\0') if a])
    except FileNotFoundError: return None


def limit_deadline(end):
    remaining = (min(END, end) - dt.datetime.now(dt.timezone.utc)).total_seconds()
    need(remaining > 0, 'Original VM allowance exhausted.')
    def expired(signum, frame): raise TimeoutError('Original absolute VM allowance reached.')
    signal.signal(signal.SIGALRM, expired); signal.setitimer(signal.ITIMER_REAL, remaining)


def prepare(args):
    core = module(pin(args.core_module), 'candidate_core_prepare')
    need(sha(args.core_module) == CORE_SHA, 'Exact public-proved core module required.')
    snapshot = load(args.snapshot_manifest); binding, queue, _, _ = core.frozen(verify(snapshot['binding']))
    proof = load(args.core_public_proof)
    need(proof['status'] == 'PASS_PUBLIC_DIRECT_PACKET_PROOF' and proof['patient_data_used'] is False and proof['producer']['sha256'] == CORE_SHA
         and proof['snapshot_manifest']['sha256'] == sha(args.snapshot_manifest), 'Exact core public proof/snapshot required.')
    names = ['direct_public_context_adapter.py', 'spliceai_vcf_codec.py', 'spliceai_packet_context.py', 'embed_vcf_evidence.py', 'export_tissue_protein_context.py', 'join_combined_annotation_master.py']
    package = args.context_module.parent; code = [pin(package / name) for name in names]
    accepted = Path(binding['accepted_filter']['path']); source = pin(accepted / 'cohort.filtered.vcf.gz')
    deadline = queue.budget_deadline(binding['budget']['vm_created_utc'], binding['budget']['hourly_usd'], binding['budget']['planned_vm_usd'])
    need(deadline == END, 'Unchanged original allowance required.')
    old = process_identity(args.superseded_pid)
    need(old and any(Path(v).name in ('parallel_annotation.py', 'direct_fully_annotated_vcf.py', 'direct_fully_annotated_vcf_v2.py') for v in old['argv']) and '--execute' in old['argv'], 'Record exact running producer; this script stops nothing.')
    need(not args.output.exists(), 'A new canonical final VCF root is required.')
    rule = dict(schema_version=1, status='CANDIDATE_DIRECT_CORE33_V2_NO_HANDOFF', root_reviewed=False, driver=pin(Path(__file__)), core=pin(args.core_module),
        snapshot_manifest=pin(args.snapshot_manifest), core_public_proof=pin(args.core_public_proof), context_module=pin(args.context_module),
        codec_module=pin(package / 'embed_vcf_evidence.py'), context_code=code, context_index=pin(args.context_index), context_data=load(args.context_index)['index'],
        accepted_filter_vcf=source,
        accepted_filter_receipt=pin(accepted / 'workflow-receipt.json'), superseded_process=old, output=str(args.output.resolve()),
        state_directory=str(args.state_directory.resolve()), workers=24, worker_cpu_ids=sorted(os.sched_getaffinity(0))[2:26], memory_allocation_ceiling_bytes=48 * 1024**3, expected_rows_per_batch=50000, deadline_utc=deadline.isoformat(),
        min_free_bytes=150 * 1024**3, max_output_bytes=48 * 1024**3, max_part_plaintext_bytes=64 * 1024**2,
        no_legacy_receipt_created=True, no_primary_process_actuation=True, no_network=True)
    args.rule.parent.mkdir(mode=0o700, parents=True, exist_ok=True); save(args.rule, rule)
    return dict(status=rule['status'], rule=str(args.rule), rule_sha256=sha(args.rule), no_process_stopped=True)


def initialize_final_worker():
    CORE_INIT()
    C.G['context'] = CONTEXT.Context(RULE['context_index']['path'], verify_index=False)


def prepare_resume(args):
    need(sha(args.parent_rule) == args.parent_reviewed_rule_sha256, 'Exact ROOT-reviewed parent rule required.')
    parent = load(args.parent_rule)
    need(parent.get('root_reviewed') is True and parent.get('status') == 'ROOT_REVIEWED_DIRECT_CORE33_V2_RULE'
         and parent['driver']['sha256'] == PARENT_DRIVER_SHA, 'Only the reviewed original24 producer may supply cached batches.')
    verify(parent['driver'])
    old_work = args.old_work.resolve(strict=True)
    need(old_work.parent == Path(parent['state_directory']).resolve() and old_work.name.startswith('work-'), 'Exact parent-owned work directory required.')
    state_path = old_work.parent / 'status.json'; state = load(state_path)
    need(type(state.get('pid')) is int and process_identity(state['pid']) is None, 'Old producer must already be stopped by ROOT; this script stops nothing.')
    need(state['expected_rows'] == 48529940 and state['batches'] == 969, 'Exact observed original cohort plan required.')
    need(not Path(parent['output']).exists(), 'Refusing existing canonical output.')
    rule = dict(parent, status='CANDIDATE_DIRECT_CORE33_V3_RESUME_NO_HANDOFF', root_reviewed=False, driver=pin(Path(__file__)),
        parent_reviewed_rule=pin(args.parent_rule), parent_driver=parent['driver'], resume_work=str(old_work),
        resume_parent_state=pin(state_path), resume_public_proof=pin(args.resume_public_proof), resume_cache_audit=pin(args.resume_cache_audit),
        superseded_process=dict(pid=state['pid'], historical_identity_basis='ROOT-reviewed parent rule and completed owned batch artifacts; process already absent'),
        state_directory=str(args.state_directory.resolve()), expected_rows_per_batch=5000, expected_total_records=907845,
        expected_total_rows=48529940, reused_body_tokens_must_remain_byte_identical=True)
    need(args.state_directory.resolve() != old_work.parent, 'Resume requires separate new state directory.')
    save(args.rule, rule)
    return dict(status=rule['status'], rule=str(args.rule), rule_sha256=sha(args.rule))


def complete_cache(directory, job, max_part_bytes):
    """Accept only fully committed owned batch artifacts, never a partial body."""
    names = ('body.vcf.gz', 'batch-proof.json', 'full-vector-proof.json', 'annotation_counts.json', 'consequence_counts.json')
    if not all((directory / name).is_file() for name in names): return None
    need(not directory.is_symlink() and all(not (directory / name).is_symlink() for name in names), 'Cache symlink rejected.')
    proof, full, counts, consequences = [load(directory / name) for name in names[1:]]
    need(all(proof[k] == job[k] for k in ('id', 'first_i', 'last_i')), 'Cached source ordinals changed.')
    need(proof['row_count'] == job['expected_rows'] and proof['record_count'] == job['records'], 'Cached source row/record counts changed.')
    need(counts['records'] == job['records'] and counts['transcript_sample_rows'] == job['expected_rows'], 'Cached counts incomplete.')
    need(all(re.fullmatch('[0-9a-f]{64}', value) for value in (proof['canonical_row_sha256'], full['canonical_row_sha256'])), 'Invalid cached vector commitment.')
    need(type(full['plaintext_bytes']) is int and 0 < full['plaintext_bytes'] <= max_part_bytes, 'Cached plaintext cap changed.')
    need(isinstance(consequences.get('definitions'), dict) and isinstance(consequences.get('samples'), dict)
         and isinstance(consequences.get('cohort'), dict) and isinstance(consequences.get('schema_version'), int), 'Cached consequence counters incomplete.')
    need(Path(proof['body']['path']).resolve() == (directory / 'body.vcf.gz').resolve(), 'Cache body provenance path changed.')
    body = verify(proof['body'])
    with body.open('rb') as stream:
        stream.seek(-len(C.BGZF_EOF), 2); need(stream.read() == C.BGZF_EOF, 'Cached BGZF EOF missing.')
    return dict(proof=proof, counts=counts, consequences=consequences, full_proof=full,
        cache_artifact_pins={name: pin(directory / name) for name in names[1:]})


def resume_partition(source, old_work, work, samples, rows_per_batch, max_part_bytes):
    """Verify exact original records; reuse committed bodies and split only missing jobs."""
    old_work = Path(old_work).resolve(strict=True)
    need(not old_work.is_symlink(), 'Resume work symlink rejected.')
    inputs = sorted((old_work / 'inputs').glob('batch-*.vcf.gz'))
    need(inputs, 'No owned original balanced inputs.')
    headers, entries, origins, ordinal, total_records, reused_plaintext = [], [], [], 1, 0, 0
    with C.text(source) as source_stream:
        source_records = iter(source_stream)
        first_record = None
        for line in source_records:
            if line.startswith('#'): headers.append(line)
            else: first_record = line; break
        import itertools
        source_records = itertools.chain([first_record] if first_record is not None else [], source_records)
        for old_id, path in enumerate(inputs):
            need(path.name == 'batch-' + str(old_id).zfill(5) + '.vcf.gz' and not path.is_symlink(), 'Original input sequence changed.')
            records, rows, observed_headers = 0, 0, []
            with C.text(path) as stream:
                for line in stream:
                    if line.startswith('#'):
                        need(records == 0, 'Old input header after records.'); observed_headers.append(line); continue
                    need(line == next(source_records, None), 'Old input record differs from pinned completed VEP source.')
                    columns = line.rstrip('\r\n').split('\t')
                    need(len(columns) == 9 + len(samples), 'Original four-sample input contract changed.')
                    info = dict((v.split('=', 1) + [''])[:2] for v in columns[7].split(';'))
                    rows += len(samples) * (len(info['CSQ'].split(',')) if info.get('CSQ') else 1); records += 1
            need(observed_headers == headers and records > 0, 'Old input headers/records incomplete.')
            old_job = dict(id=old_id, input=str(path), first_i=ordinal, last_i=ordinal + rows - 1, records=records, expected_rows=rows)
            cached = complete_cache(old_work / ('result-' + str(old_id).zfill(5)), old_job, max_part_bytes)
            origin = dict(old_job=old_job, input_pin=pin(path), reuse=cached is not None, new_batch_ids=[])
            if cached is not None:
                new_id = len(entries); cached['proof'] = dict(cached['proof'], id=new_id)
                entries.append(dict(cached=cached)); origin['new_batch_ids'].append(new_id)
                origin['cache_artifact_pins'] = cached['cache_artifact_pins']; origin['cached_body'] = cached['proof']['body']
                reused_plaintext += cached['full_proof']['plaintext_bytes']
            else:
                _, split_jobs, split_rows = C.partition(path, work / ('resume-input-' + str(old_id).zfill(5)), samples, rows_per_batch)
                need(split_rows == rows and sum(j['records'] for j in split_jobs) == records, 'Missing-input repartition lost rows/records.')
                for job in split_jobs:
                    job.update(id=len(entries), first_i=job['first_i'] + ordinal - 1, last_i=job['last_i'] + ordinal - 1)
                    entries.append(dict(job=job)); origin['new_batch_ids'].append(job['id'])
            origins.append(origin); ordinal += rows; total_records += records
        need(next(source_records, None) is None, 'Old inputs omit source VEP records.')
    return headers, entries, ordinal - 1, dict(source_origin_mapping=origins, source_record_count=total_records,
        reused_batches=sum(v['reuse'] for v in origins), missing_original_batches=sum(not v['reuse'] for v in origins),
        reused_rows=sum(v['old_job']['expected_rows'] for v in origins if v['reuse']), reused_plaintext_bytes=reused_plaintext,
        regenerated_tasks=sum('job' in e for e in entries), original_source_record_bytes_verified=True, cached_bodies_unchanged=True)


class ExtendedSink:
    def __init__(self, source, output, fields, first_i):
        CORE_SINK.__init__(self, source, output, fields, first_i)
        self.full_digest = C.VectorDigest(); self.extra_fields = C.G['extra_fields']; self.part_plaintext_bytes = 0

    def writerow(self, row):
        ordinal = self.ordinal
        CORE_SINK.writerow(self, row)
        extra = C.G['context'].augment(dict(row), ordinal)
        need(set(extra) == set(self.extra_fields) and all(isinstance(v, str) for v in extra.values()), 'Unknown/missing/non-string context field.')
        self.rows[-1]['values'] += [extra[field] for field in self.extra_fields]
        self.full_digest.row(self.rows[-1]['values'])

    def flush(self):
        if self.current is None: return
        if self.record_count % 1000 == 0:
            need(dt.datetime.now(dt.timezone.utc) < dt.datetime.fromisoformat(RULE['deadline_utc']), 'Original deadline reached.')
            need(shutil.disk_usage(RULE['state_directory']).free >= RULE['min_free_bytes'], 'Direct producer disk reserve reached.')
        columns = self.current.rstrip('\r\n').split('\t')
        packet = dict(v=1, allele=[columns[0], int(columns[1]), columns[3], columns[4]], scope='cohort', sample=None,
            rows=self.rows, membership=self.membership, evidence=C.G['context'].allele_evidence(columns[0], columns[1], columns[3], columns[4]))
        token = CODEC.encode_packet(packet)
        need(CODEC.decode_packet(token) == packet, 'Final packet roundtrip lost an ordered field vector.')
        readable = CODEC.rebuild_readable(packet, C.G['full_fields'])
        need(set(readable) <= set(CODEC.INFO_IDS) and 'WGS_EVIDENCE' not in readable, 'Unknown readable INFO field.')
        previous = {v.split('=', 1)[0] for v in columns[7].split(';')}
        need(not previous.intersection(CODEC.INFO_IDS), 'Overlay namespace collision in accepted VEP component.')
        additions = ['WGS_EVIDENCE=' + token] + [name + '=' + str(value) for name, value in readable.items()]
        columns[7] = ('' if columns[7] == '.' else columns[7] + ';') + ';'.join(additions)
        ending = '\r\n' if self.current.endswith('\r\n') else '\n'
        encoded = ('\t'.join(columns) + ending).encode('utf-8')
        need(CODEC.strip_overlay(encoded) == self.current.encode('utf-8'), 'Original VCF record/INFO/GT changed.')
        self.part_plaintext_bytes += len(encoded)
        need(self.part_plaintext_bytes <= RULE['max_part_plaintext_bytes'], 'One balanced batch exceeds its reviewed plaintext ceiling; no allele was split or dropped.')
        with C.G['output_plaintext_budget'].get_lock():
            need(C.G['output_plaintext_budget'].value + len(encoded) <= RULE['max_output_bytes'], 'Reviewed total packet output ceiling reached.')
            C.G['output_plaintext_budget'].value += len(encoded)
        self.bgzip.stdin.write(encoded); self.record_count += 1
        self.current, self.rows = None, []

    def close(self):
        CORE_SINK.close(self)
        save(Path(self.output_file.name).parent / 'full-vector-proof.json', dict(canonical_row_sha256=self.full_digest.digest.hexdigest(), plaintext_bytes=self.part_plaintext_bytes))


def install(core, context, codec, rule, stages, snapshot, manifest, manifest_path, config, alpha_inputs, work, public_alpha=None):
    global C, CONTEXT, CODEC, RULE, CORE_SINK, CORE_INIT
    C, CONTEXT, CODEC, RULE = core, context, codec, rule
    CORE_SINK, CORE_INIT = core.PacketSink, core.initialize_worker
    core.PacketSink = type('InlineContextPacketSink', (ExtendedSink, CORE_SINK), {})
    fields = list(stages.FIELDS); extra = list(context.APPENDED)
    need(len(set(fields + extra)) == len(fields + extra), 'Context field namespace collision.')
    core.G.update(stages=stages, snapshot=snapshot['snapshot']['path'], alpha_source_rows=snapshot['alpha_source_rows'], work=str(work),
        config=config, manifest=manifest, manifest_path=manifest_path, alpha_inputs=alpha_inputs, genes=stages.gene_evidence(manifest, manifest_path),
        hpo=stages.hpo_map(manifest, manifest_path), bgzip=shutil.which('bgzip'), tabix=shutil.which('tabix'), full_fields=fields + extra, extra_fields=extra,
        output_plaintext_budget=multiprocessing.Value('Q', 0))
    if public_alpha: core.G['public_alpha'] = public_alpha
    need(core.G['bgzip'] and core.G['tabix'], 'Existing bgzip/tabix required.')
    return fields, extra


def schema(fields, proofs, rows, original=None):
    original = fields if original is None else original
    return dict(v=1, codec=C.CODEC, scope='cohort', sample=None, original_fields=original, appended_fields=fields[len(original):],
        full_fields=fields, row_count=rows, ordinal_in_vcf_order=True, vector_digest_contract=C.DIGEST_CONTRACT,
        canonical_header_sha256=hashlib.sha256(C.canonical_row_bytes(fields)).hexdigest(), batches=[C.portable_batch(p) for p in proofs],
        ordered_batch_merkle_root=C.merkle_digest(fields, proofs), producer_kind='DIRECT_FROZEN_ANNOTATION_ROW_V1_NO_LEGACY_06_07_RECEIPTS')


def assemble(headers, results, full_schema, context, output):
    bundles = context.header_bundles()
    new_headers = CODEC.headers(full_schema, CODEC.encode_packet(bundles['gene_context']), CODEC.encode_packet(bundles['public_reference']))
    lines = []
    for line in headers:
        if line.startswith('#CHROM'): lines.extend(new_headers)
        lines.append(line.encode('utf-8'))
    header = output.parent / 'header.temporary.bgzf'
    with header.open('xb') as stream: subprocess.run([C.G['bgzip'], '-c', '-@', '1'], input=b''.join(lines), stdout=stream, check=True)
    with output.open('xb') as out:
        for path in [header] + [verify(r['proof']['body']) for r in results]:
            with path.open('rb') as stream:
                size = path.stat().st_size; stream.seek(size - len(C.BGZF_EOF))
                need(stream.read() == C.BGZF_EOF, 'Missing BGZF EOF block.'); stream.seek(0); remaining = size - len(C.BGZF_EOF)
                while remaining:
                    block = stream.read(min(8 * 1024**2, remaining)); need(block, 'Truncated BGZF body.')
                    need(out.tell() + len(block) <= RULE['max_output_bytes'] and shutil.disk_usage(output.parent).free >= RULE['min_free_bytes'], 'Assembly output ceiling or disk reserve reached.')
                    out.write(block); remaining -= len(block)
        out.write(C.BGZF_EOF)
    header.unlink(); subprocess.run([C.G['tabix'], '-p', 'vcf', str(output)], check=True, stdout=subprocess.DEVNULL)
    added = set(new_headers); restored = []
    with C.text(output) as stream:
        for line in stream:
            raw = line.encode('utf-8')
            if raw not in added: restored.append(raw)
            if line.startswith('#CHROM'): break
    need(b''.join(restored) == ''.join(headers).encode('utf-8'), 'Original header bytes changed.')


def execute(rule, public=None):
    global C
    C = module(rule['core'], 'proved_direct_core'); need(sha(rule['core']['path']) == CORE_SHA, 'Exact public-proved core required.')
    snapshot = load(verify(rule['snapshot_manifest'])); binding_path = verify(snapshot['binding'])
    binding, queue, stages, common = C.frozen(binding_path)
    for item in rule['context_code']: verify(item)
    context_module = module(rule['context_module'], 'direct_context_inline'); codec = module(rule['codec_module'], 'direct_codec_inline')
    verify(rule['context_index']); verify(rule['context_data']); verify(snapshot['snapshot']); verify(rule['core_public_proof'])
    end = dt.datetime.fromisoformat(rule['deadline_utc']); limit_deadline(end)
    if public is None:
        current = process_identity(rule['superseded_process']['pid'])
        need(current is None or current == rule['superseded_process'], 'Superseded PID/start-time/argv identity changed.')
        need(current is None, 'Existing primary still runs; ROOT must perform a separately reviewed handoff. This script stops nothing.')
    output_root = Path(rule['output']); need(not output_root.exists(), 'Refusing replacement of canonical final VCF.')
    state = Path(rule['state_directory']); state.mkdir(mode=0o700, parents=True, exist_ok=True)
    import fcntl
    with (state / 'direct-producer.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        def status(value, **extra):
            temporary = state / 'status.tmp'
            temporary.write_text(json.dumps(dict(status=value, pid=os.getpid(), updated_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                deadline_utc=end.isoformat(), no_legacy_receipt_created=True, **extra), indent=2) + '\n'); os.replace(temporary, state / 'status.json')
        status('VERIFYING_DIRECT_COMPONENTS_NO_LEGACY_06_07')
        context = context_module.Context(rule['context_index']['path'], verify_index=False)
        cfg = load(binding['annotation_config']['path']); mp = Path(cfg['annotation']['resource_manifest'])
        manifest = stages.load_resources(mp, verify_tree=False, local_adapters=cfg['annotation'].get('local_adapters'))
        source = verify(snapshot['annotated_vcf']) if public is None else Path(public['source'])
        alpha_path = Path(binding['alpha']['path']); need(sha(alpha_path) == binding['alpha']['sha256'], 'Alpha manifest changed.')
        alpha_manifest = load(alpha_path); alpha_files = {name: dict(matches=verify(item)) for name, item in snapshot['alpha_matches'].items()}
        alpha_inputs = (alpha_manifest, alpha_files, alpha_path); public_alpha = None
        if public: alpha_inputs = public['alpha_inputs']; public_alpha = (alpha_inputs[1], alpha_inputs[0])
        samples = stages.vcf_samples(source); need(len(samples) == 4, 'Exactly four source samples required.')
        if public: need(all(sample.startswith('PUBLIC_') for sample in samples), 'Public synthetic controls only.')
        work = state / ('work-' + uuid.uuid4().hex); work.mkdir(mode=0o700)
        fields, appended = install(C, context_module, codec, rule, stages, snapshot, manifest, mp, cfg['annotation'], alpha_inputs, work, public_alpha)
        import resource
        resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
        need(shutil.disk_usage(state).free >= rule['min_free_bytes'], 'Direct producer disk reserve reached.')
        status('PLANNING_BALANCED_ALLELE_PRESERVING_BATCHES')
        resume_audit = None
        if rule.get('resume_work'):
            headers, entries, rows, resume_audit = resume_partition(source, rule['resume_work'], work, samples,
                rule['expected_rows_per_batch'], rule['max_part_plaintext_bytes'])
            jobs = [v['job'] for v in entries if 'job' in v]
            need(rows == rule['expected_total_rows'] and resume_audit['source_record_count'] == rule['expected_total_records'], 'Exact whole-source totals changed.')
            if public is None:
                need(resume_audit == load(verify(rule['resume_cache_audit'])), 'Revalidated cache source-origin inventory differs from ROOT-reviewed audit.')
            need(resume_audit['reused_plaintext_bytes'] <= rule['max_output_bytes'], 'Cached plaintext exceeds total ceiling.')
            C.G['output_plaintext_budget'].value = resume_audit['reused_plaintext_bytes']
            save(work / 'resume-audit.json', resume_audit)
        else:
            headers, jobs, rows = C.partition(source, work / 'inputs', samples, rule['expected_rows_per_batch'])
            entries = [dict(job=job) for job in jobs]
        status('DIRECT_FROZEN_ROWS_INLINE_CONTEXT_RUNNING', batches=len(entries), expected_rows=rows, workers=rule['workers'],
            reused_batches=resume_audit['reused_batches'] if resume_audit else 0,
            reused_rows=resume_audit['reused_rows'] if resume_audit else 0, regenerated_tasks=len(jobs))
        need(len(rule['worker_cpu_ids']) == 24 and set(rule['worker_cpu_ids']) <= os.sched_getaffinity(0), 'Reviewed24-CPU allocation unavailable.')
        os.sched_setaffinity(0, set(rule['worker_cpu_ids']))
        started = time.monotonic()
        with multiprocessing.get_context('fork').Pool(rule['workers'], initializer=initialize_final_worker) as pool: produced = pool.map(C.produce_batch, jobs)
        produced = iter(produced)
        results = [entry['cached'] if 'cached' in entry else next(produced) for entry in entries]
        seconds = time.monotonic() - started; counts, consequences = C.merged_counts(results)
        planned_records = resume_audit['source_record_count'] if resume_audit else sum(job['records'] for job in jobs)
        need(counts['records'] == planned_records and counts['transcript_sample_rows'] == rows, 'A planned record/row count changed.')
        expected_i = 1
        for result in results:
            proof = result['proof']
            need(proof['first_i'] == expected_i and proof['last_i'] - proof['first_i'] + 1 == proof['row_count'], 'Resumed batch ordinal gap/overlap.')
            expected_i = proof['last_i'] + 1
        need(expected_i - 1 == rows, 'Resumed total ordinal count differs.')
        prefix_proofs = [r['proof'] for r in results]
        full_proofs = [dict(p, canonical_row_sha256=load(Path(p['body']['path']).parent / 'full-vector-proof.json')['canonical_row_sha256']) for p in prefix_proofs]
        core_schema = schema(fields, prefix_proofs, rows); full_schema = schema(fields + appended, full_proofs, rows, fields)
        staging = output_root.parent / ('.' + output_root.name + '.direct-' + uuid.uuid4().hex); staging.mkdir(mode=0o700)
        output = staging / 'cohort.core-context.vcf.gz'; status('ASSEMBLING_AND_INDEXING_DIRECT_CORE33_VCF', rows=rows)
        assemble(headers, results, full_schema, context, output); context.close()
        need(output.stat().st_size <= rule['max_output_bytes'], 'Direct VCF output ceiling exceeded.')
        source_pins = dict(binding=snapshot['binding'], config=pin(binding['annotation_config']['path']), resources=pin(mp), snapshot_database=snapshot['snapshot'],
            snapshot_manifest=rule['snapshot_manifest'], completed_vep_vcf=pin(source), completed_vep_index=snapshot['index'] if public is None else pin(str(source) + '.tbi'),
            vep_component_producer=snapshot['source_component'], filter_receipt=rule['accepted_filter_receipt'], filter_vcf=rule['accepted_filter_vcf'], alpha_manifest=pin(alpha_path),
            alpha_canonical_matches=pin(alpha_inputs[1]['canonical']['matches']), alpha_isoform_matches=pin(alpha_inputs[1]['isoforms']['matches']), core_public_proof=rule['core_public_proof'], core_module=rule['core'])
        code = dict(final_driver=sha(Path(__file__)), core_module=CORE_SHA, frozen_package_manifest=binding['code_package']['sha256'], context_adapter=rule['context_module']['sha256'], codec=rule['codec_module']['sha256'])
        provenance_pin = None
        if resume_audit is not None:
            code['resume_parent_driver'] = rule['parent_driver']['sha256']
            save(staging / 'resume-provenance.json', dict(parent_reviewed_rule=rule['parent_reviewed_rule'],
                parent_driver=rule['parent_driver'], parent_state=rule['resume_parent_state'], **resume_audit))
            provenance_pin = pin(staging / 'resume-provenance.json'); provenance_pin['path'] = str(output_root / 'resume-provenance.json')
        save(staging / 'core-manifest.json', dict(schema_version=1, status='DIRECT_PACKET_BASELINE_PASS', no_legacy_receipt_created=True, schema=core_schema,
            producer_code_sha256=code, source_pins=source_pins, annotation_counts=counts, consequence_counts=consequences,
            resume_provenance=provenance_pin, synthetic_validation=public is not None))
        save(staging / 'annotation_counts.json', counts); save(staging / 'consequence_counts.json', consequences)
        def final_pin(path):
            item = pin(path); item['path'] = str(output_root / path.name); return item
        final = dict(schema_version=1, status='INTERMEDIATE_CORE33_PASS', final_delivery_complete=False,
            spliceai_storage='PENDING_SEPARATE_ALLELE_INFO_OVERLAY', enriched_vcf=final_pin(output), index=final_pin(Path(str(output) + '.tbi')),
            core_manifest=final_pin(staging / 'core-manifest.json'), core_schema=core_schema, schema=full_schema, samples=samples, producer_code_sha256=code,
            context_index=rule['context_index'], context_data=rule['context_data'],
            resume_provenance=provenance_pin,
            original_records_fields_preserved=True, packet_decode_verified=True, direct_core_prefix_verified=True, no_legacy_receipt_created=True,
            original_header_bytes_sha256=hashlib.sha256(''.join(headers).encode('utf-8')).hexdigest(),
            no_network=True, synthetic_validation=public is not None, row_count=rows, record_count=counts['records'], worker_wall_seconds=seconds,
            max_batch_plaintext_bytes=max(load(Path(p['body']['path']).parent / 'full-vector-proof.json')['plaintext_bytes'] for p in prefix_proofs),
            created_utc=dt.datetime.now(dt.timezone.utc).isoformat())
        save(staging / 'manifest.json', final)
        artifacts = {}
        for path in staging.iterdir():
            if path.is_file(): item = pin(path); artifacts[path.name] = {k: item[k] for k in ('bytes', 'sha256')}
        save(staging / 'workflow-receipt.json', dict(schema_version=1, status='PASS', stage='direct-core-context-v1',
            identity=hashlib.sha256(C.canonical(dict(rule=rule, source_pins=source_pins, code=code))).hexdigest(), artifacts=artifacts,
            metadata=dict(new_direct_producer=True, no_legacy_receipt_created=True)))
        common.promote_exclusive(staging, output_root)
        status(final['status'], manifest_sha256=sha(output_root / 'manifest.json'), rows=rows)
        return final


def run(args):
    need(sha(args.rule) == args.reviewed_rule_sha256, 'Exact reviewed direct producer rule SHA256 required.')
    rule = load(args.rule); verify(rule['driver']); need(rule['driver']['sha256'] == sha(Path(__file__)), 'Driver changed.')
    need(rule.get('root_reviewed') is True and rule.get('status') == 'ROOT_REVIEWED_DIRECT_CORE33_V2_RULE', 'Exact ROOT-reviewed status/boolean required.')
    need(rule['workers'] == 24 and len(rule['worker_cpu_ids']) == 24 and rule['memory_allocation_ceiling_bytes'] <= 48 * 1024**3 and rule['expected_rows_per_batch'] == 5000 and rule['min_free_bytes'] >= 150 * 1024**3
         and rule['max_part_plaintext_bytes'] <= 64 * 1024**2 and dt.datetime.fromisoformat(rule['deadline_utc']) == END, 'Bounded reviewed direct execution required.')
    parent = load(verify(rule['parent_reviewed_rule']))
    need(parent.get('root_reviewed') is True and parent.get('status') == 'ROOT_REVIEWED_DIRECT_CORE33_V2_RULE'
         and parent['driver']['sha256'] == PARENT_DRIVER_SHA and parent['driver'] == rule['parent_driver'], 'Reviewed frozen parent identity changed.')
    for key in ('core','snapshot_manifest','core_public_proof','context_module','codec_module','context_code','context_index','context_data',
                'accepted_filter_vcf','accepted_filter_receipt','output','workers','worker_cpu_ids','memory_allocation_ceiling_bytes',
                'deadline_utc','min_free_bytes','max_output_bytes','max_part_plaintext_bytes'):
        need(rule[key] == parent[key], 'Resume altered reviewed scientific input or resource cap: ' + key)
    proof = load(verify(rule['resume_public_proof']))
    need(proof.get('status') == 'PASS_PUBLIC_MIXED_CACHE_RESUME_PROOF' and proof.get('producer_sha256') == sha(Path(__file__))
         and proof.get('patient_data_used') is False and proof.get('all_ordered_vectors_exact') is True
         and proof.get('counters_exact') is True and proof.get('old_cache_files_unchanged') is True, 'Exact mixed-cache public proof required.')
    verify(rule['resume_parent_state']); audit = load(verify(rule['resume_cache_audit']))
    need(audit.get('original_source_record_bytes_verified') is True and audit.get('cached_bodies_unchanged') is True
         and audit.get('source_record_count') == rule['expected_total_records'], 'Exact resume audit required.')
    if not args.execute: return dict(status='PLAN_ONLY_PRIMARY_UNCHANGED', rule_sha256=sha(args.rule), output=rule['output'])
    return execute(rule)


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('prepare')
    for name in ('core-module', 'snapshot-manifest', 'core-public-proof', 'context-module', 'context-index', 'output', 'state-directory', 'rule'):
        prep.add_argument('--' + name, type=Path, required=True)
    prep.add_argument('--superseded-pid', type=int, required=True)
    resume = sub.add_parser('prepare-resume')
    for name in ('parent-rule','old-work','resume-public-proof','resume-cache-audit','state-directory','rule'):
        resume.add_argument('--' + name, type=Path, required=True)
    resume.add_argument('--parent-reviewed-rule-sha256', required=True)
    action = sub.add_parser('run'); action.add_argument('--rule', type=Path, required=True); action.add_argument('--reviewed-rule-sha256', required=True); action.add_argument('--execute', action='store_true')
    args = parser.parse_args(); os.umask(0o077)
    print(json.dumps({'prepare': prepare, 'prepare-resume': prepare_resume, 'run': run}[args.command](args), indent=2))


if __name__ == '__main__': main()
