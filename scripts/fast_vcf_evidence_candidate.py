#!/usr/bin/env python3
"""Candidate direct VCF packet producer; no legacy 06/07 completion is claimed.

The frozen make_tables loop, annotation_row, quality gates and counters execute
unchanged. Only its CSV sinks and already completed database initialization are
adapted in isolated worker processes. Default commands are plans, not execution.
"""
import argparse
import base64
import contextlib
import csv
import datetime as dt
import gzip
import hashlib
import io
import itertools
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import signal
import sqlite3
import subprocess
import sys
import time
import uuid
import zlib

CODEC = 'json-utf8+zlib+base64url-nopad'
DIGEST_CONTRACT = 'ordered_batch_canonical_TSV_utf8_LF_sha256_v1'
BGZF_EOF = bytes.fromhex('1f8b08040000000000ff0600424302001b0003000000000000000000')
MAX_PACKET = 64 * 1024**2
G = {}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def pin(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def verify(item):
    path = Path(item['path'])
    need(not path.is_symlink() and path.is_file() and path.stat().st_size == item['bytes']
         and sha(path) == item['sha256'], 'A pinned candidate input changed.')
    return path


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def text(path):
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if str(path).endswith('.gz') else Path(path).open(encoding='utf-8', newline='')


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(',', ':')).encode()


def canonical_row_bytes(values):
    """String field vector -> the exact UTF-8/LF frozen TSV row representation."""
    stream = io.StringIO(newline='')
    csv.writer(stream, delimiter='\t', lineterminator='\n').writerow(values)
    return stream.getvalue().encode('utf-8')


def portable_batch(proof):
    return {key: proof[key] for key in ('id', 'first_i', 'last_i', 'row_count', 'record_count', 'canonical_row_sha256')}


def batch_proof_digest(proof):
    """Commit one portable batch descriptor, without deployment paths/timing."""
    return hashlib.sha256(canonical(portable_batch(proof))).hexdigest()


def merkle_digest(fields, ordered_proofs):
    """Ordered one-level Merkle commitment; explicitly not a global TSV SHA."""
    value = dict(contract=DIGEST_CONTRACT,
                 canonical_header_sha256=hashlib.sha256(canonical_row_bytes(fields)).hexdigest(),
                 batches=[portable_batch(proof) for proof in ordered_proofs])
    return hashlib.sha256(canonical(value)).hexdigest()


def encode(value):
    raw = canonical(value)
    need(len(raw) <= MAX_PACKET, 'One whole allele exceeds the packet memory ceiling; no rows were dropped.')
    return base64.urlsafe_b64encode(zlib.compress(raw, 6)).decode().rstrip('=')


def decode(token):
    decoder = zlib.decompressobj()
    raw = decoder.decompress(base64.b64decode(token + '=' * (-len(token) % 4), altchars=b'-_', validate=True), MAX_PACKET + 1)
    need(len(raw) <= MAX_PACKET and decoder.eof and not decoder.unconsumed_tail and not decoder.unused_data,
         'Invalid, oversized or truncated packet.')
    return json.loads(raw)


class VectorDigest:
    def __init__(self):
        self.digest = hashlib.sha256()
        self.writer = csv.writer(self, delimiter='\t', lineterminator='\n')

    def write(self, value):
        self.digest.update(value.encode('utf-8'))
        return len(value)

    def row(self, values):
        self.writer.writerow(values)


def frozen(binding_path):
    binding = load(binding_path)
    package = Path(binding['code_package']['path'])
    need(sha(package) == binding['code_package']['sha256'], 'Frozen code package manifest changed.')
    for item in load(package)['files']:
        path = package.parent / item['path']
        need(path.resolve(strict=True).is_relative_to(package.parent.resolve()) and not path.is_symlink()
             and path.stat().st_size == item['bytes'] and sha(path) == item['sha256'], 'A frozen helper changed.')
    sys.path.insert(0, str(package.parent / 'scripts'))
    import annotation_queue as queue
    import annotation_stages as stages
    import workflow_common as common
    queue.verify_bindings(binding)
    return binding, queue, stages, common


def readonly_database(path):
    connection = sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro&immutable=1', uri=True)
    connection.execute('PRAGMA query_only=ON')
    connection.execute('PRAGMA cache_size=-16384')
    return connection


def deadline_guard(binding, queue):
    budget = binding['budget']
    end = queue.budget_deadline(budget['vm_created_utc'], budget['hourly_usd'], budget['planned_vm_usd'])
    seconds = min(1800, (end - dt.datetime.now(dt.timezone.utc)).total_seconds())
    need(seconds > 0, 'Original VM allowance exhausted.')
    def expired(signum, frame):
        raise TimeoutError('Candidate proof stopped at its bound or unchanged VM allowance.')
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    return end


def compression_limit():
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (128 * 1024**2, 128 * 1024**2))


def snapshot(args):
    binding, queue, _, _ = frozen(args.bindings)
    source = args.source_dir.resolve(strict=True)
    need(source.parent == Path(binding['output_root']) / 'stages' and source.name.startswith('06_annotate.partial-'),
         'Only the exact primary annotation partial may supply a completed component.')
    producer = load(source / 'parallel-annotation-producer.json')
    need(producer.get('scientific_fields_unchanged') is True, 'Completed ordered VEP component proof required.')
    if not args.execute:
        return dict(status='SNAPSHOT_PLAN_ONLY', no_primary_mutation=True)
    deadline_guard(binding, queue)
    need(args.output.resolve().is_relative_to((Path(binding['data_root']) / 'private').resolve()), 'Keep the snapshot private.')
    args.output.mkdir(mode=0o700, parents=True, exist_ok=False)
    source_db = source / 'local-adapter-index.sqlite'
    before = source_db.stat()
    target = args.output / 'annotation-index.snapshot.sqlite'
    # SQLite backup reads a consistent committed database. The live source is
    # opened read-only, never hardlinked, copied with raw cp, or modified.
    original = sqlite3.connect(source_db.resolve().as_uri() + '?mode=ro', uri=True)
    original.execute('PRAGMA query_only=ON')
    with contextlib.closing(sqlite3.connect(target)) as destination:
        original.backup(destination, pages=4096, sleep=0.01)
        need(destination.execute('PRAGMA quick_check').fetchone() == ('ok',), 'Snapshot SQLite integrity check failed.')
        alpha_rows = destination.execute('SELECT count(*) FROM alpha').fetchone()[0]
        evidence_rows = destination.execute('SELECT count(*) FROM evidence').fetchone()[0]
    original.close()
    after = source_db.stat()
    need((before.st_ino, before.st_size, before.st_mtime_ns) == (after.st_ino, after.st_size, after.st_mtime_ns),
         'Live database changed during snapshot; candidate retained for inspection only.')
    result = dict(schema_version=1, status='CONSISTENT_READONLY_SNAPSHOT_CANDIDATE',
        binding=pin(args.bindings), frozen_package=binding['code_package'], snapshot=pin(target),
        source_component=pin(source / 'parallel-annotation-producer.json'),
        annotated_vcf=pin(source / 'annotated.vcf.gz'), index=pin(source / 'annotated.vcf.gz.tbi'),
        alpha_matches={kind: pin(source / ('alphamissense-' + kind + '.matches.tsv')) for kind in ('canonical', 'isoforms')},
        alpha_manifest=binding['alpha'], alpha_source_rows=alpha_rows, evidence_rows=evidence_rows,
        source_database_stat_unchanged=True, source_open_mode='SQLite mode=ro; consistent backup API',
        original_stage06_complete=False, no_legacy_receipt_created=True, patient_database_keep_private=True,
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat())
    save(args.output / 'snapshot-manifest.json', result)
    return result


def partition(source, directory, samples, rows_per_batch):
    """Balanced expected row work; an allele and all its samples stay together."""
    directory.mkdir()
    headers, jobs, output, records, rows, ordinal = [], [], None, 0, 0, 1
    def finish():
        nonlocal output, records, rows, ordinal
        if output:
            output.close()
            jobs[-1].update(records=records, expected_rows=rows, last_i=ordinal + rows - 1)
            ordinal += rows
            records, rows = 0, 0
            output = None
    with text(source) as stream:
        for line in stream:
            if line.startswith('#'):
                need(not jobs, 'Header after input records.')
                headers.append(line)
                continue
            fields = line.rstrip('\r\n').split('\t')
            need(len(fields) == 9 + len(samples), 'Four-sample VCF column contract changed.')
            info = dict((entry.split('=', 1) + [''])[:2] for entry in fields[7].split(';'))
            need('WGS_EVIDENCE' not in info, 'Core input already has packets.')
            expected = len(samples) * (len(info['CSQ'].split(',')) if info.get('CSQ') else 1)
            if output and rows >= rows_per_batch:
                finish()
            if output is None:
                path = directory / ('batch-' + str(len(jobs)).zfill(5) + '.vcf.gz')
                output = gzip.open(path, 'xt', encoding='utf-8', newline='', compresslevel=1)
                output.writelines(headers)
                jobs.append(dict(id=len(jobs), input=str(path), first_i=ordinal))
            output.write(line)
            records += 1
            rows += expected
    finish()
    need(jobs, 'No input records.')
    return headers, jobs, ordinal - 1


class PacketSink:
    def __init__(self, source, output, fields, first_i):
        self.source = text(source)
        self.raw_records = (line for line in self.source if not line.startswith('#'))
        self.output_file = Path(output).open('xb')
        self.bgzip = subprocess.Popen([G['bgzip'], '-c', '-@', '1'], stdin=subprocess.PIPE,
                                      stdout=self.output_file, stderr=subprocess.PIPE, preexec_fn=compression_limit)
        self.fields, self.ordinal = fields, first_i
        self.digest = VectorDigest()
        self.current, self.membership, self.rows = None, None, []
        self.record_count, self.row_count = 0, 0

    def begin_record(self, membership):
        self.flush()
        self.current = next(self.raw_records, None)
        need(self.current is not None, 'Frozen loop produced an extra record.')
        self.membership, self.rows = membership, []

    def writerow(self, row):
        need(self.current is not None and row['variant_id'] == self.membership['variant_id'], 'Frozen row/allele order changed.')
        values = ['' if row.get(field, '') is None else str(row.get(field, '')) for field in self.fields]
        self.digest.row(values)
        self.rows.append(dict(i=self.ordinal, values=values))
        self.ordinal += 1
        self.row_count += 1

    def flush(self):
        if self.current is None:
            return
        columns = self.current.rstrip('\r\n').split('\t')
        packet = dict(v=1, allele=[columns[0], int(columns[1]), columns[3], columns[4]],
                      scope='cohort', sample=None, rows=self.rows, membership=self.membership)
        token = encode(packet)
        # Decode the persisted representation before writing. Every complete
        # ordered vector and explicit membership must survive the packet codec.
        need(decode(token) == packet, 'Packet roundtrip changed a field vector.')
        original_info = columns[7]
        columns[7] = ('' if original_info == '.' else original_info + ';') + 'WGS_EVIDENCE=' + token
        ending = '\r\n' if self.current.endswith('\r\n') else '\n'
        self.bgzip.stdin.write(('\t'.join(columns) + ending).encode('utf-8'))
        columns[7] = original_info
        need('\t'.join(columns) + ending == self.current, 'An original VCF field changed.')
        self.record_count += 1
        self.current, self.rows = None, []

    def close(self):
        self.flush()
        need(next(self.raw_records, None) is None, 'Frozen loop omitted an input record.')
        self.source.close()
        self.bgzip.stdin.close()
        errors = self.bgzip.stderr.read()
        need(self.bgzip.wait() == 0, 'Isolated BGZF compression failed: ' + errors[:200].decode(errors='replace'))
        self.output_file.close()

    def abort(self):
        if self.bgzip.poll() is None:
            self.bgzip.terminate()
            self.bgzip.wait(timeout=10)
        self.source.close()
        self.output_file.close()


def initialize_worker():
    # A worker and its bgzip child share one CPU. These limits affect only the
    # new candidate processes, never the existing analysis.
    import resource
    available = sorted(os.sched_getaffinity(0))
    number = multiprocessing.current_process()._identity[-1] - 1
    os.sched_setaffinity(0, {available[number % len(available)]})
    resource.setrlimit(resource.RLIMIT_AS, (1536 * 1024**2, 1536 * 1024**2))
    m = G['stages']; alpha_module = m.alpha_module(); original_lookup = alpha_module.Lookup
    public_alpha = None
    if G.get('public_alpha'):
        files, manifest = G['public_alpha']
        public_alpha = original_lookup(sqlite3.connect(':memory:'), files, manifest)
    def database(manifest, manifest_path, destination):
        os.link(G['snapshot'], destination)  # Only an owned link to the new backup.
        return readonly_database(destination)
    def lookup(db, files, manifest):
        if public_alpha:
            return public_alpha
        value = object.__new__(original_lookup)
        value.db, value.manifest, value.source_rows = db, manifest, G['alpha_source_rows']
        return value
    m.make_adapter_database = database
    alpha_module.Lookup = lookup
    m.gene_evidence = lambda *_: G['genes']
    m.hpo_map = lambda *_: G['hpo']


def produce_batch(job):
    m = G['stages']; directory = Path(G['work']) / ('result-' + str(job['id']).zfill(5))
    directory.mkdir()
    sink = PacketSink(job['input'], directory / 'body.vcf.gz', m.FIELDS, job['first_i'])
    class Membership:
        writerow = lambda self, row: sink.begin_record(row)
    def writer(stack, path, fields):
        if path.name == 'transcript_samples.tsv.gz':
            need(fields == m.FIELDS, 'Frozen output field order changed.')
            return sink
        need(path.name == 'intersections.tsv.gz', 'Unexpected CSV sink; refusing silent omission.')
        return Membership()
    m.writer = writer
    started = time.monotonic()
    try:
        counts = m.make_tables(Path(job['input']), directory, G['config'], G['manifest'], G['manifest_path'], G['alpha_inputs'])
        sink.close()
    except BaseException:
        sink.abort()
        raise
    need(sink.row_count == job['expected_rows'] and sink.record_count == job['records']
         and sink.ordinal == job['last_i'] + 1, 'Batch record/row/ordinal count differs.')
    proof = dict(id=job['id'], first_i=job['first_i'], last_i=job['last_i'], row_count=sink.row_count,
                 record_count=sink.record_count, canonical_row_sha256=sink.digest.digest.hexdigest(),
                 body=pin(directory / 'body.vcf.gz'), seconds=time.monotonic() - started)
    save(directory / 'batch-proof.json', proof)
    save(directory / 'annotation_counts.json', counts)
    return dict(proof=proof, counts=counts, consequences=load(directory / 'consequence_counts.json'))


def merged_counts(results):
    first = results[0]['counts']
    counts = dict(first)
    for name in ('records', 'transcript_sample_rows', 'missing_csq', 'nonsequence_records_preserved'):
        counts[name] = sum(result['counts'][name] for result in results)
    counts['alphamissense_transcript_sample_statuses'] = {}
    for result in results:
        value = result['counts']
        need(value['samples'] == first['samples'] and value['alphamissense_matched_source_rows'] == first['alphamissense_matched_source_rows'],
             'Batch sample/Alpha source context differs; source_rows must not be summed.')
        for key, count in value['alphamissense_transcript_sample_statuses'].items():
            counts['alphamissense_transcript_sample_statuses'][key] = counts['alphamissense_transcript_sample_statuses'].get(key, 0) + count
    def add(left, right):
        if isinstance(right, int):
            return left + right
        value = dict(left)
        for key, item in right.items():
            value[key] = add(value[key], item) if key in value else item
        return value
    consequences = dict(results[0]['consequences'])
    for result in results[1:]:
        right = result['consequences']
        need(right['definitions'] == consequences['definitions'] and right['schema_version'] == consequences['schema_version'], 'Counter definitions changed.')
        for key in ('cohort', 'samples'):
            consequences[key] = add(consequences[key], right[key])
    return counts, consequences


def assemble(headers, results, output, fields, rows):
    proofs = [result['proof'] for result in results]
    expected = 1
    for proof in proofs:
        need(proof['first_i'] == expected and proof['last_i'] - proof['first_i'] + 1 == proof['row_count'], 'Batch ordinal gap/overlap.')
        expected = proof['last_i'] + 1
    need(expected - 1 == rows, 'Total ordinal count differs.')
    portable = [portable_batch(proof) for proof in proofs]
    header_digest = VectorDigest(); header_digest.row(fields)
    schema = dict(v=1, codec=CODEC, scope='cohort', sample=None, original_fields=fields,
        appended_fields=[], full_fields=fields, row_count=rows, ordinal_in_vcf_order=True,
        vector_digest_contract=DIGEST_CONTRACT, canonical_header_sha256=header_digest.digest.hexdigest(),
        batches=portable, ordered_batch_merkle_root=merkle_digest(fields, proofs),
        producer_kind='DIRECT_FROZEN_ANNOTATION_ROW_V1_NO_LEGACY_06_07_RECEIPTS')
    header = []
    for line in headers:
        if line.startswith('#CHROM'):
            header += ['##INFO=<ID=WGS_EVIDENCE,Number=1,Type=String,Description="Lossless frozen base annotation row vectors; decode WGS_EVIDENCE_SCHEMA">\n',
                       '##WGS_EVIDENCE_SCHEMA=<Version=1,Codec="' + CODEC + '",Data="' + encode(schema) + '">\n']
        header.append(line)
    header_file = output.parent / 'packet-header.bgzf'
    with header_file.open('xb') as stream:
        subprocess.run([G['bgzip'], '-c', '-@', '1'], input=''.join(header).encode(), stdout=stream, check=True)
    with output.open('xb') as out:
        for path in [header_file] + [Path(p['body']['path']) for p in proofs]:
            with path.open('rb') as stream:
                size = path.stat().st_size
                stream.seek(size - len(BGZF_EOF)); need(stream.read() == BGZF_EOF, 'Missing BGZF EOF block.')
                stream.seek(0); remaining = size - len(BGZF_EOF)
                while remaining:
                    block = stream.read(min(8 * 1024**2, remaining)); need(block, 'Truncated BGZF batch.')
                    out.write(block); remaining -= len(block)
        out.write(BGZF_EOF)
    subprocess.run([G['tabix'], '-p', 'vcf', str(output)], check=True, stdout=subprocess.DEVNULL)
    return schema


def prove(args):
    snapshot_manifest = load(args.snapshot_manifest)
    need(snapshot_manifest['status'] == 'CONSISTENT_READONLY_SNAPSHOT_CANDIDATE', 'Consistent candidate snapshot required.')
    binding_path = verify(snapshot_manifest['binding'])
    binding, queue, m, common = frozen(binding_path)
    public = load(args.public_proof / 'public-scatter-proof.json')
    need(public.get('status') == 'PASS' and public.get('patient_data_used') is False
         and public['binding_sha256'] == sha(binding_path), 'Matching sealed public controls required.')
    serial = args.public_proof / 'serial'
    for name in ('annotated.vcf.gz', 'transcript_samples.tsv.gz', 'annotation_counts.json', 'consequence_counts.json'):
        item = public['artifacts']['serial/' + name]
        need(sha(serial / name) == item['sha256'] and (serial / name).stat().st_size == item['bytes'], 'A sealed public baseline changed.')
    if not args.execute:
        return dict(status='PUBLIC_DIRECT_PACKET_PROOF_PLAN_ONLY', no_primary_mutation=True)
    deadline = deadline_guard(binding, queue)
    args.output.mkdir(mode=0o700, parents=True, exist_ok=False)
    snapshot_file = verify(snapshot_manifest['snapshot'])
    cfg = load(binding['annotation_config']['path']); manifest_path = Path(cfg['annotation']['resource_manifest'])
    manifest = m.load_resources(manifest_path, verify_tree=False, local_adapters=cfg['annotation'].get('local_adapters'))
    alpha_dir = args.output / 'public-alpha'; alpha_dir.mkdir()
    cfg['annotation']['input_vcf'] = str(serial / 'annotated.vcf.gz')
    cfg['output_root'] = str(args.output)
    alpha_inputs = m.alpha_module().prepare_matches(cfg, alpha_dir,
        common.Runtime(cfg, True, alpha_dir, 'annotate'), serial / 'annotated.vcf.gz')
    samples = m.vcf_samples(serial / 'annotated.vcf.gz')
    need(len(samples) == 4 and all(s.startswith('PUBLIC_') for s in samples), 'Public synthetic controls only.')
    headers, jobs, rows = partition(serial / 'annotated.vcf.gz', args.output / 'inputs', samples, args.rows_per_batch)
    need(len(jobs) >= 2, 'Public proof must exercise multiple balanced batches.')
    G.update(stages=m, snapshot=str(snapshot_file), alpha_source_rows=snapshot_manifest['alpha_source_rows'],
        work=str(args.output), config=cfg['annotation'], manifest=manifest, manifest_path=manifest_path,
        alpha_inputs=alpha_inputs, public_alpha=(alpha_inputs[1], alpha_inputs[0]),
        genes=m.gene_evidence(manifest, manifest_path), hpo=m.hpo_map(manifest, manifest_path),
        bgzip=shutil.which('bgzip'), tabix=shutil.which('tabix'))
    need(G['bgzip'] and G['tabix'], 'Existing bgzip/tabix required; the candidate installs nothing.')
    import resource
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
    started = time.monotonic()
    with multiprocessing.get_context('fork').Pool(args.workers, initializer=initialize_worker) as pool:
        results = pool.map(produce_batch, jobs)
    worker_seconds = time.monotonic() - started
    counts, consequences = merged_counts(results)
    need(counts == load(serial / 'annotation_counts.json'), 'Direct packet annotation counts differ from sealed make_tables.')
    need(consequences == load(serial / 'consequence_counts.json'), 'Direct packet consequence counters differ from sealed make_tables.')
    output = args.output / 'public.core-evidence.vcf.gz'
    schema = assemble(headers, results, output, m.FIELDS, rows)
    expected_ordinal, records = 1, 0
    with gzip.open(serial / 'transcript_samples.tsv.gz', 'rt', encoding='utf-8', newline='') as baseline, text(output) as vcf:
        expected = csv.reader(baseline, delimiter='\t'); need(next(expected) == m.FIELDS, 'Sealed baseline field schema differs.')
        for line in vcf:
            if line.startswith('#'):continue
            records += 1
            info = dict(v.split('=', 1) for v in line.rstrip('\r\n').split('\t')[7].split(';') if '=' in v)
            packet = decode(info['WGS_EVIDENCE'])
            for row in packet['rows']:
                need(row['i'] == expected_ordinal and row['values'] == next(expected, None), 'An ordered complete field vector differs from sealed make_tables.')
                expected_ordinal += 1
        need(next(expected, None) is None and expected_ordinal - 1 == rows and records == 20, 'Public full-vector/record count mismatch.')
    for chrom in sorted({line.split('\t', 1)[0] for line in text(serial / 'annotated.vcf.gz') if not line.startswith('#')}):
        subprocess.run([G['tabix'], str(output), chrom], stdout=subprocess.DEVNULL, check=True)
    save(args.output / 'annotation_counts.json', counts)
    save(args.output / 'consequence_counts.json', consequences)
    proof = dict(schema_version=1, status='PASS_PUBLIC_DIRECT_PACKET_PROOF', patient_data_used=False,
        public_variants=20, full_ordered_field_vectors=rows, batches=len(jobs), workers=args.workers,
        one_cpu_per_worker=True, worker_address_limit_bytes=1536 * 1024**2,
        parent_plus_worker_allocation_ceiling_bytes=16 * 1024**3, compression_child_limit_bytes=128 * 1024**2,
        deadline_utc=deadline.isoformat(), balanced_target_rows=args.rows_per_batch,
        direct_worker_wall_seconds=worker_seconds, every_field_value_and_order_equals_sealed_make_tables=True,
        annotation_counts_equal=True, consequence_counters_equal=True, bgzip_TBI_random_access=True,
        baseline_proof=pin(args.public_proof / 'public-scatter-proof.json'), snapshot_manifest=pin(args.snapshot_manifest),
        producer=pin(Path(__file__)), frozen_package=binding['code_package'], schema=schema,
        output=pin(output), index=pin(Path(str(output) + '.tbi')),
        no_legacy_stage06_or_07_receipt_created=True, no_patient_pipeline_or_controller_actuation=True)
    save(args.output / 'public-direct-proof.json', proof)
    return proof


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    snap = sub.add_parser('snapshot')
    snap.add_argument('--bindings', type=Path, required=True)
    snap.add_argument('--source-dir', type=Path, required=True)
    snap.add_argument('--output', type=Path, required=True)
    snap.add_argument('--execute', action='store_true')
    proof = sub.add_parser('prove')
    proof.add_argument('--snapshot-manifest', type=Path, required=True)
    proof.add_argument('--public-proof', type=Path, required=True)
    proof.add_argument('--output', type=Path, required=True)
    proof.add_argument('--workers', type=int, default=8, choices=range(1, 9))
    proof.add_argument('--rows-per-batch', type=int, default=250000)
    proof.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    need(args.command != 'prove' or 1 <= args.rows_per_batch <= 250000, 'Bounded balanced row target required.')
    os.umask(0o077)
    print(json.dumps({'snapshot': snapshot, 'prove': prove}[args.command](args), indent=2))


if __name__ == '__main__':
    main()
