#!/usr/bin/env python3
"""Reviewed offline supplementary annotation, portable VCF evidence and functional views."""
import argparse
import copy
import csv
import datetime as dt
import fcntl
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

END = dt.datetime.fromisoformat('2026-11-06T06:30:00+00:00')
TAG = 'WGS_SUPPLEMENTARY_EVIDENCE'
SCHEMA = 'WGS_SUPPLEMENTARY_SCHEMA'


def need(value, message):
    if not value:
        raise ValueError(message)


def now():
    return dt.datetime.now(dt.timezone.utc)


def safe(path):
    path = Path(path).absolute()
    need(not any(p.is_symlink() for p in [path, *path.parents]), 'Symlink source/output refused')
    return path


def sha(path):
    value = hashlib.sha256()
    with safe(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            value.update(block)
    return value.hexdigest()


def pin(path):
    path = safe(path)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def verify(value):
    path = safe(value['path'])
    need(path.stat().st_size == value['bytes'] and sha(path) == value['sha256'], 'Pinned artifact changed')
    return path


def load(path):
    return json.loads(safe(path).read_text(encoding='utf-8-sig'))


def write(path, value, replace=False):
    path = safe(path)
    if replace:
        temp = path.with_suffix(path.suffix + '.tmp')
        temp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
        os.replace(temp, path)
    else:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(value, stream, indent=2, ensure_ascii=False)
            stream.write('\n')


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def configure(base, source, output, forks=8):
    cfg = copy.deepcopy(base)
    cfg['output_root'] = str(output)
    cfg['annotation'].update(input_vcf=str(source), forks=forks, memory_gib=32,
                             timeout_seconds=7200, require_alphamissense=True)
    cfg.setdefault('resources', {})['min_free_gib'] = 150
    # Original budget stays byte-for-byte equivalent in the copied configuration.
    return cfg


def upstream_pin(value, root, filename):
    artifacts = value.get('artifacts', {})
    if isinstance(artifacts, dict):
        candidates = [dict(item, path=item.get('path', str(root / name)))
                      for name, item in artifacts.items() if Path(name).name == filename]
    elif isinstance(artifacts, list):
        candidates = [item for item in artifacts if Path(item.get('path', '')).name == filename]
    else:
        raise ValueError('Unsupported upstream artifact schema')
    need(len(candidates) == 1, 'Exactly one completed upstream artifact required: ' + filename)
    item = candidates[0]
    path = Path(item['path'])
    path = path if path.is_absolute() else root / path
    need(path.resolve().is_relative_to(root.resolve()), 'Upstream artifact escapes reviewed root')
    item = dict(item, path=str(path))
    return dict(path=item['path'], bytes=item['bytes'], sha256=item['sha256'])


def checked_upstream(profile):
    root = safe(profile['upstream_output'])
    receipt_path, manifest_path = root / 'workflow-receipt.json', root / 'manifest.json'
    if not receipt_path.exists() or not manifest_path.exists():
        return None
    receipt, manifest = load(receipt_path), load(manifest_path)
    need(receipt.get('status') == 'PASS' and receipt.get('stage') == profile['upstream_stage'], 'Upstream completion stage/status differs')
    need(manifest.get('status') == profile['upstream_manifest_status'], 'Upstream final manifest not PASS')
    need(manifest.get('profile_sha256') == profile['upstream_profile_identity'], 'Caller canonical profile lineage differs')
    need(manifest.get('producer_sha256') == profile['upstream_executor']['sha256'], 'Caller producer lineage differs')
    stored_profile = manifest.get('source_profile', {})
    need(safe(stored_profile['path']).resolve().is_relative_to(root.resolve()), 'Stored caller profile escapes reviewed output')
    verify(stored_profile)
    need(load(stored_profile['path']) == load(profile['upstream_profile']['path']), 'Stored caller profile differs from reviewed source')
    verify(upstream_pin(receipt, root, 'manifest.json'))
    source = upstream_pin(manifest, root, profile['upstream_vcf_name'])
    index = upstream_pin(manifest, root, profile['upstream_vcf_name'] + '.tbi')
    for item in [source, index]:
        verify(item)
        need(upstream_pin(receipt, root, Path(item['path']).name) == item, 'Caller manifest/receipt artifact differs')
    return dict(vcf=source, index=index, manifest=pin(manifest_path), receipt=pin(receipt_path))


def text(path):
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if str(path).endswith('.gz') else Path(path).open(encoding='utf-8', newline='')


def header(stream):
    lines = []
    for line in stream:
        need(line.startswith('#'), 'VCF sample header absent')
        lines.append(line.rstrip('\r\n'))
        if line.startswith('#CHROM'):
            return lines
    raise ValueError('VCF header absent')


def key(fields):
    return fields[0].removeprefix('chr'), int(fields[1]), fields[3], fields[4]


def row_key(row):
    return row['chrom'].removeprefix('chr'), int(row['pos']), row['ref'], row['alt']


def info_tokens(value):
    return [] if value == '.' else value.split(';')


def packet_header(codec, fields, scope='cohort', sample=None):
    schema = dict(v=1, fields=fields, codec=codec.CODEC, scope=scope,
                  rows='Ordered global ordinal i plus complete string values vector', sample=sample)
    return '##' + SCHEMA + '=<Version=1,Data="' + codec.encode_packet(schema) + '">'


def embed(source, annotated, table, output, codec, guard):
    import pysam
    csv.field_size_limit(64 * 1024**2)
    rows_total = records = 0
    original_digest, decoded_digest = hashlib.sha256(), hashlib.sha256()
    with text(source) as src, text(annotated) as ann, text(table) as tab, pysam.BGZFile(str(output), 'wb') as out:
        original_header, annotation_header = header(src), header(ann)
        need(original_header[-1] == annotation_header[-1], 'Annotation sample/header order changed')
        need(not any('ID=' + TAG + ',' in line or line.startswith('##' + SCHEMA + '=') for line in original_header), 'Supplementary namespace already present')
        reader = csv.DictReader(tab, delimiter='\t')
        fields = reader.fieldnames
        need(fields and len(set(fields)) == len(fields) and 'sample' in fields, 'Unique complete master fields required')
        original_digest.update(codec.canonical_tsv(fields))
        out_header = list(original_header[:-1])
        out_header += [line for line in annotation_header[:-1] if line not in original_header and not line.startswith('##fileformat=')]
        out_header += ['##INFO=<ID=' + TAG + ',Number=1,Type=String,Description="Lossless ordered supplementary annotation rows: JSON UTF8 zlib base64url; bounded64MiB plaintext; schema in ' + SCHEMA + '">', packet_header(codec, fields), original_header[-1]]
        out.write(('\n'.join(out_header) + '\n').encode())
        pending = next(reader, None)
        for raw_source in src:
            guard()
            raw_annotation = next(ann, None)
            need(raw_annotation is not None, 'Annotation lost a source record')
            original, added = raw_source.rstrip('\r\n').split('\t'), raw_annotation.rstrip('\r\n').split('\t')
            need(key(original) == key(added) and original[8:] == added[8:], 'Annotation changed allele/GT/FORMAT')
            old_tokens = info_tokens(original[7]); names = {t.split('=', 1)[0]: t for t in old_tokens}
            extra = []
            for token in info_tokens(added[7]):
                name = token.split('=', 1)[0]
                if name in names:
                    need(token == names[name], 'Annotation changed original INFO value')
                else:
                    extra.append(token)
            vectors, plain_bytes = [], 0
            while pending is not None and row_key(pending) == key(original):
                values = [pending[f] for f in fields]
                need(all(isinstance(value, str) for value in values), 'Master vector contains missing/non-string fields')
                plain_bytes += len(json.dumps(values, ensure_ascii=False).encode())
                need(plain_bytes <= 64 * 1024**2, 'One allele exceeds64MiB; no truncation')
                vectors.append(dict(i=rows_total, values=values))
                original_digest.update(codec.canonical_tsv(values)); rows_total += 1
                pending = next(reader, None)
            need(vectors, 'Every source allele requires actual annotation rows')
            packet = dict(v=1, scope='cohort', allele=list(key(original)), rows=vectors)
            token = codec.encode_packet(packet)
            need(codec.decode_packet(token) == packet, 'Evidence codec roundtrip differs')
            final = list(original)
            final[7] = ';'.join(old_tokens + extra + [TAG + '=' + token])
            out.write(('\t'.join(final) + '\n').encode()); records += 1
        need(next(ann, None) is None and pending is None, 'Unexpected extra annotation records/rows')
    guard(); pysam.tabix_index(str(output), preset='vcf'); guard()
    # Read persisted BGZF and reconstruct the actual master used by all views.
    extracted = output.parent / 'master.from-vcf.tsv.gz'
    with text(output) as stream, text(source) as check_source, gzip.open(extracted, 'wt', encoding='utf-8', newline='') as master:
        persisted_header = header(stream)
        source_header = header(check_source)
        need(persisted_header[:len(original_header)-1] == original_header[:-1], 'Original header sequence changed')
        writer = csv.writer(master, delimiter='\t', lineterminator='\n'); writer.writerow(fields)
        decoded_digest.update(codec.canonical_tsv(fields)); expected_i = 0; persisted_records = 0
        for line in stream:
            guard(); record = line.rstrip('\r\n').split('\t')
            original_line = next(check_source, None)
            need(original_line is not None, 'Persisted VCF has an extra record')
            original_record = original_line.rstrip('\r\n').split('\t')
            need(record[:7] == original_record[:7] and record[8:] == original_record[8:], 'Persisted original record fields/GT/FORMAT differ')
            original_info = info_tokens(original_record[7])
            need(info_tokens(record[7])[:len(original_info)] == original_info, 'Persisted original INFO tokens/order differ')
            encoded = [t.split('=', 1)[1] for t in record[7].split(';') if t.startswith(TAG + '=')]
            need(len(encoded) == 1, 'Exactly one complete persisted evidence packet required')
            packet = codec.decode_packet(encoded[0])
            need(tuple(packet['allele']) == key(record), 'Persisted packet allele differs')
            for row in packet['rows']:
                need(row['i'] == expected_i and len(row['values']) == len(fields), 'Row ordinal/field vector differs')
                writer.writerow(row['values']); decoded_digest.update(codec.canonical_tsv(row['values'])); expected_i += 1
            persisted_records += 1
        need(next(check_source, None) is None and expected_i == rows_total and persisted_records == records and original_digest.digest() == decoded_digest.digest(), 'Persisted complete master vector digest differs')
    return dict(records=records, rows=rows_total, fields=fields, ordered_vector_sha256=original_digest.hexdigest(),
                actual_persisted_packet_decode_verified=True, original_GT_FORMAT_INFO_header_preserved=True,
                samples=original_header[-1].split('\t')[9:], extracted_master=pin(extracted))


def prune_sample(path, codec, fields, guard, all_samples=()):
    import pysam
    temp = path.with_name(path.name + '.private-rewrite.gz')
    with text(path) as stream, pysam.BGZFile(str(temp), 'wb') as output:
        headers = header(stream); samples = headers[-1].split('\t')[9:]
        need(len(samples) == 1, 'Private subset requires exactly one sample')
        sample = samples[0]; sample_index = fields.index('sample')
        # No cohort row metadata survives in a private sample schema.
        public_prefixes = ('##fileformat=', '##contig=', '##INFO=', '##FORMAT=', '##FILTER=', '##ALT=', '##reference=', '##VEP=')
        headers = [line for line in headers[:-1] if line.startswith(public_prefixes)
                   and not any(other in line for other in all_samples if other != sample)] + [packet_header(codec, fields, 'sample', sample), headers[-1]]
        output.write(('\n'.join(headers) + '\n').encode())
        for line in stream:
            guard(); row = line.rstrip('\r\n').split('\t'); tokens = info_tokens(row[7])
            encoded = [t.split('=', 1)[1] for t in tokens if t.startswith(TAG + '=')]
            need(len(encoded) == 1, 'Sample missing portable evidence')
            packet = codec.decode_packet(encoded[0])
            selected = [r for r in packet['rows'] if r['values'][sample_index] == sample]
            need(selected and all(r['values'][sample_index] == sample for r in selected), 'Private rows absent or mixed')
            clean = dict(v=1, scope='sample', sample=sample, allele=packet['allele'], rows=selected)
            row[7] = ';'.join([t for t in tokens if not t.startswith(TAG + '=')] + [TAG + '=' + codec.encode_packet(clean)])
            output.write(('\t'.join(row) + '\n').encode())
    os.replace(temp, path); Path(str(path) + '.tbi').unlink(); pysam.tabix_index(str(path), preset='vcf')
    with text(path) as stream:
        header(stream)
        for line in stream:
            token = next(t.split('=', 1)[1] for t in line.split('\t')[7].split(';') if t.startswith(TAG + '='))
            packet = codec.decode_packet(token)
            need(packet['scope'] == 'sample' and packet['sample'] == sample and all(r['values'][sample_index] == sample for r in packet['rows']), 'Persisted private packet leaks another sample')


def finalize_views(folder, provenance, guard):
    """Declare only the actual post-pruning bytes; the old engine stays unchanged."""
    folder = safe(folder)
    manifest = load(folder/'manifest.json')
    need(manifest.get('status') == 'PASS', 'Completed engine views required before private packet finalization')
    manifest.update(source_authority='ACTUAL_SUPPLEMENTARY_ANNOTATION_PASS_AND_VCF_DECODED_MASTER',
                    packet_sample_privacy_verified=True, no_legacy06_07_authority=True,
                    supporting_provenance=provenance)
    artifacts = {}
    for path in sorted(folder.rglob('*')):
        if path.is_file() and path.name not in ['manifest.json', 'workflow-receipt.json']:
            guard(); item = pin(path)
            artifacts[str(path.relative_to(folder))] = dict(bytes=item['bytes'], sha256=item['sha256'])
    manifest['artifacts'] = artifacts
    write(folder/'manifest.json', manifest, replace=True)
    receipt_artifacts = {}
    for path in sorted(folder.rglob('*')):
        if path.is_file() and path.name != 'workflow-receipt.json':
            guard(); receipt_artifacts[str(path.relative_to(folder))] = pin(path)
    receipt = dict(schema_version=1, status='PASS', stage='supplementary-functional-views-v1',
                   supporting_provenance=provenance, packet_sample_privacy_verified=True,
                   no_legacy06_07_authority=True, artifacts=receipt_artifacts)
    write(folder/'workflow-receipt.json', receipt, replace=True)
    for name, item in manifest['artifacts'].items():
        guard(); actual = pin(folder/name)
        need(actual['bytes'] == item['bytes'] and actual['sha256'] == item['sha256']
             and actual == receipt['artifacts'][name], 'Post-pruning manifest/receipt artifact differs')
    need(receipt['artifacts']['manifest.json'] == pin(folder/'manifest.json'), 'Final manifest not pinned by receipt')
    return dict(manifest=pin(folder/'manifest.json'), receipt=pin(folder/'workflow-receipt.json'),
                artifacts_verified=len(manifest['artifacts']))


def prepare(args):
    deadline = dt.datetime.fromisoformat(args.deadline_utc.replace('Z', '+00:00'))
    need(deadline.tzinfo is not None and now() < deadline <= END, 'Explicit user-authorized deadline must be<=06:30UTC')
    runtime = safe(args.runtime_scripts)
    modules = {name: pin(runtime / name) for name in ['annotation_stages.py', 'workflow_common.py', 'alphamissense_adapter.py']}
    cfg = load(args.base_config)
    caller_profile = load(args.upstream_profile)
    caller_identity = hashlib.sha256(json.dumps(caller_profile, sort_keys=True).encode()).hexdigest()
    verify(caller_profile['executor'])
    alpha = Path(cfg['annotation'].get('alphamissense_manifest') or Path(cfg['data_root']) / 'resources/annotation-supplements/alphamissense-2023-v2-grch38-v3/resource-manifest.json')
    need(not safe(args.output).exists() and not safe(args.views_output).exists(), 'Fresh output roots required')
    profile = dict(schema_version=1, status='CANDIDATE_SUPPLEMENTARY_ANNOTATION_NOT_RUNNING', root_reviewed=False,
                   producer=pin(__file__), base_config=pin(args.base_config), runtime_modules=modules,
                   codec=pin(args.codec), views_engine=pin(args.views_engine), resource_manifest=pin(cfg['annotation']['resource_manifest']), alphamissense_manifest=pin(alpha),
                   upstream_profile=pin(args.upstream_profile), upstream_profile_identity=caller_identity,
                   upstream_executor=caller_profile['executor'], upstream_output=str(safe(args.upstream_output)),
                   upstream_stage='xy-germline-v1', upstream_manifest_status='PASS_XY_GERMLINE_CALLING', upstream_vcf_name='cohort.xy.soft-filtered.vcf.gz',
                   output=str(safe(args.output)), views_output=str(safe(args.views_output)),
                   state_directory=str(safe(args.profile).parent / 'supplementary-annotation-state-v1'),
                   deadline_utc=deadline.isoformat(), forks=8, memory_gib=32, min_free_gib=150,
                   scope='USER_AUTHORIZED_SUPPLEMENTARY_XY_ANNOTATION; no clinical karyotype/automaticACMG; original autosomal outputs untouched')
    safe(args.profile).parent.mkdir(parents=True, exist_ok=True)
    write(args.profile, profile); print(json.dumps(dict(status=profile['status'], profile=pin(args.profile)), indent=2))


def run(args):
    need(sha(args.profile) == args.reviewed_profile_sha256, 'Exact reviewed profile SHA required')
    p = load(args.profile)
    if not args.execute:
        print('PLAN_ONLY_NO_ANNOTATION'); return
    need(p['root_reviewed'] is True and p['status'] == 'ROOT_REVIEWED_SUPPLEMENTARY_ANNOTATION', 'ROOT-reviewed exact execution profile required')
    deadline = dt.datetime.fromisoformat(p['deadline_utc'])
    need(now() < deadline <= END and p['forks'] == 8 and p['memory_gib'] == 32 and p['min_free_gib'] >= 150, 'Deadline/resource caps changed')
    need(verify(p['producer']) == Path(__file__).absolute(), 'Producer path differs')
    for key_name in ['base_config', 'codec', 'views_engine', 'resource_manifest', 'alphamissense_manifest', 'upstream_profile', 'upstream_executor']:
        verify(p[key_name])
    for value in p['runtime_modules'].values():
        verify(value)
    state = safe(p['state_directory']); state.mkdir(mode=0o700, parents=True, exist_ok=True)
    lock = (state / '.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    def progress(phase, **counts):
        write(state / 'status.json', dict(status=phase, pid=os.getpid(), updated_utc=now().isoformat(), deadline_utc=p['deadline_utc'], **counts), replace=True)
    def guard():
        need(now() < deadline, 'Supplementary annotation deadline reached')
        need(shutil.disk_usage(state).free >= p['min_free_gib'] * 1024**3, 'Disk reserve150GiB reached')
    def expired(*unused):
        raise TimeoutError('User-authorized supplementary annotation deadline reached')
    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, (deadline-now()).total_seconds())
    os.sched_setaffinity(0, set(sorted(os.sched_getaffinity(0))[:8]))
    resource.setrlimit(resource.RLIMIT_AS, (32*1024**3, 32*1024**3))
    common = module(p['runtime_modules']['workflow_common.py']['path'], '_supplementary_runtime')
    annotation = module(p['runtime_modules']['annotation_stages.py']['path'], '_supplementary_annotation')
    codec = module(p['codec']['path'], '_supplementary_evidence_codec')
    class DeadlineRuntime(common.Runtime):
        def remaining_seconds(self):
            guard(); return (deadline-now()).total_seconds()
        def docker(self, image, arguments, **kwargs):
            need(kwargs.get('cpus', 4) <= 8 and kwargs.get('memory_gib', 16) <= 32 and not kwargs.get('gpu'), 'Annotation runtime exceeds reviewed limits')
            return super().docker(image, arguments, **kwargs)
        def run(self, argv, logname='command', capture=False, cwd=None, stdout=None, stderr=None, timeout=None):
            allowed = min(self.remaining_seconds(), timeout or 7200)
            self.commands.append(list(map(str, argv)))
            container = argv[argv.index('--name') + 1] if '--name' in argv else None
            handles = []
            proc = None
            try:
                log = (self.workdir / (str(logname) + '.log')).open('ab'); handles.append(log)
                if isinstance(stdout, (str, Path)):
                    output = Path(stdout).open('wb'); handles.append(output)
                else:
                    output = stdout if stdout is not None else subprocess.PIPE if capture else log
                if isinstance(stderr, (str, Path)):
                    errors = Path(stderr).open('ab'); handles.append(errors)
                else:
                    errors = stderr if stderr is not None else log
                proc = subprocess.Popen(list(map(str, argv)), cwd=cwd, stdin=subprocess.DEVNULL, stdout=output, stderr=errors, start_new_session=True)
                data, _ = proc.communicate(timeout=allowed)
                need(proc.returncode == 0, 'Offline annotation command failed; inspect private log')
                guard(); return data.decode() if capture and data else ''
            except BaseException:
                if proc is not None and proc.poll() is None:
                    os.killpg(proc.pid, signal.SIGKILL); proc.wait(timeout=5)
                if container:
                    subprocess.run(self.docker_prefix() + ['rm', '-f', container], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=10)
                raise
            finally:
                for stream in handles:
                    stream.close()
    output, views_output = safe(p['output']), safe(p['views_output'])
    need(not output.exists() and not views_output.exists(), 'Existing supplementary output refused')
    try:
        source = checked_upstream(p)
        while source is None:
            guard(); progress('WAITING_FOR_XY_CALLING_PASS'); time.sleep(min(15, max(0, (deadline-now()).total_seconds())))
            source = checked_upstream(p)
        workdir = output.with_name(output.name + '.unpromoted-' + str(os.getpid()))
        workdir.mkdir(mode=0o700, parents=True, exist_ok=False)
        cfg = configure(load(p['base_config']['path']), source['vcf']['path'], output)
        runtime = DeadlineRuntime(cfg, True, workdir, 'annotate')
        progress('RUNNING_OFFLINE_VEP_AND_FULL_RESOURCE_VERIFICATION')
        result = annotation.annotate(cfg, workdir, runtime)
        progress('EMBEDDING_AND_DECODING_COMPLETE_TRANSCRIPT_EVIDENCE')
        portable = workdir / 'cohort.supplementary.annotated.vcf.gz'
        proof = embed(source['vcf']['path'], workdir/'annotated.vcf.gz', workdir/'transcript_samples.tsv.gz', portable, codec, guard)
        proof['extracted_master']['path'] = str(output/'master.from-vcf.tsv.gz')
        manifest = dict(schema_version=1, status='PASS_SUPPLEMENTARY_ANNOTATION', source=source,
                        profile_sha256=sha(args.profile), producer=p['producer'], modules=p['runtime_modules'], codec=p['codec'],
                        packet_proof=proof, annotation=result['metadata'], no_patient_variant_remote_queries=True,
                        original_budget_config_unchanged=True, standalone_stage_not_legacy06_07=True)
        write(workdir/'manifest.json', manifest)
        artifacts = {f.name: pin(f) for f in workdir.iterdir() if f.is_file()}
        # Relocate artifact paths in the promotion receipt; bytes/hashes stay exact.
        artifacts = {name: dict(item, path=str(output/name)) for name, item in artifacts.items()}
        write(workdir/'workflow-receipt.json', dict(schema_version=1, status='PASS', stage='supplementary-annotation-v1', artifacts=artifacts))
        need(not output.exists(), 'Output appeared during annotation'); os.rename(workdir, output)
        guard(); progress('EXPORTING_ALL_FUNCTIONAL_COHORT_AND_SAMPLE_VIEWS')
        engine = module(p['views_engine']['path'], '_supplementary_views')
        engine.ABSOLUTE_DEADLINE = deadline
        def actual_annotation_inputs(unused, view_guard):
            receipt = load(output/'workflow-receipt.json')
            need(receipt['status'] == 'PASS' and receipt['stage'] == 'supplementary-annotation-v1', 'Actual supplementary annotation PASS required')
            sources = dict(vcf=output/portable.name, vcf_index=output/(portable.name+'.tbi'), master=output/'master.from-vcf.tsv.gz')
            pins = {str(output/'workflow-receipt.json'): pin(output/'workflow-receipt.json')}
            for path in sources.values():
                item = receipt['artifacts'][path.name]; verify(item); pins[str(path)] = item
            return sources, pins
        engine.checked_inputs = actual_annotation_inputs
        engine.export(output, views_output, p['runtime_modules']['annotation_stages.py']['path'], p['runtime_modules']['annotation_stages.py']['sha256'], deadline=deadline)
        progress('PRUNING_SAMPLE_PACKETS_TO_SELECTED_SAMPLE_ONLY')
        for path in (views_output/'vcf').glob('*.vcf.gz'):
            prune_sample(path, codec, proof['fields'], guard, proof['samples'])
        finalize_views(views_output, dict(producer=p['producer'], views_engine=p['views_engine'],
                       reviewed_profile=pin(args.profile), actual_annotation_receipt=pin(output/'workflow-receipt.json'),
                       evidence_master_origin='master.from-vcf.tsv.gz decoded from actual persisted supplementary VCF'), guard)
        progress('PASS_SUPPLEMENTARY_ANNOTATION_AND_ALL_FUNCTIONAL_VIEWS', vcf_records=proof['records'], transcript_rows=proof['rows'])
    except BaseException as error:
        progress('FAILED_OR_BOUNDED_STOP', error_type=type(error).__name__)
        raise


def main():
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('prepare')
    for name in ['base-config','runtime-scripts','codec','views-engine','upstream-profile','upstream-output','output','views-output','profile']:
        prep.add_argument('--'+name, type=Path, required=True)
    prep.add_argument('--deadline-utc', required=True)
    action = sub.add_parser('run'); action.add_argument('--profile', type=Path, required=True)
    action.add_argument('--reviewed-profile-sha256', required=True); action.add_argument('--execute', action='store_true')
    args = parser.parse_args(); {'prepare': prepare, 'run': run}[args.command](args)


if __name__ == '__main__':
    main()
