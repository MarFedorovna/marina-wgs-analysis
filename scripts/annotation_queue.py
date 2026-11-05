#!/usr/bin/env python3
"""Detached, hash-bound continuation of accepted stages 5 -> public controls -> 6/7."""
import argparse
import copy
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import time
import uuid

CODE = Path(__file__).resolve().parents[1]
PACKAGE_FILES = [
    'scripts/workflow.py', 'scripts/workflow_common.py', 'scripts/calling_stages.py',
    'scripts/joint_scatter.py', 'scripts/configure_workflow.py', 'scripts/annotation_stages.py',
    'scripts/alphamissense_adapter.py', 'scripts/plot_report.py', 'scripts/export_report_docx.py',
    'scripts/validate_annotation_controls.py', 'scripts/annotation_queue.py',
    'config/evidence-sources.example.json', 'config/annotation-resources.example.json',
    'config/annotation-base-resources.example.json', 'requirements-report.txt',
]
STOP = False


def need(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 ** 2), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def encoded(value):
    return (json.dumps(value, indent=2, ensure_ascii=False) + '\n').encode()


def exclusive(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    need(not path.is_symlink() and (not path.exists() or path.read_bytes() == content), 'Refusing to overwrite a changed sealed artifact.')
    if not path.exists():
        with path.open('xb') as stream:
            stream.write(content)


def write_status(path, value):
    temporary = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    temporary.write_bytes(encoded(value))
    os.replace(temporary, path)


def budget_deadline(created, hourly, cap):
    origin = dt.datetime.fromisoformat(created.replace('Z', '+00:00'))
    need(origin.tzinfo is not None and 0 < hourly and 0 < cap <= 40, 'A timezone-qualified actual VM start and unchanged allowance <=$40 are required.')
    return origin.astimezone(dt.timezone.utc) + dt.timedelta(hours=cap / hourly)


def verify_stage(directory, expected_receipt_sha=None):
    receipt_path = directory / 'workflow-receipt.json'
    if not receipt_path.is_file():
        return False
    if expected_receipt_sha:
        need(digest(receipt_path) == expected_receipt_sha, 'Accepted upstream receipt changed.')
    receipt = load(receipt_path)
    need(receipt.get('status') == 'PASS' and receipt.get('artifacts'), 'Stage lacks a nonempty PASS receipt.')
    for name, item in receipt['artifacts'].items():
        relative = Path(name)
        path = directory / relative
        need(not relative.is_absolute() and '..' not in relative.parts and not path.is_symlink()
             and path.resolve(strict=True).is_relative_to(directory.resolve())
             and path.stat().st_size == item['bytes'] and digest(path) == item['sha256'],
             'Completed stage artifact differs from its accepted receipt.')
    return True


def prepare(args):
    os.umask(0o077)
    base_path = args.base_config.resolve(strict=True)
    base = load(base_path)
    root = Path(base['data_root']).resolve(strict=True)
    target = args.runtime_root.resolve()
    expected_filter = Path(base['output_root']) / 'stages/05_filter'
    need(verify_stage(expected_filter), 'Only an accepted PASS stage5 can be queued.')
    need(digest(Path(base['reference']['fasta'])) == base['reference']['sha256'], 'Fixed analysis FASTA changed.')
    target.mkdir(parents=True, exist_ok=True)
    code_files = []
    for relative in PACKAGE_FILES:
        source = CODE / relative
        need(source.is_file() and not source.is_symlink(), 'Missing reviewed continuation source: ' + relative)
        data = source.read_bytes()
        exclusive(target / relative, data)
        code_files.append(dict(path=relative, bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    package = dict(schema_version=1, files=code_files, purpose='reviewed research annotation continuation; no patient inputs in package')
    package_path = target / 'code-manifest.json'
    exclusive(package_path, encoded(package))
    package_sha = digest(package_path)
    alpha_path = root / 'resources/annotation-supplements/alphamissense-2023-v2-grch38-v3/resource-manifest.json'
    need(load(alpha_path).get('status') == 'PASS', 'Mandatory AlphaMissense supplement is not ready.')
    resources = Path(base['annotation']['resource_manifest'])
    final_recipe = CODE / 'config/annotation-resources.example.json'
    recipe = load(final_recipe)
    need(resources.parent.name == recipe['bundle_id'], 'Resource target differs from the reviewed v2 recipe.')
    annotation = copy.deepcopy(base)
    annotation['tools']['report_python'] = str(root / 'tools/report-venv/bin/python')
    annotation['annotation'].update(alphamissense_manifest=str(alpha_path), alphamissense_manifest_sha256=digest(alpha_path),
        require_alphamissense=True, disease_context={'mondo_ids': ['MONDO:0005045'],
            'description': 'HCM research context explicitly requested; not a patient diagnosis'}, pathway_mode='disabled')
    annotation['budget'].update(vm_created_utc=args.created_utc, hourly_usd=2.39, planned_vm_usd=40, user_budget_usd=50)
    annotation['continuation_provenance'] = dict(base_config_sha256=digest(base_path),
        accepted_filter_receipt_sha256=digest(expected_filter / 'workflow-receipt.json'),
        code_package_sha256=package_sha, upstream_stages_not_rerun=True, clinical_diagnosis=False)
    annotation_path = root / 'private/annotation-v1.json'
    exclusive(annotation_path, encoded(annotation))
    external = {}
    for name in ('workflow.py', 'workflow_common.py', 'calling_stages.py', 'configure_workflow.py', 'joint_scatter.py'):
        path = args.upstream_code_root.resolve(strict=True) / 'scripts' / name
        need(path.is_file(), 'Missing frozen upstream source file.')
        external[str(path)] = digest(path)
    bindings = dict(schema_version=1, status='SEALED_FOR_ROOT_REVIEW', data_root=str(root),
        base_config={'path': str(base_path), 'sha256': digest(base_path)},
        annotation_config={'path': str(annotation_path), 'sha256': digest(annotation_path)},
        accepted_filter={'path': str(expected_filter), 'receipt_sha256': digest(expected_filter / 'workflow-receipt.json')},
        reference=base['reference'], code_package={'path': str(package_path), 'sha256': package_sha},
        external_upstream_code=external, alpha={'path': str(alpha_path), 'sha256': digest(alpha_path)},
        resources={'path': str(resources), 'expected_recipe_sha256': digest(final_recipe), 'vep_image': base['tools']['vep_image']},
        budget={'vm_created_utc': args.created_utc, 'hourly_usd': 2.39, 'planned_vm_usd': 40},
        output_root=base['output_root'], approval_file=str(root / 'private/annotation-queue-approved.json'),
        status_file=str(root / 'private/annotation-queue-status.json'),
        log_dir=str(root / 'logs/annotation-continuation'),
        allowed_changes=['annotation only', 'tools.report_python', 'VM creation budget anchor', 'continuation provenance'],
        patient_diagnosis_automatic=False)
    bindings_path = root / 'private/annotation-queue-bindings.json'
    exclusive(bindings_path, encoded(bindings))
    print(json.dumps({'status': 'SEALED_FOR_ROOT_REVIEW', 'bindings': str(bindings_path), 'bindings_sha256': digest(bindings_path),
        'package_sha256': package_sha, 'annotation_config_sha256': digest(annotation_path), 'upstream_stages_unchanged': True}))


def verify_bindings(bindings):
    for key in ('base_config', 'annotation_config', 'alpha'):
        item = bindings[key]
        need(digest(Path(item['path'])) == item['sha256'], 'A sealed config/Alpha source changed: ' + key)
    package_path = Path(bindings['code_package']['path'])
    need(digest(package_path) == bindings['code_package']['sha256'], 'Continuation code manifest changed.')
    for item in load(package_path)['files']:
        path = package_path.parent / item['path']
        need(path.resolve(strict=True).is_relative_to(package_path.parent.resolve()) and not path.is_symlink()
             and path.stat().st_size == item['bytes'] and digest(path) == item['sha256'], 'A sealed continuation source changed.')
    for filename, expected in bindings['external_upstream_code'].items():
        need(digest(Path(filename)) == expected, 'A frozen upstream code binding changed.')
    need(digest(Path(bindings['accepted_filter']['path']) / 'workflow-receipt.json') == bindings['accepted_filter']['receipt_sha256'],
         'Accepted stage5 receipt changed.')


def stop_child(child, log_path, cfg):
    if child and child.poll() is None:
        os.killpg(child.pid, signal.SIGINT)
        try:
            child.wait(timeout=25)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait(timeout=3)
    # Commands are printed by this exact child before execution. Never remove another job's container.
    if log_path.is_file():
        for line in log_path.read_text(errors='replace').splitlines():
            try:
                argv = json.loads(line).get('command', [])
            except (ValueError, AttributeError):
                continue
            if '--name' in argv and cfg['tools']['vep_image'] in argv:
                name = argv[argv.index('--name') + 1]
                if re.fullmatch('marina-wf-[a-f0-9]{12}', name):
                    for prefix in (['docker'], ['sudo', '-n', 'docker']):
                        result = subprocess.run(prefix + ['rm', '-f', name], capture_output=True, timeout=30)
                        if result.returncode == 0:
                            break


def execute_queue(args):
    global STOP
    os.umask(0o077)
    path = args.bindings.resolve(strict=True)
    binding_sha = digest(path)
    bindings = load(path)
    need(bindings.get('schema_version') == 1, 'Unsupported queue bindings.')
    state_path = Path(bindings['status_file'])
    if args.status:
        print(json.dumps(load(state_path) if state_path.is_file() else {'status': 'NOT_STARTED'}, indent=2))
        return
    verify_bindings(bindings)
    cfg = load(bindings['annotation_config']['path'])
    end = budget_deadline(**dict(created=bindings['budget']['vm_created_utc'], hourly=bindings['budget']['hourly_usd'], cap=bindings['budget']['planned_vm_usd']))
    monotonic_end = time.monotonic() + max(0, (end - dt.datetime.now(dt.timezone.utc)).total_seconds())
    code_root = Path(bindings['code_package']['path']).parent
    logs = Path(bindings['log_dir'])
    logs.mkdir(parents=True, exist_ok=True)
    child = None
    child_log = logs / 'queue-child.log'
    def remaining():
        return min((end - dt.datetime.now(dt.timezone.utc)).total_seconds(), monotonic_end - time.monotonic())
    def status(value, **extra):
        result = dict(schema_version=1, status=value, binding_sha256=binding_sha,
            updated_utc=dt.datetime.now(dt.timezone.utc).isoformat(), deadline_utc=end.isoformat(),
            output_root=bindings['output_root'], child_pid=child.pid if child and child.poll() is None else None,
            clinical_diagnosis_automatic=False, **extra)
        write_status(state_path, result)
        print(json.dumps(result), flush=True)
    def run(argv, logname):
        nonlocal child, child_log
        need(remaining() > 0, 'Existing VM allowance exhausted.')
        verify_bindings(bindings)
        child_log = logs / logname
        with child_log.open('a') as stream:
            child = subprocess.Popen(list(map(str, argv)), stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        while child.poll() is None:
            if STOP or remaining() <= 0:
                stop_child(child, child_log, cfg)
                raise TimeoutError('Queue stopped at cancellation or unchanged VM allowance deadline.')
            time.sleep(min(5, max(0, remaining())))
        need(child.returncode == 0, 'A continuation step failed; see its private log. No following patient stage was started.')
        child = None
    with (Path(bindings['data_root']) / 'private/.annotation-continuation.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            if not args.execute:
                status('PLAN_ONLY', patient_stages_started=False)
                return
            while True:
                if STOP or remaining() <= 0:
                    status('STOPPED_VM_ALLOWANCE_OR_CANCEL', billing_continues_until_VM_termination=True)
                    return
                approval = Path(bindings['approval_file'])
                if not approval.is_file() or load(approval).get('binding_sha256') != binding_sha:
                    status('WAITING_ROOT_REVIEW_OF_EXACT_SEALED_BINDING')
                elif not Path(bindings['resources']['path']).is_file():
                    status('WAITING_VEP_RESOURCES', AlphaMissense='ready_sealed')
                else:
                    resource = load(bindings['resources']['path'])
                    need(resource.get('status') == 'PASS' and resource.get('config_sha256') == bindings['resources']['expected_recipe_sha256']
                         and resource['vep']['image'] == bindings['resources']['vep_image'], 'VEP resource manifest differs from the approved recipe.')
                    break
                time.sleep(min(60, max(0, remaining())))
            verify_bindings(bindings)
            need(verify_stage(Path(bindings['accepted_filter']['path']), bindings['accepted_filter']['receipt_sha256']), 'Accepted stage5 no longer verifies.')
            need(digest(Path(bindings['reference']['fasta'])) == bindings['reference']['sha256'], 'Analysis reference changed.')
            status('PUBLIC_20_VARIANT_CONTROLS_RUNNING')
            run(['python3', code_root / 'scripts/validate_annotation_controls.py', '--config', bindings['annotation_config']['path'],
                 '--output-root', Path(bindings['data_root']) / 'resources/annotation-validation', '--package-sha', bindings['code_package']['sha256'], '--execute'], 'public-controls.log')
            status('PUBLIC_CONTROLS_PASS_PATIENT_ANNOTATION_RUNNING')
            run(['python3', code_root / 'scripts/workflow.py', 'annotate', '--root', bindings['data_root'], '--config', bindings['annotation_config']['path'],
                 '--execute'], '06-annotate.log')
            need(verify_stage(Path(bindings['output_root']) / 'stages/06_annotate'), 'Annotation did not produce a valid PASS receipt.')
            status('ANNOTATION_PASS_REPORT_RUNNING')
            run(['python3', code_root / 'scripts/workflow.py', 'report', '--root', bindings['data_root'], '--config', bindings['annotation_config']['path'],
                 '--execute'], '07-report.log')
            need(verify_stage(Path(bindings['output_root']) / 'stages/07_report'), 'Report did not produce a valid PASS receipt.')
            status('RESEARCH_REPORT_COMPLETE', report=str(Path(bindings['output_root']) / 'stages/07_report/report.html'),
                   docx=str(Path(bindings['output_root']) / 'stages/07_report/research_report.docx'), manual_ACMG_review_required=True)
        except BaseException as error:
            stop_child(child, child_log, cfg)
            status('FAILED_OR_BUDGET_STOPPED', error_type=type(error).__name__, inspect_private_logs=True)
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    setup = sub.add_parser('prepare')
    setup.add_argument('--base-config', type=Path, required=True)
    setup.add_argument('--runtime-root', type=Path, required=True)
    setup.add_argument('--upstream-code-root', type=Path, required=True)
    setup.add_argument('--created-utc', required=True)
    wait = sub.add_parser('run')
    wait.add_argument('--bindings', type=Path, default=Path('/srv/marina/private/annotation-queue-bindings.json'))
    wait.add_argument('--execute', action='store_true')
    wait.add_argument('--status', action='store_true')
    args = parser.parse_args()
    def stop(signum, frame):
        global STOP
        STOP = True
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    prepare(args) if args.action == 'prepare' else execute_queue(args)
