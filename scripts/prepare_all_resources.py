#!/usr/bin/env python3
"""One server entry point for automatic public resources and explicit remaining inventory."""
import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

CODE = Path(__file__).resolve().parents[1]


def need(ok, message):
    if not ok:
        raise ValueError(message)


def load(path, default=None):
    return json.loads(Path(path).read_text()) if Path(path).is_file() else (default or {})


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 ** 2), b''):
            value.update(block)
    return value.hexdigest()


def check_package():
    manifest = CODE / 'bootstrap-code-manifest.json'
    if manifest.is_file():
        for item in load(manifest)['files']:
            relative = Path(item['path'])
            path = CODE / relative
            need(not relative.is_absolute() and '..' not in relative.parts and not path.is_symlink()
                 and path.resolve(strict=True).is_relative_to(CODE) and path.stat().st_size == item['bytes']
                 and digest(path) == item['sha256'], 'Reviewed bootstrap source changed.')
    return digest(manifest) if manifest.is_file() else None


def lock_busy(path):
    if not path.exists():
        return False
    need(not path.is_symlink(), 'Refusing resource lock symlink.')
    with path.open('r+') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(stream, fcntl.LOCK_UN)
    return False


def source_items(recipe):
    cfg = load(recipe)
    prior = source_items(recipe.parent / cfg['base_config']) if cfg.get('base_config') else []
    return prior + [item for item in cfg['resources'] if item.get('enabled', True)]


def small_state(root, item, verify=False):
    key = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()
    folder = root / 'resources/annotation/.public-source-cache' / key
    receipt = load(folder / 'resource.json')
    path = folder / item['destination']
    if not receipt:
        return {'status': 'pending', 'version': item['version']}
    need(receipt.get('status') == 'available' and receipt.get('name') == item['name'] and receipt.get('version') == item['version']
         and path.is_file() and not path.is_symlink() and path.stat().st_size == receipt['bytes'],
         'Existing small resource receipt/file differs: ' + item['name'])
    if verify:
        need(digest(path) == receipt['sha256'], 'Existing small resource SHA mismatch: ' + item['name'])
    return {'status': 'ready_receipt', 'version': receipt['version'], 'sha256': receipt['sha256'], 'bytes': receipt['bytes']}


def inventory(root, recipe, verify=False):
    cfg = load(recipe)
    policy = load(CODE / 'config/resources-inventory.example.json')
    final = load(root / 'resources/annotation' / cfg['bundle_id'] / 'resource-manifest.json')
    rows = []
    for item in source_items(recipe):
        if item['name'] != 'vep_cache':
            rows.append(dict(resource=item['name'], mode='auto', **small_state(root, item, verify)))
    alpha_lock = CODE / 'config/alphamissense-resources.lock.json'
    alpha_recipe = load(alpha_lock)
    alpha_path = root / 'resources/annotation-supplements' / alpha_recipe['bundle_id'] / 'resource-manifest.json'
    alpha = load(alpha_path)
    if alpha:
        need(alpha.get('status') == 'PASS' and alpha.get('lock_sha256') == digest(alpha_lock), 'Existing AlphaMissense recipe differs.')
        for artifact in alpha['artifacts']:
            relative = Path(artifact['path'])
            path = alpha_path.parent / relative
            need(not relative.is_absolute() and '..' not in relative.parts and not path.is_symlink()
                 and path.resolve(strict=True).is_relative_to(alpha_path.parent.resolve())
                 and path.stat().st_size == artifact['bytes'], 'Existing AlphaMissense file differs.')
            if verify:
                need(digest(path) == artifact['sha256'], 'Existing AlphaMissense artifact SHA mismatch.')
    rows.append(dict(resource='alphamissense', mode='auto', status='ready_receipt' if alpha else 'pending', version=alpha.get('dataset_version')))
    for name, item in policy['cache_resources'].items():
        component = final.get('resources', {}).get(item['manifest_key'], {}) if item.get('manifest_key') else {}
        metadata = [line for line in final.get('cache_info', '').splitlines() if item.get('metadata_label', '\0').lower() in line.lower()]
        rows.append(dict(resource=name, mode='cache', status=component.get('status', 'cache_metadata_present' if metadata else 'waiting_verified_cache_info'),
            version=component.get('version'), cache_info_lines=metadata, note=item['note']))
    rows.append(dict(resource='vep_transcript_consequences_SO', mode='cache',
        status='ready_for_calculation' if final.get('status') == 'PASS' else 'waiting_verified_cache_and_public_controls',
        note=policy['computed_consequences']))
    for mode in ['manual', 'access_required', 'optional']:
        rows += [dict(resource=name, mode=mode, status='not_installed_by_this_command', note=note) for name, note in policy[mode].items()]
    known = {row['resource'].casefold() for row in rows}
    known |= {'hpo', 'clingen_gene_validity', 'ncbi_refseq_gene_mane', 'gnomad'}
    for item in load(CODE / 'config/evidence-sources.example.json')['resources']:
        if item['id'].casefold() not in known:
            rows.append(dict(resource=item['id'], mode='manual', status=policy['other_registry_default'], scope=item.get('scope'),
                note='Registry requirement acknowledged; not a claim of installation or full source content.'))
    return rows


def budget(root, created=None):
    state = load(root / 'private/databases-supervisor-status.json')
    settings = load(root / 'private/server-budget.json')
    original = state.get('identity', {}).get('created_utc') or settings.get('vm_created_utc')
    prior = state.get('identity') or settings
    need(float(prior.get('hourly_usd', 2.39)) == 2.39 and float(prior.get('planned_vm_usd', 40)) == 40,
         'This reviewed recipe is for $2.39/hour and the existing $40 allowance; review rates for a different VM.')
    text = created or original
    need(text, 'Actual VM creation UTC is required: --created-utc or existing private budget/supervisor receipt.')
    start = dt.datetime.fromisoformat(text.replace('Z', '+00:00'))
    need(start.tzinfo is not None, 'VM creation timestamp requires a timezone.')
    if original:
        need(start == dt.datetime.fromisoformat(original.replace('Z', '+00:00')), 'Refusing to reset the existing VM creation anchor.')
    now = dt.datetime.now(dt.timezone.utc)
    need(start <= now, 'VM creation timestamp cannot be in the future.')
    deadline = start + dt.timedelta(hours=40 / 2.39)
    need(deadline > now, 'Original $40 VM allowance is exhausted; no new download started. Billing must be stopped separately.')
    return start.isoformat(), deadline


def run_bounded(argv, deadline):
    remaining = (deadline - dt.datetime.now(dt.timezone.utc)).total_seconds()
    need(remaining > 0, 'Original VM allowance deadline reached.')
    print('RUN ' + json.dumps(list(map(str, argv))), flush=True)
    child = subprocess.Popen(list(map(str, argv)), stdin=subprocess.DEVNULL, start_new_session=True)
    try:
        need(child.wait(timeout=remaining) == 0, 'Resource preparation failed; partials preserved, inspect log.')
    except BaseException:
        if child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        raise


def show(result):
    print('СОСТОЯНИЕ: ' + result['status'], flush=True)
    for row in result['resources']:
        if row['mode'] in ('auto', 'cache'):
            print(f"{row['mode']:5} {row['resource']:30} {row['status']} {row.get('version') or ''}", flush=True)
    for mode in ('manual', 'access_required', 'optional'):
        names = ', '.join(row['resource'] for row in result['resources'] if row['mode'] == mode)
        print(mode + ' — НЕ устанавливается автоматически: ' + names, flush=True)
    print('Логи: ' + ', '.join(result['logs']), flush=True)
    print('Готовность ресурса не означает завершённую аннотацию или клиническую классификацию.', flush=True)


def prepare(root, execute=False, created=None, status=False):
    root = root.resolve()
    need(root.is_absolute() and len(root.parts) >= 3 and root.is_dir(), 'Use an existing dedicated absolute Linux data root.')
    package_sha = check_package()
    recipe = CODE / 'config/annotation-resources.example.json'
    state = load(root / 'private/databases-supervisor-status.json')
    rows = inventory(root, recipe, verify=False)
    supervisor_active = lock_busy(root / 'private/.databases-supervisor.lock')
    result = dict(status='EXISTING_SUPERVISOR_ACTIVE' if supervisor_active else 'RESOURCE_PREPARATION_PLAN',
        execute_requested=execute, no_patient_queries=True, no_analysis_started=True, source_package_sha256=package_sha,
        supervisor=state.get('status'), installer_pid=state.get('current_pid'),
        resources=rows, logs=[str(root / 'logs/annotation-databases.log'), str(root / 'logs/annotation-supervisor.log')],
        note='auto means permitted locked installer; manual/access_required/optional entries are not installed')
    show(result)
    if not execute or status:
        return result
    created, deadline = budget(root, created)
    (root / 'private').mkdir(exist_ok=True)
    bootstrap_lock = root / 'private/.all-resources-bootstrap.lock'
    need(not bootstrap_lock.is_symlink(), 'Refusing bootstrap lock symlink.')
    with bootstrap_lock.open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            result['status'] = 'ANOTHER_UNIFIED_BOOTSTRAP_ACTIVE'
            return result
        # Recheck after acquiring our lock; never duplicate an existing supervisor.
        supervisor_active = lock_busy(root / 'private/.databases-supervisor.lock')
        pending = any(row['mode'] == 'auto' and row['status'] == 'pending' for row in rows)
        if not pending and supervisor_active:
            result['status'] = 'ALREADY_RUNNING_REUSING_READY_SMALLS_AND_ALPHA'
            print(result['status'] + '; logs=' + ', '.join(result['logs']), flush=True)
            return result
        inventory(root, recipe, verify=True)
        if pending:
            run_bounded([sys.executable, CODE / 'scripts/prepare_annotation_resources.py', '--config', recipe,
                '--root', root / 'resources/annotation', '--small-only', '--execute'], deadline)
            run_bounded([sys.executable, CODE / 'scripts/prepare_alphamissense.py', '--lock', CODE / 'config/alphamissense-resources.lock.json',
                '--root', root, '--execute'], deadline)
        import supervise_databases as supervisor
        target = root / 'resources/annotation' / load(recipe)['bundle_id'] / 'resource-manifest.json'
        if target.exists() and supervisor.verified_ready(target, recipe):
            result['status'] = 'AUTOMATIC_RESOURCES_READY_MANUAL_SOURCES_SEPARATE'
            print(result['status'], flush=True)
            return result
        if lock_busy(root / 'private/.databases-supervisor.lock'):
            result['status'] = 'ALREADY_RUNNING'
            return result
        state = load(root / 'private/databases-supervisor-status.json')
        old = state.get('identity', {})
        installer = Path(old.get('installer') or CODE / 'publicstages/01_prepare_databases.sh')
        if old:
            need(old.get('recipe_sha256') == digest(recipe) and old.get('installer_sha256') == digest(installer)
                 and digest(installer) == digest(CODE / 'publicstages/01_prepare_databases.sh')
                 and digest(installer.parents[1] / 'scripts/prepare_annotation_resources.py') == digest(CODE / 'scripts/prepare_annotation_resources.py'),
                 'Existing supervisor uses a different reviewed installer/recipe; refusing replacement.')
        args = [sys.executable, '-u', str(CODE / 'scripts/supervise_databases.py'), '--root', str(root),
            '--installer', str(installer), '--config', str(recipe), '--created-utc', created,
            '--planned-usd', '40', '--hourly-usd', '2.39', '--workers', '32', '--retry-until-deadline']
        pid = state.get('current_pid')
        if pid and supervisor.process_token(pid, installer):
            args += ['--adopt-pid', str(pid)]
        (root / 'logs').mkdir(exist_ok=True)
        with (root / 'logs/annotation-supervisor.log').open('a') as log:
            child = subprocess.Popen(args, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        time.sleep(1)
        need(child.poll() is None or child.returncode == 0, 'Supervisor failed to start; inspect annotation-supervisor.log.')
        result.update(status='BUDGETED_SUPERVISOR_STARTED', supervisor_pid=child.pid, deadline_utc=deadline.isoformat())
        print(result['status'] + f" PID={child.pid}; deadline={deadline.isoformat()}", flush=True)
        return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--created-utc', help='Actual original VM creation timestamp, required on a new installation.')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--execute', action='store_true')
    mode.add_argument('--dry-run', action='store_true')
    mode.add_argument('--status', action='store_true')
    args = parser.parse_args()
    try:
        prepare(args.root, args.execute, args.created_utc, args.status)
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print('RESOURCE BOOTSTRAP FAILED: ' + str(error), file=sys.stderr)
        raise SystemExit(1)
