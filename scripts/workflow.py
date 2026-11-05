#!/usr/bin/env python3
"""Seven numbered, explicit stages. Default is a plan; --execute performs a stage."""
import argparse
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from workflow_common import Runtime, STAGES, json_write, require, sha256

CODE = Path(__file__).resolve().parent.parent
ALIASES = {'1': 'prepare', '2': 'align', '3': 'call', '4': 'joint', '5': 'filter', '6': 'annotate', '7': 'report'}


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def default_config(args):
    return args.config or args.root / 'private' / ('workflow-' + args.label + '.json')


def base_config(root):
    cfg = {'data_root': str(root), 'output_root': str(root / 'outputs' / 'preparation'),
           'resources': {'max_hours': 6, 'min_free_gib': 100}, 'budget': {'hourly_usd': 2.39,
           'planned_vm_usd': 40, 'user_budget_usd': 50}}
    budget = root / 'private/server-budget.json'
    if budget.is_file():
        cfg['budget'].update(load(budget))
    return cfg


def prepare(args):
    cfg = base_config(args.root)
    runtime = Runtime(cfg, args.execute, args.root / 'logs/workflow-preparation', 'prepare')
    resource_config = args.resources_config or CODE / 'config/annotation-resources.example.json'
    resource_spec = load(resource_config)
    resource_root = args.root / 'resources/annotation'
    manifest = args.resource_manifest or resource_root / resource_spec['bundle_id'] / 'resource-manifest.json'
    print('STAGE 1: install/check tools and offline databases; no patient variant calls.', flush=True)
    command = [sys.executable, str(CODE / 'scripts/prepare_annotation_resources.py'),
               '--config', str(resource_config), '--root', str(resource_root)]
    if args.execute:
        command += ['--execute', '--pull-image']
        runtime.check_space()
        # Respect an already running preparation process instead of duplicating large downloads.
        lock_path = resource_root / '.annotation-prep.lock'
        if lock_path.exists() and not manifest.exists():
            with lock_path.open('a') as lock:
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    fcntl.flock(lock, fcntl.LOCK_UN)
                except BlockingIOError:
                    require(args.wait, 'Resource installation is already running. Retry stage1 with --wait.')
                    while not manifest.exists():
                        runtime.remaining_seconds()
                        print('WAITING: offline database installation is still running.', flush=True)
                        time.sleep(min(30, runtime.remaining_seconds()))
        runtime.run(command, 'annotation-resources')
        from configure_workflow import GATK_IMAGE
        runtime.run(runtime.docker_prefix() + ['pull', GATK_IMAGE], 'gatk-image-install', timeout=1800)
        runtime.run(runtime.docker(GATK_IMAGE, ['gatk', '--version'], cpus=2, memory_gib=8), 'gatk-version')
        command = [sys.executable, str(CODE / 'scripts/configure_workflow.py'), '--root', str(args.root),
                   '--resource-manifest', str(manifest), '--label', args.label, '--execute']
        if args.callsets_manifest:
            command += ['--callsets-manifest', str(args.callsets_manifest)]
        runtime.run(command, 'private-workflow-configuration')
    else:
        print('PLAN:', ' '.join(command))
        print('The prepared configuration will use verified BAMs, fixed software and explicit intervals.')
    return manifest


def align(args):
    root = args.root
    cfg = base_config(root)
    runtime = Runtime(cfg, args.execute, root / 'logs/workflow-alignment', 'align')
    print('STAGE 2: preserve the first verified whole-reference BAM and align the remaining samples.', flush=True)
    configure = [sys.executable, str(CODE / 'scripts/configure_remaining_align.py'), '--root', str(root)]
    runner = [sys.executable, str(CODE / 'scripts/run_parabricks.py'), '--config',
              str(root / 'private/remaining-align.json'), '--stage', 'align']
    if args.execute:
        configure.append('--execute'); runner.append('--execute')
        queue_lock = root / 'outputs/remaining-align/.gpu-queue.lock'
        if queue_lock.exists():
            with queue_lock.open('a') as lock:
                while True:
                    try:
                        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        fcntl.flock(lock, fcntl.LOCK_UN)
                        break
                    except BlockingIOError:
                        require(args.wait, 'The remaining alignment queue is already running. Use --wait to wait for it.')
                        print('WAITING: the existing alignment queue is running; it will not be duplicated.', flush=True)
                        time.sleep(min(30, runtime.remaining_seconds()))
        runtime.run(configure, 'remaining-config')
        runtime.run(runner, 'remaining-alignment')
    else:
        print('PLAN:', ' '.join(configure)); print('PLAN:', ' '.join(runner + ['--dry-run']))


def dependencies(stage, cfg, runtime):
    paths = [Path(cfg['samples_manifest']), Path(cfg['alignment_manifest'])]
    if stage == 'call':
        import csv
        with Path(cfg['alignment_manifest']).open(newline='') as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                paths += [Path(row['receipt']), Path(row['bam'])]
        paths += [Path(cfg['reference']['fasta'])]
        paths += [Path(item['intervals']) for item in cfg['callsets']]
    else:
        previous = {'joint': 'call', 'filter': 'joint', 'annotate': 'filter', 'report': 'annotate'}[stage]
        previous_dir = runtime.stage_path(previous)
        receipt = previous_dir / 'workflow-receipt.json'
        require(receipt.is_file(), 'The preceding numbered stage is not completed: ' + previous)
        saved = load(receipt)
        require(saved.get('status') == 'PASS', 'Previous stage does not have a successful receipt.')
        # Verify the actual predecessor artifacts, not only its status JSON.
        for name, expected in saved['artifacts'].items():
            path = (previous_dir / name).resolve(strict=True)
            require(path.is_relative_to(previous_dir.resolve()) and path.stat().st_size == expected['bytes']
                    and runtime.file_hash(path) == expected['sha256'], 'A preceding stage artifact changed.')
        paths.append(receipt)
    if stage in ('annotate', 'report'):
        paths.append(Path(cfg['annotation']['resource_manifest']))
    return paths


def downstream(args, stage, cfg, shared=None):
    runtime = shared.child(Path(cfg['output_root']), stage) if shared else Runtime(cfg, args.execute, stage=stage)
    print('STAGE ' + STAGES[stage] + ': ' + ('EXECUTE' if args.execute else 'PLAN'), flush=True)
    print('Analysis scope:', json.dumps(cfg.get('scope', {}), ensure_ascii=False), flush=True)
    if not args.execute:
        print('No patient calls, filtering or annotation started. Use --execute for this numbered stage.')
        return
    require(cfg.get('schema_version') == 2, 'Expected a prepared schema_version2 workflow configuration.')
    paths = dependencies(stage, cfg, runtime)
    if stage in ('call', 'joint', 'filter'):
        from calling_stages import run_stage
    else:
        from annotation_stages import run_stage
    def produce(directory):
        child = runtime.child(directory, stage)
        return run_stage(stage, cfg, directory, child)
    result = runtime.run_stage(stage, paths, {'stage': stage}, produce)
    json_write(Path(cfg['data_root']) / 'private/workflow-status.json', {
        'status': 'PASS', 'stage': STAGES[stage], 'output': str(result),
        'updated_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'scope': cfg.get('scope', {}), 'clinical_classification_complete': False})


def status(args):
    config = default_config(args)
    result = {'configuration_ready': config.is_file(), 'stages': {},
              'note': 'Completion is technical; inspect scope, coverage and manual evidence requirements.'}
    if config.is_file():
        cfg = load(config)
        result['scope'] = cfg.get('scope')
        for stage in ('call', 'joint', 'filter', 'annotate', 'report'):
            path = Path(cfg['output_root']) / 'stages' / STAGES[stage] / 'workflow-receipt.json'
            result['stages'][STAGES[stage]] = load(path).get('status') if path.is_file() else 'NOT_COMPLETED'
    result['alignment_process_note'] = 'The original alignment command has its own receipts under outputs.'
    print(json.dumps(result, indent=2, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=tuple(ALIASES) + tuple(STAGES) + ('all', 'status'))
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--config', type=Path)
    parser.add_argument('--label', default='v1')
    parser.add_argument('--resources-config', type=Path)
    parser.add_argument('--resource-manifest', type=Path)
    parser.add_argument('--callsets-manifest', type=Path)
    parser.add_argument('--wait', action='store_true')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--execute', action='store_true')
    mode.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    args.root = args.root.resolve()
    require(args.root.is_absolute() and len(args.root.parts) >= 3, 'Use a dedicated absolute project root.')
    os.umask(0o077)
    stage = ALIASES.get(args.stage, args.stage)
    if stage == 'status':
        status(args); return
    if stage in ('prepare', 'align'):
        getattr(sys.modules[__name__], stage)(args); return
    config_path = default_config(args)
    require(config_path.is_file(), 'Prepare the private workflow plan first with stage1; annotation resources must be ready.')
    cfg = load(config_path)
    require(Path(cfg['data_root']).resolve() == args.root, 'Configuration belongs to a different data root.')
    if not args.execute:
        for item in ('call', 'joint', 'filter', 'annotate', 'report') if stage == 'all' else (stage,):
            downstream(args, item, cfg)
        return
    private = args.root / 'private'
    private.mkdir(parents=True, exist_ok=True)
    with (private / '.downstream-workflow.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Another downstream workflow is already running; use status.') from None
        shared = Runtime(cfg, True)
        for item in ('call', 'joint', 'filter', 'annotate', 'report') if stage == 'all' else (stage,):
            shared.remaining_seconds()
            downstream(args, item, cfg, shared)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, subprocess.SubprocessError) as error:
        print('WORKFLOW FAILED:', str(error), file=sys.stderr)
        raise SystemExit(1)
