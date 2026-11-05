#!/usr/bin/env python3
"""Bounded, explicitly reviewed, local-only postannotation sidecar continuation."""
import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

ABSOLUTE_DEADLINE = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def validate(rule, rule_path, reviewed_rule_sha256, current_time=None):
    if sha(rule_path) != reviewed_rule_sha256:
        raise ValueError('ROOT-reviewed rule SHA differs')
    if rule.get('review_status') != 'ROOT_REVIEWED_EXACT_LOCAL_ONLY_SIDECAR_RULE':
        raise ValueError('ROOT-reviewed rule required before execution')
    expected_keys = {'schema_version', 'review_status', 'reviewed_utc', 'code_sha256', 'context_index',
                     'context_index_sha256', 'context_manifest_sha256', 'workflow_root', 'output', 'state_directory',
                     'deadline_utc', 'interval_seconds', 'max_hours', 'no_network', 'no_primary_workflow_mutation'}
    if set(rule) != expected_keys or rule['schema_version'] != 1:
        raise ValueError('Unknown or missing rule fields')
    if rule['no_network'] is not True or rule['no_primary_workflow_mutation'] is not True:
        raise ValueError('Local-only sidecar contract required')
    if not 10 <= rule['interval_seconds'] <= 60 or not 0 < rule['max_hours'] <= 8:
        raise ValueError('Bounded wait required')
    deadline = dt.datetime.fromisoformat(rule['deadline_utc'])
    current_time = current_time or dt.datetime.now(dt.timezone.utc)
    if deadline > ABSOLUTE_DEADLINE or deadline <= current_time:
        raise ValueError('Original absolute budget deadline exceeded or already reached')
    scripts = Path(__file__).resolve().parent
    expected_scripts = {str(Path(__file__).resolve()), str(scripts / 'export_tissue_protein_context.py')}
    if set(rule['code_sha256']) != expected_scripts:
        raise ValueError('Only exact reviewed controller/exporter pair accepted')
    for path, expected in rule['code_sha256'].items():
        if Path(path).is_symlink() or sha(path) != expected:
            raise ValueError('Reviewed code digest differs: ' + path)
    context = Path(rule['context_index'])
    if context.is_symlink() or sha(context) != rule['context_index_sha256'] or sha(context.parent / 'manifest.json') != rule['context_manifest_sha256']:
        raise ValueError('Public context digest differs')
    workflow = Path(rule['workflow_root']).resolve()
    output = Path(rule['output']).resolve()
    state = Path(rule['state_directory']).resolve()
    if output == workflow or workflow in output.parents or state == workflow or workflow in state.parents:
        raise ValueError('State/output must be outside primary workflow')
    if state == output or state in output.parents or output in state.parents:
        raise ValueError('State and immutable output directories must be separate')
    if output.exists():
        raise ValueError('Sidecar output already exists; refusing overwrite')
    return min(deadline, current_time + dt.timedelta(hours=rule['max_hours']))


def readiness(workflow):
    try:
        receipts = [load(Path(workflow) / 'stages' / stage / 'workflow-receipt.json') for stage in ('06_annotate', '07_report')]
    except FileNotFoundError:
        return False
    if any(str(r.get('status', '')).startswith('FAIL') for r in receipts):
        raise ValueError('Upstream stage failed; no automatic restart')
    return all(r.get('status') == 'PASS' for r in receipts)


def run(rule_path, reviewed_rule_sha256, execute=False):
    rule = load(rule_path)
    state = Path(rule['state_directory'])
    if not execute:
        return load(state / 'status.json') if (state / 'status.json').exists() else dict(status='NOT_STARTED_ROOT_REVIEW_REQUIRED')
    deadline = validate(rule, rule_path, reviewed_rule_sha256)
    os.umask(0o077); state.mkdir(parents=True, mode=0o700, exist_ok=True)
    import fcntl
    with (state / 'continuation.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('This sidecar continuation already owns its lock')
        def status(value, **extra):
            doc = dict(status=value, pid=os.getpid(), rule_sha256=reviewed_rule_sha256,
                       updated_utc=dt.datetime.now(dt.timezone.utc).isoformat(), deadline_utc=deadline.isoformat(), **extra)
            path = state / 'status.json'; temp = state / 'status.tmp'
            if path.is_symlink() or temp.is_symlink(): raise ValueError('Symlink state refused')
            temp.write_text(json.dumps(doc, indent=2) + '\n'); temp.replace(path)
            return doc
        try:
            while not readiness(rule['workflow_root']):
                remaining = (deadline - dt.datetime.now(dt.timezone.utc)).total_seconds()
                if remaining <= 0:
                    return status('STOPPED_DEADLINE_NO_EXPORT_NO_PRIMARY_ACTUATION')
                status('WAITING_STAGE6_AND_STAGE7_PASS_NO_PRIMARY_ACTUATION')
                time.sleep(min(rule['interval_seconds'], remaining))
            # Validate pins again at point of use after potentially hours of waiting.
            validate(rule, rule_path, reviewed_rule_sha256)
            status('EXPORTING_SEPARATE_LOCAL_RESEARCH_CONTEXT')
            helper = Path(__file__).resolve().parent / 'export_tissue_protein_context.py'
            spec = importlib.util.spec_from_file_location('reviewed_context_exporter', helper)
            module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
            result = module.enrich(rule['workflow_root'], rule['context_index'], rule['output'], deadline.isoformat())
            return status('PASS_SEPARATE_CONTEXT_EXPORT_PRIMARY_UNCHANGED', result=result)
        except BaseException as error:
            status('FAILED_NO_AUTOMATIC_RETRY_PRIMARY_UNCHANGED', error=str(error)[:400])
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rule', type=Path, required=True)
    parser.add_argument('--reviewed-rule-sha256')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    print(json.dumps(run(args.rule, args.reviewed_rule_sha256, args.execute), indent=2))
