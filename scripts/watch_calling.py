#!/usr/bin/env python3
"""Read-only live status for stages 3-5; no sample identifiers are displayed."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import subprocess
import time


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def snapshot(root, label):
    private = root / 'private'
    cfg_path = private / ('workflow-' + label + '.json')
    result = {'utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'configuration_ready': cfg_path.is_file()}
    if not cfg_path.is_file():
        result['status'] = 'WAITING_FOR_CONFIGURATION'
        return result
    cfg = load(cfg_path)
    out = Path(cfg['output_root'])
    state_path = private / 'calling-3to5-state.json'
    state = load(state_path) if state_path.is_file() else {}
    result.update(status=state.get('status', 'NOT_STARTED'), stage=state.get('stage'),
                  sample_count=state.get('sample_count'), scope='explicit configuration callsets',
                  covered_callsets=[c['name'] for c in cfg.get('callsets', [])])
    pid = state.get('pid')
    alive = False
    if isinstance(pid, int):
        try:
            os.kill(pid, 0)
            alive = True
        except OSError:
            pass
    result['controller_alive'] = alive
    if result['status'] == 'RUNNING' and not alive:
        result['status'] = 'STOPPED_UNEXPECTEDLY'
    checkpoints = out / 'checkpoints' / '03_call' / 'gvcf'
    result['completed_gvcfs'] = sum(load(p).get('status') == 'PASS'
                                   for p in checkpoints.glob('*/*/workflow-receipt.json'))
    result['completed_stages'] = [stage for stage in ('03_call', '04_joint', '05_filter')
                                  if (out / 'stages' / stage / 'workflow-receipt.json').is_file()
                                  and load(out / 'stages' / stage / 'workflow-receipt.json').get('status') == 'PASS']
    logs = (list(checkpoints.glob('*/*.partial-*/haplotypecaller.log'))
            + list(checkpoints.glob('*/*.partial-*/parabricks.log')))
    if logs:
        latest = max(logs, key=lambda p: p.stat().st_mtime)
        result['HC_log_age_seconds'] = round(max(0, time.time() - latest.stat().st_mtime), 1)
        # Only parse progress counters. Never echo arbitrary patient paths/log lines.
        tail = latest.read_text(errors='replace').splitlines()[-100:]
        result['HC_progress'] = 'Running; waiting for the next periodic progress counter.'
        for line in reversed(tail):
            match = re.search(r'\b(chr[^:\s]+):[0-9,]+\s+([0-9.]+)\s+([0-9,]+)\s+([0-9.]+)', line)
            if match:
                result['HC_progress'] = {'contig': match[1], 'elapsed_minutes': float(match[2]),
                                         'regions_processed': int(match[3].replace(',', '')),
                                         'regions_per_minute': float(match[4])}
                break
            if 'Total time taken:' in line:
                result['HC_progress'] = 'GPU caller finished; output/index checks are in progress.'
                break
    if result.get('stage') == '4':
        result['current_task'] = 'Joint cohort genotyping (CPU); GPU may be idle.'
        scatter_root = out / 'checkpoints' / '04_joint' / 'scatter'
        result['joint_chromosomes_completed'] = sum(load(p).get('status') == 'PASS'
            for p in scatter_root.glob('*/workflow-receipt.json'))
        result['joint_chromosome_jobs_active'] = sum(not (p / 'failure.json').is_file()
            for p in scatter_root.glob('*.partial-*') if p.is_dir())
        result['joint_workers_limit'] = state.get('joint_workers')
    elif result.get('stage') == '5':
        result['current_task'] = 'Normalization and soft filtering (CPU); GPU may be idle.'
    try:
        text = subprocess.check_output(['nvidia-smi', '--query-gpu=utilization.gpu,memory.used,memory.total',
                '--format=csv,noheader,nounits'], text=True, stderr=subprocess.DEVNULL, timeout=10).strip().splitlines()
        gpu = []
        for line in text:
            values = [float(v.strip()) for v in line.split(',')]
            gpu.append({'utilization_percent': values[0], 'used_MiB': values[1], 'total_MiB': values[2]})
        result['GPU_live'] = gpu
    except (OSError, ValueError, subprocess.SubprocessError):
        result['GPU_live'] = 'unavailable'
    result['annotation_resource_manifest_available'] = Path(cfg['annotation']['resource_manifest']).is_file() if 'annotation' in cfg else False
    result['annotation_started_by_this_controller'] = False
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('/srv/marina'))
    parser.add_argument('--label', default='v1')
    parser.add_argument('--follow', action='store_true', help='Refresh until PASS/FAIL; Ctrl-C stops the watcher only.')
    parser.add_argument('--interval', type=float, default=15)
    parser.add_argument('--json', action='store_true')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,60}', args.label) or not 5 <= args.interval <= 300:
        parser.error('Use a safe label and a 5..300 second interval.')
    try:
        while True:
            try:
                result = snapshot(args.root, args.label)
            except (OSError, KeyError, ValueError):
                result = {'status': 'STATUS_TEMPORARILY_UNAVAILABLE'}
            if args.json:
                print(json.dumps(result), flush=True)
            else:
                print('UTC ' + result.get('utc', '') + ' | ' + result['status']
                      + ' | stage ' + str(result.get('stage', '-'))
                      + ' | gVCF ' + str(result.get('completed_gvcfs', 0)) + '/' + str(result.get('sample_count') or '?'), flush=True)
                if 'GPU_live' in result:
                    print('GPU live: ' + json.dumps(result['GPU_live']), flush=True)
                if 'HC_progress' in result:
                    print('HC: ' + json.dumps(result['HC_progress']), flush=True)
                if 'current_task' in result:
                    print(result['current_task'], flush=True)
                if 'joint_chromosomes_completed' in result:
                    print('Joint chromosomes completed: ' + str(result['joint_chromosomes_completed'])
                          + '; active jobs: ' + str(result['joint_chromosome_jobs_active'])
                          + '; worker limit: ' + str(result.get('joint_workers_limit') or '?'), flush=True)
            if not args.follow or result['status'] in ('PASS', 'FAIL', 'STOPPED_UNEXPECTEDLY'):
                return 0 if result['status'] not in ('FAIL', 'STOPPED_UNEXPECTEDLY') else 1
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print('Watcher stopped. The analysis keeps running.', flush=True)
        return 0


if __name__ == '__main__':
    raise SystemExit(main())
