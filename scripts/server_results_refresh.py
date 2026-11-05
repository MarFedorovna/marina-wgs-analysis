#!/usr/bin/env python3
"""Maintain the receipt-backed server RESULTS view independently of the laptop."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--remote-root', required=True)
    parser.add_argument('--workspace-root', required=True)
    parser.add_argument('--jupyter-root', required=True)
    parser.add_argument('--interval', type=int, default=60)
    parser.add_argument('--max-hours', type=float, default=8)
    args = parser.parse_args()
    if args.interval < 60 or not 0 < args.max_hours <= 12:
        raise ValueError('Use interval >=60 and bounded max-hours <=12')
    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo / 'transfer'))
    from prepare_native_plan import HELPER
    os.nice(10)
    deadline = time.time() + args.max_hours * 3600
    cfg = {'root': args.remote_root, 'workspace': args.workspace_root, 'jupyter': args.jupyter_root}
    while time.time() < deadline:
        result = subprocess.run([sys.executable, '-c', HELPER, json.dumps(cfg)], text=True,
                                capture_output=True, timeout=min(120, max(1, deadline-time.time())))
        if result.returncode:
            print('Server view refresh error: ' + result.stderr, flush=True)
        else:
            plan = json.loads(result.stdout)
            print(json.dumps({'ServerView': plan['ServerView'], 'CompletedArtifacts': len(plan['Files'])}), flush=True)
        time.sleep(min(args.interval, max(0, deadline-time.time())))


if __name__ == '__main__':
    main()
