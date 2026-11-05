#!/usr/bin/env python3
"""Persist the user's local-backup deadline and stop only this project's watchers.

Does not stop the remote server or control browsers. Native-download UI is paused
by its authorized UI owner. Existing files, partials, plans and indices are kept.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import signal
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--deadline-utc', required=True)
    p.add_argument('--backup-root', required=True)
    p.add_argument('--state-root', required=True)
    p.add_argument('--execute', action='store_true')
    args = p.parse_args()
    due = dt.datetime.fromisoformat(args.deadline_utc.replace('Z', '+00:00'))
    if due.tzinfo is None:
        raise ValueError('Absolute UTC deadline required')
    root, state = Path(args.backup_root).resolve(strict=True), Path(args.state_root).resolve(strict=True)
    if not args.execute:
        print('DRY RUN: no pause marker or process changes')
        return
    session = {'Status': 'ACTIVE_UNTIL_DEADLINE', 'DeadlineUtc': due.isoformat(), 'ServerWorkUnaffected': True, 'RequiresUserResumeAfterDeadline': True}
    (state / 'local-download-session.json').write_text(json.dumps(session, indent=2), encoding='utf-8')
    while time.time() < due.timestamp():
        time.sleep(min(1, max(0, due.timestamp() - time.time())))
    flag = root / 'Settings' / 'LOCAL_DOWNLOADS_PAUSED.json'
    flag.parent.mkdir(parents=True, exist_ok=True)
    payload = {'Status': 'PAUSED', 'DeadlineUtc': due.isoformat(),
               'PausedUtc': dt.datetime.now(dt.timezone.utc).isoformat(),
               'PartialsPreserved': True, 'ServerWorkUnaffected': True, 'RequiresUserResume': True,
               'Reason': 'Explicit user request: local downloads allowed only until deadline'}
    try:
        with flag.open('x', encoding='utf-8') as stream:
            json.dump(payload, stream, indent=2)
    except FileExistsError:
        pass
    scripts = Path(__file__).resolve().parent
    allowed = {str(scripts / 'watch_native_imports.py'), str(scripts / 'prepare_native_plan.py')}
    stopped = 0
    for proc in Path('/proc').iterdir():
        if not proc.name.isdigit():
            continue
        try:
            argv = (proc / 'cmdline').read_bytes().split(b'\0')
            paths = [a.decode(errors='replace') for a in argv]
            if not any(a in allowed for a in paths) or '--execute' not in paths:
                continue
            # Exact project script paths identify these owned local controllers.
            os.kill(int(proc.name), signal.SIGTERM)
            stopped += 1
        except (OSError, ValueError):
            pass
    payload['StoppedLocalControllers'] = stopped
    (state / 'local-download-session.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    (root / 'native-watch-status.json').write_text(json.dumps(payload, indent=2), encoding='utf-8')
    print(f'PAUSED at user deadline; stopped {stopped} owned controllers; partials preserved; server unaffected', flush=True)


if __name__ == '__main__':
    main()
