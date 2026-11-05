#!/usr/bin/env python3
"""Private, resumable stage execution shared by the downstream workflow."""
import datetime as dt
import ctypes
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import signal
import subprocess
import time
import uuid

STAGES = {'prepare': '01_prepare', 'align': '02_align', 'call': '03_call',
          'joint': '04_joint', 'filter': '05_filter', 'annotate': '06_annotate', 'report': '07_report'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), default=str).encode()).hexdigest()


def safe_key(value):
    path = PurePosixPath(str(value))
    require(not path.is_absolute() and path.parts and all(
        re.fullmatch(r'[A-Za-z0-9_.-]{1,160}', part) and part not in ('.', '..') for part in path.parts),
        'Invalid stage/checkpoint name.')
    return Path(*path.parts)


def json_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name('.' + path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temporary.open('x', encoding='utf-8') as stream:
            json.dump(value, stream, indent=2, ensure_ascii=False, default=str)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def artifact_manifest(directory):
    directory = Path(directory).resolve()
    result = {}
    for path in sorted(directory.rglob('*')):
        require(not path.is_symlink(), 'Stage outputs cannot contain symlinks.')
        if path.is_file() and path != directory / 'workflow-receipt.json':
            result[path.relative_to(directory).as_posix()] = {
                'bytes': path.stat().st_size, 'sha256': sha256(path)}
    require(result, 'Stage produced no output artifacts.')
    return result


def promote_exclusive(source, target):
    libc = ctypes.CDLL(None, use_errno=True)
    require(hasattr(libc, 'renameat2'), 'Atomic no-replace promotion requires Linux renameat2.')
    rename = libc.renameat2
    rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(source), -100, os.fsencode(target), 1):
        number = ctypes.get_errno()
        raise OSError(number, os.strerror(number), str(target))


class Runtime:
    def __init__(self, cfg, execute=False, workdir=None, stage=None, started=None):
        self.cfg = cfg
        self.execute = bool(execute)
        self.dry_run = not self.execute
        self.output_root = Path(cfg['output_root']).resolve()
        self.workdir = Path(workdir).resolve() if workdir else self.output_root
        self.stage = stage or 'prepare'
        self.started = started if started is not None else time.monotonic()
        self.commands = []
        self._hashes = {}
        self._docker = None
        self._containers = set()
        self._resources = None
        self.code_root = Path(__file__).resolve().parent.parent

    def child(self, directory, stage=None):
        child = Runtime(self.cfg, self.execute, directory, stage or self.stage, self.started)
        child._hashes = self._hashes
        child._containers = self._containers
        child._docker = self._docker
        return child

    def file_hash(self, path):
        path = Path(path).resolve(strict=True)
        st = path.stat()
        key = (str(path), st.st_size, st.st_mtime_ns, st.st_ino)
        if key not in self._hashes:
            self._hashes[key] = sha256(path)
        return self._hashes[key]

    def input_identity(self, paths):
        result = {}
        for source in paths:
            path = Path(source).resolve(strict=True)
            if path.is_dir():
                files = [p for p in sorted(path.rglob('*')) if p.is_file()]
                result[str(path)] = {p.relative_to(path).as_posix(): self.file_hash(p) for p in files}
            else:
                result[str(path)] = self.file_hash(path)
        return result

    def code_identity(self):
        names = ['workflow_common.py', 'workflow.py']
        names += ['calling_stages.py'] if self.stage in ('call', 'joint', 'filter') else ['annotation_stages.py']
        if self.stage == 'joint' and (self.code_root / 'scripts/joint_scatter.py').is_file():
            names.append('joint_scatter.py')
        if self.stage in ('annotate', 'report') and (self.code_root / 'scripts/alphamissense_adapter.py').is_file():
            names.append('alphamissense_adapter.py')
        if self.stage == 'report' and (self.code_root / 'scripts/export_report_docx.py').is_file():
            names.append('export_report_docx.py')
        return {name: self.file_hash(self.code_root / 'scripts' / name) for name in sorted(names)}

    def remaining_seconds(self):
        limits = self.cfg.get('resources', {})
        remaining = float(limits.get('max_hours', 6)) * 3600 - (time.monotonic() - self.started)
        budget = self.cfg.get('budget', {})
        if budget.get('vm_created_utc'):
            created = dt.datetime.fromisoformat(budget['vm_created_utc'].replace('Z', '+00:00'))
            require(created.tzinfo is not None, 'VM budget start needs a timezone.')
            hours = (dt.datetime.now(dt.timezone.utc) - created).total_seconds() / 3600
            rate = float(budget['hourly_usd'])
            cap = float(budget.get('planned_vm_usd', 40))
            require(rate > 0 and 0 < cap <= float(budget.get('user_budget_usd', 50)), 'Invalid cost allocation.')
            remaining = min(remaining, (cap / rate - hours) * 3600)
        require(remaining > 0, 'Workflow time/cost allocation exhausted. VM billing must be ended separately.')
        return remaining

    def check_space(self):
        parent = self.output_root
        while not parent.exists():
            parent = parent.parent
        reserve = int(self.cfg.get('resources', {}).get('min_free_gib', 100)) * 2**30
        require(shutil.disk_usage(parent).free >= reserve, 'Insufficient disk reserve for the next stage.')

    def docker_prefix(self):
        if self._docker is None:
            for prefix in (['docker'], ['sudo', '-n', 'docker']):
                if shutil.which(prefix[0]) and subprocess.run(prefix + ['info'], stdin=subprocess.DEVNULL,
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=30).returncode == 0:
                    self._docker = prefix
                    break
            require(self._docker is not None, 'Docker access is unavailable; run preparation first.')
        return self._docker

    def docker(self, image, args, mounts=None, gpu=False, cpus=4, memory_gib=16):
        require(re.fullmatch(r'[A-Za-z0-9./:_-]+@sha256:[0-9a-f]{64}', image or ''),
                'Containers must use an immutable image digest.')
        require(0 < float(cpus) <= 30 and 0 < float(memory_gib) <= 180, 'Invalid container resource limits.')
        prefix = self.docker_prefix() if self.execute else ['docker']
        name = 'marina-wf-' + uuid.uuid4().hex[:12]
        command = prefix + ['run', '--rm', '--pull=never', '--name', name, '--network=none',
            '--cpus', str(cpus), '--memory', str(memory_gib) + 'g', '--user', f'{os.getuid()}:{os.getgid()}']
        if gpu:
            command += ['--gpus', 'device=' + str(self.cfg.get('resources', {}).get('gpu_device', 0))]
        for mount in mounts or []:
            source = Path(mount['source']).resolve()
            target = PurePosixPath(str(mount['target']))
            require(target.is_absolute() and '..' not in target.parts and ',' not in str(source), 'Invalid container mount.')
            spec = f'type=bind,src={source},dst={target}'
            if mount.get('readonly', False):
                spec += ',readonly'
            command += ['--mount', spec]
        command += [image] + list(map(str, args))
        return command

    def run(self, argv, logname='command', capture=False, cwd=None, stdout=None, stderr=None, timeout=None):
        argv = list(map(str, argv))
        self.commands.append(argv)
        print(json.dumps({'stage': self.stage, 'command': argv}), flush=True)
        if not self.execute:
            return ''
        self.workdir.mkdir(parents=True, exist_ok=True)
        logpath = Path(logname)
        if not logpath.is_absolute():
            logpath = self.workdir / (str(logname) if str(logname).endswith('.log') else str(logname) + '.log')
        logpath.parent.mkdir(parents=True, exist_ok=True)
        allowed = self.remaining_seconds()
        if timeout is not None:
            allowed = min(allowed, float(timeout))
        container = argv[argv.index('--name') + 1] if '--name' in argv and '--pull=never' in argv else None
        if container:
            self._containers.add(container)
        handles = []
        proc = None
        try:
            log = logpath.open('ab'); handles.append(log)
            if isinstance(stdout, (str, Path)):
                output = Path(stdout).open('wb'); handles.append(output)
            else:
                output = stdout if stdout is not None else (subprocess.PIPE if capture else log)
            if isinstance(stderr, (str, Path)):
                errors = Path(stderr).open('ab'); handles.append(errors)
            else:
                errors = stderr if stderr is not None else log
            proc = subprocess.Popen(argv, cwd=cwd, stdin=subprocess.DEVNULL, stdout=output,
                                    stderr=errors, start_new_session=True)
            captured, _ = proc.communicate(timeout=allowed)
            if captured:
                log.write(captured)
            if proc.returncode:
                raise subprocess.CalledProcessError(proc.returncode, argv)
            return captured.decode('utf-8') if capture and captured else ''
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            if proc and proc.poll() is None:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
            if container:
                subprocess.run(self.docker_prefix() + ['stop', '--time', '30', container],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=45)
            raise
        finally:
            for handle in handles:
                handle.close()
            if container:
                self._containers.discard(container)

    def require_resource(self, name):
        manifest_path = Path(self.cfg['annotation']['resource_manifest']).resolve(strict=True)
        if self._resources is None:
            self._resources = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
        entry = self._resources['resources'].get(name)
        require(entry and str(entry.get('status', '')).lower() in ('ready', 'pass', 'available'),
                'Required offline resource is not available: ' + name)
        path = Path(entry['path'])
        if not path.is_absolute():
            path = manifest_path.parent / path
        path = path.resolve(strict=True)
        if path.is_file() and entry.get('sha256'):
            require(self.file_hash(path) == entry['sha256'], 'Resource checksum changed: ' + name)
        return path

    def stage_path(self, stage):
        require(stage in STAGES, 'Unknown workflow stage.')
        return self.output_root / 'stages' / STAGES[stage]

    def _valid(self, directory, identity):
        receipt_path = directory / 'workflow-receipt.json'
        require(receipt_path.is_file(), 'Existing output has no completed receipt; refusing overwrite.')
        receipt = json.loads(receipt_path.read_text())
        require(receipt.get('status') == 'PASS' and receipt.get('identity') == identity,
                'Completed output has different inputs, configuration or code; use a new output root.')
        for name, expected in receipt['artifacts'].items():
            path = (directory / name).resolve(strict=True)
            require(path.is_relative_to(directory.resolve()) and path.stat().st_size == expected['bytes']
                    and self.file_hash(path) == expected['sha256'], 'A completed output artifact changed.')
        return receipt

    def checkpoint(self, key, input_paths, parameters, producer):
        target = self.output_root / 'checkpoints' / STAGES[self.stage] / safe_key(key)
        return self._checkpoint_at(target, key, input_paths, parameters, producer)

    def run_stage(self, stage, input_paths, parameters, producer):
        require(self.stage == stage, 'Runtime stage must match the requested stage.')
        target = self.stage_path(stage)
        return self._checkpoint_at(target, STAGES[stage], input_paths, parameters, producer)

    def _checkpoint_at(self, target, label, input_paths, parameters, producer):
        identity = canonical_hash({'stage': label, 'inputs': self.input_identity(input_paths),
            'parameters': parameters, 'config': self.cfg, 'code': self.code_identity()})
        if target.exists():
            self._valid(target, identity)
            print('VERIFIED COMPLETE: ' + str(label), flush=True)
            return target
        if not self.execute:
            print('PLANNED STAGE: ' + str(label), flush=True)
            return target
        self.remaining_seconds()
        self.check_space()
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_name(target.name + '.partial-' + uuid.uuid4().hex[:12])
        partial.mkdir(mode=0o700)
        previous_workdir = self.workdir
        self.workdir = partial
        started = dt.datetime.now(dt.timezone.utc).isoformat()
        try:
            result = producer(partial) or {}
            artifacts = artifact_manifest(partial)
            receipt = {'schema_version': 1, 'status': 'PASS', 'stage': str(label), 'identity': identity,
                'started_utc': started, 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
                'scope': self.cfg.get('scope', {}), 'artifacts': artifacts, 'metadata': result,
                'note': 'Technical stage completion is not a clinical classification.'}
            json_write(partial / 'workflow-receipt.json', receipt)
            require(not target.exists(), 'Destination appeared during execution; refusing overwrite.')
            promote_exclusive(partial, target)
            print('COMPLETED: ' + str(label), flush=True)
            return target
        except BaseException as error:
            json_write(partial / 'failure.json', {'status': 'FAIL', 'stage': str(label),
                'error_type': type(error).__name__, 'started_utc': started,
                'failed_utc': dt.datetime.now(dt.timezone.utc).isoformat()})
            raise
        finally:
            self.workdir = previous_workdir
