#!/usr/bin/env python3
"""Install an isolated, pinned Mutalyzer proof using only official public fixtures.

No patient input; never installs into an annotation runtime. Network is used only
for PyPI package files and the explicit upstream commit fixture allowlist.
"""
import argparse
import datetime as dt
import hashlib
import gzip
import csv
import json
from pathlib import Path
import subprocess
import sys
import urllib.request
import venv

VERSION = '3.1.1'
COMMIT = '171f36ada65b228c443b359925305bd85e9fa40e'
REFERENCE_IDS = ['NM_003002.4', 'NC_012920.1']
DEADLINE = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')


def get(url, limit=2_000_000):
    with urllib.request.urlopen(url, timeout=30) as response:
        value = response.read(limit + 1)
    if len(value) > limit:
        raise ValueError('Public metadata/fixture exceeds fixed byte ceiling')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--lock', type=Path, default=Path(__file__).resolve().parents[1] / 'config/offline-hgvs-mutalyzer3.1.1-py310.lock.txt')
    parser.add_argument('--mane-summary', type=Path, help='Existing public MANE1.4 summary; adds two fixed public HCM transcript controls')
    args = parser.parse_args()
    if dt.datetime.now(dt.timezone.utc) >= DEADLINE:
        raise SystemExit('Original absolute VM deadline reached')
    if not args.lock.is_file():
        raise SystemExit('Checked dependency hash lock is required')
    root = args.root.resolve()
    if root.exists() and any(root.iterdir()):
        raise SystemExit('Use a new empty isolated directory; refusing overwrite')
    root.mkdir(parents=True, exist_ok=True)
    source = root / 'public-reference-fixtures'
    source.mkdir()
    manifest = []
    reference_ids = list(REFERENCE_IDS)
    for reference in REFERENCE_IDS:
        for suffix in ['fasta', 'gff3']:
            name = reference + '.' + suffix
            url = f'https://raw.githubusercontent.com/mutalyzer/mutalyzer/{COMMIT}/tests/data/{name}'
            payload = get(url)
            (source / name).write_bytes(payload)
            manifest.append(dict(reference_id=reference, name=name, url=url,
                                 bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(),
                                 role='upstream public test fixture; not project MANE/GENCODE cache'))
    mane_pairs = []
    if args.mane_summary:
        expected = {'MYBPC3': 'NM_000256.3', 'MYH7': 'NM_000257.4'}
        with gzip.open(args.mane_summary, 'rt') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row.get('symbol') in expected and row.get('MANE_status') == 'MANE Select':
                    if row.get('RefSeq_nuc') != expected[row['symbol']]:
                        raise RuntimeError('MANE transcript version changed; review required')
                    mane_pairs.append({key: row[key] for key in ['symbol', 'RefSeq_nuc', 'RefSeq_prot', 'Ensembl_nuc', 'Ensembl_prot', 'GRCh38_chr', 'chr_strand']})
        if len(mane_pairs) != 2:
            raise RuntimeError('Expected two public HCM MANE Select rows')
        for reference in expected.values():
            for suffix in ['fasta', 'gff3']:
                url = ('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&id='
                       + reference + '&report=' + suffix + '&retmode=text')
                payload = get(url)
                if (suffix == 'fasta' and not payload.startswith(b'>')) or (suffix == 'gff3' and b'##gff-version 3' not in payload[:100]):
                    raise RuntimeError('Unexpected NCBI public reference content')
                name = reference + '.' + suffix
                (source / name).write_bytes(payload)
                manifest.append(dict(reference_id=reference, name=name, url=url,
                                     bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest(),
                                     role='NCBI public versioned HCM transcript; matched to existing MANE1.4 summary'))
            reference_ids.append(reference)
    for name in ['LICENSE', 'tests/commons.py', 'tests/conftest.py']:
        url = f'https://raw.githubusercontent.com/mutalyzer/mutalyzer/{COMMIT}/{name}'
        try:
            payload = get(url)
        except Exception:
            if name == 'LICENSE':
                url = f'https://raw.githubusercontent.com/mutalyzer/mutalyzer/{COMMIT}/LICENSE.md'
                payload = get(url)
            else:
                raise
        target = root / ('upstream-' + name.replace('/', '-'))
        target.write_bytes(payload)
        manifest.append(dict(name=target.name, url=url, bytes=len(payload),
                             sha256=hashlib.sha256(payload).hexdigest()))
    env = root / 'venv'
    venv.EnvBuilder(with_pip=True).create(env)
    python = env / ('Scripts/python.exe' if sys.platform == 'win32' else 'bin/python')
    pip_metadata_url = 'https://pypi.org/pypi/pip/25.3/json'
    pip_metadata = json.loads(get(pip_metadata_url))
    pip_wheel = next(entry for entry in pip_metadata['urls'] if entry['filename'].endswith('.whl'))
    payload = get(pip_wheel['url'], limit=5_000_000)
    if hashlib.sha256(payload).hexdigest() != pip_wheel['digests']['sha256']:
        raise RuntimeError('Publisher pip wheel digest mismatch')
    pip_path = root / pip_wheel['filename']
    pip_path.write_bytes(payload)
    manifest.append(dict(name=pip_path.name, url=pip_wheel['url'],
                         metadata_url=pip_metadata_url, bytes=len(payload),
                         sha256=pip_wheel['digests']['sha256']))
    subprocess.run([str(python), '-m', 'pip', 'install', '--no-deps', str(pip_path)],
                   check=True, timeout=90)
    report = root / 'pip-install-report.json'
    subprocess.run([str(python), '-m', 'pip', 'install', '--disable-pip-version-check',
                    '--index-url', 'https://pypi.org/simple', 'setuptools==80.9.0', 'wheel==0.45.1'],
                   check=True, timeout=90)
    subprocess.run([str(python), '-m', 'pip', 'install', '--disable-pip-version-check',
                    '--index-url', 'https://pypi.org/simple', '--no-build-isolation',
                    '--require-hashes', '-r', str(args.lock), '--report', str(report)],
                   check=True, timeout=240)
    installed = json.loads(report.read_text())
    lines = []
    for entry in installed['install']:
        metadata = entry['metadata']
        sha = entry['download_info']['archive_info']['hashes']['sha256']
        url = entry['download_info']['url']
        if not url.startswith('https://files.pythonhosted.org/'):
            raise RuntimeError('Unexpected package file host')
        lines.append(f"{metadata['name']}=={metadata['version']} --hash=sha256:{sha}")
    (root / 'requirements.lock.txt').write_text('\n'.join(sorted(lines)) + '\n')
    (root / 'public-source-manifest.json').write_text(json.dumps(dict(
        schema_version=1, tool='Mutalyzer', version=VERSION, upstream_commit=COMMIT,
        fixture_reference_ids=reference_ids, mane14_public_pairs=mane_pairs, created_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        files=manifest, patient_data=False, isolated=True,
        compatibility='Public reference normalization proof; selected MANE1.4 HCM IDs optionally included; genomic-to-transcript mapping and full VEP115 compatibility require separate controls'), indent=2) + '\n')
    runner = Path(__file__).with_name('run_offline_hgvs.py')
    subprocess.run([str(python), str(runner), '--root', str(root), '--public-controls'],
                   check=True, timeout=90)
    print(json.dumps({'status': 'SETUP_AND_PUBLIC_CONTROLS_PASS', 'root': str(root)}))


if __name__ == '__main__':
    main()
