#!/usr/bin/env python3
"""Official DELLY public smoke proof and read-only BAM fitcheck; no patient SV calling."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import resource
import statistics
import subprocess
import time
import urllib.request

COMMIT = '3d505adc17afccb4443cf97dd01f78afc1d5a267'
VERSION = '2.7.0'
BINARY_SHA256 = '0a2a29c6b827b130948a0431ecba64942257e187c52d2b58ade014c032dab665'
DEADLINE = dt.datetime.fromisoformat('2026-11-06T03:46:48.046025+00:00')
# MANE1.4 public gene bodies plus 50kb; these are not capture targets.
WINDOWS = {'MYBPC3_flanks50kb': 'chr11:47281406-47402702',
           'MYH7_flanks50kb': 'chr14:23362740-23485660'}


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def command(argv, timeout=45, **kwargs):
    if dt.datetime.now(dt.timezone.utc) >= DEADLINE:
        raise RuntimeError('Original VM deadline reached')
    return subprocess.run(argv, timeout=timeout, check=True, capture_output=True,
                          text=True, **kwargs)


def public_proof(root):
    root.mkdir(parents=True, exist_ok=False)
    entries = [
        ('delly', f'https://github.com/dellytools/delly/releases/download/v{VERSION}/delly-v{VERSION}-linux-amd64', 9871504, BINARY_SHA256),
        ('ref.fa', f'https://raw.githubusercontent.com/dellytools/delly/{COMMIT}/example/ref.fa', 203339, None),
        ('ref.fa.fai', f'https://raw.githubusercontent.com/dellytools/delly/{COMMIT}/example/ref.fa.fai', 18, None),
        ('sr.bam', f'https://raw.githubusercontent.com/dellytools/delly/{COMMIT}/example/sr.bam', 1381351, None),
        ('sr.bam.bai', f'https://raw.githubusercontent.com/dellytools/delly/{COMMIT}/example/sr.bam.bai', 672, None),
        ('LICENSE', f'https://raw.githubusercontent.com/dellytools/delly/{COMMIT}/LICENSE', None, None),
    ]
    files = []
    received = 0
    for name, url, size, expected in entries:
        with urllib.request.urlopen(url, timeout=30) as response:
            payload = response.read(15_000_001)
        received += len(payload)
        if len(payload) > 15_000_000 or received > 50_000_000:
            raise RuntimeError('Public download cap exceeded')
        sha = hashlib.sha256(payload).hexdigest()
        if (size is not None and len(payload) != size) or (expected and sha != expected):
            raise RuntimeError('Official source identity mismatch')
        (root / name).write_bytes(payload)
        files.append(dict(name=name, url=url, bytes=len(payload), sha256=sha,
                          checksum_provenance='publisher release asset SHA256' if expected else 'own observation over HTTPS, pinned upstream commit'))
    binary = root / 'delly'
    binary.chmod(0o700)
    version_result = command([str(binary), '--version'])
    version = (version_result.stdout + version_result.stderr).strip()
    if VERSION not in version:
        raise RuntimeError('Executable version output mismatch')
    help_result = command([str(binary), 'sr', '--help'])
    (root / 'sr-help.txt').write_text(help_result.stdout + help_result.stderr)
    command(['samtools', 'quickcheck', '-v', str(root / 'sr.bam')])
    env = {**os.environ, 'OMP_NUM_THREADS': '2'}
    started = time.monotonic()
    run = command([str(binary), 'sr', '-g', str(root / 'ref.fa'), '-o', str(root / 'public-sr.bcf'), str(root / 'sr.bam')], timeout=180, env=env)
    (root / 'public-sr.log').write_text(run.stdout + run.stderr)
    records = command(['bcftools', 'view', '-H', str(root / 'public-sr.bcf')]).stdout.splitlines()
    receipt = dict(schema_version=1, status='PUBLIC_CLI_FIXTURE_SMOKE_PASS',
                   version=version, public_files=files, received_bytes=received,
                   command=['delly', 'sr', '-g', 'ref.fa', '-o', 'public-sr.bcf', 'sr.bam'],
                   public_output_records=len(records), elapsed_seconds=round(time.monotonic()-started, 3),
                   no_patient_data=True, patient_sv_caller_launched=False,
                   golden_truth_accuracy_benchmark=False, clinical_detection_sensitivity_validated=False)
    save(root / 'public-proof-receipt.json', receipt)
    print(json.dumps(receipt))


def fitcheck(cfg_path, output):
    cfg = json.loads(cfg_path.read_text())
    output.mkdir(parents=True, exist_ok=False)
    results = []
    overall = time.monotonic()
    for sample in cfg['samples']:
        bam = Path(sample['bam'])
        before = (bam.stat().st_size, bam.stat().st_mtime_ns)
        receipt = json.loads(Path(sample['receipt']).read_text())
        if not receipt.get('completed_utc') or receipt.get('sample') != sample['sample_id']:
            raise RuntimeError('Incomplete or mismatched new BAM receipt')
        command(['samtools', 'quickcheck', '-v', str(bam)])
        windows = {}
        for name, region in WINDOWS.items():
            if time.monotonic()-overall > 900:
                raise RuntimeError('Read-only fitcheck wall-time ceiling reached')
            view = command(['samtools', 'view', str(bam), region], timeout=45).stdout
            lines = view.splitlines()
            counts = dict(records=len(lines), primary=0, paired=0, proper_pair=0,
                          duplicate=0, secondary=0, supplementary=0, mapq_lt20=0,
                          softclipped=0, with_SA=0, mapped_mate=0, positive_tlen_records=0)
            tlen = []
            for line in lines:
                fields = line.split('\t');flag=int(fields[1]);mq=int(fields[4])
                for key, bit in [('paired',1),('proper_pair',2),('duplicate',1024),('secondary',256),('supplementary',2048)]:
                    counts[key] += bool(flag & bit)
                counts['primary'] += not bool(flag & (256|2048))
                counts['mapq_lt20'] += mq < 20
                counts['softclipped'] += 'S' in fields[5]
                counts['with_SA'] += any(v.startswith('SA:Z:') for v in fields[11:])
                counts['mapped_mate'] += bool(flag & 1) and not bool(flag & 8)
                length=int(fields[8])
                if length > 0 and not flag & (256|2048|1024):tlen.append(length)
            counts['positive_tlen_records']=len(tlen)
            # This measures the two public windows in the complete indexed BAM, never a mate subset.
            depth = command(['samtools','depth','-aa','-s','-q','0','-Q','20','-G','3844','-r',region,str(bam)], timeout=60).stdout
            depths=[int(line.split('\t')[2]) for line in depth.splitlines()]
            window_length=int(region.split('-')[-1])-int(region.split(':')[-1].split('-')[0])+1
            windows[name]=dict(region=region, read_geometry=counts,
                               positive_template_length_median=statistics.median(tlen) if tlen else None,
                               depth_positions=len(depths), expected_window_positions=window_length,
                               mean_depth=sum(depths)/len(depths) if depths else None,
                               median_depth=statistics.median(depths) if depths else None,
                               fraction_ge1=sum(d>=1 for d in depths)/len(depths) if depths else None,
                               fraction_ge10=sum(d>=10 for d in depths)/len(depths) if depths else None,
                               fraction_ge20=sum(d>=20 for d in depths)/len(depths) if depths else None)
        if before != (bam.stat().st_size,bam.stat().st_mtime_ns):
            raise RuntimeError('Source BAM changed during read-only inspection')
        results.append(dict(local_sample_id=sample['sample_id'], source_bam=str(bam),
                            receipt_sha256=hashlib.sha256(Path(sample['receipt']).read_bytes()).hexdigest(),
                            unchanged_stat=True, windows=windows))
    result=dict(schema_version=1,status='READ_ONLY_FITCHECK_COMPLETE',samples=results,
                public_windows=WINDOWS,source_bam_subset_created=False,patient_sv_caller_launched=False,
                target_bed_is_capture_design=False,uniform_wgs_proven=False,
                clinical_negative_result=False,elapsed_seconds=round(time.monotonic()-overall,3),
                limitations=['Only two public windows assessed; whole-genome sensitivity unassessed',
                             'Regional template-length distribution is capture/selection biased',
                             'CDS depth does not establish intronic breakpoint coverage',
                             'Actual capture BED and matching normal reference still required for targeted CNV design'])
    save(output/'fitcheck.private.json',result)
    print(json.dumps({'status':result['status'],'samples':len(results),'elapsed_seconds':result['elapsed_seconds'],'patient_sv_caller_launched':False}))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--public-proof-root',type=Path)
    group.add_argument('--fitcheck-config',type=Path)
    p.add_argument('--private-output',type=Path)
    a=p.parse_args()
    os.umask(0o077)
    os.nice(10)
    if hasattr(os,'sched_getaffinity'):
        os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
    resource.setrlimit(resource.RLIMIT_AS,(4*2**30,4*2**30))
    if dt.datetime.now(dt.timezone.utc)>=DEADLINE:raise SystemExit('Original deadline reached')
    if a.public_proof_root:public_proof(a.public_proof_root)
    else:
        if not a.private_output:raise SystemExit('--private-output is required for fitcheck')
        fitcheck(a.fitcheck_config,a.private_output)


if __name__=='__main__':main()
