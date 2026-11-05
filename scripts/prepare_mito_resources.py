#!/usr/bin/env python3
"""Prepare pinned public mitochondrial references and a full-genome NUMT index."""
import argparse
import base64
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import urllib.request
import uuid

from calling_stages import mounts, require, sha256
from workflow_common import Runtime, json_write, promote_exclusive


def sequence(path):
    return ''.join(line.strip() for line in Path(path).read_text().splitlines() if not line.startswith('>')).upper()


def prepare(cfg, lock_path, root, execute):
    lock = json.loads(Path(lock_path).read_text())
    reference = Path(cfg['reference']['fasta']).resolve(strict=True)
    require(sha256(reference) == cfg['reference']['sha256'], 'Whole-genome reference SHA256 mismatch.')
    if not execute:
        print(json.dumps({'status': 'PLAN', 'public_assets': len(lock['objects']), 'whole_reference_sha256': cfg['reference']['sha256'],
                          'build_image_from_existing_indexes': True, 'patient_data_used': False}))
        return
    root.mkdir(parents=True, exist_ok=True)
    with (root / '.prepare.lock').open('a') as guard:
        fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest_path = root / 'resource-manifest.json'
        if manifest_path.exists():
            saved = json.loads(manifest_path.read_text())
            require(saved['status'] == 'PASS' and saved['whole_reference_sha256'] == cfg['reference']['sha256']
                    and saved['lock_sha256'] == sha256(lock_path), 'Existing mitochondrial resource identity differs.')
            for name, expected in saved['files'].items():
                require(Path(name).stat().st_size == expected['bytes'] and sha256(name) == expected['sha256'], 'Prepared MT resource changed.')
            print(json.dumps({'status': 'VERIFIED_COMPLETE', 'manifest': str(manifest_path)}))
            return
        require(shutil.disk_usage(root).free >= 150 * 2**30, 'Less than 150GiB disk reserve.')
        bundle = root / 'reference-bundle'
        bundle.mkdir(exist_ok=True)
        for item in lock['objects']:
            require(Path(item['destination']).name == item['destination'], 'Unsafe locked resource basename.')
            target = bundle / item['destination']
            if not target.exists():
                temporary = bundle / ('.download-' + uuid.uuid4().hex)
                try:
                    with urllib.request.urlopen(item['url'], timeout=120) as response, temporary.open('xb') as stream:
                        total = 0
                        for block in iter(lambda: response.read(1024**2), b''):
                            total += len(block)
                            require(total <= item['bytes'], 'Public reference exceeds locked size.')
                            stream.write(block)
                    require(temporary.stat().st_size == item['bytes'] and sha256(temporary) == item['sha256'], 'Locked reference SHA mismatch.')
                    os.link(temporary, target)
                finally:
                    temporary.unlink(missing_ok=True)
            require(target.stat().st_size == item['bytes'] and sha256(target) == item['sha256'], 'Existing locked reference differs.')
            require(base64.b64encode(hashlib.md5(target.read_bytes()).digest()).decode() == item['publisher_md5_base64'], 'Publisher MD5 mismatch.')
        normal = bundle / 'Homo_sapiens_assembly38.chrM.fasta'
        shifted = bundle / 'Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta'
        seq = sequence(normal)
        require(len(seq) == 16569 and hashlib.md5(seq.encode()).hexdigest() == lock['normal_sequence_md5'], 'Normal reference is not verified rCRS.')
        require(sequence(shifted) == seq[8000:] + seq[:8000], 'Shifted reference rotation is not exactly 8000 bases.')
        dictionary = reference.with_suffix('.dict').read_text().splitlines()
        mt = [dict(x.split(':', 1) for x in line.split('\t')[1:]) for line in dictionary
              if line.startswith('@SQ\t') and '\tSN:chrM\t' in line]
        require(len(mt) == 1 and mt[0]['M5'] == lock['normal_sequence_md5'] and mt[0]['LN'] == '16569',
                'Whole-reference chrM does not match rCRS.')
        # Verify that existing BWA indexes belong to the verified whole-reference bundle.
        original = json.loads((reference.parent / 'reference-receipt.json').read_text())
        artifacts = {x['path']: x for x in original['artifacts']}
        index_paths = [Path(str(reference) + suffix) for suffix in ('.amb', '.ann', '.bwt', '.pac', '.sa', '.alt')]
        for path in index_paths:
            record = artifacts[path.name]
            require(path.stat().st_size == record['bytes'] and sha256(path) == record['sha256'], 'Original BWA index differs from its reference receipt.')
        runtime = Runtime(cfg, True, root / 'logs', 'prepare')
        runtime.workdir.mkdir(exist_ok=True)
        helper = Path(__file__).with_name('MitoIndexImage.java').resolve(strict=True)
        final_image_dir = root / 'whole-genome-index'
        image_dir = root / ('whole-genome-index.partial-' + uuid.uuid4().hex[:10])
        index_image = image_dir / 'whole-genome.img'
        jar = '/gatk/gatk-package-4.3.0.0-local.jar'
        compile_cmd = ['javac', '-cp', jar, '-d', str(image_dir), str(helper)]
        build_cmd = ['java', '-Xmx2g', '-cp', jar + ':' + str(image_dir), 'MitoIndexImage',
                     str(reference), str(index_image), str(reference) + '.fai']
        script = shlex.join(compile_cmd) + ' && ' + shlex.join(build_cmd)
        if final_image_dir.exists():
            seal = json.loads((final_image_dir / 'image-identity.json').read_text())
            index_image = final_image_dir / 'whole-genome.img'
            require(seal['whole_reference_sha256'] == cfg['reference']['sha256']
                    and seal['helper_sha256'] == sha256(helper) and sha256(index_image) == seal['image_sha256'],
                    'Existing whole-genome image has a different identity.')
        else:
            image_dir.mkdir()
            runtime.run(runtime.docker(lock['gatk_image'], ['sh', '-c', script],
                        mounts=mounts((reference.parent, True), (helper.parent, True), (image_dir, False)),
                        gpu=False, cpus=4, memory_gib=16), 'whole-genome-image', timeout=3600)
            require(index_image.is_file() and index_image.stat().st_size > 0, 'No BWA image produced.')
            json_write(image_dir / 'image-identity.json', {'whole_reference_sha256': cfg['reference']['sha256'],
                       'helper_sha256': sha256(helper), 'image_sha256': sha256(index_image),
                       'source_index_artifacts': {p.name: sha256(p) for p in index_paths}})
            promote_exclusive(image_dir, final_image_dir)
            index_image = final_image_dir / 'whole-genome.img'
        for image in (lock['gatk_image'], lock['haplochecker_image']):
            runtime.run(runtime.docker_prefix() + ['image', 'inspect', image], 'pinned-image-check', capture=True)
        files = [p for p in bundle.iterdir() if p.is_file()] + [index_image, helper]
        result = {'schema_version': 1, 'status': 'PASS', 'lock_sha256': sha256(lock_path),
                  'gatk_image': lock['gatk_image'], 'haplochecker_image': lock['haplochecker_image'],
                  'normal_fasta': str(normal), 'shifted_fasta': str(shifted),
                  'shift_back_chain': str(bundle / 'ShiftBack.chain'),
                  'blacklist_bed': str(bundle / 'blacklist_sites.hg38.chrM.bed'),
                  'selection_intervals': str(bundle / 'chrMWithFinalNuMTs.hg38.interval_list'),
                  'whole_genome_image': str(index_image), 'whole_reference_sha256': cfg['reference']['sha256'],
                  'normal_sequence_md5': lock['normal_sequence_md5'], 'shift_bases': 8000,
                  'image_creation': 'Existing verified BWA index files via GATK BwaMemIndex API; no FASTA reindexing',
                  'files': {str(p): {'bytes': p.stat().st_size, 'sha256': sha256(p)} for p in files}}
        # The directory name remains explicit in the manifest; only this PASS receipt declares readiness.
        json_write(manifest_path, result)
        print(json.dumps({'status': 'PASS', 'manifest': str(manifest_path), 'public_assets': len(lock['objects']), 'patient_data_used': False}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-config', type=Path, required=True)
    parser.add_argument('--lock', type=Path, default=Path(__file__).resolve().parents[1] / 'config/mitochondria.resources.lock.json')
    parser.add_argument('--root', type=Path, default=Path('/srv/marina/resources/mitochondria/gatk-4.3.0.0'))
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    os.umask(0o077)
    prepare(json.loads(args.workflow_config.read_text()), args.lock.resolve(), args.root.resolve(), args.execute)


if __name__ == '__main__':
    main()
