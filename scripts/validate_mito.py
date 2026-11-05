#!/usr/bin/env python3
"""Explicit opt-in synthetic circular-rCRS pipeline smoke; never uses patient reads."""
import argparse
import csv
import datetime as dt
import json
import os
from pathlib import Path
import random
import subprocess

import mitochondria as mt
from calling_stages import require, sha256
from prepare_mito_resources import sequence
from workflow_common import json_write


def validate(source_cfg, resource_path, base):
    resources = json.loads(resource_path.read_text())
    base.mkdir(parents=True)
    (base / 'private').mkdir()
    sample = 'MITO_SYNTHETIC'
    sample_file = base / 'samples.tsv'
    sample_file.write_text('sample_id\n' + sample + '\n')
    manifest = base / 'alignments.tsv'
    cfg = {'data_root': str(base), 'output_root': str(base / 'output'), 'reference': source_cfg['reference'],
           'samples_manifest': str(sample_file), 'alignment_manifest': str(manifest), 'tools': {'bcftools': 'bcftools'},
           'resources': {'max_hours': 1, 'min_free_gib': 0}, 'budget': source_cfg['budget'],
           'mt_resource_manifest_sha256': sha256(resource_path), 'synthetic_only': True}
    runtime = mt.MitoRuntime(cfg, True, base, 'mitochondria')
    ref = sequence(resources['normal_fasta'])
    haplotype = list(ref)
    expected = {}
    for position in (100, 5000):
        alt = next(a for a in 'ACGT' if a != ref[position - 1])
        haplotype[position - 1] = alt
        expected[position] = alt
    haplotype = ''.join(haplotype)
    randomizer = random.Random(79123)
    r1, r2 = base / 'synthetic_R1.fastq', base / 'synthetic_R2.fastq'
    complement = str.maketrans('ACGTN', 'TGCAN')
    with r1.open('x') as first, r2.open('x') as second:
        for index in range(3000):
            start = randomizer.randrange(16569)
            template = haplotype if index % 2 == 0 else ref
            fragment = ''.join(template[(start + i) % 16569] for i in range(180))
            left, right = fragment[:100], fragment[-100:].translate(complement)[::-1]
            name = 'MITO_SYNTH_READ_' + str(index)
            first.write('@' + name + '\n' + left + '\n+\n' + 'I' * 100 + '\n')
            second.write('@' + name + '\n' + right + '\n+\n' + 'I' * 100 + '\n')
    original = base / 'original-alignment'
    original.mkdir()
    sam, bam = original / 'whole-reference.sam', original / 'marked.bam'
    runtime.run(['bwa', 'mem', '-Y', '-t', '2', '-R', '@RG\\tID:mitosynth\\tSM:' + sample +
                 '\\tLB:analysis_library\\tPU:analysis_unit_synthetic', cfg['reference']['fasta'], str(r1), str(r2)],
                'synthetic-whole-reference-alignment', stdout=sam)
    runtime.run(['samtools', 'sort', '-@', '2', '-o', str(bam), str(sam)], 'synthetic-sort')
    runtime.run(['samtools', 'index', str(bam)], 'synthetic-index')
    receipt = original / 'receipt.json'
    json_write(receipt, {'sample': sample, 'config': {'reference': cfg['reference']},
                'synthetic_fixture_not_patient_data': True,
                'method': 'Generated circular rCRS paired FASTQs with known100/5000 variants; BWA whole-genome alignment +samtools sort',
                'outputs': {p.name: sha256(p) for p in (bam, Path(str(bam) + '.bai'))}})
    manifest.write_text('sample_id\tbam\treceipt\n' + sample + '\t' + str(bam) + '\t' + str(receipt) + '\n')
    json_write(base / 'config.json', cfg)
    mt.run(cfg, resource_path, True, workers=1)
    target = Path(cfg['output_root']) / 'checkpoints' / sample / 'calls' / 'mitochondria.raw.vcf.gz'
    records = subprocess.check_output(['bcftools', 'query', '-f', '%POS\t%ALT\n', str(target)], text=True).splitlines()
    actual = {int(line.split('\t')[0]): line.split('\t')[1].split(',') for line in records}
    require(all(p in actual and alt in actual[p] for p, alt in expected.items()), 'Known normal/control synthetic variants were not recovered.')
    report = {'status': 'PASS', 'synthetic_only': True, 'patient_reads_used': False,
              'circular_reference_length': 16569, 'known_control_variant_recovered': 100,
              'known_normal_variant_recovered': 5000, 'caller': 'Mutect2 mitochondria-mode',
              'whole_genome_NUMT_image_used': True, 'haplochecker_executed_or_explicit_unavailable': True,
              'clinical_detection_limit_validated': False, 'fixture': str(base)}
    json_write(base / 'validation-summary.json', report)
    print(json.dumps(report), flush=True)
    return base


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workflow-config', type=Path, required=True)
    parser.add_argument('--resources', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args(); os.umask(0o077)
    require(args.execute, 'Explicit --execute required for synthetic smoke.')
    base = args.output or Path('/srv/marina/scratch') / ('mito-synthetic-' + dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    validate(json.loads(args.workflow_config.read_text()), args.resources.resolve(), base.resolve())


if __name__ == '__main__':
    main()
