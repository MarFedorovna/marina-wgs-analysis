"""Small synthetic checks; no patient data, Docker or GPU required."""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('calling_stages', Path(__file__).resolve().parents[1] / 'scripts/calling_stages.py')
stages = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stages)

HEADER = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n'
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n'
          '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH_A\tSYNTH_B\n')
SNP_INFO = 'QD=10;FS=1;SOR=1;MQ=60;MQRankSum=0;ReadPosRankSum=0'


class FakeRuntime:
    dry_run = True
    execute = False

    def __init__(self, root):
        self.root, self.commands, self.parameters = Path(root), [], []

    def run(self, argv, logname='command', capture=False, **kwargs):
        self.commands.append(list(map(str, argv)))
        return ''

    def docker(self, image, args, mounts=None, gpu=False, cpus=1, memory_gib=4):
        return ['docker', 'run', '--gpu-test' if gpu else '--cpu-test', image, *args]

    def checkpoint(self, key, input_paths, parameters, producer):
        directory = self.root / 'nested' / key
        directory.mkdir(parents=True)
        self.parameters.append(parameters)
        producer(directory)
        return directory

    def stage_path(self, stage):
        return self.root / stage


class CallingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def filter(self, records, header=HEADER, policy=None):
        source, tagged, selected = [self.root / x for x in ('source.vcf', 'tagged.vcf', 'selected.vcf')]
        source.write_text(header + ''.join(records))
        result = stages.soft_filter_vcf(source, tagged, selected, policy)
        read = lambda p: [x.split('\t') for x in p.read_text().splitlines() if not x.startswith('#')]
        return result, read(tagged), read(selected)

    def fixture(self, sex=False):
        reference = self.root / 'reference.fasta'
        reference.write_text('>chr1\n' + 'A' * 100 + '\n' + ('>chrX\n' + 'A' * 100 + '\n' if sex else ''))
        Path(str(reference) + '.fai').write_text('chr1\t100\t6\t100\t101\n' + ('chrX\t100\t113\t100\t101\n' if sex else ''))
        reference.with_suffix('.dict').write_text('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n')
        intervals = self.root / 'regions.bed'
        intervals.write_text(('chrX' if sex else 'chr1') + '\t0\t100\n')
        samples = self.root / 'samples.tsv'
        samples.write_text('sample_id\nSYNTH_A\nSYNTH_B\n')
        manifest = self.root / 'alignments.tsv'
        rows = []
        for sid in ('SYNTH_A', 'SYNTH_B'):
            directory = self.root / ('alignment-' + sid)
            directory.mkdir()
            bam = directory / 'marked.bam'
            bam.write_bytes(b'synthetic-placeholder-bam')
            Path(str(bam) + '.bai').write_bytes(b'synthetic-index')
            receipt = directory / 'receipt.json'
            receipt.write_text(json.dumps({'sample': sid, 'config': {'reference': {'sha256': stages.sha256(reference)}},
                'outputs': {p.name: stages.sha256(p) for p in (bam, Path(str(bam) + '.bai'))}}))
            rows.append({'sample_id': sid, 'bam': str(bam), 'receipt': str(receipt)})
        with manifest.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        return {'reference': {'fasta': str(reference), 'sha256': stages.sha256(reference), 'assembly_label': 'SYNTHETIC'},
                'samples_manifest': str(samples), 'alignment_manifest': str(manifest),
                'callsets': [{'name': 'synthetic', 'intervals': str(intervals), 'intervals_sha256': stages.sha256(intervals),
                             'ploidy': 2, 'samples': ['SYNTH_A', 'SYNTH_B'], 'scope': 'autosomes', 'excluded_contigs': ['chrX', 'chrY', 'chrM']}],
                'tools': {'parabricks_image': 'nvcr.io/nvidia/clara/clara-parabricks@sha256:' + 'a' * 64,
                          'gatk_image': 'broadinstitute/gatk@sha256:' + 'b' * 64, 'bcftools': 'bcftools'},
                'output_root': str(self.root / 'outputs')}

    def test_all_rows_and_gt_preserved_with_individual_ft(self):
        record = 'chr1\t10\t.\tA\tC\t100\tPASS\t' + SNP_INFO + '\tGT:DP:GQ\t0/1:20:60\t1/1:2:5\n'
        summary, tagged, selected = self.filter([record])
        self.assertEqual(len(tagged), 1)
        self.assertEqual(len(selected), 1)
        self.assertEqual(tagged[0][9], '0/1:20:60:PASS')
        self.assertEqual(tagged[0][10], '1/1:2:5:LowDP;LowGQ')
        self.assertEqual(summary['counts']['all_records'], 1)

    def test_nocall_and_missing_confidence_fail_without_masking(self):
        _, tagged, selected = self.filter(['chr1\t10\t.\tA\tC\t100\tPASS\t' + SNP_INFO + '\tGT:DP:GQ\t./.:.:.\t0/.:20:40\n'])
        self.assertEqual(selected, [])
        self.assertTrue(tagged[0][9].startswith('./.:.:.:NoCall;LowDP;LowGQ'))
        self.assertEqual(tagged[0][10], '0/.:20:40:NoCall')

    def test_snp_indel_distinct_thresholds(self):
        snp = 'chr1\t10\t.\tA\tC\t100\tPASS\t' + SNP_INFO.replace('FS=1;', 'FS=100;') + '\tGT:DP:GQ\t0/1:20:60\t0/0:20:60\n'
        indel = 'chr1\t20\t.\tA\tAC\t100\tPASS\tQD=10;FS=100;ReadPosRankSum=-10\tGT:DP:GQ\t0/1:20:60\t0/0:20:60\n'
        _, tagged, selected = self.filter([snp, indel])
        self.assertEqual(tagged[0][6], 'Fail_FS')
        self.assertEqual(tagged[1][6], 'PASS')
        self.assertEqual(len(selected), 1)

    def test_missing_site_metrics_and_existing_filter_retained(self):
        _, tagged, selected = self.filter(['chr1\t10\t.\tA\tC\t.\tOriginalFailure\tQD=10\tGT:DP:GQ\t0/1:20:60\t0/0:20:60\n'])
        self.assertEqual(selected, [])
        self.assertEqual(set(tagged[0][6].split(';')), {'OriginalFailure', 'LowQUAL', 'MissingAnnotation'})

    def test_existing_ft_and_haploid_gt_preserved(self):
        header = HEADER.replace('#CHROM', '##FORMAT=<ID=FT,Number=1,Type=String,Description="FT">\n#CHROM')
        _, tagged, selected = self.filter(['chr1\t10\t.\tA\tC\t100\tPASS\t' + SNP_INFO + '\tGT:DP:GQ:FT\t1:20:60:OriginalGT\t0:20:60:PASS\n'], header)
        self.assertEqual(tagged[0][9], '1:20:60:OriginalGT')
        self.assertEqual(tagged[0][8].count('FT'), 1)
        self.assertEqual(selected, [])

    def test_filter_threshold_equality_passes(self):
        _, tagged, selected = self.filter(['chr1\t10\t.\tA\tC\t30\tPASS\tQD=2;FS=60;SOR=3;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8\tGT:DP:GQ\t0/1:10:20\t0/0:10:20\n'])
        self.assertEqual(tagged[0][6], 'PASS')
        self.assertEqual(len(selected), 1)

    def test_scope_excluded_sample_remains_missing_never_hom_ref(self):
        source, target = self.root / 'source.vcf', self.root / 'expanded.vcf'
        source.write_text(HEADER.replace('\tSYNTH_B', '') +
            'chr1\t10\t.\tA\tC\t100\tPASS\t' + SNP_INFO + '\tGT:DP:GQ\t1:20:60\n')
        stages.expand_sample_columns(source, target, ['SYNTH_B', 'SYNTH_A'])
        record = [line for line in target.read_text().splitlines() if not line.startswith('#')][0].split('\t')
        self.assertEqual(record[9], '.:.:.')
        self.assertEqual(record[10], '1:20:60')
        self.assertIn('##marina_scope_excluded_samples=["SYNTH_B"]', target.read_text())

    def test_sex_contig_requires_explicit_sample_ploidy(self):
        cfg = self.fixture(sex=True)
        with self.assertRaisesRegex(ValueError, 'Autosome'):
            stages.context(cfg)
        cfg['callsets'][0]['scope'] = 'explicit_nonPAR_X'
        with self.assertRaisesRegex(ValueError, 'per-sample'):
            stages.context(cfg)
        cfg['callsets'][0]['per_sample_ploidy'] = {'SYNTH_A': 1, 'SYNTH_B': 2}
        stages.context(cfg)

    def test_reference_and_interval_hash_checks(self):
        cfg = self.fixture()
        Path(cfg['callsets'][0]['intervals']).write_text('chr1\t0\t90\n')
        with self.assertRaisesRegex(ValueError, 'intervals SHA256'):
            stages.context(cfg)

    def test_bam_reuse_checks_original_receipt_not_fake_import(self):
        cfg = self.fixture()
        rows = stages.read_tsv(cfg['alignment_manifest'])
        Path(rows[0]['bam']).write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'BAM source receipt SHA256'):
            stages.alignment_rows(cfg, ['SYNTH_A'])

    def test_call_constructs_gvcf_from_existing_bam(self):
        cfg = self.fixture()
        workdir = self.root / 'call'
        workdir.mkdir()
        runtime = FakeRuntime(self.root)
        result = stages.run_stage('call', cfg, workdir, runtime)
        commands = [c for c in runtime.commands if 'haplotypecaller' in c]
        self.assertEqual(len(commands), 2)
        self.assertTrue(all('--gvcf' in c and '--interval-file' in c and '--in-recal-file' not in c for c in commands))
        self.assertEqual(result['metadata']['gvcfs'], 2)
        self.assertTrue(result['metadata']['source_alignments_reused_without_receipt_changes'])

    def test_joint_partitions_mixed_ploidy_without_reference_imputation(self):
        cfg = self.fixture(sex=True)
        cfg['callsets'][0].update(scope='explicit_nonPAR_X', per_sample_ploidy={'SYNTH_A': 1, 'SYNTH_B': 2})
        runtime = FakeRuntime(self.root)
        call_dir = self.root / 'call'
        call_dir.mkdir()
        stages.run_stage('call', cfg, call_dir, runtime)
        joint_dir = self.root / 'joint'
        joint_dir.mkdir()
        result = stages.run_stage('joint', cfg, joint_dir, runtime)
        commands = [c for c in runtime.commands if 'GenotypeGVCFs' in c]
        self.assertEqual({c[c.index('--sample-ploidy') + 1] for c in commands}, {'1', '2'})
        self.assertTrue(any(c[:2] == ['bcftools', 'merge'] for c in runtime.commands))
        self.assertTrue(all('--missing-to-ref' not in c for c in runtime.commands))
        self.assertIn('cohort.raw.vcf.gz', result['outputs'])

    @unittest.skipUnless(os.environ.get('MARINA_SYNTHETIC_CPU_EXECUTE') == '1',
                         'Explicit MARINA_SYNTHETIC_CPU_EXECUTE=1 required for Docker CPU integration')
    def test_actual_runtime_gatk_joint_filter_cpu_only(self):
        """Opt-in real CLI smoke: four CPUs/16 GiB, artificial GVCFs, zero GPU access."""
        self.assertTrue(shutil.which('docker') and shutil.which('bcftools'))
        common_spec = importlib.util.spec_from_file_location('workflow_common',
            Path(__file__).resolve().parents[1] / 'scripts/workflow_common.py')
        common = importlib.util.module_from_spec(common_spec)
        common_spec.loader.exec_module(common)
        cfg = self.fixture()
        cfg['tools']['gatk_image'] = 'broadinstitute/gatk@sha256:e7996ba655225c1cde0a1faec6a113e217758310af2cf99b00d61dae8ec6e9f2'
        cfg['callsets'][0].update(scope='explicit_synthetic_mixed_ploidy',
                                  per_sample_ploidy={'SYNTH_A': 1, 'SYNTH_B': 2})
        cfg['resources'] = {'joint_cpus': 4, 'joint_memory_gib': 16, 'joint_java_heap_gib': 8,
                            'min_free_gib': 0, 'max_hours': 0.1}
        runtime = common.Runtime(cfg, execute=True, workdir=self.root, stage='joint')
        source_dir = self.root / 'synthetic-gvcfs'
        source_dir.mkdir()
        header = ('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
                  '##ALT=<ID=NON_REF,Description="Any alternative allele">\n'
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
                  '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        for key, kind, cardinality in [('GT', 'String', '1'), ('DP', 'Integer', '1'),
                                       ('AD', 'Integer', 'R'), ('GQ', 'Integer', '1'),
                                       ('MIN_DP', 'Integer', '1'), ('PL', 'Integer', 'G')]:
            header += f'##FORMAT=<ID={key},Number={cardinality},Type={kind},Description="Synthetic {key}">\n'
        rows = []
        for sid in ['SYNTH_A', 'SYNTH_B']:
            plain = source_dir / (sid + '.g.vcf')
            gvcf = Path(str(plain) + '.gz')
            gt = ('1:0,20,0:20:60:100,0,200' if sid == 'SYNTH_A'
                  else '1/1:0,20,0:20:60:200,60,0,300,300,400')
            block = '0:20:60:20:0,600' if sid == 'SYNTH_A' else '0/0:20:60:20:0,60,600'
            plain.write_text(header + '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sid + '\n'
                + 'chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=49\tGT:DP:GQ:MIN_DP:PL\t' + block + '\n'
                + 'chr1\t50\t.\tA\tC,<NON_REF>\t100\tPASS\tDP=20\tGT:AD:DP:GQ:PL\t' + gt + '\n'
                + 'chr1\t51\t.\tA\t<NON_REF>\t.\t.\tEND=100\tGT:DP:GQ:MIN_DP:PL\t' + block + '\n')
            runtime.run(['bcftools', 'view', '-Oz', '-o', str(gvcf), str(plain)], 'synthetic-compress')
            runtime.run(['bcftools', 'index', '--tbi', str(gvcf)], 'synthetic-index')
            c = cfg['callsets'][0]
            rows.append(dict(sample_id=sid, callset=c['name'], ploidy=str(c['per_sample_ploidy'][sid]), gvcf=str(gvcf),
                intervals=c['intervals'], intervals_sha256=c['intervals_sha256'], reference_sha256=cfg['reference']['sha256']))
        call_dir = runtime.stage_path('call')
        call_dir.mkdir(parents=True)
        manifest = call_dir / 'gvcfs.tsv'
        with manifest.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        child_commands = []
        for stage in ('joint', 'filter'):
            current = runtime.child(self.root, stage)
            inputs = [manifest, cfg['reference']['fasta']] + [r['gvcf'] for r in rows]
            if stage == 'filter':
                inputs += [str(runtime.stage_path('joint') / 'cohort.raw.vcf.gz')]
            producer = lambda directory, s=stage, r=current: stages.run_stage(s, cfg, directory, r)
            directory = current.run_stage(stage, inputs, {'synthetic_only': True}, producer)
            self.assertEqual(json.loads((directory / 'workflow-receipt.json').read_text())['status'], 'PASS')
            # The real shared Runtime must verify completed outputs and skip the producer.
            self.assertEqual(current.run_stage(stage, inputs, {'synthetic_only': True}, producer), directory)
            child_commands.extend(current.commands)
        raw = runtime.stage_path('joint') / 'cohort.raw.vcf.gz'
        count = runtime.run(['bcftools', 'index', '-n', str(raw)], 'synthetic-count', capture=True)
        self.assertEqual(int(count), 1)
        genotypes = runtime.run(['bcftools', 'query', '-f', '[%GT\t]\n', str(raw)], 'synthetic-ploidy', capture=True).strip().split()
        self.assertEqual(genotypes, ['1', '1/1'])
        filtered = runtime.stage_path('filter')
        policy = json.loads((filtered / 'filtering-policy.json').read_text())
        self.assertEqual(policy['counts']['all_records'], 1)
        self.assertIn('cohort.filtered.vcf.gz', [p.name for p in filtered.iterdir()])
        self.assertFalse(any('--gpus' in command for command in runtime.commands + child_commands))


if __name__ == '__main__':
    unittest.main()
