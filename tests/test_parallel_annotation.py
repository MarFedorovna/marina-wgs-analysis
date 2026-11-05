"""Synthetic-only chromosome scatter integrity tests; no Docker or private data."""
import gzip
import importlib.util
import datetime as dt
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/parallel_annotation.py'
spec = importlib.util.spec_from_file_location('parallel_annotation', SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

HEADER = ('##fileformat=VCFv4.2\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPUBLIC_A\tPUBLIC_B\tPUBLIC_C\tPUBLIC_D\n')


def record(chrom, pos=1):
    return f'{chrom}\t{pos}\t.\tA\tC\t100\tPASS\tSOURCE=public\tGT:FT\t0/1:PASS\t0/0:PASS\t./.:NoCall\t0/1:LowGQ\n'


class ParallelAnnotationTests(unittest.TestCase):
    def fixture(self, root):
        source = root / 'input.vcf'
        source.write_text(HEADER + record('chr2') + record('chr2', 2) + record('chr1'))
        shards = m.partition(source, root / 'shards')
        for shard in shards:
            output = Path(shard['input']).with_name(shard['chromosome'] + '.annotated.vcf.gz')
            with gzip.open(shard['input'], 'rt') as inp, gzip.open(output, 'wt') as out:
                for line in inp:
                    out.write(line if line.startswith('#') else line.replace('SOURCE=public', 'SOURCE=public;CSQ=C|missense_variant'))
            shard['output'] = str(output)
        return source, shards

    def test_exact_input_block_order_and_all_sample_states_survive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, shards = self.fixture(root)
            out = root / 'merged.vcf'
            self.assertEqual(m.merge_shards(source, shards, out), 3)
            self.assertEqual([s['chromosome'] for s in shards], ['chr2', 'chr1'])
            self.assertIn('./.:NoCall\t0/1:LowGQ', out.read_text())
            self.assertEqual(out.read_text().count('CSQ=C|missense_variant'), 3)

    def test_ref_alt_filter_format_genotype_and_info_changes_fail(self):
        replacements = [('\tA\tC\t', '\tA\tG\t'), ('\t100\tPASS\t', '\t100\tLowQUAL\t'),
                        ('GT:FT', 'GT'), ('./.:NoCall', '0/0:PASS'), ('SOURCE=public;', '')]
        for old, new in replacements:
            with self.subTest(old=old), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source, shards = self.fixture(root)
                path = shards[0]['output']
                with gzip.open(path, 'rt') as stream:
                    text = stream.read()
                with gzip.open(path, 'wt') as stream:
                    stream.write(text.replace(old, new))
                with self.assertRaises(ValueError):
                    m.merge_shards(source, shards, root / 'merged.vcf')

    def test_missing_extra_reordered_records_and_schema_changes_fail(self):
        for transform in (lambda lines: lines[:-1], lambda lines: lines + [lines[-1]],
                          lambda lines: lines[:-2] + list(reversed(lines[-2:])),
                          lambda lines: ['##INFO=<ID=CSQ,Number=.,Type=String,Description="different">\n'] + lines):
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                source, shards = self.fixture(root)
                path = shards[0]['output']
                with gzip.open(path, 'rt') as stream:
                    lines = stream.readlines()
                with gzip.open(path, 'wt') as stream:
                    stream.writelines(transform(lines))
                with self.assertRaises(ValueError):
                    m.merge_shards(source, shards, root / 'merged.vcf')

    def test_repeated_contig_and_non_autosomes_refused(self):
        for text in (record('chr1') + record('chr2') + record('chr1', 2), record('chrX'), record('chrM')):
            with tempfile.TemporaryDirectory() as tmp:
                source = Path(tmp) / 'input.vcf'
                source.write_text(HEADER + text)
                with self.assertRaises(ValueError):
                    m.partition(source, Path(tmp) / 'shards')

    def test_only_explicit_execution_header_provenance_is_normalized(self):
        self.assertEqual(m.header_key('##VEP="v115.2" time="one" cache="115"\n'),
                         m.header_key('##VEP="v115.2" time="two" cache="115"\n'))
        self.assertNotEqual(m.header_key('##VEP="v115.2" time="one" cache="115"\n'),
                            m.header_key('##VEP="v115.2" time="one" cache="116"\n'))
        self.assertNotEqual(m.header_key('##INFO=<ID=CSQ,Description="a">\n'),
                            m.header_key('##INFO=<ID=CSQ,Description="b">\n'))

    def handoff_fixture(self, root, completed):
        rule = dict(producer_sha256='reviewed', producer_path=str(SCRIPT.resolve()),
                    binding_path='binding', binding_sha256='reviewed', proof_path='proof', proof_sha256='reviewed',
                    deadline_utc='future', old_processes=[dict(pid=101, start_ticks='old', argv=['controller']),
                                                        dict(pid=102, start_ticks='old', argv=['queue'])])
        binding = dict(annotation_config={'path': 'config'}, output_root=str(root / 'primary'),
                       budget={'vm_created_utc': 'created', 'hourly_usd': 2.39, 'planned_vm_usd': 40})
        deadline = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=1)
        rule['deadline_utc'] = deadline.isoformat()
        queue = SimpleNamespace(budget_deadline=lambda *a: deadline,
                                verify_stage=lambda path: path.name in completed)
        args = SimpleNamespace(rule=root / 'rule.json', reviewed_rule_sha256='reviewed', execute=True)
        return args, rule, binding, queue

    def test_completed_primary_annotation_or_report_is_never_stopped(self):
        for completed in ({'06_annotate'}, {'06_annotate', '07_report'}):
            with self.subTest(completed=completed), tempfile.TemporaryDirectory() as tmp:
                args, rule, binding, queue = self.handoff_fixture(Path(tmp), completed)
                with patch.object(m, 'sha', return_value='reviewed'), \
                     patch.object(m, 'load', side_effect=lambda p: rule if p == args.rule else {}), \
                     patch.object(m, 'frozen', return_value=(binding, queue, None, None, None)), \
                     patch.object(m, 'verify_proof'), patch.object(m.os, 'kill') as kill:
                    self.assertIn('COMPLETE', m.run(args)['status'])
                    kill.assert_not_called()

    def test_pid_reuse_refuses_handoff_before_any_signal(self):
        with tempfile.TemporaryDirectory() as tmp:
            args, rule, binding, queue = self.handoff_fixture(Path(tmp), set())
            fake_fcntl = SimpleNamespace(LOCK_EX=1, LOCK_NB=2, flock=lambda *a: None)
            with patch.object(m, 'sha', return_value='reviewed'), \
                 patch.object(m, 'load', side_effect=lambda p: rule if p == args.rule else {}), \
                 patch.object(m, 'frozen', return_value=(binding, queue, None, None, None)), \
                 patch.object(m, 'verify_proof'), patch.object(m, 'guard_deadline'), \
                 patch.object(m, 'process_identity', return_value=dict(pid=101, start_ticks='REUSED', argv=['other'])), \
                 patch.dict('sys.modules', fcntl=fake_fcntl), patch.object(m.os, 'kill') as kill:
                with self.assertRaisesRegex(ValueError, 'identity changed'):
                    m.run(args)
                kill.assert_not_called()


if __name__ == '__main__':
    unittest.main()
