"""Synthetic-only packet, row-order, counter and read-only database checks."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/fast_vcf_evidence_candidate.py'
spec = importlib.util.spec_from_file_location('direct_packet_candidate', SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

SAMPLES = ['PUBLIC_A', 'PUBLIC_B', 'PUBLIC_C', 'PUBLIC_D']
HEADER = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(SAMPLES) + '\n'


def record(position, transcripts):
    return 'chr1\t' + str(position) + '\t.\tA\tC\t100\tPASS\tCSQ=' + ','.join(['C|missense_variant'] * transcripts) + '\tGT:FT\t0/1:PASS\t0/0:PASS\t./.:NoCall\t0/1:LowGQ\n'


class DirectPacketTests(unittest.TestCase):
    def test_balanced_batches_preserve_whole_allele_and_global_order(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); source = root / 'input.vcf'
            source.write_text(HEADER + record(1, 2) + record(2, 1) + record(3, 4), newline='')
            headers, jobs, rows = m.partition(source, root / 'jobs', SAMPLES, 8)
            self.assertEqual(rows, 28)
            self.assertEqual([(v['first_i'], v['last_i']) for v in jobs], [(1, 8), (9, 28)])
            self.assertEqual([v['records'] for v in jobs], [1, 2])
            joined = ''
            for job in jobs:
                with gzip.open(job['input'], 'rt') as stream:
                    joined += ''.join(line for line in stream if not line.startswith('#'))
            self.assertEqual(joined, record(1, 2) + record(2, 1) + record(3, 4))

    def test_oversize_single_allele_is_not_split_to_meet_batch_target(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); source = root / 'input.vcf'
            source.write_text(HEADER + record(1, 20) + record(2, 1), newline='')
            _, jobs, rows = m.partition(source, root / 'jobs', SAMPLES, 8)
            self.assertEqual([v['expected_rows'] for v in jobs], [80, 4])
            self.assertEqual(rows, 84)

    def test_codec_preserves_ordered_unicode_escaped_and_missing_strings(self):
        value = {'v': 1, 'rows': [{'i': 1, 'values': ['', '.', './.', '0|1', 'a\tb', '"quoted"', 'строка', 'one\ntwo']}],
                 'membership': {'unknown_or_filtered_samples': 'PUBLIC_C&PUBLIC_D'}}
        self.assertEqual(m.decode(m.encode(value)), value)
        with self.assertRaises(Exception):
            m.decode('bad-token!')

    def test_incremental_vector_digest_equals_canonical_tsv_without_loss(self):
        rows = [['PUBLIC_A', '1', 'a\tb'], ['PUBLIC_B', '2', '"quoted"'], ['PUBLIC_C', '', 'строка']]
        digest = m.VectorDigest()
        for row in rows:digest.row(row)
        self.assertEqual(digest.digest.hexdigest(), hashlib.sha256(b''.join(m.canonical_row_bytes(row) for row in rows)).hexdigest())

    def test_merkle_binds_field_order_batch_order_and_vector_digest(self):
        proof = dict(id=0, first_i=1, last_i=4, row_count=4, record_count=1, canonical_row_sha256='a' * 64)
        second = dict(proof, id=1, first_i=5, last_i=8)
        expected = m.merkle_digest(['sample', 'GT'], [proof, second])
        self.assertNotEqual(expected, m.merkle_digest(['GT', 'sample'], [proof, second]))
        self.assertNotEqual(expected, m.merkle_digest(['sample', 'GT'], [second, proof]))
        self.assertNotEqual(expected, m.merkle_digest(['sample', 'GT'], [dict(proof, canonical_row_sha256='b' * 64), second]))

    def test_counter_merge_sums_rows_but_never_duplicate_alpha_source_rows(self):
        counts = dict(records=1, transcript_sample_rows=4, missing_csq=0, nonsequence_records_preserved=0,
                      samples=SAMPLES, alphamissense_matched_source_rows=17,
                      alphamissense_transcript_sample_statuses={'not_covered_in_source; not_benign': 4})
        consequences = dict(schema_version=1, definitions={'count_unit': 'unique allele'},
                            cohort={'raw': {'unique_normalized_alleles': 1, 'nonexclusive_flags': {'missense': 1}}},
                            samples={'PUBLIC_A': {'raw': {'unique_normalized_alleles': 1}}})
        result, merged = m.merged_counts([dict(counts=counts, consequences=consequences)] * 2)
        self.assertEqual(result['records'], 2)
        self.assertEqual(result['transcript_sample_rows'], 8)
        self.assertEqual(result['alphamissense_matched_source_rows'], 17)
        self.assertEqual(merged['cohort']['raw']['unique_normalized_alleles'], 2)
        self.assertEqual(merged['definitions'], consequences['definitions'])

    def test_snapshot_database_connection_cannot_mutate_data(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'snapshot.sqlite'
            with sqlite3.connect(path) as connection:
                connection.execute('CREATE TABLE evidence (value TEXT)')
                connection.execute('INSERT INTO evidence VALUES ("public")')
            connection.close()
            before = path.read_bytes()
            connection = m.readonly_database(path)
            self.assertEqual(connection.execute('SELECT value FROM evidence').fetchone(), ('public',))
            with self.assertRaises(sqlite3.OperationalError):
                connection.execute('DELETE FROM evidence')
            connection.close()
            self.assertEqual(path.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
