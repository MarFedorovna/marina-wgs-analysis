"""An offline HTML link still resolves after files receive category folders."""
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest

source = Path(__file__).resolve().parents[1] / 'transfer/write_local_index.py'
spec = importlib.util.spec_from_file_location('local_result_index', source)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class OfflineBundleTests(unittest.TestCase):
    def test_report_link_to_categorized_vcf_preserves_one_inode(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            objects = root / 'private-storage'; objects.mkdir()
            vcf = objects / 'public-fixture'; vcf.write_bytes(b'public fixture')
            entry = {'RelativePath': 'outputs/cohort-sample-views-v1/vcf/sample_0001.annotated.vcf.gz',
                     'Sha256': 'a' * 64, 'SizeBytes': vcf.stat().st_size}
            record = m.publish_size_only(root, entry, vcf, 'synthetic_fixture')
            canonical = root / record['LocalPath']
            report = root / 'Reports/outputs/cohort-sample-views-v1/index.html'
            report.write_text('<a href="vcf/sample_0001.annotated.vcf.gz">VCF</a>')
            relative_link = report.parent / 'vcf/sample_0001.annotated.vcf.gz'
            self.assertTrue(os.path.samefile(relative_link, canonical))
            self.assertTrue(os.path.samefile(relative_link, vcf))
            self.assertEqual(record['Verification'], 'size_only')
            self.assertIsNone(record['ActualSha256'])
            table = 'outputs/cohort-sample-views-v1/tables/cohort.functional.tsv.gz'
            self.assertEqual(m.category_path(table)[0], 'Reports')
            with self.assertRaises(ValueError):
                m.publish_report_alias(root, {'RelativePath': 'outputs/cohort-sample-views-v1/../../escape'}, vcf)


if __name__ == '__main__':
    unittest.main()
