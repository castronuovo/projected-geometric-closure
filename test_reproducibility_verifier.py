"""Regression checks for complete, partial and mismatched output sets."""
import contextlib
import csv
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import verify_reproducibility as verifier


class OutputCoverageTests(unittest.TestCase):
    def invoke(self, omitted=()):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory)
            source = Path(verifier.__file__).parent
            for name in verifier.FILES:
                if name not in omitted:
                    (reference / name).write_bytes((source / name).read_bytes())
            output = io.StringIO()
            with patch('sys.argv', ['verify_reproducibility.py', directory]):
                with contextlib.redirect_stdout(output):
                    verifier.main()
            return output.getvalue()

    def test_complete_set_includes_both_new_outputs(self):
        self.assertEqual(len(verifier.FILES), 19)
        self.assertIn('All 19 scientific outputs', self.invoke())

    def test_previous_complete_set_is_reported_as_partial(self):
        message = self.invoke(verifier.V16_FILES_NOT_IN_V15_CI)
        self.assertIn('17 of 19', message)
        self.assertIn('partial coverage', message)
        self.assertNotIn('All 19', message)

    def test_legacy_set_is_reported_as_partial(self):
        missing = verifier.V16_FILES_NOT_IN_V15_CI | verifier.V15_FILES_NOT_IN_LEGACY_CI
        self.assertIn('14 of 19', self.invoke(missing))

    def test_one_missing_new_output_is_rejected(self):
        with self.assertRaises(FileNotFoundError):
            self.invoke({'causal_second_order_certified_range.csv'})

    def test_missing_inherited_output_is_rejected(self):
        with self.assertRaises(FileNotFoundError):
            self.invoke({'causal_growth_transport_summary.json'})

    def test_changed_margin_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(verifier.__file__).parent / 'causal_second_order_certified_range.csv'
            modified = Path(directory) / source.name
            with source.open(newline='') as stream:
                rows = list(csv.reader(stream))
            column = rows[0].index('margin_new')
            rows[1][column] = str(float(rows[1][column]) + 0.01)
            with modified.open('w', newline='') as stream:
                csv.writer(stream).writerows(rows)
            with self.assertRaises(AssertionError):
                verifier.compare_csv(source, modified, optimization=True)

    def test_changed_structural_metadata_is_rejected(self):
        with self.assertRaises(AssertionError):
            verifier.compare_json({'n_masses_pair_grid': 243},
                                  {'n_masses_pair_grid': 242}, optimization=True)

    def test_changed_verification_boolean_is_rejected(self):
        with self.assertRaises(AssertionError):
            verifier.compare_json({'new_allowance_bounds_all_stored_stationary_errors': True},
                                  {'new_allowance_bounds_all_stored_stationary_errors': False},
                                  optimization=True)


if __name__ == '__main__':
    unittest.main()
