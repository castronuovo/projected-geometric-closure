"""Regression checks for complete, partial and mismatched output sets."""
import contextlib
import csv
import io
import json
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

    def test_complete_set_includes_causal_amplitude_output(self):
        self.assertEqual(len(verifier.FILES), 32)
        self.assertIn('All 32 scientific outputs', self.invoke())

    def test_initial_v18_set_is_reported_as_partial(self):
        self.assertIn('27 of 32', self.invoke(verifier.V18_EXTRA_FILES_NOT_IN_INITIAL_CI))

    def test_missing_one_extra_output_is_rejected(self):
        with self.assertRaises(FileNotFoundError):
            self.invoke({'rigorous_lp_certificate_M800_m025.json'})

    def test_changed_extra_target_bound_is_rejected(self):
        source = Path(verifier.__file__).parent / 'rigorous_lp_certificate_M800_m025.json'
        reference = json.loads(source.read_text())
        modified = json.loads(source.read_text())
        modified['rows'][0]['estimated_lower_bound'] += 1e-6
        with self.assertRaises(AssertionError):
            verifier.compare_linear_program_output(reference, modified)

    def test_v17_set_is_reported_as_partial(self):
        self.assertIn('22 of 32', self.invoke(verifier.V18_FILES_NOT_IN_V17_CI))

    def test_v16_set_is_reported_as_partial(self):
        self.assertIn('19 of 32', self.invoke(
            verifier.V18_FILES_NOT_IN_V17_CI | verifier.V17_FILES_NOT_IN_V16_CI))

    def test_previous_complete_set_is_reported_as_partial(self):
        message = self.invoke(verifier.V18_FILES_NOT_IN_V17_CI
                              | verifier.V17_FILES_NOT_IN_V16_CI | verifier.V16_FILES_NOT_IN_V15_CI)
        self.assertIn('17 of 32', message)
        self.assertIn('partial coverage', message)
        self.assertNotIn('All 32', message)

    def test_legacy_set_is_reported_as_partial(self):
        missing = (verifier.V18_FILES_NOT_IN_V17_CI | verifier.V17_FILES_NOT_IN_V16_CI
                   | verifier.V16_FILES_NOT_IN_V15_CI | verifier.V15_FILES_NOT_IN_LEGACY_CI)
        self.assertIn('14 of 32', self.invoke(missing))

    def test_workflow_preserves_exact_covered_output_set(self):
        import re
        workflow = (Path(verifier.__file__).parent / '.github/workflows/reproducibility.yml').read_text()
        block = workflow.split('- name: Preserve frozen scientific outputs',1)[1].split('- name:',1)[0]
        preserved = set(re.findall(r'([A-Za-z0-9_]+\.(?:json|csv))',block))
        self.assertEqual(preserved,set(verifier.FILES))

    def test_changed_linear_program_bound_is_rejected(self):
        source = Path(verifier.__file__).parent / 'lp_certificate_results.json'
        reference = json.loads(source.read_text())
        modified = json.loads(source.read_text())
        modified['rows'][1]['LP_lower_bound_first_order'] *= 1.001
        with self.assertRaises(AssertionError):
            verifier.compare_linear_program_output(reference, modified)

    def test_linear_program_history_vector_is_not_compared(self):
        reference = {'rows': [{'LP': 1.0e-4, 'b_nodes': [0.02, 0.03]}]}
        modified = {'rows': [{'LP': 1.0e-4, 'b_nodes': [0.021, 0.029]}]}
        verifier.compare_linear_program_output(reference, modified)

    def test_changed_joint_history_prior_is_rejected(self):
        source = Path(verifier.__file__).parent / 'bounded_joint_history_results.json'
        reference = json.loads(source.read_text())
        modified = json.loads(source.read_text())
        modified['cases'][0]['H8']['conservative_log_slope_bound'] = 2.1
        with self.assertRaises(AssertionError):
            verifier.compare_bounded_joint_history(reference, modified)

    def test_changed_joint_history_snapshot_is_rejected(self):
        source = Path(verifier.__file__).parent / 'bounded_joint_history_results.json'
        reference = json.loads(source.read_text())
        modified = json.loads(source.read_text())
        modified['cases'][0]['H8']['prediction'][0] += 0.001
        with self.assertRaises(AssertionError):
            verifier.compare_bounded_joint_history(reference, modified)

    def joint_outputs(self):
        source = Path(verifier.__file__).parent / 'bounded_joint_history_results.json'
        return json.loads(source.read_text()), json.loads(source.read_text())

    def test_equivalent_joint_history_witness_is_accepted(self):
        import numpy as np
        from bounded_joint_history import restrictions
        reference, modified = self.joint_outputs()
        case = modified['cases'][0]
        fit = case['H8']
        fit['coefficients'][0] += 1e-10
        prediction, states = verifier._joint_history_prediction(tuple(fit['coefficients']))
        fit.update(restrictions(np.asarray(fit['coefficients'])))
        fit['prediction'], fit['states'] = prediction, states
        residual = [p-t for p, t in zip(prediction, case['target'])]
        fit['max_joint_residual'] = max(map(abs, residual))
        fit['max_power_residual'] = max(abs(residual[i]) for i in (0, 1, 4, 5))
        fit['max_growthrate_residual'] = max(abs(residual[i]) for i in (2, 3, 6, 7))
        verifier.compare_bounded_joint_history(reference, modified)

    def test_cross_platform_joint_endpoint_regression(self):
        fit = {'max_joint_residual': 1.6690502537830554e-6}
        tolerance = verifier._joint_fit_comparison_tolerance(fit, fit)
        verifier.compare_json(0.00394694610733332, 0.00394694603509595,
                              absolute_tolerance=tolerance)
        with self.assertRaises(AssertionError):
            verifier.compare_json(0.00394694610733332,
                                  0.00394694610733332 + 2*tolerance,
                                  absolute_tolerance=tolerance)

    def test_inconsistent_joint_coefficients_are_rejected(self):
        import numpy as np
        from bounded_joint_history import restrictions
        reference, modified = self.joint_outputs()
        fit = modified['cases'][0]['H8']
        fit['coefficients'][0] += 1e-5
        fit.update(restrictions(np.asarray(fit['coefficients'])))
        with self.assertRaisesRegex(AssertionError, 'forward_prediction'):
            verifier.compare_bounded_joint_history(reference, modified)

    def test_joint_residual_must_match_recorded_snapshots(self):
        reference, modified = self.joint_outputs()
        modified['cases'][0]['H8']['max_joint_residual'] += 1e-8
        with self.assertRaisesRegex(AssertionError, 'max_joint_residual'):
            verifier.compare_bounded_joint_history(reference, modified)

    def test_joint_tied_family_is_preserved(self):
        reference, modified = self.joint_outputs()
        modified['cases'][0]['H2']['coefficients'][0] += 1e-7
        with self.assertRaisesRegex(AssertionError, 'tied coefficients'):
            verifier.compare_bounded_joint_history(reference, modified)

    def test_nonfinite_joint_coefficient_is_rejected(self):
        reference, modified = self.joint_outputs()
        modified['cases'][0]['H8']['coefficients'][0] = float('nan')
        with self.assertRaisesRegex(AssertionError, 'coefficients invalid'):
            verifier.compare_bounded_joint_history(reference, modified)

    def test_joint_contract_retains_strict_tolerance(self):
        reference, modified = self.joint_outputs()
        modified['contract']['Omega_m0'] += 1e-8
        with self.assertRaisesRegex(AssertionError, 'contract.Omega_m0'):
            verifier.compare_bounded_joint_history(reference, modified)

    def test_changed_inverse_residue_is_rejected(self):
        source = Path(verifier.__file__).parent / 'gapped_residue_inverse_diagnostic_results.json'
        reference = json.loads(source.read_text())
        modified = json.loads(source.read_text())
        modified['rows'][0]['residue_upper'] += 0.01
        with self.assertRaises(AssertionError):
            verifier.compare_inverse_diagnostic(reference, modified)

    def test_changed_counterfit_amplitude_is_rejected(self):
        source = Path(verifier.__file__).parent / 'causal_amplitude_reachability_results.json'
        reference = json.loads(source.read_text())
        modified = json.loads(source.read_text())
        modified['exact_amplitude_bounds'][0] = -0.01
        with self.assertRaises(AssertionError):
            verifier.compare_causal_amplitude(reference, modified)

    def test_failed_counterfit_residual_is_rejected(self):
        source = Path(verifier.__file__).parent / 'causal_amplitude_reachability_results.json'
        reference = json.loads(source.read_text())
        modified = json.loads(source.read_text())
        modified['independent_Radau_max_power_residual'] = 1e-6
        with self.assertRaises(AssertionError):
            verifier.compare_causal_amplitude(reference, modified)

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
