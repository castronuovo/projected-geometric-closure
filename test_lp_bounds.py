"""Regressions for bounds and history contracts corrected during review."""
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import rigorous_lp_certificate as constant
import parent_law_comparison as parent

class BoundContractTests(unittest.TestCase):
    def test_dual_respects_known_box_extremum(self):
        # Average perturbation plus one: minimum 1-ETA at constant lower endpoint.
        for module in (constant,parent):
            J=np.ones((8,module.M+1))/(module.M+1)
            lower,db,primal,gap=module.lp_box(np.ones(8),J)
            self.assertLessEqual(lower,1-module.ETA+1e-13)
            self.assertAlmostEqual(primal,1-module.ETA,places=10)
            self.assertGreaterEqual(gap,-1e-13)

    def test_failed_lp_cannot_emit_bound(self):
        with patch.object(constant,'linprog',return_value=SimpleNamespace(success=False,message='failed')):
            with self.assertRaisesRegex(RuntimeError,'failed'):
                constant.lp_box(np.ones(8),np.zeros((8,constant.M+1)))

    def test_early_derivative_norm_not_replaced_by_final_value(self):
        ep=[(np.ones(2),np.zeros(2),np.ones(2)),
            (np.ones(2),np.zeros(2),np.zeros(2))]
        rem=constant.remainder(np.ones(2),np.zeros(2),np.ones(2),ep)
        self.assertGreater(rem[2],0)
        self.assertEqual(rem[6],0)

    def test_parent_zero_perturbation_is_exact_reference_history(self):
        expected=parent.obs(parent.bref,.15)
        actual=parent.nonlinear(np.zeros(parent.M+1),.15)
        np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-12)

if __name__=='__main__':unittest.main()
