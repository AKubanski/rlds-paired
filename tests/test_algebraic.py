"""Tests for the algebraic IDE model (core/algebraic.py)."""

import numpy as np
import pytest
from rlds_paired.constants import DR, K_BEST, OMEGA_M, RC


class TestAlgebraicLCDMLimit:
    """K=0 must reproduce exact ΛCDM."""

    def test_F_unity_at_K0(self):
        from rlds_paired.core.algebraic import F_transition

        assert F_transition(0.5, K=0.0) == 1.0

    def test_E2_normalised_at_K0(self):
        from rlds_paired.core.algebraic import E_squared_alg

        assert abs(E_squared_alg(0.0, OMEGA_M, K=0.0) - 1.0) < 1e-10

    def test_solveR_lcdm_initial_guess(self):
        from rlds_paired.core.algebraic import solve_R

        # At K=0 the fixed point should equal the ΛCDM matter fraction
        z = 1.0
        R = solve_R(z, OMEGA_M, K=0.0, Rc=RC, dR=DR)
        R_lcdm = OMEGA_M * (1 + z) ** 3 / (OMEGA_M * (1 + z) ** 3 + (1 - OMEGA_M))
        assert abs(R - R_lcdm) < 1e-10


class TestAlgebraicMonotonicity:
    """E²(z) must increase monotonically with z."""

    def test_E2_monotone_K0(self):
        from rlds_paired.core.algebraic import E_squared_alg

        z = np.linspace(0.0, 3.0, 30)
        E2 = [E_squared_alg(zi, OMEGA_M, K=0.0) for zi in z]
        assert all(E2[i] <= E2[i + 1] for i in range(len(E2) - 1))

    def test_E2_monotone_Kbest(self):
        from rlds_paired.core.algebraic import E_squared_alg

        z = np.linspace(0.0, 3.0, 30)
        E2 = [E_squared_alg(zi, OMEGA_M, K=K_BEST) for zi in z]
        assert all(E2[i] <= E2[i + 1] for i in range(len(E2) - 1))


class TestPhantomExcursion:
    """K>0 must produce w_eff < -1 at late times."""

    def test_wmin_below_minus1(self):
        from rlds_paired.core.algebraic import w_eff_alg

        z_grid = np.linspace(0.01, 2.0, 100)
        w_min = min(w_eff_alg(z, OMEGA_M, K_BEST, RC, DR) for z in z_grid)
        assert w_min < -1.0

    def test_wmin_K0_is_minus1(self):
        from rlds_paired.core.algebraic import w_eff_alg

        w = w_eff_alg(0.5, OMEGA_M, K=0.0, Rc=RC, dR=DR)
        assert abs(w + 1.0) < 0.01
