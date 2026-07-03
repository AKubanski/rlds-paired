"""Tests for the linear growth factor (core/growth.py)."""

import numpy as np
import pytest
from rlds_paired.constants import K_BEST, OMEGA_M


class TestGrowthFactorLCDM:
    def test_range_z05(self):
        from rlds_paired.core.growth import growth_factor

        fs8 = growth_factor(OMEGA_M, K=0.0, z_eval=0.5)
        assert 0.35 < fs8 < 0.55

    def test_decreases_with_z(self):
        from rlds_paired.core.growth import growth_factor

        fs8_low = growth_factor(OMEGA_M, K=0.0, z_eval=0.3)
        fs8_high = growth_factor(OMEGA_M, K=0.0, z_eval=1.0)
        assert fs8_low > fs8_high


class TestLinearScaling:
    def test_r2_above_threshold(self):
        from rlds_paired.core.growth import growth_factor
        from scipy.stats import linregress

        K_vals = np.linspace(0.02, 0.20, 10)
        fs8_0 = growth_factor(OMEGA_M, K=0.0, z_eval=0.5)
        deltas = [
            (growth_factor(OMEGA_M, K=K, z_eval=0.5) - fs8_0) / fs8_0 for K in K_vals
        ]
        _, _, r, _, _ = linregress(K_vals, deltas)
        assert r**2 > 0.95

    def test_canonical_coefficient_computed(self):
        # Canonical protocol: full coupled-ODE k-scan must give a NEGATIVE
        # coefficient in the documented band. Recomputed, never read from a
        # stored literal. n_points=12 keeps the test fast but stable.
        from rlds_paired.analysis import k_scan

        res = k_scan(n_points=12, verbose=False)
        assert res.coefficient < 0
        assert -0.21 < res.coefficient < -0.17
        assert res.r2 > 0.95
