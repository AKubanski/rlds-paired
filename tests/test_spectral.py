"""Tests for the spectral-gap validation (Paper V, Appendix C).

These assert the *honest* computed behaviour of the corrected coupled-ODE
background. At the canonical benchmark the Fenichel persistence bound (0.10) is
VIOLATED (delta_max ~ 0.245, peak in the transition window near z ~ 0.42); the
slow-manifold reduction is therefore only a marginal (~25%) correction there.
No headline number is read from a stored constant — everything is recomputed.
"""

import pytest
from rlds_paired.constants import FENICHEL_BOUND


class TestSpectralGapCanonical:
    def test_fenichel_violated(self, canonical_sg):
        # Expected: the bound does NOT hold at the canonical benchmark.
        assert not canonical_sg.fenichel_satisfied
        assert canonical_sg.delta_max > FENICHEL_BOUND

    def test_delta_max_value(self, canonical_sg):
        # Computed value ~0.245 (recomputed, not compared to a stored literal).
        assert 0.23 < canonical_sg.delta_max < 0.26

    def test_z_at_max_in_transition(self, canonical_sg):
        # Peak sits in the transition window, NOT pinned at the z~3 edge.
        # This is the regression guard against the old clipping bug.
        assert 0.35 < canonical_sg.z_at_max < 0.50


class TestSpectralGapPhysics:
    def test_decreases_with_K(self):
        # delta_max is weakly decreasing in K; both still exceed the bound.
        from rlds_paired.core.spectral import compute_spectral_gap

        sg_low = compute_spectral_gap(K=0.0106)  # canonical DR2 K
        sg_high = compute_spectral_gap(K=0.100)
        assert sg_high.delta_max < sg_low.delta_max
        assert sg_low.delta_max > FENICHEL_BOUND
        assert sg_high.delta_max > FENICHEL_BOUND

    def test_fenichel_violated_across_physical_K(self):
        # No physical coupling satisfies Fenichel; all delta_max land in band.
        from rlds_paired.core.spectral import spectral_gap_table

        rows = spectral_gap_table([0.005, 0.0106, 0.050, 0.100, 0.200])
        assert not any(r.fenichel_satisfied for r in rows)
        assert all(r.delta_max > FENICHEL_BOUND for r in rows)
        assert all(0.19 < r.delta_max < 0.26 for r in rows)

    def test_summary_contains_key_info(self, canonical_sg):
        s = canonical_sg.summary()
        assert "delta_max" in s
        assert "Fenichel" in s
