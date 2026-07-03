"""Regression tests that lock the canonical headline results and the public
array-input contract of RLDSTrajectory.

These guard the *scientific* numbers (paired coefficient, spectral gap,
attractor, present-day matter fraction) against silent regressions, and the
R/E/E2 array support against re-breaking.
"""

import numpy as np

import rlds_paired as rp
from rlds_paired.core.ode import RLDSTrajectory
from rlds_paired.reference import reference_results


def test_reference_headline_values():
    """Canonical results recomputed from the model (Papers III-V)."""
    ref = reference_results()
    # paired coefficient: converged value ~ -0.193 (DR2), growth suppression (negative)
    assert ref.paired_coefficient < 0
    assert -0.21 < ref.paired_coefficient < -0.17
    assert ref.paired_coefficient_err < 0.02
    assert ref.paired_r2 > 0.95
    # spectral gap ~ 0.244 (DR2), peak in the transition window, Fenichel (0.10) violated
    assert 0.23 < ref.delta_max < 0.26
    assert 0.2 < ref.delta_max_z < 0.7
    assert ref.fenichel_satisfied is False
    assert ref.delta_max > ref.fenichel_bound


def test_trajectory_canonical_points():
    """Background sanity: normalisation, attractor, present-day fraction."""
    t = RLDSTrajectory()
    assert abs(t.E(0.0) - 1.0) < 1e-6  # H(z=0)/H_0 = 1
    assert abs(t.R_star - rp.R_STAR) < 5e-3  # attractor ~ sqrt(K)
    assert 0.28 < t.R_today < 0.31  # present-day matter fraction
    assert t.R_today > t.R_star  # above the attractor today (recent transition)


def test_R_E_accept_arrays():
    """Public contract: R/E/E2 accept ndarray (annotation says float | ndarray)."""
    t = RLDSTrajectory()
    z = np.array([0.0, 0.5, 1.0, 2.0])
    for fn in (t.R, t.E2, t.E):
        out = fn(z)
        assert isinstance(out, np.ndarray)
        assert out.shape == z.shape
        assert np.all(np.isfinite(out))


def test_R_E_scalar_returns_float():
    """Scalar input must still return a plain Python float."""
    t = RLDSTrajectory()
    assert isinstance(t.R(0.0), float)
    assert isinstance(t.E2(0.5), float)
    assert isinstance(t.E(1.0), float)
