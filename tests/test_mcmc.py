"""Tests for the shooting-MCMC module.

The pure-logic tests (parameters, priors) always run. The reproduction tests
require the external Pantheon+ data file and are skipped when it is absent
(e.g. on CI without the data record); set $RLDS_DATA_DIR or place
Pantheon+SH0ES.dat in ~/cosmology_data to enable them.
"""

import numpy as np
import pytest

from rlds_paired import mcmc


def _have_pantheon() -> bool:
    try:
        from rlds_paired.data import load_pantheon_plus

        load_pantheon_plus()
        return True
    except Exception:
        return False


HAVE_DATA = _have_pantheon()
needs_data = pytest.mark.skipif(
    not HAVE_DATA, reason="Pantheon+ data file not available"
)


# --- pure logic (no data) --------------------------------------------------
def test_param_layout():
    assert mcmc.NDIM == 5
    assert mcmc.PARAM_NAMES == ["K", "H0", "R_c", "Delta_R", "alpha"]


def test_prior_accepts_canonical():
    # Arbitrary in-bounds probe vector [K, H0, R_c, Delta_R, alpha] used only to
    # exercise the prior's accept path; NOT the DR2 canonical parameter set
    # (see rlds_paired.constants / CANONICAL for those).
    assert mcmc.log_prior([0.012, 67.44, 0.51, 0.15, 3.02]) == 0.0


@pytest.mark.parametrize(
    "theta",
    # Each vector pushes exactly ONE parameter out of bounds; the remaining
    # entries are arbitrary in-bounds fixtures (non-canonical), not headline values.
    [
        [0.012, 67.44, 0.51, 0.15, 2.99],   # alpha below the 3.01 floor
        [0.30, 67.44, 0.51, 0.15, 3.02],    # K above upper bound
        [0.012, 80.0, 0.51, 0.15, 3.02],    # H0 above upper bound
        [0.012, 67.44, 0.51, 0.25, 3.02],   # Delta_R above upper bound
    ],
)
def test_prior_rejects_out_of_bounds(theta):
    assert mcmc.log_prior(theta) == -np.inf


def test_alpha_floor_is_phantom_only():
    # the prior forbids alpha < 3.01, i.e. w_DE > -1.003: phantom-only by prior
    lo = mcmc.PRIOR_BOUNDS["alpha"][0]
    assert lo == pytest.approx(3.01)
    assert -lo / 3 < -1.0  # w_DE strictly below -1 at the floor


# --- reproduction (needs data) --------------------------------------------
@needs_data
def test_reference_chi2_reproduces():
    chi2 = mcmc.reproduce_reference_chi2(fast_distance=False)
    assert chi2 == pytest.approx(mcmc.REFERENCE_CHI2, abs=0.01)


@needs_data
def test_fast_distance_matches_slow():
    slow = mcmc.reproduce_reference_chi2(fast_distance=False)
    fast = mcmc.reproduce_reference_chi2(fast_distance=True)
    assert abs(fast - slow) < 0.05


@needs_data
def test_best_fit_model_self_consistent():
    m = mcmc.get_selfconsistent_model(*mcmc.REFERENCE_BEST_FIT)
    assert m is not None and m._valid
    # shooting guarantees R(z=0) = Omega_m
    assert abs(m.R_today - m.Omega_m) < 1e-3
    # derived present-day matter fraction in the expected band
    assert 0.28 < m.Omega_m < 0.31
