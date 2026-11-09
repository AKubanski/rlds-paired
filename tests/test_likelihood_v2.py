"""Tests for the corrected likelihood pipeline (Stage 1, tasks 1.1-1.8)."""
import os
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from rlds_paired import likelihood_v2 as L
from rlds_paired import mcmc as M

# Planck 2018 TT,TE,EE+lowE mean point (defines the Chen+2019 prior means)
WB, WC, H0 = 0.02236, 0.1202, 67.27
OM = (WB + WC + L.OMEGA_NU) / (H0 / 100) ** 2

_DATA = Path(os.environ.get("RLDS_DATA_DIR", Path.home() / "cosmology_data"))
needs_sn = pytest.mark.skipif(not (_DATA / "Pantheon+SH0ES.dat").exists(),
                              reason="Pantheon+ data not available")
needs_cov = pytest.mark.skipif(not (_DATA / "Pantheon+SH0ES_STAT+SYS.cov").exists(),
                               reason="Pantheon+ covariance not available")


def test_radiation_depends_on_h():                                    # 1.3
    a, b = L.lcdm(0.3, 60, WB), L.lcdm(0.3, 75, WB)
    assert a.Omega_r / b.Omega_r == pytest.approx((75 / 60) ** 2)


def test_omega_b_is_independent_parameter():                          # 1.2
    bg = L.lcdm(0.3, 67, 0.0210)
    assert bg.cmb_vector()[2] == 0.0210


def test_comoving_distance_matches_quadrature():                     # grid accuracy
    bg = L.lcdm(OM, H0, WB)
    Or = bg.Omega_r
    E = lambda z: np.sqrt(OM * (1 + z) ** 3 + Or * (1 + z) ** 4 + 1 - OM - Or)
    for z in (0.3, 1.0, 2.33, 1090.0):
        ref = L.C_LIGHT / H0 * quad(lambda x: 1 / E(x), 0, z, limit=500)[0]
        assert bg.DM(z) / ref - 1 == pytest.approx(0, abs=2e-5)


def test_sound_horizon_against_class_reference():                     # 1.4
    # CLASS (v3, same cosmology as test_vs_class) comoving sound horizon at z=1091.91
    bg = L.lcdm(0.31519, 67.36, 0.02237)
    assert bg.r_s(1091.91) == pytest.approx(144.26, rel=1.5e-3)


def test_z_star_formula():                                            # 1.5
    assert L.z_star_hu_sugiyama(0.02237, 0.1430) == pytest.approx(1091.9, abs=0.2)


def test_cmb_vector_reproduces_prior_means_at_planck_point():         # 1.6
    R, lA, wb = L.lcdm(OM, H0, WB).cmb_vector()
    assert R == pytest.approx(L.CMB_MEAN[0], abs=2e-4)      # 0.04 sigma
    assert lA == pytest.approx(L.CMB_MEAN[1], abs=1e-6)     # by calibration
    assert L.chi2_cmb(L.lcdm(OM, H0, WB)) < 0.1


def test_r_drag_formula():                                            # 1.7
    # CLASS rs_d = 147.097 Mpc for (0.02237, 0.1200, mnu=0.06)
    assert L.r_drag_aubourg(0.1200 + 0.02237, 0.02237) == pytest.approx(147.097, rel=5e-4)


def test_rlds_free_reduces_to_lcdm():
    a = L.lcdm(0.31, 67.5, WB)
    b = L.rlds_free(0.31, 67.5, WB, 0.0, 0.5, 0.1, 3.0)
    for z in (0.5, 2.0, 1090.0):
        assert b.DM(z) / a.DM(z) - 1 == pytest.approx(0, abs=1e-6)
    assert b.omega_m_early == pytest.approx(a.omega_m_early, rel=1e-6)


def test_rlds_free_alpha_is_constant_w_when_K_zero():
    a = L.wcdm(0.31, 67.5, WB, -2.85 / 3)
    b = L.rlds_free(0.31, 67.5, WB, 0.0, 0.5, 0.1, 2.85)
    assert b.DM(1.0) / a.DM(1.0) - 1 == pytest.approx(0, abs=1e-6)


@pytest.mark.parametrize("K", [0.001, 0.01, 0.05])
def test_comparison_principle_no_shooting_root_at_alpha_3(K):
    """At alpha=3 and K>0 the shooting residual R(0)-Om is positive for every Om."""
    res = [M._integrate_Ru(3.0, K, Om, 0.565, 0.17)[0] - Om
           for Om in np.linspace(0.15, 0.5, 8)]
    assert all(r > 0 for r in res)


@needs_sn
def test_sn_marginalisation_invariant_to_offset():                    # 1.1
    sn = dict(L.load_sn(full_cov=False))
    bg = L.lcdm(0.31, 67.5, WB)
    c0 = L.chi2_sn(bg, sn)
    sn["mu"] = sn["mu"] + 0.37
    assert L.chi2_sn(bg, sn) == pytest.approx(c0, rel=1e-9)


@needs_sn
@needs_cov
def test_full_covariance_loaded():                                    # 1.1
    sn = L.load_sn(full_cov=True)
    assert sn["Ci"].shape == (1590, 1590)


@needs_sn
def test_old_pipeline_regression():                                   # 1.8
    assert M.reproduce_reference_chi2(fast_distance=True) == pytest.approx(
        M.REFERENCE_CHI2, abs=0.01)
