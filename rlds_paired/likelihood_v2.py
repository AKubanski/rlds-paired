"""
rlds_paired.likelihood_v2
=========================
Corrected likelihood pipeline (Stage 1 of the 2026-11 technical specification).

Differences from the original pipeline in :mod:`rlds_paired.mcmc`
-----------------------------------------------------------------
1.1  Supernovae: full STAT+SYS Pantheon+ covariance (default), analytic
     marginalisation over the absolute magnitude; d_L = (1+z_hel) D_M(z_HD).
1.2  omega_b is a free parameter (no longer tied to H0).
1.3  Radiation: Omega_r = omega_gamma (1 + 0.2271 N_eff) / h^2  (depends on h).
1.4  Sound horizon r_s(z) by direct integration of c_s / H.
1.5  z_* from the Hu & Sugiyama fitting formula.
1.6  CMB: Planck 2018 distance priors (R, l_A, omega_b) with correlations,
     Chen, Huang & Wang 2019 (arXiv:1808.05724), flat LCDM table.
1.7  r_d for BAO: Aubourg et al. 2015, eq. 16 (0.021 % accuracy for standard
     early-universe physics).

Early-universe matter density
-----------------------------
In models where matter does not scale as a^-3 at late times (RLDS), the
quantities that set the CMB and the sound horizon are the *early-time* matter
density.  We therefore use ``omega_m_early`` = rho_m(a) a^3 h^2 evaluated deep
in matter/radiation domination (a = 1e-5) in z_*, R, and r_d.  For LCDM, wCDM
and CPL this equals Omega_m h^2 exactly.

Neutrinos: one massive species, sum m_nu = 0.06 eV, treated as matter at all
times relevant here (omega_nu = 0.0107 * 0.06).  It enters only through
omega_cb = omega_m - omega_nu in the r_d formula.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp

from .data import get_desi_dr2_bao
from . import mcmc as _old

C_LIGHT = 299792.458  # km/s
T_CMB = 2.7255
OMEGA_GAMMA = 2.4728e-5 * (T_CMB / 2.7255) ** 4  # omega_gamma = Omega_gamma h^2
N_EFF = 3.046
OMEGA_R_H2 = OMEGA_GAMMA * (1.0 + 0.22710731766 * N_EFF)
SUM_MNU = 0.06
OMEGA_NU = 0.0107 * SUM_MNU

A_MIN = 1e-5                 # earliest scale factor of the tabulated background
N_GRID = 6000                # points in ln a for distance integrals

# ---------------------------------------------------------------------------
# CMB distance priors (Chen, Huang & Wang 2019, Table I, flat LCDM,
# Planck 2018 TT,TE,EE+lowE).  Order: R, l_A, omega_b.
# ---------------------------------------------------------------------------
CMB_MEAN = np.array([1.7502, 301.471, 0.02236])
CMB_SIG = np.array([0.0046, 0.0895, 0.00015])
CMB_CORR = np.array([[1.00, 0.46, -0.66],
                     [0.46, 1.00, -0.33],
                     [-0.66, -0.33, 1.00]])
CMB_COV = CMB_CORR * np.outer(CMB_SIG, CMB_SIG)
CMB_COV_INV = np.linalg.inv(CMB_COV)

# l_A convention calibration.
# Chen+2019 compute z_* and r_s(z_*) with CosmoMC conventions.  Our integral
# r_s (validated against CLASS to < 0.1 %) and the Hu-Sugiyama z_* reproduce
# their R exactly at the Planck 2018 TT,TE,EE+lowE mean point
# (omega_b=0.02236, omega_c=0.1202, H0=67.27): R = 1.75024 vs 1.7502,
# but give l_A = 302.314 vs 301.471, a constant convention offset of 0.28 %.
# Validated against CLASS (validation/lA_convention.py): the offset varies by
# < 7e-5 (< 0.2 sigma in l_A) across +-2 sigma shifts in omega_c, omega_b, H0.
# We therefore rescale l_A by the factor below, fixed at the Planck mean point.
LA_CALIBRATION = 0.9972124312921956


# ===========================================================================
# Physics helpers
# ===========================================================================
def z_star_hu_sugiyama(omega_b, omega_m):
    """Hu & Sugiyama (1996) decoupling redshift (coefficient 0.0783)."""
    g1 = 0.0783 * omega_b ** -0.238 / (1.0 + 39.5 * omega_b ** 0.763)
    g2 = 0.560 / (1.0 + 21.1 * omega_b ** 1.81)
    return 1048.0 * (1.0 + 0.00124 * omega_b ** -0.738) * (1.0 + g1 * omega_m ** g2)


def r_drag_aubourg(omega_cb, omega_b, omega_nu=OMEGA_NU):
    """Aubourg et al. (2015) eq. 16, Mpc."""
    return (55.154 * np.exp(-72.3 * (omega_nu + 0.0006) ** 2)
            / (omega_cb ** 0.25351 * omega_b ** 0.12807))


# ===========================================================================
# Background: tabulated rho_nr(a) in units of today's critical density
# ===========================================================================
@dataclass
class Background:
    """Flat background defined by the non-relativistic density rho_nr(ln a).

    rho_nr includes matter + dark energy (everything except radiation),
    in units of rho_crit,0.  ``rho_m`` is the matter part.
    """
    H0: float
    omega_b: float
    lna: np.ndarray        # increasing, from ln(A_MIN) to 0
    rho_nr: np.ndarray
    rho_m: np.ndarray
    label: str = ""

    def __post_init__(self):
        self.h = self.H0 / 100.0
        self.Omega_r = OMEGA_R_H2 / self.h ** 2
        a = np.exp(self.lna)
        self.E_grid = np.sqrt(self.rho_nr + self.Omega_r * a ** -4)
        # comoving distance chi(a) = c/H0 * int_a^1 dlna' / (a' E)
        integrand = 1.0 / (a * self.E_grid)
        cum = cumulative_trapezoid(integrand, self.lna, initial=0.0)
        self.chi_grid = (C_LIGHT / self.H0) * (cum[-1] - cum)   # chi(a)
        self.omega_m_early = float(self.rho_m[0] * a[0] ** 3 * self.h ** 2)
        self.Omega_m_today = float(self.rho_m[-1])
        # sound horizon r_s(a) = int_0^a c_s da' / (a'^2 H)
        Rb = 3.0 * self.omega_b / (4.0 * OMEGA_GAMMA) * a
        cs = C_LIGHT / np.sqrt(3.0 * (1.0 + Rb))
        rs_int = cs / (a * self.E_grid * self.H0)               # d r_s / d ln a
        # analytic piece for a < A_MIN (radiation domination, R_b -> 0)
        rs0 = C_LIGHT / np.sqrt(3.0) * a[0] / (self.H0 * np.sqrt(self.Omega_r))
        self.rs_grid = rs0 + cumulative_trapezoid(rs_int, self.lna, initial=0.0)

    # -- interpolation in ln a ---------------------------------------------
    def _interp(self, arr, z):
        return np.interp(-np.log1p(np.asarray(z, float)), self.lna, arr)

    def E(self, z):
        return self._interp(self.E_grid, z)

    def DM(self, z):
        return self._interp(self.chi_grid, z)

    def DH(self, z):
        return C_LIGHT / (self.H0 * self.E(z))

    def DV(self, z):
        z = np.asarray(z, float)
        return (z * self.DM(z) ** 2 * self.DH(z)) ** (1.0 / 3.0)

    def r_s(self, z):
        return self._interp(self.rs_grid, z)

    # -- derived CMB / BAO quantities --------------------------------------
    @property
    def z_star(self):
        return z_star_hu_sugiyama(self.omega_b, self.omega_m_early)

    @property
    def r_d(self):
        return r_drag_aubourg(self.omega_m_early - OMEGA_NU, self.omega_b)

    def cmb_vector(self):
        zs = self.z_star
        dm = self.DM(zs)
        R = np.sqrt(self.omega_m_early) * 100.0 * dm / C_LIGHT
        lA = LA_CALIBRATION * np.pi * dm / self.r_s(zs)
        return np.array([R, lA, self.omega_b])


def _grid():
    return np.linspace(np.log(A_MIN), 0.0, N_GRID)


def lcdm(Om, H0, omega_b):
    lna = _grid(); a = np.exp(lna)
    Or = OMEGA_R_H2 / (H0 / 100) ** 2
    rm = Om * a ** -3
    return Background(H0, omega_b, lna, rm + (1 - Om - Or), rm, "LCDM")


def wcdm(Om, H0, omega_b, w):
    lna = _grid(); a = np.exp(lna)
    Or = OMEGA_R_H2 / (H0 / 100) ** 2
    rm = Om * a ** -3
    return Background(H0, omega_b, lna, rm + (1 - Om - Or) * a ** (-3 * (1 + w)), rm, "wCDM")


def cpl(Om, H0, omega_b, w0, wa):
    lna = _grid(); a = np.exp(lna)
    Or = OMEGA_R_H2 / (H0 / 100) ** 2
    rm = Om * a ** -3
    de = (1 - Om - Or) * a ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
    return Background(H0, omega_b, lna, rm + de, rm, "CPL")


def rlds_free(Om, H0, omega_b, K, R_c, Delta_R, alpha):
    """RLDS ODE with Omega_m free: R(0) = Om/(1-Omega_r), integrated backward."""
    Or = OMEGA_R_H2 / (H0 / 100) ** 2
    lna = _grid()

    def ode(N, y):
        R = min(max(y[0], 1e-10), 1 - 1e-10)
        return [_old.F_rhs(R, alpha, K, R_c, Delta_R), -3.0 + alpha * (1.0 - R)]

    sol = solve_ivp(ode, [0.0, lna[0]], [Om / (1 - Or), 0.0], method="RK45",
                    rtol=1e-9, atol=1e-12, dense_output=True, max_step=0.05)
    if not sol.success:
        return None
    R, u = sol.sol(lna)
    R = np.clip(R, 1e-10, 1 - 1e-10)
    rho_nr = (1 - Or) * np.exp(u)
    return Background(H0, omega_b, lna, rho_nr, R * rho_nr, "RLDS-free")


def rlds_shooting(K, H0, omega_b, R_c, Delta_R, alpha):
    """Original self-consistent RLDS: Omega_m fixed by the shooting condition
    of :mod:`rlds_paired.mcmc` (start from LCDM at z=1100, require R(0)=Omega_m).
    The trajectory is unique, so it is rebuilt here by backward integration."""
    m = _old.ShootingRLDS(alpha, K, H0, R_c, Delta_R)
    if not m._valid:
        return None
    return rlds_free(m.Omega_m, H0, omega_b, K, R_c, Delta_R, alpha)


# ===========================================================================
# Data
# ===========================================================================
_SN = {}


def load_sn(data_dir=None, full_cov=True):
    """Pantheon+SH0ES: same selection as the original loader (1590 SNe)."""
    import os
    d = Path(data_dir or os.environ.get("RLDS_DATA_DIR") or Path.home() / "cosmology_data")
    key = (str(d), full_cov)
    if key in _SN:
        return _SN[key]
    zhd, zhel, mu, err, keep = [], [], [], [], []
    with open(d / "Pantheon+SH0ES.dat") as fh:
        fh.readline()
        for line in fh:
            p = line.split()
            z, zh, m, e = float(p[2]), float(p[6]), float(p[10]), float(p[11])
            ok = 0.01 < z < 2.5 and m > 0 and 0 < e < 10
            keep.append(ok)
            if ok:
                zhd.append(z); zhel.append(zh); mu.append(m); err.append(e)
    keep = np.array(keep)
    if full_cov:
        covf = d / "Pantheon+SH0ES_STAT+SYS.cov"
        n = int(open(covf).readline())
        C = np.loadtxt(covf, skiprows=1).reshape(n, n)[np.ix_(keep, keep)]
    else:
        C = np.diag(np.array(err) ** 2)
    Ci = np.linalg.inv(C)
    out = dict(zhd=np.array(zhd), zhel=np.array(zhel), mu=np.array(mu), Ci=Ci,
               Ci_sum=Ci.sum(), Ci_rowsum=Ci.sum(axis=1))
    _SN[key] = out
    return out


# ===========================================================================
# Likelihoods (return chi^2)
# ===========================================================================
def chi2_sn(bg, sn):
    dl = (1.0 + sn["zhel"]) * bg.DM(sn["zhd"])
    if np.any(dl <= 0) or not np.all(np.isfinite(dl)):
        return np.inf
    delta = sn["mu"] - (5.0 * np.log10(dl) + 25.0)
    A = delta @ sn["Ci"] @ delta
    B = delta @ sn["Ci_rowsum"]
    return float(A - B * B / sn["Ci_sum"])


def chi2_bao(bg):
    entries, data, _cov, cov_inv = get_desi_dr2_bao()
    rd = bg.r_d
    fn = {"DM": bg.DM, "DH": bg.DH, "DV": bg.DV}
    th = np.array([fn[q](z) / rd for z, q in entries])
    d = data - th
    return float(d @ cov_inv @ d)


def chi2_cmb(bg):
    d = bg.cmb_vector() - CMB_MEAN
    return float(d @ CMB_COV_INV @ d)


def is_physical(bg):
    """Reject backgrounds with non-finite H(z) or (almost) no early matter
    (strong RLDS coupling can drain it; omega_cb must stay positive)."""
    return (bg is not None and np.all(np.isfinite(bg.E_grid))
            and np.all(bg.rho_nr > 0) and bg.omega_m_early > OMEGA_NU + 1e-4)


def chi2_total(bg, sn=None, use=("sn", "bao", "cmb")):
    if not is_physical(bg):
        return np.inf
    c = 0.0
    if "cmb" in use:
        c += chi2_cmb(bg)
    if "bao" in use:
        c += chi2_bao(bg)
    if "sn" in use:
        c += chi2_sn(bg, sn if sn is not None else load_sn())
    return c if np.isfinite(c) else np.inf


__all__ = ["Background", "lcdm", "wcdm", "cpl", "rlds_free", "rlds_shooting",
           "load_sn", "chi2_sn", "chi2_bao", "chi2_cmb", "chi2_total",
           "z_star_hu_sugiyama", "r_drag_aubourg", "CMB_MEAN", "CMB_COV"]
