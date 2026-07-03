"""
rlds_paired.mcmc
================
Faithful reproduction of the canonical *shooting* MCMC that underlies the
RLDS cosmology results (Papers III--V).

Sampled parameters (5):  [K, H0, R_c, Delta_R, alpha].
Omega_m is **not** sampled -- for each (K, alpha, R_c, Delta_R) it is found by
shooting so that R(z=0) = Omega_m *exactly*. The posterior median of the
derived Omega_m is ~0.297.

Data (via :mod:`rlds_paired.data`):
    Pantheon+SH0ES (diagonal likelihood, analytic M-marginalisation),
    DESI DR2 BAO, Planck 2018 distance prior.

Reproduction notes
------------------
* ``EPS`` is 1e-6 here, matching the code that produced the archived chain.
  The maintained analysis branch (``rlds_paired.analysis``) uses 1e-12; the
  two agree to within sampling error (R(z=0) differs at the sub-sigma level).
* The reference run fixed **no** random seed, so only *statistical*
  reproduction is possible from scratch; the archived ``.npy`` chain is the
  canonical artefact. ``run_mcmc(seed=...)`` is provided for new runs.
* The supernova likelihood defaults to the diagonal error (exact reproduction).
  ``use_full_cov=True`` switches to a full covariance if one is supplied.
* The slow distance path (``quad`` per SN) is the original; ``fast_distance``
  builds one comoving-distance interpolator per model (~100x faster) and
  agrees with the slow path to < 1e-4 mag.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, quad, cumulative_trapezoid
from scipy.optimize import brentq, minimize_scalar

from .data import (
    load_pantheon_plus,
    get_desi_dr2_bao,
    get_planck2018_distance_prior,
    PLANCK2018_ZSTAR,
)

# ---------------------------------------------------------------------------
# CONSTANTS (match the reference run)
# ---------------------------------------------------------------------------
C_LIGHT = 299792.458            # km/s
EPS = 1e-6                      # ODE regularisation used for the archived chain
OMEGA_R = 9.15e-5              # radiation density today
OMEGA_B = 0.02237 / (67.36 / 100) ** 2   # baryon density parameter
Z_START = 1100.0               # ODE integration start (deep matter era)

PARAM_NAMES = ["K", "H0", "R_c", "Delta_R", "alpha"]
PARAM_LABELS = [r"$K$", r"$H_0$", r"$R_c$", r"$\Delta R$", r"$\alpha$"]
NDIM = 5

PRIOR_BOUNDS = {
    "K":       (0.005, 0.25),
    "H0":      (60.0, 75.0),
    "R_c":     (0.30, 0.60),
    "Delta_R": (0.03, 0.20),
    "alpha":   (3.01, 3.30),   # w_DE in [-1.10, -1.003]: phantom region by prior
}

# theta=[0.05, 67.44, 0.45, 0.12, 3.04] is the reference starting point.
THETA_START = np.array([0.05, 67.44, 0.45, 0.12, 3.04])
START_SCATTER = np.array([0.01, 0.5, 0.03, 0.02, 0.02])

# Archived best-fit (argmax log-prob) of the canonical DESI DR2 32x3000
# shooting chain, and its chi2 (diagonal SN, eps=1e-6, full DR2 BAO covariance).
# Used by the smoke test / CLI to confirm faithful reproduction of the chain.
# Note: R_c and Delta_R sit at their upper prior edges in the argmax; the
# headline canonical set uses the posterior medians (see constants.py).
REFERENCE_BEST_FIT = np.array(
    [0.00963176, 67.78986273, 0.59981095, 0.1998306, 3.01000794]
)
REFERENCE_CHI2 = 750.950730  # chi2 at REFERENCE_BEST_FIT (DR2, full cov)


# ===========================================================================
# RLDS ODE COMPONENTS
# ===========================================================================
def Psi(R, eps: float = EPS):
    return (1.0 - R) / (R + eps)


def S_switch(R, R_c, Delta_R):
    arg = np.clip((R_c - R) / Delta_R, -50.0, 50.0)
    return 0.5 * (1.0 + np.tanh(arg))


def F_rhs(R, alpha, K, R_c, Delta_R):
    R = np.clip(R, 1e-10, 1.0 - 1e-10)
    B = K * alpha
    return -alpha * R * (1.0 - R) + B * Psi(R) * S_switch(R, R_c, Delta_R)


def H_aux(R, R_c, Delta_R, eps: float = EPS):
    psi = Psi(R, eps)
    s = S_switch(R, R_c, Delta_R)
    denom = psi * s
    if denom < 1e-30:
        return 1e30
    return R * (1.0 - R) / denom


def find_attractor(K, R_c, Delta_R):
    try:
        f = lambda R: H_aux(R, R_c, Delta_R) - K
        R_lo, R_hi = 1e-4, R_c - 1e-4
        if f(R_lo) * f(R_hi) > 0:
            R_hi = min(R_c + 0.1, 0.99)
            if f(R_lo) * f(R_hi) > 0:
                return np.nan
        return brentq(f, R_lo, R_hi, xtol=1e-12)
    except Exception:
        return np.nan


def R_LCDM(z, Omega_m):
    Om_z = Omega_m * (1.0 + z) ** 3
    return Om_z / (Om_z + (1.0 - Omega_m))


def _integrate_Ru(alpha, K, Omega_m, R_c, Delta_R, z_start=Z_START):
    """Integrate coupled [R, u] from z_start to z=0; return (R0, u0, sol)."""
    N_start = -np.log(1.0 + z_start)
    R_init = R_LCDM(z_start, Omega_m)

    def ode(N, y):
        R = np.clip(y[0], 1e-10, 1.0 - 1e-10)
        dR = F_rhs(R, alpha, K, R_c, Delta_R)
        du = -3.0 + alpha * (1.0 - R)
        return [dR, du]

    sol = solve_ivp(
        ode, [N_start, 0.0], [R_init, 0.0],
        method="RK45", rtol=1e-10, atol=1e-12,
        dense_output=True, max_step=0.05,
    )
    if not sol.success:
        return np.nan, np.nan, None
    R0 = float(np.clip(sol.sol(0.0)[0], 1e-10, 1 - 1e-10))
    u0 = float(sol.sol(0.0)[1])
    return R0, u0, sol


# ===========================================================================
# SELF-CONSISTENT MODEL (shooting on Omega_m)
# ===========================================================================
class ShootingRLDS:
    """Self-consistent RLDS background; Omega_m derived so R(z=0)=Omega_m."""

    def __init__(self, alpha, K, H0, R_c, Delta_R, *, _skip_shooting=False):
        self.alpha, self.K, self.H0 = alpha, K, H0
        self.R_c, self.Delta_R = R_c, Delta_R
        self.Omega_r = OMEGA_R
        self._dc_interp = None  # lazy fast-distance interpolator
        if _skip_shooting:
            return
        Om_sc, R0 = self._shoot()
        if Om_sc is None:
            self._mark_invalid()
            return
        self._finalise(Om_sc)

    # -- shooting -----------------------------------------------------------
    def _residual(self, Om):
        R0, _, _ = _integrate_Ru(self.alpha, self.K, Om, self.R_c, self.Delta_R)
        return 1.0 if np.isnan(R0) else R0 - Om

    def _shoot(self, Om_lo=0.15, Om_hi=0.50, tol=1e-5):
        r_lo, r_hi = self._residual(Om_lo), self._residual(Om_hi)
        if np.isnan(r_lo) or np.isnan(r_hi):
            return None, None
        if r_lo * r_hi > 0:
            try:
                res = minimize_scalar(lambda Om: abs(self._residual(Om)),
                                      bounds=(Om_lo, Om_hi), method="bounded")
                if abs(self._residual(res.x)) < tol:
                    return res.x, res.x - self._residual(res.x)
            except Exception:
                pass
            return None, None
        try:
            Om = brentq(self._residual, Om_lo, Om_hi, xtol=tol, maxiter=50)
            return Om, Om + self._residual(Om)
        except Exception:
            return None, None

    def _finalise(self, Om_sc):
        _, self.u_today, self.sol = _integrate_Ru(
            self.alpha, self.K, Om_sc, self.R_c, self.Delta_R)
        self.Omega_m = Om_sc
        self.R_today = float(np.clip(self.sol.sol(0.0)[0], 1e-10, 1 - 1e-10))
        self._rho_nr_ref = Om_sc * (1 + Z_START) ** 3 + (1 - Om_sc)
        self._E2_norm = self._rho_nr_ref * np.exp(self.u_today) + self.Omega_r
        self.R_star = find_attractor(self.K, self.R_c, self.Delta_R)
        self.self_consistency_error = abs(self.R_today - Om_sc)
        self._valid = True

    def _mark_invalid(self):
        self._valid = False
        self.Omega_m = self.R_today = self.u_today = self.R_star = np.nan
        self.self_consistency_error = np.inf
        self.sol = None
        self._rho_nr_ref = self._E2_norm = 1.0

    # -- background ---------------------------------------------------------
    def R(self, z):
        if not self._valid:
            return np.nan
        return float(np.clip(self.sol.sol(-np.log(1.0 + z))[0], 1e-10, 1 - 1e-10))

    def u(self, z):
        if not self._valid:
            return 0.0
        return float(self.sol.sol(-np.log(1.0 + z))[1])

    def E_squared(self, z):
        if not self._valid:
            return 1.0
        rho_nr = self._rho_nr_ref * np.exp(self.u(z))
        return max((rho_nr + self.Omega_r * (1 + z) ** 4) / self._E2_norm, 1e-30)

    def E(self, z):
        return np.sqrt(self.E_squared(z))

    # -- distances (slow / exact) -------------------------------------------
    def comoving_distance(self, z):
        result, _ = quad(lambda zp: 1.0 / self.E(zp), 0, z,
                         limit=200, epsabs=1e-10, epsrel=1e-10)
        return (C_LIGHT / self.H0) * result

    def luminosity_distance(self, z):
        return (1 + z) * self.comoving_distance(z)

    def distance_modulus(self, z):
        dL = self.luminosity_distance(z)
        return np.inf if dL <= 0 else 5.0 * np.log10(dL) + 25.0

    def DM(self, z):
        return self.comoving_distance(z)

    def DH(self, z):
        return C_LIGHT / (self.H0 * self.E(z))

    def DV(self, z):
        DA = self.comoving_distance(z) / (1 + z)
        Hz = self.H0 * self.E(z)
        return ((1 + z) ** 2 * DA ** 2 * C_LIGHT * z / Hz) ** (1.0 / 3.0)

    # -- distances (fast interpolation) -------------------------------------
    def _build_dc_interp(self, z_max=2.6, n=2000):
        zg = np.linspace(0.0, z_max, n)
        Ng = -np.log(1.0 + zg)
        ys = self.sol.sol(Ng)                       # vectorised dense output
        rho_nr = self._rho_nr_ref * np.exp(ys[1])
        E2 = (rho_nr + self.Omega_r * (1 + zg) ** 4) / self._E2_norm
        inv_E = 1.0 / np.sqrt(np.clip(E2, 1e-30, None))
        dc = (C_LIGHT / self.H0) * cumulative_trapezoid(inv_E, zg, initial=0.0)
        self._dc_interp = (zg, dc)

    def comoving_distance_fast(self, z):
        if self._dc_interp is None:
            self._build_dc_interp()
        zg, dc = self._dc_interp
        return np.interp(z, zg, dc)

    def distance_modulus_fast(self, z):
        dL = (1.0 + np.asarray(z)) * self.comoving_distance_fast(z)
        with np.errstate(divide="ignore"):
            return 5.0 * np.log10(dL) + 25.0

    # -- diagnostics --------------------------------------------------------
    def find_z_acc(self):
        try:
            if self.R(0) > 2.0 / 3.0:
                return np.nan
            return brentq(lambda z: self.R(z) - 2.0 / 3.0, 0, 5, xtol=1e-6)
        except Exception:
            return np.nan


# ---------------------------------------------------------------------------
# fast cached self-consistent model (used inside the posterior)
# ---------------------------------------------------------------------------
_model_cache: dict = {}


def get_selfconsistent_model(K, H0, R_c, Delta_R, alpha):
    """Cached self-consistent model; returns None if shooting fails."""
    key = (round(alpha, 4), round(K, 5), round(H0, 3),
           round(R_c, 4), round(Delta_R, 4))
    if key in _model_cache:
        return _model_cache[key]
    m = ShootingRLDS(alpha, K, H0, R_c, Delta_R)
    if not m._valid or not (0.20 <= m.Omega_m <= 0.45):
        _model_cache[key] = None
        return None
    if len(_model_cache) > 300:
        _model_cache.clear()
    _model_cache[key] = m
    return m


# ---------------------------------------------------------------------------
# LCDM helpers (sound horizon for BAO)
# ---------------------------------------------------------------------------
def sound_horizon_Eh98(Omega_m, Omega_b, h):
    omega_m = Omega_m * h ** 2
    omega_b = Omega_b * h ** 2
    return 147.49 * (omega_m / 0.1415) ** (-0.255) * (omega_b / 0.02236) ** (-0.128)


# ===========================================================================
# LIKELIHOODS
# ===========================================================================
def log_likelihood_SN(model, *, use_full_cov=False, fast_distance=False,
                      cov_inv=None):
    """Pantheon+ SN likelihood with analytic absolute-magnitude marginalisation.

    Diagonal (default) reproduces the reference run exactly. If ``use_full_cov``
    and a ``cov_inv`` (inverse covariance) is supplied, the full-covariance
    marginalised form is used instead.
    """
    z_sn, mu_obs, mu_err = load_pantheon_plus()
    if fast_distance:
        mu_th = model.distance_modulus_fast(z_sn)
    else:
        mu_th = np.array([model.distance_modulus(z) for z in z_sn])
    if np.any(~np.isfinite(mu_th)):
        return -np.inf
    delta = mu_obs - mu_th

    if use_full_cov:
        if cov_inv is None:
            raise ValueError("use_full_cov=True requires cov_inv.")
        A = delta @ cov_inv @ delta
        B = np.sum(cov_inv @ delta)
        Cc = np.sum(cov_inv)
        return -0.5 * (A - B * B / Cc)

    # diagonal: analytic marginalisation over offset (S1/S2/S3 trick)
    sigma2 = mu_err ** 2
    S1 = np.sum(1.0 / sigma2)
    S2 = np.sum(delta / sigma2)
    S3 = np.sum(delta ** 2 / sigma2)
    return -0.5 * (S3 - S2 ** 2 / S1)


def log_likelihood_BAO(model):
    h = model.H0 / 100.0
    r_d = sound_horizon_Eh98(model.Omega_m, OMEGA_B, h)
    entries, data, _cov, cov_inv = get_desi_dr2_bao()
    theory = np.empty(len(entries))
    for i, (z, quantity) in enumerate(entries):
        if quantity == "DM":
            th = model.DM(z) / r_d
        elif quantity == "DH":
            th = model.DH(z) / r_d
        elif quantity == "DV":
            th = model.DV(z) / r_d
        else:  # pragma: no cover - guarded by the loader's label map
            return -np.inf
        if not np.isfinite(th) or th <= 0:
            return -np.inf
        theory[i] = th
    delta = data - theory
    return -0.5 * float(delta @ cov_inv @ delta)


def log_likelihood_CMB(model):
    h = model.H0 / 100.0
    omega_b = OMEGA_B * h ** 2
    means, cov_inv = get_planck2018_distance_prior()
    DM_star = model.DM(PLANCK2018_ZSTAR)
    R_th = np.sqrt(model.Omega_m) * model.H0 * DM_star / C_LIGHT
    if not np.isfinite(R_th):
        return -np.inf
    delta = means - np.array([R_th, omega_b])
    return -0.5 * (delta @ cov_inv @ delta)


# ===========================================================================
# PRIOR + POSTERIOR
# ===========================================================================
def log_prior(theta):
    for val, name in zip(theta, PARAM_NAMES):
        lo, hi = PRIOR_BOUNDS[name]
        if not (lo <= val <= hi):
            return -np.inf
    return 0.0


def log_posterior(theta, *, use_full_cov=False, fast_distance=False,
                  sn_cov_inv=None):
    """Full log-posterior. Defaults reproduce the archived run (diagonal, slow)."""
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    K, H0, Rc, dR, alpha = theta
    model = get_selfconsistent_model(K, H0, Rc, dR, alpha)
    if model is None:
        return -np.inf

    ll_sn = log_likelihood_SN(model, use_full_cov=use_full_cov,
                              fast_distance=fast_distance, cov_inv=sn_cov_inv)
    if not np.isfinite(ll_sn):
        return -np.inf
    ll_bao = log_likelihood_BAO(model)
    if not np.isfinite(ll_bao):
        return -np.inf
    ll_cmb = log_likelihood_CMB(model)
    if not np.isfinite(ll_cmb):
        return -np.inf
    return lp + ll_sn + ll_bao + ll_cmb


def reproduce_reference_chi2(*, fast_distance=False):
    """chi2 at :data:`REFERENCE_BEST_FIT`; equals :data:`REFERENCE_CHI2` to < 0.01.

    A fast, exact regression check that the model + likelihood pipeline still
    reproduces the canonical DESI DR2 chain. Requires the Pantheon+ data file
    (see :func:`rlds_paired.data.load_pantheon_plus`).
    """
    return -2.0 * log_posterior(REFERENCE_BEST_FIT, fast_distance=fast_distance)


# ===========================================================================
# RUNNER
# ===========================================================================
def run_mcmc(n_walkers=32, n_steps=3000, n_burn=200, *, seed=None,
             output_dir="mcmc_output_shooting", fast_distance=True,
             use_full_cov=False, sn_cov_inv=None, progress=True):
    """Run the shooting MCMC and save the canonical ``.npy`` artefacts.

    Defaults match the reference run dimensions (32 walkers x 3000 steps,
    burn-in 200). ``fast_distance=True`` is recommended for new runs (the
    slow path is ~100x slower); set it False to follow the original code path.
    Set ``seed`` for a reproducible chain (the reference run used none).
    """
    import emcee

    rng = np.random.default_rng(seed)
    out = Path(output_dir)
    out.mkdir(exist_ok=True, parents=True)

    def lnpost(th):
        return log_posterior(th, use_full_cov=use_full_cov,
                             fast_distance=fast_distance, sn_cov_inv=sn_cov_inv)

    # initialise walkers near the reference start, clipped into the prior box
    pos = THETA_START + START_SCATTER * rng.standard_normal((n_walkers, NDIM))
    for i, name in enumerate(PARAM_NAMES):
        lo, hi = PRIOR_BOUNDS[name]
        pos[:, i] = np.clip(pos[:, i], lo + 1e-5, hi - 1e-5)
    for i in range(n_walkers):
        if not np.isfinite(lnpost(pos[i])):
            pos[i] = np.clip(
                THETA_START + 1e-3 * rng.standard_normal(NDIM),
                [PRIOR_BOUNDS[n][0] + 1e-5 for n in PARAM_NAMES],
                [PRIOR_BOUNDS[n][1] - 1e-5 for n in PARAM_NAMES],
            )

    _model_cache.clear()
    sampler = emcee.EnsembleSampler(n_walkers, NDIM, lnpost)
    sampler.run_mcmc(pos, n_steps, progress=progress)

    chains = sampler.get_chain()
    logprob = sampler.get_log_prob()
    flat = sampler.get_chain(discard=n_burn, flat=True)
    flat_lp = sampler.get_log_prob(discard=n_burn, flat=True)

    np.save(out / "chains_raw.npy", chains)
    np.save(out / "logprob_raw.npy", logprob)
    np.save(out / "flat_samples.npy", flat)
    np.save(out / "flat_logprob.npy", flat_lp)

    n_om = min(len(flat), 2000)
    Om = np.full(n_om, np.nan)
    for j in range(n_om):
        m = get_selfconsistent_model(*[flat[j, 0], flat[j, 1], flat[j, 2],
                                       flat[j, 3], flat[j, 4]])
        if m is not None:
            Om[j] = m.Omega_m
    np.save(out / "omega_m_derived.npy", Om)

    idx = int(np.argmax(flat_lp))
    summary = {
        "param_names": PARAM_NAMES,
        "theta_best": flat[idx],
        "chi2_best": float(-2 * flat_lp[idx]),
        "acceptance": float(np.mean(sampler.acceptance_fraction)),
        "seed": seed,
        "n_walkers": n_walkers,
        "n_steps": n_steps,
        "n_burn": n_burn,
    }
    np.save(out / "summary.npy", summary)
    return sampler, summary


__all__ = [
    "ShootingRLDS",
    "get_selfconsistent_model",
    "log_prior",
    "log_posterior",
    "log_likelihood_SN",
    "log_likelihood_BAO",
    "log_likelihood_CMB",
    "run_mcmc",
    "reproduce_reference_chi2",
    "REFERENCE_BEST_FIT",
    "REFERENCE_CHI2",
    "PARAM_NAMES",
    "PARAM_LABELS",
    "PRIOR_BOUNDS",
    "NDIM",
    "EPS",
]
