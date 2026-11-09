"""
Controlled comparison for the old Delta chi2 ~ +36 of Papers III-IV (review item 1):
change ONE element at a time -- likelihood (old / new) x priors (old / new) --
with the same optimiser everywhere, for the shooting version and for LCDM.

  L_old : rlds_paired.mcmc likelihood of Papers III-IV (CMB: R only at fixed z*,
          omega_b tied to H0 through fixed Omega_b; r_d from the EH98-type fit with
          today's Omega_m; Pantheon+ with diagonal errors).
  L_new : rlds_paired.likelihood_v2 (CMB R, l_A, omega_b with early omega_m;
          Aubourg r_d; Pantheon+ full covariance).
  P_old : K in [0.005, 0.25], alpha in [3.01, 3.30] (mcmc.PRIOR_BOUNDS).
  P_new : K in [0.0005, 0.25], alpha in [3.0, 3.5] (derived).
  Common: H0 [60, 75], R_c [0.3, 0.6], dR [0.03, 0.2], Om [0.2, 0.45];
          omega_b [0.019, 0.026] in L_new only (tied in L_old).

Optimiser: Nelder-Mead in the (K, H0, [wb], Om, Rc, dR) chart with alpha from the
shooting condition (as in profile_K.py), several starts. The archived best fit of the
old chain is also evaluated (optimiser effect).
Usage: controlled_comparison.py old|new
"""
import json, sys
import numpy as np
from scipy.optimize import minimize
from rlds_paired import mcmc as M, likelihood_v2 as L
import profile_K as P

lik = sys.argv[1]
SN_NEW = L.load_sn() if lik == "new" else None


# ---------------- models in each likelihood ----------------
def old_model(K, H0, Om, Rc, dR, alpha):
    m = M.ShootingRLDS(alpha, K, H0, Rc, dR, _skip_shooting=True)
    m._finalise(Om)
    return m if m._valid else None


def chi2_old(model):
    try:
        ll = (M.log_likelihood_SN(model, fast_distance=True) + M.log_likelihood_BAO(model)
              + M.log_likelihood_CMB(model))
    except Exception:
        return np.inf
    return -2 * ll if np.isfinite(ll) else np.inf


def chi2_shoot(theta, prior):
    if lik == "old":
        K, H0, Om, Rc, dR = theta; wb = None
    else:
        K, H0, wb, Om, Rc, dR = theta
    a = P.alpha_from(K, Om, Rc, dR)
    if a is None:
        return 1e10, None
    a_lo, a_hi = (3.01, 3.30) if prior == "old" else (3.0, 3.5)
    if not a_lo <= a <= a_hi:
        return 1e10, a
    if lik == "old":
        m = old_model(K, H0, Om, Rc, dR, a)
        c = chi2_old(m) if m is not None else np.inf
    else:
        c = L.chi2_total(L.rlds_free(Om, H0, wb, K, Rc, dR, a), SN_NEW, ("cmb", "bao", "sn"))
    return (c if np.isfinite(c) else 1e10), a


def chi2_lcdm(theta):
    if lik == "old":
        Om, H0 = theta
        if not (0.2 <= Om <= 0.45 and 60 <= H0 <= 75):
            return 1e10
        m = old_model(0.0, H0, Om, 0.5, 0.1, 3.0)
        c = chi2_old(m) if m is not None else np.inf
    else:
        Om, H0, wb = theta
        if not (0.2 <= Om <= 0.45 and 60 <= H0 <= 75 and 0.019 <= wb <= 0.026):
            return 1e10
        c = L.chi2_total(L.lcdm(Om, H0, wb), SN_NEW, ("cmb", "bao", "sn"))
    return c if np.isfinite(c) else 1e10


# ---------------- bounds / scales ----------------
if lik == "old":
    NAMES = ["K", "H0", "Om", "Rc", "dR"]
    SC = np.array([0.005, 0.5, 0.01, 0.05, 0.03])
    BLO = {"old": [0.005, 60, 0.2, 0.3, 0.03], "new": [0.0005, 60, 0.2, 0.3, 0.03]}
    BHI = [0.25, 75, 0.45, 0.6, 0.2]
else:
    NAMES = ["K", "H0", "wb", "Om", "Rc", "dR"]
    SC = np.array([0.005, 0.5, 0.0002, 0.01, 0.05, 0.03])
    BLO = {"old": [0.005, 60, 0.019, 0.2, 0.3, 0.03], "new": [0.0005, 60, 0.019, 0.2, 0.3, 0.03]}
    BHI = [0.25, 75, 0.026, 0.45, 0.6, 0.2]


def fit_shoot(prior, starts):
    lo, hi = np.array(BLO[prior]), np.array(BHI)
    def f(y):
        x = y * SC
        if np.any(x < lo) or np.any(x > hi):
            return 1e10
        return chi2_shoot(x, prior)[0]
    best = None
    for s in starts:
        r = minimize(f, np.clip(s, lo + 1e-6, hi - 1e-6) / SC, method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-4, maxfev=1500, adaptive=True))
        print(f"  {prior} start {np.round(s, 4).tolist()} -> {r.fun:.3f}", flush=True)
        if best is None or r.fun < best.fun:
            best = r
    x = best.x * SC
    c, a = chi2_shoot(x, prior)
    return dict(chi2=float(c), alpha=a, **dict(zip(NAMES, map(float, x))))


# archived best fit of the old chain: K, H0, Rc, dR, alpha
KA, HA, RcA, dRA, aA = M.REFERENCE_BEST_FIT
OmA = M.ShootingRLDS(aA, KA, HA, RcA, dRA).Omega_m
out = dict(likelihood=lik)

# LCDM
l0 = [0.30, 68.0] if lik == "old" else [0.30, 68.0, 0.0224]
r = minimize(chi2_lcdm, l0, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-5, maxfev=3000))
out["lcdm"] = dict(chi2=float(r.fun), params=list(map(float, r.x)))
print("LCDM", out["lcdm"], flush=True)

# archived point evaluated in this likelihood (alpha re-derived from Om to stay on the ridge)
arch = [KA, HA, OmA, RcA, dRA] if lik == "old" else [KA, HA, 0.02237, OmA, RcA, dRA]
out["archived_point"] = dict(chi2=float(chi2_shoot(np.array(arch), "old")[0]), theta=list(map(float, arch)))
print("archived point", out["archived_point"], flush=True)

new_best = json.load(open("validation/profile_K_cmb_bao_sn.json"))
nb = min(new_best, key=lambda r: r["chi2"])
st_new = [nb["K"], nb["H0"], nb["Om"], nb["Rc"], nb["dR"]] if lik == "old" else [nb["K"], nb["H0"], nb["wb"], nb["Om"], nb["Rc"], nb["dR"]]
st_small = [0.01, 68.0, 0.30, 0.55, 0.15] if lik == "old" else [0.01, 68.0, 0.0224, 0.30, 0.55, 0.15]
for prior in ("old", "new"):
    starts = [np.array(arch), np.array(st_new), np.array(st_small)]
    res = fit_shoot(prior, starts)
    res["dchi2_vs_lcdm"] = res["chi2"] - out["lcdm"]["chi2"]
    out[f"shoot_prior_{prior}"] = res
    print(f"SHOOT prior={prior}", json.dumps(res), flush=True)
json.dump(out, open(f"validation/controlled_{lik}.json", "w"), indent=1)
