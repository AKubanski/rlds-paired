"""
Profile likelihood chi2(K) for the self-consistent (shooting) RLDS model.
Replaces the MCMC on K (analysis-plan amendment of 2026-11-05): K is degenerate
with alpha, Omega_m and R_c (Paper VI, Sec. 4.2), so a Bayesian bound on K would
be prior-dominated; the profile likelihood is not.

Reparametrisation: on the shooting ridge alpha is a steep function of the other
parameters, so we sample Omega_m instead and solve the shooting condition
R(z=0; alpha, K, Omega_m, Rc, dR) = Omega_m for alpha (monotone, brentq).
This is the same model as likelihood_v2.rlds_shooting, only a better chart.

For each fixed K: minimise chi2 over (H0, omega_b, Omega_m, Rc, dR) with the
pre-registered priors of mcmc_v2 (alpha in [3.0, 3.5] is derived).
Usage: profile_K.py cmb,bao,sn [--om-prior MEAN,SIGMA] [--rc-max X]
  --om-prior : Gaussian prior on Omega_m today (independent probe, e.g. low-z
               cluster f_gas, Mantz et al. 2022: 0.260 +/- 0.040)
  --rc-max   : widen the upper prior edge of R_c (default 0.6, pre-registered)
"""
import json, sys, time
import numpy as np
from scipy.optimize import brentq, minimize
from rlds_paired import mcmc as M, likelihood_v2 as L

NAMES = ["H0", "wb", "Om", "Rc", "dR"]
LO = np.array([60, 0.019, 0.20, 0.30, 0.03])
HI = np.array([75, 0.026, 0.45, 0.60, 0.20])
SCALE = np.array([0.5, 0.0002, 0.01, 0.05, 0.03])
K_GRID = [0.002, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05,
          0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.25]
A_LO, A_HI = 3.0, 3.5
OM_PRIOR = []          # [mean, sigma] or empty


def alpha_from(K, Om, Rc, dR):
    f = lambda a: M._integrate_Ru(a, K, Om, Rc, dR)[0] - Om
    flo, fhi = f(A_LO), f(A_HI)
    if not (np.isfinite(flo) and np.isfinite(fhi)) or flo * fhi > 0:
        return None
    return brentq(f, A_LO, A_HI, xtol=1e-9)


def make_chi2(K, sn, use):
    def chi2(y):
        x = y * SCALE
        if np.any(x < LO) or np.any(x > HI):
            return 1e10
        H0, wb, Om, Rc, dR = x
        a = alpha_from(K, Om, Rc, dR)
        if a is None:
            return 1e10
        c = L.chi2_total(L.rlds_free(Om, H0, wb, K, Rc, dR, a), sn, use)
        if OM_PRIOR:
            c += ((Om - OM_PRIOR[0]) / OM_PRIOR[1]) ** 2
        return c if np.isfinite(c) else 1e10
    return chi2


def fit(K, x0, sn, use):
    f = make_chi2(K, sn, use)
    best = None
    for start in x0:
        y0 = np.clip(start, LO + 1e-6, HI - 1e-6) / SCALE
        r = minimize(f, y0, method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-3, maxfev=500, adaptive=True))
        if best is None or r.fun < best.fun:
            best = r
    x = best.x * SCALE
    return float(best.fun), x, alpha_from(K, *x[2:])


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("use")
    ap.add_argument("--om-prior", default=None); ap.add_argument("--rc-max", type=float, default=0.6)
    args = ap.parse_args()
    use = tuple(args.use.split(","))
    HI[3] = args.rc_max
    if args.om_prior:
        OM_PRIOR[:] = [float(v) for v in args.om_prior.split(",")]
    sn = L.load_sn() if "sn" in use else None
    tag = "_".join(use) + (f"_om{OM_PRIOR[0]:g}-{OM_PRIOR[1]:g}" if OM_PRIOR else "") \
        + (f"_rc{args.rc_max:g}" if args.rc_max != 0.6 else "")
    default = np.array([68.2, 0.02252, 0.325, 0.513, 0.05])
    prev = default.copy()
    rows = []
    for K in K_GRID:
        t = time.time()
        starts = [prev, np.array([68.0, 0.0225, prev[2], 0.40, 0.12])]
        c, x, a = fit(K, starts, sn, use)
        prev = x
        chi2_prior = ((x[2] - OM_PRIOR[0]) / OM_PRIOR[1]) ** 2 if OM_PRIOR else 0.0
        row = dict(K=K, chi2=c, chi2_om_prior=float(chi2_prior), alpha=a, **dict(zip(NAMES, map(float, x))),
                   seconds=round(time.time() - t))
        rows.append(row)
        print(json.dumps(row), flush=True)
        json.dump(rows, open(f"validation/profile_K_{tag}.json", "w"), indent=1)
