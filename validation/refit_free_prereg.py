"""Re-optimisation of RLDS-free with the PRE-REGISTERED priors (review of draft v4):
the reviewer found chi2 = 1415.0386 (dchi2 = -6.32) at a point below the MCMC maximum (-5.70).
Evaluate that point, then optimise from it and from other starts (Nelder-Mead, scaled)."""
import json, numpy as np
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
NAMES = ["Om", "H0", "wb", "K", "Rc", "dR", "alpha"]
LO = np.array([0.20, 60, 0.019, 0.0, 0.30, 0.03, 2.5]); HI = np.array([0.45, 75, 0.026, 0.25, 0.60, 0.20, 3.5])
SC = np.array([0.01, 0.5, 0.0002, 0.01, 0.05, 0.02, 0.02]); LCDM = 1421.359008662079
sn = L.load_sn(); use = ("cmb", "bao", "sn")
def chi2(x):
    if np.any(x < LO) or np.any(x > HI): return 1e10
    try:
        c = L.chi2_total(L.rlds_free(*x), sn, use)
    except Exception:
        return 1e10
    return c if np.isfinite(c) else 1e10
rev = np.array([0.44999999999954526, 67.66520543096937, 0.022487543321884865, 0.18915279481079547, 0.5949279502080446, 0.1953969566469472, 3.277824564260472])
print("reviewer point chi2 =", chi2(rev), flush=True)
starts = [rev, np.array([0.4439, 67.4535, 0.0225, 0.207, 0.5414, 0.1806, 3.1959]),
          np.array([0.43, 67.6, 0.0225, 0.12, 0.59, 0.05, 3.3]), np.array([0.40, 67.8, 0.0225, 0.08, 0.58, 0.10, 3.15])]
out = []
for s in starts:
    f = lambda y: chi2(y * SC)
    r = minimize(f, np.clip(s, LO + 1e-7, HI - 1e-7) / SC, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-5, maxfev=4000, adaptive=True))
    r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxfev=3000, adaptive=True))
    x = r.x * SC
    row = dict(start=list(map(float, s)), chi2=float(r.fun), dchi2=float(r.fun - LCDM), params=dict(zip(NAMES, map(float, x))),
               at_edge=[n for n, v, lo, hi in zip(NAMES, x, LO, HI) if min(v - lo, hi - v) < 1e-3 * (hi - lo)])
    print(json.dumps(row), flush=True); out.append(row)
best = min(out, key=lambda r: r["chi2"])
bg = L.rlds_free(*[best["params"][k] for k in NAMES])
best["dOm"] = float(bg.Omega_m_today - bg.omega_m_early / bg.h ** 2)
json.dump(dict(runs=out, best=best, reviewer_point_chi2=float(chi2(rev))), open("validation/refit_free_prereg.json", "w"), indent=1)
print("BEST", json.dumps(best))
