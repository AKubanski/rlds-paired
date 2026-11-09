"""
Full optimisation of RLDS-free with the WIDE priors used in the scan
(Rc <= 0.9, dR >= 0.01, alpha in [2.5, 3.5]), CMB+BAO+SN, with and without
the f_gas prior. Starts: the best scan hits (scan_hits_chi2.json) and the
maximum-likelihood point of the pre-registered MCMC. Nelder-Mead in scaled
coordinates, then a restart from the best point.
"""
import json
import numpy as np
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
import derived_observables as D

NAMES = ["Om", "H0", "wb", "K", "Rc", "dR", "alpha"]
LO = np.array([0.20, 60, 0.019, 0.0, 0.30, 0.01, 2.5])
HI = np.array([0.45, 75, 0.026, 0.25, 0.90, 0.20, 3.5])
SC = np.array([0.01, 0.5, 0.0002, 0.01, 0.05, 0.02, 0.02])
LCDM, LCDM_F = 1421.359, 1422.423
sn = L.load_sn(); use = ("cmb", "bao", "sn")


def chi2(x, fgas):
    if np.any(x < LO) or np.any(x > HI):
        return 1e10
    bg = L.rlds_free(*x)
    v = L.chi2_total(bg, sn, use) if bg is not None else np.inf
    if not np.isfinite(v):
        return 1e10
    return v + (((x[0] - 0.26) / 0.04) ** 2 if fgas else 0.0)


hits = sorted(json.load(open("validation/scan_hits_chi2.json")), key=lambda r: r["dchi2"])
starts = [np.array([h[k] for k in ["Om", "H0", "wb", "K", "Rc", "dR", "alpha"]]) for h in hits[:4]]
starts.append(np.array([0.4439, 67.4535, 0.0225, 0.207, 0.5414, 0.1806, 3.1959]))   # MCMC max (pre-registered)
starts.append(np.array([0.33, 68.0, 0.0225, 0.03, 0.75, 0.03, 3.03]))               # moderate start
res = {}
for fgas in (False, True):
    best = None
    for s in starts:
        f = lambda y: chi2(y * SC, fgas)
        r = minimize(f, np.clip(s, LO + 1e-6, HI - 1e-6) / SC, method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-4, maxfev=3000, adaptive=True))
        r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-5, maxfev=2000, adaptive=True))
        print("fgas" if fgas else "none", "start", np.round(s, 3).tolist(), "->", round(r.fun - (LCDM_F if fgas else LCDM), 3), flush=True)
        if best is None or r.fun < best.fun:
            best = r
    x = best.x * SC
    p = dict(zip(NAMES, map(float, x)))
    d = D.derived(dict(p, Om=p["Om"]))
    res["fgas" if fgas else "none"] = dict(params=p, chi2=float(best.fun),
                                           dchi2=float(best.fun - (LCDM_F if fgas else LCDM)),
                                           at_edge={n: bool(abs(v - lo) < 1e-3 * (hi - lo) or abs(v - hi) < 1e-3 * (hi - lo))
                                                    for n, v, lo, hi in zip(NAMES, x, LO, HI)},
                                           derived=d)
    print(json.dumps(res["fgas" if fgas else "none"]), flush=True)
json.dump(res, open("validation/fit_free_wide.json", "w"), indent=1)
