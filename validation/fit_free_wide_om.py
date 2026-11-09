"""
Edge check (Om <= 0.60) for the full optimisation of RLDS-free with the WIDE priors used in the scan
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



HI[0] = 0.60
best = json.load(open("validation/fit_free_wide.json"))["none"]["params"]
starts = [np.array([best[k] for k in NAMES]), np.array([0.50, 67.5, 0.0225, 0.15, 0.75, 0.03, 3.4])]
out = []
for s in starts:
    f = lambda y: chi2(y * SC, False)
    r = minimize(f, np.clip(s, LO + 1e-6, HI - 1e-6) / SC, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-4, maxfev=3000, adaptive=True))
    r = minimize(f, r.x, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-5, maxfev=2000, adaptive=True))
    x = r.x * SC; p = dict(zip(NAMES, map(float, x))); d = D.derived(p)
    row = dict(params=p, dchi2=float(r.fun - LCDM), dOm=d["dOm"], w0=d["w0_fit"], wa=d["wa_fit"], Om_early=d["Om_early"])
    print(json.dumps(row), flush=True); out.append(row)
json.dump(out, open("validation/fit_free_wide_om060.json", "w"), indent=1)
