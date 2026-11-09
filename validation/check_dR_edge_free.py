"""Edge check for profile_dOm_free.py (f_gas prior): dR sits at its lower bound 0.01 along the
profile. Re-optimise two profile points (dOm = 0.01 and 0.08) with dR >= 0.002."""
import json, numpy as np
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
NAMES = ["Om", "H0", "wb", "K", "Rc", "dR", "alpha"]
LO = np.array([0.20, 60, 0.019, 0.0, 0.30, 0.002, 2.5]); HI = np.array([0.60, 75, 0.026, 0.25, 0.90, 0.20, 3.5])
SC = np.array([0.01, 0.5, 0.0002, 0.01, 0.05, 0.02, 0.02])
sn = L.load_sn(); use = ("cmb", "bao", "sn")
rows = {round(r["target"], 3): r for r in json.load(open("validation/profile_dOm_free_om0.26-0.04_ommax0.6.json"))}
out = []
for t in (0.01, 0.08):
    r0 = rows[t]
    def f(y):
        x = y * SC
        if np.any(x < LO) or np.any(x > HI): return 1e10
        try:
            bg = L.rlds_free(*x); c = L.chi2_total(bg, sn, use)
        except Exception:
            return 1e10
        if not np.isfinite(c): return 1e10
        d = bg.Omega_m_today - bg.omega_m_early / bg.h ** 2
        return c + ((x[0] - 0.26) / 0.04) ** 2 + ((d - t) / 0.001) ** 2
    x0 = np.array([r0[k] for k in NAMES]); x0[5] = 0.005
    r = minimize(f, x0 / SC, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-4, maxfev=2500, adaptive=True))
    x = r.x * SC
    res = dict(target=t, chi2_dR001=r0["chi2"], chi2_dR0002=float(r.fun), dR=float(x[5]), improvement=float(r0["chi2"] - r.fun))
    print(json.dumps(res), flush=True); out.append(res)
json.dump(out, open("validation/check_dR_edge_free.json", "w"), indent=1)
