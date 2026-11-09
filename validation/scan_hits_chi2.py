"""
Likelihood of the scan points that land inside the CPL contours (scan_w0wa.py).
For each hit with CMB-consistent early matter (0.29 < Om_early < 0.33) the
shape parameters (Om, K, Rc, dR, alpha) are kept fixed and chi2(CMB+BAO+SN)
is minimised over (H0, omega_b) only -- a partial profile, i.e. an upper bound
on the best chi2 reachable at that point of the (w0, wa) projection.
Reported relative to LCDM, without and with the f_gas Omega_m prior.
"""
import json
import numpy as np
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.path import Path
import emcee
from rlds_paired import likelihood_v2 as L

LCDM, LCDM_F = 1421.359, 1422.423
b = emcee.backends.HDFBackend("validation/chains/CPL.h5", read_only=True)
c = b.get_chain(discard=1000, thin=10, flat=True)
H, xe, ye = np.histogram2d(c[:, 3], c[:, 4], bins=40)
Hs = np.sort(H.ravel())[::-1]; cs = np.cumsum(Hs) / Hs.sum()
lev = [Hs[np.searchsorted(cs, q)] for q in (0.95, 0.68)]
cn = plt.contour(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, levels=lev)
paths = {q: [Path(s) for s in segs if len(s) > 10] for q, segs in zip((95, 68), cn.allsegs)}
inside = lambda q, p: any(P.contains_point(p) for P in paths[q])

sn = L.load_sn()
use = ("cmb", "bao", "sn")
out = []
for kind, seed in [("free", 1), ("shooting", 2)]:
    pts = [x for x in json.load(open(f"validation/scan_w0wa_{kind}_{seed}.json"))
           if x["rhoDE_min"] > 0 and 0.29 < x["Om_early"] < 0.33
           and inside(95, (x["w0_fit"], x["wa_fit"]))]
    for x in pts:
        def f(y):
            H0, wb = y
            if not (60 < H0 < 75 and 0.019 < wb < 0.026):
                return 1e10
            bg = L.rlds_free(x["Om"], H0, wb, x["K"], x["Rc"], x["dR"], x["alpha"])
            v = L.chi2_total(bg, sn, use) if bg is not None else np.inf
            return v if np.isfinite(v) else 1e10
        r = minimize(f, [68.0, 0.0225], method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-3, maxfev=300, initial_simplex=[[68, 0.0225], [69, 0.0225], [68, 0.0228]]))
        pri = ((x["Om"] - 0.26) / 0.04) ** 2
        row = dict(kind=kind, in68=inside(68, (x["w0_fit"], x["wa_fit"])), w0=x["w0_fit"], wa=x["wa_fit"],
                   dOm=x["dOm"], Om=x["Om"], Om_early=x["Om_early"], K=x["K"], Rc=x["Rc"], dR=x["dR"], alpha=x["alpha"],
                   H0=float(r.x[0]), wb=float(r.x[1]), chi2=float(r.fun),
                   dchi2=float(r.fun - LCDM), dchi2_fgas=float(r.fun + pri - LCDM_F))
        out.append(row)
        print(json.dumps(row), flush=True)
json.dump(out, open("validation/scan_hits_chi2.json", "w"), indent=1)
