"""
Profile likelihood in dOm = Om_today - omega_m(early)/h^2 for the FREE-Om version of
RLDS (alpha is a free parameter, not fixed by the shooting condition), CMB+BAO+SN,
with the widened priors of the scan: Om in [0.2, OM_MAX], H0 [60, 75], wb [0.019, 0.026],
K [0, 0.25], Rc [0.3, 0.9], dR [0.01, 0.2], alpha [2.5, 3.5].
Optional Gaussian prior on Om_today (cluster f_gas, Mantz+2022: 0.260 +/- 0.040).

For each target dOm: minimise chi2 (+ prior) + ((dOm - target)/0.001)^2 over all seven
parameters (Nelder-Mead, scaled coordinates). Starts: points of the pre-registered
free-model MCMC (thinned), the widened-prior best fits, the scan hits, and the solution
at the previous target (continuation).
Note: dOm >= 0 for K >= 0; dOm -> 0 is the wCDM limit (w = -alpha/3), not LCDM.
Usage: profile_dOm_free.py [--om-prior 0.260,0.040] [--om-max 0.60]
"""
import argparse, json, os, time
import numpy as np
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L

ap = argparse.ArgumentParser()
ap.add_argument("--om-prior", default=None)
ap.add_argument("--om-max", type=float, default=0.60)
ap.add_argument("--targets", default="0.002,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.10,0.12,0.14,0.17")
args = ap.parse_args()
prior = [float(v) for v in args.om_prior.split(",")] if args.om_prior else None
TARGETS = [float(v) for v in args.targets.split(",")]
NAMES = ["Om", "H0", "wb", "K", "Rc", "dR", "alpha"]
LO = np.array([0.20, 60, 0.019, 0.0, 0.30, 0.01, 2.5])
HI = np.array([args.om_max, 75, 0.026, 0.25, 0.90, 0.20, 3.5])
SC = np.array([0.01, 0.5, 0.0002, 0.01, 0.05, 0.02, 0.02])
PEN = 0.001
sn = L.load_sn(); use = ("cmb", "bao", "sn")
tag = "free" + (f"_om{prior[0]:g}-{prior[1]:g}" if prior else "") + f"_ommax{args.om_max:g}"


def evaluate(x):
    if np.any(x < LO) or np.any(x > HI):
        return None
    try:
        bg = L.rlds_free(*x)
        c = L.chi2_total(bg, sn, use)
    except Exception:
        return None
    if not np.isfinite(c):
        return None
    if prior:
        c += ((x[0] - prior[0]) / prior[1]) ** 2
    return c, bg.Omega_m_today - bg.omega_m_early / bg.h ** 2


# ---------- start pool ----------
pool_file = "validation/profile_dOm_free_pool.json"
if os.path.exists(pool_file):
    pool = json.load(open(pool_file))
else:
    import emcee
    b = emcee.backends.HDFBackend("validation/chains/RLDS-free.h5", read_only=True)
    ch = b.get_chain(discard=3000, thin=20, flat=True)
    rng = np.random.default_rng(1)
    pool = []
    for x in ch[rng.choice(len(ch), 1500, replace=False)]:
        try:
            bg = L.rlds_free(*x); c = L.chi2_total(bg, sn, use)
        except Exception:
            continue
        if np.isfinite(c):
            pool.append(dict(x=list(map(float, x)), chi2=float(c), dOm=float(bg.Omega_m_today - bg.omega_m_early / bg.h ** 2)))
    for key, f in [("none", "validation/fit_free_wide.json"), ("fgas", "validation/fit_free_wide.json")]:
        p = json.load(open(f))[key]["params"]
        x = [p[k] for k in NAMES]; bg = L.rlds_free(*x)
        pool.append(dict(x=x, chi2=float(L.chi2_total(bg, sn, use)), dOm=float(bg.Omega_m_today - bg.omega_m_early / bg.h ** 2)))
    for r in json.load(open("validation/fit_free_wide_om060.json")):
        p = r["params"]; x = [p[k] for k in NAMES]; bg = L.rlds_free(*x)
        pool.append(dict(x=x, chi2=float(L.chi2_total(bg, sn, use)), dOm=float(bg.Omega_m_today - bg.omega_m_early / bg.h ** 2)))
    for r in json.load(open("validation/scan_hits_chi2.json")):
        if r["kind"] != "free":
            continue
        x = [r["Om"], r["H0"], r["wb"], r["K"], r["Rc"], r["dR"], r["alpha"]]
        pool.append(dict(x=x, chi2=r["chi2"], dOm=r["dOm"]))
    json.dump(pool, open(pool_file, "w"))
print("start pool", len(pool), flush=True)

score = lambda p: p["chi2"] + (((p["x"][0] - prior[0]) / prior[1]) ** 2 if prior else 0.0)
rows, prev = [], None
for t in TARGETS:
    t0 = time.time()
    inside = [p for p in pool if abs(p["dOm"] - t) < 0.015 and np.all(np.array(p["x"]) <= HI) and np.all(np.array(p["x"]) >= LO)]
    starts = [np.array(p["x"]) for p in sorted(inside, key=score)[:2]]
    if not starts:
        starts = [np.array(min(pool, key=lambda p: abs(p["dOm"] - t))["x"])]
    if prev is not None:
        starts.append(prev)
    def f(y):
        r = evaluate(y * SC)
        return 1e10 if r is None else r[0] + ((r[1] - t) / PEN) ** 2
    best = None
    for s in starts:
        r = minimize(f, np.clip(s, LO + 1e-7, HI - 1e-7) / SC, method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-4, maxfev=2500, adaptive=True))
        if best is None or r.fun < best.fun:
            best = r
    x = best.x * SC
    c, d = evaluate(x)
    prev = x
    row = dict(target=t, dOm=float(d), chi2=float(c), **dict(zip(NAMES, map(float, x))),
               at_edge=[n for n, v, lo, hi in zip(NAMES, x, LO, HI) if min(v - lo, hi - v) < 1e-3 * (hi - lo)],
               n_starts=len(starts), seconds=round(time.time() - t0))
    rows.append(row)
    print(json.dumps(row), flush=True)
    json.dump(rows, open(f"validation/profile_dOm_{tag}.json", "w"), indent=1)
