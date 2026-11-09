"""
Profile likelihood in the derived quantity dOm = Om_today - omega_m(early)/h^2
for RLDS-shooting (the single combination the background data constrain, see
derived_observables.py). For each target dOm: minimise chi2 over
(K, H0, wb, Om, Rc, dR), alpha from the shooting condition, with a stiff
penalty ((dOm - target)/0.001)^2 holding dOm at the target.
Starts: the best K-profile points of the same data combination near the target.
Usage: profile_dOm.py cmb,bao,sn [--om-prior MEAN,SIGMA]
"""
import argparse, json, time
import numpy as np
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
import profile_K as P

NAMES = ["K", "H0", "wb", "Om", "Rc", "dR"]
LO = np.array([0.0005, 60, 0.019, 0.20, 0.30, 0.03])
HI = np.array([0.25, 75, 0.026, 0.45, 0.60, 0.20])
SCALE = np.array([0.01, 0.5, 0.0002, 0.01, 0.05, 0.03])
TARGETS = [0.002, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10]
PEN = 0.001

ap = argparse.ArgumentParser(); ap.add_argument("use"); ap.add_argument("--om-prior", default=None)
ap.add_argument("--rc-max", type=float, default=0.6); ap.add_argument("--dr-min", type=float, default=0.03)
ap.add_argument("--targets", default=None)
args = ap.parse_args()
HI[4] = args.rc_max; LO[5] = args.dr_min
if args.targets: TARGETS = [float(v) for v in args.targets.split(",")]
wide = args.rc_max != 0.6 or args.dr_min != 0.03
use = tuple(args.use.split(","))
prior = [float(v) for v in args.om_prior.split(",")] if args.om_prior else None
sn = L.load_sn() if "sn" in use else None
src = "_".join(use) + (f"_om{prior[0]:g}-{prior[1]:g}" if prior else "")
pool = json.load(open("validation/derived_observables.json"))
pts = pool[src]
extra = pool["_".join(use) + "_rc0.8"] if wide else []
out_tag = src + (f"_rc{args.rc_max:g}_dr{args.dr_min:g}" if wide else "") + (f"_t{TARGETS[0]:g}" if args.targets else "")


def evaluate(x):
    K, H0, wb, Om, Rc, dR = x
    a = P.alpha_from(K, Om, Rc, dR)
    if a is None:
        return None
    bg = L.rlds_free(Om, H0, wb, K, Rc, dR, a)
    c = L.chi2_total(bg, sn, use)
    if prior:
        c += ((Om - prior[0]) / prior[1]) ** 2
    return c, bg.Omega_m_today - bg.omega_m_early / bg.h ** 2, a


def objective(target):
    def f(y):
        x = y * SCALE
        if np.any(x < LO) or np.any(x > HI):
            return 1e10
        r = evaluate(x)
        if r is None or not np.isfinite(r[0]):
            return 1e10
        return r[0] + ((r[1] - target) / PEN) ** 2
    return f


rows = []
for t in TARGETS:
    t0 = time.time()
    key = lambda p: (abs(p["dOm"] - t) > 0.012, p["chi2"])
    near = sorted(pts, key=key)[:2] if not wide else [sorted(pts, key=key)[0], sorted(extra, key=key)[0]]
    if wide and rows:
        near.append(dict(rows[-1]))
    best = None
    for p in near:
        x0 = np.array([p["K"], p["H0"], p["wb"], p["Om"] + (t - p["dOm"]), p["Rc"], p["dR"]])
        y0 = np.clip(x0, LO + 1e-6, HI - 1e-6) / SCALE
        r = minimize(objective(t), y0, method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-3, maxfev=700, adaptive=True))
        if best is None or r.fun < best.fun:
            best = r
    x = best.x * SCALE
    c, d, a = evaluate(x)
    row = dict(target=t, dOm=d, chi2=c, alpha=a, **dict(zip(NAMES, map(float, x))),
               seconds=round(time.time() - t0))
    rows.append(row)
    print(json.dumps(row), flush=True)
    json.dump(rows, open(f"validation/profile_dOm_{out_tag}.json", "w"), indent=1)
