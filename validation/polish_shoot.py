"""Final polishing of the RLDS-shooting best fits listed in Table 2 of Paper VI.

Near its minimum the shooting likelihood is flat in K and dR, and the best fits sit
at the R_c prior edge. Nelder-Mead with a hard wall at the prior edges stalls there,
so minima found by different optimisations differ by a few hundredths in chi2: the
K profile on its grid (profile_K.py), the controlled comparison of Appendix C
(controlled_comparison.py new) and fit_shoot_wide.py. This script starts from the
lowest of these points inside each prior set and restarts Nelder-Mead with a fresh
simplex, evaluating chi2 at the point clipped to the prior box plus a quadratic
penalty outside it (so the simplex can slide along an edge), until two consecutive
restarts improve chi2 by less than 2e-4. Table 2 lists the resulting values.

Usage : polish_shoot.py prereg|rc0.8|wide
Output: validation/polish_shoot_<name>.json
"""
import json, sys, numpy as np
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
import profile_K as P
import derived_observables as D

NAMES = ["K", "H0", "wb", "Om", "Rc", "dR"]
CONF = {"prereg": (0.60, 0.03), "rc0.8": (0.80, 0.03), "wide": (0.90, 0.01)}  # (Rc max, dR min)
name = sys.argv[1]; rc_max, dr_min = CONF[name]
LO = np.array([0.0005, 60, 0.019, 0.20, 0.30, dr_min]); HI = np.array([0.25, 75, 0.026, 0.45, rc_max, 0.20])
SC = np.array([0.005, 0.5, 0.0002, 0.01, 0.05, 0.03]); LCDM = 1421.359
sn = L.load_sn(); use = ("cmb", "bao", "sn")


def chi2(x):
    K, H0, wb, Om, Rc, dR = x
    a = P.alpha_from(K, Om, Rc, dR)
    if a is None:
        return 1e10
    v = L.chi2_total(L.rlds_free(Om, H0, wb, K, Rc, dR, a), sn, use)
    return v if np.isfinite(v) else 1e10


def f(y):
    x = y * SC; xc = np.clip(x, LO, HI)
    return chi2(xc) + 1e3 * float(np.sum(((x - xc) / SC) ** 2))


def vec(r):
    p = r.get("params", r)
    return np.array([p[k] for k in NAMES], dtype=float)


# candidate starts: every earlier minimum that lies inside this prior set
cand = []
for fn in ["validation/profile_K_cmb_bao_sn.json", "validation/profile_K_cmb_bao_sn_rc0.8.json"]:
    cand += [vec(r) for r in json.load(open(fn))]
c = json.load(open("validation/controlled_new.json"))["shoot_prior_new"]
cand.append(np.array([c[k] for k in NAMES]))
cand += [vec(r) for r in json.load(open("validation/fit_shoot_wide.json"))]
cand = [x for x in cand if np.all(x >= LO - 1e-9) and np.all(x <= HI + 1e-9)]
vals = [chi2(np.clip(x, LO, HI)) for x in cand]
i0 = int(np.argmin(vals)); x = np.clip(cand[i0], LO, HI); best = vals[i0]; start = best
print(f"{name}: best earlier minimum dchi2 = {best - LCDM:.4f}", flush=True)
small = 0
for it in range(12):
    r = minimize(f, x / SC, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-5, maxfev=1200, adaptive=True))
    xn = np.clip(r.x * SC, LO, HI); cn = chi2(xn); gain = best - cn
    if cn < best:
        best, x = cn, xn
    print(f"  restart {it}: dchi2 = {best - LCDM:.4f} (gain {gain:+.4f}, nfev {r.nfev})", flush=True)
    small = small + 1 if gain < 2e-4 else 0
    if small == 2:
        break
a = P.alpha_from(x[0], x[3], x[4], x[5])
p = dict(zip(NAMES, map(float, x))); p["alpha"] = float(a)
d = D.derived(dict(p))
row = dict(prior=name, Rc_max=rc_max, dR_min=dr_min, params=p, chi2=float(best), dchi2=float(best - LCDM),
           dchi2_start=float(start - LCDM), at_edge=[k for k, v, lo, hi in zip(NAMES, x, LO, HI) if v <= lo + 1e-6 or v >= hi - 1e-6],
           dOm=d["dOm"], w0=d["w0_fit"], wa=d["wa_fit"], Om_early=d["Om_early"])
print(json.dumps(row), flush=True)
json.dump(row, open(f"validation/polish_shoot_{name}.json", "w"), indent=1)
