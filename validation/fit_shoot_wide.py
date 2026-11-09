"""Global optimisation of RLDS-shooting with the widened priors (Rc <= 0.9, dR >= 0.01),
CMB+BAO+SN, no f_gas prior. alpha from the shooting condition (profile_K.alpha_from)."""
import json, numpy as np
from scipy.optimize import minimize
from rlds_paired import likelihood_v2 as L
import profile_K as P
import derived_observables as D
NAMES = ["K", "H0", "wb", "Om", "Rc", "dR"]
LO = np.array([0.0005, 60, 0.019, 0.20, 0.30, 0.01]); HI = np.array([0.25, 75, 0.026, 0.45, 0.90, 0.20])
SC = np.array([0.01, 0.5, 0.0002, 0.01, 0.05, 0.02]); LCDM = 1421.359
sn = L.load_sn(); use = ("cmb", "bao", "sn")
def chi2(x):
    if np.any(x < LO) or np.any(x > HI): return 1e10
    K, H0, wb, Om, Rc, dR = x
    a = P.alpha_from(K, Om, Rc, dR)
    if a is None: return 1e10
    v = L.chi2_total(L.rlds_free(Om, H0, wb, K, Rc, dR, a), sn, use)
    return v if np.isfinite(v) else 1e10
rc8 = {r["K"]: r for r in json.load(open("validation/profile_K_cmb_bao_sn_rc0.8.json"))}[0.08]
starts = [np.array([rc8["K"], rc8["H0"], rc8["wb"], rc8["Om"], rc8["Rc"], rc8["dR"]]),
          np.array([0.12, 67.3, 0.0225, 0.42, 0.75, 0.02])]
out = []
for s in starts:
    f = lambda y: chi2(y * SC)
    r = minimize(f, np.clip(s, LO + 1e-6, HI - 1e-6) / SC, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-4, maxfev=2500, adaptive=True))
    x = r.x * SC; a = P.alpha_from(x[0], x[3], x[4], x[5])
    p = dict(zip(NAMES, map(float, x))); p["alpha"] = a
    d = D.derived(dict(Om=p["Om"], H0=p["H0"], wb=p["wb"], K=p["K"], Rc=p["Rc"], dR=p["dR"], alpha=a))
    row = dict(params=p, dchi2=float(r.fun - LCDM), dOm=d["dOm"], w0=d["w0_fit"], wa=d["wa_fit"], Om_early=d["Om_early"])
    print(json.dumps(row), flush=True); out.append(row)
json.dump(out, open("validation/fit_shoot_wide.json", "w"), indent=1)
