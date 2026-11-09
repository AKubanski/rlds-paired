"""Small numerical checks quoted in RLDS VI (Sec. 2.2, 2.3, App. A):
 (a) Eq. (5): analytic w_eff(0) vs numerical w_eff(0) at all K-profile points;
 (b) regulariser: epsilon = 1e-6 vs 1e-12 at the shooting and free best fits;
 (c) shooting ridge: change of the derived Om for a change of alpha by 0.005 at K = 0.01.
Output: validation/checks_misc.json"""
import json
import numpy as np
from rlds_paired import mcmc as M, likelihood_v2 as L
import profile_K as P

out = {}
# (a)
d = json.load(open("validation/derived_observables.json"))
err = []
for rows in d.values():
    for r in rows:
        Or = L.OMEGA_R_H2 / (r["H0"] / 100) ** 2
        pred = -(r["alpha"] / 3) * (1 - r["dOm"] / (1 - Or - r["Om_early"]))
        err.append(abs(pred - r["w_z0"]))
out["eq5_points"] = len(err); out["eq5_max_abs_diff"] = float(max(err))
# (b)
sn = L.load_sn(); use = ("cmb", "bao", "sn")
shoot = min(json.load(open("validation/profile_K_cmb_bao_sn.json")), key=lambda r: r["chi2"])
free = json.load(open("validation/fit_free_wide.json"))["none"]["params"]
def chis():
    a = P.alpha_from(shoot["K"], shoot["Om"], shoot["Rc"], shoot["dR"])
    c1 = L.chi2_total(L.rlds_free(shoot["Om"], shoot["H0"], shoot["wb"], shoot["K"], shoot["Rc"], shoot["dR"], a), sn, use)
    c2 = L.chi2_total(L.rlds_free(*[free[k] for k in ["Om", "H0", "wb", "K", "Rc", "dR", "alpha"]]), sn, use)
    return c1, c2
c6 = chis()
M.Psi = lambda R, eps=1e-12: (1.0 - R) / (R + eps)
c12 = chis()
M.Psi = lambda R, eps=M.EPS: (1.0 - R) / (R + eps)
out["eps_dchi2_shooting"] = float(c12[0] - c6[0]); out["eps_dchi2_free"] = float(c12[1] - c6[1])
# (c)
for Rc, dR in ((0.565, 0.17), (0.513, 0.03)):
    oms = [M.ShootingRLDS(a, 0.01, 68.0, Rc, dR).Omega_m for a in (3.005, 3.010)]
    out[f"ridge_K0.01_Rc{Rc}_dR{dR}_alpha3.005_3.010_Om"] = list(map(float, oms))
json.dump(out, open("validation/checks_misc.json", "w"), indent=1)
print(json.dumps(out, indent=1))
