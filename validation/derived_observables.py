"""
Derived, data-facing quantities for every point of the K profiles.

For each profile point (K and the optimised H0, wb, Om, Rc, dR, alpha):
  Om_early  = omega_m(early)/h^2  -- Omega_m today if matter were conserved
  dOm       = Om_today - Om_early -- matter created by the coupling since z~1100
  rho_DE^eff(a) = rho_nr(a) - Om_early a^-3   (what an observer who assumes
                 conserved matter calls dark energy)
  w_eff(a)  = -1 - (1/3) dln rho_DE^eff / dln a
  (w0, wa)  = least-squares CPL fit of w_eff(a) on 0 <= z <= 2.5
"""
import json, glob
import numpy as np
from rlds_paired import mcmc as M, likelihood_v2 as L
from scipy.optimize import brentq

LCDM = {"cmb_bao_sn": 1421.359, "cmb_bao_sn_om0.26-0.04": 1422.423,
        "cmb_bao_sn_om0.26-0.08": 1421.627, "cmb_bao_sn_rc0.8": 1421.359,
        "cmb_bao": None}


def derived(p):
    bg = L.rlds_free(p["Om"], p["H0"], p["wb"], p["K"], p["Rc"], p["dR"], p["alpha"])
    a = np.exp(bg.lna)
    Om_early = bg.omega_m_early / bg.h ** 2
    rde = bg.rho_nr - Om_early * a ** -3
    sel = a > 1 / 3.5
    w = -1 - np.gradient(np.log(np.abs(rde[sel])), bg.lna[sel]) / 3
    aa = a[sel]
    zsel = aa >= 1 / 3.5
    A = np.vstack([np.ones(zsel.sum()), 1 - aa[zsel]]).T
    w0, wa = np.linalg.lstsq(A, w[zsel], rcond=None)[0]
    w_at = lambda z: float(np.interp(1 / (1 + z), aa, w))
    return dict(Om_today=bg.Omega_m_today, Om_early=Om_early,
                dOm=bg.Omega_m_today - Om_early, rhoDE_min=float(rde[sel].min()),
                w0_fit=float(w0), wa_fit=float(wa),
                w_z0=w_at(0), w_z05=w_at(0.5), w_z1=w_at(1.0), w_z2=w_at(2.0))


if __name__ == "__main__":
    out = {}
    for f in sorted(glob.glob("validation/profile_K_*.json")):
        tag = f.split("profile_K_")[1][:-5]
        rows = json.load(open(f))
        for r in rows:
            r.update(derived(r))
            if LCDM.get(tag) is not None:
                r["dchi2_vs_lcdm"] = r["chi2"] - LCDM[tag]
        out[tag] = rows
        print(tag)
        for r in rows:
            print("  K=%.3f dchi2=%6.2f Om=%.3f Om_early=%.3f dOm=%+.3f w0=%.3f wa=%+.3f w(z=0,0.5,1,2)=%.3f %.3f %.3f %.3f"
                  % (r["K"], r.get("dchi2_vs_lcdm", np.nan), r["Om_today"], r["Om_early"], r["dOm"],
                     r["w0_fit"], r["wa_fit"], r["w_z0"], r["w_z05"], r["w_z1"], r["w_z2"]))
    json.dump(out, open("validation/derived_observables.json", "w"), indent=1)
