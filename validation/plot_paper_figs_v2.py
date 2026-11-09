"""Figs. 3 and 4 of Paper VI, v2: both closures shown and labelled explicitly.
Fig. 3 left : profile likelihood of dOm -- shooting (solid) and free (dashed), without (blue)
              and with (red) the f_gas prior; 95% upper limits with the prior.
Fig. 3 right: effective (w0, wa) of the profile minima of both closures vs the CPL posterior.
Fig. 4      : reachable region from the uniform scan; black = profile minima of the closure
              shown in that panel (shooting left, free right)."""
import json, os, numpy as np, emcee, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from rlds_paired import likelihood_v2 as L
import derived_observables as D

LCDM, LCDM_F = 1421.359, 1422.423
b = emcee.backends.HDFBackend("validation/chains/CPL.h5", read_only=True)
c = b.get_chain(discard=1000, thin=10, flat=True)
H, xe, ye = np.histogram2d(c[:, 3], c[:, 4], bins=40)
Hs = np.sort(H.ravel())[::-1]; cs = np.cumsum(Hs) / Hs.sum(); lev = [Hs[np.searchsorted(cs, q)] for q in (0.95, 0.68)]
X, Y = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])


def good(rows):
    return [r for r in rows if abs(r["dOm"] - r["target"]) < 1e-3]


def bound95(rows, ref):
    """dOm where the profile (relative to its minimum) crosses 4 on the high side."""
    d = np.array([r["dOm"] for r in rows]); y = np.array([r["chi2"] for r in rows]); y = y - y.min()
    i0 = int(np.argmin(y))
    for i in range(i0 + 1, len(d)):
        if y[i] >= 4:
            return d[i - 1] + (4 - y[i - 1]) * (d[i] - d[i - 1]) / (y[i] - y[i - 1])
    return None


def wfit(r, free):
    p = dict(Om=r["Om"], H0=r["H0"], wb=r["wb"], K=r["K"], Rc=r["Rc"], dR=r["dR"], alpha=r["alpha"])
    d = D.derived(p)
    return d["w0_fit"], d["wa_fit"]


sets = {
    "shoot": good(json.load(open("validation/profile_dOm_cmb_bao_sn.json"))),
    "shoot_f": good(json.load(open("validation/profile_dOm_cmb_bao_sn_om0.26-0.04.json"))),
}
for k, f in [("free", "validation/profile_dOm_free_ommax0.6.json"), ("free_f", "validation/profile_dOm_free_om0.26-0.04_ommax0.6.json")]:
    if os.path.exists(f):
        sets[k] = good(json.load(open(f)))

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
style = {"shoot": ("C0", "-", "shooting, no $\\Omega_m$ prior", LCDM), "shoot_f": ("C3", "-", "shooting, $f_{gas}$ prior", LCDM_F),
         "free": ("C0", "--", "free $\\Omega_m$, no $\\Omega_m$ prior", LCDM), "free_f": ("C3", "--", "free $\\Omega_m$, $f_{gas}$ prior", LCDM_F)}
out = {}
for k, rows in sets.items():
    col, ls, lab, ref = style[k]
    rows = sorted(rows, key=lambda r: r["dOm"])
    ax[0].plot([r["dOm"] for r in rows], [r["chi2"] - ref for r in rows], ls, color=col, marker="o" if ls == "-" else "s", ms=4, label=lab)
    if k.endswith("_f"):
        bnd = bound95(rows, ref); out[k] = bnd
        if bnd:
            ax[0].axvline(bnd, color=col, ls=ls, lw=1)
    free = k.startswith("free")
    w = [wfit(r, free) for r in rows]
    ax[1].plot([a for a, _ in w], [b_ for _, b_ in w], ls, color=col, marker="o" if not free else "s", ms=4, label=lab)
kp = json.load(open("validation/derived_observables.json"))["cmb_bao_sn"]
ax[0].plot([r["dOm"] for r in kp], [r["dchi2_vs_lcdm"] for r in kp], "o", color="C0", alpha=.25, ms=5, label="shooting, K-profile points")
ax[0].axhline(0, c="k", lw=.8)
ax[0].set_xlabel(r"$\Delta\Omega_m=\Omega_{m,0}-\omega_m^{\rm early}/h^2$")
ax[0].set_ylabel(r"$\chi^2_{\rm RLDS}-\chi^2_{\Lambda{\rm CDM}}$ (same prior)")
ax[0].legend(fontsize=7.5); ax[0].grid(alpha=.3); ax[0].set_title(r"Profile likelihood of $\Delta\Omega_m$; vertical lines: 95% limits with $f_{gas}$", fontsize=10)
ax[1].contour(X, Y, H.T, levels=lev, colors="k"); ax[1].plot([], [], "k-", label="CPL posterior (68/95%)")
ax[1].plot(-1, 0, "k*", ms=11, label=r"$\Lambda$CDM")
ax[1].set_xlabel(r"effective $w_0$"); ax[1].set_ylabel(r"effective $w_a$"); ax[1].grid(alpha=.3)
ax[1].legend(fontsize=7, loc="lower left"); ax[1].set_title(r"Profile minima of both closures in the $w_0$–$w_a$ plane", fontsize=10)
plt.tight_layout(); plt.savefig("validation/fig3_dOm_w0wa_v2.png", dpi=150)
json.dump(out, open("validation/fig3_bounds.json", "w"), indent=1)
print("95% bounds with f_gas:", out)

# ---- Fig. 4 ----
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.9), sharey=True)
for a, (k, s, t, prof, free) in zip(ax, [("shooting", 2, "RLDS, shooting closure", "shoot", False), ("free", 1, r"RLDS, free $\Omega_m$", "free", True)]):
    r = [x for x in json.load(open(f"validation/scan_w0wa_{k}_{s}.json")) if x["rhoDE_min"] > 0]
    ok = [x for x in r if 0.29 < x["Om_early"] < 0.33]; no = [x for x in r if not 0.29 < x["Om_early"] < 0.33]
    a.scatter([x["w0_fit"] for x in no], [x["wa_fit"] for x in no], s=3, c="0.82", label=r"other $\Omega_e$")
    sc = a.scatter([x["w0_fit"] for x in ok], [x["wa_fit"] for x in ok], c=[x["dOm"] for x in ok], s=9, cmap="viridis", vmin=0, vmax=0.15,
                   label=r"$0.29<\Omega_e<0.33$ (colour: $\Delta\Omega_m$)")
    a.contour(X, Y, H.T, levels=lev, colors="r")
    if prof in sets:
        rows = sorted(sets[prof], key=lambda q: q["dOm"]); w = [wfit(q, free) for q in rows]
        a.plot([u for u, _ in w], [v for _, v in w], "k.-", ms=4, label=("shooting" if not free else "free-$\\Omega_m$") + r" $\Delta\Omega_m$-profile minima")
    a.plot(-1, 0, "k*", ms=12); a.set_title(f"{t}\n{len(r)} solutions, {len(ok)} with $0.29<\\Omega_e<0.33$", fontsize=10)
    a.set_xlabel(r"effective $w_0$ (CPL fit, $0\leq z\leq2.5$)"); a.grid(alpha=.3); a.set_xlim(-1.2, -0.3)
    a.legend(fontsize=6.5, loc="lower left")
ax[0].set_ylabel(r"effective $w_a$"); ax[0].plot([], [], "r-", label="CPL posterior, 68/95%")
fig.colorbar(sc, ax=ax, label=r"$\Delta\Omega_m$", shrink=0.9)
plt.savefig("validation/fig4_scan_w0wa_v2.png", dpi=150, bbox_inches="tight")
