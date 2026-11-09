"""English figures for Paper VI (Figs. 2 and 3)."""
import json, numpy as np, emcee, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
LCDM = {"cmb_bao_sn": 1421.359, "cmb_bao_sn_om0.26-0.04": 1422.423,
        "cmb_bao_sn_om0.26-0.08": 1421.627, "cmb_bao_sn_rc0.8": 1421.359}
lab = {"cmb_bao_sn": r"no $\Omega_m$ prior ($R_c\leq0.6$)", "cmb_bao_sn_rc0.8": r"no $\Omega_m$ prior, $R_c\leq0.8$",
       "cmb_bao_sn_om0.26-0.08": r"$f_{gas}$: $\Omega_m=0.26\pm0.08$", "cmb_bao_sn_om0.26-0.04": r"$f_{gas}$: $\Omega_m=0.26\pm0.04$"}
# Fig 2
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for i, (t, l) in enumerate(lab.items()):
    r = json.load(open(f"validation/profile_K_{t}.json")); K = [x["K"] for x in r]
    ax[0].plot(K, [x["chi2"] - LCDM[t] for x in r], "o-", c=f"C{i}", label=l)
    ax[1].plot(K, [x["Om"] for x in r], "o-", c=f"C{i}")
    ax[2].plot(K, [x["Rc"] for x in r], "o-", c=f"C{i}")
ax[0].axhline(0, c="k", lw=.8); ax[0].set_ylabel(r"$\chi^2_{\rm RLDS}-\chi^2_{\Lambda{\rm CDM}}$ (same prior)")
ax[1].set_ylabel(r"$\Omega_{m,0}$ at profile minimum"); ax[2].set_ylabel(r"$R_c$ at profile minimum")
for a in ax: a.set_xlabel("K"); a.grid(alpha=.3)
ax[0].legend(fontsize=8); ax[0].set_title("Profile of K, CMB+BAO+SN")
ax[1].set_title(r"K – $\Omega_m$ degeneracy"); ax[2].set_title(r"K – $R_c$ degeneracy")
plt.tight_layout(); plt.savefig("validation/fig2_profile_K_priors.png", dpi=150)
# Fig 3
d = json.load(open("validation/derived_observables.json"))
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for t, l, c in [("cmb_bao_sn", r"no $\Omega_m$ prior", "C0"), ("cmb_bao_sn_om0.26-0.04", r"$f_{gas}$: $\Omega_m=0.26\pm0.04$", "C3")]:
    ax[0].plot([x["dOm"] for x in d[t]], [x["dchi2_vs_lcdm"] for x in d[t]], "o", c=c, alpha=.25, ms=4)
    p = [x for x in json.load(open(f"validation/profile_dOm_{t}.json")) if abs(x["dOm"] - x["target"]) < 1e-3]
    ax[0].plot([x["dOm"] for x in p], [x["chi2"] - LCDM[t] for x in p], "-s", c=c, label=l + " (profile)")
ax[0].axhline(0, c="k", lw=.8); ax[0].axvline(0.0733, c="C3", ls="--", lw=1)
ax[0].text(0.0745, 3.2, "95% upper limit\n(f_gas prior)", color="C3", fontsize=8)
ax[0].set_xlabel(r"$\Delta\Omega_m=\Omega_{m,0}-\omega_m^{\rm early}/h^2$")
ax[0].set_ylabel(r"$\chi^2_{\rm RLDS}-\chi^2_{\Lambda{\rm CDM}}$ (same prior)")
ax[0].legend(fontsize=8, loc="lower left"); ax[0].grid(alpha=.3)
ax[0].set_title(r"Profile likelihood of $\Delta\Omega_m$ (dots: K-profile points)")
for i, t in enumerate(lab):
    ax[1].plot([x["w0_fit"] for x in d[t]], [x["wa_fit"] for x in d[t]], "o", c=f"C{i}", ms=4, label="RLDS, " + lab[t])
b = emcee.backends.HDFBackend("validation/chains/CPL.h5", read_only=True)
c = b.get_chain(discard=1000, thin=10, flat=True)
H, xe, ye = np.histogram2d(c[:, 3], c[:, 4], bins=40)
Hs = np.sort(H.ravel())[::-1]; cs = np.cumsum(Hs) / Hs.sum()
lev = [Hs[np.searchsorted(cs, q)] for q in (0.95, 0.68)]
ax[1].contour(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, levels=lev, colors="k")
ax[1].plot([], [], "k-", label="CPL posterior (68/95%)"); ax[1].plot(-1, 0, "k*", ms=12, label=r"$\Lambda$CDM")
ax[1].set_xlabel(r"effective $w_0$"); ax[1].set_ylabel(r"effective $w_a$"); ax[1].grid(alpha=.3)
ax[1].legend(fontsize=7, loc="lower left"); ax[1].set_title(r"RLDS solutions in the $w_0$–$w_a$ plane")
plt.tight_layout(); plt.savefig("validation/fig3_dOm_w0wa.png", dpi=150)
