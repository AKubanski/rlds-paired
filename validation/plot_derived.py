"""Diagnostic working plot, not a paper figure: Delta chi2 against dOm and the effective
(w0, wa) of all K-profile points. Paper figures: plot_profile_K.py (Fig. 1),
plot_paper_figs.py (Fig. 2), plot_paper_figs_v2.py (Figs. 3-4)."""
import json, numpy as np, emcee, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = json.load(open("validation/derived_observables.json"))
lab = {"cmb_bao_sn": "no Ωm prior (Rc ≤ 0.6)", "cmb_bao_sn_rc0.8": "no Ωm prior, Rc ≤ 0.8",
       "cmb_bao_sn_om0.26-0.08": "f_gas: Ωm = 0.26 ± 0.08", "cmb_bao_sn_om0.26-0.04": "f_gas: Ωm = 0.26 ± 0.04"}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for i, (t, l) in enumerate(lab.items()):
    r = d[t]
    ax[0].plot([x["dOm"] for x in r], [x["dchi2_vs_lcdm"] for x in r], "o", c=f"C{i}", label=l)
    ax[1].plot([x["w0_fit"] for x in r], [x["wa_fit"] for x in r], "o", c=f"C{i}", ms=4, label="RLDS: " + l)
ax[0].axhline(0, c="k", lw=.8); ax[0].set_xlabel("ΔΩm = Ωm(today) − ωm(early)/h²  (matter created by the coupling)")
ax[0].set_ylabel("χ²(RLDS) − χ²(ΛCDM, same prior)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[0].set_title("Δχ² against ΔΩm: lower envelope = profile in ΔΩm")
b = emcee.backends.HDFBackend("validation/chains/CPL.h5", read_only=True)
c = b.get_chain(discard=1000, thin=10, flat=True)
H, xe, ye = np.histogram2d(c[:, 3], c[:, 4], bins=40)
Hs = np.sort(H.ravel())[::-1]; cs = np.cumsum(Hs) / Hs.sum()
lev = [Hs[np.searchsorted(cs, q)] for q in (0.95, 0.68)]
ax[1].contour(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, levels=lev, colors="k")
ax[1].plot([], [], "k-", label="CPL, our MCMC (68/95%)")
ax[1].plot(-1, 0, "k*", ms=12, label="ΛCDM")
ax[1].set_xlabel("w0 (effective)"); ax[1].set_ylabel("wa (effective)"); ax[1].grid(alpha=.3)
ax[1].legend(fontsize=7, loc="lower left"); ax[1].set_title("What RLDS can mimic in the w0–wa plane")
plt.tight_layout(); plt.savefig("validation/derived_observables.png", dpi=150)
