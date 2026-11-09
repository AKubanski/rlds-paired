"""Diagnostic working plot, not a paper figure: K profiles with the f_gas Omega_m priors
and with R_c <= 0.8, and the K-Omega_m and K-R_c degeneracies. The paper's Fig. 2 is
made by plot_paper_figs.py."""
import json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
LCDM = {"cmb_bao_sn": 1421.359, "cmb_bao_sn_om0.26-0.04": 1422.423,
        "cmb_bao_sn_om0.26-0.08": 1421.627, "cmb_bao_sn_rc0.8": 1421.359}
lab = {"cmb_bao_sn": "no Ωm prior (Rc ≤ 0.6)", "cmb_bao_sn_rc0.8": "no Ωm prior, Rc ≤ 0.8",
       "cmb_bao_sn_om0.26-0.08": "f_gas: Ωm = 0.26 ± 0.08", "cmb_bao_sn_om0.26-0.04": "f_gas: Ωm = 0.26 ± 0.04"}
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for i, (t, l) in enumerate(lab.items()):
    r = json.load(open(f"validation/profile_K_{t}.json")); K = [x["K"] for x in r]
    ax[0].plot(K, [x["chi2"] - LCDM[t] for x in r], "o-", c=f"C{i}", label=l)
    ax[1].plot(K, [x["Om"] for x in r], "o-", c=f"C{i}", label=l)
    ax[2].plot(K, [x["Rc"] for x in r], "o-", c=f"C{i}", label=l)
ax[0].axhline(0, c="k", lw=.8); ax[0].set_ylabel("χ²(RLDS) − χ²(ΛCDM, same prior)")
ax[1].set_ylabel("Ωm today at the minimum"); ax[2].set_ylabel("Rc at the minimum")
for a in ax: a.set_xlabel("K"); a.grid(alpha=.3)
ax[0].legend(fontsize=8); ax[0].set_title("Profile of K, CMB+BAO+SN")
ax[1].set_title("K–Ωm degeneracy"); ax[2].set_title("K–Rc degeneracy")
plt.tight_layout(); plt.savefig("validation/profile_K_priors.png", dpi=150)
