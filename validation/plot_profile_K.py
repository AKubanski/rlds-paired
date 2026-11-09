import json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for f, lab, c in [("cmb_bao_sn", "CMB+BAO+SN", "C0"), ("cmb_bao", "CMB+BAO", "C1")]:
    r = json.load(open(f"validation/profile_K_{f}.json")); m = min(x["chi2"] for x in r)
    K = [x["K"] for x in r]
    ax[0].plot(K, [x["chi2"] - m for x in r], "o-", c=c, label=lab)
    ax[1].plot(K, [x["Om"] for x in r], "o-", c=c, label=lab)
ax[0].axhline(4, ls="--", c="gray"); ax[0].set_xlabel("K"); ax[0].set_ylabel(r"$\Delta\chi^2$ (profile)")
ax[1].axhline(0.45, ls=":", c="gray"); ax[1].set_xlabel("K"); ax[1].set_ylabel(r"$\Omega_m$ today at profile minimum")
for a in ax: a.legend(); a.grid(alpha=.3)
ax[0].set_title("RLDS-shooting: profile of K"); ax[1].set_title("Degeneracy K – $\Omega_m$")
plt.tight_layout(); plt.savefig("validation/profile_K.png", dpi=150)
