import json, numpy as np, emcee, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
b = emcee.backends.HDFBackend("validation/chains/CPL.h5", read_only=True)
c = b.get_chain(discard=1000, thin=10, flat=True)
H, xe, ye = np.histogram2d(c[:, 3], c[:, 4], bins=40)
Hs = np.sort(H.ravel())[::-1]; cs = np.cumsum(Hs) / Hs.sum(); lev = [Hs[np.searchsorted(cs, q)] for q in (0.95, 0.68)]
prof = json.load(open("validation/derived_observables.json"))["cmb_bao_sn"]
for a, (k, s, t) in zip(ax, [("shooting", 2, "RLDS-shooting"), ("free", 1, "RLDS-free")]):
    r = [x for x in json.load(open(f"validation/scan_w0wa_{k}_{s}.json")) if x["rhoDE_min"] > 0]
    sc = a.scatter([x["w0_fit"] for x in r], [x["wa_fit"] for x in r], c=[x["dOm"] for x in r], s=5, cmap="viridis", vmin=0, vmax=0.25)
    a.contour(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, levels=lev, colors="r")
    a.plot([x["w0_fit"] for x in prof], [x["wa_fit"] for x in prof], "k.-", ms=4, label="profile minima (best fits)")
    a.plot(-1, 0, "k*", ms=12); a.set_title(f"{t}: reachable region (uniform prior scan)")
    a.set_xlabel(r"effective $w_0$"); a.grid(alpha=.3)
ax[0].set_ylabel(r"effective $w_a$"); ax[0].plot([], [], "r-", label="CPL posterior (68/95%)"); ax[0].legend(fontsize=8, loc="lower left")
fig.colorbar(sc, ax=ax, label=r"$\Delta\Omega_m$")
plt.savefig("validation/scan_w0wa.png", dpi=150, bbox_inches="tight")
