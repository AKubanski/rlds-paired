"""Fig. 4 of Paper VI: reachable effective (w0, wa) region of RLDS (uniform scan of the wide priors)."""
import json, numpy as np, emcee, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
b = emcee.backends.HDFBackend("validation/chains/CPL.h5", read_only=True)
c = b.get_chain(discard=1000, thin=10, flat=True)
H, xe, ye = np.histogram2d(c[:, 3], c[:, 4], bins=40)
Hs = np.sort(H.ravel())[::-1]; cs = np.cumsum(Hs) / Hs.sum(); lev = [Hs[np.searchsorted(cs, q)] for q in (0.95, 0.68)]
prof = json.load(open("validation/derived_observables.json"))["cmb_bao_sn"]
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.9), sharey=True)
for a, (k, s, t) in zip(ax, [("shooting", 2, "RLDS, shooting closure"), ("free", 1, r"RLDS, free $\Omega_m$")]):
    r = [x for x in json.load(open(f"validation/scan_w0wa_{k}_{s}.json")) if x["rhoDE_min"] > 0]
    ok = [x for x in r if 0.29 < x["Om_early"] < 0.33]; no = [x for x in r if not 0.29 < x["Om_early"] < 0.33]
    a.scatter([x["w0_fit"] for x in no], [x["wa_fit"] for x in no], s=3, c="0.82", label=r"other $\Omega_e$")
    sc = a.scatter([x["w0_fit"] for x in ok], [x["wa_fit"] for x in ok], c=[x["dOm"] for x in ok], s=9, cmap="viridis", vmin=0, vmax=0.15,
                   label=r"$0.29<\Omega_e<0.33$ (coloured by $\Delta\Omega_m$)")
    a.contour(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), H.T, levels=lev, colors="r")
    a.plot([x["w0_fit"] for x in prof], [x["wa_fit"] for x in prof], "k.-", ms=4, label="profile minima")
    a.plot(-1, 0, "k*", ms=12); a.set_title(f"{t}\n{len(r)} solutions, {len(ok)} with $0.29<\\Omega_e<0.33$", fontsize=10)
    a.set_xlabel(r"effective $w_0$ (CPL fit, $0\leq z\leq2.5$)"); a.grid(alpha=.3); a.set_xlim(-1.2, -0.3)
ax[0].set_ylabel(r"effective $w_a$"); ax[0].plot([], [], "r-", label="CPL posterior, 68/95%")
ax[0].legend(fontsize=7, loc="lower left")
fig.colorbar(sc, ax=ax, label=r"$\Delta\Omega_m$", shrink=0.9)
plt.savefig("validation/fig4_scan_w0wa.png", dpi=150, bbox_inches="tight")
