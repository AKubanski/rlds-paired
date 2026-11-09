"""Summarise HDF5 chains: tau, convergence, posterior mean ± 68 %, 95 % upper limits."""
import sys, numpy as np, emcee
from rlds_paired.mcmc_v2 import MODELS
for model in sys.argv[1:]:
    b = emcee.backends.HDFBackend(f"validation/chains/{model}.h5", read_only=True)
    n = b.iteration
    tau = b.get_autocorr_time(tol=0)
    burn = int(5 * np.max(tau)); thin = max(1, int(np.min(tau) / 2))
    flat = b.get_chain(discard=burn, thin=thin, flat=True)
    print(f"\n{model}: steps={n}, max tau={np.max(tau):.0f}, steps/tau={n/np.max(tau):.0f} "
          f"({'OK' if n/np.max(tau) >= 50 else 'NOT converged by 50-tau rule'}), "
          f"acc={np.mean(b.accepted/n):.2f}, samples={len(flat)}")
    for i, p in enumerate(MODELS[model][1]):
        x = flat[:, i]; lo, md, hi = np.percentile(x, [16, 50, 84])
        print(f"  {p:6s} mean={x.mean():.4f}  68%: [{lo:.4f}, {hi:.4f}]  95% upper={np.percentile(x,95):.4f}  tau={tau[i]:.0f}")
