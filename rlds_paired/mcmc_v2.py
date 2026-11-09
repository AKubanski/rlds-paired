"""
rlds_paired.mcmc_v2
===================
Resumable emcee runs on the corrected likelihood (likelihood_v2).
Priors are the pre-registered ones from the technical specification, section 5.1.
"""
from __future__ import annotations
import numpy as np
from . import likelihood_v2 as L

MODELS = {
    "LCDM": (L.lcdm, ["Om", "H0", "wb"],
             [(0.2, 0.45), (60, 75), (0.019, 0.026)]),
    "wCDM": (L.wcdm, ["Om", "H0", "wb", "w"],
             [(0.2, 0.45), (60, 75), (0.019, 0.026), (-1.5, -0.5)]),
    "CPL": (L.cpl, ["Om", "H0", "wb", "w0", "wa"],
            [(0.2, 0.45), (60, 75), (0.019, 0.026), (-2, 0), (-3, 2)]),
    "RLDS-free": (L.rlds_free, ["Om", "H0", "wb", "K", "Rc", "dR", "alpha"],
                  [(0.2, 0.45), (60, 75), (0.019, 0.026), (0, 0.25),
                   (0.3, 0.6), (0.03, 0.2), (2.5, 3.5)]),
}

_SN = None


def log_prob(theta, model, use):
    global _SN
    fn, _, bounds = MODELS[model]
    for v, (lo, hi) in zip(theta, bounds):
        if not lo <= v <= hi:
            return -np.inf
    if model == "CPL" and theta[3] + theta[4] >= 0:      # early matter domination
        return -np.inf
    if "sn" in use and _SN is None:
        _SN = L.load_sn()
    try:
        c = L.chi2_total(fn(*theta), _SN, use)
    except Exception:
        return -np.inf
    return -0.5 * c if np.isfinite(c) else -np.inf


def run(model, use, start, scatter, backend_file, n_steps, n_walkers=32,
        seed=0, resume=True, processes=2, progress=False):
    import emcee
    from multiprocessing import Pool
    ndim = len(MODELS[model][1])
    import os
    exists = os.path.exists(backend_file)
    backend = emcee.backends.HDFBackend(backend_file)
    fresh = not (resume and exists and backend.iteration > 0)
    if fresh:
        backend.reset(n_walkers, ndim)
        rng = np.random.default_rng(seed)
        lo = np.array([b[0] for b in MODELS[model][2]])
        hi = np.array([b[1] for b in MODELS[model][2]])
        p0 = np.clip(np.asarray(start) + np.asarray(scatter) * rng.standard_normal((n_walkers, ndim)),
                     lo + 1e-6, hi - 1e-6)
        if model == "CPL":
            bad = p0[:, 3] + p0[:, 4] >= 0
            p0[bad, 4] = -p0[bad, 3] - 0.1
    else:
        p0 = None
    with Pool(processes) as pool:
        s = emcee.EnsembleSampler(n_walkers, ndim, log_prob, args=(model, use),
                                  pool=pool, backend=backend)
        s.run_mcmc(p0, n_steps, progress=progress)
    return backend
