"""Usage: run_mcmc.py MODEL N_STEPS   (data: CMB+BAO+SN; resumable)"""
import sys, json, numpy as np, emcee
from rlds_paired import mcmc_v2 as MC
model, n = sys.argv[1], int(sys.argv[2])
START = {"LCDM": ([0.301, 68.4, 0.02253], [0.003, 0.3, 0.0001]),
         "wCDM": ([0.305, 67.9, 0.02257, -0.977], [0.003, 0.3, 0.0001, 0.01]),
         "CPL": ([0.311, 67.6, 0.02245, -0.85, -0.52], [0.003, 0.3, 0.0001, 0.03, 0.1]),
         "RLDS-free": ([0.33, 68.0, 0.0225, 0.04, 0.45, 0.11, 3.03], [0.02, 0.3, 0.0001, 0.02, 0.08, 0.04, 0.03])}
f = f"validation/chains/{model}.h5"
import os; os.makedirs("validation/chains", exist_ok=True)
b = MC.run(model, ("cmb", "bao", "sn"), *START[model], f, n, seed=20261104)
ch = b.get_chain()
try: tau = b.get_autocorr_time(tol=0)
except Exception: tau = np.full(ch.shape[2], np.nan)
print(json.dumps(dict(model=model, steps=int(b.iteration), tau=[round(float(t), 1) for t in tau],
      steps_over_tau=round(float(b.iteration / np.nanmax(tau)), 1),
      acc=round(float(np.mean(b.accepted / b.iteration)), 3))))
