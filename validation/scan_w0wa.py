"""
Reachable region of RLDS in the effective (w0, wa) plane -- not only the
profile minima. Uniform random samples over the WIDE priors
(Om 0.2-0.45, K 0-0.25, Rc 0.3-0.9, dR 0.01-0.2; alpha 2.5-3.5 for the free
model, alpha from the shooting condition for the shooting model).
w_eff is defined as in derived_observables.py. H0 and omega_b are fixed
(68, 0.0225): w_eff depends on them only through radiation, negligible at z<2.5.
No likelihood is evaluated: this maps what the model CAN produce.
Usage: scan_w0wa.py free|shooting N SEED
"""
import json, sys
import numpy as np
from rlds_paired import likelihood_v2 as L
import profile_K as P
from derived_observables import derived

kind, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rng = np.random.default_rng(seed)
out = []
for i in range(N):
    Om, K = rng.uniform(0.2, 0.45), rng.uniform(0, 0.25)
    Rc, dR = rng.uniform(0.3, 0.9), rng.uniform(0.01, 0.2)
    if kind == "free":
        a = rng.uniform(2.5, 3.5)
    else:
        a = P.alpha_from(K, Om, Rc, dR)
        if a is None:
            continue
    try:
        d = derived(dict(Om=Om, H0=68.0, wb=0.0225, K=K, Rc=Rc, dR=dR, alpha=a))
    except Exception:
        continue
    if not np.isfinite(d["w0_fit"]):
        continue
    out.append(dict(Om=Om, K=K, Rc=Rc, dR=dR, alpha=a, **d))
json.dump(out, open(f"validation/scan_w0wa_{kind}_{seed}.json", "w"))
print(kind, len(out), "valid of", N)
