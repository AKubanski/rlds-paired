"""
Canonical benchmark *inputs* for the RLDS series (Papers I–V).

This module holds only **input parameters** and **closed-form** quantities
(flow parameters, densities, the analytic attractor R* = sqrt(K), the bare
equation of state w = -alpha/3, fiducial sigma_8, etc.).

It deliberately does **not** store computed headline *results* (the paired
coefficient, its error / r^2, or the spectral-gap delta_max). Those are
recomputed from the model code on demand in :mod:`rlds_paired.reference`
(see :func:`rlds_paired.reference.reference_results`), so a stale literal can
never silently ship. Several previously hardcoded results did not survive
recomputation against the corrected coupled-ODE background.

References
----------
Paper III : https://doi.org/10.5281/zenodo.19636371
Paper IV  : https://doi.org/10.5281/zenodo.19640061
Paper V   : Zenodo (Aleksander Kubanski, 2026)
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# RLDS flow parameters
# ---------------------------------------------------------------------------

#: Coupling constant (DESI DR2 shooting posterior median)
K_BEST: float = 0.0106

#: Dilution rate  alpha = 3|w_DE|  (DR2 median; railed at the prior floor 3.01)
ALPHA: float = 3.0102

#: Switching threshold (matter fraction; DR2 shooting median)
RC: float = 0.565

#: Switching width (DR2 shooting median; canonical)
DR: float = 0.170

#: ODE regularisation epsilon
EPSILON: float = 1e-12

# ---------------------------------------------------------------------------
# Derived quantities (Paper IV)
# ---------------------------------------------------------------------------

#: Matter density (RLDS shooting self-consistency; DR2)
OMEGA_M: float = 0.297

#: Radiation density (photons + neutrinos)
OMEGA_R: float = 9.15e-5

#: Attractor position R* = sqrt(K) at canonical parameters
R_STAR: float = 0.1029

#: Redshift of cosmic acceleration onset, defined as R(z_acc) = 2/3 on the
#: DR2 ODE background (same threshold formula as Paper IV and mcmc.find_z_acc).
#: Computed via brentq on RLDSTrajectory.R at canonical K=0.0106 -> 0.7144,
#: rounded to 0.714.  Documentation-only; not used in computations or tests.
Z_ACC: float = 0.714

#: Bare dark-energy equation of state  w = -alpha/3
W_DE: float = -1.003

#: Integration start redshift
Z_INIT: float = 1100.0

# ---------------------------------------------------------------------------
# Paired signature / spectral gap
# ---------------------------------------------------------------------------
# NOTE: the paired coefficient (slope of Delta(f*sigma8)/f*sigma8 = c * K, its
# std-err and r^2) and the spectral-gap delta_max are RESULTS, not inputs, and
# are intentionally not stored here. Import them from rlds_paired.reference:
#     from rlds_paired.reference import reference_results
#     r = reference_results()      # r.paired_coefficient, r.delta_max, ...
# Only the Fenichel *criterion* lives here, as a fixed threshold to test against.

#: Fenichel persistence bound (a criterion delta_max is tested against, not a result)
FENICHEL_BOUND: float = 0.10

# ---------------------------------------------------------------------------
# Observational / fiducial values
# ---------------------------------------------------------------------------

#: Fiducial sigma_8 at z=0
SIGMA8: float = 0.811

#: Planck 2018 matter density (external reference)
OMEGA_M_PLANCK: float = 0.315
