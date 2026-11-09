"""
rlds-paired: RLDS Paired Dark-Sector Signature (Papers III–V)
=============================================================

Python library for reproducing the numerical results of the RLDS
cosmology series (Papers III–V) by Aleksander Kubanski.

The central result is the **paired growth–expansion signature**:

    Δ(f·σ₈) / (f·σ₈)_ΛCDM  =  c · K + O(K²)
    c  ≈  −0.19   (negative: growth suppression; computed, not stored —
                   see ``rlds_paired.reference.reference_results``)

``rlds-paired`` is the observational–cosmological branch of the RLDS SDK
family.  The mathematical core is in ``rlds-mat``; the microphysical
descent (Paper V) is validated here via the spectral-gap module.

Based on the research published in:

    Kubanski, A. (2026). *Paired observational signature of the RLDS
    dark sector: growth suppression and phantom excursion.*
    Zenodo. DOI: 10.5281/zenodo.19636371  [Paper III]

    Kubanski, A. (2026). *Cosmological consequences under RLDS closure.*
    Zenodo. DOI: 10.5281/zenodo.19640061  [Paper IV]

    Kubanski, A. (2026). *Microphysical descent to RLDS.*
    Zenodo. [Paper V]

Quick start
-----------
>>> import rlds_paired
>>>
>>> # Spectral gap validation (Paper V, Appendix C)
>>> sg = rlds_paired.compute_spectral_gap()
>>> print(sg.summary())
>>>
>>> # ODE trajectory (Paper IV)
>>> traj = rlds_paired.RLDSTrajectory()
>>> print(traj.R_star, traj.E(0.5))
>>>
>>> # All canonical benchmark values
>>> print(rlds_paired.CANONICAL)

Import cost
-----------
Importing ``rlds_paired`` is lightweight: only constants and exceptions are
loaded eagerly.  The physics modules (numpy/scipy) and the analysis module
(matplotlib) are imported lazily on first use via PEP 562 ``__getattr__``.

Homepage  : https://research.kubanski.pro
Repository: https://github.com/AKubanski/rlds-paired
Author    : Aleksander Kubanski <research@kubanski.pro>
License   : MIT
"""

from __future__ import annotations

from ._version import __version__
from .exceptions import (
    RLDSPairedError,
    ConvergenceError,
    PhysicsError,
    ConfigurationError,
)

__author__ = "Aleksander Kubanski"
__email__ = "research@kubanski.pro"
__license__ = "MIT"

# ---------------------------------------------------------------------------
# Eager imports — constants only (no heavy deps)
# ---------------------------------------------------------------------------
from .constants import (
    K_BEST,
    ALPHA,
    RC,
    DR,
    EPSILON,
    OMEGA_M,
    OMEGA_R,
    R_STAR,
    W_DE,
    Z_ACC,
    FENICHEL_BOUND,
    SIGMA8,
    Z_INIT,
)

# ---------------------------------------------------------------------------
# Lazy attribute access (PEP 562) for scipy/matplotlib-heavy modules.
#
# ``RLDSTrajectory``, ``compute_spectral_gap``, ``k_scan``, etc.
# are only imported when first accessed as ``rlds_paired.<Name>``.
# ---------------------------------------------------------------------------
_LAZY: dict[str, tuple[str, str]] = {
    # core.algebraic
    "F_transition": ("core.algebraic", "F_transition"),
    "solve_R": ("core.algebraic", "solve_R"),
    "E_squared_alg": ("core.algebraic", "E_squared_alg"),
    "E_alg": ("core.algebraic", "E_alg"),
    "w_eff_alg": ("core.algebraic", "w_eff_alg"),
    # core.ode
    "RLDSTrajectory": ("core.ode", "RLDSTrajectory"),
    "flow": ("core.ode", "flow"),
    "psi": ("core.ode", "psi"),
    "switching": ("core.ode", "switching"),
    # core.growth
    "growth_factor": ("core.growth", "growth_factor"),
    # core.spectral
    "SpectralGapResult": ("core.spectral", "SpectralGapResult"),
    "compute_spectral_gap": ("core.spectral", "compute_spectral_gap"),
    "spectral_gap_table": ("core.spectral", "spectral_gap_table"),
    # analysis
    "KScanResult": ("analysis", "KScanResult"),
    "k_scan": ("analysis", "k_scan"),
    "FSigma8Result": ("analysis", "FSigma8Result"),
    "fsigma8_comparison": ("analysis", "fsigma8_comparison"),
    # reference (recomputed headline results)
    "reference_results": ("reference", "reference_results"),
    "ReferenceResults": ("reference", "ReferenceResults"),
    # mcmc (shooting MCMC; run_mcmc needs the optional 'emcee' extra)
    "ShootingRLDS": ("mcmc", "ShootingRLDS"),
    "run_mcmc": ("mcmc", "run_mcmc"),
    "log_posterior": ("mcmc", "log_posterior"),
    "get_selfconsistent_model": ("mcmc", "get_selfconsistent_model"),
    "reproduce_reference_chi2": ("mcmc", "reproduce_reference_chi2"),
}


def __getattr__(name: str):
    """PEP 562 lazy loader for heavy submodule objects."""
    target = _LAZY.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_path, attr_name = target
    import importlib

    module = importlib.import_module(f".{module_path}", __name__)
    value = getattr(module, attr_name)
    # Cache on the module so subsequent lookups skip __getattr__.
    globals()[name] = value
    return value


def __dir__():
    """Include lazily-loaded names in dir(rlds_paired)."""
    return sorted(list(globals().keys()) + list(_LAZY.keys()))


# ---------------------------------------------------------------------------
# Convenience dict of canonical benchmark *inputs* (closed-form only).
# Computed results (paired coefficient, delta_max) are not here by design;
# get them from rlds_paired.reference.reference_results().
# ---------------------------------------------------------------------------
CANONICAL: dict = {
    "K": K_BEST,
    "alpha": ALPHA,
    "Rc": RC,
    "dR": DR,
    "Omega_m": OMEGA_M,
    "R_star": R_STAR,
    "w_DE": W_DE,
    "z_acc": Z_ACC,
    "fenichel": FENICHEL_BOUND,
}
