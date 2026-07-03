"""
Validated reference results for the RLDS canonical benchmark (Papers III-V).

These numbers are **not stored as literals**. Each is recomputed from the model
code on demand via the documented protocol below, so the package can never
silently ship a stale headline number. (Several previously hardcoded values did
not survive recomputation against the corrected background; this module exists
so that cannot recur.)

Usage
-----
>>> from rlds_paired.reference import reference_results
>>> r = reference_results()          # cached: one k-scan + one ODE pass
>>> print(r.summary())

or a standalone report::

    python -m rlds_paired.reference

Canonical protocol
------------------
Parameters (``rlds_paired.constants``): ``alpha=ALPHA``, ``Rc=RC``, ``dR=DR``,
``Omega_m=OMEGA_M``; canonical coupling ``K=K_BEST``.

Paired coefficient ``c``
    Slope of  ``Delta(f*sigma8)/f*sigma8 = c * K``  from :func:`analysis.k_scan`
    on the **full coupled-ODE background**, with the baseline taken at ``K=0``
    on the *same* ODE path (this removes the ~1% ODE-vs-analytic offset),
    ``z_eval = 0.5``, K-scan over ``[0.01, 0.40]``. The sign is **negative**:
    the dark-sector energy exchange ``Q > 0`` drains the scalar, so structure
    growth is *suppressed* (an algebraic-background model gives the wrong,
    positive sign -- see ``core/growth.py``).

Spectral gap (Fenichel scale separation)
    ``delta_max = max_{z < 3} |dR/dN| / 3``  from
    :func:`core.spectral.compute_spectral_gap` at ``K = K_BEST``. The Fenichel
    persistence bound (0.10) is a **criterion that is reported pass/fail -- it
    is not assumed to hold**. At the canonical benchmark it is *violated*
    (``delta_max ~ 0.25``), with the peak in the transition window
    (``z ~ 0.42``, i.e. the observable epoch sampled by DESI/Euclid f*sigma8).
    The slow-manifold reduction therefore supplies a *marginal* (~25%)
    correction there rather than an asymptotically exact one.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from rlds_paired.constants import ALPHA, DR, FENICHEL_BOUND, K_BEST, OMEGA_M, RC


@dataclass(frozen=True)
class ReferenceResults:
    """Computed canonical headline results (see module docstring for protocol)."""

    # paired signature
    paired_coefficient: float
    paired_coefficient_err: float
    paired_r2: float
    # spectral gap
    delta_max: float
    delta_max_z: float
    fenichel_bound: float
    fenichel_satisfied: bool
    # protocol echo (so a printed report is self-describing)
    K: float
    alpha: float
    Rc: float
    dR: float
    Omega_m: float

    def summary(self) -> str:
        fen = "SATISFIED" if self.fenichel_satisfied else "VIOLATED"
        return (
            "RLDS canonical reference (recomputed from code)\n"
            f"  parameters : K={self.K}  alpha={self.alpha}  "
            f"Rc={self.Rc}  dR={self.dR}  Omega_m={self.Omega_m}\n"
            f"  paired c   : {self.paired_coefficient:+.4f} "
            f"+/- {self.paired_coefficient_err:.4f}   (r2 = {self.paired_r2:.4f})\n"
            f"             : sign negative => growth suppression\n"
            f"  delta_max  : {self.delta_max:.4f}  at z = {self.delta_max_z:.3f} "
            f"(transition / observable epoch)\n"
            f"  Fenichel   : bound delta < {self.fenichel_bound}: {fen} "
            f"(reduction is marginal, ~{round(100*self.delta_max)}% in the transition)\n"
        )


@lru_cache(maxsize=1)
def reference_results() -> ReferenceResults:
    """Recompute the canonical headline results. Cached (one pass per process)."""
    # Imported lazily so importing this module stays cheap and avoids cycles.
    from rlds_paired.analysis import k_scan
    from rlds_paired.core.spectral import compute_spectral_gap

    ks = k_scan(verbose=False)  # canonical protocol (defaults)
    sg = compute_spectral_gap(K=K_BEST)  # canonical coupling

    return ReferenceResults(
        paired_coefficient=float(ks.coefficient),
        paired_coefficient_err=float(ks.std_err),
        paired_r2=float(ks.r2),
        delta_max=float(sg.delta_max),
        delta_max_z=float(sg.z_at_max),
        fenichel_bound=float(FENICHEL_BOUND),
        fenichel_satisfied=bool(sg.fenichel_satisfied),
        K=float(K_BEST),
        alpha=float(ALPHA),
        Rc=float(RC),
        dR=float(DR),
        Omega_m=float(OMEGA_M),
    )


if __name__ == "__main__":
    print(reference_results().summary())
