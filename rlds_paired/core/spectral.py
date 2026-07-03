"""
Spectral gap validation for the RLDS slow-manifold reduction (Paper V, Appendix C).

Computes the scale-separation parameter:
    delta(N) = |dR/dN| / 3

The Fenichel persistence theorem would require delta_max < 0.10. This is treated
here as a **criterion that is tested and reported pass/fail, not assumed**. At
the canonical benchmark (K ~ 0.0106) it is VIOLATED: delta_max ~ 0.244, with the
peak near z ~ 0.40 (the transition / observable epoch). The slow-manifold
reduction therefore supplies a *marginal* (~25%) correction there rather than an
asymptotically exact one. The computed value is not stored; see
:func:`rlds_paired.reference.reference_results`.

Reference: Paper V, Section 5.1 and Appendix C.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rlds_paired.constants import (
    ALPHA,
    DR,
    EPSILON,
    FENICHEL_BOUND,
    K_BEST,
    OMEGA_M,
    RC,
)
from rlds_paired.core.ode import RLDSTrajectory, flow


@dataclass
class SpectralGapResult:
    """Result of a spectral gap computation.

    Attributes
    ----------
    K:
        Coupling constant.
    delta_max:
        Maximum delta(N) = |dR/dN|/3 over z < 3.
    delta_mean:
        Mean delta(N) over z < 3.
    z_at_max:
        Redshift at which delta_max is attained.
    fenichel_satisfied:
        True if delta_max < FENICHEL_BOUND (0.10).
    N_arr, R_arr, delta_arr:
        Full trajectory arrays.
    """

    K: float
    delta_max: float
    delta_mean: float
    z_at_max: float
    fenichel_satisfied: bool
    N_arr: np.ndarray
    R_arr: np.ndarray
    delta_arr: np.ndarray

    def summary(self) -> str:
        status = "SATISFIED ✓" if self.fenichel_satisfied else "VIOLATED ✗"
        return (
            f"Spectral gap  (K = {self.K:.3f})\n"
            f"  delta_max  = {self.delta_max:.4f}  at z = {self.z_at_max:.3f}\n"
            f"  delta_mean = {self.delta_mean:.4f}\n"
            f"  alpha*K    = {ALPHA * self.K:.4f}  (Paper V leading estimate)\n"
            f"  Fenichel bound delta < {FENICHEL_BOUND}: {status}\n"
        )


def compute_spectral_gap(
    K: float = K_BEST,
    alpha: float = ALPHA,
    Rc: float = RC,
    dr: float = DR,
    omega_m: float = OMEGA_M,
    eps: float = EPSILON,
    z_max: float = 3.0,
) -> SpectralGapResult:
    """Compute delta(N) = |dR/dN|/3 along the RLDS trajectory.

    Parameters
    ----------
    K:
        Coupling constant.
    alpha, Rc, dr, omega_m, eps:
        Model parameters.
    z_max:
        Maximum redshift for statistics.

    Returns
    -------
    SpectralGapResult
    """
    traj = RLDSTrajectory(K=K, alpha=alpha, Rc=Rc, dr=dr, omega_m=omega_m, eps=eps)
    dRdN = np.array([flow(r, K, alpha, Rc, dr, eps) for r in traj.R_arr])
    delta = np.abs(dRdN) / 3.0

    mask = traj.z_arr < z_max
    delta_late = delta[mask]
    z_late = traj.z_arr[mask]
    idx = int(np.argmax(delta_late))

    return SpectralGapResult(
        K=K,
        delta_max=float(delta_late[idx]),
        delta_mean=float(np.mean(delta_late)),
        z_at_max=float(z_late[idx]),
        fenichel_satisfied=bool(delta_late[idx] < FENICHEL_BOUND),
        N_arr=traj.N_arr,
        R_arr=traj.R_arr,
        delta_arr=delta,
    )


def spectral_gap_table(
    K_values: list[float] | None = None,
) -> list[SpectralGapResult]:
    """Table of delta_max across coupling values (Paper V, Appendix C).

    Parameters
    ----------
    K_values:
        List of K to evaluate.  Defaults to Paper V canonical table.

    Returns
    -------
    list[SpectralGapResult]
    """
    if K_values is None:
        K_values = [0.005, 0.0106, 0.016, 0.030, 0.050, 0.100, 0.200]
    return [compute_spectral_gap(K=K) for K in K_values]
