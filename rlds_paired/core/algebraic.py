"""
Algebraic IDE expansion history (Paper III model).

E^2(z) = [Omega_m*(1+z)^3 + Omega_r*(1+z)^4 + Omega_DE0*F(R)] / norm
F(R)   = 1 + K * tanh((R_c - R) / Delta_R)

R is solved self-consistently via fixed-point iteration.
K=0 recovers standard ΛCDM exactly.

This algebraic background reproduces the expansion history. NOTE: it is *not*
used for the paired-signature coefficient — that ansatz makes rho_DE rise toward
low z and yields the wrong (positive) growth sign. The paired coefficient is
computed on the full coupled-ODE background (``core/growth.py``), giving c < 0
(≈ −0.19, computed; see ``rlds_paired.reference``).

Reference: Paper III, Appendix B (algebraic approximation).
"""

from __future__ import annotations

import numpy as np

from rlds_paired.constants import DR, OMEGA_M, OMEGA_R, RC


def F_transition(
    R: float | np.ndarray,
    K: float,
    Rc: float = RC,
    dR: float = DR,
) -> float | np.ndarray:
    """DE modification factor  F(R) = 1 + K * tanh((R_c - R) / Delta_R).

    Sign convention: when R < R_c (matter subdominant), tanh is positive,
    boosting rho_DE.  This produces w_eff < -1 (phantom crossing) and
    matches the RLDS switching function S(R_c - R) in the ODE formulation.

    Parameters
    ----------
    R:
        Matter fraction in (0, 1).
    K:
        Coupling constant.  K=0 → F=1 → exact ΛCDM.
    Rc, dR:
        Switching threshold and width.

    Returns
    -------
    float | np.ndarray
        F(R) ≥ 0.
    """
    return 1.0 + K * np.tanh((Rc - R) / dR)


def solve_R(
    z: float,
    omega_m: float = OMEGA_M,
    K: float = 0.0,
    Rc: float = RC,
    dR: float = DR,
    tol: float = 1e-12,
    maxiter: int = 200,
) -> float:
    """Matter fraction R(z) from fixed-point iteration.

    Solves:  R = Omega_m*(1+z)^3 / [Omega_m*(1+z)^3 + Omega_DE0*F(R)]

    Parameters
    ----------
    z:
        Redshift.
    omega_m:
        Present-day matter density.
    K, Rc, dR:
        Model parameters.
    tol:
        Convergence tolerance.
    maxiter:
        Maximum iterations.

    Returns
    -------
    float
        Self-consistent R(z).
    """
    omega_de0 = 1.0 - omega_m
    om_z = omega_m * (1.0 + z) ** 3
    R = om_z / (om_z + omega_de0)  # ΛCDM initial guess

    for _ in range(maxiter):
        F = F_transition(R, K, Rc, dR)
        R_new = om_z / (om_z + omega_de0 * F)
        if abs(R_new - R) < tol:
            return float(R_new)
        R = R_new

    return float(R)


def E_squared_alg(
    z: float,
    omega_m: float = OMEGA_M,
    K: float = 0.0,
    Rc: float = RC,
    dR: float = DR,
    omega_r: float = OMEGA_R,
) -> float:
    """Normalised squared Hubble rate E^2(z) = H^2(z)/H_0^2 (algebraic model).

    Parameters
    ----------
    z:
        Redshift.
    omega_m:
        Present-day matter density parameter.
    K, Rc, dR:
        Model parameters.  K=0 → exact ΛCDM.
    omega_r:
        Radiation density parameter.

    Returns
    -------
    float
        E^2(z) ≥ 1e-30.
    """
    omega_de0 = 1.0 - omega_m

    def _rho(zz: float) -> float:
        R = solve_R(zz, omega_m, K, Rc, dR)
        F = F_transition(R, K, Rc, dR)
        return omega_m * (1.0 + zz) ** 3 + omega_r * (1.0 + zz) ** 4 + omega_de0 * F

    numerator = _rho(z)
    denominator = _rho(0.0)
    return max(numerator / denominator, 1e-30)


def E_alg(
    z: float,
    omega_m: float = OMEGA_M,
    K: float = 0.0,
    Rc: float = RC,
    dR: float = DR,
) -> float:
    """Normalised Hubble rate E(z) = H(z)/H_0 (algebraic model)."""
    return float(np.sqrt(E_squared_alg(z, omega_m, K, Rc, dR)))


def w_eff_alg(
    z: float,
    omega_m: float = OMEGA_M,
    K: float = 0.0,
    Rc: float = RC,
    dR: float = DR,
    dz: float = 0.005,
) -> float:
    """Effective dark-energy equation of state from rho_DE evolution.

    w_eff(z) = -1 + (1+z)/(3*rho_DE) * d(rho_DE)/dz

    Parameters
    ----------
    z:
        Redshift.
    omega_m, K, Rc, dR:
        Model parameters.
    dz:
        Finite-difference step.

    Returns
    -------
    float
        w_eff(z).
    """
    omega_de0 = 1.0 - omega_m

    def rho_de(zz: float) -> float:
        R = solve_R(zz, omega_m, K, Rc, dR)
        return omega_de0 * F_transition(R, K, Rc, dR)

    rho_c = rho_de(z)
    if rho_c < 1e-30:
        return -1.0

    rho_p = rho_de(z + dz)
    rho_m = rho_de(z - dz) if z > dz else rho_de(z)
    drho = (rho_p - rho_m) / (2.0 * dz) if z > dz else (rho_p - rho_c) / dz

    return -1.0 + (1.0 + z) / (3.0 * rho_c) * drho
