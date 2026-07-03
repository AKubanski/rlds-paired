"""
Full two-phase RLDS ODE trajectory (Paper IV).

Solves the coupled system:
    dR/dN = -alpha*R*(1-R) + K*alpha*Psi(R)*S(R_c-R)
    du/dN = -3 + alpha*(1-R)

where N = ln(a), R = rho_DM/rho_nr, u = ln(rho_nr/rho_nr_ref).

Reference: Paper IV, Section 2.3 and Appendix B.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from rlds_paired.constants import (
    ALPHA,
    DR,
    EPSILON,
    K_BEST,
    OMEGA_M,
    OMEGA_R,
    RC,
    Z_INIT,
)

# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------


def psi(R: float | np.ndarray, eps: float = EPSILON) -> float | np.ndarray:
    """Singular restoring profile  Psi(R) = (1-R) / (R + eps)."""
    return (1.0 - R) / (R + eps)


def switching(y: float | np.ndarray, dR: float = DR) -> float | np.ndarray:
    """Smooth switching function  S(y) = (1 + tanh(y / Delta_R)) / 2."""
    return 0.5 * (1.0 + np.tanh(y / dR))


def flow(
    R: float,
    K: float,
    alpha: float = ALPHA,
    Rc: float = RC,
    dR: float = DR,
    eps: float = EPSILON,
) -> float:
    """RLDS reduced flow  F(R) = dR/dN."""
    # Clip only against the hard domain edges {0, 1}, NOT against the
    # regularisation eps. With eps=1e-6 the initial condition at z=1100
    # (1-R0 = 1.78e-9) was clamped to 1e-6 -- ~560x too large -- which
    # accelerated matter dilution and pushed the whole transition out to
    # z ~ 7-9 (nonphysical). 1e-12 leaves the true trajectory intact.
    R = float(np.clip(R, 1e-12, 1.0 - 1e-12))
    return -alpha * R * (1.0 - R) + K * alpha * psi(R, eps) * switching(Rc - R, dR)


# ---------------------------------------------------------------------------
# Two-phase ODE solver
# ---------------------------------------------------------------------------


class RLDSTrajectory:
    """Numerical solution of the coupled [R, u] ODE system.

    Parameters
    ----------
    K:
        Coupling constant.
    alpha:
        Dilution rate.
    Rc, dr, omega_m, eps:
        Model parameters.
    n_points:
        Number of output points.
    rtol, atol:
        ODE solver tolerances.

    Attributes
    ----------
    N_arr:
        E-fold time array.
    R_arr:
        Matter fraction R(N).
    u_arr:
        Log-density u(N).
    E2_arr:
        Squared Hubble rate E^2(z).
    z_arr:
        Redshift array.
    """

    def __init__(
        self,
        K: float = K_BEST,
        alpha: float = ALPHA,
        Rc: float = RC,
        dr: float = DR,
        omega_m: float = OMEGA_M,
        eps: float = EPSILON,
        n_points: int = 5000,
        rtol: float = 1e-10,
        atol: float = 1e-12,
    ) -> None:
        self.K = K
        self.alpha = alpha
        self.Rc = Rc
        self.dr = dr
        self.omega_m = omega_m
        self.eps = eps
        self._n = n_points
        self._rtol = rtol
        self._atol = atol
        self._solve()

    def _solve(self) -> None:
        N0 = -np.log(1.0 + Z_INIT)
        R0 = (
            self.omega_m
            * (1 + Z_INIT) ** 3
            / (self.omega_m * (1 + Z_INIT) ** 3 + (1.0 - self.omega_m))
        )

        def rhs(N: float, y: list) -> list:
            R, u = y
            R = float(np.clip(R, 1e-12, 1.0 - 1e-12))  # hard edges only (see flow())
            dR = flow(R, self.K, self.alpha, self.Rc, self.dr, self.eps)
            du = -3.0 + self.alpha * (1.0 - R)
            return [dR, du]

        sol = solve_ivp(
            rhs,
            [N0, 0.0],
            [R0, 0.0],
            method="RK45",
            rtol=self._rtol,
            atol=self._atol,
            max_step=0.05,
            dense_output=True,
        )
        if not sol.success:
            raise RuntimeError(f"ODE failed: {sol.message}")

        N_arr = np.linspace(N0, 0.0, self._n)
        y_arr = sol.sol(N_arr)
        R_arr = np.clip(y_arr[0], self.eps, 1.0 - self.eps)
        u_arr = y_arr[1]
        z_arr = np.exp(-N_arr) - 1.0

        rho_nr_ref = self.omega_m * (1 + Z_INIT) ** 3 + (1.0 - self.omega_m)
        E2_norm = rho_nr_ref * np.exp(u_arr[-1]) + OMEGA_R
        E2_arr = (rho_nr_ref * np.exp(u_arr) + OMEGA_R * (1 + z_arr) ** 4) / E2_norm

        self.N_arr = N_arr
        self.R_arr = R_arr
        self.u_arr = u_arr
        self.E2_arr = E2_arr
        self.z_arr = z_arr

    def R(self, z: float | np.ndarray) -> float | np.ndarray:
        """R(z) by interpolation. Scalar in -> float; array in -> ndarray."""
        out = interp1d(
            self.z_arr[::-1],
            self.R_arr[::-1],
            bounds_error=False,
            fill_value="extrapolate",
        )(z)
        return float(out) if np.ndim(z) == 0 else out

    def E2(self, z: float | np.ndarray) -> float | np.ndarray:
        """E^2(z) by interpolation. Scalar in -> float; array in -> ndarray."""
        out = np.maximum(
            interp1d(
                self.z_arr[::-1],
                self.E2_arr[::-1],
                bounds_error=False,
                fill_value="extrapolate",
            )(z),
            1e-30,
        )
        return float(out) if np.ndim(z) == 0 else out

    def E(self, z: float | np.ndarray) -> float | np.ndarray:
        """E(z) = H(z)/H_0. Scalar in -> float; array in -> ndarray."""
        e2 = self.E2(z)
        return float(np.sqrt(e2)) if np.ndim(z) == 0 else np.sqrt(e2)

    @property
    def R_today(self) -> float:
        """Present-day matter fraction R(z=0) (endpoint of the integration)."""
        return float(self.R_arr[-1])

    @property
    def R_star(self) -> float:
        """Asymptotic attractor R* = stable fixed point of the reduced flow F(R)=0.

        This is the FUTURE (N -> +inf) value, not R(z=0). In the attractor
        branch (R < Rc so the switch S -> 1) the reduced flow is
        F(R) = alpha*(1-R)*(K - R^2)/(R+eps), whose physical root is R* ~ sqrt(K).
        We locate it by bracketed root-finding so it stays correct if the
        switching/regularisation details change.

        Note: with the clip fix, R(z=0) (= ``R_today`` ~ 0.297 at the canonical
        benchmark) and the attractor (~ sqrt(K) ~ 0.103) are distinct. They
        coincided only under the old clipping bug, which drove the trajectory
        into the attractor regime far too early (z ~ 7-9).
        """
        from scipy.optimize import brentq

        f = lambda R: flow(R, self.K, self.alpha, self.Rc, self.dr, self.eps)
        lo, hi = 1e-8, max(self.Rc - 1e-6, 1e-6)
        try:
            if f(lo) * f(hi) < 0.0:
                return float(brentq(f, lo, hi, xtol=1e-12))
        except (ValueError, RuntimeError):
            pass
        # Fallback: analytic attractor of the S->1 branch.
        return float(np.sqrt(max(self.K, 0.0)))

    @property
    def wmin(self) -> float:
        """Minimum w_eff over z in [0.01, 2.0]."""
        from rlds_paired.core.algebraic import w_eff_alg

        z_grid = np.linspace(0.01, 2.0, 300)
        return float(
            min(w_eff_alg(z, self.omega_m, self.K, self.Rc, self.dr) for z in z_grid)
        )
