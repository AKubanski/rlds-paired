"""
Linear growth factor on the FULL coupled-ODE RLDS background (Paper III, Sec 2.5 / App B.4).

This replaces the previous *algebraic-background* growth model, which gave the WRONG
SIGN for the paired coefficient (+0.57, spurious enhancement). Root cause: the algebraic
ansatz F(R) = 1 + K*tanh((Rc - R)/dR) makes rho_DE RISE toward low z, whereas the true
coupled dark sector (energy exchange Q > 0 drains the scalar) makes rho_DE FALL. The two
backgrounds therefore yield opposite-sign growth responses. The full coupled ODE gives
c < 0 (suppression), in agreement with Paper III's stated full-ODE result.

PREREQUISITE: the clip fix in core/ode.py (clip lower bound 1e-12, not 1e-6). Without it
the RLDS trajectory itself is corrupted (transition pushed to z ~ 7-9).

Standard sub-horizon growth equation (Paper III Eq. 8):
    D'' + (3/a + E'/E) D' - (3/2) * Omega_m0 / (a^5 E^2) * D = 0
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import odeint

from rlds_paired.constants import DR, OMEGA_M, RC, SIGMA8
from rlds_paired.core.ode import RLDSTrajectory

_A = np.linspace(0.02, 1.0, 3000)
_ZA = 1.0 / _A - 1.0


def growth_factor(
    omega_m: float = OMEGA_M,
    K: float = 0.0,
    Rc: float = RC,
    dR: float = DR,
    sigma8: float = SIGMA8,
    z_eval: float = 0.5,
) -> float:
    """f*sigma8(z_eval) on the full coupled-ODE RLDS background. K=0 -> baseline."""
    traj = RLDSTrajectory(K=K, Rc=Rc, dr=dR, omega_m=omega_m)
    E2 = np.maximum(np.interp(_ZA, traj.z_arr[::-1], traj.E2_arr[::-1]), 1e-30)
    E = np.sqrt(E2)
    dEda = np.gradient(E, _A)
    cDp = -(3.0 / _A + dEda / E)
    cD = 1.5 * omega_m / (_A**5 * E2)

    def _ode(y, a):
        D, Dp = y
        return [Dp, np.interp(a, _A, cDp) * Dp + np.interp(a, _A, cD) * D]

    sol = odeint(_ode, [0.02, 1.0], _A, rtol=1e-9, atol=1e-11)
    D = sol[:, 0] / sol[-1, 0]
    Dp = sol[:, 1] / sol[-1, 0]
    fs8 = _A * Dp / D * sigma8 * D
    return float(np.interp(z_eval, _ZA[::-1], fs8[::-1]))
