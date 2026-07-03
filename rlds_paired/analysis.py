"""
High-level analysis functions for the RLDS paired signature (Paper III).

Quick start
-----------
>>> from rlds_paired.analysis import k_scan, fsigma8_comparison
>>> result = k_scan()
>>> print(result.summary())
>>> result.save_figure("k_scan.png")

Paired signature
----------------
``k_scan`` measures the linear scaling  Delta(f*sigma8)/f*sigma8 = c * K
on the **full coupled-ODE background** (``core.growth.growth_factor``), with
the baseline taken at ``K=0`` on the *same* ODE path. The coefficient is
**negative** (growth suppression): the dark-sector energy exchange ``Q > 0``
drains the scalar. (An algebraic-background model gives the spurious *positive*
sign; that bug is described in ``core/growth.py``.)

The canonical value is not stored here; recompute it via
:func:`rlds_paired.reference.reference_results`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

from rlds_paired.constants import DR, K_BEST, OMEGA_M, RC, SIGMA8
from rlds_paired.core.algebraic import w_eff_alg
from rlds_paired.core.growth import growth_factor

# ---------------------------------------------------------------------------
# K-scan
# ---------------------------------------------------------------------------


@dataclass
class KScanResult:
    """Result of a K-scan for the paired signature coefficient.

    Attributes
    ----------
    K_values:
        Coupling constants used.
    delta_fs8:
        Fractional deviation  Delta(f*sigma8) / f*sigma8_ref  at each K.
    coefficient:
        Slope  c  from linear regression  delta = c * K.
    std_err:
        Standard error of c.
    r2:
        Goodness of linear fit.
    w_mins:
        w_min(K) values (phantom deviation).
    """

    K_values: np.ndarray
    delta_fs8: np.ndarray
    w_mins: np.ndarray
    coefficient: float
    std_err: float
    r2: float
    dR: float
    Omega_m: float
    z_eval: float

    def summary(self) -> str:
        return (
            f"K-scan paired signature  (full coupled-ODE background)\n"
            f"  Delta_R = {self.dR},  Omega_m = {self.Omega_m},  z = {self.z_eval}\n"
            f"  c  =  {self.coefficient:+.3f} ± {self.std_err:.3f}\n"
            f"  r² =  {self.r2:.2f}\n"
            f"  sign negative => growth suppression (Q>0 drains the scalar)\n"
        )

    def save_figure(self, path: str = "k_scan.png") -> None:
        """Save K-scan figure with two panels."""
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        valid = np.isfinite(self.delta_fs8)
        K_fit = np.linspace(0, self.K_values[-1] * 1.05, 100)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

        # Panel A: delta vs K
        ax1.scatter(
            self.K_values[valid],
            self.delta_fs8[valid],
            s=40,
            color="#2166ac",
            zorder=5,
            label="Full coupled-ODE model",
        )
        y_fit = self.coefficient * K_fit
        ax1.plot(
            K_fit,
            y_fit,
            "b--",
            lw=2,
            label=f"Linear fit:  c = {self.coefficient:.3f} ± {self.std_err:.3f}",
        )
        ax1.axhline(0, color="gray", lw=0.7, ls=":")
        ax1.set_xlabel("$K$", fontsize=13)
        ax1.set_ylabel(r"$\Delta(f\sigma_8)\,/\,f\sigma_8$", fontsize=13)
        ax1.set_title("(A) Paired-signature K-scan", fontsize=12)
        ax1.legend(fontsize=10)
        ax1.grid(alpha=0.2)

        # Panel B: w_min + 1 vs K
        valid_w = np.isfinite(self.w_mins)
        ax2.scatter(
            self.K_values[valid_w],
            self.w_mins[valid_w] + 1,
            s=40,
            color="#d73027",
            zorder=5,
        )
        ax2.axhline(0, color="gray", lw=0.8, ls="--", label="$w=-1$ (ΛCDM)")
        ax2.set_xlabel("$K$", fontsize=13)
        ax2.set_ylabel(r"$w_{\min}+1$", fontsize=13)
        ax2.set_title("(B) Phantom deviation vs K", fontsize=12)
        ax2.legend(fontsize=10)
        ax2.grid(alpha=0.2)

        fig.suptitle(
            f"RLDS Paired Signature  "
            f"($\\Omega_m={self.Omega_m}$, $R_c={RC}$, $\\Delta R={self.dR}$)",
            fontsize=13,
        )
        fig.tight_layout()
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved → {path}")


def k_scan(
    dR: float = DR,
    Omega_m: float = OMEGA_M,
    Rc: float = RC,
    K_min: float = 0.01,
    K_max: float = 0.40,
    n_points: int = 40,
    z_eval: float = 0.5,
    sigma8: float = SIGMA8,
    verbose: bool = True,
) -> KScanResult:
    """Compute the paired-signature K-scan.

    Calculates  Delta(f*sigma8) / f*sigma8_LCDM  across a range of K,
    fits a linear trend, and returns the slope c.

    Also computes w_min(K) to characterise the phantom excursion side of
    the paired signature.

    Parameters
    ----------
    dR, Omega_m, Rc:
        Model parameters.
    K_min, K_max, n_points:
        Scan range.
    z_eval:
        Redshift at which f*sigma8 is evaluated.
    sigma8:
        Fiducial sigma_8.
    verbose:
        Print progress.

    Returns
    -------
    KScanResult

    Notes
    -----
    The regression is performed by ``scipy.stats.linregress``, which fits the
    full model ``delta = slope * K + intercept`` (free intercept).  The
    reported slope ``c`` and ``r2`` are from that full model.  The intercept
    is approximately ``-0.012`` and arises from the ~1% offset between the
    ``K=0`` ODE baseline and pure LCDM ("ODE-vs-analytic offset" documented
    in :mod:`rlds_paired.reference`).  The formula ``delta = c * K`` (through
    origin) is a leading-order approximation; at ``K = K_BEST = 0.0106`` the
    direct ODE value (``~-0.0042``) differs from ``c * K`` (``~-0.0021``) by
    roughly 2x due to nonlinearity at the low-K edge of the scan range.
    """
    K_values = np.linspace(K_min, K_max, n_points)
    delta_arr = np.full(n_points, np.nan)
    wmin_arr = np.full(n_points, np.nan)

    fs8_ref = growth_factor(Omega_m, K=0.0, Rc=Rc, dR=dR, sigma8=sigma8, z_eval=z_eval)

    z_wgrid = np.linspace(0.01, 2.0, 200)

    for i, K in enumerate(K_values):
        try:
            fs8_k = growth_factor(
                Omega_m, K=K, Rc=Rc, dR=dR, sigma8=sigma8, z_eval=z_eval
            )
            delta_arr[i] = (fs8_k - fs8_ref) / fs8_ref

            w_arr = [w_eff_alg(z, Omega_m, K, Rc, dR) for z in z_wgrid]
            wmin_arr[i] = float(np.min(w_arr))

        except Exception:
            pass

        if verbose and (i + 1) % 10 == 0:
            print(f"  K={K:.3f}: delta={delta_arr[i]:+.5f}  w_min={wmin_arr[i]:.4f}")

    valid = np.isfinite(delta_arr)
    slope, _, r, _, se = linregress(K_values[valid], delta_arr[valid])

    return KScanResult(
        K_values=K_values,
        delta_fs8=delta_arr,
        w_mins=wmin_arr,
        coefficient=float(slope),
        std_err=float(se),
        r2=float(r**2),
        dR=dR,
        Omega_m=Omega_m,
        z_eval=z_eval,
    )


# ---------------------------------------------------------------------------
# f*sigma8 comparison
# ---------------------------------------------------------------------------


@dataclass
class FSigma8Result:
    """RLDS vs ΛCDM f*sigma8(z) curves.

    Attributes
    ----------
    z_grid:
        Redshift array.
    fs8_rlds:
        f*sigma8(z) for RLDS (K > 0).
    fs8_lcdm:
        f*sigma8(z) for K=0 (ΛCDM limit).
    K, Omega_m:
        Parameters used.
    """

    z_grid: np.ndarray
    fs8_rlds: np.ndarray
    fs8_lcdm: np.ndarray
    K: float
    Omega_m: float

    def save_figure(self, path: str = "fsigma8_comparison.png") -> None:
        """Save f*sigma8(z) comparison figure."""
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(self.z_grid, self.fs8_lcdm, "k--", lw=2, label="ΛCDM  ($K=0$)")
        ax.plot(
            self.z_grid,
            self.fs8_rlds,
            color="#2166ac",
            lw=2.5,
            label=f"RLDS  ($K={self.K:.3f}$,  $\\Omega_m={self.Omega_m:.3f}$)",
        )
        ax.set_xlabel("Redshift $z$", fontsize=13)
        ax.set_ylabel(r"$f\sigma_8(z)$", fontsize=13)
        ax.set_title("Growth rate: RLDS vs ΛCDM", fontsize=13)
        ax.legend(fontsize=12)
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved → {path}")


def fsigma8_comparison(
    K: float = K_BEST,
    Omega_m: float = OMEGA_M,
    Rc: float = RC,
    dR: float = DR,
    sigma8: float = SIGMA8,
    z_min: float = 0.05,
    z_max: float = 2.0,
    n_points: int = 60,
) -> FSigma8Result:
    """Compute RLDS and ΛCDM f*sigma8(z) curves (Paper III, Figure 4).

    Parameters
    ----------
    K, Omega_m, Rc, dR:
        Model parameters.
    sigma8:
        Fiducial sigma_8.
    z_min, z_max, n_points:
        Redshift grid for output.

    Returns
    -------
    FSigma8Result
    """
    z_grid = np.linspace(z_min, z_max, n_points)

    fs8_lcdm = np.array(
        [
            growth_factor(Omega_m, K=0.0, Rc=Rc, dR=dR, sigma8=sigma8, z_eval=z)
            for z in z_grid
        ]
    )
    fs8_rlds = np.array(
        [
            growth_factor(Omega_m, K=K, Rc=Rc, dR=dR, sigma8=sigma8, z_eval=z)
            for z in z_grid
        ]
    )

    return FSigma8Result(
        z_grid=z_grid,
        fs8_rlds=fs8_rlds,
        fs8_lcdm=fs8_lcdm,
        K=K,
        Omega_m=Omega_m,
    )
