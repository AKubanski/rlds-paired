"""
rlds_paired.data
================
Single provenance point for every dataset that enters the paired-signature
likelihood (Papers III--V).

Three sources are used, exactly as in the canonical shooting MCMC:

1. Pantheon+SH0ES supernovae  -- external file, loaded here (with SHA-256 check).
2. DESI DR2 BAO               -- 11 distance-ratio points, tabulated below.
3. Planck 2018 distance prior -- CMB shift parameter R and omega_b, below.

The supernova likelihood uses the SH0ES distance modulus column with its
DIAGONAL error and an analytic marginalisation over the absolute magnitude
(the S1/S2/S3 trick in ``mcmc.log_likelihood_SN``). The full systematic
covariance is *not* used; this is a deliberate simplification of the
reference run and is documented as such.

Reproduction note
-----------------
The archived chain was produced with these exact inputs. ``load_pantheon_plus``
verifies the file by SHA-256 so that a re-run uses byte-identical data.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# 1. PANTHEON+SH0ES
# ---------------------------------------------------------------------------
#
# Source : Brout et al. (2022), ApJ 938, 110; Scolnic et al. (2022).
#          Public release: github.com/PantheonPlusSH0ES/DataRelease
# File   : Pantheon+SH0ES.dat  (1701 SNe; 1590 pass the 0.01 < z < 2.5 cut)
# Columns used (0-indexed, whitespace-separated):
#     [2]  zHD                  -- Hubble-diagram redshift
#     [10] MU_SH0ES             -- SH0ES-calibrated distance modulus
#     [11] MU_SH0ES_ERR_DIAG    -- diagonal distance-modulus error
# Selection: 0.01 < z < 2.5, MU_SH0ES > 0, 0 < err < 10.
#
PANTHEON_FILENAMES = ("Pantheon+SH0ES.dat", "Pantheon_SH0ES.dat")
PANTHEON_SHA256 = (
    "1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8"
)
PANTHEON_N_EXPECTED = 1590  # SNe after the selection cut

_Z_MIN, _Z_MAX = 0.01, 2.5
_sn_cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}


def _resolve_pantheon_path(data_dir: str | Path | None) -> Path:
    """Find the Pantheon file, tolerating the '+' vs '_' filename variants.

    Resolution order: explicit ``data_dir`` arg, then ``$RLDS_DATA_DIR``,
    then ``~/cosmology_data``.
    """
    import os

    if data_dir is None:
        data_dir = os.environ.get("RLDS_DATA_DIR") or (Path.home() / "cosmology_data")
    data_dir = Path(data_dir)
    # If a file was passed directly, use it.
    if data_dir.is_file():
        return data_dir
    for name in PANTHEON_FILENAMES:
        candidate = data_dir / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Pantheon file not found in {data_dir} "
        f"(looked for {', '.join(PANTHEON_FILENAMES)})."
    )


def sha256_of(path: str | Path) -> str:
    """Return the SHA-256 hex digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify_pantheon(path: str | Path, *, strict: bool = False) -> bool:
    """Check the Pantheon file against the recorded SHA-256.

    Returns True on match. On mismatch: raise if ``strict``, else warn-return
    False (a different but valid release may still reproduce the science to
    within sampling error; bitwise reproduction requires the exact file).
    """
    digest = sha256_of(path)
    if digest == PANTHEON_SHA256:
        return True
    msg = (
        f"Pantheon SHA-256 mismatch:\n  expected {PANTHEON_SHA256}\n  got      "
        f"{digest}\nBitwise reproduction of the archived chain requires the "
        f"exact reference file."
    )
    if strict:
        raise ValueError(msg)
    import warnings

    warnings.warn(msg, RuntimeWarning, stacklevel=2)
    return False


def load_pantheon_plus(
    data_dir: str | Path | None = None, *, verify: bool = True
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load (z, mu, mu_err) for the SNe passing the reference selection.

    Parameters
    ----------
    data_dir : path to the directory containing the file, or the file itself.
               Defaults to ``~/cosmology_data``.
    verify   : if True, check the SHA-256 (non-strict warning on mismatch).

    Returns
    -------
    (z, mu, mu_err) : three 1-D float arrays of equal length.
    """
    path = _resolve_pantheon_path(data_dir)
    key = str(path.resolve())
    if key in _sn_cache:
        return _sn_cache[key]

    if verify:
        verify_pantheon(path, strict=False)

    rows = []
    with open(path, "r") as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.split()
            if len(parts) < 12:
                continue
            try:
                z, mu, err = float(parts[2]), float(parts[10]), float(parts[11])
            except ValueError:
                continue
            if _Z_MIN < z < _Z_MAX and mu > 0 and 0 < err < 10:
                rows.append((z, mu, err))

    arr = np.asarray(rows, dtype=float)
    if arr.shape[0] != PANTHEON_N_EXPECTED:
        import warnings

        warnings.warn(
            f"Pantheon selection returned {arr.shape[0]} SNe "
            f"(reference run used {PANTHEON_N_EXPECTED}).",
            RuntimeWarning,
            stacklevel=2,
        )
    out = (arr[:, 0], arr[:, 1], arr[:, 2])
    _sn_cache[key] = out
    return out


# ---------------------------------------------------------------------------
# 2. DESI DR2 BAO
# ---------------------------------------------------------------------------
#
# Source : DESI Collaboration (2025), "DESI DR2 Results II: Measurements of
# Baryon Acoustic Oscillations and Cosmological Constraints" (arXiv:2503.14738),
# combined galaxy + quasar + Lyman-alpha BAO (GCcomb). Distributed via the
# CobayaSampler/bao_data repository as
# desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt; both files are
# bundled verbatim under desi_dr2/ and loaded in their native row order so the
# covariance indices stay aligned with the data vector.
#
# Thirteen distance ratios at seven effective redshifts: BGS DV/rd at z=0.295;
# DM/rd and DH/rd for LRG1 (0.510), LRG2 (0.706), LRG3+ELG1 (0.934),
# ELG2 (1.321), QSO (1.484); and Lya (2.330), given as DH/rd then DM/rd in the
# DESI ordering. The 13x13 covariance is block-diagonal: DM and DH correlate
# within a tracer, while the BGS DV point and different redshift bins are
# mutually uncorrelated.
#
_DESI_DR2_DIR = Path(__file__).resolve().parent / "desi_dr2"
_DESI_DR2_MEAN_FILE = _DESI_DR2_DIR / "desi_gaussian_bao_ALL_GCcomb_mean.txt"
_DESI_DR2_COV_FILE = _DESI_DR2_DIR / "desi_gaussian_bao_ALL_GCcomb_cov.txt"

# DESI quantity labels mapped to the short forms used by the model interface.
_DESI_QUANTITY_MAP = {"DV_over_rs": "DV", "DM_over_rs": "DM", "DH_over_rs": "DH"}

_desi_dr2_cache: tuple | None = None


def get_desi_dr2_bao() -> tuple[
    list[tuple[float, str]], np.ndarray, np.ndarray, np.ndarray
]:
    """Load the DESI DR2 BAO data vector and full (block-diagonal) covariance.

    Returns a 4-tuple ``(entries, data, cov, cov_inv)`` where ``entries`` is a
    list of ``(z, quantity)`` pairs (quantity in ``{"DV", "DM", "DH"}``) in the
    native DESI row order, ``data`` is the (13,) array of measured distance
    ratios, and ``cov`` / ``cov_inv`` are the (13, 13) covariance and its
    inverse. The rows of ``cov`` align one-to-one with ``entries`` and
    ``data``. The result is cached after the first call.
    """
    global _desi_dr2_cache
    if _desi_dr2_cache is not None:
        return _desi_dr2_cache

    entries: list[tuple[float, str]] = []
    values: list[float] = []
    with open(_DESI_DR2_MEAN_FILE, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            z_str, val_str, label = line.split()
            quantity = _DESI_QUANTITY_MAP.get(label)
            if quantity is None:
                raise ValueError(f"Unknown DESI DR2 quantity label: {label!r}")
            entries.append((float(z_str), quantity))
            values.append(float(val_str))

    data = np.asarray(values, dtype=float)
    cov = np.loadtxt(_DESI_DR2_COV_FILE, dtype=float)
    n = len(entries)
    if cov.shape != (n, n):
        raise ValueError(
            f"DESI DR2 covariance shape {cov.shape} does not match {n} points"
        )
    cov_inv = np.linalg.inv(cov)
    _desi_dr2_cache = (entries, data, cov, cov_inv)
    return _desi_dr2_cache


# ---------------------------------------------------------------------------
# 3. PLANCK 2018 DISTANCE PRIOR
# ---------------------------------------------------------------------------
#
# Source : Planck 2018 (Aghanim et al. 2020); compressed CMB likelihood.
# Two-vector [R, omega_b] with diagonal errors, where R is the CMB shift
# parameter R = sqrt(Omega_m) * H0 * D_M(z_*) / c, evaluated at z_* = 1089.9.
#
PLANCK2018_MEANS = np.array([1.7488, 0.02237])      # [R, omega_b]
PLANCK2018_ERRORS = np.array([0.0074, 0.00015])     # 1-sigma, diagonal
PLANCK2018_ZSTAR = 1089.9                            # recombination redshift


def get_planck2018_distance_prior() -> tuple[np.ndarray, np.ndarray]:
    """Return (means, inverse_covariance) for the [R, omega_b] distance prior."""
    cov_inv = np.linalg.inv(np.diag(PLANCK2018_ERRORS ** 2))
    return PLANCK2018_MEANS, cov_inv


__all__ = [
    "load_pantheon_plus",
    "verify_pantheon",
    "sha256_of",
    "PANTHEON_SHA256",
    "PANTHEON_N_EXPECTED",
    "get_desi_dr2_bao",
    "get_planck2018_distance_prior",
    "PLANCK2018_MEANS",
    "PLANCK2018_ERRORS",
    "PLANCK2018_ZSTAR",
]
