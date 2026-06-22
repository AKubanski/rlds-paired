# Data sources

Every dataset that enters the shooting-MCMC likelihood (`rlds_paired.mcmc`,
via `rlds_paired.data`). Two of the three are tabulated directly in the code;
only the supernova catalogue is an external file.

---

## 1. Pantheon+SH0ES (Type Ia supernovae)

- **Reference:** Brout et al. (2022), *The Pantheon+ Analysis: Cosmological
  Constraints*, ApJ 938, 110; Scolnic et al. (2022), ApJ 938, 113.
- **Release:** public, https://github.com/PantheonPlusSH0ES/DataRelease
- **File:** `Pantheon+SH0ES.dat`
- **SHA-256:** `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8`
- **Rows:** 1701 SNe in the file; **1590** pass the selection below.
- **Columns used** (0-indexed, whitespace-separated):
  - `[2]  zHD`               — Hubble-diagram redshift
  - `[10] MU_SH0ES`          — SH0ES-calibrated distance modulus
  - `[11] MU_SH0ES_ERR_DIAG` — diagonal distance-modulus error
- **Selection:** `0.01 < z < 2.5`, `MU_SH0ES > 0`, `0 < err < 10`.
- **Likelihood:** diagonal error with analytic marginalisation over the
  absolute-magnitude offset (the S1/S2/S3 form in `mcmc.log_likelihood_SN`).
  The full systematic covariance is **not** used in the reference run; a
  full-covariance path is available via `use_full_cov=True` (supply the
  inverse covariance), but it is **not** what produced the archived chain.
- **Placement:** `~/cosmology_data/Pantheon+SH0ES.dat`, or set `$RLDS_DATA_DIR`.
  The `_` vs `+` filename variants are both accepted.

> Bitwise reproduction of the archived chain requires this exact file;
> `rlds_paired.data.verify_pantheon` checks the SHA-256.

---

## 2. DESI DR2 BAO

- **Reference:** DESI Collaboration (2025), *DESI DR2 Results II: Measurements
  of Baryon Acoustic Oscillations and Cosmological Constraints*
  (arXiv:2503.14738).
- **Form:** 13 distance ratios at seven effective redshifts (BGS `DV/rd` at
  z=0.295; `DM/rd` and `DH/rd` for LRG1 0.510, LRG2 0.706, LRG3+ELG1 0.934,
  ELG2 1.321, QSO 1.484; Lya 2.330), with the **full 13x13 covariance**
  (block-diagonal: `DM`-`DH` correlated within a tracer; the BGS `DV` point and
  different redshift bins uncorrelated).
- **Files:** the official combined data vector and covariance,
  `desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt`, distributed via the
  CobayaSampler/bao_data repository, are bundled verbatim under
  `rlds_paired/desi_dr2/` and loaded in their native row order so the
  covariance indices stay aligned with the data vector.
- **Location:** `get_desi_dr2_bao()` returns `(entries, data, cov, cov_inv)`.

---

## 3. Planck 2018 distance prior

- **Reference:** Planck Collaboration (2020), Aghanim et al., A&A 641, A6
  (compressed CMB / distance-prior summary).
- **Form:** two-vector `[R, omega_b]` with diagonal errors, where `R` is the
  CMB shift parameter `R = sqrt(Omega_m) H0 D_M(z*) / c` at `z* = 1089.9`.
  - means  = `[1.7488, 0.02237]`
  - errors = `[0.0074, 0.00015]`
- **Location:** `rlds_paired.data.PLANCK2018_MEANS / PLANCK2018_ERRORS`
  (returned by `get_planck2018_distance_prior()`).

---

### Notes on redistribution

Pantheon+SH0ES is publicly released for research use; only the data vector and
diagonal error enter here. The DESI DR2 BAO data vector and covariance are the
collaboration's publicly released files (via CobayaSampler/bao_data); the
Planck distance prior is a summary statistic reproduced from the cited
publication.
