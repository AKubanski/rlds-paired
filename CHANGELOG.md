# Changelog

> **Versioning reset (2026-06-21).** The earlier `rlds-paired` Zenodo records
> (v1.0.0 `…20771955`, v1.1.0 `…20781102`, concept `…20771954`) were **deleted**
> to avoid confusion. The DESI DR2 release will be published **fresh as
> v1.0.0** with a new concept DOI. The development history below is retained as
> an internal record; the public v1.0.0 release notes will be consolidated into
> a single entry at release time.

## [Unreleased] — DESI DR2 + fresh v1.0.0 (1.0.0.dev0)

- **BAO data migrated from DESI DR1 to DESI DR2** (arXiv:2503.14738). The
  reference BAO likelihood previously used the 11-point DESI DR1 (2024) data
  vector with diagonal errors (mislabeled internally as "DR2"); it now uses the
  official 13-point DESI DR2 combined data vector **with the full 13x13
  block-diagonal covariance** (adds the QSO bin; `DM`-`DH` correlations within
  each tracer). The official files
  `desi_gaussian_bao_ALL_GCcomb_{mean,cov}.txt` (CobayaSampler/bao_data) are
  bundled verbatim under `rlds_paired/desi_dr2/` and loaded in native order.
  - `data.get_desi_dr2_bao()` now returns `(entries, data, cov, cov_inv)`
    instead of a list of `(z, obs, val, sigma)` tuples; the `DESI_DR2_BAO`
    constant is removed.
  - `mcmc.log_likelihood_BAO` now evaluates `chi2 = d^T C^-1 d` over the full
    covariance.
  - **Canonical numbers will change**: a fresh shooting MCMC on DR2 is required.
    Until that chain is regenerated, `REFERENCE_BEST_FIT` is the *former* DR1
    best-fit retained as a pipeline-drift anchor, and `REFERENCE_CHI2` is the
    DR2 chi2 at that point (753.91, provisional). The canonical DR2 best-fit
    and chi2 replace both once the re-run completes.
  - The `[mcmc]` extra now also declares `tqdm` (emcee's progress bar, used by
    `run_mcmc(progress=True)`); it was previously missing.
  - All 40 tests pass.
- **Canonical numbers synced to the DR2 chain.** `constants.py` and the MCMC
  reference now hold the DR2 posterior medians: `K=0.0106`, `α=3.0102` (still
  **100% railed** at the prior floor 3.01 — phantom remains prior-driven),
  `R_c=0.565`, `ΔR=0.170`, `Ω_m=0.297`, `R*=√K=0.1029` (~10% matter / 90% DE),
  with recomputed `c=−0.193±0.005` (r²=0.97), `δ_max=0.244 @ z≈0.40`,
  `R(z=0)=0.298`, `z_acc≈0.71`. `REFERENCE_BEST_FIT` is the DR2 argmax,
  `REFERENCE_CHI2=750.95`.
- **Fix:** `RLDSTrajectory`'s default `K` was a hardcoded `0.012` (DR1) while
  its other defaults already read from `constants`; it now reads `K_BEST`, so
  the canonical trajectory tracks the convention sheet. `cli verify` and the
  test bands for paired `c` and `δ_max` were updated to the DR2 values.

## [Unreleased] — canonical parameter sync

  - **Release hardening** (pre-1.0.0): added CITATION.cff, .zenodo.json, SECURITY.md,
    CODE_OF_CONDUCT.md, RELEASING.md and scripts/bump_version.py; release workflow moved to
    a manual `workflow_dispatch` pipeline (test + build gated before tag/Release) with PyPI
    Trusted Publishing (OIDC, no stored token); RLDSTrajectory.R/E/E2 now correctly accept
    ndarray input; MANIFEST.in now ships tests/ (incl. conftest.py) so sdist tests pass;
    removed the non-functional `spectral --figure` argument; added tests/test_reference.py
    (headline values c, delta_max, R*, R(z=0)) and tests/test_cli.py (30 tests total); the CLI forces UTF-8 stdout/stderr so the ✓/Greek glyphs no
    longer crash on Windows cp1251/cp1252 consoles or pipes; the scan CLI prints a coarse-grid note when n<40; code
    formatted with Black.

- Switching width `DR` set to the shooting best-fit **0.15** (was 0.12), matching the
  locked canonical set used in Papers III/IV/V (Rc=0.51, dR=0.15, Omega_m=0.296).
  Recomputed headline: paired c = -0.148 +/- 0.004 (r2=0.97); delta_max = 0.245
  (Fenichel violated). All 30 tests pass.

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.1.0] — 2026-06

### Added — shooting MCMC (Papers III–V)

- `rlds_paired/mcmc.py`: faithful reproduction of the canonical *shooting* MCMC.
  Five sampled parameters `[K, H0, R_c, Delta_R, alpha]`; `Omega_m` is **derived**
  by shooting so `R(z=0) = Omega_m` exactly (posterior median 0.296). Provides
  `ShootingRLDS`, `log_prior/log_likelihood_*/log_posterior`, and
  `run_mcmc(seed=...)`. Uses `eps=1e-6` (chain fidelity); `analysis` keeps
  `1e-12` (sub-σ difference, documented).
  - Diagonal Pantheon+ likelihood with analytic absolute-magnitude
    marginalisation (reference run); optional `use_full_cov=True`.
  - Optional `fast_distance` path (one comoving-distance interpolator per
    model, ~100× faster, agrees with the `quad` path to < 2e-5 mag).
  - `reproduce_reference_chi2()` recomputes χ² at the archived best-fit; it
    equals the recorded **727.665** to < 0.01 (exact regression check).
- `rlds_paired/data.py`: single data-provenance module — Pantheon+SH0ES
  (SHA-256 verified loader; `$RLDS_DATA_DIR` honoured), DESI DR2 BAO, and the
  Planck 2018 distance prior.
- CLI `mcmc` subcommand: `rlds-paired mcmc` (reproduce best-fit χ²) and
  `rlds-paired mcmc --run [--steps --walkers --seed --out]` (new chain).
- Optional dependency extra: `pip install "rlds-paired[mcmc]"` (emcee).
- `tests/test_mcmc.py`: prior logic (always) + χ²-reproduction, fast-vs-slow,
  and self-consistency (skipped when the Pantheon file is absent).
- `SOURCES.md` (data provenance) and `chains/METADATA.md` (chain passport).

### Notes

- The `alpha` posterior is 100 % railed against the prior floor 3.01: the
  phantom excursion is prior-dependent and the data are consistent with
  `w = −1`. Documented in `chains/METADATA.md`.
- The concept DOI `10.5281/zenodo.20771954` is unchanged; this release adds
  the MCMC code referenced by Papers III/IV.

---

## [1.0.0] — 2026-06

First public release.

### Added

**Core physics**
- `core/algebraic.py`: Algebraic IDE expansion history — `E_squared_alg`,
  `solve_R`, `F_transition`, `w_eff_alg`.  K=0 recovers exact ΛCDM.
- `core/ode.py`: Full two-phase coupled [R, u] ODE trajectory (`RLDSTrajectory`),
  `flow`, `psi`, `switching`.
- `core/growth.py`: Linear growth factor `growth_factor()` on the **full
  coupled-ODE** RLDS background (Paper III, Sec 2.5 / App B.4).
- `core/spectral.py`: Pointwise spectral-gap validation δ(N) = |dR/dN|/3
  (`compute_spectral_gap`, `spectral_gap_table`, `SpectralGapResult`).

**Analysis**
- `analysis.py`: `k_scan()` — K-scan of the paired-signature coefficient
  (linear scaling, negative slope = growth suppression).
- `analysis.py`: `fsigma8_comparison()` — RLDS vs ΛCDM f·σ₈(z) curves.

**Infrastructure**
- `constants.py`: Canonical benchmark *inputs* (Papers III–V). Computed
  results are not stored here — see `reference.py`.
- `reference.py`: `reference_results()` — recomputes the headline numbers
  (paired coefficient, δ_max) from the model on demand, so a stale literal
  can never ship.
- `exceptions.py`: `RLDSPairedError`, `ConvergenceError`, `PhysicsError`,
  `ConfigurationError`.
- `cli.py`: `rlds-paired verify | scan | spectral | fsigma8`.
- Lazy imports via PEP 562 `__getattr__` — fast `import rlds_paired`.
- Full test suite: 30 tests across 6 modules, all passing.
- GitHub Actions CI (Ubuntu, macOS, Windows × Python 3.10–3.12).

### Fixed (pre-release)

- **ODE initial-condition clipping** (`core/ode.py`): the slow-manifold clip
  lower bound was `1e-6`, ~560× larger than the physical `1 − R(z=1100)`. It
  clamped the high-z initial condition and pushed the matter→DE transition to
  z ≈ 7–9. Lowered to `1e-12` at both flow sites. R(z=0) now ≈ 0.297 (was
  0.110); the R=2/3 acceleration crossing now lands at z ≈ 0.72 (was ≈ 9).
- **Paired-signature sign** (`core/growth.py`): growth was being computed on
  the *algebraic* background F = 1 + K·tanh(...), where ρ_DE rises toward low
  z — giving a spurious **positive** coefficient (+0.57). Switched to the full
  coupled-ODE background (Q > 0 drains the scalar), yielding the correct
  **negative** sign.
- `R_star` now returns the genuine attractor (stable fixed point of the reduced
  flow, ≈ √K); added `R_today` = R(z=0). Pre-fix these coincided only because
  of the clipping bug.

### Canonical results (Paper III–V, recomputed from code)

| Quantity | Value |
|---|---|
| Paired coefficient *c* | ≈ −0.14 (negative: growth suppression) |
| Attractor *R** | 0.109 (= √K, K = 0.012) |
| Present-day matter fraction *R(z=0)* | 0.297 |
| Spectral gap *δ_max* | ≈ 0.245 at z ≈ 0.42 — Fenichel bound (0.10) **violated**; reduction marginal (~25%) |
