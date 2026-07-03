# Reproducing the RLDS figures (Papers III–V)

This file maps each headline figure / result of the RLDS observational series
to the exact entry point in `rlds-paired` that reproduces it, on the **DESI DR2
canonical frame** (`K=0.0106`, `α=3.0102`, `R_c=0.565`, `ΔR=0.170`,
`Ω_m=0.297`). All headline numbers are recomputed from code; none are read from
a stored literal (see `rlds_paired/reference.py`).

## Paper III — paired growth–expansion signature

- **Paired coefficient `c` and the K-scan:** `rlds_paired.analysis.k_scan()`
  - Slope of `Δ(f·σ₈)/(f·σ₈)_ΛCDM = c·K` on the full coupled-ODE background.
  - Canonical K-scan defaults (`K ∈ [0.01, 0.40]`, `n = 40`, `z_eval = 0.5`)
    reproduce `c ≈ −0.1936 ± 0.0053` (`r² ≈ 0.97`). A coarser grid will not
    reproduce this slope — keep the defaults.
  - **Regression note:** `scipy.linregress` fits the full model
    `δ = slope·K + intercept`; the intercept is approximately `−0.012` and
    reflects the ~1% offset between the `K=0` ODE baseline and pure ΛCDM.
    The reported `r² = 0.972` is for the full model with intercept, not for
    the zero-intercept approximation `δ = c·K`. At `K = K_BEST = 0.0106`,
    the direct ODE value (`≈ −0.0042`) differs from the linear prediction
    `c·K` (`≈ −0.0021`) by ~2×, consistent with the O(K²) nonlinearity
    at the low-K edge of the scan.
- **f·σ₈ comparison curves:** `rlds_paired.analysis.fsigma8_comparison()`.

## Paper IV — cosmological consequences (RLDS closure)

- **Attractor `R*`, present-day `R(z=0)`, expansion `E(z)`:**
  `rlds_paired.core.ode.RLDSTrajectory`
  (`.R_star ≈ 0.1027`, `.R_today ≈ 0.2974`).
- **`w_eff(z)`, algebraic `E²(z)`, acceleration onset `z_acc`:**
  `rlds_paired.core.algebraic` (`w_eff_alg`, `E_squared_alg`, `solve_R`,
  `F_transition`); `z_acc = 0.714` (threshold `R = 2/3`, reproduced via
  `scipy.optimize.brentq` on `RLDSTrajectory.R` at canonical `K=0.0106`
  → exact value 0.7144; see `constants.Z_ACC`).

## Paper V — microphysical descent

- **Figures 1–4 (microphysical validation):**
  `scripts/paper_v/fig5_generate.py`
  - Run: `python scripts/paper_v/fig5_generate.py --out outputs/paper_v`
  - See `scripts/paper_v/README_PAPER_V_REPRODUCTION.md` for outputs and
    expected key values.
- **Figure 5 (spectral-gap / Fenichel scale separation, `δ_max`):**
  `rlds_paired.core.spectral.compute_spectral_gap()`
  - At the canonical benchmark, `δ_max ≈ 0.244` at `z ≈ 0.40`; the Fenichel
    bound (0.10) is reported **violated** (the reduction is a ~24% correction
    in the transition window).

> **Note.** Paper V Figure 5 (the spectral gap / `δ_max`) is reproduced by the
> core library through `rlds_paired.core.spectral`, **not** by
> `scripts/paper_v/fig5_generate.py`. The script `fig5_generate.py` is
> responsible only for the Paper V **microphysical Figures 1–4**.

## Model comparison (`Δχ²`, `ΔBIC`)

- **Reference values:** `data/reference/model_comparison_reference.json`
  - The DR2 comparison reports the RLDS shooting form disfavored vs ΛCDM:
    `Δχ² ≈ +36`, `ΔBIC ≈ +58` (sign convention `Δ = RLDS − LCDM`).
  - These are **stored reference values, not recomputed by the SDK**: the
    comparison is supernova-inclusive and needs the external Pantheon+ /
    Pantheon+SH0ES dataset, which is not redistributed here. The JSON makes the
    claim traceable, not re-derivable from the archive alone. See the README
    section "Model-comparison reference values and external data policy".

## Quick command reference

```bash
# Papers III/IV headline numbers (recomputed):
python -m rlds_paired.reference

# Paper III K-scan / paired coefficient:
rlds-paired scan --n 40

# Paper V Figure 5 (spectral gap):
rlds-paired spectral --table

# Paper V Figures 1–4 (microphysical validation):
python scripts/paper_v/fig5_generate.py --out outputs/paper_v
```
