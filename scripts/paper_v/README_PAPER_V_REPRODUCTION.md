# Paper V — microphysical reproduction (Figures 1–4)

This directory holds the official reproduction script for the **microphysical
validation figures** of RLDS Paper V (*Microphysical descent to RLDS*).

The script `fig5_generate.py` regenerates Paper V Figures **1–4** from a
single-field thawing-scalar microphysical model and compares its background
matter-fraction trajectory `R(N)` against the reduced RLDS flow evaluated on
the **DESI DR2 canonical frame** (imported from `rlds_paired.constants` and the
`rlds_paired.core.ode.RLDSTrajectory` solver). It is the Paper V counterpart to
the core-library reproduction paths used for Papers III and IV (see the
repository-level `REPRODUCE_FIGURES.md`).

> **Scope.** This script covers Paper V Figures **1–4** (microphysical
> validation). Paper V **Figure 5** (spectral-gap / `δ_max`) is reproduced by
> the core library via `rlds_paired.core.spectral.compute_spectral_gap()`, not
> by this script.

## What it reproduces

| Figure | File | Content |
|---|---|---|
| Fig. 1 | `fig5_trajectory.{pdf,png}` | `R_micro(N)` vs reduced `R_RLDS(N)`, with residual `ΔR` |
| Fig. 2 | `fig5_restoring.{pdf,png}` | reconstructed micro current vs the RLDS restoring term (log-y) |
| Fig. 3 | `fig5_attractor.{pdf,png}` | late-time convergence `\|R−R*\|` and the `γ` relaxation guide |
| Fig. 4 | `fig5_switching.{pdf,png}` | projected switching profile and coupling activation `β_eff(N)` |

Each figure ships a matching CSV of the underlying curve data, and a single
machine-readable `fig5_summary.json` of the key numbers.

**Honest status (unchanged from the paper).** The script reproduces the RLDS
background trajectory `R(N)` and the structural form of the switching, but the
reconstructed microscopic current `Q_micro` sits orders of magnitude below the
RLDS restoring amplitude — the micro–meso gap discussed in Paper V (Section 9 /
10, Prop. 10.1, and the Fig. 2 caption). The figures are a *structural*
validation, not a claim of amplitude-level matching.

## Dependencies

- `numpy`, `scipy`, `matplotlib` (all are already core dependencies of
  `rlds-paired`, so a normal install is sufficient)
- the `rlds_paired` package itself, importable on `PYTHONPATH`

```bash
pip install -e .          # from the repository root; pulls numpy/scipy/matplotlib
```

## Run

From the repository root:

```bash
python scripts/paper_v/fig5_generate.py --out outputs/paper_v
```

`--out` is optional; it defaults to `outputs/paper_v` (relative to the current
working directory) and the directory is created automatically. There is no
hardcoded absolute output path.

## Expected output files

In the chosen output directory:

```
fig5_trajectory.pdf      fig5_trajectory.png
fig5_restoring.pdf       fig5_restoring.png
fig5_attractor.pdf       fig5_attractor.png
fig5_switching.pdf       fig5_switching.png
fig5_data_trajectory.csv
fig5_data_restoring.csv
fig5_data_convergence.csv
fig5_data_switching.csv
fig5_summary.json
```

(4 PDFs, 4 PNGs, 4 CSVs, 1 JSON.)

## Expected key values (DR2 canon)

These are the values printed to the console and written to `fig5_summary.json`:

| Key | Expected |
|---|---|
| `lambda` | `0.5` |
| `beta_star` | `6e-4` |
| `dphi` | `0.3` |
| `phi_c` | `0.85` |
| `V0` | ≈ `1.886e-5` |
| `R0` (= R(z=0)) | ≈ `0.2974` |
| `R_star` | ≈ `0.1027` |
| `gamma_RLDS` | ≈ `5.40` |
| `max_abs_delta_R_z_0_10` | ≈ `0.0070` (< 0.01) |
| `beta_eff_peak` | ≈ `0.0010` |
| `beta_eff_peak_N` | ≈ `1.62` |
| `beta_eff_peak_R` | ≈ `0.0046` |
| `Q_rec_peak` | ≈ `5e-6` |
| `RLDS_term_max` | ≈ `3.15` |
| `micro_meso_gap` | ≈ `O(100)` (≈ 123) |

`V0` is obtained by shooting so that the microphysical `R(z=0) = 0.2974` (DR2);
the RLDS branch is taken from the `rlds-paired` SDK at the DR2 canonical
parameters (`K=0.0106`, `α=3.0102`, `R_c=0.565`, `ΔR=0.170`, `Ω_m=0.297`).

## JSON key ↔ Paper V symbol

The keys in `fig5_summary.json` use ASCII identifiers; the worked example in
Paper V (Section 6) uses mathematical notation. The mapping is one-to-one:

| JSON key | Paper V symbol |
|---|---|
| `lambda` | λ |
| `beta_star` | β_* |
| `dphi` | Δφ |
| `phi_c` | φ_c |
| `V0` | V₀ |
| `R0` | R(z=0) |
| `R_star` | R* |
| `gamma_RLDS` | γ_RLDS |
| `max_abs_delta_R_z_0_10` | max\|ΔR\| over z ∈ [0, 10] |
| `beta_eff_peak` | β_eff (peak value) |
| `beta_eff_peak_N` | N at the β_eff peak |
| `beta_eff_peak_R` | R at the β_eff peak |
| `Q_rec_peak` | Q_rec (peak value) |
| `RLDS_term_max` | max of the RLDS restoring term |
| `micro_meso_gap` | R_c^eff / R_c^inst |

This lets an SDK↔TEX consistency check match field names to symbols
mechanically, without inferring the λ → `lambda` style transliteration.

## Reproducibility notes

- The recipe constants (`lambda`, `beta_*`, `dphi`, `phi_c`) and the `V0`
  shooting target are **not** changed by this packaging; they are identical to
  the original generation script.
- The RLDS branch is always recomputed from `rlds_paired` at import time, so the
  Paper V figures track the same canonical frame as Papers III–IV. No headline
  number is read from a stored literal here.
