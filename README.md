# rlds-paired

**Regime-Limited Dynamical Systems — Paired Dark-Sector Signature**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org)
[![CI](https://github.com/AKubanski/rlds-paired/actions/workflows/ci.yml/badge.svg)](https://github.com/AKubanski/rlds-paired/actions/workflows/ci.yml)

A Python library for reproducing the numerical results of the RLDS
cosmology series (Papers III–V).  `rlds-paired` is the observational–cosmological
branch of the RLDS SDK family.

The central result is the **paired growth–expansion signature**:

```
Δ(f·σ₈) / (f·σ₈)_ΛCDM  =  c · K + O(K²)
c  ≈  −0.19      (negative: growth suppression)
```

The coefficient is **computed from the model code, not stored as a literal** —
recompute it with `rlds_paired.reference.reference_results()`.

> **Status after Paper VI (November 2026).** The observational conclusions of
> Papers III–IV — the shooting fit worse than ΛCDM by Δχ² ≈ +36, the "canonical"
> frame `K = 0.0106`, `α = 3.0102`, `Ω_m = 0.297`, and the paired coefficient as a
> physical prediction — are superseded by **RLDS VI**, the closing analysis
> against DESI DR2 with a corrected, validated likelihood
> (`rlds_paired/likelihood_v2.py`). The legacy likelihood in
> `rlds_paired/mcmc.py` and the canonical frame below are kept unchanged for
> reproducibility of Papers III–V. See [`REPRODUCE_VI.md`](REPRODUCE_VI.md)
> for every number, table and figure of Paper VI.

## Theory references

This library implements the methods from:

1. Kubanski, A. (2026). *Paired observational signature of the RLDS dark sector:
   growth suppression and phantom excursion.*
   Zenodo. [doi:10.5281/zenodo.19636371](https://doi.org/10.5281/zenodo.19636371) [Paper III]

2. Kubanski, A. (2026). *Cosmological consequences under RLDS closure.*
   Zenodo. [doi:10.5281/zenodo.19640061](https://doi.org/10.5281/zenodo.19640061) [Paper IV]

3. Kubanski, A. (2026). *Microphysical descent to RLDS.*
   Zenodo. [doi:10.5281/zenodo.19640630](https://doi.org/10.5281/zenodo.19640630) [Paper V]

4. Kubanski, A. (2026). *RLDS VI: What the Expansion History Tests in a
   Regime-Limited Interacting Dark Sector — A Closing Analysis against DESI DR2.*
   Zenodo (DOI to be added on publication). [Paper VI]

The mathematical core (Papers I–II) is implemented in
[rlds-mat](https://github.com/AKubanski/rlds-mat)
([doi:10.5281/zenodo.20342965](https://doi.org/10.5281/zenodo.20342965)).

## Features

- **Algebraic IDE model**: E²(z) with tanh-transition, self-consistent R(z),
  phantom w_eff(z).
- **Full ODE trajectory**: Coupled [R, u] two-phase system (Paper IV), global
  attractor R* ≈ 0.103.
- **Growth factor**: Linear growth f·σ₈(z) on the full coupled-ODE RLDS
  background (Paper III).
- **K-scan**: Paired-signature coefficient scan; linear scaling with a
  **negative** slope (growth suppression).
- **Spectral gap**: δ(N) = |dR/dN|/3 validation (Paper V, Appendix C). The
  Fenichel persistence bound (0.10) is reported pass/fail; at the canonical
  benchmark it is **violated** (δ_max ≈ 0.245), so the slow-manifold reduction
  is a marginal (~25%) correction in the transition window.
- **CLI**: `rlds-paired verify | scan | spectral | fsigma8`.
- **Canonical inputs**: Paper III–V benchmark *parameters* in one place;
  computed results live in `rlds_paired.reference`.

## Installation

```bash
pip install rlds-paired

# Development install
pip install -e ".[dev]"
```

## Quick start

### From Python

```python
import rlds_paired

# Canonical benchmark *inputs* (computed results live in rlds_paired.reference)
print(rlds_paired.CANONICAL)
# {'K': 0.0106, 'alpha': 3.0102, 'Rc': 0.565, 'Omega_m': 0.297, ...}

# ODE trajectory (Paper IV)
traj = rlds_paired.RLDSTrajectory()
print(f"R* = {traj.R_star:.4f}")        # 0.1027  (attractor, ≈ sqrt(K))
print(f"R(z=0) = {traj.R_today:.4f}")   # 0.2974  (present-day matter fraction)
print(f"E(z=0.5) = {traj.E(0.5):.4f}")

# Spectral gap validation (Paper V, Appendix C)
sg = rlds_paired.compute_spectral_gap()
print(sg.summary())
# delta_max ≈ 0.245 at z ≈ 0.42  →  Fenichel bound (0.10) VIOLATED (expected)

# Recomputed headline results (paired coefficient, delta_max)
print(rlds_paired.reference_results().summary())

# Growth factor (Paper III)
fs8 = rlds_paired.growth_factor(rlds_paired.OMEGA_M, K=rlds_paired.K_BEST, z_eval=0.5)
print(f"f·σ₈(z=0.5) = {fs8:.4f}")

# K-scan
result = rlds_paired.k_scan(n_points=20)
print(result.summary())

# f·σ₈ comparison figure
fs8_res = rlds_paired.fsigma8_comparison()
fs8_res.save_figure("fsigma8.png")
```

### From the command line

```bash
# Verify installation
rlds-paired verify

# Spectral gap table (Paper V, Appendix C)
rlds-paired spectral --table

# K-scan (20 points)
rlds-paired scan --n 20

# f·σ₈ comparison figure
rlds-paired fsigma8 --figure fsigma8.png
```

## Reproducing the figures

A full map from each paper's figures/results to the exact entry point is in
[`REPRODUCE_FIGURES.md`](REPRODUCE_FIGURES.md). In short: Paper III's paired
coefficient and K-scan come from `rlds_paired.analysis.k_scan()`; Paper IV's
`R*`, `z_acc`, and `w_eff` from `rlds_paired.core.ode` /
`rlds_paired.core.algebraic`; and Paper V's spectral gap (`δ_max`, Figure 5)
from `rlds_paired.core.spectral.compute_spectral_gap()`.

### Paper V microphysical reproduction

Paper V's **microphysical validation figures (Figures 1–4)** are regenerated by
a dedicated script bundled in the source archive:

```bash
python scripts/paper_v/fig5_generate.py --out outputs/paper_v
```

It writes four figures (PDF + PNG), four CSVs, and a machine-readable
`fig5_summary.json`. The recipe (`λ=0.5`, `β_*=6e-4`, `Δφ=0.3`, `φ_c=0.85`,
`V₀≈1.886e-5` by shooting) and the RLDS branch are evaluated on the DR2
canonical frame. See
[`scripts/paper_v/README_PAPER_V_REPRODUCTION.md`](scripts/paper_v/README_PAPER_V_REPRODUCTION.md)
for outputs and expected key values. (Paper V's Figure 5 / `δ_max` is reproduced
by `rlds_paired.core.spectral`, not by this script.)

## Model-comparison reference values and external data policy

The DR2 model comparison reports the RLDS shooting form as disfavored relative
to ΛCDM, with `Δχ² ≈ +36` and `ΔBIC ≈ +58`. The sign convention is:

```
Δχ²  = χ²_RLDS  − χ²_LCDM
ΔBIC = BIC_RLDS − BIC_LCDM
```

so **positive `Δχ²` / `ΔBIC` means the RLDS fit is disfavored** relative to
ΛCDM for the corresponding comparison.

The **raw Pantheon+ / Pantheon+SH0ES dataset is not redistributed** in this
archive: it is an external third-party dataset. A full likelihood rerun requires
the user to obtain it from the official source and place it under the documented
local data path. What this archive *does* include from the original run is the
set of derived artifacts: archived chains, summary arrays,
`chains/omega_m_derived.npy`, and a reference model-comparison file:

```
data/reference/model_comparison_reference.json
```

These make the paper-level `Δχ²` / `ΔBIC` claims **traceable**, not fully
re-derivable from the archive alone — the comparison is supernova-inclusive, so
the likelihood-level values are stored reference numbers from the papers rather
than quantities the SDK recomputes at run time. The model-level results (paired
coefficient, spectral gap, attractor, present-day matter fraction) remain fully
reproducible from the code, as described above.

## Canonical benchmark values

| Quantity | Value | Source |
|---|---|---|
| Coupling *K* | 0.0106 | Paper IV, DR2 shooting MCMC |
| Dilution rate *α* | 3.0102 | Paper IV (DR2; railed at prior floor 3.01) |
| Switching threshold *R_c* | 0.565 | Paper IV |
| Switching width *ΔR* | 0.170 | Paper IV, DR2 shooting best-fit |
| Derived *Ω_m* | 0.297 | Paper IV, self-consistency |
| Attractor *R** | 0.1027 (≈ √K) | Paper IV |
| Present-day matter fraction *R(z=0)* | 0.2974 | computed |
| Equation of state *w_DE* | −1.003 | Paper IV |
| Acceleration onset *z_acc* | 0.7144 | Paper IV, *q(z)=0* on DR2 background |
| Paired coefficient *c* | ≈ −0.1936 ± 0.005 (computed; r² ≈ 0.97) | Paper III |
| Spectral gap *δ_max* | ≈ 0.245 — Fenichel (0.10) violated | Paper V, Appendix C |

## Package structure

```
rlds_paired/
├── _version.py       ← single source of truth for version
├── constants.py      ← all canonical benchmark values
├── exceptions.py     ← custom exception hierarchy
├── analysis.py       ← k_scan(), fsigma8_comparison()
├── cli.py            ← rlds-paired command-line tool
└── core/
    ├── algebraic.py  ← algebraic IDE model (E², R, F, w_eff)
    ├── ode.py        ← two-phase coupled [R, u] ODE (Paper IV)
    ├── growth.py     ← linear growth factor (Paper III)
    └── spectral.py   ← spectral gap δ(N) = |dR/dN|/3 (Paper V)
```

## SDK family

| Library | Papers | Description |
|---|---|---|
| [rlds-mat](https://github.com/AKubanski/rlds-mat) | I–II | Mathematical core: ODE analysis, attractor proofs |
| **rlds-paired** | III–V | Observational: paired signature, growth, spectral gap |

## Author

**Aleksander Kubanski**
Independent Researcher, Poland
ORCID: [0009-0006-7239-5074](https://orcid.org/0009-0006-7239-5074)
Research site: [research.kubanski.pro](https://research.kubanski.pro)

## License

MIT License. See [LICENSE](LICENSE).
