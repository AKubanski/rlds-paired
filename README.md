# rlds-paired

**Regime-Limited Dynamical Systems — Paired Dark-Sector Signature**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org)
[![CI](https://github.com/Askuba77/rlds-paired/actions/workflows/ci.yml/badge.svg)](https://github.com/Askuba77/rlds-paired/actions/workflows/ci.yml)

A Python library for reproducing the numerical results of the RLDS
cosmology series (Papers III–V).  `rlds-paired` is the observational–cosmological
branch of the RLDS SDK family.

The central result is the **paired growth–expansion signature**:

```
Δ(f·σ₈) / (f·σ₈)_ΛCDM  =  c · K + O(K²)
c  ≈  −0.15      (negative: growth suppression)
```

The coefficient is **computed from the model code, not stored as a literal** —
recompute it with `rlds_paired.reference.reference_results()`.

## Theory references

This library implements the methods from:

1. Kubanski, A. (2026). *Paired observational signature of the RLDS dark sector:
   growth suppression and phantom excursion.*
   Zenodo. [doi:10.5281/zenodo.19636371](https://doi.org/10.5281/zenodo.19636371) [Paper III]

2. Kubanski, A. (2026). *Cosmological consequences under RLDS closure.*
   Zenodo. [doi:10.5281/zenodo.19640061](https://doi.org/10.5281/zenodo.19640061) [Paper IV]

3. Kubanski, A. (2026). *Microphysical descent to RLDS.*
   Zenodo (forthcoming). [Paper V]

The mathematical core (Papers I–II) is implemented in
[rlds-mat](https://github.com/AKubanski/rlds-mat).

## Features

- **Algebraic IDE model**: E²(z) with tanh-transition, self-consistent R(z),
  phantom w_eff(z).
- **Full ODE trajectory**: Coupled [R, u] two-phase system (Paper IV), global
  attractor R* = 0.109.
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
# {'K': 0.012, 'alpha': 3.01, 'Rc': 0.51, 'Omega_m': 0.296, ...}

# ODE trajectory (Paper IV)
traj = rlds_paired.RLDSTrajectory()
print(f"R* = {traj.R_star:.4f}")        # 0.1095  (attractor, sqrt(K))
print(f"R(z=0) = {traj.R_today:.4f}")   # 0.2969  (present-day matter fraction)
print(f"E(z=0.5) = {traj.E(0.5):.4f}")

# Spectral gap validation (Paper V, Appendix C)
sg = rlds_paired.compute_spectral_gap()
print(sg.summary())
# delta_max ≈ 0.245 at z ≈ 0.42  →  Fenichel bound (0.10) VIOLATED (expected)

# Recomputed headline results (paired coefficient, delta_max)
print(rlds_paired.reference_results().summary())

# Growth factor (Paper III)
fs8 = rlds_paired.growth_factor(rlds_paired.OMEGA_M, K=0.012, z_eval=0.5)
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

## Canonical benchmark values

| Quantity | Value | Source |
|---|---|---|
| Coupling *K* | 0.012 | Paper IV, shooting MCMC |
| Dilution rate *α* | 3.010 | Paper IV |
| Switching threshold *R_c* | 0.51 | Paper IV |
| Switching width *ΔR* | 0.15 | Paper IV, shooting best-fit |
| Derived *Ω_m* | 0.296 | Paper IV, self-consistency |
| Attractor *R** | 0.109 (= √K) | Paper IV |
| Present-day matter fraction *R(z=0)* | 0.297 | computed |
| Equation of state *w_DE* | −1.003 | Paper IV |
| Acceleration onset *z_acc* | 0.705 | Paper IV |
| Paired coefficient *c* | ≈ −0.15 (computed) | Paper III |
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
