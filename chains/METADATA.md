# Archived chain — metadata

Canonical shooting MCMC for the RLDS cosmology series (Papers III–V),
on **DESI DR2** BAO (full covariance). The six `.npy` files are archived in
this directory.

## Sampler

- **Engine:** `emcee` (affine-invariant ensemble; stretch move).
- **Walkers × steps:** 32 × 3000.  **Burn-in:** 200.  **Flat samples:** 89 600.
- **Acceptance:** 0.366.
- **RNG seed:** *not fixed* in the reference run → only *statistical*
  reproduction from scratch. These archived arrays are the canonical artefact;
  use `rlds-paired mcmc --run --seed N` for a fresh, seeded chain.
- **ODE regularisation:** `eps = 1e-6` (matches `rlds_paired.mcmc.EPS`; the
  maintained `analysis` branch uses `1e-12`, a sub-σ difference).

## Parameters (5 sampled)

| col | name | prior | posterior median |
|----:|------|-------|------------------|
| 0 | `K`        | [0.005, 0.25] | 0.0106 |
| 1 | `H0`       | [60, 75]      | 67.76  |
| 2 | `R_c`      | [0.30, 0.60]  | 0.565  |
| 3 | `Delta_R`  | [0.03, 0.20]  | 0.170  |
| 4 | `alpha`    | [3.01, 3.30]  | 3.0102 |

`Omega_m` is **derived** (shooting: `R(z=0) = Omega_m`), not sampled —
posterior median **0.297 ± 0.007**.

The medians above are the locked **canonical set** (Rc=0.565, ΔR=0.170,
Ω_m=0.297) used throughout Papers III/IV/V. Attractor `R* = √K ≈ 0.103`.

## Honest note on the phantom

`alpha` is **100 % railed against the prior floor 3.01** (min 3.0100,
median 3.0102, max 3.0127) → `w_DE ≈ −1.003`. The data pull toward `w = −1`;
the phantom excursion is a consequence of the prior floor, not a data
preference. This is the "prior-dependent phantom" of Paper IV.

## Best fit (argmax log-prob)

`theta_best = [K=0.00963, H0=67.790, R_c=0.600, ΔR=0.200, α=3.0100]`,
`chi2 = 750.95`. Note this point sits at the prior-box edges and differs from
the medians; `rlds-paired mcmc` recomputes this χ² as an exact regression check.

## Likelihood

Pantheon+SH0ES (diagonal, analytic M-marginalisation) ⊕ DESI DR2 BAO
(13-point combined data vector, full block covariance) ⊕ Planck 2018 distance
prior. See `SOURCES.md`.

## Files

| file | shape | contents |
|------|-------|----------|
| `chains_raw.npy`     | (3000, 32, 5) | full chain (steps, walkers, params) |
| `logprob_raw.npy`    | (3000, 32)    | log-posterior per step/walker |
| `flat_samples.npy`   | (89600, 5)    | post-burn-in, flattened |
| `flat_logprob.npy`   | (89600,)      | log-posterior of flat samples |
| `omega_m_derived.npy`| (2000,)       | derived `Omega_m` per (sub)sample |
| `summary.npy`        | () dict       | best-fit, χ², acceptance, settings |

## Reproduce

```bash
pip install "rlds-paired[mcmc]"
export RLDS_DATA_DIR=/path/to/cosmology_data    # holds Pantheon+SH0ES.dat
rlds-paired mcmc                 # checks chi2 at the archived best-fit (== 750.95)
rlds-paired mcmc --run --seed 0  # run a fresh seeded chain
```
