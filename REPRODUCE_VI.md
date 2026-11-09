# Reproducing RLDS VI (closing analysis against DESI DR2)

This file maps every number, table and figure of

> A. Kubanski, *RLDS VI: What the Expansion History Tests in a Regime-Limited
> Interacting Dark Sector — A Closing Analysis against DESI DR2*, Version 1
> (November 2026)

to the script that produces it and the output file it writes. All outputs listed
here are committed in `validation/`, so every number can be checked without
re-running anything; the commands regenerate them.

## Setup

```bash
pip install -e ".[dev,mcmc]" h5py   # numpy, scipy, matplotlib, emcee, pytest; h5py for the chains
pip install classy                  # only for the CLASS cross-checks (Sec. 3.1)
export RLDS_DATA_DIR=/path/to/data  # Pantheon+SH0ES.dat and Pantheon+SH0ES_STAT+SYS.cov
export PYTHONPATH=.:validation      # run every command from the repository root
python -m pytest -q                 # 55 tests (15 for likelihood_v2)
```

Data: Pantheon+SH0ES distances and STAT+SYS covariance (Scolnic et al. 2022,
public data release; SHA-256 of the `.dat` file used:
`1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8`), DESI DR2 BAO
(bundled in `rlds_paired/desi_dr2/`), Planck 2018 distance priors of Chen,
Huang & Wang (2019) (hard-coded in `rlds_paired/likelihood_v2.py`).

Run times are for one CPU core. The final runs and checks used Python 3.13, numpy 2.5,
scipy 1.18, matplotlib 3.11 and emcee 3.1.6.

## Likelihood (Sec. 3)

| Paper item | Script / module | Output |
|---|---|---|
| Corrected likelihood (CMB R, l_A, omega_b with early omega_m; Aubourg r_d; full SN covariance) | `rlds_paired/likelihood_v2.py` | — |
| Model constructors (LCDM, wCDM, CPL, RLDS free, RLDS shooting) | `rlds_paired/likelihood_v2.py` | — |
| Tests 1.1–1.8 (radiation, omega_b, distances, r_s, z*, CMB vector, r_d, limits, comparison principle, SN offset invariance, regression of the old pipeline) | `tests/test_likelihood_v2.py` | pytest |
| Distances and r_d against CLASS (Sec. 3.1) | `validation/check_vs_class.py` | stdout |
| l_A calibration factor at four LCDM control cosmologies | `validation/lA_convention.py` | stdout |
| l_A calibration at the Planck-mean LCDM point and four wCDM / CPL / phantom cosmologies | `validation/lA_convention_de.py` | stdout (spread 5.6e-5 over the five points) |

## Validation against DESI DR2 (Sec. 3.3, Table 1)

| Paper item | Command | Output |
|---|---|---|
| Best fits of LCDM, wCDM, CPL for DESI+CMB and DESI+CMB+Pantheon+ | `python validation/fit_models.py LCDM,wCDM,CPL` (~10 min) | `validation/fit_LCDM_wCDM_CPL.json` |
| Table 1 with differences in sigma | — | `validation/VALIDATION_DESI_DR2.md` |
| Control MCMC chains (LCDM, wCDM, CPL; >= 50 tau) | `python validation/run_mcmc.py MODEL N_STEPS`, then `python validation/analyze_chains.py LCDM wCDM CPL` | `validation/chains/*.h5`, `validation/MCMC_RUN1.txt` |

Note: `fit_models.py` uses `H0 in [55, 80]`; all optima are interior, so this does
not affect any number.

Note on the chains: `run_mcmc.py` resumes an existing chain file and adds N_STEPS to it.
The committed chains have 3000 (LCDM, wCDM), 4000 (CPL) and 11 000 (RLDS-free) steps, so
running the command on them extends them. To rerun a chain from scratch, move the committed
file away first, e.g. `mv validation/chains/CPL.h5 validation/chains/CPL_archived.h5`, then
run `python validation/run_mcmc.py CPL 4000`. The initial walker positions are seeded
(seed 20261104) but the emcee proposals are not, so a rerun agrees with the archived chain
statistically, not bit for bit.

## Model comparison (Sec. 4.1, Table 2)

| Row of Table 2 | Command | Output |
|---|---|---|
| wCDM, CPL (and LCDM chi2 = 1421.36) | `fit_models.py` (above) | `validation/fit_LCDM_wCDM_CPL.json` |
| RLDS shooting, pre-registered: -4.80 | `python validation/polish_shoot.py prereg` (~15 min), restarting from the minima of `profile_K.py cmb,bao,sn` (~1.5 h) and `controlled_comparison.py new` | `validation/polish_shoot_prereg.json` |
| RLDS shooting, R_c <= 0.8: -4.98 | `python validation/polish_shoot.py rc0.8` (~15 min), restarting from the minimum of `profile_K.py cmb,bao,sn --rc-max 0.8` | `validation/polish_shoot_rc0.8.json` |
| RLDS shooting, R_c <= 0.9, dR >= 0.01: -5.01 | `python validation/polish_shoot.py wide` (~10 min), restarting from the result of `fit_shoot_wide.py` (~25 min) | `validation/polish_shoot_wide.json` |
| RLDS free, pre-registered: -6.36 | `python validation/refit_free_prereg.py` (~10 min) | `validation/refit_free_prereg.json` |
| RLDS free, widened: -7.97 (four of six starts) | `python validation/fit_free_wide.py` (~40 min) | `validation/fit_free_wide.json` (`none`) |
| RLDS free, Om <= 0.60: -8.53 | `python validation/fit_free_wide_om.py` | `validation/fit_free_wide_om060.json` |
| RLDS free, widened, with f_gas prior: -0.80 | `fit_free_wide.py` (above) | `validation/fit_free_wide.json` (`fgas`) |
| MCMC of the free version (11 000 steps, not converged) | `python validation/run_mcmc.py RLDS-free 11000` | `validation/chains/RLDS-free.h5`, `validation/MCMC_RUN1.txt` |

The shooting likelihood is flat in K and dR near its minimum, which lies at a prior edge. Earlier runs
stopped slightly higher: -4.77 at the minimum of the K profile on its grid, -4.79 in the controlled
comparison (Table 4). `polish_shoot.py` restarts Nelder-Mead from the lowest of them with soft prior
walls until two consecutive restarts gain less than 2e-4 (Appendix A of the paper).

## Profile of K and degeneracies (Sec. 4.2, Figs. 1–2)

| Paper item | Command | Output |
|---|---|---|
| Profile of K, CMB+BAO+SN and CMB+BAO | `python validation/profile_K.py cmb,bao,sn` and `... cmb,bao` | `validation/profile_K_cmb_bao_sn.json`, `validation/profile_K_cmb_bao.json` |
| Profiles with f_gas prior (0.26 +/- 0.04, +/- 0.08) and R_c <= 0.8 | `profile_K.py cmb,bao,sn --om-prior 0.260,0.040` (and `0.260,0.080`; `--rc-max 0.8`) | `validation/profile_K_cmb_bao_sn_om0.26-0.04.json`, `..._om0.26-0.08.json`, `..._rc0.8.json` |
| Fig. 1 | `python validation/plot_profile_K.py` | `validation/profile_K.png` |
| Fig. 2 | `python validation/plot_paper_figs.py` | `validation/fig2_profile_K_priors.png` |

## Created matter dOm (Sec. 4.3, Eqs. 6–7, Fig. 3)

| Paper item | Command | Output |
|---|---|---|
| Derived quantities (Om_early, dOm, effective w0/wa) of all K-profile points | `python validation/derived_observables.py` | `validation/derived_observables.json` |
| Eq. (5) against numerical w_eff(0) at 80 points (max diff 5.0e-4); epsilon 1e-6 vs 1e-12 (dchi2 ~ 1e-6); shooting ridge at K = 0.01 | `python validation/checks_misc.py` | `validation/checks_misc.json` |
| Shooting dOm profile, no prior | `python validation/profile_dOm.py cmb,bao,sn` | `validation/profile_dOm_cmb_bao_sn.json` |
| Shooting dOm profile, f_gas prior: 0.018 +/- 0.028, <= 0.074 (one-sided 0.065) | `python validation/profile_dOm.py cmb,bao,sn --om-prior 0.260,0.040` | `validation/profile_dOm_cmb_bao_sn_om0.26-0.04.json` |
| Same with R_c <= 0.9, dR >= 0.01: 0.022 +/- 0.026, <= 0.074 | `profile_dOm.py ... --rc-max 0.9 --dr-min 0.01 --targets ...` (two runs) | `validation/profile_dOm_cmb_bao_sn_om0.26-0.04_rc0.9_dr0.01_t*.json` |
| Free-version dOm profile, no prior (up to 0.17, -8.5) | `python validation/profile_dOm_free.py` (~45 min) | `validation/profile_dOm_free_ommax0.6.json` |
| Free-version dOm profile, f_gas prior: 0.010 +/- 0.035, <= 0.080 (one-sided 0.068) | `python validation/profile_dOm_free.py --om-prior 0.260,0.040` | `validation/profile_dOm_free_om0.26-0.04_ommax0.6.json` |
| dR edge check (dR >= 0.002): change ~1e-3 | `python validation/check_dR_edge_free.py` | `validation/check_dR_edge_free.json` |
| Fig. 3 and the 95% limits drawn in it | `python validation/plot_paper_figs_v2.py` | `validation/fig3_dOm_w0wa_v2.png`, `validation/fig3_bounds.json` |

## Evolving-dark-energy signal (Sec. 4.4, Fig. 4)

| Paper item | Command | Output |
|---|---|---|
| Uniform scan of the widened priors (free 4000, shooting 2500 draws) | `python validation/scan_w0wa.py free 4000 1`, `... shooting 2500 2` | `validation/scan_w0wa_free_1.json`, `validation/scan_w0wa_shooting_2.json` |
| Likelihood of the scan points inside the CPL contours (six 68% hits: -5.9 ... +6) | `python validation/scan_hits_chi2.py` | `validation/scan_hits_chi2.json` |
| Fig. 4 | `python validation/plot_paper_figs_v2.py` | `validation/fig4_scan_w0wa_v2.png` |

## Controlled comparison (Appendix C, Table 4)

| Paper item | Command | Output |
|---|---|---|
| Old likelihood (Papers III–IV) with old and new priors | `python validation/controlled_comparison.py old` | `validation/controlled_old.json` |
| New likelihood with old and new priors | `python validation/controlled_comparison.py new` | `validation/controlled_new.json` |

The old likelihood is the unchanged `rlds_paired/mcmc.py`; the test
`test_old_pipeline_regression` checks that it still reproduces the archived
chi2 = 750.95.

## Figures

All figures are regenerated from the committed JSON outputs and chains; the PNG files themselves
are not tracked (`.gitignore`). The curves and the limits drawn in them are reproduced exactly; the
appearance of the PNG files (fonts, sizes) depends on the matplotlib version. `plot_derived.py` and
`plot_profile_K2.py` make diagnostic working plots that are not figures of the paper.

## Superseded intermediate files

`validation/SUMMARY_STAGE1.txt` and `validation/fit_RLDS-free_RLDS-shooting.json`
contain the first-stage RLDS fits (shooting -1.36 at a local minimum, free
-4.75). They are kept as a record; the values used in the paper are those listed
above. The canonical frame of `REPRODUCE_FIGURES.md` (K = 0.0106, alpha = 3.0102)
belongs to Papers III–V and is superseded by Paper VI.
