# Contributing to rlds-paired

Thank you for your interest in contributing to `rlds-paired` — the
observational–cosmological SDK for the RLDS framework.

## Quick start

```bash
git clone https://github.com/AKubanski/rlds-paired.git
cd rlds-paired
pip install -e ".[dev]"
pytest tests/
```

If `pytest` reports all tests passing, your setup is correct.

## Scope of the library

`rlds-paired` implements the observational and cosmological results of the
RLDS framework as published in:

- *Paper III: Paired observational signature of the RLDS dark sector*
  (Zenodo, [doi:10.5281/zenodo.19636371](https://doi.org/10.5281/zenodo.19636371))
- *Paper IV: Cosmological consequences under RLDS closure*
  (Zenodo, [doi:10.5281/zenodo.19640061](https://doi.org/10.5281/zenodo.19640061))
- *Paper V: Microphysical descent to RLDS*
  (Zenodo, [doi:10.5281/zenodo.19640630](https://doi.org/10.5281/zenodo.19640630))

The library focuses on: the algebraic IDE expansion history, the two-phase
ODE trajectory, linear growth factor, paired-signature K-scan, and spectral-gap
validation.  Domain-general dynamical systems analysis belongs in
[rlds-mat](https://github.com/AKubanski/rlds-mat).

## Running the tests

```bash
pytest tests/ -v                    # all tests, verbose
pytest tests/test_spectral.py -v   # specific file
pytest tests/ --cov=rlds_paired    # with coverage
```

## Code style

- Formatter: `black` (line length 88)
- All public functions must have NumPy-style docstrings in English
- No Russian comments in source files
- Type hints on all public API functions

## Reporting issues

Please use the GitHub issue templates:
- **Bug report**: for unexpected behaviour or wrong numerical results
- **Feature request**: for new functionality

## Physics conventions

- N = ln(a): e-fold time (a = scale factor)
- R = ρ_DM / ρ_nr: matter fraction in the dark sector
- K: dimensionless coupling constant (K=0 → exact ΛCDM limit)
- All canonical parameters are in `rlds_paired/constants.py`
