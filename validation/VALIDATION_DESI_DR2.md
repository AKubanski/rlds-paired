# Stage 1.9 — validation against DESI DR2 (arXiv:2503.14738v2, Tables V–VI)

Ours: best-fit (MAP) with `likelihood_v2` (Planck distance priors R, l_A, omega_b;
DESI DR2 BAO full cov; Pantheon+SH0ES full STAT+SYS cov). DESI: posterior means ± 68 %
with the full CMB likelihood (Table V).

| Model | Data | Param | DESI DR2 | Ours | Diff / σ |
|---|---|---|---|---|---|
| ΛCDM | DESI+CMB | Ωm | 0.3027 ± 0.0036 | 0.3000 | −0.75 |
| ΛCDM | DESI+CMB | H0 | 68.17 ± 0.28 | 68.48 | +1.1 |
| wCDM | DESI+CMB | w | −1.055 ± 0.036 | −1.026 | +0.8 |
| wCDM | DESI+CMB | Ωm | 0.2927 ± 0.0073 | 0.2957 | +0.4 |
| wCDM | DESI+CMB+Pantheon+ | w | −0.995 ± 0.023 | −0.977 | +0.8 |
| wCDM | DESI+CMB+Pantheon+ | Ωm | 0.3047 ± 0.0051 | 0.3047 | 0.0 |
| CPL | DESI+CMB | w0 | −0.42 ± 0.21 | −0.43 | −0.05 |
| CPL | DESI+CMB | wa | −1.75 ± 0.58 | −1.68 | +0.1 |
| CPL | DESI+CMB+Pantheon+ | w0 | −0.838 ± 0.055 | −0.852 | −0.25 |
| CPL | DESI+CMB+Pantheon+ | wa | −0.62 (+0.22/−0.19) | −0.525 | +0.45 |

Δχ²(CPL − ΛCDM), DESI+CMB: DESI −12.5 (full CMB), **−8.0 with compressed (θ*, ω_b, ω_bc) CMB**;
ours −8.1 (compressed distance priors) — matches the compressed-CMB value.
DESI+CMB+Pantheon+: DESI −10.7 (full CMB); ours −7.7. The ~30 % lower Δχ² is the expected
loss of information from a compressed CMB likelihood (same ratio as for DESI+CMB).

**Verdict: Stage 1.9 PASSED** — all parameters within 1.1 σ (best-fit vs posterior mean).
**Limitation to state in papers:** compressed CMB underestimates Δχ² for dynamical-DE models by ~30 %.
