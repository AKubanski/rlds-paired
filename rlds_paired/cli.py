"""CLI for rlds-paired.

Usage
-----
rlds-paired verify
rlds-paired scan [--n 20]
rlds-paired spectral [--table]
rlds-paired fsigma8 [--figure path]
"""

from __future__ import annotations
import argparse, sys


def _force_utf8_output():
    """Force stdout/stderr to UTF-8 so the check/Greek glyphs (✓ α δ ²) never
    crash on a Windows console or pipe whose locale codepage (cp1251/cp1252)
    cannot encode them. No-op where reconfigure is unavailable."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8")
        except (AttributeError, ValueError, OSError):
            pass


def _h(t):
    print(f"\n{'='*62}\n  {t}\n{'='*62}")


def _ok(m):
    print(f"  \033[92m✓\033[0m  {m}")


def _fail(m):
    print(f"  \033[91m✗\033[0m  {m}")


def cmd_verify(_args):
    _h("rlds-paired — Verification")
    import rlds_paired as rp
    from rlds_paired.reference import reference_results

    _ok(f"rlds-paired {rp.__version__}")
    _ok(f"K_BEST={rp.K_BEST}  ALPHA={rp.ALPHA}  RC={rp.RC}")

    t = rp.RLDSTrajectory()
    ok1 = abs(t.E(0) - 1.0) < 1e-6
    (_ok if ok1 else _fail)(f"E(z=0) = {t.E(0):.6f}  (must be 1.0)")

    ok2 = abs(t.R_star - rp.R_STAR) < 0.005
    (_ok if ok2 else _fail)(f"R* (attractor) = {t.R_star:.4f}  (expected {rp.R_STAR})")

    ok_today = 0.28 < t.R_today < 0.31
    (_ok if ok_today else _fail)(
        f"R(z=0) = {t.R_today:.4f}  (present-day matter fraction, expect ~{rp.OMEGA_M})"
    )

    ref = reference_results()

    # Fenichel is REPORTED, not gated: at the canonical benchmark it is violated by
    # design (delta_max ~ 0.244, peak in the transition window). What we *do* gate on
    # is that delta_max sits in the expected band and its peak is in the transition
    # (a guard against the old clipping bug, which pinned the peak at the z~3 edge).
    sg = rp.compute_spectral_gap()
    fen = "SATISFIED" if sg.fenichel_satisfied else "VIOLATED (expected)"
    _ok(
        f"δ_max = {sg.delta_max:.4f} at z = {sg.z_at_max:.3f}   Fenichel < {rp.FENICHEL_BOUND}: {fen}"
    )
    ok3 = (0.20 < sg.delta_max < 0.30) and (0.30 < sg.z_at_max < 0.55)
    (_ok if ok3 else _fail)(
        f"δ_max in band [0.20,0.30] and peak in transition: {'OK' if ok3 else 'FAIL'}"
    )

    ok4 = (ref.paired_coefficient < 0) and (-0.21 < ref.paired_coefficient < -0.17)
    (_ok if ok4 else _fail)(
        f"paired c = {ref.paired_coefficient:+.3f} ± {ref.paired_coefficient_err:.3f} "
        f"(r²={ref.paired_r2:.2f})  — negative, in band [-0.21,-0.17]: {'OK' if ok4 else 'FAIL'}"
    )

    print()
    if ok1 and ok2 and ok_today and ok3 and ok4:
        print("  All checks passed.\n")
        return 0
    print("  Some checks failed.\n")
    return 1


def cmd_scan(args):
    from rlds_paired.analysis import k_scan

    _h("K-scan — paired signature coefficient")
    result = k_scan(n_points=args.n, verbose=True)
    print()
    print(result.summary())
    if args.n < 40:
        print(
            f"\n  note: coarse grid (n={args.n}) — the converged value at the default "
            f"n≥40 is c ≈ -0.193 (matches `rlds-paired verify` and Papers III–V)."
        )
    result.save_figure("k_scan.png")
    return 0


def cmd_spectral(args):
    from rlds_paired.core.spectral import compute_spectral_gap, spectral_gap_table
    from rlds_paired.constants import ALPHA, FENICHEL_BOUND, K_BEST

    _h("Spectral gap validation — Paper V, Appendix C")
    if args.table:
        rows = spectral_gap_table()
        print(f"\n  {'K':>8} | {'αK':>8} | {'δ_max':>10} | Fenichel OK?")
        print("  " + "-" * 50)
        for r in rows:
            print(
                f"  {r.K:>8.3f} | {ALPHA*r.K:>8.4f} | {r.delta_max:>10.4f} | {'YES ✓' if r.fenichel_satisfied else 'NO ✗'}"
            )
    else:
        sg = compute_spectral_gap()
        print()
        print(sg.summary())
    return 0


def cmd_fsigma8(args):
    from rlds_paired.analysis import fsigma8_comparison

    _h("f*sigma8 comparison — RLDS vs ΛCDM")
    fs8 = fsigma8_comparison()
    path = args.figure or "fsigma8_comparison.png"
    fs8.save_figure(path)
    return 0


def cmd_mcmc(args):
    _h("Shooting MCMC — Papers III–V")
    from rlds_paired import mcmc
    from rlds_paired.data import load_pantheon_plus

    try:
        z, _, _ = load_pantheon_plus()
        _ok(f"Pantheon+ loaded: {len(z)} SNe")
    except FileNotFoundError as exc:
        _fail(
            f"Pantheon+ data not found ({exc}). "
            f"Set $RLDS_DATA_DIR or place Pantheon+SH0ES.dat in ~/cosmology_data."
        )
        return 1

    if args.run:
        try:
            import emcee  # noqa: F401
        except ImportError:
            _fail("`--run` needs emcee:  pip install 'rlds-paired[mcmc]'")
            return 1
        _ok(f"running MCMC: {args.walkers} walkers x {args.steps} steps (seed={args.seed})")
        _, summ = mcmc.run_mcmc(
            n_walkers=args.walkers, n_steps=args.steps, seed=args.seed,
            output_dir=args.out, fast_distance=True, progress=True,
        )
        _ok(
            f"done: acceptance={summ['acceptance']:.3f}, "
            f"chi2_best={summ['chi2_best']:.1f}, saved to {args.out}/"
        )
        return 0

    # default: reproduce the archived best-fit chi2 (fast, exact regression)
    chi2 = mcmc.reproduce_reference_chi2(fast_distance=False)
    ok = abs(chi2 - mcmc.REFERENCE_CHI2) < 0.01
    (_ok if ok else _fail)(
        f"reference chi2 = {chi2:.4f}  (archived {mcmc.REFERENCE_CHI2:.4f}, "
        f"Δ={chi2 - mcmc.REFERENCE_CHI2:+.4f})"
    )
    print()
    return 0 if ok else 1


def main(argv=None):
    _force_utf8_output()
    p = argparse.ArgumentParser(
        prog="rlds-paired",
        description="RLDS paired dark-sector signature — Papers III–V",
    )
    sub = p.add_subparsers(dest="command", metavar="COMMAND")
    sub.add_parser("verify", help="Check installation and canonical benchmarks.")
    ps = sub.add_parser("scan", help="K-scan of paired signature coefficient.")
    ps.add_argument("--n", type=int, default=40, help="Number of K points (default 40)")
    psp = sub.add_parser("spectral", help="Spectral gap validation (Paper V).")
    psp.add_argument("--table", action="store_true", help="Full table across K values.")
    pf = sub.add_parser("fsigma8", help="f*sigma8 comparison figure.")
    pf.add_argument("--figure", type=str, default=None)
    pm = sub.add_parser(
        "mcmc", help="Reproduce (default) or run the shooting MCMC (Papers III–V)."
    )
    pm.add_argument("--run", action="store_true", help="Run a new chain (needs emcee).")
    pm.add_argument("--steps", type=int, default=3000, help="MCMC steps (default 3000).")
    pm.add_argument("--walkers", type=int, default=32, help="Walkers (default 32).")
    pm.add_argument("--seed", type=int, default=None, help="RNG seed (default none).")
    pm.add_argument("--out", type=str, default="mcmc_output_shooting",
                    help="Output directory for chain artefacts.")
    args = p.parse_args(argv)
    if not args.command:
        p.print_help()
        return 0
    return {
        "verify": cmd_verify,
        "scan": cmd_scan,
        "spectral": cmd_spectral,
        "fsigma8": cmd_fsigma8,
        "mcmc": cmd_mcmc,
    }[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
