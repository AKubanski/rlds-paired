"""
Basic usage examples for rlds-paired.

Run:
    python examples/basic_usage.py
"""

import rlds_paired


def main():
    print("=" * 60)
    print("rlds-paired  —  basic usage examples")
    print("=" * 60)

    # ── 1. Canonical benchmark values ──────────────────────────────
    print("\n1. Canonical benchmark values:")
    for k, v in rlds_paired.CANONICAL.items():
        print(f"   {k:15s} = {v}")

    # ── 2. ODE trajectory ──────────────────────────────────────────
    print("\n2. RLDS ODE trajectory (Paper IV):")
    traj = rlds_paired.RLDSTrajectory()
    print(f"   E(z=0)   = {traj.E(0.0):.6f}  (must be 1.0)")
    print(f"   E(z=0.5) = {traj.E(0.5):.4f}")
    print(f"   R*       = {traj.R_star:.4f}  (attractor, stable fixed point)")
    print(f"   R(z=0)   = {traj.R_today:.4f}  (present-day matter fraction)")

    # ── 3. Spectral gap ────────────────────────────────────────────
    print("\n3. Spectral gap validation (Paper V, Appendix C):")
    sg = rlds_paired.compute_spectral_gap()
    print(sg.summary())

    # ── 4. Growth factor ───────────────────────────────────────────
    print("4. Growth factor (Paper III):")
    fs8_lcdm = rlds_paired.growth_factor(rlds_paired.OMEGA_M, K=0.0, z_eval=0.5)
    fs8_rlds = rlds_paired.growth_factor(
        rlds_paired.OMEGA_M, K=rlds_paired.K_BEST, z_eval=0.5
    )
    print(f"   f·σ₈(z=0.5, K=0)     = {fs8_lcdm:.4f}  (ΛCDM limit)")
    print(f"   f·σ₈(z=0.5, K={rlds_paired.K_BEST}) = {fs8_rlds:.4f}")
    print(
        f"   Δ/f·σ₈               = {(fs8_rlds-fs8_lcdm)/fs8_lcdm:+.5f}  (single point)"
    )

    # ── 5. Recomputed headline results ─────────────────────────────
    print("\n5. Canonical headline results (recomputed from code):")
    print(rlds_paired.reference_results().summary())

    print("Done.")


if __name__ == "__main__":
    main()
