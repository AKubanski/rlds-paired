"""Tests for the two-phase RLDS ODE trajectory (core/ode.py)."""

import numpy as np
import pytest
from rlds_paired.constants import OMEGA_M, R_STAR


class TestRLDSTrajectoryBasic:
    def test_E_today(self, canonical_traj):
        assert abs(canonical_traj.E(0.0) - 1.0) < 1e-6

    def test_attractor(self, canonical_traj):
        assert abs(canonical_traj.R_star - R_STAR) < 0.005

    def test_present_day_R(self, canonical_traj):
        # R(z=0) is the present-day matter fraction (~Omega_m), distinct from
        # the attractor R* ~ sqrt(K). They coincided pre-fix only via the bug.
        assert 0.28 < canonical_traj.R_today < 0.31

    def test_E_monotone(self, canonical_traj):
        z = np.linspace(0.0, 2.0, 20)
        E = [canonical_traj.E(zi) for zi in z]
        assert all(E[i] <= E[i + 1] for i in range(len(E) - 1))

    def test_R_positive(self, canonical_traj):
        assert all(r > 0 for r in canonical_traj.R_arr)

    def test_R_below_one(self, canonical_traj):
        assert all(r < 1.0 for r in canonical_traj.R_arr)
