"""Shared pytest fixtures for rlds-paired tests."""

import pytest


@pytest.fixture(scope="session")
def canonical_traj():
    """RLDS trajectory at canonical parameters — computed once per session."""
    from rlds_paired.core.ode import RLDSTrajectory

    return RLDSTrajectory()


@pytest.fixture(scope="session")
def canonical_sg():
    """Spectral gap at canonical parameters (default K=K_BEST=0.0106) — computed once per session."""
    from rlds_paired.core.spectral import compute_spectral_gap

    return compute_spectral_gap()
