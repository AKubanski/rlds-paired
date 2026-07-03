"""
Custom exceptions for rlds-paired.

All public exceptions inherit from :class:`RLDSPairedError` so that callers
can catch the entire family with a single ``except RLDSPairedError`` clause.
"""

from __future__ import annotations


class RLDSPairedError(Exception):
    """Base exception for all rlds-paired errors."""


class ConvergenceError(RLDSPairedError):
    """Raised when a numerical solver fails to converge."""


class PhysicsError(RLDSPairedError):
    """Raised when physical constraints are violated (e.g. K < 0)."""


class ConfigurationError(RLDSPairedError):
    """Raised when parameters are inconsistent or out of range."""
