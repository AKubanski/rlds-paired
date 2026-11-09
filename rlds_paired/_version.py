"""
Single source of truth for the package version.

This module is intentionally tiny and dependency-free: it imports nothing.
That lets tooling (release workflows, build backends, packagers) read the
version with a simple regex or a lightweight import, without pulling in the
runtime dependencies (numpy, scipy, matplotlib).

``pyproject.toml`` reads ``__version__`` from here via
``[tool.setuptools.dynamic]``, and ``rlds_paired/__init__.py`` re-exports it.
"""

__version__ = "1.1.0"
