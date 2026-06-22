"""
rlds-paired — setup.py shim
===========================

All project metadata is in ``pyproject.toml`` (PEP 621). This file exists
solely as a fallback for environments or tooling that still expect a
``setup.py`` to be present (some IDEs, legacy CI configurations,
``python setup.py develop``, etc.).

For modern usage:

    pip install .              # build and install from pyproject.toml
    pip install -e .           # editable install
    pip install -e ".[dev]"    # editable install with dev extras
    python -m build            # build sdist + wheel
"""

from setuptools import setup

setup()
