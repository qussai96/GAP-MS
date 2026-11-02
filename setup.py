"""Compatibility shim; v1 metadata is declared in pyproject.toml.

The gapms-v1 distribution depends on the separately installed legacy gapms 0.3.0
and does not ship or replace its command, model, or namespace.
"""
from setuptools import setup
setup()
