# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**polartools** is a Python package for reading and processing data from the APS Polar beamline (4-ID-D at Argonne National Laboratory). It handles X-ray absorption (XAS/XMCD), diffraction, and spectroscopic data analysis. See [ARCHITECTURE.md](ARCHITECTURE.md) for the module map, catalog-backend details, and test-data fixtures.

## Commands

### Install (development)
```bash
pip install -e .[dev]
```

### Run tests
```bash
# Full suite with coverage (config lives in pyproject.toml, not .coveragerc):
pytest --cov=polartools --cov-report=xml

# Single file / single test:
pytest polartools/tests/test_absorption.py
pytest polartools/tests/test_absorption.py::test_function_name
```
No manual `databroker-unpack` step is needed — `polartools/tests/conftest.py` registers the packed test catalog automatically on import.

### Lint and format
Uses `ruff` (not flake8/black) for both linting and formatting; config in `pyproject.toml` (`[tool.ruff]`).
```bash
ruff check .
ruff format --check .   # add --diff to preview, drop --check to apply
```

## Versioning
Version is derived from git tags via `hatch-vcs` (see `[tool.hatch.version]` in `pyproject.toml`). Do not manually edit `polartools/_version.py`.

## Known issues
Don't rely on a hardcoded list here — it goes stale fast. Check current open issues with:
```bash
gh issue list --state open
```
