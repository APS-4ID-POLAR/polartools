# Installation

## From PyPI

```bash
pip install polartools
```

## From conda-forge

If you already manage your scientific Python stack with conda or mamba:

```bash
mamba install -c conda-forge polartools
```

## From source

For development, clone and install in editable mode with the dev extras:

```bash
git clone https://github.com/APS-4ID-POLAR/polartools
cd polartools
pip install -e .[dev,docs]
```

Or with conda for the heavier native dependencies (Qt, HDF5 plugins):

```bash
git clone https://github.com/APS-4ID-POLAR/polartools
cd polartools
mamba env create -f environment.yml
mamba activate polartools-dev
```

## Verifying

```bash
python -c "import polartools; print(polartools.__version__)"
```

Should print a version string like `0.5.5` (or `0.5.5.postN+gHASH` for a
development checkout).

## Optional: enable pre-commit hooks

If you'll be contributing, install the pre-commit hooks so ruff runs on every
commit:

```bash
pre-commit install
```
