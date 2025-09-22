# Getting Started

## Prerequisites

- Python 3.9 or newer.
- A virtual environment to isolate dependencies is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
```

## Install the package

Install the project in editable mode so that local code changes are reflected immediately:

```bash
python -m pip install -e .
```

## Run the test suite

Pytest expects the `src/` directory on `PYTHONPATH` so the package can be located:

```bash
PYTHONPATH=src pytest -s tests
```

Use `-k` to narrow to a single test module when iterating, for example:

```bash
PYTHONPATH=src pytest -s tests/test_dipole.py
```

