# biliresp

Utilities for parsing RESP output (`resp.out`) and ESP grid (`esp.xyz`) files. The package supplies:

- A parser (`resp.ParseRespDotOut`) for extracting RESP frames and ESP grids.
- A linear ESP charge fitting implementation (`linearESPcharges.linear`).
- Dipole post-processing helpers (`dipole.three_dipoles_last_frame`).
- Command-line entry points in `scripts/` for quick comparisons.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

## Run the test suite

```bash
PYTHONPATH=src pytest -s tests
```

Use `-k` to narrow to an individual module while iterating (for example `-k test_dipole`).

## Command-line entry points

All commands assume the sample RESP outputs in `data/raw/` and 78 atoms; adjust to your system as needed.

```bash
# Compare RESP ESP charges to fitted charges
python scripts/compare_charges.py data/raw/resp.out data/raw/esp.xyz 78 --frame -1 --solver explicit

# Print QM, ESP, and fitted dipoles for a frame
python scripts/print_dipoles.py data/raw/resp.out data/raw/esp.xyz 78 --frame -1 --solver explicit
```

Both scripts accept `--help` for a summary of arguments. Swap `--solver kkt` to use the block KKT solver.

## Documentation site

MkDocs config lives in `mkdocs.yml`. Install MkDocs (and optionally `mkdocs-material`) and serve locally:

```bash
python -m pip install mkdocs mkdocs-material
mkdocs serve
```

Navigate to the printed URL to browse the generated documentation, which covers project setup and the two workflows above.

## Packaging and distribution

The project already exposes the package via `pyproject.toml`. Once you create a GitHub repository:

1. Commit the sources, tests, and documentation.
2. Push to GitHub.
3. (Optional) Enable GitHub Pages and publish the MkDocs site via `mkdocs gh-deploy`.

