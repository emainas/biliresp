# biliresp

<p>
  <img src="docs/img/profile.png" alt="Electrostatic potential for biliverdin" width="200">
</p>

<p align="center">
  <a href="https://github.com/emainas/biliresp/actions">
    <img src="https://img.shields.io/badge/status-alpha-orange" alt="Status: alpha">
  </a>
</p>

> **Status:** pre-release, under active development. Interfaces may change without notice.

Utilities for parsing electrostatic potential output (In Terachem this is included in `resp.out`) and ESP grid (In Terachem this is outputed as `esp.xyz`) files. The package supplies:

- A parser (`resp.ParseRespDotOut`) for extracting RESP frames and ESP grids from an ab initio Molecular Dynamics trajectory or QM/MM trajectory (or a single conformer calculation can be used).
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
