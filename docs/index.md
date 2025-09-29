# biliresp

Utilities for analyzing electrostatic potential (ESP) outputs. The project includes a RESP parser for Terachem, a linear-ESP (raw, unrestrained) charge constraint optimization solver, examples for validating dipoles, and a mass-weighted center-of-mass calculator that can operate on RESP or xyz geometries.

Use the navigation to find quick-start installation instructions and focused walkthroughs for the two shipped command-line entry points:

- Fitting charges with the linear ESP solver.
- Cross-validating dipoles between RESP and fitted charges (and reconciling centers of mass) using the helper in `scripts/print_dipoles.py` or the reference implementation in `tests/test_dipole.py`.
