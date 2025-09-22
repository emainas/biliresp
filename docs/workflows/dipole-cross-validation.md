# Dipole Cross-Validation

`dipole.three_dipoles_last_frame` compares three sets of dipole moments for the final frame of a RESP run:

1. **QM dipole** reported directly by Terachem.
2. **ESP unrestrained dipole** computed from RESP's own ESP charges.
3. **Optimized dipole** generated from your fitted charges.

The function parses the last "ESP unrestrained charges" block, the center of mass, and the dipole magnitude/vector from `resp.out`, then computes dipoles in Debye via:

```python
μ = Σ q_A (R_A - R_COM)
```

with coordinates expressed in bohr. It returns a dictionary containing the vectors, magnitudes, and the deltas relative to the QM reference.

## Script entry point

`scripts/print_dipoles.py` bundles the full workflow:

```bash
python scripts/print_dipoles.py data/raw/resp.out data/raw/esp.xyz 78 --frame -1 --solver explicit
```

- The script builds the linear system, fits charges (`explicit` or `kkt`), and prints the QM, ESP, and optimized dipole magnitudes along with delta vectors.
- `--frame` chooses which RESP frame to analyze (`-1` = last); `--solver` swaps between the two solvers exported from `linearESPcharges`.

## Typical output

```
QM (from resp.out log)
  vector (Debye): [-0.123 ...]
  |μ| (Debye):   2.345678

ESP unrestrained charges
  vector (Debye): [-0.121 ...]
  |μ| (Debye):   2.346100
  Δ vector vs QM (Debye): [0.002 ...]
  Δ|μ| vs QM (Debye):     0.000422

Fitted charges (explicit)
  vector (Debye): [-0.120 ...]
  |μ| (Debye):   2.345900
  Δ vector vs QM (Debye): [0.003 ...]
  Δ|μ| vs QM (Debye):     0.000222
```

Use this readout to sanity-check that your fitted charges reproduce the QM dipole within an acceptable tolerance before exporting them to downstream workflows.

## Programmatic use

```python
from dipole import three_dipoles_last_frame

# Assume you already prepared A, V, Q, resp_charges, coords_bohr via prepare_linear_system(..., return_positions=True)
dipoles = three_dipoles_last_frame("data/raw/resp.out", coords_bohr, fitted_result["q"])
print(dipoles["optimized_dipole_mag_D"], dipoles["delta_opt_vs_qm_mag_D"])
```

