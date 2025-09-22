# Linear ESP Charge Solver

The linear solver fits electrostatic potential (ESP) charges that reproduce grid values exported by RESP/TeraChem. Everything lives in `linearESPcharges.linear` and is backed by numpy.

## Pipeline

1. **Parse Terachem RESP output** with `ParseRespDotOut` to obtain atomic positions and RESP ESP charges for each frame.
2. **Read ESP grid** points from `esp.xyz`.
3. **Build the design matrix** `A` where `A[i, j] = 1 / r_ij` for grid point `i` and atom `j`.
4. **Solve** the constrained optimization problem `A q ≈ V` subject to `Σ q = Q` using one of two solvers:
   - `explicit_solution`: closed-form projection that first finds the unconstrained least squares solution and then enforces the total charge (Lagrange multiplier method).
   - `KKTblock_solution`: solves the block matrix from the Karush–Kuhn–Tucker conditions in a single shot.

Both solvers accept an optional `ridge` hyper-parameter if you want to add a small diagonal Tikhonov term for numerical stability.

## Script entry point

The `scripts/compare_charges.py` wrapper prepares the system and prints per-atom differences between RESP unrestrained charges and the fitted charges:

```bash
python scripts/compare_charges.py data/raw/resp.out data/raw/esp.xyz 78 --frame -1 --solver explicit
```

- `resp.out`: RESP log file containing the ESP unrestrained block.
- `esp.xyz`: grid potentials from TeraChem.
- `78`: number of atoms in the RESP job.
- `--frame`: zero-based frame index (use `-1` for the last frame).
- `--solver`: either `explicit` (default) or `kkt`.

The output lists each atom index with the RESP charge, the fitted charge, and the difference. A footer prints charge conservation, RMSE, and RRMS metrics so you can quickly gauge the fit quality.

## Programmatic use

You can access the components directly:

```python
from linearESPcharges.linear import prepare_linear_system, explicit_solution

A, V, Q, resp_charges = prepare_linear_system("data/raw/resp.out", "data/raw/esp.xyz", 78, frame_index=-1)
solver = explicit_solution(ridge=0.0)
result = solver.fit(A, V, Q)
print(result["rmse"], result["sum_q"])
```

The returned dictionary includes the fitted charges `q`, intermediate matrices, RMSE/RRMS, and the enforced total charge.

