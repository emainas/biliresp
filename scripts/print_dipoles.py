from __future__ import annotations

import argparse
from pathlib import Path

from dipole import three_dipoles_for_frame
from linearESPcharges.linear import explicit_solution, prepare_linear_system


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print QM, ESP-unrestrained, and fitted dipole moments for the selected frame."
    )
    parser.add_argument("resp_out", type=Path)
    parser.add_argument("esp_xyz", type=Path)
    parser.add_argument("geom_xyz", type=Path)
    parser.add_argument("n_atoms", type=int)
    parser.add_argument("--frame", type=int, default=-1, help="Frame index (default: last)")
    args = parser.parse_args()

    A, V, Q, resp_charges, coords_bohr = prepare_linear_system(
        args.resp_out,
        args.esp_xyz,
        args.n_atoms,
        frame_index=args.frame,
        return_positions=True,
    )

    solver = explicit_solution()
    fit_result = solver.fit(A, V, Q)

    dipoles = three_dipoles_for_frame(
        args.resp_out,
        args.geom_xyz,
        coords_bohr,
        fit_result["q"],
        frame_index=args.frame,
    )

    print("QM (from resp.out log)")
    print("  vector (Debye):", dipoles["qm_dipole_vec_D"])
    print("  |μ| (Debye):   {:.6f}".format(dipoles["qm_dipole_mag_D"]))
    print()
    print("Terachem charges (RESP log)")
    print("  vector (Debye):", dipoles["terachem_dipole_vec_D"])
    print("  |μ| (Debye):   {:.6f}".format(dipoles["terachem_dipole_mag_D"]))
    print("  Δ vector vs QM (Debye):", dipoles["delta_terachem_vs_qm_vec_D"])
    print("  Δ|μ| vs QM (Debye):     {:.6f}".format(dipoles["delta_terachem_vs_qm_mag_D"]))
    print()
    print("Lagrange multiplier fit (explicit)")
    print("  vector (Debye):", dipoles["lagrange_dipole_vec_D"])
    print("  |μ| (Debye):   {:.6f}".format(dipoles["lagrange_dipole_mag_D"]))
    print("  Δ vector vs QM (Debye):", dipoles["delta_lagrange_vs_qm_vec_D"])
    print("  Δ|μ| vs QM (Debye):     {:.6f}".format(dipoles["delta_lagrange_vs_qm_mag_D"]))


if __name__ == "__main__":
    main()
