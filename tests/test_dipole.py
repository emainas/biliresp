from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dipole import center_of_mass_bohr_from_xyz, three_dipoles_for_frame
from dipole.dipole import BOHR_PER_ANG
from resp import ParseDotXYZ
from linearESPcharges.linear import explicit_solution, prepare_linear_system

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"
RESP_OUT = DATA_DIR / "resp.out"
ESP_XYZ = DATA_DIR / "esp.xyz"
GEOM_XYZ = DATA_DIR / "1.pose.xyz"
NUMBER_OF_ATOMS = 78


def test_three_dipoles_for_frame():
    A, V, Q, resp_charges, coords_bohr = prepare_linear_system(
        RESP_OUT,
        ESP_XYZ,
        NUMBER_OF_ATOMS,
        frame_index=-1,
        return_positions=True,
    )

    solver = explicit_solution(ridge=0.0)
    res = solver.fit(A, V, Q)

    frame_index = -1
    dipoles = three_dipoles_for_frame(
        RESP_OUT,
        GEOM_XYZ,
        coords_bohr,
        res["q"],
        frame_index=frame_index,
    )

    # Basic sanity: magnitudes and vectors exist
    for key in (
        "qm_dipole_vec_D",
        "qm_dipole_mag_D",
        "terachem_dipole_vec_D",
        "terachem_dipole_mag_D",
        "lagrange_dipole_vec_D",
        "lagrange_dipole_mag_D",
    ):
        assert key in dipoles

    # Optimized dipole should be extremely close to ESP-unrestrained dipole
    np.testing.assert_allclose(
        dipoles["lagrange_dipole_vec_D"],
        dipoles["terachem_dipole_vec_D"],
        atol=5e-4,
    )
    assert dipoles["lagrange_dipole_mag_D"] == pytest.approx(
        dipoles["terachem_dipole_mag_D"], abs=5e-4
    )

    # Both classical dipoles should be within ~0.1 Debye of the QM dipole
    assert dipoles["terachem_dipole_mag_D"] == pytest.approx(
        dipoles["qm_dipole_mag_D"], abs=1e-1
    )
    assert dipoles["lagrange_dipole_mag_D"] == pytest.approx(
        dipoles["qm_dipole_mag_D"], abs=1e-1
    )

    # Ensure COM estimates are close (RESP log vs mass-weighted)
    com_resp = dipoles["COM_bohr_resp"]
    com_mass = dipoles["COM_bohr_mass"]
    np.testing.assert_allclose(com_resp, com_mass, atol=1e-2)

    print("Mass-weighted COM (angstrom):", com_mass / BOHR_PER_ANG)

    print("QM dipole vector (Debye):", dipoles["qm_dipole_vec_D"])
    print("QM |μ| (Debye): {:.6f}".format(dipoles["qm_dipole_mag_D"]))
    print("Terachem dipole vector (Debye):", dipoles["terachem_dipole_vec_D"])
    print("Terachem |μ| (Debye): {:.6f}".format(dipoles["terachem_dipole_mag_D"]))
    print("Lagrange dipole vector (Debye):", dipoles["lagrange_dipole_vec_D"])
    print("Lagrange |μ| (Debye): {:.6f}".format(dipoles["lagrange_dipole_mag_D"]))


def test_center_of_mass_bohr_from_xyz_matches_length():
    # Default path uses coordinates from the xyz file
    com_bohr_default = center_of_mass_bohr_from_xyz(GEOM_XYZ, frame_index=-1)
    assert com_bohr_default.shape == (3,)

    # Custom coordinates in Angstrom should yield the same center
    geom_frames = ParseDotXYZ(GEOM_XYZ).elements()
    coords_ang = np.asarray(geom_frames[0].coordinates, dtype=float)
    com_bohr_custom = center_of_mass_bohr_from_xyz(
        GEOM_XYZ,
        frame_index=-1,
        coords=coords_ang,
        coords_unit="ang",
    )
    np.testing.assert_allclose(com_bohr_default, com_bohr_custom, atol=1e-12)
