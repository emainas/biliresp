from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dipole import three_dipoles_last_frame
from linearESPcharges.linear import explicit_solution, prepare_linear_system

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"
RESP_OUT = DATA_DIR / "resp.out"
ESP_XYZ = DATA_DIR / "esp.xyz"
NUMBER_OF_ATOMS = 78


def test_three_dipoles_last_frame():
    A, V, Q, resp_charges, coords_bohr = prepare_linear_system(
        RESP_OUT,
        ESP_XYZ,
        NUMBER_OF_ATOMS,
        frame_index=-1,
        return_positions=True,
    )

    solver = explicit_solution(ridge=0.0)
    res = solver.fit(A, V, Q)

    dipoles = three_dipoles_last_frame(RESP_OUT, coords_bohr, res["q"])

    # Basic sanity: magnitudes and vectors exist
    for key in (
        "qm_dipole_vec_D",
        "qm_dipole_mag_D",
        "esp_unrestrained_dipole_vec_D",
        "esp_unrestrained_dipole_mag_D",
        "optimized_dipole_vec_D",
        "optimized_dipole_mag_D",
    ):
        assert key in dipoles

    # Optimized dipole should be extremely close to ESP-unrestrained dipole
    np.testing.assert_allclose(
        dipoles["optimized_dipole_vec_D"],
        dipoles["esp_unrestrained_dipole_vec_D"],
        atol=5e-4,
    )
    assert dipoles["optimized_dipole_mag_D"] == pytest.approx(
        dipoles["esp_unrestrained_dipole_mag_D"], abs=5e-4
    )

    # Both classical dipoles should be within ~0.1 Debye of the QM dipole
    assert dipoles["esp_unrestrained_dipole_mag_D"] == pytest.approx(
        dipoles["qm_dipole_mag_D"], abs=1e-1
    )
    assert dipoles["optimized_dipole_mag_D"] == pytest.approx(
        dipoles["qm_dipole_mag_D"], abs=1e-1
    )

    #print("QM dipole vector (Debye):", dipoles["qm_dipole_vec_D"])
    print("QM |μ| (Debye): {:.6f}".format(dipoles["qm_dipole_mag_D"]))
    #print("ESP-unrestrained vector (Debye):", dipoles["esp_unrestrained_dipole_vec_D"])
    print("ESP |μ| (Debye): {:.6f}".format(dipoles["esp_unrestrained_dipole_mag_D"]))
    #print("Δ vector (ESP-QM):", dipoles["delta_esp_vs_qm_vec_D"])
    #print("Δ|μ| (ESP-QM): {:.6f}".format(dipoles["delta_esp_vs_qm_mag_D"]))
    #print("Optimized vector (Debye):", dipoles["optimized_dipole_vec_D"])
    print("Optimized |μ| (Debye): {:.6f}".format(dipoles["optimized_dipole_mag_D"]))
    #print("Δ vector (Opt-QM):", dipoles["delta_opt_vs_qm_vec_D"])
    #print("Δ|μ| (Opt-QM): {:.6f}".format(dipoles["delta_opt_vs_qm_mag_D"]))
