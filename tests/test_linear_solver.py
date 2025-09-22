from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from linearESPcharges.linear import (
    KKTblock_solution,
    explicit_solution,
    prepare_linear_system,
)


DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"
RESP_OUT = DATA_DIR / "resp.out"
ESP_XYZ = DATA_DIR / "esp.xyz"
NUMBER_OF_ATOMS = 78


def test_linear_solvers_agree_and_reduce_residual():
    A, V, Q, resp_charges = prepare_linear_system(
        RESP_OUT,
        ESP_XYZ,
        NUMBER_OF_ATOMS,
        frame_index=-1,
    )

    expl = explicit_solution(ridge=0.0)
    kkt = KKTblock_solution(ridge=0.0)

    res_expl = expl.fit(A, V, Q)
    res_kkt = kkt.fit(A, V, Q)

    # sanity: last unrestrained RESP charge first entry should match known value
    assert resp_charges[0] == pytest.approx(-0.517177, abs=1e-6)

    # Solvers should agree to numerical tolerance
    np.testing.assert_allclose(res_expl["q"], res_kkt["q"], rtol=0.0, atol=1e-10)
    np.testing.assert_allclose(res_expl["sum_q"], Q, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(res_kkt["sum_q"], Q, rtol=0.0, atol=1e-12)

    baseline_residual = A @ resp_charges - V
    baseline_rmse = float(np.sqrt(np.mean(baseline_residual**2)))
    baseline_rrms = float(np.sqrt(np.mean(baseline_residual**2) / np.mean(V**2)))
    assert res_expl["rmse"] <= baseline_rmse
    assert res_kkt["rmse"] <= baseline_rmse

    header = (
        "atom  resp_charge   explicit      kkt     explicit-resp     kkt-resp"
    )
    print(header)
    for idx in range(78):
        print(
            f"{idx:4d}  {resp_charges[idx]:+11.6f}  {res_expl['q'][idx]:+11.6f}"
            f"  {res_kkt['q'][idx]:+11.6f}  {res_expl['q'][idx]-resp_charges[idx]:+14.6f}"
            f"  {res_kkt['q'][idx]-resp_charges[idx]:+11.6f}"
        )

    print(
        f"\nΣq (explicit) = {res_expl['sum_q']:.12f}, RMSE = {res_expl['rmse']:.6e}, RRMS = {res_expl['rrms']:.6e}"
    )
    print(
        f"Σq (KKT)      = {res_kkt['sum_q']:.12f}, RMSE = {res_kkt['rmse']:.6e}, RRMS = {res_kkt['rrms']:.6e}"
    )
    reported_rrms = _read_unrestrained_rrms(RESP_OUT)
    print(
        f"Resp baseline Σq = {Q:.12f}, RMSE (calc) = {baseline_rmse:.6e}, RRMS (calc) = {baseline_rrms:.6e}, RRMS (reported) = {reported_rrms:.6e}"
    )


def _read_unrestrained_rrms(resp_path: Path) -> float:
    lines = resp_path.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if 'ESP unrestrained charges:' in line]
    if not starts:
        raise ValueError('ESP unrestrained charges block not found in resp.out')
    start = starts[-1]
    reported = None
    for line in lines[start + 1 :]:
        if 'ESP restrained charges:' in line or 'ESP unrestrained charges:' in line:
            break
        if 'Quality of fit (RRMS):' in line:
            reported = float(line.split(':')[-1])
            break
    if reported is None:
        raise ValueError('Could not find Quality of fit (RRMS) for unrestrained block in resp.out')
    return reported
