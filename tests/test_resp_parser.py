from __future__ import annotations

from pathlib import Path

import pytest

from resp import ParseRespDotOut


DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "raw"


def test_extract_frames_with_esp():
    parser = ParseRespDotOut(DATA_DIR / "resp.out", 78)
    assert parser.success_check()

    esp_frames, resp_frames = parser.extract_frames_with_esp()

    assert esp_frames, "Expected at least one ESP frame"
    assert resp_frames, "Expected at least one RESP frame"

    first_esp = esp_frames[0]
    first_resp = resp_frames[0]

    assert len(first_esp.positions) == 78
    assert len(first_esp.esp_charges) == 78
    assert len(first_esp.exposure_fractions) == 78
    assert len(first_resp.charges) == 78

    assert first_esp.positions[-1] == pytest.approx((29.349334, 44.452024, 48.860755))
    assert first_esp.esp_charges[-1] == pytest.approx(0.439305, abs=1e-6)
    assert first_resp.charges[-1] == pytest.approx(0.44013, abs=1e-6)

    assert sum(first_esp.esp_charges) == pytest.approx(1.0, abs=1e-5)
    assert sum(first_resp.charges) == pytest.approx(1.0, abs=1e-5)


def test_extract_esp_from_espxyz():
    parser = ParseRespDotOut(DATA_DIR / "resp.out", 78)
    frames = parser.extract_esp_from_espxyz(DATA_DIR / "esp.xyz")

    assert frames, "Expected ESP grid frames from esp.xyz"

    first_frame = frames[0]
    assert len(first_frame.coordinates) == len(first_frame.potentials)
    assert len(first_frame.coordinates) == 11486
    assert first_frame.coordinates[0] == pytest.approx((21.7902980, 23.20555, 23.4869096))
    assert first_frame.potentials[0] == pytest.approx(0.164433438524, abs=1e-12)
