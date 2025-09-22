import re
import numpy as np
from pathlib import Path
from typing import Dict, Tuple

BOHR_PER_ANG = 1.8897261254578281
DEBYE_PER_E_BOHR = 2.541746

# --- helpers ---
def _dipole_from_charges(q: np.ndarray, R_bohr: np.ndarray, origin_bohr: np.ndarray) -> Tuple[np.ndarray, float]:
    """μ = Σ q_A (R_A - origin). Returns (vec_D, mag_D)."""
    Rc = R_bohr - origin_bohr[None, :]
    mu_vec_e_bohr = (q[:, None] * Rc).sum(axis=0)
    mu_vec_D = DEBYE_PER_E_BOHR * mu_vec_e_bohr
    return mu_vec_D, float(np.linalg.norm(mu_vec_D))

def _parse_last_resp_bits(resp_out_path: str | Path):
    """
    Parse from resp.out (LAST frame only):
      - ESP unrestrained charges block: (symbols, coords_bohr, charges)
      - CENTER OF MASS line (Å -> bohr)
      - DIPOLE MOMENT line (Debye)
    """
    text = Path(resp_out_path).read_text(errors="ignore")
    lines = text.splitlines()

    # Last ESP unrestrained block
    blocks, cap, cur = [], False, []
    for line in lines:
        if re.search(r"ESP\s+unrestrained\s+charges", line, flags=re.I):
            if cur: blocks.append(cur); cur = []
            cap = True; continue
        if cap:
            if re.search(r"ESP\s+restrained\s+charges", line, flags=re.I):
                if cur: blocks.append(cur); cur = []
                cap = False; continue
            if re.match(r"\s*-{3,}\s*$", line) or re.search(r"^\s*Atom\s+X\s+Y\s+Z\s+Charge\s+Exposure", line, flags=re.I):
                continue
            s = line.strip()
            if not s:
                if cur: blocks.append(cur); cur = []
                cap = False; continue
            parts = s.split()
            if len(parts) >= 6:
                sym = parts[0]
                x, y, z = map(float, parts[1:4])  # bohr (in resp.out table)
                q = float(parts[4])               # column 5 = charge
                cur.append((sym, x, y, z, q))
    if cur: blocks.append(cur)
    if not blocks:
        raise RuntimeError("No 'ESP unrestrained charges' blocks found in resp.out")
    last = blocks[-1]
    symbols_last = [t[0] for t in last]
    coords_last_bohr = np.array([[t[1], t[2], t[3]] for t in last], float)
    charges_last = np.array([t[4] for t in last], float)

    # Last COM (Å -> bohr)
    com_match = None
    for m in re.finditer(r"CENTER OF MASS:\s*\{([^}]+)\}\s*ANGS", text, flags=re.I):
        com_match = m
    if not com_match:
        raise RuntimeError("CENTER OF MASS not found in resp.out")
    com_vals_ang = list(map(float, com_match.group(1).replace(",", " ").split()))
    if len(com_vals_ang) != 3:
        raise RuntimeError("CENTER OF MASS line did not parse 3 numbers")
    com_bohr = np.array(com_vals_ang, float) * BOHR_PER_ANG

    # Last QM dipole (Debye)
    dip_match = None
    dip_pattern = re.compile(
        r"DIPOLE MOMENT:\s*\{([^}]+)\}\s*\(\|D\|\s*=\s*([0-9.Ee+-]+)\)\s*DEBYE",
        re.I,
    )
    for m in dip_pattern.finditer(text):
        dip_match = m
    if not dip_match:
        raise RuntimeError("DIPOLE MOMENT not found in resp.out")
    qm_vec_vals = list(map(float, dip_match.group(1).replace(",", " ").split()))
    qm_dipole_vec_D = np.array(qm_vec_vals, float)
    qm_dipole_mag_D = float(dip_match.group(2))

    return {
        "symbols_last": symbols_last,
        "coords_last_bohr": coords_last_bohr,
        "charges_esp_unrestrained": charges_last,
        "com_bohr": com_bohr,
        "qm_dipole_vec_D": qm_dipole_vec_D,
        "qm_dipole_mag_D": qm_dipole_mag_D,
    }

def three_dipoles_last_frame(resp_out_path: str | Path,
                             R_bohr_lastframe: np.ndarray,
                             q_opt: np.ndarray) -> Dict[str, np.ndarray | float]:
    """
    Compute and return the three dipoles for the LAST frame:

    Inputs:
      - resp_out_path: path to resp.out (must contain the COM & dipole lines shown)
      - R_bohr_lastframe: (N,3) coordinates in BOHR **from the LAST frame in resp.out**
      - q_opt: (N,) your optimized charges (unitless)

    Returns dict with vectors & magnitudes (Debye) for:
      - 'qm_dipole_vec_D', 'qm_dipole_mag_D'
      - 'esp_unrestrained_dipole_vec_D', 'esp_unrestrained_dipole_mag_D'
      - 'optimized_dipole_vec_D', 'optimized_dipole_mag_D'
      plus deltas vs QM.
    """
    parsed = _parse_last_resp_bits(resp_out_path)

    # 1) QM dipole (already Debye)
    qm_vec = parsed["qm_dipole_vec_D"]
    qm_mag = parsed["qm_dipole_mag_D"]

    # 2) Dipole from Terachem's ESP-unrestrained charges (use coords and COM from that last block)
    q_tc = parsed["charges_esp_unrestrained"]
    R_tc = parsed["coords_last_bohr"]
    COM = parsed["com_bohr"]
    esp_vec, esp_mag = _dipole_from_charges(q_tc, R_tc, COM)

    # 3) Dipole from YOUR optimized charges using YOUR last-frame coordinates
    opt_vec, opt_mag = _dipole_from_charges(q_opt, R_bohr_lastframe, COM)

    return {
        # QM
        "qm_dipole_vec_D": qm_vec,
        "qm_dipole_mag_D": qm_mag,
        # ESP (Terachem unrestrained)
        "esp_unrestrained_dipole_vec_D": esp_vec,
        "esp_unrestrained_dipole_mag_D": esp_mag,
        # Yours
        "optimized_dipole_vec_D": opt_vec,
        "optimized_dipole_mag_D": opt_mag,
        # Differences vs QM
        "delta_esp_vs_qm_vec_D": esp_vec - qm_vec,
        "delta_esp_vs_qm_mag_D": esp_mag - qm_mag,
        "delta_opt_vs_qm_vec_D": opt_vec - qm_vec,
        "delta_opt_vs_qm_mag_D": opt_mag - qm_mag,
        # For reference
        "COM_bohr": COM,
    }
