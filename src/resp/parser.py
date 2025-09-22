from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple


@dataclass
class ESPFrame:
    positions: List[Tuple[float, float, float]]
    esp_charges: List[float]
    exposure_fractions: List[float]


@dataclass
class RESPFrame:
    charges: List[float]


@dataclass
class ESPGridFrame:
    coordinates: List[Tuple[float, float, float]]
    potentials: List[float]


class ParseRespDotOut:
    def __init__(self, filename: Path | str, number_of_atoms: int) -> None:
        self.file = Path(filename)
        self.number_of_atoms = number_of_atoms

    # made this function here so that we can look through the file make sure the job is finished just like in resp
    def success_check(self) -> bool:
        success_string = "| Job finished:"
        with self.file.open() as f:
            for line in f:
                if line.strip().startswith(success_string):
                    return True
        return False

    def extract_frames_with_esp(self) -> Tuple[List[ESPFrame], List[RESPFrame]]:
        lines = self._read_lines()

        unrestrained_header = "ESP unrestrained charges:"
        restrained_header = "ESP restrained charges:"

        esp_indices = [i for i, line in enumerate(lines) if line.strip() == unrestrained_header]
        resp_indices = [i for i, line in enumerate(lines) if line.strip() == restrained_header]

        esp_frames: List[ESPFrame] = []
        for idx in esp_indices:
            start = idx + 3
            data_lines = lines[start : start + self.number_of_atoms]
            positions: List[Tuple[float, float, float]] = []
            esp_charges: List[float] = []
            exposures: List[float] = []
            for line in data_lines:
                parts = line.split()
                positions.append((float(parts[1]), float(parts[2]), float(parts[3])))
                esp_charges.append(float(parts[4]))
                exposures.append(float(parts[5]))
            esp_frames.append(ESPFrame(positions, esp_charges, exposures))

        resp_frames: List[RESPFrame] = []
        for idx in resp_indices:
            start = idx + 3
            data_lines = lines[start : start + self.number_of_atoms]
            charges = [float(line.split()[4]) for line in data_lines]
            resp_frames.append(RESPFrame(charges))

        return esp_frames, resp_frames

    def extract_esp_from_espxyz(self, filename: Path | str) -> List[ESPGridFrame]:
        frames: List[ESPGridFrame] = []
        with Path(filename).open() as f:
            lines = [line.strip() for line in f if line.strip()]
        idx = 0
        total_lines = len(lines)
        while idx < total_lines:
            natoms = int(lines[idx])
            idx += 1  # skip atom count
            idx += 1  # skip comment line
            coords: List[Tuple[float, float, float]] = []
            potentials: List[float] = []
            for _ in range(natoms):
                parts = lines[idx].split()
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
                potentials.append(float(parts[4]))
                idx += 1
            frames.append(ESPGridFrame(coords, potentials))
        return frames

    def _read_lines(self) -> List[str]:
        with self.file.open() as f:
            return [line.rstrip("\n") for line in f]
