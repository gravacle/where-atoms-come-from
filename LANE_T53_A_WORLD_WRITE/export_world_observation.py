#!/usr/bin/env python3
"""Normalize the five VSM pairs into the closed world-observation input contract."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
MODEL = HERE.parent / "model"
if str(MODEL) not in sys.path:
    sys.path.insert(0, str(MODEL))

from lakeshore_vsm import WORLD_OBSERVATION_COLUMNS, normalized_world_observation_csv


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    raw = HERE / "raw"
    source = json.loads((raw / "SOURCE.json").read_text(encoding="utf-8"))
    frozen = {row["file"]: row["sha256"] for row in source["selected_files"]}
    files = sorted(raw.glob("*.csv"))
    if {path.name for path in files} != set(frozen):
        raise SystemExit("REFUSED: source membership differs from SOURCE.json")
    for path in files:
        if digest(path) != frozen[path.name]:
            raise SystemExit(f"REFUSED: source hash mismatch for {path.name}")

    grouped: dict[str, dict[str, Path]] = {}
    for path in files:
        if "-hys-dcd-forc1 " not in path.name:
            raise SystemExit(f"REFUSED: unrecognized source filename {path.name}")
        sample, kind = path.name.split("-hys-dcd-forc1 ", 1)
        grouped.setdefault(sample, {})[kind] = path

    data_path = HERE / "world_observation.csv"
    data_path.write_bytes(normalized_world_observation_csv(files))
    columns = list(WORLD_OBSERVATION_COLUMNS)

    protocol = HERE / "OBSERVATION_PROTOCOL.md"
    source_artifact_paths = [HERE / "raw" / "SOURCE.json", Path(__file__).resolve()] + files
    manifest = {
        "schema": "WAC_WORLD_OBSERVATION_V001",
        "surface_id": "IODP-U1537-MAGNETIC-REMANENCE",
        "sample_ids": sorted(grouped),
        "run_id": "U1537-VSM-HYS-DCD-RETROSPECTIVE-V001",
        "instrument_id": "LAKE-SHORE-8600-VSM-EM7-CSB",
        "source_uri": "https://doi.org/10.5281/zenodo.14564186",
        "evidence_class": "ACTUAL_SURFACE_MEASUREMENT",
        "observation_scope": "WRITE_POST_ONLY",
        "protocol_timing": "RETROSPECTIVE_POST_OUTCOME",
        "protocol_frozen_at_utc": "2026-08-22T08:16:07Z",
        "protocol_path": protocol.relative_to(HERE).as_posix(),
        "protocol_sha256": digest(protocol),
        "data_file": data_path.relative_to(HERE).as_posix(),
        "data_sha256": digest(data_path),
        "columns": columns,
        "units": {"time": "ordinal", "coordinate": "T", "value": "A m^2"},
        "extra_columns": [
            {"name": "writer_field_T", "type": "float", "unit": "T"},
            {"name": "source_time_s", "type": "float", "unit": "s"},
            {"name": "protocol_segment", "type": "string", "unit": "1"},
        ],
        "controls": ["OPPOSITE_WRITE", "LOW_PERTURBATION"],
        "normalization": {
            "adapter_id": "LAKESHORE_VSM_WRITER_OFF_V001",
            "adapter_sha256": digest(MODEL / "lakeshore_vsm.py"),
            "source_paths": [path.relative_to(HERE).as_posix() for path in files],
        },
        "source_artifacts": [
            {
                "path": path.relative_to(HERE).as_posix(),
                "sha256": digest(path),
                "media_type": (
                    "text/csv" if path.suffix == ".csv"
                    else "application/json" if path.suffix == ".json"
                    else "text/x-python"
                ),
            }
            for path in sorted(source_artifact_paths)
        ],
    }
    manifest_path = HERE / "world_observation.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"manifest": str(manifest_path), "rows": 760}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
