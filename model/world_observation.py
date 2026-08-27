"""Closed, measurement-only input contract for physical record-surface observations.

This module validates custody, units, provenance, and data coverage.  It deliberately does
not decide whether a record formed or whether a general theory is true.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any


SCHEMA = "WAC_WORLD_OBSERVATION_V001"
CERTIFICATE_SCHEMA = "WAC_WORLD_OBSERVATION_CERTIFICATE_V001"
CORE_COLUMNS = (
    "row_id",
    "record_id",
    "event_id",
    "stage",
    "role",
    "time",
    "coordinate",
    "value",
)
MANIFEST_KEYS = {
    "schema",
    "surface_id",
    "sample_ids",
    "run_id",
    "instrument_id",
    "source_uri",
    "evidence_class",
    "observation_scope",
    "protocol_timing",
    "protocol_frozen_at_utc",
    "protocol_path",
    "protocol_sha256",
    "data_file",
    "data_sha256",
    "columns",
    "units",
    "extra_columns",
    "controls",
    "normalization",
    "source_artifacts",
}
STAGES = {"BEFORE", "WRITE", "HOLD", "POST_WRITE", "READ", "RATE"}
ROLES = {"WRITER_CONTROL", "RECORD_READOUT", "ENVIRONMENT", "CONTROL_READOUT"}
EVIDENCE_CLASSES = {"ACTUAL_SURFACE_MEASUREMENT", "SYNTHETIC_TEST_ONLY"}
OBSERVATION_SCOPES = {
    "WRITE_POST_ONLY",
    "FULL_WRITE_HOLD_READ",
    "RATE_STABILITY_ONLY",
    "CONFIGURATION_ONLY",
}
PROTOCOL_TIMINGS = {"PROSPECTIVE_PRE_OUTCOME", "RETROSPECTIVE_POST_OUTCOME"}
CONTROLS = {
    "NO_WRITE",
    "LABEL_SHUFFLE",
    "KNOWN_STATE_COPY",
    "SHORT_WRITE",
    "ERASE",
    "LOW_PERTURBATION",
    "OPPOSITE_WRITE",
    "NONE",
}
EXTRA_TYPES = {"string", "float", "int", "bool"}
SHA256 = re.compile(r"^[0-9a-f]{64}$")
UTC = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
PROTOCOL_MARKER = "# WAC WORLD OBSERVATION PROTOCOL V001"


class ObservationRefusal(ValueError):
    """The supplied observation violates the closed input contract."""


@dataclass(frozen=True)
class WorldObservation:
    manifest_path: Path
    manifest_sha256: str
    manifest: dict[str, Any]
    rows: tuple[dict[str, Any], ...]
    source_hashes: tuple[tuple[str, str], ...]

    def certificate(self) -> dict[str, Any]:
        times = [row["time"] for row in self.rows]
        coordinates = [row["coordinate"] for row in self.rows]
        values = [row["value"] for row in self.rows]
        records = sorted({row["record_id"] for row in self.rows})
        events = {(row["record_id"], row["event_id"]) for row in self.rows}
        stages = sorted({row["stage"] for row in self.rows})
        roles = sorted({row["role"] for row in self.rows})
        classification = {
            "WRITE_POST_ONLY": "CONFIGURATION_EVIDENCE_ONLY",
            "FULL_WRITE_HOLD_READ": "FULL_PROCESS_COVERAGE_ONLY",
            "RATE_STABILITY_ONLY": "RATE_STABILITY_EVIDENCE_ONLY",
            "CONFIGURATION_ONLY": "CONFIGURATION_EVIDENCE_ONLY",
        }[self.manifest["observation_scope"]]
        return {
            "schema": CERTIFICATE_SCHEMA,
            "input_contract_sha256": sha256_file(Path(__file__).resolve()),
            "manifest_sha256": self.manifest_sha256,
            "data_sha256": self.manifest["data_sha256"],
            "surface_id": self.manifest["surface_id"],
            "sample_ids": self.manifest["sample_ids"],
            "run_id": self.manifest["run_id"],
            "instrument_id": self.manifest["instrument_id"],
            "source_uri": self.manifest["source_uri"],
            "evidence_class": self.manifest["evidence_class"],
            "custody_validation": "CONTENT_HASHED_NOT_PHYSICAL_ORIGIN_AUTHENTICATED",
            "independent_reproduction_attested": False,
            "observation_scope": self.manifest["observation_scope"],
            "protocol_timing": self.manifest["protocol_timing"],
            "protocol_sha256": self.manifest["protocol_sha256"],
            "normalization": {
                "adapter_id": self.manifest["normalization"]["adapter_id"],
                "adapter_sha256": self.manifest["normalization"]["adapter_sha256"],
                "normalized_source_join": "VERIFIED",
            },
            "source_hashes": [
                {"path": path, "sha256": digest}
                for path, digest in self.source_hashes
            ],
            "coverage": {
                "row_count": len(self.rows),
                "record_ids": records,
                "event_count": len(events),
                "stages": stages,
                "roles": roles,
                "controls": self.manifest["controls"],
                "time_range": [min(times), max(times)],
                "coordinate_range": [min(coordinates), max(coordinates)],
                "value_range": [min(values), max(values)],
            },
            "scope_classification": classification,
            "scientific_verdict": "NONE_NOT_SCORED",
            "record_formation_proof_authorized": False,
            "universal_claim_authorized": False,
            "public_urm_registration_authorized": False,
        }


def _refuse(message: str) -> None:
    raise ObservationRefusal("WORLD OBSERVATION REFUSES: " + message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            _refuse(f"duplicate JSON member name {key!r}")
        value[key] = item
    return value


def _file_identity(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino


def _normalization_adapter(adapter_id: str):
    if adapter_id == "DIRECT_CSV_COPY_TEST_ONLY_V001":
        def direct(paths: list[Path]) -> bytes:
            if len(paths) != 1:
                _refuse("DIRECT_CSV_COPY_TEST_ONLY_V001 requires exactly one source")
            return paths[0].read_bytes()
        return Path(__file__).resolve(), direct, True
    if adapter_id == "LAKESHORE_VSM_WRITER_OFF_V001":
        from lakeshore_vsm import normalized_world_observation_csv
        import lakeshore_vsm
        return Path(lakeshore_vsm.__file__).resolve(), normalized_world_observation_csv, False
    _refuse(f"normalization adapter {adapter_id!r} is not registered")


def _nonempty(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        _refuse(f"{field} must be a nonempty string")
    return value


def _relative_file(root: Path, value: Any, field: str) -> Path:
    raw = _nonempty(value, field)
    candidate = Path(raw)
    if candidate.is_absolute() or ".." in candidate.parts:
        _refuse(f"{field} must be a safe relative path")
    path = root / candidate
    if not path.is_file():
        _refuse(f"{field} does not name a regular file")
    try:
        resolved = path.resolve(strict=True)
        resolved.relative_to(root.resolve(strict=True))
    except (OSError, ValueError):
        _refuse(f"{field} escapes the observation bundle through a symlink")
    return resolved


def _digest(value: Any, field: str) -> str:
    if not isinstance(value, str) or not SHA256.fullmatch(value):
        _refuse(f"{field} must be a lowercase SHA-256")
    return value


def _parse_extra(value: str, kind: str, field: str) -> Any:
    try:
        if kind == "string":
            if not value:
                _refuse(f"{field} is empty")
            return value
        if kind == "float":
            parsed = float(value)
            if not math.isfinite(parsed):
                _refuse(f"{field} is non-finite")
            return parsed
        if kind == "int":
            return int(value)
        if kind == "bool":
            if value not in {"true", "false"}:
                _refuse(f"{field} must be true or false")
            return value == "true"
    except ValueError as exc:
        _refuse(f"{field} does not match declared type {kind}: {exc}")
    _refuse(f"unsupported extra-column type {kind}")


def load_world_observation(manifest_path: str | Path) -> WorldObservation:
    """Load and validate one closed observation bundle without scientific scoring."""
    manifest_path = Path(manifest_path)
    if manifest_path.is_symlink():
        _refuse("submitted manifest must not be a symlink")
    if manifest_path.parent.is_symlink():
        _refuse("submitted manifest bundle directory must not be a symlink")
    try:
        manifest_bytes = manifest_path.read_bytes()
        manifest = json.loads(manifest_bytes, object_pairs_hook=_unique_json_object)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        _refuse(f"manifest is not readable JSON: {exc}")
    if not isinstance(manifest, dict):
        _refuse("manifest root must be an object")
    if set(manifest) != MANIFEST_KEYS:
        missing = sorted(MANIFEST_KEYS - set(manifest))
        extra = sorted(set(manifest) - MANIFEST_KEYS)
        _refuse(f"manifest key closure failed; missing={missing}, extra={extra}")
    if manifest["schema"] != SCHEMA:
        _refuse(f"unsupported schema {manifest['schema']!r}")

    for field in ("surface_id", "run_id", "instrument_id", "source_uri"):
        _nonempty(manifest[field], field)
    if not manifest["source_uri"].startswith(("https://", "doi:", "urn:")):
        _refuse("source_uri must use https, doi, or urn")
    samples = manifest["sample_ids"]
    if (not isinstance(samples, list) or not samples
            or any(not isinstance(item, str) or not item.strip() for item in samples)
            or len(samples) != len(set(samples))):
        _refuse("sample_ids must be a nonempty unique string list")
    if manifest["evidence_class"] not in EVIDENCE_CLASSES:
        _refuse("evidence_class is not registered")
    if manifest["observation_scope"] not in OBSERVATION_SCOPES:
        _refuse("observation_scope is not registered")
    if manifest["protocol_timing"] not in PROTOCOL_TIMINGS:
        _refuse("protocol_timing is not registered")
    if not isinstance(manifest["protocol_frozen_at_utc"], str) or not UTC.fullmatch(
        manifest["protocol_frozen_at_utc"]
    ):
        _refuse("protocol_frozen_at_utc must be second-resolution UTC")
    try:
        datetime.strptime(manifest["protocol_frozen_at_utc"], "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        _refuse("protocol_frozen_at_utc is not a real UTC date")

    try:
        resolved_manifest = manifest_path.resolve(strict=True)
        root = manifest_path.absolute().parent.resolve(strict=True)
        resolved_manifest.relative_to(root)
    except (OSError, ValueError):
        _refuse("submitted manifest escapes its lexical bundle root")
    protocol_path = _relative_file(root, manifest["protocol_path"], "protocol_path")
    protocol_digest = _digest(manifest["protocol_sha256"], "protocol_sha256")
    if sha256_file(protocol_path) != protocol_digest:
        _refuse("protocol_sha256 mismatch")
    try:
        protocol_first = protocol_path.read_text(encoding="utf-8").splitlines()[0]
    except (UnicodeDecodeError, IndexError):
        _refuse("protocol is not a marked UTF-8 protocol document")
    if protocol_first != PROTOCOL_MARKER:
        _refuse("protocol document marker is not registered")
    data_path = _relative_file(root, manifest["data_file"], "data_file")
    data_digest = _digest(manifest["data_sha256"], "data_sha256")
    if sha256_file(data_path) != data_digest:
        _refuse("data_sha256 mismatch")
    manifest_identity = _file_identity(resolved_manifest)
    protocol_identity = _file_identity(protocol_path)
    data_identity = _file_identity(data_path)
    if len({manifest_identity, protocol_identity, data_identity}) != 3:
        _refuse("manifest, protocol, and data must be distinct physical files")

    extras = manifest["extra_columns"]
    if not isinstance(extras, list):
        _refuse("extra_columns must be a list")
    extra_names: list[str] = []
    extra_types: dict[str, str] = {}
    for index, item in enumerate(extras):
        if not isinstance(item, dict) or set(item) != {"name", "type", "unit"}:
            _refuse(f"extra_columns[{index}] must have exact name/type/unit keys")
        name = _nonempty(item["name"], f"extra_columns[{index}].name")
        if name in CORE_COLUMNS or name in extra_names:
            _refuse(f"duplicate or reserved extra column {name}")
        if item["type"] not in EXTRA_TYPES:
            _refuse(f"unsupported type for extra column {name}")
        _nonempty(item["unit"], f"extra_columns[{index}].unit")
        extra_names.append(name)
        extra_types[name] = item["type"]
    expected_columns = list(CORE_COLUMNS) + extra_names
    if manifest["columns"] != expected_columns:
        _refuse("columns must equal the core prefix plus declared extras in order")
    units = manifest["units"]
    if not isinstance(units, dict) or set(units) != {"time", "coordinate", "value"}:
        _refuse("units must contain exact time/coordinate/value keys")
    for field, value in units.items():
        _nonempty(value, f"units.{field}")

    controls = manifest["controls"]
    if (not isinstance(controls, list) or not controls
            or len(controls) != len(set(controls))
            or any(item not in CONTROLS for item in controls)):
        _refuse("controls must be a nonempty unique registered list")
    if "NONE" in controls and len(controls) != 1:
        _refuse("NONE cannot be combined with another control")

    artifacts = manifest["source_artifacts"]
    if not isinstance(artifacts, list) or not artifacts:
        _refuse("source_artifacts must be nonempty")
    source_hashes: list[tuple[str, str]] = []
    seen_artifacts: set[str] = set()
    seen_artifact_files: set[Path] = set()
    seen_artifact_identities: set[tuple[int, int]] = set()
    source_path_map: dict[str, Path] = {}
    for index, item in enumerate(artifacts):
        if not isinstance(item, dict) or set(item) != {"path", "sha256", "media_type"}:
            _refuse(f"source_artifacts[{index}] has an open or incomplete schema")
        relative = _nonempty(item["path"], f"source_artifacts[{index}].path")
        if relative in seen_artifacts:
            _refuse(f"duplicate source artifact {relative}")
        seen_artifacts.add(relative)
        path = _relative_file(root, relative, f"source_artifacts[{index}].path")
        if path in seen_artifact_files:
            _refuse(f"duplicate resolved source artifact {relative}")
        seen_artifact_files.add(path)
        identity = _file_identity(path)
        if identity in seen_artifact_identities:
            _refuse(f"duplicate physical source artifact {relative}")
        if identity in {manifest_identity, protocol_identity, data_identity}:
            _refuse(f"source artifact {relative} aliases a manifest role file")
        seen_artifact_identities.add(identity)
        expected = _digest(item["sha256"], f"source_artifacts[{index}].sha256")
        _nonempty(item["media_type"], f"source_artifacts[{index}].media_type")
        if sha256_file(path) != expected:
            _refuse(f"source artifact hash mismatch for {relative}")
        source_hashes.append((relative, expected))
        source_path_map[relative] = path

    normalization = manifest["normalization"]
    if not isinstance(normalization, dict) or set(normalization) != {
        "adapter_id", "adapter_sha256", "source_paths"
    }:
        _refuse("normalization must have exact adapter_id/adapter_sha256/source_paths keys")
    adapter_id = _nonempty(normalization["adapter_id"], "normalization.adapter_id")
    adapter_path, adapter, test_only = _normalization_adapter(adapter_id)
    if test_only and manifest["evidence_class"] != "SYNTHETIC_TEST_ONLY":
        _refuse("test-only normalization cannot certify actual-surface evidence")
    adapter_digest = _digest(normalization["adapter_sha256"], "normalization.adapter_sha256")
    if sha256_file(adapter_path) != adapter_digest:
        _refuse("normalization adapter hash mismatch")
    normalization_sources = normalization["source_paths"]
    if (not isinstance(normalization_sources, list) or not normalization_sources
            or len(normalization_sources) != len(set(normalization_sources))
            or any(not isinstance(item, str) for item in normalization_sources)):
        _refuse("normalization.source_paths must be a nonempty unique string list")
    if any(item not in source_path_map for item in normalization_sources):
        _refuse("normalization source is absent from source_artifacts")
    try:
        rebuilt_data = adapter([source_path_map[item] for item in normalization_sources])
    except ObservationRefusal:
        raise
    except Exception as exc:
        _refuse(f"registered normalization adapter failed: {type(exc).__name__}: {exc}")
    if rebuilt_data != data_path.read_bytes():
        _refuse("normalized CSV is not the registered adapter output for the frozen sources")

    rows: list[dict[str, Any]] = []
    row_ids: set[str] = set()
    last_time_by_event: dict[tuple[str, str], float] = {}
    try:
        with data_path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, strict=True)
            if reader.fieldnames != expected_columns:
                _refuse("CSV header differs from closed columns")
            for line_number, raw in enumerate(reader, start=2):
                if None in raw:
                    _refuse(f"CSV row {line_number} has excess cells")
                if any(value is None for value in raw.values()):
                    _refuse(f"CSV row {line_number} has missing cells")
                row_id = _nonempty(raw["row_id"], f"row {line_number} row_id")
                if row_id in row_ids:
                    _refuse(f"duplicate row_id {row_id}")
                row_ids.add(row_id)
                record_id = _nonempty(raw["record_id"], f"row {line_number} record_id")
                if record_id not in samples:
                    _refuse(f"row {line_number} record_id is outside sample_ids")
                event_id = _nonempty(raw["event_id"], f"row {line_number} event_id")
                stage = raw["stage"]
                role = raw["role"]
                if stage not in STAGES:
                    _refuse(f"row {line_number} stage is not registered")
                if role not in ROLES:
                    _refuse(f"row {line_number} role is not registered")
                numeric: dict[str, float] = {}
                for field in ("time", "coordinate", "value"):
                    try:
                        numeric[field] = float(raw[field])
                    except ValueError:
                        _refuse(f"row {line_number} {field} is nonnumeric")
                    if not math.isfinite(numeric[field]):
                        _refuse(f"row {line_number} {field} is non-finite")
                event_key = (record_id, event_id)
                if numeric["time"] < last_time_by_event.get(event_key, -math.inf):
                    _refuse(f"row {line_number} time decreases within record/event")
                last_time_by_event[event_key] = numeric["time"]
                parsed: dict[str, Any] = {
                    "row_id": row_id,
                    "record_id": record_id,
                    "event_id": event_id,
                    "stage": stage,
                    "role": role,
                    **numeric,
                }
                for name, kind in extra_types.items():
                    parsed[name] = _parse_extra(raw[name], kind, f"row {line_number} {name}")
                rows.append(parsed)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        _refuse(f"CSV is unreadable: {exc}")
    if not rows:
        _refuse("CSV contains no observations")

    observed_records = {row["record_id"] for row in rows}
    if observed_records != set(samples):
        _refuse("every declared sample must have at least one observation row")
    stages_by_record = {
        sample: {row["stage"] for row in rows if row["record_id"] == sample}
        for sample in samples
    }
    scope = manifest["observation_scope"]
    if scope == "WRITE_POST_ONLY":
        if any("POST_WRITE" not in stages for stages in stages_by_record.values()):
            _refuse("WRITE_POST_ONLY requires a POST_WRITE row for every sample")
        if any("BEFORE" in stages for stages in stages_by_record.values()):
            _refuse("WRITE_POST_ONLY excludes BEFORE rows")
    elif scope == "FULL_WRITE_HOLD_READ":
        required = {"BEFORE", "WRITE", "HOLD", "READ"}
        trajectories: dict[tuple[str, str], list[dict[str, Any]]] = {}
        for row in rows:
            trajectories.setdefault((row["record_id"], row["event_id"]), []).append(row)
        complete_samples: set[str] = set()
        for (record_id, _), trajectory in trajectories.items():
            stages = {row["stage"] for row in trajectory}
            if not required.issubset(stages):
                continue
            write_rows = [row for row in trajectory if row["stage"] == "WRITE"]
            read_rows = [row for row in trajectory if row["stage"] == "READ"]
            if not any(row["role"] == "WRITER_CONTROL" for row in write_rows):
                _refuse("complete trajectory lacks a WRITER_CONTROL row at WRITE")
            if not any(
                row["role"] in {"RECORD_READOUT", "CONTROL_READOUT"}
                for row in read_rows
            ):
                _refuse("complete trajectory lacks a readout row at READ")
            bounds = {
                stage: (
                    min(row["time"] for row in trajectory if row["stage"] == stage),
                    max(row["time"] for row in trajectory if row["stage"] == stage),
                )
                for stage in required
            }
            if not (
                bounds["BEFORE"][1] <= bounds["WRITE"][0]
                <= bounds["WRITE"][1] <= bounds["HOLD"][0]
                <= bounds["HOLD"][1] <= bounds["READ"][0]
            ):
                _refuse("complete trajectory stages are not time ordered")
            complete_samples.add(record_id)
        if complete_samples != set(samples):
            _refuse("FULL_WRITE_HOLD_READ lacks a complete joined trajectory for every sample")
    elif scope == "RATE_STABILITY_ONLY" and any(
        "RATE" not in stages for stages in stages_by_record.values()
    ):
        _refuse("RATE_STABILITY_ONLY requires a RATE row for every sample")

    return WorldObservation(
        manifest_path=manifest_path.resolve(),
        manifest_sha256=sha256_bytes(manifest_bytes),
        manifest=manifest,
        rows=tuple(rows),
        source_hashes=tuple(source_hashes),
    )


def certificate_json(observation: WorldObservation) -> str:
    return json.dumps(observation.certificate(), indent=2, sort_keys=True) + "\n"
