#!/usr/bin/env python3
"""Independent, default-refuted verifier for the T-53A observation ingest path.

This verifier owns no scientific scoring rule.  It checks the frozen implementation,
reconstructs the normalized table without importing ``lakeshore_vsm``, and attacks the
generic contract with deliberately malformed or falsely upgraded bundles.
"""

from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Callable


VERIFY = Path(__file__).resolve().parent
LANE = VERIFY.parent
REPO = LANE.parent
MODEL = REPO / "model"
RAW = LANE / "raw"

FROZEN_HASHES = {
    "model/world_observation.py": "842451df14b6cf5677c99e103753daf5312c02aded93e775e4545c5c9b8afa99",
    "model/run_world_observation.py": "7a10ad8ee91b5283c00f42492eac7dada2df2a434610f2f467abfe9b1c0fb922",
    "model/checks_world_observation.py": "e56ba0eb31c3ad2e0f1b0bd773e85c3d6db45d243f9c53518277c1c372e8267b",
    "LANE_T53_A_WORLD_WRITE/OBSERVATION_PROTOCOL.md": "926be7fe16940d7874a05e1ef1bdd8d0251c3b0c91ae39b1c7f314a0542c5060",
    "LANE_T53_A_WORLD_WRITE/export_world_observation.py": "0f2cd50064c3f9bab0787edf030347ed082f3939ddc8df6416b713aec623b5f0",
    "LANE_T53_A_WORLD_WRITE/world_observation.csv": "82fca51d5f7923107763b01340a85842d5a06d770ecc88e8e518a696a2e3a891",
    "LANE_T53_A_WORLD_WRITE/world_observation.json": "01a9221c5a549401eddc8312a36d2e8cfb5daf09c5115de79565d9fbbdf203f5",
    "LANE_T53_A_WORLD_WRITE/world_observation_certificate.json": "497b3f875e9aaafe5e317c5b4e7d5cde81483eedaae914fa7f093d3a3ae1bdae",
}

CORE = [
    "row_id",
    "record_id",
    "event_id",
    "stage",
    "role",
    "time",
    "coordinate",
    "value",
]
PAIR_TOKEN = "-hys-dcd-forc1 "


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_json(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def freeze_check() -> dict[str, Any]:
    actual = {name: digest(REPO / name) for name in FROZEN_HASHES}
    return {
        "expected": FROZEN_HASHES,
        "actual": actual,
        "pass": actual == FROZEN_HASHES,
    }


def import_target():
    sys.path.insert(0, str(MODEL))
    try:
        import world_observation  # type: ignore
    finally:
        sys.path.pop(0)
    return world_observation


TARGET = import_target()


def run_official_paths() -> dict[str, Any]:
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    check = subprocess.run(
        [sys.executable, "-B", str(MODEL / "checks_world_observation.py")],
        cwd=REPO,
        env=environment,
        capture_output=True,
        check=False,
    )
    command = [
        sys.executable,
        "-B",
        str(MODEL / "run_world_observation.py"),
        str(LANE / "world_observation.json"),
    ]
    renders = [
        subprocess.run(
            command,
            cwd=REPO,
            env=environment,
            capture_output=True,
            check=False,
        )
        for _ in range(3)
    ]
    stored = (LANE / "world_observation_certificate.json").read_bytes()
    outputs = [result.stdout for result in renders]
    parsed = json.loads(outputs[0]) if outputs and renders[0].returncode == 0 else {}
    safety = {
        "scope_classification": parsed.get("scope_classification"),
        "scientific_verdict": parsed.get("scientific_verdict"),
        "record_formation_proof_authorized": parsed.get(
            "record_formation_proof_authorized"
        ),
        "universal_claim_authorized": parsed.get("universal_claim_authorized"),
        "public_urm_registration_authorized": parsed.get(
            "public_urm_registration_authorized"
        ),
        "independent_reproduction_attested": parsed.get(
            "independent_reproduction_attested"
        ),
        "custody_validation": parsed.get("custody_validation"),
        "row_count": parsed.get("coverage", {}).get("row_count"),
    }
    expected_safety = {
        "scope_classification": "CONFIGURATION_EVIDENCE_ONLY",
        "scientific_verdict": "NONE_NOT_SCORED",
        "record_formation_proof_authorized": False,
        "universal_claim_authorized": False,
        "public_urm_registration_authorized": False,
        "independent_reproduction_attested": False,
        "custody_validation": "CONTENT_HASHED_NOT_PHYSICAL_ORIGIN_AUTHENTICATED",
        "row_count": 760,
    }
    return {
        "local_checks_returncode": check.returncode,
        "local_checks_stdout": check.stdout.decode("utf-8", errors="replace").strip(),
        "local_checks_stderr": check.stderr.decode("utf-8", errors="replace").strip(),
        "cli_returncodes": [result.returncode for result in renders],
        "three_runs_byte_identical": len(set(outputs)) == 1,
        "stored_certificate_byte_identical": bool(outputs) and outputs[0] == stored,
        "render_sha256": hashlib.sha256(outputs[0]).hexdigest() if outputs else None,
        "safety_fields": safety,
        "safety_fields_pass": safety == expected_safety,
        "pass": (
            check.returncode == 0
            and check.stdout.strip() == b"WORLD_OBSERVATION_CHECKS: 16/16 PASS"
            and all(result.returncode == 0 for result in renders)
            and len(set(outputs)) == 1
            and outputs[0] == stored
            and safety == expected_safety
        ),
    }


def one_header(lines: list[str], prefix: str) -> str:
    hits = [line.split(":", 1)[1].strip() for line in lines if line.startswith(prefix + ":")]
    if len(hits) != 1:
        raise AssertionError(f"expected exactly one {prefix!r}, got {len(hits)}")
    return hits[0]


def parse_raw_independently(path: Path) -> dict[str, Any]:
    """Parse only the documented ASCII portions; this does not import project adapters."""
    payload = path.read_bytes()
    lines = payload.decode("latin-1").splitlines()
    table_indices = [i for i, line in enumerate(lines) if line.startswith("Step,Iteration,")]
    if len(table_indices) != 1:
        raise AssertionError(f"{path.name}: expected one data table")
    table = table_indices[0]
    fields: list[float] = []
    moments: list[float] = []
    times: list[float] = []
    for line_number, row in enumerate(csv.reader(lines[table + 1 :]), start=table + 2):
        if not row or all(not item.strip() for item in row):
            continue
        if len(row) < 8:
            raise AssertionError(f"{path.name}:{line_number}: short source row")
        values = [float(row[index]) for index in (3, 4, 5)]
        if not all(math.isfinite(item) for item in values):
            raise AssertionError(f"{path.name}:{line_number}: nonfinite source value")
        if row[6].strip() != "GOOD" or row[7].strip() != "GOOD":
            raise AssertionError(f"{path.name}:{line_number}: non-GOOD source status")
        fields.append(values[0])
        moments.append(values[1])
        times.append(values[2])
    if len(fields) < 3 or any(right < left for left, right in zip(times, times[1:])):
        raise AssertionError(f"{path.name}: insufficient or nonmonotone trace")
    start = datetime.strptime(one_header(lines, "START TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    finish = datetime.strptime(one_header(lines, "FINISH TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    if finish < start:
        raise AssertionError(f"{path.name}: finish precedes start")
    if "#HYSTERESIS MEASUREMENT" in lines:
        kind = "hys"
        writer_off = None
    elif "#REMANENCE CURVES MEASUREMENT" in lines:
        kind = "dcd"
        flag = one_header(lines, "Measure moment at applied fields")
        if flag != "False":
            raise AssertionError(f"{path.name}: DCD is not a writer-off readout")
        writer_off = True
    else:
        raise AssertionError(f"{path.name}: unsupported measurement")
    return {
        "path": path,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "kind": kind,
        "fields": fields,
        "moments": moments,
        "times": times,
        "start": start,
        "finish": finish,
        "magnet": one_header(lines, "Magnet"),
        "instrument_sample_id": one_header(lines, "ID"),
        "writer_off": writer_off,
    }


def strict_interpolant(x: list[float], y: list[float]) -> tuple[float, float, int]:
    crossings = [index for index in range(len(x) - 1) if x[index] * x[index + 1] < 0]
    if len(crossings) != 1:
        raise AssertionError(f"expected one strict crossing, got {len(crossings)}")
    index = crossings[0]
    fraction = -x[index] / (x[index + 1] - x[index])
    return y[index] + fraction * (y[index + 1] - y[index]), fraction, index


def expected_normalized_rows(sample: str, hys: dict[str, Any], dcd: dict[str, Any]):
    turn = min(range(len(hys["fields"])), key=hys["fields"].__getitem__)
    branches = [
        (
            "positive-history",
            hys["fields"][: turn + 1],
            hys["moments"][: turn + 1],
            hys["times"][: turn + 1],
            max(hys["fields"]),
        ),
        (
            "negative-history",
            hys["fields"][turn:],
            hys["moments"][turn:],
            hys["times"][turn:],
            min(hys["fields"]),
        ),
    ]
    rows: list[dict[str, Any]] = []
    ordinal = 0
    for label, fields, moments, times, writer_field in branches:
        moment, fraction, crossing = strict_interpolant(fields, moments)
        source_time = times[crossing] + fraction * (times[crossing + 1] - times[crossing])
        rows.append(
            {
                "row_id": f"{sample}:hys:{label}",
                "record_id": sample,
                "event_id": f"hys:{label}",
                "stage": "POST_WRITE",
                "role": "RECORD_READOUT",
                "time": float(ordinal),
                "coordinate": 0.0,
                "value": moment,
                "writer_field_T": writer_field,
                "source_time_s": source_time,
                "protocol_segment": "HYSTERESIS_ZERO_FIELD_INTERPOLANT",
            }
        )
        ordinal += 1
    for index, (field, moment, source_time) in enumerate(
        zip(dcd["fields"], dcd["moments"], dcd["times"])
    ):
        rows.append(
            {
                "row_id": f"{sample}:dcd:{index:03d}",
                "record_id": sample,
                "event_id": f"dcd:{index:03d}",
                "stage": "POST_WRITE",
                "role": "RECORD_READOUT",
                "time": float(ordinal),
                "coordinate": 0.0,
                "value": moment,
                "writer_field_T": field,
                "source_time_s": source_time,
                "protocol_segment": "DCD_WRITER_OFF_READOUT",
            }
        )
        ordinal += 1
    return rows


def close_float(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=2e-15, abs_tol=1e-20)


def verify_actual_join() -> dict[str, Any]:
    source = json.loads((RAW / "SOURCE.json").read_text(encoding="utf-8"))
    selected = {item["file"]: item for item in source["selected_files"]}
    manifest = json.loads((LANE / "world_observation.json").read_text(encoding="utf-8"))
    artifacts = {item["path"]: item for item in manifest["source_artifacts"]}
    raw_files = sorted(RAW.glob("*.csv"))
    source_membership = {path.name for path in raw_files} == set(selected)
    raw_hashes_pass = all(
        digest(path) == selected[path.name]["sha256"]
        and digest(path) == artifacts[f"raw/{path.name}"]["sha256"]
        and path.stat().st_size == selected[path.name]["bytes"]
        for path in raw_files
    )
    traces = {path.name: parse_raw_independently(path) for path in raw_files}
    grouped: dict[str, dict[str, dict[str, Any]]] = {}
    for name, trace in traces.items():
        sample, remainder = name.split(PAIR_TOKEN, 1)
        kind = "hys" if remainder.startswith("Step 1 ") else "dcd"
        grouped.setdefault(sample, {})[kind] = trace
    pair_pass = (
        sorted(grouped) == manifest["sample_ids"]
        and all(set(pair) == {"hys", "dcd"} for pair in grouped.values())
    )
    expected: list[dict[str, Any]] = []
    for sample in sorted(grouped):
        expected.extend(expected_normalized_rows(sample, grouped[sample]["hys"], grouped[sample]["dcd"]))
    with (LANE / "world_observation.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        actual = list(csv.DictReader(stream))
    exact_fields = ["row_id", "record_id", "event_id", "stage", "role", "protocol_segment"]
    numeric_fields = ["time", "coordinate", "value", "writer_field_T", "source_time_s"]
    mismatches: list[str] = []
    if len(actual) != len(expected):
        mismatches.append(f"row count: actual={len(actual)} expected={len(expected)}")
    for index, (actual_row, expected_row) in enumerate(zip(actual, expected)):
        for field in exact_fields:
            if actual_row[field] != expected_row[field]:
                mismatches.append(f"row {index} {field}")
        for field in numeric_fields:
            try:
                agrees = close_float(float(actual_row[field]), float(expected_row[field]))
            except ValueError:
                agrees = False
            if not agrees:
                mismatches.append(f"row {index} {field}")
        if len(mismatches) >= 20:
            break
    frozen = datetime.strptime(manifest["protocol_frozen_at_utc"], "%Y-%m-%dT%H:%M:%SZ")
    all_outcomes_precede_protocol = all(trace["finish"] < frozen for trace in traces.values())
    protocol_text = (LANE / "OBSERVATION_PROTOCOL.md").read_text(encoding="utf-8")
    protocol_words = " ".join(protocol_text.split())
    protocol_honesty = all(
        phrase in protocol_words
        for phrase in (
            "retrospective, measurement-only, no scientific scoring rule",
            "WRITE_POST_ONLY",
            "neither this normalization nor a successful ingestion certificate is a record-formation proof",
        )
    )
    raw_semantics = all(
        trace["magnet"] == "EM7-CSB"
        and trace["instrument_sample_id"] == "Sample1"
        and (trace["kind"] != "dcd" or trace["writer_off"] is True)
        for trace in traces.values()
    )
    return {
        "raw_file_count": len(raw_files),
        "sample_count": len(grouped),
        "expected_normalized_rows": len(expected),
        "actual_normalized_rows": len(actual),
        "source_membership_pass": source_membership,
        "source_and_manifest_hashes_pass": raw_hashes_pass,
        "five_complete_pairs_pass": pair_pass,
        "independent_row_join_pass": not mismatches,
        "first_mismatches": mismatches,
        "all_raw_outcomes_precede_protocol_freeze": all_outcomes_precede_protocol,
        "retrospective_label_honest": (
            manifest["protocol_timing"] == "RETROSPECTIVE_POST_OUTCOME"
            and all_outcomes_precede_protocol
        ),
        "writer_off_and_instrument_headers_pass": raw_semantics,
        "protocol_nonproof_boundary_pass": protocol_honesty,
        "pass": all(
            (
                source_membership,
                raw_hashes_pass,
                pair_pass,
                not mismatches,
                len(actual) == 760,
                all_outcomes_precede_protocol,
                raw_semantics,
                protocol_honesty,
            )
        ),
    }


def write_synthetic_bundle(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    protocol = root / "PROTOCOL.txt"
    source = root / "raw.dat"
    data = root / "observations.csv"
    protocol.write_text("synthetic protocol fixture\n", encoding="utf-8")
    source.write_bytes(b"synthetic raw fixture\n")
    columns = CORE + ["quality"]
    with data.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(["r1", "sample-a", "event-a", "POST_WRITE", "RECORD_READOUT", "0", "0", "2", "1"])
        writer.writerow(["r2", "sample-a", "event-a", "POST_WRITE", "RECORD_READOUT", "1", "0", "-2", "1"])
    manifest = {
        "schema": "WAC_WORLD_OBSERVATION_V001",
        "surface_id": "synthetic-surface",
        "sample_ids": ["sample-a"],
        "run_id": "synthetic-run",
        "instrument_id": "synthetic-instrument",
        "source_uri": "urn:wac:synthetic-fixture",
        "evidence_class": "SYNTHETIC_TEST_ONLY",
        "observation_scope": "WRITE_POST_ONLY",
        "protocol_timing": "PROSPECTIVE_PRE_OUTCOME",
        "protocol_frozen_at_utc": "2026-08-22T00:00:00Z",
        "protocol_path": protocol.name,
        "protocol_sha256": digest(protocol),
        "data_file": data.name,
        "data_sha256": digest(data),
        "columns": columns,
        "units": {"time": "s", "coordinate": "1", "value": "1"},
        "extra_columns": [{"name": "quality", "type": "float", "unit": "1"}],
        "controls": ["OPPOSITE_WRITE"],
        "source_artifacts": [
            {"path": source.name, "sha256": digest(source), "media_type": "application/octet-stream"}
        ],
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(stable_json(manifest), encoding="utf-8")
    return manifest_path


def load_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def save_manifest(path: Path, manifest: dict[str, Any]) -> None:
    path.write_text(stable_json(manifest), encoding="utf-8")


def update_data_hash(manifest_path: Path) -> None:
    manifest = load_manifest(manifest_path)
    manifest["data_sha256"] = digest(manifest_path.parent / manifest["data_file"])
    save_manifest(manifest_path, manifest)


def target_accepts(path: Path) -> tuple[bool, str]:
    try:
        observation = TARGET.load_world_observation(path)
    except TARGET.ObservationRefusal as exc:
        return False, str(exc)
    return True, TARGET.certificate_json(observation)


def edit_manifest(path: Path, mutation: Callable[[dict[str, Any]], None]) -> None:
    value = load_manifest(path)
    mutation(value)
    save_manifest(path, value)


def replace_data(path: Path, rows: list[list[str]], extra: tuple[str, str, str] | None = None) -> None:
    manifest = load_manifest(path)
    if extra is None:
        columns = CORE
        manifest["extra_columns"] = []
    else:
        name, kind, unit = extra
        columns = CORE + [name]
        manifest["extra_columns"] = [{"name": name, "type": kind, "unit": unit}]
    manifest["columns"] = columns
    data = path.parent / manifest["data_file"]
    with data.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    manifest["data_sha256"] = digest(data)
    save_manifest(path, manifest)


def mutation_suite() -> dict[str, Any]:
    scratch = VERIFY / "_scratch"
    scratch.mkdir(exist_ok=True)
    cases: list[dict[str, Any]] = []

    def execute(name: str, expected: str, mutation: Callable[[Path], Path | None]) -> None:
        with tempfile.TemporaryDirectory(prefix=name + "-", dir=scratch) as temporary:
            root = Path(temporary)
            manifest_path = write_synthetic_bundle(root / "bundle")
            returned = mutation(manifest_path)
            target = returned if returned is not None else manifest_path
            accepted, detail = target_accepts(target)
            observed = "ACCEPT" if accepted else "REFUSE"
            cases.append(
                {
                    "name": name,
                    "expected": expected,
                    "observed": observed,
                    "expectation_met": observed == expected,
                    "detail": detail[:500],
                }
            )

    execute("valid_baseline", "ACCEPT", lambda path: None)
    execute("open_outer_manifest", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("extra", 1)), None)[1])
    execute("open_source_schema", "REFUSE", lambda path: (edit_manifest(path, lambda value: value["source_artifacts"][0].__setitem__("extra", 1)), None)[1])
    execute("open_extra_column_schema", "REFUSE", lambda path: (edit_manifest(path, lambda value: value["extra_columns"][0].__setitem__("extra", 1)), None)[1])
    execute("open_units_schema", "REFUSE", lambda path: (edit_manifest(path, lambda value: value["units"].__setitem__("extra", "1")), None)[1])
    execute("parent_path_escape", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("protocol_path", "../PROTOCOL.txt")), None)[1])
    execute("absolute_path_escape", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("protocol_path", str((path.parent / "PROTOCOL.txt").resolve()))), None)[1])

    def outside_symlink(path: Path) -> None:
        outside = path.parent.parent / "outside.dat"
        outside.write_bytes(b"outside\n")
        link = path.parent / "outside-link"
        link.symlink_to(outside)
        edit_manifest(
            path,
            lambda value: value["source_artifacts"].__setitem__(
                0,
                {"path": link.name, "sha256": digest(outside), "media_type": "application/octet-stream"},
            ),
        )

    execute("source_symlink_escape", "REFUSE", outside_symlink)

    def nested_symlink(path: Path) -> None:
        outside = path.parent.parent / "outside-dir"
        outside.mkdir()
        payload = outside / "payload.dat"
        payload.write_bytes(b"outside nested\n")
        (path.parent / "alias-dir").symlink_to(outside, target_is_directory=True)
        edit_manifest(
            path,
            lambda value: value["source_artifacts"].__setitem__(
                0,
                {"path": "alias-dir/payload.dat", "sha256": digest(payload), "media_type": "application/octet-stream"},
            ),
        )

    execute("nested_directory_symlink_escape", "REFUSE", nested_symlink)

    def duplicate_symlink(path: Path) -> None:
        alias = path.parent / "raw-alias.dat"
        alias.symlink_to(path.parent / "raw.dat")
        manifest = load_manifest(path)
        manifest["source_artifacts"].append(
            {"path": alias.name, "sha256": digest(alias), "media_type": "application/octet-stream"}
        )
        save_manifest(path, manifest)

    execute("duplicate_source_symlink_alias", "REFUSE", duplicate_symlink)
    execute("data_hash_drift", "REFUSE", lambda path: ((path.parent / "observations.csv").write_text("changed\n", encoding="utf-8"), None)[1])
    execute("protocol_hash_drift", "REFUSE", lambda path: ((path.parent / "PROTOCOL.txt").write_text("changed\n", encoding="utf-8"), None)[1])
    execute("source_hash_drift", "REFUSE", lambda path: ((path.parent / "raw.dat").write_bytes(b"changed\n"), None)[1])

    def wrong_header(path: Path) -> None:
        data = path.parent / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("quality", "wrong"), encoding="utf-8")
        update_data_hash(path)

    execute("wrong_csv_header", "REFUSE", wrong_header)

    def excess_cell(path: Path) -> None:
        data = path.parent / "observations.csv"
        lines = data.read_text(encoding="utf-8").splitlines()
        lines[1] += ",excess"
        data.write_text("\n".join(lines) + "\n", encoding="utf-8")
        update_data_hash(path)

    execute("csv_excess_cell", "REFUSE", excess_cell)

    def missing_cell(path: Path) -> None:
        data = path.parent / "observations.csv"
        lines = data.read_text(encoding="utf-8").splitlines()
        lines[1] = lines[1].rsplit(",", 1)[0]
        data.write_text("\n".join(lines) + "\n", encoding="utf-8")
        update_data_hash(path)

    execute("csv_missing_cell", "REFUSE", missing_cell)

    def nonfinite(path: Path, token: str, column: int) -> None:
        data = path.parent / "observations.csv"
        rows = list(csv.reader(data.read_text(encoding="utf-8").splitlines()))
        rows[1][column] = token
        with data.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream, lineterminator="\n").writerows(rows)
        update_data_hash(path)

    execute("nan_core", "REFUSE", lambda path: nonfinite(path, "nan", 7))
    execute("positive_inf_core", "REFUSE", lambda path: nonfinite(path, "inf", 5))
    execute("negative_inf_extra", "REFUSE", lambda path: nonfinite(path, "-inf", 8))

    def duplicate_row(path: Path) -> None:
        data = path.parent / "observations.csv"
        lines = data.read_text(encoding="utf-8").splitlines()
        lines[2] = "r1," + lines[2].split(",", 1)[1]
        data.write_text("\n".join(lines) + "\n", encoding="utf-8")
        update_data_hash(path)

    execute("duplicate_row_id", "REFUSE", duplicate_row)

    def decreasing_event_time(path: Path) -> None:
        data = path.parent / "observations.csv"
        rows = list(csv.reader(data.read_text(encoding="utf-8").splitlines()))
        rows[1][5] = "2"
        rows[2][5] = "1"
        with data.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream, lineterminator="\n").writerows(rows)
        update_data_hash(path)

    execute("decreasing_time_same_record_event", "REFUSE", decreasing_event_time)
    execute("impossible_utc_date", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("protocol_frozen_at_utc", "2026-02-30T00:00:00Z")), None)[1])
    execute("unregistered_uri_scheme", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("source_uri", "file:///tmp/fake")), None)[1])
    execute("unregistered_scope", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("observation_scope", "ALL_PROVED")), None)[1])
    execute("full_scope_missing_stages", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ")), None)[1])
    execute("none_control_combined", "REFUSE", lambda path: (edit_manifest(path, lambda value: value.__setitem__("controls", ["NONE", "OPPOSITE_WRITE"])), None)[1])

    def duplicate_json_key(path: Path) -> None:
        text = path.read_text(encoding="utf-8")
        path.write_text(text.replace("{\n", "{\n  \"schema\": \"ATTACKER_VALUE\",\n", 1), encoding="utf-8")

    execute("duplicate_outer_json_key", "REFUSE", duplicate_json_key)

    def duplicate_nested_json_key(path: Path) -> None:
        text = path.read_text(encoding="utf-8")
        needle = '"units": {\n'
        path.write_text(text.replace(needle, needle + '    "time": "ATTACKER_VALUE",\n', 1), encoding="utf-8")

    execute("duplicate_nested_json_key", "REFUSE", duplicate_nested_json_key)

    def malformed_unclosed_quote(path: Path) -> None:
        manifest = load_manifest(path)
        manifest["extra_columns"] = [{"name": "note", "type": "string", "unit": "1"}]
        manifest["columns"] = CORE + ["note"]
        data = path.parent / "observations.csv"
        data.write_text(
            ",".join(CORE + ["note"]) + "\n"
            + "r1,sample-a,event-a,POST_WRITE,RECORD_READOUT,0,0,1,\"unterminated\n",
            encoding="utf-8",
        )
        manifest["data_sha256"] = digest(data)
        save_manifest(path, manifest)

    execute("malformed_csv_unclosed_quote", "REFUSE", malformed_unclosed_quote)

    def fragmented_full_scope(path: Path) -> None:
        rows = [
            ["r1", "sample-a", "event-before", "BEFORE", "ENVIRONMENT", "0", "0", "0"],
            ["r2", "sample-a", "event-write", "WRITE", "RECORD_READOUT", "0", "0", "0"],
            ["r3", "sample-b", "event-hold", "HOLD", "WRITER_CONTROL", "0", "0", "0"],
            ["r4", "sample-b", "event-read", "READ", "ENVIRONMENT", "0", "0", "0"],
        ]
        replace_data(path, rows)
        edit_manifest(
            path,
            lambda value: (
                value.__setitem__("sample_ids", ["sample-a", "sample-b"]),
                value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"),
            ),
        )

    execute("fragmented_false_full_scope_upgrade", "REFUSE", fragmented_full_scope)

    def normalized_trace_drift(path: Path) -> None:
        # Hashes remain self-consistent, but the normalized value no longer derives from raw.dat.
        data = path.parent / "observations.csv"
        rows = list(csv.reader(data.read_text(encoding="utf-8").splitlines()))
        rows[1][7] = "999999"
        with data.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream, lineterminator="\n").writerows(rows)
        update_data_hash(path)

    execute("normalized_value_not_joined_to_source", "REFUSE", normalized_trace_drift)

    def fake_actual(path: Path) -> None:
        edit_manifest(
            path,
            lambda value: (
                value.__setitem__("evidence_class", "ACTUAL_SURFACE_MEASUREMENT"),
                value.__setitem__("source_uri", "https://example.invalid/not-an-observation"),
            ),
        )

    execute("self_asserted_actual_origin", "ACCEPT", fake_actual)

    def fake_prospective(path: Path) -> None:
        edit_manifest(
            path,
            lambda value: (
                value.__setitem__("protocol_timing", "PROSPECTIVE_PRE_OUTCOME"),
                value.__setitem__("protocol_frozen_at_utc", "2000-01-01T00:00:00Z"),
            ),
        )

    execute("self_asserted_prospective_timing", "ACCEPT", fake_prospective)

    def manifest_symlink_escape(path: Path) -> Path:
        nominal = path.parent.parent / "nominal"
        nominal.mkdir()
        link = nominal / "manifest.json"
        link.symlink_to(path)
        return link

    execute("manifest_symlink_re_roots_bundle", "REFUSE", manifest_symlink_escape)

    def duplicate_hardlink(path: Path) -> None:
        original = path.parent / "raw.dat"
        alias = path.parent / "raw-hardlink.dat"
        os.link(original, alias)
        manifest = load_manifest(path)
        manifest["source_artifacts"].append(
            {"path": alias.name, "sha256": digest(alias), "media_type": "application/octet-stream"}
        )
        save_manifest(path, manifest)

    execute("duplicate_source_hardlink_alias", "REFUSE", duplicate_hardlink)

    def data_is_protocol(path: Path) -> None:
        manifest = load_manifest(path)
        data = path.parent / manifest["data_file"]
        manifest["protocol_path"] = manifest["data_file"]
        manifest["protocol_sha256"] = digest(data)
        save_manifest(path, manifest)

    execute("csv_data_accepted_as_protocol", "REFUSE", data_is_protocol)

    expected_failures = [case for case in cases if not case["expectation_met"]]
    contained_metadata_boundaries = {
        "self_asserted_actual_origin": "Certificate explicitly disclaims physical-origin authentication.",
        "self_asserted_prospective_timing": "Timing is copied metadata, not independently attested.",
    }
    return {
        "case_count": len(cases),
        "expectations_met": len(cases) - len(expected_failures),
        "unexpected_accept_or_refuse_count": len(expected_failures),
        "unexpected_cases": [case["name"] for case in expected_failures],
        "contained_metadata_boundaries": contained_metadata_boundaries,
        "cases": cases,
        "pass": not expected_failures,
    }


def main() -> int:
    initial_freeze = freeze_check()
    official = run_official_paths() if initial_freeze["pass"] else {"pass": False, "skipped": "freeze mismatch"}
    join = verify_actual_join() if initial_freeze["pass"] else {"pass": False, "skipped": "freeze mismatch"}
    mutations = mutation_suite() if initial_freeze["pass"] else {"pass": False, "skipped": "freeze mismatch"}
    final_freeze = freeze_check()
    vulnerabilities = mutations.get("unexpected_cases", [])
    actual_bundle_pass = all(
        item.get("pass", False) for item in (initial_freeze, official, join, final_freeze)
    )
    verdict = "NOT_REFUTED" if actual_bundle_pass and not vulnerabilities else "REFUTED"
    result = {
        "schema": "WAC_T53A_WORLD_INGEST_VERIFY_V001",
        "verdict": verdict,
        "default": "REFUTED",
        "actual_760_row_bundle_pass": actual_bundle_pass,
        "generic_contract_adversarial_pass": not vulnerabilities,
        "scientific_claim_scored": False,
        "record_formation_proof_authorized": False,
        "universal_claim_authorized": False,
        "public_urm_registration_authorized": False,
        "initial_freeze": initial_freeze,
        "official_paths": official,
        "independent_actual_trace_join": join,
        "mutations": mutations,
        "final_freeze": final_freeze,
        "refuting_findings": vulnerabilities,
    }
    (VERIFY / "RESULT.json").write_text(stable_json(result), encoding="utf-8")
    lines = [
        f"T53A WORLD INGEST VERIFIER: {verdict}",
        f"actual_760_row_bundle_pass={str(actual_bundle_pass).lower()}",
        f"generic_contract_adversarial_pass={str(not vulnerabilities).lower()}",
        f"mutation_cases={mutations.get('case_count', 0)}",
        f"unexpected_cases={','.join(vulnerabilities) if vulnerabilities else 'NONE'}",
        "scientific_claim_scored=false",
        "record_formation_proof_authorized=false",
        "universal_claim_authorized=false",
        "public_urm_registration_authorized=false",
    ]
    (VERIFY / "RESULT.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0 if verdict == "NOT_REFUTED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
