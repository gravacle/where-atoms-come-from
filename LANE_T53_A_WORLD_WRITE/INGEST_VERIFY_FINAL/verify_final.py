#!/usr/bin/env python3
"""Final scientific provenance/schema QA for the T-53A world-observation bundle.

Default verdict is REFUTED.  The verifier reconstructs the 760-row table independently,
checks the registered normalization, and submits closed-schema/custody substitutions to the
current input contract.  It performs no scientific scoring.
"""

from __future__ import annotations

import csv
from datetime import datetime
import hashlib
import io
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
    "model/world_observation.py": "41249b88d24772b46e2ef7e357e62f1906555510518b00c24c0f87127d0fa716",
    "model/run_world_observation.py": "7a10ad8ee91b5283c00f42492eac7dada2df2a434610f2f467abfe9b1c0fb922",
    "model/checks_world_observation.py": "cb76bf865868d3ebc79833dc44b34810461597668c16d1f9e86c80d01b29ef2b",
    "model/lakeshore_vsm.py": "9e0e26fb168ed4d1cb155f0af1b04b8e0a3aa6a31aae2292d80a169ca87eb097",
    "model/checks_lakeshore_vsm.py": "3be19c33b071bc3fcf8c5b2d7aa43f487f2a8909f38f4e1a2513a42cc5ce2519",
    "LANE_T53_A_WORLD_WRITE/OBSERVATION_PROTOCOL.md": "1d32ac55feced8e15d1183e8652dfe7f7d74d19aa1ac0b68fc8b288ae12c6b70",
    "LANE_T53_A_WORLD_WRITE/export_world_observation.py": "74f2952cedeef3c01b776ca3fa3632e56941d85513f585218f76fb12108bfc54",
    "LANE_T53_A_WORLD_WRITE/raw/SOURCE.json": "76d214c1e0cad7a3416d039fa92d0872565f736f885df5b9014678f9719a3e90",
    "LANE_T53_A_WORLD_WRITE/world_observation.csv": "82fca51d5f7923107763b01340a85842d5a06d770ecc88e8e518a696a2e3a891",
    "LANE_T53_A_WORLD_WRITE/world_observation.json": "c388f7cef2d53242be1bd9f151543953f40319aeb10bb6be3a85656d566df253",
    "LANE_T53_A_WORLD_WRITE/world_observation_certificate.json": "99f470574454fd8c4efdd0235a67956fb3d641ca0577d41cc1633a13f2212a93",
}

CORE = [
    "row_id", "record_id", "event_id", "stage", "role", "time", "coordinate", "value"
]
EXTRAS = ["writer_field_T", "source_time_s", "protocol_segment"]
PAIR_TOKEN = "-hys-dcd-forc1 "


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino


def render_json(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def freeze() -> dict[str, Any]:
    actual = {name: sha256(REPO / name) for name in FROZEN_HASHES}
    return {"expected": FROZEN_HASHES, "actual": actual, "pass": actual == FROZEN_HASHES}


def import_current_modules():
    sys.path.insert(0, str(MODEL))
    try:
        import world_observation  # type: ignore
        import lakeshore_vsm  # type: ignore
    finally:
        sys.path.pop(0)
    return world_observation, lakeshore_vsm


CONTRACT, LAKESHORE = import_current_modules()


def contract_result(path: Path) -> tuple[str, str, dict[str, Any] | None]:
    try:
        observation = CONTRACT.load_world_observation(path)
    except CONTRACT.ObservationRefusal as exc:
        return "REFUSE", str(exc), None
    rendered = CONTRACT.certificate_json(observation)
    return "ACCEPT", rendered, json.loads(rendered)


def run_official_checks() -> dict[str, Any]:
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    expected = {
        "checks_world_observation.py": "WORLD_OBSERVATION_CHECKS: 27/27 PASS",
        "checks_lakeshore_vsm.py": "LAKESHORE_VSM_CHECKS: 9/9 PASS",
    }
    checks = []
    for name in expected:
        process = subprocess.run(
            [sys.executable, "-B", str(MODEL / name)], cwd=REPO, env=environment,
            capture_output=True, check=False,
        )
        checks.append(
            {
                "name": name,
                "returncode": process.returncode,
                "stdout": process.stdout.decode("utf-8", errors="replace").strip(),
                "stderr": process.stderr.decode("utf-8", errors="replace").strip(),
            }
        )
    command = [
        sys.executable, "-B", str(MODEL / "run_world_observation.py"),
        str(LANE / "world_observation.json"),
    ]
    runs = [
        subprocess.run(command, cwd=REPO, env=environment, capture_output=True, check=False)
        for _ in range(3)
    ]
    outputs = [run.stdout for run in runs]
    stored = (LANE / "world_observation_certificate.json").read_bytes()
    certificate = json.loads(outputs[0]) if runs[0].returncode == 0 else {}
    safety_expected = {
        "row_count": 760,
        "scope_classification": "CONFIGURATION_EVIDENCE_ONLY",
        "scientific_verdict": "NONE_NOT_SCORED",
        "normalized_source_join": "VERIFIED",
        "custody_validation": "CONTENT_HASHED_NOT_PHYSICAL_ORIGIN_AUTHENTICATED",
        "independent_reproduction_attested": False,
        "record_formation_proof_authorized": False,
        "universal_claim_authorized": False,
        "public_urm_registration_authorized": False,
    }
    safety_actual = {
        "row_count": certificate.get("coverage", {}).get("row_count"),
        "scope_classification": certificate.get("scope_classification"),
        "scientific_verdict": certificate.get("scientific_verdict"),
        "normalized_source_join": certificate.get("normalization", {}).get(
            "normalized_source_join"
        ),
        "custody_validation": certificate.get("custody_validation"),
        "independent_reproduction_attested": certificate.get(
            "independent_reproduction_attested"
        ),
        "record_formation_proof_authorized": certificate.get(
            "record_formation_proof_authorized"
        ),
        "universal_claim_authorized": certificate.get("universal_claim_authorized"),
        "public_urm_registration_authorized": certificate.get(
            "public_urm_registration_authorized"
        ),
    }
    checks_pass = all(
        check["returncode"] == 0 and check["stdout"] == expected[check["name"]]
        for check in checks
    )
    return {
        "checks": checks,
        "cli_returncodes": [run.returncode for run in runs],
        "three_cli_runs_byte_identical": len(set(outputs)) == 1,
        "stored_certificate_byte_identical": outputs[0] == stored,
        "certificate_sha256": hashlib.sha256(outputs[0]).hexdigest(),
        "safety_fields": safety_actual,
        "safety_fields_pass": safety_actual == safety_expected,
        "pass": (
            checks_pass
            and all(run.returncode == 0 for run in runs)
            and len(set(outputs)) == 1
            and outputs[0] == stored
            and safety_actual == safety_expected
        ),
    }


def header_value(lines: list[str], prefix: str) -> str:
    values = [line.split(":", 1)[1].strip() for line in lines if line.startswith(prefix + ":")]
    if len(values) != 1:
        raise AssertionError(f"expected one {prefix!r}; found {len(values)}")
    return values[0]


def independent_trace(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    lines = payload.decode("latin-1").splitlines()
    tables = [index for index, line in enumerate(lines) if line.startswith("Step,Iteration,")]
    if len(tables) != 1:
        raise AssertionError(f"{path.name}: table count {len(tables)}")
    expected_header = (
        "Step,Iteration,Segment,Field (µ0H) [T],Moment (m) [A·m²],"
        "Time Stamp [s],,Field Status,Moment (m) Status,"
    )
    # The Latin-1 decode of the source bytes is stable; require the whole instrument header.
    if lines[tables[0]] != expected_header:
        raise AssertionError(f"{path.name}: data-table header differs")
    try:
        parsed = list(csv.reader(lines[tables[0] + 1 :], strict=True))
    except csv.Error as exc:
        raise AssertionError(f"{path.name}: malformed CSV: {exc}") from exc
    fields: list[float] = []
    moments: list[float] = []
    times: list[float] = []
    for line_number, row in enumerate(parsed, start=tables[0] + 2):
        if not row or all(not cell.strip() for cell in row):
            continue
        if len(row) != 9 or row[8] != "":
            raise AssertionError(f"{path.name}:{line_number}: source cell closure failed")
        values = [float(row[index]) for index in (3, 4, 5)]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError(f"{path.name}:{line_number}: nonfinite value")
        if row[6].strip() != "GOOD" or row[7].strip() != "GOOD":
            raise AssertionError(f"{path.name}:{line_number}: non-GOOD status")
        fields.append(values[0])
        moments.append(values[1])
        times.append(values[2])
    if len(fields) < 3 or any(right < left for left, right in zip(times, times[1:])):
        raise AssertionError(f"{path.name}: insufficient/nonmonotone trace")
    has_hys = "#HYSTERESIS MEASUREMENT" in lines
    has_dcd = "#REMANENCE CURVES MEASUREMENT" in lines
    if has_hys == has_dcd:
        raise AssertionError(f"{path.name}: ambiguous measurement kind")
    kind = "hys" if has_hys else "dcd"
    writer_off = None if has_hys else header_value(
        lines, "Measure moment at applied fields"
    ) == "False"
    start = datetime.strptime(header_value(lines, "START TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    finish = datetime.strptime(header_value(lines, "FINISH TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    if finish < start:
        raise AssertionError(f"{path.name}: finish precedes start")
    return {
        "path": path,
        "kind": kind,
        "writer_off": writer_off,
        "field": fields,
        "moment": moments,
        "time": times,
        "magnet": header_value(lines, "Magnet"),
        "instrument_sample": header_value(lines, "ID"),
        "start": start,
        "finish": finish,
    }


def zero_interpolant(x: list[float], y: list[float]) -> tuple[float, float, int]:
    crossings = [index for index in range(len(x) - 1) if x[index] * x[index + 1] < 0]
    if len(crossings) != 1:
        raise AssertionError(f"strict zero crossing count {len(crossings)}")
    index = crossings[0]
    fraction = -x[index] / (x[index + 1] - x[index])
    return y[index] + fraction * (y[index + 1] - y[index]), fraction, index


def independent_normalization(paths: list[Path]) -> bytes:
    traces = {path.name: independent_trace(path) for path in paths}
    grouped: dict[str, dict[str, dict[str, Any]]] = {}
    for name, trace in traces.items():
        sample, suffix = name.split(PAIR_TOKEN, 1)
        kind = {
            "Step 1 Hysteresis Measurement.csv": "hys",
            "Step 2 Remanence Curves.csv": "dcd",
        }.get(suffix)
        if kind is None or trace["kind"] != kind:
            raise AssertionError(f"{name}: filename/measurement mismatch")
        grouped.setdefault(sample, {})[kind] = trace
    if len(grouped) != 5 or any(set(pair) != {"hys", "dcd"} for pair in grouped.values()):
        raise AssertionError("expected five complete HYS/DCD pairs")
    rows: list[list[Any]] = []
    for sample in sorted(grouped):
        hys = grouped[sample]["hys"]
        dcd = grouped[sample]["dcd"]
        if hys["magnet"] != "EM7-CSB" or dcd["magnet"] != "EM7-CSB":
            raise AssertionError(f"{sample}: instrument mismatch")
        if not dcd["writer_off"]:
            raise AssertionError(f"{sample}: DCD is not writer-off")
        turn = min(range(len(hys["field"])), key=hys["field"].__getitem__)
        branches = [
            ("positive-history", slice(0, turn + 1), max(hys["field"])),
            ("negative-history", slice(turn, len(hys["field"])), min(hys["field"])),
        ]
        ordinal = 0
        for label, branch, writer_field in branches:
            fields = hys["field"][branch]
            moments = hys["moment"][branch]
            times = hys["time"][branch]
            moment, fraction, crossing = zero_interpolant(fields, moments)
            source_time = times[crossing] + fraction * (times[crossing + 1] - times[crossing])
            rows.append([
                f"{sample}:hys:{label}", sample, f"hys:{label}", "POST_WRITE",
                "RECORD_READOUT", ordinal, 0.0, moment, writer_field, source_time,
                "HYSTERESIS_ZERO_FIELD_INTERPOLANT",
            ])
            ordinal += 1
        for index, (field, moment, source_time) in enumerate(
            zip(dcd["field"], dcd["moment"], dcd["time"])
        ):
            rows.append([
                f"{sample}:dcd:{index:03d}", sample, f"dcd:{index:03d}",
                "POST_WRITE", "RECORD_READOUT", ordinal, 0.0, moment, field,
                source_time, "DCD_WRITER_OFF_READOUT",
            ])
            ordinal += 1
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(CORE + EXTRAS)
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def actual_data_integrity() -> dict[str, Any]:
    manifest_path = LANE / "world_observation.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    declared = manifest["normalization"]["source_paths"]
    expected = sorted(path.relative_to(LANE).as_posix() for path in RAW.glob("*.csv"))
    paths = [LANE / relative for relative in declared]
    registered = LAKESHORE.normalized_world_observation_csv(paths)
    independent = independent_normalization(paths)
    stored = (LANE / manifest["data_file"]).read_bytes()
    artifacts = {item["path"]: item for item in manifest["source_artifacts"]}
    source_registry = json.loads((RAW / "SOURCE.json").read_text(encoding="utf-8"))
    registered_hashes = {item["file"]: item["sha256"] for item in source_registry["selected_files"]}
    registry_pass = (
        set(registered_hashes) == {path.name for path in paths}
        and all(registered_hashes[path.name] == sha256(path) for path in paths)
    )
    all_roles = [
        manifest_path,
        LANE / manifest["protocol_path"],
        LANE / manifest["data_file"],
        *[LANE / item["path"] for item in manifest["source_artifacts"]],
    ]
    role_ids = [identity(path.resolve()) for path in all_roles]
    protocol_words = " ".join(
        (LANE / manifest["protocol_path"]).read_text(encoding="utf-8").split()
    )
    boundary_pass = all(
        phrase in protocol_words
        for phrase in (
            "retrospective, measurement-only, no scientific scoring rule",
            "WRITE_POST_ONLY",
            "neither this normalization nor a successful ingestion certificate is a record-formation proof",
        )
    )
    facts = {
        "declared_source_count": len(declared),
        "declared_sources_are_exact_ten_raw_paths": declared == expected,
        "all_declared_sources_are_hashed_artifacts": all(path in artifacts for path in declared),
        "source_registry_hashes_pass": registry_pass,
        "registered_normalization_exact_bytes": registered == stored,
        "independent_normalization_exact_bytes": independent == stored,
        "registered_and_independent_byte_identical": registered == independent,
        "normalized_sha256": hashlib.sha256(registered).hexdigest(),
        "all_file_role_identities_disjoint": len(role_ids) == len(set(role_ids)),
        "adapter_hash_bound": (
            manifest["normalization"]["adapter_sha256"] == sha256(MODEL / "lakeshore_vsm.py")
        ),
        "protocol_nonproof_boundary_pass": boundary_pass,
    }
    facts["pass"] = all(
        (
            facts["declared_source_count"] == 10,
            facts["declared_sources_are_exact_ten_raw_paths"],
            facts["all_declared_sources_are_hashed_artifacts"],
            facts["source_registry_hashes_pass"],
            facts["registered_normalization_exact_bytes"],
            facts["independent_normalization_exact_bytes"],
            facts["registered_and_independent_byte_identical"],
            facts["all_file_role_identities_disjoint"],
            facts["adapter_hash_bound"],
            facts["protocol_nonproof_boundary_pass"],
        )
    )
    return facts


def write_fixture(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    protocol = root / "PROTOCOL.txt"
    data = root / "observations.csv"
    source = root / "raw.csv"
    protocol.write_text(
        "# WAC WORLD OBSERVATION PROTOCOL V001\n\nfinal provenance fixture\n",
        encoding="utf-8",
    )
    columns = CORE + ["quality"]
    with data.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerow(["r1", "sample-a", "event-a", "POST_WRITE", "RECORD_READOUT", 0, 0, 2, 1])
        writer.writerow(["r2", "sample-a", "event-a", "POST_WRITE", "RECORD_READOUT", 1, 0, -2, 1])
    source.write_bytes(data.read_bytes())
    manifest = {
        "schema": "WAC_WORLD_OBSERVATION_V001",
        "surface_id": "synthetic-surface",
        "sample_ids": ["sample-a"],
        "run_id": "final-input-fixture",
        "instrument_id": "synthetic-instrument",
        "source_uri": "urn:wac:final-input-fixture",
        "evidence_class": "SYNTHETIC_TEST_ONLY",
        "observation_scope": "WRITE_POST_ONLY",
        "protocol_timing": "PROSPECTIVE_PRE_OUTCOME",
        "protocol_frozen_at_utc": "2026-08-22T00:00:00Z",
        "protocol_path": protocol.name,
        "protocol_sha256": sha256(protocol),
        "data_file": data.name,
        "data_sha256": sha256(data),
        "columns": columns,
        "units": {"time": "s", "coordinate": "1", "value": "1"},
        "extra_columns": [{"name": "quality", "type": "float", "unit": "1"}],
        "controls": ["OPPOSITE_WRITE"],
        "normalization": {
            "adapter_id": "DIRECT_CSV_COPY_TEST_ONLY_V001",
            "adapter_sha256": sha256(MODEL / "world_observation.py"),
            "source_paths": [source.name],
        },
        "source_artifacts": [
            {"path": source.name, "sha256": sha256(source), "media_type": "text/csv"}
        ],
    }
    path = root / "manifest.json"
    path.write_text(render_json(manifest), encoding="utf-8")
    return path


def read_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_manifest(path: Path, value: dict[str, Any]) -> None:
    path.write_text(render_json(value), encoding="utf-8")


def change_manifest(path: Path, operation: Callable[[dict[str, Any]], None]) -> None:
    value = read_manifest(path)
    operation(value)
    write_manifest(path, value)


def sync_direct_source(path: Path) -> None:
    value = read_manifest(path)
    data = path.parent / value["data_file"]
    source_name = value["normalization"]["source_paths"][0]
    source = path.parent / source_name
    source.write_bytes(data.read_bytes())
    value["data_sha256"] = sha256(data)
    for artifact in value["source_artifacts"]:
        if artifact["path"] == source_name:
            artifact["sha256"] = sha256(source)
    write_manifest(path, value)


def replace_rows(path: Path, rows: list[list[Any]]) -> None:
    value = read_manifest(path)
    data = path.parent / value["data_file"]
    with data.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(value["columns"])
        writer.writerows(rows)
    sync_direct_source(path)


def clone_actual(root: Path) -> Path:
    root.mkdir(parents=True)
    (root / "raw").mkdir()
    for name in (
        "OBSERVATION_PROTOCOL.md", "export_world_observation.py",
        "world_observation.csv", "world_observation.json",
    ):
        shutil.copy2(LANE / name, root / name)
    for source in RAW.iterdir():
        if source.is_file():
            shutil.copy2(source, root / "raw" / source.name)
    return root / "world_observation.json"


def update_actual_source_registry(manifest_path: Path, relative: str) -> None:
    value = read_manifest(manifest_path)
    source_path = manifest_path.parent / relative
    source_digest = sha256(source_path)
    for artifact in value["source_artifacts"]:
        if artifact["path"] == relative:
            artifact["sha256"] = source_digest
    registry_path = manifest_path.parent / "raw" / "SOURCE.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    for item in registry["selected_files"]:
        if item["file"] == source_path.name:
            payload = source_path.read_bytes()
            item["sha256"] = source_digest
            item["bytes"] = len(payload)
            item["md5_repository"] = hashlib.md5(payload).hexdigest()
    registry_path.write_text(render_json(registry), encoding="utf-8")
    for artifact in value["source_artifacts"]:
        if artifact["path"] == "raw/SOURCE.json":
            artifact["sha256"] = sha256(registry_path)
    write_manifest(manifest_path, value)


def validation_cases() -> dict[str, Any]:
    scratch = VERIFY / "_scratch"
    scratch.mkdir(exist_ok=True)
    cases: list[dict[str, Any]] = []

    def run_case(
        name: str,
        group: str,
        expected: str,
        operation: Callable[[Path], Path | None],
        *,
        actual_bundle: bool = False,
    ) -> None:
        with tempfile.TemporaryDirectory(prefix=name + "-", dir=scratch) as temporary:
            directory = Path(temporary)
            manifest = (
                clone_actual(directory / "bundle")
                if actual_bundle
                else write_fixture(directory / "bundle")
            )
            returned = operation(manifest)
            target = returned if returned is not None else manifest
            observed, detail, certificate = contract_result(target)
            cases.append(
                {
                    "name": name,
                    "group": group,
                    "expected": expected,
                    "observed": observed,
                    "expectation_met": observed == expected,
                    "detail": detail[:600],
                    "scope_classification": (
                        certificate.get("scope_classification") if certificate else None
                    ),
                }
            )

    def manifest_edit(operation: Callable[[dict[str, Any]], None]):
        return lambda path: (change_manifest(path, operation), None)[1]

    run_case("valid_baseline", "control", "ACCEPT", lambda path: None)
    run_case("open_outer_manifest", "closure", "REFUSE", manifest_edit(lambda value: value.__setitem__("unexpected", True)))
    run_case("open_source_artifact_schema", "closure", "REFUSE", manifest_edit(lambda value: value["source_artifacts"][0].__setitem__("unexpected", True)))
    run_case("open_extra_column_schema", "closure", "REFUSE", manifest_edit(lambda value: value["extra_columns"][0].__setitem__("unexpected", True)))
    run_case("open_units_schema", "closure", "REFUSE", manifest_edit(lambda value: value["units"].__setitem__("unexpected", "1")))
    run_case("open_normalization_schema", "closure", "REFUSE", manifest_edit(lambda value: value["normalization"].__setitem__("unexpected", True)))
    run_case("parent_path_component", "path", "REFUSE", manifest_edit(lambda value: value.__setitem__("protocol_path", "../PROTOCOL.txt")))
    run_case("absolute_protocol_path", "path", "REFUSE", manifest_edit(lambda value: value.__setitem__("protocol_path", str((Path("/") / "tmp" / "PROTOCOL.txt")))))

    def source_link_outside(path: Path) -> None:
        outside = path.parent.parent / "outside.csv"
        outside.write_bytes(b"outside\n")
        link = path.parent / "outside-link.csv"
        link.symlink_to(outside)
        change_manifest(
            path,
            lambda value: value["source_artifacts"].__setitem__(
                0, {"path": link.name, "sha256": sha256(outside), "media_type": "text/csv"}
            ),
        )

    run_case("source_file_symlink_outside", "path", "REFUSE", source_link_outside)

    def nested_link_outside(path: Path) -> None:
        outside = path.parent.parent / "external-directory"
        outside.mkdir()
        payload = outside / "source.csv"
        payload.write_bytes(b"external\n")
        (path.parent / "linked-directory").symlink_to(outside, target_is_directory=True)
        change_manifest(
            path,
            lambda value: value["source_artifacts"].__setitem__(
                0,
                {"path": "linked-directory/source.csv", "sha256": sha256(payload), "media_type": "text/csv"},
            ),
        )

    run_case("nested_directory_symlink_outside", "path", "REFUSE", nested_link_outside)

    def duplicate_symlink_source(path: Path) -> None:
        alias = path.parent / "raw-alias.csv"
        alias.symlink_to(path.parent / "raw.csv")
        value = read_manifest(path)
        value["source_artifacts"].append(
            {"path": alias.name, "sha256": sha256(alias), "media_type": "text/csv"}
        )
        write_manifest(path, value)

    run_case("duplicate_source_symlink_alias", "identity", "REFUSE", duplicate_symlink_source)

    run_case("data_hash_drift", "hash", "REFUSE", lambda path: ((path.parent / "observations.csv").write_text("changed\n", encoding="utf-8"), None)[1])
    run_case("protocol_hash_drift", "hash", "REFUSE", lambda path: ((path.parent / "PROTOCOL.txt").write_text("changed\n", encoding="utf-8"), None)[1])
    run_case("source_hash_drift", "hash", "REFUSE", lambda path: ((path.parent / "raw.csv").write_text("changed\n", encoding="utf-8"), None)[1])

    def edit_data(path: Path, edit: Callable[[list[list[str]]], None]) -> None:
        data = path.parent / "observations.csv"
        rows = list(csv.reader(data.read_text(encoding="utf-8").splitlines()))
        edit(rows)
        with data.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream, lineterminator="\n").writerows(rows)
        sync_direct_source(path)

    run_case("wrong_csv_header", "csv", "REFUSE", lambda path: edit_data(path, lambda rows: rows[0].__setitem__(-1, "wrong")))
    run_case("csv_excess_cell", "csv", "REFUSE", lambda path: edit_data(path, lambda rows: rows[1].append("extra")))
    run_case("csv_missing_cell", "csv", "REFUSE", lambda path: edit_data(path, lambda rows: rows[1].pop()))
    run_case("nan_core_value", "numeric", "REFUSE", lambda path: edit_data(path, lambda rows: rows[1].__setitem__(7, "nan")))
    run_case("positive_infinite_time", "numeric", "REFUSE", lambda path: edit_data(path, lambda rows: rows[1].__setitem__(5, "inf")))
    run_case("negative_infinite_extra", "numeric", "REFUSE", lambda path: edit_data(path, lambda rows: rows[1].__setitem__(8, "-inf")))
    run_case("duplicate_row_id", "identity", "REFUSE", lambda path: edit_data(path, lambda rows: rows[2].__setitem__(0, rows[1][0])))
    run_case("decreasing_time_same_record_event", "time", "REFUSE", lambda path: edit_data(path, lambda rows: (rows[1].__setitem__(5, "2"), rows[2].__setitem__(5, "1"))))
    run_case("impossible_utc_date", "time", "REFUSE", manifest_edit(lambda value: value.__setitem__("protocol_frozen_at_utc", "2026-02-30T00:00:00Z")))
    run_case("unregistered_uri_scheme", "schema", "REFUSE", manifest_edit(lambda value: value.__setitem__("source_uri", "file:///not-a-source")))
    run_case("unregistered_scope", "schema", "REFUSE", manifest_edit(lambda value: value.__setitem__("observation_scope", "ALL_PROVED")))
    run_case("full_scope_missing_stages", "trajectory", "REFUSE", manifest_edit(lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ")))
    run_case("none_control_combined", "schema", "REFUSE", manifest_edit(lambda value: value.__setitem__("controls", ["NONE", "OPPOSITE_WRITE"])))

    def duplicate_outer(path: Path) -> None:
        text = path.read_text(encoding="utf-8")
        path.write_text(text.replace("{\n", '{\n  "schema": "DUPLICATE",\n', 1), encoding="utf-8")

    run_case("duplicate_outer_json_member", "json", "REFUSE", duplicate_outer)

    def duplicate_nested(path: Path) -> None:
        text = path.read_text(encoding="utf-8")
        marker = '"units": {\n'
        path.write_text(text.replace(marker, marker + '    "time": "DUPLICATE",\n', 1), encoding="utf-8")

    run_case("duplicate_nested_json_member", "json", "REFUSE", duplicate_nested)

    def malformed_normalized(path: Path) -> None:
        value = read_manifest(path)
        value["extra_columns"] = [{"name": "note", "type": "string", "unit": "1"}]
        value["columns"] = CORE + ["note"]
        data = path.parent / value["data_file"]
        data.write_text(
            ",".join(value["columns"]) + "\n"
            + 'r1,sample-a,event-a,POST_WRITE,RECORD_READOUT,0,0,1,"unterminated\n',
            encoding="utf-8",
        )
        write_manifest(path, value)
        sync_direct_source(path)

    run_case("malformed_normalized_csv_quote", "csv", "REFUSE", malformed_normalized)

    def fragmented_full(path: Path) -> None:
        replace_rows(
            path,
            [
                ["f1", "sample-a", "before-event", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1],
                ["f2", "sample-a", "write-event", "WRITE", "WRITER_CONTROL", 1, 0, 1, 1],
                ["f3", "sample-a", "hold-event", "HOLD", "ENVIRONMENT", 2, 0, 1, 1],
                ["f4", "sample-a", "read-event", "READ", "RECORD_READOUT", 3, 0, 1, 1],
            ],
        )
        change_manifest(path, lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"))

    run_case("fragmented_full_trajectory", "trajectory", "REFUSE", fragmented_full)

    def normalized_drift(path: Path) -> None:
        value = read_manifest(path)
        data = path.parent / value["data_file"]
        rows = list(csv.reader(data.read_text(encoding="utf-8").splitlines()))
        rows[1][7] = "999999"
        with data.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream, lineterminator="\n").writerows(rows)
        value["data_sha256"] = sha256(data)
        write_manifest(path, value)

    run_case("normalized_value_without_source_derivation", "normalization", "REFUSE", normalized_drift)
    run_case(
        "test_adapter_cannot_self_declare_actual",
        "declared_boundary",
        "REFUSE",
        manifest_edit(
            lambda value: value.__setitem__(
                "evidence_class", "ACTUAL_SURFACE_MEASUREMENT"
            )
        ),
    )

    run_case("self_declared_prospective_timing_boundary", "declared_boundary", "ACCEPT", manifest_edit(lambda value: (value.__setitem__("protocol_timing", "PROSPECTIVE_PRE_OUTCOME"), value.__setitem__("protocol_frozen_at_utc", "2000-01-01T00:00:00Z"))))

    def manifest_file_link(path: Path) -> Path:
        nominal = path.parent.parent / "nominal"
        nominal.mkdir()
        link = nominal / "manifest.json"
        link.symlink_to(path)
        return link

    run_case("linked_manifest_file", "path", "REFUSE", manifest_file_link)

    def duplicate_hardlink(path: Path) -> None:
        source = path.parent / "raw.csv"
        alias = path.parent / "raw-hardlink.csv"
        os.link(source, alias)
        value = read_manifest(path)
        value["source_artifacts"].append(
            {"path": alias.name, "sha256": sha256(alias), "media_type": "text/csv"}
        )
        write_manifest(path, value)

    run_case("duplicate_source_hardlink", "identity", "REFUSE", duplicate_hardlink)

    def protocol_as_data(path: Path) -> None:
        value = read_manifest(path)
        protocol = path.parent / value["protocol_path"]
        value["data_file"] = value["protocol_path"]
        value["data_sha256"] = sha256(protocol)
        write_manifest(path, value)

    run_case("protocol_and_data_same_file", "identity", "REFUSE", protocol_as_data)

    def invalid_protocol_marker(path: Path) -> None:
        protocol = path.parent / "PROTOCOL.txt"
        protocol.write_text("# UNREGISTERED PROTOCOL\n", encoding="utf-8")
        change_manifest(path, lambda value: value.__setitem__("protocol_sha256", sha256(protocol)))

    run_case("invalid_protocol_marker", "protocol", "REFUSE", invalid_protocol_marker)

    def valid_joined_full(path: Path) -> None:
        replace_rows(
            path,
            [
                ["j1", "sample-a", "trial-1", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1],
                ["j2", "sample-a", "trial-1", "WRITE", "WRITER_CONTROL", 1, 0, 1, 1],
                ["j3", "sample-a", "trial-1", "HOLD", "ENVIRONMENT", 2, 0, 1, 1],
                ["j4", "sample-a", "trial-1", "READ", "RECORD_READOUT", 3, 0, 1, 1],
            ],
        )
        change_manifest(path, lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"))

    run_case("valid_joined_full_trajectory", "trajectory", "ACCEPT", valid_joined_full)

    def missing_writer(path: Path) -> None:
        valid_joined_full(path)
        data = path.parent / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("WRITE,WRITER_CONTROL", "WRITE,ENVIRONMENT"), encoding="utf-8")
        sync_direct_source(path)

    run_case("joined_full_missing_writer_role", "trajectory", "REFUSE", missing_writer)

    def missing_read(path: Path) -> None:
        valid_joined_full(path)
        data = path.parent / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("READ,RECORD_READOUT", "READ,ENVIRONMENT"), encoding="utf-8")
        sync_direct_source(path)

    run_case("joined_full_missing_readout_role", "trajectory", "REFUSE", missing_read)

    def out_of_order(path: Path) -> None:
        valid_joined_full(path)
        data = path.parent / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("HOLD,ENVIRONMENT,2", "HOLD,ENVIRONMENT,4"), encoding="utf-8")
        sync_direct_source(path)

    run_case("joined_full_time_order_violation", "trajectory", "REFUSE", out_of_order)

    def incomplete_second_sample(path: Path) -> None:
        valid_joined_full(path)
        value = read_manifest(path)
        value["sample_ids"] = ["sample-a", "sample-b"]
        data = path.parent / "observations.csv"
        with data.open("a", encoding="utf-8", newline="") as stream:
            csv.writer(stream, lineterminator="\n").writerow(
                ["b1", "sample-b", "trial-b", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1]
            )
        write_manifest(path, value)
        sync_direct_source(path)

    run_case("full_second_sample_incomplete", "trajectory", "REFUSE", incomplete_second_sample)
    run_case("normalization_adapter_hash_drift", "normalization", "REFUSE", manifest_edit(lambda value: value["normalization"].__setitem__("adapter_sha256", "0" * 64)))
    run_case("unregistered_normalization_adapter", "normalization", "REFUSE", manifest_edit(lambda value: value["normalization"].__setitem__("adapter_id", "UNREGISTERED_ADAPTER")))
    run_case("normalization_source_not_in_artifacts", "normalization", "REFUSE", manifest_edit(lambda value: value["normalization"].__setitem__("source_paths", ["missing.csv"])))

    def source_alias_data(path: Path) -> None:
        value = read_manifest(path)
        data = path.parent / value["data_file"]
        alias = path.parent / "source-data-link.csv"
        os.link(data, alias)
        value["source_artifacts"] = [
            {"path": alias.name, "sha256": sha256(alias), "media_type": "text/csv"}
        ]
        value["normalization"]["source_paths"] = [alias.name]
        write_manifest(path, value)

    run_case("source_and_data_hardlink_identity", "identity", "REFUSE", source_alias_data)

    def linked_parent(path: Path) -> Path:
        alias = path.parent.parent / "linked-bundle"
        alias.symlink_to(path.parent, target_is_directory=True)
        return alias / path.name

    run_case("linked_parent_directory", "second_round", "REFUSE", linked_parent)

    def malformed_raw(path: Path) -> None:
        value = read_manifest(path)
        relative = value["normalization"]["source_paths"][0]
        source = path.parent / relative
        lines = source.read_bytes().splitlines()
        lines[-1] += b'"unterminated'
        source.write_bytes(b"\r\n".join(lines) + b"\r\n")
        update_actual_source_registry(path, relative)

    run_case("malformed_raw_lakeshore_csv", "second_round", "REFUSE", malformed_raw, actual_bundle=True)
    run_case("declared_sample_without_rows", "second_round", "REFUSE", manifest_edit(lambda value: value.__setitem__("sample_ids", ["sample-a", "ghost-sample"])))

    def rate_missing_sample(path: Path) -> None:
        replace_rows(
            path,
            [["rate-a", "sample-a", "rate-event", "RATE", "RECORD_READOUT", 0, 0, 1, 1]],
        )
        change_manifest(
            path,
            lambda value: (
                value.__setitem__("sample_ids", ["sample-a", "sample-b"]),
                value.__setitem__("observation_scope", "RATE_STABILITY_ONLY"),
            ),
        )

    run_case("rate_scope_missing_sample", "second_round", "REFUSE", rate_missing_sample)

    def rate_missing_for_observed_sample(path: Path) -> None:
        replace_rows(
            path,
            [
                ["rate-a", "sample-a", "rate-event", "RATE", "RECORD_READOUT", 0, 0, 1, 1],
                ["post-b", "sample-b", "post-event", "POST_WRITE", "RECORD_READOUT", 0, 0, 1, 1],
            ],
        )
        change_manifest(
            path,
            lambda value: (
                value.__setitem__("sample_ids", ["sample-a", "sample-b"]),
                value.__setitem__("observation_scope", "RATE_STABILITY_ONLY"),
            ),
        )

    run_case("rate_scope_observed_sample_lacks_rate", "independent", "REFUSE", rate_missing_for_observed_sample)

    def raw_extra_cell(path: Path) -> None:
        value = read_manifest(path)
        relative = value["normalization"]["source_paths"][0]
        source = path.parent / relative
        lines = source.read_bytes().splitlines()
        lines[-1] += b"EXTRA"
        source.write_bytes(b"\r\n".join(lines) + b"\r\n")
        update_actual_source_registry(path, relative)

    run_case("raw_lakeshore_nonempty_trailing_cell", "independent", "REFUSE", raw_extra_cell, actual_bundle=True)

    def raw_header_drift(path: Path) -> None:
        value = read_manifest(path)
        relative = value["normalization"]["source_paths"][0]
        source = path.parent / relative
        payload = source.read_bytes().replace(
            b"Moment (m) Status,\r\n", b"Moment (m) Status,Unexpected Column\r\n", 1
        )
        source.write_bytes(payload)
        update_actual_source_registry(path, relative)

    run_case("raw_lakeshore_table_header_drift", "independent", "REFUSE", raw_header_drift, actual_bundle=True)

    failures = [case for case in cases if not case["expectation_met"]]
    groups: dict[str, dict[str, Any]] = {}
    for group in sorted({case["group"] for case in cases}):
        members = [case for case in cases if case["group"] == group]
        groups[group] = {
            "count": len(members),
            "pass": all(case["expectation_met"] for case in members),
        }
    return {
        "case_count": len(cases),
        "expectations_met": len(cases) - len(failures),
        "unexpected_case_count": len(failures),
        "first_surviving_case": failures[0]["name"] if failures else None,
        "unexpected_cases": [case["name"] for case in failures],
        "groups": groups,
        "cases": cases,
        "pass": not failures,
    }


def write_d24(result: dict[str, Any]) -> None:
    cases = result["validation_cases"]
    unexpected = cases["unexpected_cases"]
    lines = [
        "# D24 — Final T-53A scientific data-integrity recheck",
        "",
        "Posture: **REFUTED by default**. This is scientific provenance/schema QA, not a",
        "scientific verdict and not a record-formation proof gate.",
        "",
        "## Frozen implementation",
        "",
    ]
    for name, value in sorted(FROZEN_HASHES.items()):
        lines.append(f"- `{name}`: `{value}`")
    lines.extend(
        [
            "",
            "## Actual-bundle reconstruction",
            "",
            f"- Actual bundle pass: `{str(result['actual_760_row_bundle_pass']).lower()}`.",
            "- Ten declared Lake Shore sources and five sample pairs are reconstructed by an",
            "  independent strict parser and by the registered adapter.",
            "- Both outputs must match the stored 760-row CSV byte for byte.",
            "- Protocol, data, manifest, and source artifact file identities must be disjoint.",
            "- Certificate remains `CONFIGURATION_EVIDENCE_ONLY`, `NONE_NOT_SCORED`, and",
            "  authorizes no proof, universal claim, independent reproduction, or public URM registration.",
            "",
            "## Verifier correction log",
            "",
            "- Initial run stopped before issuing a verdict because the independent parser's",
            "  literal Latin-1 table header mistakenly expected double-encoded unit symbols.",
            "- The verifier literal was corrected to the source bytes' `µ` and `·` characters;",
            "  no model, source, normalized data, manifest, or certificate file was changed.",
            "",
            "## Input-validation cases",
            "",
            f"- Cases: `{cases['case_count']}`.",
            f"- Expectations met: `{cases['expectations_met']}`.",
            f"- First surviving case: `{cases['first_surviving_case'] or 'NONE'}`.",
            "",
        ]
    )
    if unexpected:
        lines.extend(["## Unresolved findings", ""])
        case_map = {case["name"]: case for case in cases["cases"]}
        for name in unexpected:
            case = case_map[name]
            lines.extend(
                [
                    f"### {name}",
                    "",
                    f"- Expected `{case['expected']}`; observed `{case['observed']}`.",
                    f"- Group: `{case['group']}`.",
                    "- The default-refuted verifier therefore does not clear the reusable input contract.",
                    "",
                ]
            )
    else:
        lines.extend(
            [
                "## Disposition",
                "",
                "No declared malformed/substitution case survived. The input contract is",
                "`NOT_REFUTED` at this verifier boundary. This says only that the frozen actual",
                "bundle and the declared local provenance/schema cases pass; it adds no scientific",
                "measurement and carries no record-formation or universality authorization.",
                "",
            ]
        )
    (VERIFY / "D24.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    initial = freeze()
    official = run_official_checks() if initial["pass"] else {"pass": False, "skipped": "freeze mismatch"}
    integrity = actual_data_integrity() if initial["pass"] else {"pass": False, "skipped": "freeze mismatch"}
    cases = validation_cases() if initial["pass"] else {"pass": False, "skipped": "freeze mismatch"}
    final = freeze()
    actual_pass = all(item.get("pass", False) for item in (initial, official, integrity, final))
    case_pass = cases.get("pass", False)
    verdict = "NOT_REFUTED" if actual_pass and case_pass else "REFUTED"
    result = {
        "schema": "WAC_T53A_WORLD_INGEST_FINAL_QA_V001",
        "default": "REFUTED",
        "verdict": verdict,
        "verifier_source_sha256": sha256(Path(__file__).resolve()),
        "actual_760_row_bundle_pass": actual_pass,
        "input_validation_pass": case_pass,
        "first_surviving_case": cases.get("first_surviving_case"),
        "scientific_claim_scored": False,
        "record_formation_proof_authorized": False,
        "universal_claim_authorized": False,
        "public_urm_registration_authorized": False,
        "initial_freeze": initial,
        "official_checks": official,
        "actual_data_integrity": integrity,
        "validation_cases": cases,
        "final_freeze": final,
    }
    (VERIFY / "RESULT.json").write_text(render_json(result), encoding="utf-8")
    lines = [
        f"T53A WORLD INGEST FINAL QA: {verdict}",
        f"actual_760_row_bundle_pass={str(actual_pass).lower()}",
        f"input_validation_pass={str(case_pass).lower()}",
        f"case_count={cases.get('case_count', 0)}",
        f"first_surviving_case={cases.get('first_surviving_case') or 'NONE'}",
        f"unexpected_cases={','.join(cases.get('unexpected_cases', [])) or 'NONE'}",
        "scientific_claim_scored=false",
        "record_formation_proof_authorized=false",
        "universal_claim_authorized=false",
        "public_urm_registration_authorized=false",
    ]
    (VERIFY / "RESULT.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_d24(result)
    print("\n".join(lines))
    return 0 if verdict == "NOT_REFUTED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
