#!/usr/bin/env python3
"""Fresh closure audit of the frozen T-53A world-observation input contract."""

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


HERE = Path(__file__).resolve().parent
LANE = HERE.parent
REPO = LANE.parent
MODEL = REPO / "model"
RAW = LANE / "raw"
CORE = ["row_id", "record_id", "event_id", "stage", "role", "time", "coordinate", "value"]
EXTRAS = ["writer_field_T", "source_time_s", "protocol_segment"]
PAIR = "-hys-dcd-forc1 "
PINS = {
    "model/world_observation.py": "8e511fad0862afa4d48be24b29feca688c4817418e6a83c614d44ca8e760b81f",
    "model/run_world_observation.py": "7a10ad8ee91b5283c00f42492eac7dada2df2a434610f2f467abfe9b1c0fb922",
    "model/checks_world_observation.py": "76a96bda992bfcf83345190175eff8b9f45030caefd58c92c3f3c3ccad88746d",
    "model/lakeshore_vsm.py": "68221ff4deab5442370e4be57db44916810dd534676d6ead39f378e905370bf1",
    "model/checks_lakeshore_vsm.py": "ee2b79c12921123ad08e0f8d19cbfd8884e129743692ee67e1630b1565de8ffa",
    "LANE_T53_A_WORLD_WRITE/OBSERVATION_PROTOCOL.md": "1d32ac55feced8e15d1183e8652dfe7f7d74d19aa1ac0b68fc8b288ae12c6b70",
    "LANE_T53_A_WORLD_WRITE/export_world_observation.py": "74f2952cedeef3c01b776ca3fa3632e56941d85513f585218f76fb12108bfc54",
    "LANE_T53_A_WORLD_WRITE/raw/SOURCE.json": "76d214c1e0cad7a3416d039fa92d0872565f736f885df5b9014678f9719a3e90",
    "LANE_T53_A_WORLD_WRITE/world_observation.csv": "82fca51d5f7923107763b01340a85842d5a06d770ecc88e8e518a696a2e3a891",
    "LANE_T53_A_WORLD_WRITE/world_observation.json": "15352eaa9c6363c5ff16c02c7c75887941999c9ab914c758c4218ff25b7e1d20",
    "LANE_T53_A_WORLD_WRITE/world_observation_certificate.json": "04beada93f26d2c7f4d3e7b88f016e4b4f27c70e38d8f8812f7e352c3eba0faf",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ident(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino


def json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def pin_check() -> dict[str, Any]:
    actual = {name: digest(REPO / name) for name in PINS}
    return {"expected": PINS, "actual": actual, "pass": actual == PINS}


sys.path.insert(0, str(MODEL))
try:
    import world_observation as contract  # type: ignore
    import lakeshore_vsm as lake  # type: ignore
finally:
    sys.path.pop(0)


def load(path: Path) -> tuple[str, str, dict[str, Any] | None]:
    try:
        observation = contract.load_world_observation(path)
    except contract.ObservationRefusal as exc:
        return "REFUSE", str(exc), None
    rendered = contract.certificate_json(observation)
    return "ACCEPT", rendered, json.loads(rendered)


def official() -> dict[str, Any]:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    expected = {
        "checks_world_observation.py": "WORLD_OBSERVATION_CHECKS: 28/28 PASS",
        "checks_lakeshore_vsm.py": "LAKESHORE_VSM_CHECKS: 11/11 PASS",
    }
    checks = []
    for name, line in expected.items():
        run = subprocess.run(
            [sys.executable, "-B", str(MODEL / name)], cwd=REPO, env=env,
            capture_output=True, check=False,
        )
        checks.append({
            "name": name, "returncode": run.returncode,
            "stdout": run.stdout.decode(errors="replace").strip(),
            "stderr": run.stderr.decode(errors="replace").strip(),
            "expected": line,
        })
    command = [sys.executable, "-B", str(MODEL / "run_world_observation.py"), str(LANE / "world_observation.json")]
    runs = [subprocess.run(command, cwd=REPO, env=env, capture_output=True, check=False) for _ in range(3)]
    outputs = [run.stdout for run in runs]
    stored = (LANE / "world_observation_certificate.json").read_bytes()
    cert = json.loads(outputs[0]) if runs[0].returncode == 0 else {}
    safety = {
        "row_count": cert.get("coverage", {}).get("row_count"),
        "scope_classification": cert.get("scope_classification"),
        "scientific_verdict": cert.get("scientific_verdict"),
        "normalized_source_join": cert.get("normalization", {}).get("normalized_source_join"),
        "custody_validation": cert.get("custody_validation"),
        "independent_reproduction_attested": cert.get("independent_reproduction_attested"),
        "record_formation_proof_authorized": cert.get("record_formation_proof_authorized"),
        "universal_claim_authorized": cert.get("universal_claim_authorized"),
        "public_urm_registration_authorized": cert.get("public_urm_registration_authorized"),
    }
    wanted = {
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
    passed = (
        all(item["returncode"] == 0 and item["stdout"] == item["expected"] for item in checks)
        and all(run.returncode == 0 for run in runs)
        and len(set(outputs)) == 1 and outputs[0] == stored and safety == wanted
    )
    return {
        "checks": checks, "cli_returncodes": [run.returncode for run in runs],
        "three_runs_identical": len(set(outputs)) == 1,
        "stored_certificate_identical": outputs[0] == stored,
        "certificate_sha256": hashlib.sha256(outputs[0]).hexdigest(),
        "safety": safety, "safety_pass": safety == wanted, "pass": passed,
    }


def one(lines: list[str], prefix: str) -> str:
    hits = [line.split(":", 1)[1].strip() for line in lines if line.startswith(prefix + ":")]
    if len(hits) != 1:
        raise AssertionError(f"{prefix}: {len(hits)} occurrences")
    return hits[0]


def parse_raw(path: Path) -> dict[str, Any]:
    lines = path.read_bytes().decode("latin-1").splitlines()
    header = "Step,Iteration,Segment,Field (µ0H) [T],Moment (m) [A·m²],Time Stamp [s],,Field Status,Moment (m) Status,"
    tables = [i for i, line in enumerate(lines) if line == header]
    if len(tables) != 1:
        raise AssertionError(f"{path.name}: table closure")
    try:
        rows = list(csv.reader(lines[tables[0] + 1:], strict=True))
    except csv.Error as exc:
        raise AssertionError(f"{path.name}: malformed source CSV") from exc
    field: list[float] = []
    moment: list[float] = []
    time: list[float] = []
    for row in rows:
        if not row or all(not cell.strip() for cell in row):
            continue
        if len(row) != 9 or row[8] != "" or row[6].strip() != "GOOD" or row[7].strip() != "GOOD":
            raise AssertionError(f"{path.name}: data row closure")
        triple = [float(row[i]) for i in (3, 4, 5)]
        if not all(math.isfinite(item) for item in triple):
            raise AssertionError(f"{path.name}: nonfinite")
        field.append(triple[0]); moment.append(triple[1]); time.append(triple[2])
    if any(b < a for a, b in zip(time, time[1:])):
        raise AssertionError(f"{path.name}: time order")
    hys = "#HYSTERESIS MEASUREMENT" in lines
    dcd = "#REMANENCE CURVES MEASUREMENT" in lines
    if hys == dcd:
        raise AssertionError(f"{path.name}: measurement identity")
    start = datetime.strptime(one(lines, "START TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    finish = datetime.strptime(one(lines, "FINISH TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    return {
        "kind": "hys" if hys else "dcd", "field": field, "moment": moment, "time": time,
        "writer_off": None if hys else one(lines, "Measure moment at applied fields") == "False",
        "magnet": one(lines, "Magnet"), "sample": one(lines, "ID"),
        "start": start, "finish": finish,
    }


def crossing(x: list[float], y: list[float]) -> tuple[float, float, int]:
    found = [i for i in range(len(x) - 1) if x[i] * x[i + 1] < 0]
    if len(found) != 1:
        raise AssertionError(f"zero crossings: {len(found)}")
    i = found[0]
    fraction = -x[i] / (x[i + 1] - x[i])
    return y[i] + fraction * (y[i + 1] - y[i]), fraction, i


def independent_csv(paths: list[Path]) -> bytes:
    traces = {path.name: parse_raw(path) for path in paths}
    groups: dict[str, dict[str, dict[str, Any]]] = {}
    suffixes = {"Step 1 Hysteresis Measurement.csv": "hys", "Step 2 Remanence Curves.csv": "dcd"}
    for name, trace in traces.items():
        sample, suffix = name.split(PAIR, 1)
        kind = suffixes.get(suffix)
        if kind is None or trace["kind"] != kind:
            raise AssertionError(f"{name}: filename/header mismatch")
        groups.setdefault(sample, {})[kind] = trace
    if len(groups) != 5 or any(set(pair) != {"hys", "dcd"} for pair in groups.values()):
        raise AssertionError("five complete pairs required")
    result: list[list[Any]] = []
    for sample in sorted(groups):
        h = groups[sample]["hys"]; d = groups[sample]["dcd"]
        if h["magnet"] != "EM7-CSB" or d["magnet"] != "EM7-CSB" or not d["writer_off"]:
            raise AssertionError(f"{sample}: instrument/writer-off mismatch")
        turn = min(range(len(h["field"])), key=h["field"].__getitem__)
        branches = [("positive-history", slice(0, turn + 1), max(h["field"])), ("negative-history", slice(turn, len(h["field"])), min(h["field"]))]
        ordinal = 0
        for label, branch, writer_field in branches:
            fields = h["field"][branch]; moments = h["moment"][branch]; times = h["time"][branch]
            value, fraction, i = crossing(fields, moments)
            source_time = times[i] + fraction * (times[i + 1] - times[i])
            result.append([f"{sample}:hys:{label}", sample, f"hys:{label}", "POST_WRITE", "RECORD_READOUT", ordinal, 0.0, value, writer_field, source_time, "HYSTERESIS_ZERO_FIELD_INTERPOLANT"])
            ordinal += 1
        for i, (writer_field, value, source_time) in enumerate(zip(d["field"], d["moment"], d["time"])):
            result.append([f"{sample}:dcd:{i:03d}", sample, f"dcd:{i:03d}", "POST_WRITE", "RECORD_READOUT", ordinal, 0.0, value, writer_field, source_time, "DCD_WRITER_OFF_READOUT"])
            ordinal += 1
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n"); writer.writerow(CORE + EXTRAS); writer.writerows(result)
    return stream.getvalue().encode()


def actual_bundle() -> dict[str, Any]:
    manifest_path = LANE / "world_observation.json"
    manifest = json.loads(manifest_path.read_text())
    declared = manifest["normalization"]["source_paths"]
    expected = sorted(path.relative_to(LANE).as_posix() for path in RAW.glob("*.csv"))
    paths = [LANE / name for name in declared]
    registered = lake.normalized_world_observation_csv(paths)
    independent = independent_csv(paths)
    stored = (LANE / manifest["data_file"]).read_bytes()
    artifacts = {item["path"]: item for item in manifest["source_artifacts"]}
    registry = json.loads((RAW / "SOURCE.json").read_text())
    registered_hashes = {item["file"]: item["sha256"] for item in registry["selected_files"]}
    roles = [manifest_path, LANE / manifest["protocol_path"], LANE / manifest["data_file"], *[LANE / item["path"] for item in manifest["source_artifacts"]]]
    facts = {
        "exact_ten_sources": declared == expected and len(declared) == 10,
        "sources_are_artifacts": all(name in artifacts for name in declared),
        "source_registry_matches": set(registered_hashes) == {path.name for path in paths} and all(registered_hashes[path.name] == digest(path) for path in paths),
        "registered_exact_bytes": registered == stored,
        "independent_exact_bytes": independent == stored,
        "registered_equals_independent": registered == independent,
        "all_role_identities_disjoint": len({ident(path.resolve()) for path in roles}) == len(roles),
        "adapter_hash_bound": manifest["normalization"]["adapter_sha256"] == digest(MODEL / "lakeshore_vsm.py"),
        "data_sha256": hashlib.sha256(stored).hexdigest(),
    }
    facts["pass"] = all(value for key, value in facts.items() if key not in {"data_sha256", "pass"})
    return facts


def fixture(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    protocol = root / "PROTOCOL.txt"; data = root / "observations.csv"; source = root / "raw.csv"
    protocol.write_text("# WAC WORLD OBSERVATION PROTOCOL V001\n\nclosure fixture\n")
    columns = CORE + ["quality"]
    with data.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n"); writer.writerow(columns)
        writer.writerow(["r1", "sample-a", "event-a", "POST_WRITE", "RECORD_READOUT", 0, 0, 2, 1])
        writer.writerow(["r2", "sample-a", "event-a", "POST_WRITE", "RECORD_READOUT", 1, 0, -2, 1])
    source.write_bytes(data.read_bytes())
    value = {
        "schema": "WAC_WORLD_OBSERVATION_V001", "surface_id": "synthetic-surface",
        "sample_ids": ["sample-a"], "run_id": "closure-fixture", "instrument_id": "synthetic-instrument",
        "source_uri": "urn:wac:closure-fixture", "evidence_class": "SYNTHETIC_TEST_ONLY",
        "observation_scope": "WRITE_POST_ONLY", "protocol_timing": "PROSPECTIVE_PRE_OUTCOME",
        "protocol_frozen_at_utc": "2026-08-22T00:00:00Z", "protocol_path": protocol.name,
        "protocol_sha256": digest(protocol), "data_file": data.name, "data_sha256": digest(data),
        "columns": columns, "units": {"time": "s", "coordinate": "1", "value": "1"},
        "extra_columns": [{"name": "quality", "type": "float", "unit": "1"}],
        "controls": ["OPPOSITE_WRITE"],
        "normalization": {"adapter_id": "DIRECT_CSV_COPY_TEST_ONLY_V001", "adapter_sha256": digest(MODEL / "world_observation.py"), "source_paths": [source.name]},
        "source_artifacts": [{"path": source.name, "sha256": digest(source), "media_type": "text/csv"}],
    }
    path = root / "manifest.json"; path.write_text(json_text(value)); return path


def read_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_manifest(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json_text(value))


def mutate(path: Path, operation: Callable[[dict[str, Any]], None]) -> None:
    value = read_manifest(path); operation(value); write_manifest(path, value)


def sync(path: Path) -> None:
    value = read_manifest(path); data = path.parent / value["data_file"]
    source_name = value["normalization"]["source_paths"][0]; source = path.parent / source_name
    source.write_bytes(data.read_bytes()); value["data_sha256"] = digest(data)
    for item in value["source_artifacts"]:
        if item["path"] == source_name: item["sha256"] = digest(source)
    write_manifest(path, value)


def rows(path: Path, content: list[list[Any]]) -> None:
    value = read_manifest(path); data = path.parent / value["data_file"]
    with data.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n"); writer.writerow(value["columns"]); writer.writerows(content)
    sync(path)


def clone_actual(root: Path) -> Path:
    root.mkdir(parents=True); (root / "raw").mkdir()
    for name in ("OBSERVATION_PROTOCOL.md", "export_world_observation.py", "world_observation.csv", "world_observation.json"):
        shutil.copy2(LANE / name, root / name)
    for path in RAW.iterdir():
        if path.is_file(): shutil.copy2(path, root / "raw" / path.name)
    return root / "world_observation.json"


def refresh_actual_source(manifest_path: Path, relative: str) -> None:
    value = read_manifest(manifest_path); source = manifest_path.parent / relative; payload = source.read_bytes()
    for item in value["source_artifacts"]:
        if item["path"] == relative: item["sha256"] = digest(source)
    registry_path = manifest_path.parent / "raw" / "SOURCE.json"; registry = json.loads(registry_path.read_text())
    for item in registry["selected_files"]:
        if item["file"] == source.name:
            item["sha256"] = digest(source); item["bytes"] = len(payload); item["md5_repository"] = hashlib.md5(payload).hexdigest()
    registry_path.write_text(json_text(registry))
    for item in value["source_artifacts"]:
        if item["path"] == "raw/SOURCE.json": item["sha256"] = digest(registry_path)
    write_manifest(manifest_path, value)


def replay() -> dict[str, Any]:
    scratch = HERE / "_scratch"; scratch.mkdir(exist_ok=True)
    results: list[dict[str, Any]] = []

    def case(name: str, expected: str, operation: Callable[[Path], Path | None], group: str = "prior54", actual: bool = False) -> None:
        with tempfile.TemporaryDirectory(prefix=name + "-", dir=scratch) as temporary:
            base = Path(temporary); manifest = clone_actual(base / "bundle") if actual else fixture(base / "bundle")
            returned = operation(manifest); target = returned if returned is not None else manifest
            observed, detail, certificate = load(target)
            results.append({"name": name, "group": group, "expected": expected, "observed": observed, "expectation_met": expected == observed, "detail": detail[:600], "scope_classification": certificate.get("scope_classification") if certificate else None})

    def edit(operation: Callable[[dict[str, Any]], None]):
        return lambda path: (mutate(path, operation), None)[1]

    def edit_csv(operation: Callable[[list[list[str]]], None]):
        def apply(path: Path) -> None:
            data = path.parent / "observations.csv"; content = list(csv.reader(data.read_text().splitlines())); operation(content)
            with data.open("w", newline="") as stream: csv.writer(stream, lineterminator="\n").writerows(content)
            sync(path)
        return apply

    case("valid_baseline", "ACCEPT", lambda path: None)
    case("open_outer_manifest", "REFUSE", edit(lambda v: v.__setitem__("unexpected", True)))
    case("open_source_artifact_schema", "REFUSE", edit(lambda v: v["source_artifacts"][0].__setitem__("unexpected", True)))
    case("open_extra_column_schema", "REFUSE", edit(lambda v: v["extra_columns"][0].__setitem__("unexpected", True)))
    case("open_units_schema", "REFUSE", edit(lambda v: v["units"].__setitem__("unexpected", "1")))
    case("open_normalization_schema", "REFUSE", edit(lambda v: v["normalization"].__setitem__("unexpected", True)))
    case("parent_path_component", "REFUSE", edit(lambda v: v.__setitem__("protocol_path", "../PROTOCOL.txt")))
    case("absolute_protocol_path", "REFUSE", edit(lambda v: v.__setitem__("protocol_path", "/tmp/PROTOCOL.txt")))

    def outside(path: Path) -> None:
        payload = path.parent.parent / "outside.csv"; payload.write_bytes(b"outside\n"); link = path.parent / "outside-link.csv"; link.symlink_to(payload)
        mutate(path, lambda v: v["source_artifacts"].__setitem__(0, {"path": link.name, "sha256": digest(payload), "media_type": "text/csv"}))
    case("source_file_symlink_outside", "REFUSE", outside)

    def nested(path: Path) -> None:
        directory = path.parent.parent / "external"; directory.mkdir(); payload = directory / "source.csv"; payload.write_bytes(b"outside\n")
        (path.parent / "linked").symlink_to(directory, target_is_directory=True)
        mutate(path, lambda v: v["source_artifacts"].__setitem__(0, {"path": "linked/source.csv", "sha256": digest(payload), "media_type": "text/csv"}))
    case("nested_directory_symlink_outside", "REFUSE", nested)

    def symlink_alias(path: Path) -> None:
        alias = path.parent / "raw-alias.csv"; alias.symlink_to(path.parent / "raw.csv"); value = read_manifest(path)
        value["source_artifacts"].append({"path": alias.name, "sha256": digest(alias), "media_type": "text/csv"}); write_manifest(path, value)
    case("duplicate_source_symlink_alias", "REFUSE", symlink_alias)
    case("data_hash_drift", "REFUSE", lambda path: ((path.parent / "observations.csv").write_text("changed\n"), None)[1])
    case("protocol_hash_drift", "REFUSE", lambda path: ((path.parent / "PROTOCOL.txt").write_text("changed\n"), None)[1])
    case("source_hash_drift", "REFUSE", lambda path: ((path.parent / "raw.csv").write_text("changed\n"), None)[1])
    case("wrong_csv_header", "REFUSE", edit_csv(lambda r: r[0].__setitem__(-1, "wrong")))
    case("csv_excess_cell", "REFUSE", edit_csv(lambda r: r[1].append("extra")))
    case("csv_missing_cell", "REFUSE", edit_csv(lambda r: r[1].pop()))
    case("nan_core_value", "REFUSE", edit_csv(lambda r: r[1].__setitem__(7, "nan")))
    case("positive_infinite_time", "REFUSE", edit_csv(lambda r: r[1].__setitem__(5, "inf")))
    case("negative_infinite_extra", "REFUSE", edit_csv(lambda r: r[1].__setitem__(8, "-inf")))
    case("duplicate_row_id", "REFUSE", edit_csv(lambda r: r[2].__setitem__(0, r[1][0])))
    case("decreasing_time_same_record_event", "REFUSE", edit_csv(lambda r: (r[1].__setitem__(5, "2"), r[2].__setitem__(5, "1"))))
    case("impossible_utc_date", "REFUSE", edit(lambda v: v.__setitem__("protocol_frozen_at_utc", "2026-02-30T00:00:00Z")))
    case("unregistered_uri_scheme", "REFUSE", edit(lambda v: v.__setitem__("source_uri", "file:///not-source")))
    case("unregistered_scope", "REFUSE", edit(lambda v: v.__setitem__("observation_scope", "ALL_PROVED")))
    case("full_scope_missing_stages", "REFUSE", edit(lambda v: v.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ")))
    case("none_control_combined", "REFUSE", edit(lambda v: v.__setitem__("controls", ["NONE", "OPPOSITE_WRITE"])))

    def dup_outer(path: Path) -> None: path.write_text(path.read_text().replace("{\n", '{\n  "schema": "DUPLICATE",\n', 1))
    case("duplicate_outer_json_member", "REFUSE", dup_outer)
    def dup_nested(path: Path) -> None: path.write_text(path.read_text().replace('"units": {\n', '"units": {\n    "time": "DUPLICATE",\n', 1))
    case("duplicate_nested_json_member", "REFUSE", dup_nested)

    def bad_quote(path: Path) -> None:
        value = read_manifest(path); value["extra_columns"] = [{"name": "note", "type": "string", "unit": "1"}]; value["columns"] = CORE + ["note"]
        (path.parent / "observations.csv").write_text(",".join(value["columns"]) + '\nr1,sample-a,event-a,POST_WRITE,RECORD_READOUT,0,0,1,"unterminated\n'); write_manifest(path, value); sync(path)
    case("malformed_normalized_csv_quote", "REFUSE", bad_quote)

    def fragmented(path: Path) -> None:
        rows(path, [["f1", "sample-a", "e1", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1], ["f2", "sample-a", "e2", "WRITE", "WRITER_CONTROL", 1, 0, 1, 1], ["f3", "sample-a", "e3", "HOLD", "ENVIRONMENT", 2, 0, 1, 1], ["f4", "sample-a", "e4", "READ", "RECORD_READOUT", 3, 0, 1, 1]])
        mutate(path, lambda v: v.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"))
    case("fragmented_full_trajectory", "REFUSE", fragmented)

    def normalized_drift(path: Path) -> None:
        value = read_manifest(path); data = path.parent / value["data_file"]; content = list(csv.reader(data.read_text().splitlines())); content[1][7] = "999999"
        with data.open("w", newline="") as stream: csv.writer(stream, lineterminator="\n").writerows(content)
        value["data_sha256"] = digest(data); write_manifest(path, value)
    case("normalized_value_without_source_derivation", "REFUSE", normalized_drift)
    case("test_adapter_cannot_self_declare_actual", "REFUSE", edit(lambda v: v.__setitem__("evidence_class", "ACTUAL_SURFACE_MEASUREMENT")))
    case("self_declared_prospective_timing_boundary", "ACCEPT", edit(lambda v: (v.__setitem__("protocol_timing", "PROSPECTIVE_PRE_OUTCOME"), v.__setitem__("protocol_frozen_at_utc", "2000-01-01T00:00:00Z"))))

    def manifest_link(path: Path) -> Path:
        nominal = path.parent.parent / "nominal"; nominal.mkdir(); link = nominal / "manifest.json"; link.symlink_to(path); return link
    case("linked_manifest_file", "REFUSE", manifest_link)

    def hardlink_source(path: Path) -> None:
        source = path.parent / "raw.csv"; alias = path.parent / "raw-hard.csv"; os.link(source, alias); value = read_manifest(path)
        value["source_artifacts"].append({"path": alias.name, "sha256": digest(alias), "media_type": "text/csv"}); write_manifest(path, value)
    case("duplicate_source_hardlink", "REFUSE", hardlink_source)

    def same_data_protocol(path: Path) -> None:
        value = read_manifest(path); protocol = path.parent / value["protocol_path"]; value["data_file"] = value["protocol_path"]; value["data_sha256"] = digest(protocol); write_manifest(path, value)
    case("protocol_and_data_same_file", "REFUSE", same_data_protocol)

    def bad_protocol(path: Path) -> None:
        protocol = path.parent / "PROTOCOL.txt"; protocol.write_text("# OTHER\n"); mutate(path, lambda v: v.__setitem__("protocol_sha256", digest(protocol)))
    case("invalid_protocol_marker", "REFUSE", bad_protocol)

    def valid_full(path: Path) -> None:
        rows(path, [["j1", "sample-a", "trial", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1], ["j2", "sample-a", "trial", "WRITE", "WRITER_CONTROL", 1, 0, 1, 1], ["j3", "sample-a", "trial", "HOLD", "ENVIRONMENT", 2, 0, 1, 1], ["j4", "sample-a", "trial", "READ", "RECORD_READOUT", 3, 0, 1, 1]])
        mutate(path, lambda v: v.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"))
    case("valid_joined_full_trajectory", "ACCEPT", valid_full)
    def no_writer(path: Path) -> None: valid_full(path); data = path.parent / "observations.csv"; data.write_text(data.read_text().replace("WRITE,WRITER_CONTROL", "WRITE,ENVIRONMENT")); sync(path)
    case("joined_full_missing_writer_role", "REFUSE", no_writer)
    def no_read(path: Path) -> None: valid_full(path); data = path.parent / "observations.csv"; data.write_text(data.read_text().replace("READ,RECORD_READOUT", "READ,ENVIRONMENT")); sync(path)
    case("joined_full_missing_readout_role", "REFUSE", no_read)
    def bad_order(path: Path) -> None: valid_full(path); data = path.parent / "observations.csv"; data.write_text(data.read_text().replace("HOLD,ENVIRONMENT,2", "HOLD,ENVIRONMENT,4")); sync(path)
    case("joined_full_time_order_violation", "REFUSE", bad_order)

    def incomplete_second(path: Path) -> None:
        valid_full(path); value = read_manifest(path); value["sample_ids"] = ["sample-a", "sample-b"]; data = path.parent / "observations.csv"
        with data.open("a", newline="") as stream: csv.writer(stream, lineterminator="\n").writerow(["b1", "sample-b", "trial-b", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1])
        write_manifest(path, value); sync(path)
    case("full_second_sample_incomplete", "REFUSE", incomplete_second)
    case("normalization_adapter_hash_drift", "REFUSE", edit(lambda v: v["normalization"].__setitem__("adapter_sha256", "0" * 64)))
    case("unregistered_normalization_adapter", "REFUSE", edit(lambda v: v["normalization"].__setitem__("adapter_id", "UNREGISTERED")))
    case("normalization_source_not_in_artifacts", "REFUSE", edit(lambda v: v["normalization"].__setitem__("source_paths", ["missing.csv"])))

    def source_data_alias(path: Path) -> None:
        value = read_manifest(path); data = path.parent / value["data_file"]; alias = path.parent / "source-data.csv"; os.link(data, alias)
        value["source_artifacts"] = [{"path": alias.name, "sha256": digest(alias), "media_type": "text/csv"}]; value["normalization"]["source_paths"] = [alias.name]; write_manifest(path, value)
    case("source_and_data_hardlink_identity", "REFUSE", source_data_alias)

    def parent_link(path: Path) -> Path:
        alias = path.parent.parent / "linked-bundle"; alias.symlink_to(path.parent, target_is_directory=True); return alias / path.name
    case("linked_parent_directory", "REFUSE", parent_link)

    def raw_change(path: Path, kind: str) -> None:
        value = read_manifest(path); relative = value["normalization"]["source_paths"][0]; source = path.parent / relative; lines = source.read_bytes().splitlines()
        if kind == "quote": lines[-1] += b'"unterminated'
        elif kind == "cell": lines[-1] += b"EXTRA"
        else: source.write_bytes(source.read_bytes().replace(b"Moment (m) Status,\r\n", b"Moment (m) Status,Unexpected Column\r\n", 1)); refresh_actual_source(path, relative); return
        source.write_bytes(b"\r\n".join(lines) + b"\r\n"); refresh_actual_source(path, relative)
    case("malformed_raw_lakeshore_csv", "REFUSE", lambda p: raw_change(p, "quote"), actual=True)
    case("declared_sample_without_rows", "REFUSE", edit(lambda v: v.__setitem__("sample_ids", ["sample-a", "ghost"])))

    def missing_rate(path: Path, observed: bool) -> None:
        content = [["a", "sample-a", "rate", "RATE", "RECORD_READOUT", 0, 0, 1, 1]]
        if observed: content.append(["b", "sample-b", "post", "POST_WRITE", "RECORD_READOUT", 0, 0, 1, 1])
        rows(path, content); mutate(path, lambda v: (v.__setitem__("sample_ids", ["sample-a", "sample-b"]), v.__setitem__("observation_scope", "RATE_STABILITY_ONLY")))
    case("rate_scope_missing_sample", "REFUSE", lambda p: missing_rate(p, False))
    case("rate_scope_observed_sample_lacks_rate", "REFUSE", lambda p: missing_rate(p, True))
    case("raw_lakeshore_nonempty_trailing_cell", "REFUSE", lambda p: raw_change(p, "cell"), actual=True)
    case("raw_lakeshore_table_header_drift", "REFUSE", lambda p: raw_change(p, "header"), actual=True)

    prior = [item for item in results if item["group"] == "prior54"]
    failed = [item for item in results if not item["expectation_met"]]
    return {
        "case_count": len(results), "prior54_count": len(prior),
        "prior54_all_pass": len(prior) == 54 and all(item["expectation_met"] for item in prior),
        "three_prior_survivors_now_refuse": all(next(item for item in results if item["name"] == name)["observed"] == "REFUSE" for name in ("rate_scope_observed_sample_lacks_rate", "raw_lakeshore_nonempty_trailing_cell", "raw_lakeshore_table_header_drift")),
        "first_failed_prior_case": failed[0]["name"] if failed else None,
        "unexpected_cases": [item["name"] for item in failed], "cases": results, "pass": not failed,
    }


def d24(result: dict[str, Any]) -> None:
    cases = result["cases"]; lines = [
        "# D24 — T-53A input-contract closure audit", "",
        "Posture: **REFUTED by default**. This is local scientific provenance/schema QA;",
        "it is not a scientific verdict and authorizes no record-formation or universal claim.", "",
        "## Results", "",
        f"- Actual 760-row bundle pass: `{str(result['actual_bundle_pass']).lower()}`.",
        f"- Prior 54 expectations all pass: `{str(cases['prior54_all_pass']).lower()}`.",
        f"- Three prior surviving cases now refuse: `{str(cases['three_prior_survivors_now_refuse']).lower()}`.",
        f"- First failed prior case: `{cases['first_failed_prior_case'] or 'NONE'}`.", "",
    ]
    if cases["unexpected_cases"]:
        lines += ["## Failed fixed expectations", ""]
        by_name = {item["name"]: item for item in cases["cases"]}
        for name in cases["unexpected_cases"]:
            item = by_name[name]; lines += [f"### {name}", "", f"- Expected `{item['expected']}`; observed `{item['observed']}`.", ""]
    else:
        lines += ["## Disposition", "", "All fixed prior expectations passed. The frozen contract is `NOT_REFUTED` at this V1 closure boundary only.", ""]
    (HERE / "D24.md").write_text("\n".join(lines))


def main() -> int:
    before = pin_check(); off = official() if before["pass"] else {"pass": False}; actual = actual_bundle() if before["pass"] else {"pass": False}; cases = replay() if before["pass"] else {"pass": False}; after = pin_check()
    actual_pass = all(item.get("pass", False) for item in (before, off, actual, after)); input_pass = cases.get("pass", False)
    verdict = "NOT_REFUTED" if actual_pass and input_pass else "REFUTED"
    result = {"schema": "WAC_T53A_INGEST_CLOSURE_V001", "default": "REFUTED", "verdict": verdict, "verifier_sha256": digest(Path(__file__).resolve()), "actual_bundle_pass": actual_pass, "input_closure_pass": input_pass, "first_failed_prior_case": cases.get("first_failed_prior_case"), "scientific_claim_scored": False, "record_formation_proof_authorized": False, "universal_claim_authorized": False, "public_urm_registration_authorized": False, "initial_pins": before, "official": off, "actual": actual, "cases": cases, "final_pins": after}
    (HERE / "RESULT.json").write_text(json_text(result))
    summary = [f"T53A INGEST CLOSURE: {verdict}", f"actual_bundle_pass={str(actual_pass).lower()}", f"prior54_all_pass={str(cases.get('prior54_all_pass', False)).lower()}", f"three_prior_survivors_now_refuse={str(cases.get('three_prior_survivors_now_refuse', False)).lower()}", f"input_closure_pass={str(input_pass).lower()}", f"first_failed_prior_case={cases.get('first_failed_prior_case') or 'NONE'}", f"unexpected_cases={','.join(cases.get('unexpected_cases', [])) or 'NONE'}", "scientific_claim_scored=false", "record_formation_proof_authorized=false", "universal_claim_authorized=false", "public_urm_registration_authorized=false"]
    (HERE / "RESULT.txt").write_text("\n".join(summary) + "\n"); d24(result); print("\n".join(summary)); return 0 if verdict == "NOT_REFUTED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
