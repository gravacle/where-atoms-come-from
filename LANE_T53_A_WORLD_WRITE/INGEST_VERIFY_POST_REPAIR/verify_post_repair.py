#!/usr/bin/env python3
"""Default-refuted post-repair verifier for the world-observation ingest contract.

All generated fixtures live below this verifier directory.  The script does not import or
execute the sealed pre-repair verifier.
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
    "model/world_observation.py": "5663f739bf1a61e4fd3399c06e08921aabaf74283bdcf5685f4d7f92c1fa50ab",
    "model/run_world_observation.py": "7a10ad8ee91b5283c00f42492eac7dada2df2a434610f2f467abfe9b1c0fb922",
    "model/checks_world_observation.py": "0324825cffd2f79071b98acc071a1e10b6dfa6aa8d57721d707dbe087baafcf9",
    "model/lakeshore_vsm.py": "6079260bd2ba6eb947370e855fea01c3a529c60c1791e663fab01e3f6527549b",
    "LANE_T53_A_WORLD_WRITE/OBSERVATION_PROTOCOL.md": "1d32ac55feced8e15d1183e8652dfe7f7d74d19aa1ac0b68fc8b288ae12c6b70",
    "LANE_T53_A_WORLD_WRITE/export_world_observation.py": "74f2952cedeef3c01b776ca3fa3632e56941d85513f585218f76fb12108bfc54",
    "LANE_T53_A_WORLD_WRITE/raw/SOURCE.json": "76d214c1e0cad7a3416d039fa92d0872565f736f885df5b9014678f9719a3e90",
    "LANE_T53_A_WORLD_WRITE/world_observation.csv": "82fca51d5f7923107763b01340a85842d5a06d770ecc88e8e518a696a2e3a891",
    "LANE_T53_A_WORLD_WRITE/world_observation.json": "79e4500cdd7cea4fee35c5f78859fb2988246017a6e567d7c2fdbb4f84219eaa",
    "LANE_T53_A_WORLD_WRITE/world_observation_certificate.json": "d4565bc8f27a5e4af798575b4810b02f0dba67387fb0b9a2ee0088b72819f597",
}

CORE = [
    "row_id", "record_id", "event_id", "stage", "role", "time", "coordinate", "value"
]
PAIR_TOKEN = "-hys-dcd-forc1 "


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def file_identity(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino


def stable_json(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True) + "\n"


def freeze_check() -> dict[str, Any]:
    actual = {name: digest(REPO / name) for name in FROZEN_HASHES}
    return {"expected": FROZEN_HASHES, "actual": actual, "pass": actual == FROZEN_HASHES}


def import_targets():
    sys.path.insert(0, str(MODEL))
    try:
        import world_observation  # type: ignore
        import lakeshore_vsm  # type: ignore
    finally:
        sys.path.pop(0)
    return world_observation, lakeshore_vsm


CONTRACT, REGISTERED_ADAPTER = import_targets()


def load_target(path: Path) -> tuple[bool, str, dict[str, Any] | None]:
    try:
        observation = CONTRACT.load_world_observation(path)
    except CONTRACT.ObservationRefusal as exc:
        return False, str(exc), None
    rendered = CONTRACT.certificate_json(observation)
    return True, rendered, json.loads(rendered)


def official_and_actual() -> dict[str, Any]:
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    checks = []
    for name in ("checks_world_observation.py", "checks_lakeshore_vsm.py"):
        result = subprocess.run(
            [sys.executable, "-B", str(MODEL / name)], cwd=REPO, env=environment,
            capture_output=True, check=False,
        )
        checks.append(
            {
                "name": name,
                "returncode": result.returncode,
                "stdout": result.stdout.decode("utf-8", errors="replace").strip(),
                "stderr": result.stderr.decode("utf-8", errors="replace").strip(),
            }
        )
    command = [
        sys.executable, "-B", str(MODEL / "run_world_observation.py"),
        str(LANE / "world_observation.json"),
    ]
    runs = [
        subprocess.run(
            command, cwd=REPO, env=environment, capture_output=True, check=False
        )
        for _ in range(3)
    ]
    outputs = [run.stdout for run in runs]
    stored = (LANE / "world_observation_certificate.json").read_bytes()
    certificate = json.loads(outputs[0]) if runs[0].returncode == 0 else {}
    expected_safety = {
        "row_count": 760,
        "scope_classification": "CONFIGURATION_EVIDENCE_ONLY",
        "scientific_verdict": "NONE_NOT_SCORED",
        "record_formation_proof_authorized": False,
        "universal_claim_authorized": False,
        "public_urm_registration_authorized": False,
        "independent_reproduction_attested": False,
        "custody_validation": "CONTENT_HASHED_NOT_PHYSICAL_ORIGIN_AUTHENTICATED",
        "normalized_source_join": "VERIFIED",
    }
    observed_safety = {
        "row_count": certificate.get("coverage", {}).get("row_count"),
        "scope_classification": certificate.get("scope_classification"),
        "scientific_verdict": certificate.get("scientific_verdict"),
        "record_formation_proof_authorized": certificate.get(
            "record_formation_proof_authorized"
        ),
        "universal_claim_authorized": certificate.get("universal_claim_authorized"),
        "public_urm_registration_authorized": certificate.get(
            "public_urm_registration_authorized"
        ),
        "independent_reproduction_attested": certificate.get(
            "independent_reproduction_attested"
        ),
        "custody_validation": certificate.get("custody_validation"),
        "normalized_source_join": certificate.get("normalization", {}).get(
            "normalized_source_join"
        ),
    }
    expected_check_lines = {
        "checks_world_observation.py": "WORLD_OBSERVATION_CHECKS: 24/24 PASS",
        "checks_lakeshore_vsm.py": "LAKESHORE_VSM_CHECKS: 8/8 PASS",
    }
    checks_pass = all(
        item["returncode"] == 0 and item["stdout"] == expected_check_lines[item["name"]]
        for item in checks
    )
    return {
        "checks": checks,
        "cli_returncodes": [run.returncode for run in runs],
        "three_cli_runs_byte_identical": len(set(outputs)) == 1,
        "stored_certificate_byte_identical": outputs[0] == stored,
        "certificate_sha256": hashlib.sha256(outputs[0]).hexdigest(),
        "safety_fields": observed_safety,
        "safety_fields_pass": observed_safety == expected_safety,
        "pass": (
            checks_pass
            and all(run.returncode == 0 for run in runs)
            and len(set(outputs)) == 1
            and outputs[0] == stored
            and observed_safety == expected_safety
        ),
    }


def one_value(lines: list[str], prefix: str) -> str:
    matches = [line.split(":", 1)[1].strip() for line in lines if line.startswith(prefix + ":")]
    if len(matches) != 1:
        raise AssertionError(f"expected one {prefix!r}, found {len(matches)}")
    return matches[0]


def parse_trace_independently(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    lines = payload.decode("latin-1").splitlines()
    tables = [index for index, line in enumerate(lines) if line.startswith("Step,Iteration,")]
    if len(tables) != 1:
        raise AssertionError(f"{path.name}: data table count {len(tables)}")
    fields: list[float] = []
    moments: list[float] = []
    times: list[float] = []
    reader = csv.reader(lines[tables[0] + 1 :], strict=True)
    for number, row in enumerate(reader, start=tables[0] + 2):
        if not row or all(not cell.strip() for cell in row):
            continue
        if len(row) not in {8, 9}:
            raise AssertionError(f"{path.name}:{number}: unexpected source width {len(row)}")
        if len(row) == 9 and row[8] != "":
            raise AssertionError(f"{path.name}:{number}: nonempty trailing source cell")
        values = [float(row[index]) for index in (3, 4, 5)]
        if not all(math.isfinite(value) for value in values):
            raise AssertionError(f"{path.name}:{number}: nonfinite source value")
        if row[6].strip() != "GOOD" or row[7].strip() != "GOOD":
            raise AssertionError(f"{path.name}:{number}: bad instrument status")
        fields.append(values[0])
        moments.append(values[1])
        times.append(values[2])
    if len(fields) < 3 or any(right < left for left, right in zip(times, times[1:])):
        raise AssertionError(f"{path.name}: insufficient/nonmonotone source")
    if "#HYSTERESIS MEASUREMENT" in lines and "#REMANENCE CURVES MEASUREMENT" not in lines:
        kind = "hys"
        writer_off = None
    elif "#REMANENCE CURVES MEASUREMENT" in lines and "#HYSTERESIS MEASUREMENT" not in lines:
        kind = "dcd"
        writer_off = one_value(lines, "Measure moment at applied fields") == "False"
    else:
        raise AssertionError(f"{path.name}: ambiguous measurement kind")
    start = datetime.strptime(one_value(lines, "START TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    finish = datetime.strptime(one_value(lines, "FINISH TIME"), "%m/%d/%Y %I:%M:%S.%f %p")
    if finish < start:
        raise AssertionError(f"{path.name}: negative duration")
    return {
        "path": path,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "kind": kind,
        "writer_off": writer_off,
        "fields": fields,
        "moments": moments,
        "times": times,
        "magnet": one_value(lines, "Magnet"),
        "sample": one_value(lines, "ID"),
        "start": start,
        "finish": finish,
    }


def strict_interpolate(x: list[float], y: list[float]) -> tuple[float, float, int]:
    brackets = [index for index in range(len(x) - 1) if x[index] * x[index + 1] < 0]
    if len(brackets) != 1:
        raise AssertionError(f"strict zero bracket count {len(brackets)}")
    index = brackets[0]
    fraction = -x[index] / (x[index + 1] - x[index])
    value = y[index] + fraction * (y[index + 1] - y[index])
    return value, fraction, index


def independent_canonical(paths: list[Path]) -> bytes:
    traces = {path.name: parse_trace_independently(path) for path in paths}
    grouped: dict[str, dict[str, dict[str, Any]]] = {}
    for name, trace in traces.items():
        sample, remainder = name.split(PAIR_TOKEN, 1)
        kind = (
            "hys" if remainder == "Step 1 Hysteresis Measurement.csv"
            else "dcd" if remainder == "Step 2 Remanence Curves.csv"
            else "unknown"
        )
        if kind == "unknown" or trace["kind"] != kind:
            raise AssertionError(f"filename/header kind mismatch: {name}")
        grouped.setdefault(sample, {})[kind] = trace
    if len(grouped) != 5 or any(set(pair) != {"hys", "dcd"} for pair in grouped.values()):
        raise AssertionError("expected five complete source pairs")
    rows: list[list[Any]] = []
    for sample in sorted(grouped):
        hys, dcd = grouped[sample]["hys"], grouped[sample]["dcd"]
        if hys["magnet"] != dcd["magnet"] or hys["magnet"] != "EM7-CSB":
            raise AssertionError(f"cross-instrument pair: {sample}")
        if not dcd["writer_off"]:
            raise AssertionError(f"non-writer-off DCD: {sample}")
        turn = min(range(len(hys["fields"])), key=hys["fields"].__getitem__)
        branches = [
            ("positive-history", slice(0, turn + 1), max(hys["fields"])),
            ("negative-history", slice(turn, len(hys["fields"])), min(hys["fields"])),
        ]
        ordinal = 0
        for label, branch, writer_field in branches:
            fields = hys["fields"][branch]
            moments = hys["moments"][branch]
            times = hys["times"][branch]
            moment, fraction, crossing = strict_interpolate(fields, moments)
            source_time = times[crossing] + fraction * (times[crossing + 1] - times[crossing])
            rows.append([
                f"{sample}:hys:{label}", sample, f"hys:{label}", "POST_WRITE",
                "RECORD_READOUT", ordinal, 0.0, moment, writer_field, source_time,
                "HYSTERESIS_ZERO_FIELD_INTERPOLANT",
            ])
            ordinal += 1
        for index, (field, moment, source_time) in enumerate(
            zip(dcd["fields"], dcd["moments"], dcd["times"])
        ):
            rows.append([
                f"{sample}:dcd:{index:03d}", sample, f"dcd:{index:03d}",
                "POST_WRITE", "RECORD_READOUT", ordinal, 0.0, moment, field,
                source_time, "DCD_WRITER_OFF_READOUT",
            ])
            ordinal += 1
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow([
        *CORE, "writer_field_T", "source_time_s", "protocol_segment"
    ])
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def actual_adapter_and_identity() -> dict[str, Any]:
    manifest_path = LANE / "world_observation.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    declared_sources = manifest["normalization"]["source_paths"]
    expected_sources = sorted(
        path.relative_to(LANE).as_posix() for path in RAW.glob("*.csv")
    )
    paths = [LANE / relative for relative in declared_sources]
    registered = REGISTERED_ADAPTER.normalized_world_observation_csv(paths)
    independent = independent_canonical(paths)
    actual = (LANE / manifest["data_file"]).read_bytes()
    artifact_paths = [LANE / item["path"] for item in manifest["source_artifacts"]]
    role_paths = [manifest_path, LANE / manifest["protocol_path"], LANE / manifest["data_file"]]
    identities = [file_identity(path.resolve()) for path in role_paths + artifact_paths]
    sources_in_artifacts = all(
        relative in {item["path"] for item in manifest["source_artifacts"]}
        for relative in declared_sources
    )
    source_registry = json.loads((RAW / "SOURCE.json").read_text(encoding="utf-8"))
    registry_hashes = {item["file"]: item["sha256"] for item in source_registry["selected_files"]}
    registry_pass = all(
        digest(path) == registry_hashes[path.name]
        for path in paths
    ) and set(registry_hashes) == {path.name for path in paths}
    return {
        "declared_source_count": len(declared_sources),
        "declared_sources_exact_ten": declared_sources == expected_sources,
        "declared_sources_all_hashed_artifacts": sources_in_artifacts,
        "registered_adapter_exact_bytes": registered == actual,
        "independent_adapter_exact_bytes": independent == actual,
        "registered_equals_independent": registered == independent,
        "normalized_sha256": hashlib.sha256(registered).hexdigest(),
        "source_registry_hashes_pass": registry_pass,
        "role_and_artifact_file_count": len(identities),
        "all_file_role_identities_disjoint": len(set(identities)) == len(identities),
        "adapter_hash_bound": (
            manifest["normalization"]["adapter_sha256"]
            == digest(MODEL / "lakeshore_vsm.py")
        ),
        "pass": all(
            (
                declared_sources == expected_sources,
                sources_in_artifacts,
                registered == actual,
                independent == actual,
                registered == independent,
                registry_pass,
                len(set(identities)) == len(identities),
                manifest["normalization"]["adapter_sha256"]
                == digest(MODEL / "lakeshore_vsm.py"),
            )
        ),
    }


def write_synthetic_bundle(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    protocol = root / "PROTOCOL.txt"
    data = root / "observations.csv"
    source = root / "raw.csv"
    protocol.write_text(
        "# WAC WORLD OBSERVATION PROTOCOL V001\n\npost-repair adversarial fixture\n",
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
        "run_id": "synthetic-run",
        "instrument_id": "synthetic-instrument",
        "source_uri": "urn:wac:post-repair-fixture",
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
        "normalization": {
            "adapter_id": "DIRECT_CSV_COPY_TEST_ONLY_V001",
            "adapter_sha256": digest(MODEL / "world_observation.py"),
            "source_paths": [source.name],
        },
        "source_artifacts": [
            {"path": source.name, "sha256": digest(source), "media_type": "text/csv"}
        ],
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(stable_json(manifest), encoding="utf-8")
    return manifest_path


def read_manifest(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_manifest(path: Path, value: dict[str, Any]) -> None:
    path.write_text(stable_json(value), encoding="utf-8")


def mutate_manifest(path: Path, mutation: Callable[[dict[str, Any]], None]) -> None:
    value = read_manifest(path)
    mutation(value)
    write_manifest(path, value)


def sync_direct(path: Path) -> None:
    value = read_manifest(path)
    data = path.parent / value["data_file"]
    source_name = value["normalization"]["source_paths"][0]
    source = path.parent / source_name
    source.write_bytes(data.read_bytes())
    value["data_sha256"] = digest(data)
    for artifact in value["source_artifacts"]:
        if artifact["path"] == source_name:
            artifact["sha256"] = digest(source)
    write_manifest(path, value)


def replace_synthetic_data(path: Path, rows: list[list[Any]]) -> None:
    value = read_manifest(path)
    data = path.parent / value["data_file"]
    with data.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(value["columns"])
        writer.writerows(rows)
    sync_direct(path)


def clone_actual(root: Path) -> Path:
    root.mkdir(parents=True)
    (root / "raw").mkdir()
    for name in (
        "OBSERVATION_PROTOCOL.md", "export_world_observation.py",
        "world_observation.csv", "world_observation.json",
    ):
        shutil.copy2(LANE / name, root / name)
    for path in RAW.iterdir():
        if path.is_file():
            shutil.copy2(path, root / "raw" / path.name)
    return root / "world_observation.json"


def mutation_suite() -> dict[str, Any]:
    scratch = VERIFY / "_scratch"
    scratch.mkdir(exist_ok=True)
    cases: list[dict[str, Any]] = []

    def execute(
        name: str,
        family: str,
        expected: str,
        mutation: Callable[[Path], Path | None],
        actual_clone: bool = False,
    ) -> None:
        with tempfile.TemporaryDirectory(prefix=name + "-", dir=scratch) as temporary:
            temp = Path(temporary)
            manifest = clone_actual(temp / "bundle") if actual_clone else write_synthetic_bundle(temp / "bundle")
            returned = mutation(manifest)
            target = returned if returned is not None else manifest
            accepted, detail, certificate = load_target(target)
            observed = "ACCEPT" if accepted else "REFUSE"
            cases.append(
                {
                    "name": name,
                    "family": family,
                    "expected": expected,
                    "observed": observed,
                    "expectation_met": observed == expected,
                    "detail": detail[:700],
                    "scope_classification": (
                        certificate.get("scope_classification") if certificate else None
                    ),
                }
            )

    execute("valid_synthetic", "control", "ACCEPT", lambda path: None)

    # The eight substitutions accepted by the sealed pre-repair verifier.
    def duplicate_outer(path: Path) -> None:
        text = path.read_text(encoding="utf-8")
        path.write_text(text.replace("{\n", '{\n  "schema": "ATTACKER",\n', 1), encoding="utf-8")

    execute("duplicate_outer_json_key", "prior_eight", "REFUSE", duplicate_outer)

    def duplicate_nested(path: Path) -> None:
        text = path.read_text(encoding="utf-8")
        needle = '"units": {\n'
        path.write_text(text.replace(needle, needle + '    "time": "ATTACKER",\n', 1), encoding="utf-8")

    execute("duplicate_nested_json_key", "prior_eight", "REFUSE", duplicate_nested)

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
        sync_direct(path)

    execute("malformed_csv_unclosed_quote", "prior_eight", "REFUSE", malformed_normalized)

    def fragmented_full(path: Path) -> None:
        replace_synthetic_data(
            path,
            [
                ["f1", "sample-a", "before-event", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1],
                ["f2", "sample-a", "write-event", "WRITE", "WRITER_CONTROL", 1, 0, 1, 1],
                ["f3", "sample-a", "hold-event", "HOLD", "ENVIRONMENT", 2, 0, 1, 1],
                ["f4", "sample-a", "read-event", "READ", "RECORD_READOUT", 3, 0, 1, 1],
            ],
        )
        mutate_manifest(path, lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"))

    execute("fragmented_false_full_scope_upgrade", "prior_eight", "REFUSE", fragmented_full)

    def normalized_drift(path: Path) -> None:
        value = read_manifest(path)
        data = path.parent / value["data_file"]
        rows = list(csv.reader(data.read_text(encoding="utf-8").splitlines()))
        rows[1][7] = "999999"
        with data.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream, lineterminator="\n").writerows(rows)
        value["data_sha256"] = digest(data)
        write_manifest(path, value)

    execute("normalized_value_not_joined_to_source", "prior_eight", "REFUSE", normalized_drift)

    def manifest_link(path: Path) -> Path:
        nominal = path.parent.parent / "nominal"
        nominal.mkdir()
        link = nominal / "manifest.json"
        link.symlink_to(path)
        return link

    execute("manifest_symlink_re_roots_bundle", "prior_eight", "REFUSE", manifest_link)

    def duplicate_hardlink(path: Path) -> None:
        source = path.parent / "raw.csv"
        alias = path.parent / "raw-hardlink.csv"
        os.link(source, alias)
        value = read_manifest(path)
        value["source_artifacts"].append(
            {"path": alias.name, "sha256": digest(alias), "media_type": "text/csv"}
        )
        write_manifest(path, value)

    execute("duplicate_source_hardlink_alias", "prior_eight", "REFUSE", duplicate_hardlink)

    def data_as_protocol(path: Path) -> None:
        value = read_manifest(path)
        protocol = path.parent / value["protocol_path"]
        value["data_file"] = value["protocol_path"]
        value["data_sha256"] = digest(protocol)
        write_manifest(path, value)

    execute("protocol_data_same_identity", "prior_eight", "REFUSE", data_as_protocol)

    # Positive and negative FULL trajectory checks.
    def valid_full(path: Path) -> None:
        replace_synthetic_data(
            path,
            [
                ["j1", "sample-a", "trial-1", "BEFORE", "RECORD_READOUT", 0, 0, 0, 1],
                ["j2", "sample-a", "trial-1", "WRITE", "WRITER_CONTROL", 1, 0, 1, 1],
                ["j3", "sample-a", "trial-1", "HOLD", "ENVIRONMENT", 2, 0, 1, 1],
                ["j4", "sample-a", "trial-1", "READ", "RECORD_READOUT", 3, 0, 1, 1],
            ],
        )
        mutate_manifest(path, lambda value: value.__setitem__("observation_scope", "FULL_WRITE_HOLD_READ"))

    execute("valid_joined_full_trajectory", "trajectory", "ACCEPT", valid_full)

    def joined_missing_writer(path: Path) -> None:
        valid_full(path)
        data = path.parent / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("WRITE,WRITER_CONTROL", "WRITE,ENVIRONMENT"), encoding="utf-8")
        sync_direct(path)

    execute("joined_full_missing_writer_role", "trajectory", "REFUSE", joined_missing_writer)

    def joined_missing_read(path: Path) -> None:
        valid_full(path)
        data = path.parent / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("READ,RECORD_READOUT", "READ,ENVIRONMENT"), encoding="utf-8")
        sync_direct(path)

    execute("joined_full_missing_readout_role", "trajectory", "REFUSE", joined_missing_read)

    def joined_out_of_order(path: Path) -> None:
        valid_full(path)
        data = path.parent / "observations.csv"
        data.write_text(data.read_text(encoding="utf-8").replace("HOLD,ENVIRONMENT,2", "HOLD,ENVIRONMENT,4"), encoding="utf-8")
        sync_direct(path)

    execute("joined_full_time_order_violation", "trajectory", "REFUSE", joined_out_of_order)

    def second_sample_incomplete(path: Path) -> None:
        valid_full(path)
        mutate_manifest(path, lambda value: value.__setitem__("sample_ids", ["sample-a", "sample-b"]))

    execute("full_second_sample_missing", "trajectory", "REFUSE", second_sample_incomplete)

    # Adapter/normalization registry attacks.
    execute(
        "normalization_adapter_hash_drift", "normalization", "REFUSE",
        lambda path: (mutate_manifest(path, lambda value: value["normalization"].__setitem__("adapter_sha256", "0" * 64)), None)[1],
    )
    execute(
        "unregistered_normalization_adapter", "normalization", "REFUSE",
        lambda path: (mutate_manifest(path, lambda value: value["normalization"].__setitem__("adapter_id", "SHELL_COMMAND_V001")), None)[1],
    )
    execute(
        "normalization_source_absent_from_artifacts", "normalization", "REFUSE",
        lambda path: (mutate_manifest(path, lambda value: value["normalization"].__setitem__("source_paths", ["missing.csv"])), None)[1],
    )
    execute(
        "test_adapter_claims_actual", "normalization", "REFUSE",
        lambda path: (mutate_manifest(path, lambda value: value.__setitem__("evidence_class", "ACTUAL_SURFACE_MEASUREMENT")), None)[1],
    )

    def source_data_hardlink(path: Path) -> None:
        value = read_manifest(path)
        data = path.parent / value["data_file"]
        source = path.parent / "data-source-alias.csv"
        os.link(data, source)
        value["source_artifacts"] = [
            {"path": source.name, "sha256": digest(source), "media_type": "text/csv"}
        ]
        value["normalization"]["source_paths"] = [source.name]
        write_manifest(path, value)

    execute("source_aliases_data_identity", "identity", "REFUSE", source_data_hardlink)

    def protocol_data_hardlink(path: Path) -> None:
        value = read_manifest(path)
        protocol = path.parent / value["protocol_path"]
        alias = path.parent / "data-protocol-hardlink.txt"
        os.link(protocol, alias)
        value["data_file"] = alias.name
        value["data_sha256"] = digest(alias)
        write_manifest(path, value)

    execute("data_hardlinks_protocol_identity", "identity", "REFUSE", protocol_data_hardlink)

    # Independent root-boundary variant: final manifest is regular, lexical parent is a symlink.
    def parent_directory_link(path: Path) -> Path:
        nominal = path.parent.parent / "nominal-bundle"
        nominal.symlink_to(path.parent, target_is_directory=True)
        return nominal / path.name

    execute("manifest_parent_directory_symlink_reroot", "independent", "REFUSE", parent_directory_link)

    # Independent registered-adapter mutation: malformed final raw row in one of ten files.
    def malformed_lakeshore_source(path: Path) -> None:
        value = read_manifest(path)
        relative = value["normalization"]["source_paths"][0]
        raw_path = path.parent / relative
        lines = raw_path.read_bytes().splitlines()
        lines[-1] = lines[-1] + b'"unterminated'
        raw_path.write_bytes(b"\r\n".join(lines) + b"\r\n")
        new_digest = digest(raw_path)
        for artifact in value["source_artifacts"]:
            if artifact["path"] == relative:
                artifact["sha256"] = new_digest
        # Keep the included registry internally self-consistent on SHA/size/MD5.
        registry_path = path.parent / "raw" / "SOURCE.json"
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
        for item in registry["selected_files"]:
            if item["file"] == raw_path.name:
                item["sha256"] = new_digest
                item["bytes"] = raw_path.stat().st_size
                item["md5_repository"] = hashlib.md5(raw_path.read_bytes()).hexdigest()
        registry_path.write_text(stable_json(registry), encoding="utf-8")
        for artifact in value["source_artifacts"]:
            if artifact["path"] == "raw/SOURCE.json":
                artifact["sha256"] = digest(registry_path)
        write_manifest(path, value)

    execute(
        "malformed_lakeshore_raw_unclosed_quote",
        "independent",
        "REFUSE",
        malformed_lakeshore_source,
        actual_clone=True,
    )

    # Independent coverage mutation: a declared sample with no rows must not disappear silently.
    execute(
        "declared_sample_without_observation_rows",
        "independent",
        "REFUSE",
        lambda path: (mutate_manifest(path, lambda value: value.__setitem__("sample_ids", ["sample-a", "ghost-sample"])), None)[1],
    )

    # Independent RATE scope fragmentation across two declared samples.
    def fragmented_rate(path: Path) -> None:
        replace_synthetic_data(
            path,
            [["rate-a", "sample-a", "rate-event", "RATE", "RECORD_READOUT", 0, 0, 1, 1]],
        )
        mutate_manifest(
            path,
            lambda value: (
                value.__setitem__("sample_ids", ["sample-a", "sample-b"]),
                value.__setitem__("observation_scope", "RATE_STABILITY_ONLY"),
            ),
        )

    execute("rate_scope_missing_declared_sample", "independent", "REFUSE", fragmented_rate)

    failures = [case for case in cases if not case["expectation_met"]]
    prior = [case for case in cases if case["family"] == "prior_eight"]
    trajectory = [case for case in cases if case["family"] == "trajectory"]
    return {
        "case_count": len(cases),
        "prior_eight_count": len(prior),
        "prior_eight_all_repaired": all(case["expectation_met"] for case in prior),
        "full_trajectory_contract_pass": all(case["expectation_met"] for case in trajectory),
        "unexpected_case_count": len(failures),
        "first_surviving_attack": failures[0]["name"] if failures else None,
        "unexpected_cases": [case["name"] for case in failures],
        "cases": cases,
        "pass": not failures,
    }


def main() -> int:
    initial = freeze_check()
    official = official_and_actual() if initial["pass"] else {"pass": False, "skipped": "freeze drift"}
    adapter = actual_adapter_and_identity() if initial["pass"] else {"pass": False, "skipped": "freeze drift"}
    mutations = mutation_suite() if initial["pass"] else {"pass": False, "skipped": "freeze drift"}
    final = freeze_check()
    actual_pass = all(item.get("pass", False) for item in (initial, official, adapter, final))
    attack_pass = mutations.get("pass", False)
    verdict = "NOT_REFUTED" if actual_pass and attack_pass else "REFUTED"
    result = {
        "schema": "WAC_T53A_WORLD_INGEST_POST_REPAIR_VERIFY_V001",
        "default": "REFUTED",
        "verdict": verdict,
        "actual_760_row_bundle_pass": actual_pass,
        "generic_contract_adversarial_pass": attack_pass,
        "first_surviving_attack": mutations.get("first_surviving_attack"),
        "scientific_claim_scored": False,
        "record_formation_proof_authorized": False,
        "universal_claim_authorized": False,
        "public_urm_registration_authorized": False,
        "initial_freeze": initial,
        "official_and_actual": official,
        "actual_adapter_and_identity": adapter,
        "mutations": mutations,
        "final_freeze": final,
    }
    (VERIFY / "RESULT.json").write_text(stable_json(result), encoding="utf-8")
    lines = [
        f"T53A WORLD INGEST POST-REPAIR VERIFIER: {verdict}",
        f"actual_760_row_bundle_pass={str(actual_pass).lower()}",
        f"prior_eight_all_repaired={str(mutations.get('prior_eight_all_repaired', False)).lower()}",
        f"full_trajectory_contract_pass={str(mutations.get('full_trajectory_contract_pass', False)).lower()}",
        f"generic_contract_adversarial_pass={str(attack_pass).lower()}",
        f"first_surviving_attack={mutations.get('first_surviving_attack') or 'NONE'}",
        f"unexpected_cases={','.join(mutations.get('unexpected_cases', [])) or 'NONE'}",
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
