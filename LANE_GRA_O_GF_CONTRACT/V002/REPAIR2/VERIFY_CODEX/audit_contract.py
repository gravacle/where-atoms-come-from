#!/usr/bin/env python3
"""Independent bounded QA of V002 Repair2.

The script writes nothing.  It verifies the three subject manifests, runs the fixed
54-case suite twice, checks the typed formation-observation joins, reproduces the
Repair1 scalar escape and Repair2 refusal, and then tests whether the Repair2
compatibility projection preserves an invalid original formation/source parent
join.  The last probe is the first scientific-state inconsistency found.
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Callable


HERE = Path(__file__).resolve().parent
REPAIR2_DIR = HERE.parent
V002_DIR = REPAIR2_DIR.parent
REPAIR1_DIR = V002_DIR / "REPAIR1"
LANE_DIR = V002_DIR.parent
ROOT = LANE_DIR.parent
for path in (REPAIR2_DIR, REPAIR1_DIR, V002_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import mutation_suite_v002 as sealed_suite  # noqa: E402
import mutation_suite_v002_repair2 as suite2  # noqa: E402
import synthetic_positive_fixture_v002 as fixture  # noqa: E402
import validator_v002 as base  # noqa: E402
import validator_v002_repair1 as repair1  # noqa: E402
import validator_v002_repair2 as repair2  # noqa: E402


PROMOTED = ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_manifest(path: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for number, raw in enumerate(path.read_text().splitlines(), 1):
        if not raw.strip():
            continue
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", raw)
        if not match:
            rows.append({"line": number, "parse": False})
            continue
        expected, relative = match.groups()
        target = ROOT / relative
        actual = sha256(target) if target.is_file() else None
        rows.append({
            "line": number,
            "path": relative,
            "parse": True,
            "exists": target.is_file(),
            "expected": expected,
            "actual": actual,
            "match": expected == actual,
        })
    parsed_paths = [row.get("path") for row in rows if row.get("parse")]
    return {
        "result": "PASS" if rows and all(row.get("match") for row in rows) and len(parsed_paths) == len(set(parsed_paths)) else "FAIL",
        "matched_entries": sum(row.get("match") is True for row in rows),
        "total_entries": len(rows),
        "unique_paths": len(parsed_paths) == len(set(parsed_paths)),
        "rows": rows,
    }


def seven_no_proof(result: dict[str, Any]) -> bool:
    outputs = result.get("authoritative_proof_outputs", {})
    return set(outputs) == set(base.CLAIMS) and all(value == "NO_PROOF_OUTPUT" for value in outputs.values())


def promoted(result: dict[str, Any]) -> bool:
    outputs = result.get("authoritative_proof_outputs", {})
    return all(outputs.get(claim) == "PASSES_DECLARED_DOMAIN" for claim in PROMOTED)


def refused_with(result: dict[str, Any], diagnostic: str) -> bool:
    return (
        result.get("accepted") is False
        and result.get("actual_platform_present") is False
        and diagnostic in json.dumps(result, sort_keys=True)
        and seven_no_proof(result)
    )


def observation_candidate(mutator: Callable[[list[dict[str, Any]]], None]) -> dict[str, Any]:
    candidate = suite2.typed_platform_candidate()
    observations = copy.deepcopy(suite2.typed_observations())
    mutator(observations)
    suite2.change_observations(candidate, observations)
    return candidate


def observation_join_audit() -> dict[str, Any]:
    positive = repair2.evaluate_instance(suite2.typed_platform_candidate(), ROOT)
    rows: list[dict[str, Any]] = []

    def run(name: str, expected: str, mutator: Callable[[list[dict[str, Any]]], None]) -> None:
        result = repair2.evaluate_instance(observation_candidate(mutator), ROOT)
        rows.append({"case": name, "expected": expected, "pass": refused_with(result, expected)})

    run("event-stability", "FORMATION_LIFECYCLE_CROSS_EVENT", lambda rows: rows[1].__setitem__("event_id", "OTHER-EVENT"))
    run("surface-membership", "FORMATION_OBSERVATION_SURFACE_UNKNOWN", lambda rows: rows[1].__setitem__("surface_id", "UNKNOWN-SURFACE"))
    run("surface-stability", "FORMATION_LIFECYCLE_CROSS_SURFACE", lambda rows: rows[1].__setitem__("surface_id", "SURFACE-1"))
    run("stage-role", "FORMATION_OBSERVATION_STAGE_ROLE_INVALID", lambda rows: rows[1].__setitem__("role", "FORMATION"))
    run("timezone-aware-time", "FORMATION_OBSERVATION_TIME_INVALID", lambda rows: rows[1].__setitem__("time", "2026-08-22T09:10:00"))
    run("strict-time-order", "FORMATION_LIFECYCLE_TIME_ORDER_INVALID", lambda rows: rows[1].__setitem__("time", rows[0]["time"]))
    run("certified-value", "FORMATION_OBSERVATION_VALUE_NOT_CERTIFIED", lambda rows: rows[1].__setitem__("value", False))
    run("value-type", "FORMATION_OBSERVATION_VALUE_NONFINITE_OR_INVALID", lambda rows: rows[1].__setitem__("value", "true"))
    run("unit-stability", "FORMATION_LIFECYCLE_UNIT_UNSTABLE", lambda rows: rows[1].__setitem__("unit", "other-unit"))
    run("source-membership", "FORMATION_OBSERVATION_SOURCE_JOIN_INVALID", lambda rows: rows[1].__setitem__("source_artifact_id", "UNJOINED.SOURCE"))
    run("source-stability", "FORMATION_LIFECYCLE_SOURCE_UNSTABLE", lambda rows: rows[1].__setitem__("source_artifact_id", "UNJOINED.SOURCE"))
    run("lifecycle-coverage", "FORMATION_LIFECYCLE_COVERAGE_INVALID", lambda rows: rows.pop())
    run("lifecycle-cardinality", "FORMATION_LIFECYCLE_CARDINALITY_INVALID", lambda rows: rows.append({**rows[-1], "predicate_id": "REAL.PERSISTENCE.EXTRA", "time": "2026-08-22T09:30:00Z"}))
    run("predicate-uniqueness", "FORMATION_PREDICATE_ID_DUPLICATE", lambda rows: rows[1].__setitem__("predicate_id", rows[0]["predicate_id"]))
    run("reproducibility", "FORMATION_OBSERVATION_NOT_REPRODUCIBLE", lambda rows: rows[1].__setitem__("reproducible", False))

    return {
        "typed_positive_control": {
            "pass": (
                positive.get("accepted") is True
                and positive.get("actual_platform_present") is True
                and positive.get("formation_observation_repair2", {}).get("valid") is True
                and promoted(positive)
            ),
            "accepted": positive.get("accepted"),
            "actual_platform_present": positive.get("actual_platform_present"),
        },
        "result": "PASS" if rows and all(row["pass"] for row in rows) else "FAIL",
        "case_count": len(rows),
        "rows": rows,
    }


def scalar_regression() -> dict[str, Any]:
    candidate = suite2.typed_platform_candidate()
    suite2.change_observations(candidate, [0])
    legacy = repair1.evaluate_instance(copy.deepcopy(candidate), ROOT)
    repaired = repair2.evaluate_instance(copy.deepcopy(candidate), ROOT)
    reproduced = (
        legacy.get("accepted") is True
        and legacy.get("actual_platform_present") is True
        and promoted(legacy)
        and refused_with(repaired, "FORMATION_OBSERVATION_MEMBER_NOT_OBJECT")
    )
    return {
        "result": "PASS" if reproduced else "FAIL",
        "mutation": "REAL.PLATFORM.FORMATION.observations=[0]",
        "repair1": {
            "accepted": legacy.get("accepted"),
            "actual_platform_present": legacy.get("actual_platform_present"),
            "custody_disposition": legacy.get("custody", {}).get("disposition"),
            "promoted_claims": [claim for claim in PROMOTED if legacy.get("authoritative_proof_outputs", {}).get(claim) == "PASSES_DECLARED_DOMAIN"],
        },
        "repair2": {
            "accepted": repaired.get("accepted"),
            "actual_platform_present": repaired.get("actual_platform_present"),
            "custody_disposition": repaired.get("custody", {}).get("disposition"),
            "member_diagnostic_present": "FORMATION_OBSERVATION_MEMBER_NOT_OBJECT" in json.dumps(repaired, sort_keys=True),
            "seven_no_proof_outputs": seven_no_proof(repaired),
        },
    }


def nonpromotion_controls() -> dict[str, Any]:
    synthetic = repair2.evaluate_instance(fixture.build_fixture(), ROOT)

    mode_candidate = suite2.typed_platform_candidate()
    mode_candidate["package_mode"] = "SYNTHETIC_TEST"
    synthetic_mode = repair2.evaluate_instance(mode_candidate, ROOT)

    flag_candidate = suite2.typed_platform_candidate()
    platform = sealed_suite.payload(flag_candidate, "REAL.PLATFORM.INSTANCE")
    platform["synthetic_only"] = True
    sealed_suite.replace_payload(flag_candidate, "REAL.PLATFORM.INSTANCE", platform)
    synthetic_flag = repair2.evaluate_instance(flag_candidate, ROOT)

    incomplete = observation_candidate(lambda rows: rows.pop())
    incomplete_result = repair2.evaluate_instance(incomplete, ROOT)

    rows = [
        {
            "case": "sealed-synthetic-fixture",
            "pass": synthetic.get("actual_platform_present") is False and seven_no_proof(synthetic),
        },
        {
            "case": "platform-under-synthetic-package-mode",
            "pass": refused_with(synthetic_mode, "PLATFORM_NOT_ACTUAL_SCIENTIFIC"),
        },
        {
            "case": "platform-synthetic-only-flag",
            "pass": refused_with(synthetic_flag, "PLATFORM_NOT_ACTUAL_SCIENTIFIC"),
        },
        {
            "case": "incomplete-lifecycle",
            "pass": refused_with(incomplete_result, "FORMATION_LIFECYCLE_COVERAGE_INVALID"),
        },
    ]
    return {"result": "PASS" if all(row["pass"] for row in rows) else "FAIL", "rows": rows}


def nongeometric_boundary() -> dict[str, Any]:
    candidate = suite2.typed_platform_candidate()
    decoded_payloads: list[Any] = []
    for artifact in candidate["artifacts"]:
        if artifact.get("storage") != "INLINE_JSON":
            continue
        decoded_payloads.append(base.strict_json_loads(__import__("base64").b64decode(artifact["payload_b64"], validate=True)))

    keys: set[str] = set()

    def walk(value: Any) -> None:
        if isinstance(value, dict):
            keys.update(str(key).lower() for key in value)
            for child in value.values():
                walk(child)
        elif isinstance(value, list):
            for child in value:
                walk(child)

    walk(decoded_payloads)
    gravitational_source_keys = {
        "geometry", "spacetime_geometry", "metric_tensor", "spacetime_metric",
        "stress_energy", "stress-energy", "stress_energy_tensor",
    }
    result = repair2.evaluate_instance(candidate, ROOT)
    absent = sorted(gravitational_source_keys & keys)
    passed = not absent and all(result.get("authoritative_proof_outputs", {}).get(claim) == "PASSES_DECLARED_DOMAIN" for claim in ("GF0", "GF1"))
    return {
        "result": "PASS" if passed else "FAIL",
        "gravitational_geometry_or_source_keys_present": absent,
        "GF0": result.get("authoritative_proof_outputs", {}).get("GF0"),
        "GF1": result.get("authoritative_proof_outputs", {}).get("GF1"),
        "note": "The transport W calibration is not a spacetime/gravitational metric requirement.",
    }


def parent_join_counterexample() -> dict[str, Any]:
    candidate = suite2.typed_platform_candidate()
    formation = next(
        artifact for artifact in candidate["artifacts"]
        if artifact.get("artifact_id") == suite2.FORMATION_ID
    )
    parent = next(
        item for item in formation["parents"]
        if item.get("artifact_id") == "REAL.PLATFORM.SOURCE"
    )
    correct_sha256 = parent["sha256"]
    corrupt_sha256 = "0" * 64 if correct_sha256 != "0" * 64 else "1" * 64
    parent["sha256"] = corrupt_sha256

    before = repair1.evaluate_instance(copy.deepcopy(candidate), ROOT)
    after = repair2.evaluate_instance(copy.deepcopy(candidate), ROOT)
    mismatch = "ARTIFACT_PARENT_HASH_MISMATCH:REAL.PLATFORM.FORMATION:REAL.PLATFORM.SOURCE"
    reproduced = (
        refused_with(before, mismatch)
        and after.get("accepted") is True
        and after.get("custody", {}).get("disposition") == "QUALIFIED"
        and after.get("actual_platform_present") is True
        and after.get("formation_observation_repair2", {}).get("valid") is True
        and promoted(after)
    )
    return {
        "id": "FORMATION_SOURCE_PARENT_HASH_MISMATCH_LAUNDERED",
        "reproduced": reproduced,
        "mutation": {
            "artifact_id": suite2.FORMATION_ID,
            "parent_artifact_id": "REAL.PLATFORM.SOURCE",
            "field": "parents[].sha256",
            "correct_sha256": correct_sha256,
            "submitted_sha256": corrupt_sha256,
            "payload_or_artifact_digest_changed": False,
        },
        "repair1_original_candidate": {
            "accepted": before.get("accepted"),
            "custody_disposition": before.get("custody", {}).get("disposition"),
            "actual_platform_present": before.get("actual_platform_present"),
            "mismatch_diagnostic_present": mismatch in json.dumps(before, sort_keys=True),
            "seven_no_proof_outputs": seven_no_proof(before),
        },
        "repair2_same_submitted_candidate": {
            "accepted": after.get("accepted"),
            "custody_disposition": after.get("custody", {}).get("disposition"),
            "actual_platform_present": after.get("actual_platform_present"),
            "formation_observation_repair2_valid": after.get("formation_observation_repair2", {}).get("valid"),
            "promoted_claims": [claim for claim in PROMOTED if after.get("authoritative_proof_outputs", {}).get(claim) == "PASSES_DECLARED_DOMAIN"],
            "UGE": after.get("authoritative_proof_outputs", {}).get("UGE"),
            "custody_errors": after.get("custody", {}).get("errors", []),
        },
        "cause": {
            "ignored_original_artifact_errors": "validate_actual_formation_observations records artifact_errors but does not fail on this parent mismatch",
            "projection_rewrite": "_project_for_repair1 rewrites every resolved parents[].sha256 to the registry digest before Repair1 evaluates the copy",
            "subject_lines": "validator_v002_repair2.py:72,195-202,226-228",
        },
    }


def main() -> int:
    manifests = {
        "v002": audit_manifest(V002_DIR / "MANIFEST_V002.sha256"),
        "repair1": audit_manifest(REPAIR1_DIR / "MANIFEST_V002_REPAIR1.sha256"),
        "repair2": audit_manifest(REPAIR2_DIR / "MANIFEST_V002_REPAIR2.sha256"),
    }
    suite_run_1 = suite2.run_suite()
    suite_run_2 = suite2.run_suite()
    fixed_suite = {
        "result": "PASS" if (
            suite_run_1.get("overall_result") == "PASS"
            and suite_run_2.get("overall_result") == "PASS"
            and suite_run_1.get("case_count") == suite_run_2.get("case_count") == 54
            and suite_run_1.get("failed_count") == suite_run_2.get("failed_count") == 0
            and suite_run_1.get("result_payload_sha256") == suite_run_2.get("result_payload_sha256")
        ) else "FAIL",
        "runs": 2,
        "case_count_each_run": [suite_run_1.get("case_count"), suite_run_2.get("case_count")],
        "failed_count_each_run": [suite_run_1.get("failed_count"), suite_run_2.get("failed_count")],
        "result_payload_sha256_each_run": [suite_run_1.get("result_payload_sha256"), suite_run_2.get("result_payload_sha256")],
        "typed_lifecycle_positive_control_each_run": [suite_run_1.get("typed_lifecycle_positive_control"), suite_run_2.get("typed_lifecycle_positive_control")],
    }
    scalar = scalar_regression()
    joins = observation_join_audit()
    nonpromotion = nonpromotion_controls()
    nongeometric = nongeometric_boundary()
    inconsistency = parent_join_counterexample()
    controls_pass = (
        all(value["result"] == "PASS" for value in manifests.values())
        and fixed_suite["result"] == "PASS"
        and scalar["result"] == "PASS"
        and joins["result"] == "PASS"
        and joins["typed_positive_control"]["pass"]
        and nonpromotion["result"] == "PASS"
        and nongeometric["result"] == "PASS"
    )
    result = {
        "audit_id": "GRA-O-GF-CONTRACT-V002-REPAIR2-VERIFY-CODEX",
        "audit_date": "2026-08-22",
        "subject": "GRA-O-GF-CONTRACT-V002-REPAIR2",
        "overall_result": "REFUTED" if inconsistency["reproduced"] else "PASS" if controls_pass else "AUDIT_CONTROL_FAILURE",
        "scientific_result": "NO_PROOF_OUTPUT",
        "manifest_checks": manifests,
        "fixed_suite": fixed_suite,
        "repair1_scalar_inconsistency_and_repair2_refusal": scalar,
        "typed_observation_join_checks": joins,
        "synthetic_and_incomplete_nonpromotion_controls_before_counterexample": nonpromotion,
        "GF0_GF1_nongeometric_boundary": nongeometric,
        "first_scientific_state_inconsistency": inconsistency,
        "subject_files_modified": False,
        "audit_files_confined_to": "LANE_GRA_O_GF_CONTRACT/V002/REPAIR2/VERIFY_CODEX",
    }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if result["overall_result"] in {"PASS", "REFUTED"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
