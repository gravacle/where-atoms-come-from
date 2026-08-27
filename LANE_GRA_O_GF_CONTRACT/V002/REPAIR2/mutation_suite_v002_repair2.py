#!/usr/bin/env python3
"""Prior 47 cases plus seven exact Repair2 observation regressions."""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Callable


REPAIR2_DIR = Path(__file__).resolve().parent
REPAIR1_DIR = REPAIR2_DIR.parent / "REPAIR1"
V002_DIR = REPAIR2_DIR.parent
for path in (REPAIR1_DIR, V002_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import mutation_suite_v002 as sealed_suite  # noqa: E402
import mutation_suite_v002_repair1 as repair1_suite  # noqa: E402
import validator_v002 as base  # noqa: E402
import validator_v002_repair1 as repair1  # noqa: E402
import validator_v002_repair2 as repair2  # noqa: E402


WORKSPACE = Path(__file__).resolve().parents[3]
FORMATION_ID = "REAL.PLATFORM.FORMATION"


def typed_observations() -> list[dict[str, Any]]:
    rows = []
    for index, (stage, role) in enumerate(repair2.LIFECYCLE):
        rows.append({
            "predicate_id": f"REAL.{stage}", "event_id": "REAL-EVENT-001",
            "surface_id": "SURFACE-0", "stage": stage, "role": role,
            "time": f"2026-08-22T09:{index * 10:02d}:00Z", "value": True,
            "unit": "certified_boolean", "source_artifact_id": "REAL.PLATFORM.SOURCE",
            "reproducible": True,
        })
    return rows


def typed_platform_candidate() -> dict[str, Any]:
    candidate = repair1_suite.joined_platform_candidate()
    evidence = sealed_suite.payload(candidate, FORMATION_ID)
    evidence["observations"] = typed_observations()
    repair1_suite.replace_payload(candidate, FORMATION_ID, evidence)
    return candidate


def change_observations(candidate: dict[str, Any], observations: list[Any]) -> None:
    evidence = sealed_suite.payload(candidate, FORMATION_ID)
    evidence["observations"] = observations
    repair1_suite.replace_payload(candidate, FORMATION_ID, evidence)


def change_nonfinite(candidate: dict[str, Any]) -> None:
    evidence = sealed_suite.payload(candidate, FORMATION_ID)
    evidence["observations"][0]["value"] = float("nan")
    raw = json.dumps(evidence, sort_keys=True, separators=(",", ":"), allow_nan=True).encode()
    artifact = sealed_suite.artifact(candidate, FORMATION_ID)
    artifact["payload_b64"] = base64.b64encode(raw).decode()
    artifact["sha256"] = hashlib.sha256(raw).hexdigest()
    artifact["byte_length"] = len(raw)
    registry = {item["artifact_id"]: item for item in candidate["artifacts"]}
    for item in candidate["artifacts"]:
        for parent in item["parents"]:
            if parent["artifact_id"] in registry:
                parent["sha256"] = registry[parent["artifact_id"]]["sha256"]


def mutate_scalar(candidate: dict[str, Any]) -> None:
    change_observations(candidate, [0])


def mutate_mixed(candidate: dict[str, Any]) -> None:
    observations = typed_observations()
    observations.append(0)
    change_observations(candidate, observations)


def mutate_empty_dict(candidate: dict[str, Any]) -> None:
    observations = typed_observations()
    observations[0] = {}
    change_observations(candidate, observations)


def mutate_unknown_key(candidate: dict[str, Any]) -> None:
    observations = typed_observations()
    observations[0]["unexpected"] = "bypass"
    change_observations(candidate, observations)


def mutate_missing_identity(candidate: dict[str, Any]) -> None:
    observations = typed_observations()
    observations[0]["event_id"] = ""
    change_observations(candidate, observations)


def mutate_cross_event(candidate: dict[str, Any]) -> None:
    observations = typed_observations()
    observations[1]["event_id"] = "REAL-EVENT-OTHER"
    change_observations(candidate, observations)


Case = tuple[str, Callable[[dict[str, Any]], None], str]
CASES: list[Case] = [
    ("REPAIR2-exact-scalar-list", mutate_scalar, "FORMATION_OBSERVATION_MEMBER_NOT_OBJECT"),
    ("REPAIR2-mixed-dict-scalar", mutate_mixed, "FORMATION_OBSERVATION_MEMBER_NOT_OBJECT"),
    ("REPAIR2-empty-dict", mutate_empty_dict, "FORMATION_OBSERVATION_KEYS_INVALID"),
    ("REPAIR2-unknown-observation-key", mutate_unknown_key, "FORMATION_OBSERVATION_KEYS_INVALID"),
    ("REPAIR2-missing-identity", mutate_missing_identity, "FORMATION_OBSERVATION_IDENTITY_MISSING"),
    ("REPAIR2-nonfinite-value", change_nonfinite, "FORMATION_OBSERVATION_NONFINITE_OR_INVALID_JSON"),
    ("REPAIR2-cross-event-lifecycle", mutate_cross_event, "FORMATION_LIFECYCLE_CROSS_EVENT"),
]


def no_proof(result: dict[str, Any]) -> bool:
    outputs = result.get("authoritative_proof_outputs", {})
    return bool(outputs) and all(value == "NO_PROOF_OUTPUT" for value in outputs.values())


def run_suite() -> dict[str, Any]:
    original_repair1_evaluate = repair1.evaluate_instance
    repair1.evaluate_instance = repair2.evaluate_instance
    try:
        prior = repair1_suite.run_suite()
    finally:
        repair1.evaluate_instance = original_repair1_evaluate

    positive = typed_platform_candidate()
    positive_result = repair2.evaluate_instance(copy.deepcopy(positive), WORKSPACE)
    positive_ok = (
        positive_result.get("accepted") is True
        and positive_result.get("actual_platform_present") is True
        and positive_result.get("formation_observation_repair2", {}).get("valid") is True
        and all(positive_result.get("authoritative_proof_outputs", {}).get(claim) == "PASSES_DECLARED_DOMAIN" for claim in ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2"))
    )

    rows: list[dict[str, Any]] = []
    for name, mutation, expected in CASES:
        candidate = typed_platform_candidate()
        mutation(candidate)
        legacy = repair1.evaluate_instance(copy.deepcopy(candidate), WORKSPACE) if name == "REPAIR2-exact-scalar-list" else None
        result = repair2.evaluate_instance(candidate, WORKSPACE)
        text = json.dumps(result, sort_keys=True)
        legacy_escape = bool(
            legacy and legacy.get("accepted") is True and legacy.get("actual_platform_present") is True
            and all(legacy.get("authoritative_proof_outputs", {}).get(claim) == "PASSES_DECLARED_DOMAIN" for claim in ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2"))
        )
        passed = (
            positive_ok and result.get("accepted") is False
            and result.get("actual_platform_present") is False
            and expected in text and no_proof(result)
            and (legacy_escape if name == "REPAIR2-exact-scalar-list" else True)
        )
        row = {"case": name, "pass": passed, "expected": expected}
        if legacy is not None:
            row["repair1_escape_reproduced"] = legacy_escape
        rows.append(row)

    failed = list(prior.get("failed_cases", [])) + [row["case"] for row in rows if not row["pass"]]
    result = {
        "suite_id": "GRA-O-GF-CONTRACT-V002-REPAIR2-MUTATIONS",
        "overall_result": "PASS" if prior.get("overall_result") == "PASS" and not failed else "FAIL",
        "case_count": prior.get("case_count", 0) + len(rows),
        "prior_repair1_cases": prior.get("case_count"),
        "repair2_observation_cases": len(rows),
        "typed_lifecycle_positive_control": positive_ok,
        "failed_count": len(failed), "failed_cases": failed,
        "prior_repair1_result_payload_sha256": prior.get("result_payload_sha256"),
        "repair2_cases": rows,
    }
    result["result_payload_sha256"] = hashlib.sha256(base.canonical_bytes(result)).hexdigest()
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    result = run_suite()
    print(json.dumps(result, sort_keys=True, indent=None if args.compact else 2))
    return 0 if result["overall_result"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
