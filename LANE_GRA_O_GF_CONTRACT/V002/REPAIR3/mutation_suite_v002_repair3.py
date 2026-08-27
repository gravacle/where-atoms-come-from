#!/usr/bin/env python3
"""Inherited 54 cases plus five immutable-custody Repair3 regressions."""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Callable


REPAIR3_DIR = Path(__file__).resolve().parent
REPAIR2_DIR = REPAIR3_DIR.parent / "REPAIR2"
REPAIR1_DIR = REPAIR3_DIR.parent / "REPAIR1"
V002_DIR = REPAIR3_DIR.parent
for path in (REPAIR2_DIR, REPAIR1_DIR, V002_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import mutation_suite_v002 as sealed_suite  # noqa: E402
import mutation_suite_v002_repair1 as repair1_suite  # noqa: E402
import mutation_suite_v002_repair2 as repair2_suite  # noqa: E402
import validator_v002 as base  # noqa: E402
import validator_v002_repair2 as repair2  # noqa: E402
import validator_v002_repair3 as repair3  # noqa: E402


WORKSPACE = Path(__file__).resolve().parents[3]
FORMATION_ID = repair2_suite.FORMATION_ID
SOURCE_ID = "REAL.PLATFORM.SOURCE"
PROMOTED_CLAIMS = ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2")


def _formation(candidate: dict[str, Any]) -> dict[str, Any]:
    return sealed_suite.artifact(candidate, FORMATION_ID)


def mutate_corrupted_parent_digest(candidate: dict[str, Any]) -> None:
    parent = next(item for item in _formation(candidate)["parents"] if item["artifact_id"] == SOURCE_ID)
    parent["sha256"] = "0" * 64


def mutate_wrong_child_digest(candidate: dict[str, Any]) -> None:
    _formation(candidate)["sha256"] = "f" * 64


def mutate_unknown_parent(candidate: dict[str, Any]) -> None:
    parent = next(item for item in _formation(candidate)["parents"] if item["artifact_id"] == SOURCE_ID)
    parent["artifact_id"] = "UNKNOWN.PARENT"
    parent["sha256"] = "0" * 64


def mutate_parent_cycle(candidate: dict[str, Any]) -> None:
    child = _formation(candidate)
    child["parents"].append({"artifact_id": FORMATION_ID, "sha256": child["sha256"]})


def mutate_orphan(candidate: dict[str, Any]) -> None:
    content = b"unbound measured bytes used only for orphan regression"
    repair1_suite.add_inline(candidate, "REPAIR3.ORPHAN.SOURCE", "SOURCE_DATA", "GENERIC", {
        "source_id": "REPAIR3.ORPHAN.SOURCE",
        "content_b64": base64.b64encode(content).decode(),
        "content_sha256": hashlib.sha256(content).hexdigest(),
    })


Case = tuple[str, Callable[[dict[str, Any]], None], str]
CASES: list[Case] = [
    ("REPAIR3-exact-corrupted-parent-digest", mutate_corrupted_parent_digest, "ARTIFACT_PARENT_HASH_MISMATCH"),
    ("REPAIR3-wrong-child-digest", mutate_wrong_child_digest, "ARTIFACT_HASH_MISMATCH:REAL.PLATFORM.FORMATION"),
    ("REPAIR3-unknown-parent", mutate_unknown_parent, "ARTIFACT_PARENT_UNRESOLVED"),
    ("REPAIR3-parent-cycle", mutate_parent_cycle, "ARTIFACT_DAG_SELF_LOOP"),
    ("REPAIR3-orphan", mutate_orphan, "ARTIFACT_ORPHAN:REPAIR3.ORPHAN.SOURCE"),
]


def no_proof(result: dict[str, Any]) -> bool:
    outputs = result.get("authoritative_proof_outputs", {})
    return bool(outputs) and all(value == "NO_PROOF_OUTPUT" for value in outputs.values())


def promoted(result: dict[str, Any]) -> bool:
    return all(
        result.get("authoritative_proof_outputs", {}).get(claim) == "PASSES_DECLARED_DOMAIN"
        for claim in PROMOTED_CLAIMS
    )


def run_prior_54_through_repair3() -> dict[str, Any]:
    """Exercise the sealed Repair2 inventory with Repair3 as its current evaluator."""
    original = repair2.evaluate_instance
    repair2.evaluate_instance = repair3.evaluate_instance
    try:
        return repair2_suite.run_suite()
    finally:
        repair2.evaluate_instance = original


def run_suite() -> dict[str, Any]:
    prior = run_prior_54_through_repair3()
    prior_rows = {row["case"]: row for row in prior.get("repair2_cases", [])}
    scalar_confirmed = prior_rows.get("REPAIR2-exact-scalar-list", {}).get("pass") is True
    mixed_confirmed = prior_rows.get("REPAIR2-mixed-dict-scalar", {}).get("pass") is True

    positive = repair2_suite.typed_platform_candidate()
    positive_before = base.canonical_bytes(positive)
    positive_result = repair3.evaluate_instance(positive, WORKSPACE)
    typed_positive = (
        positive_result.get("accepted") is True
        and positive_result.get("actual_platform_present") is True
        and positive_result.get("formation_observation_repair3", {}).get("valid") is True
        and promoted(positive_result)
        and base.canonical_bytes(positive) == positive_before
    )

    rows: list[dict[str, Any]] = []
    exact_laundering_reproduced = False
    for name, mutation, expected in CASES:
        candidate = repair2_suite.typed_platform_candidate()
        mutation(candidate)
        candidate_before = base.canonical_bytes(candidate)
        legacy = repair2.evaluate_instance(copy.deepcopy(candidate), WORKSPACE) if name == "REPAIR3-exact-corrupted-parent-digest" else None
        result = repair3.evaluate_instance(candidate, WORKSPACE)
        result_text = json.dumps(result, sort_keys=True, allow_nan=False)
        immutable = base.canonical_bytes(candidate) == candidate_before
        legacy_escape = bool(
            legacy
            and legacy.get("accepted") is True
            and legacy.get("actual_platform_present") is True
            and promoted(legacy)
        )
        if legacy is not None:
            exact_laundering_reproduced = legacy_escape
        passed = (
            typed_positive
            and result.get("accepted") is False
            and result.get("actual_platform_present") is False
            and expected in result_text
            and no_proof(result)
            and immutable
            and (legacy_escape if legacy is not None else True)
        )
        row: dict[str, Any] = {
            "case": name,
            "pass": passed,
            "expected": expected,
            "input_byte_identity_preserved": immutable,
        }
        if legacy is not None:
            row["repair2_laundering_reproduced"] = legacy_escape
        rows.append(row)

    failed = list(prior.get("failed_cases", [])) + [row["case"] for row in rows if not row["pass"]]
    controls_ok = scalar_confirmed and mixed_confirmed and typed_positive and exact_laundering_reproduced
    if not controls_ok:
        failed.append("REPAIR3-required-controls")
    failed = sorted(set(failed))
    result = {
        "suite_id": "GRA-O-GF-CONTRACT-V002-REPAIR3-MUTATIONS",
        "overall_result": "PASS" if prior.get("overall_result") == "PASS" and not failed else "FAIL",
        "case_count": prior.get("case_count", 0) + len(rows),
        "prior_repair2_cases": prior.get("case_count"),
        "repair3_custody_cases": len(rows),
        "typed_lifecycle_positive_control": typed_positive,
        "scalar_list_regression_confirmed": scalar_confirmed,
        "mixed_list_regression_confirmed": mixed_confirmed,
        "repair2_exact_parent_laundering_reproduced": exact_laundering_reproduced,
        "failed_count": len(failed),
        "failed_cases": failed,
        "prior_repair2_result_payload_sha256": prior.get("result_payload_sha256"),
        "repair3_cases": rows,
    }
    result["result_payload_sha256"] = hashlib.sha256(base.canonical_bytes(result)).hexdigest()
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    result = run_suite()
    print(json.dumps(result, sort_keys=True, indent=None if args.compact else 2, allow_nan=False))
    return 0 if result["overall_result"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
