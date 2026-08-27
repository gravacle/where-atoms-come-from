#!/usr/bin/env python3
"""Sealed V002 suite plus the four narrow PLATFORM Repair1 regressions."""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import json
import sys
from pathlib import Path
from typing import Any


REPAIR_DIR = Path(__file__).resolve().parent
PARENT = REPAIR_DIR.parent
if str(PARENT) not in sys.path:
    sys.path.insert(0, str(PARENT))

import mutation_suite_v002 as sealed_suite  # noqa: E402
import synthetic_positive_fixture_v002 as fixture  # noqa: E402
import validator_v002 as base  # noqa: E402
import validator_v002_repair1 as repair  # noqa: E402


WORKSPACE = Path(__file__).resolve().parents[3]


def add_inline(candidate: dict[str, Any], aid: str, kind: str, role: str, payload: dict[str, Any]) -> str:
    registry = {item["artifact_id"]: item for item in candidate["artifacts"]}
    refs = sorted(base.referenced_artifact_ids(payload))
    raw = base.canonical_bytes(payload)
    artifact = {
        "artifact_id": aid, "kind": kind, "role": role,
        "sha256": hashlib.sha256(raw).hexdigest(), "byte_length": len(raw),
        "storage": "INLINE_JSON", "locator": None,
        "payload_b64": base64.b64encode(raw).decode(),
        "parents": [{"artifact_id": ref, "sha256": registry.get(ref, {"sha256": "0" * 64})["sha256"]} for ref in refs],
    }
    candidate["artifacts"].append(artifact)
    return aid


def replace_payload(candidate: dict[str, Any], aid: str, payload: dict[str, Any]) -> None:
    sealed_suite.replace_payload(candidate, aid, payload)


def change_payload(candidate: dict[str, Any], aid: str, updates: dict[str, Any]) -> None:
    value = sealed_suite.payload(candidate, aid)
    value.update(updates)
    replace_payload(candidate, aid, value)


def exact_counterexample() -> dict[str, Any]:
    candidate = fixture.build_fixture()
    candidate["package_mode"] = "SCIENTIFIC"
    platform = {
        "platform_id": "", "surface_ids": [], "synthetic_only": False,
        "platform_map_artifact_ids": [], "freeze_time": "",
    }
    add_inline(candidate, "PLATFORM.ATTACK.EXACT", "PLATFORM_INSTANTIATION", "DEVELOPMENT", platform)
    candidate["bindings"]["platform_instantiation_artifact_id"] = "PLATFORM.ATTACK.EXACT"
    return candidate


def source_leaf_ids(candidate: dict[str, Any], roots: set[str]) -> list[str]:
    registry = {item["artifact_id"]: item for item in candidate["artifacts"]}
    parents = {aid: {parent["artifact_id"] for parent in item["parents"]} for aid, item in registry.items()}
    leaves: set[str] = set()
    visited: set[str] = set()
    stack = list(roots)
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        if registry.get(current, {}).get("kind") == "SOURCE_DATA":
            leaves.add(current)
        stack.extend(parents.get(current, ()))
    return sorted(leaves)


def joined_platform_candidate() -> dict[str, Any]:
    candidate = fixture.build_fixture()
    candidate["package_mode"] = "SCIENTIFIC"
    package_id = candidate["package_id"]
    bindings = candidate["bindings"]
    platform_id = "REAL-PLATFORM-REPAIR1"

    change_payload(candidate, "DILATION.MEASURE", {
        "source_id": "DILATION.MEASURE",
        "content_b64": base64.b64encode(b"measured full exchange stream").decode(),
        "content_sha256": hashlib.sha256(b"measured full exchange stream").hexdigest(),
    })
    change_payload(candidate, "DILATION.FULL", {"surface_id": "SURFACE-0"})
    change_payload(candidate, "TRANSPORT.V004", {"surface_id": "SURFACE-0"})

    real_source = add_inline(candidate, "REAL.PLATFORM.SOURCE", "SOURCE_DATA", "DEVELOPMENT", {
        "source_id": "REAL.PLATFORM.SOURCE",
        "content_b64": base64.b64encode(b"measured physical platform bytes").decode(),
        "content_sha256": hashlib.sha256(b"measured physical platform bytes").hexdigest(),
    })
    formation = add_inline(candidate, "REAL.PLATFORM.FORMATION", "GATE_EVIDENCE", "DEVELOPMENT", {
        "gate_id": "FORM.ALLOW0",
        "observations": [{"predicate_id": "REAL.FORMATION", "value": True, "reproducible": True}],
        "source_artifact_ids": [real_source], "taxonomy_match": "KNOWN", "reproducible": True,
    })
    roots = set(bindings["gamma_process_artifact_ids"]) | {
        bindings["seed_definition_artifact_id"], bindings["dilation_artifact_id"],
        bindings["transport_artifact_id"], formation,
    }
    custody_sources = source_leaf_ids(candidate, roots)
    custody = add_inline(candidate, "REAL.PLATFORM.CUSTODY", "PLATFORM_SOURCE_CUSTODY", "DEVELOPMENT", {
        "platform_id": platform_id, "package_id": package_id,
        "source_artifact_ids": custody_sources, "acquisition_ids": ["REAL-ACQUISITION"],
        "specimen_ids": ["REAL-SPECIMEN"], "independent_unit_ids": ["REAL-UNIT"],
        "provenance": "content-addressed measured platform acquisition", "license_id": "REAL-DATA-LICENSE",
        "synthetic_test": False,
    })
    process_by_surface = {
        sealed_suite.payload(candidate, process_id)["surface_id"]: process_id
        for process_id in bindings["gamma_process_artifact_ids"]
    }
    map_ids: list[str] = []
    for index, (surface, process_id) in enumerate(sorted(process_by_surface.items())):
        map_ids.append(add_inline(candidate, f"REAL.PLATFORM.MAP.{index}", "PLATFORM_MAP", "DEVELOPMENT", {
            "platform_id": platform_id, "package_id": package_id, "surface_id": surface,
            "process_artifact_id": process_id,
            "seed_definition_artifact_id": bindings["seed_definition_artifact_id"],
            "dilation_artifact_id": bindings["dilation_artifact_id"],
            "transport_artifact_id": bindings["transport_artifact_id"],
            "formation_evidence_artifact_ids": [formation],
            "source_custody_artifact_ids": [custody],
        }))
    platform = {
        "platform_id": platform_id, "package_id": package_id,
        "surface_ids": sorted(process_by_surface), "synthetic_only": False,
        "platform_map_artifact_ids": map_ids, "formation_evidence_artifact_ids": [formation],
        "gamma_process_artifact_ids": list(bindings["gamma_process_artifact_ids"]),
        "seed_definition_artifact_id": bindings["seed_definition_artifact_id"],
        "dilation_artifact_id": bindings["dilation_artifact_id"],
        "transport_artifact_id": bindings["transport_artifact_id"],
        "source_custody_artifact_ids": [custody], "freeze_time": "2026-08-22T09:45:00Z",
    }
    platform_aid = add_inline(candidate, "REAL.PLATFORM.INSTANCE", "PLATFORM_INSTANTIATION", "DEVELOPMENT", platform)
    candidate["bindings"]["platform_instantiation_artifact_id"] = platform_aid
    return candidate


def repair_case(name: str, candidate: dict[str, Any], expected: str) -> dict[str, Any]:
    result = repair.evaluate_instance(candidate, WORKSPACE)
    text = json.dumps(result, sort_keys=True)
    proof = result.get("authoritative_proof_outputs", {})
    passed = (
        result.get("accepted") is False
        and result.get("actual_platform_present") is False
        and expected in text
        and proof and all(value == "NO_PROOF_OUTPUT" for value in proof.values())
    )
    return {"case": name, "pass": passed, "expected": expected}


def run_suite() -> dict[str, Any]:
    original_evaluate = base.evaluate_instance
    base.evaluate_instance = repair.evaluate_instance
    try:
        sealed = sealed_suite.run_suite()
    finally:
        base.evaluate_instance = original_evaluate

    rows: list[dict[str, Any]] = []
    exact = exact_counterexample()
    legacy = repair.evaluate_legacy_exact(copy.deepcopy(exact), WORKSPACE)
    repaired = repair.evaluate_instance(copy.deepcopy(exact), WORKSPACE)
    legacy_promoted = (
        legacy.get("accepted") is True and legacy.get("actual_platform_present") is True
        and all(legacy.get("authoritative_proof_outputs", {}).get(claim) == "PASSES_DECLARED_DOMAIN" for claim in ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2"))
    )
    exact_pass = (
        legacy_promoted and repaired.get("accepted") is False
        and repaired.get("actual_platform_present") is False
        and "PLATFORM_INSTANCE_SCHEMA_INVALID" in json.dumps(repaired)
        and all(value == "NO_PROOF_OUTPUT" for value in repaired.get("authoritative_proof_outputs", {}).values())
    )
    rows.append({"case": "REPAIR1-exact-accepted-counterexample", "pass": exact_pass, "legacy_refutation_reproduced": legacy_promoted, "expected": "PLATFORM_INSTANCE_SCHEMA_INVALID"})

    empty = joined_platform_candidate()
    change_payload(empty, "REAL.PLATFORM.INSTANCE", {
        **sealed_suite.payload(empty, "REAL.PLATFORM.INSTANCE"),
        "platform_id": "", "surface_ids": [], "platform_map_artifact_ids": [],
    })
    rows.append(repair_case("REPAIR1-empty-platform-fields", empty, "PLATFORM_ID_EMPTY_OR_INVALID"))

    unknown = joined_platform_candidate()
    unknown_payload = sealed_suite.payload(unknown, "REAL.PLATFORM.INSTANCE")
    unknown_payload["surface_ids"] = unknown_payload["surface_ids"] + ["SURFACE-UNKNOWN"]
    replace_payload(unknown, "REAL.PLATFORM.INSTANCE", unknown_payload)
    rows.append(repair_case("REPAIR1-unknown-surface", unknown, "PLATFORM_SURFACE_UNKNOWN_OR_INCOMPLETE"))

    cross = joined_platform_candidate()
    change_payload(cross, "REAL.PLATFORM.INSTANCE", {"package_id": "OTHER.PACKAGE"})
    rows.append(repair_case("REPAIR1-cross-package", cross, "PLATFORM_CROSS_PACKAGE_BINDING"))

    failed = list(sealed.get("failed_cases", [])) + [row["case"] for row in rows if not row["pass"]]
    result = {
        "suite_id": "GRA-O-GF-CONTRACT-V002-REPAIR1-MUTATIONS",
        "overall_result": "PASS" if sealed.get("overall_result") == "PASS" and not failed else "FAIL",
        "case_count": sealed.get("case_count", 0) + len(rows),
        "sealed_v002_cases": sealed.get("case_count"),
        "repair1_platform_cases": len(rows),
        "failed_count": len(failed), "failed_cases": failed,
        "sealed_v002_result_payload_sha256": sealed.get("result_payload_sha256"),
        "repair1_cases": rows,
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
