#!/usr/bin/env python3
"""Bounded read-only structural checks for the V002 lane."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import synthetic_positive_fixture_v002 as fixture
import validator_v002 as v


WORKSPACE = Path(__file__).resolve().parents[2]
LANE = Path(__file__).resolve().parent


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(path: Path) -> tuple[int, int, list[str]]:
    matched = 0
    rows = 0
    failures: list[str] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        rows += 1
        expected, relative = line.split(None, 1)
        target = WORKSPACE / relative.strip()
        if target.is_file() and digest(target) == expected:
            matched += 1
        else:
            failures.append(relative.strip())
    return matched, rows, failures


def local_refs(value: Any) -> list[str]:
    refs: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "$ref" and isinstance(child, str) and child.startswith("#/$defs/"):
                refs.append(child.removeprefix("#/$defs/"))
            refs.extend(local_refs(child))
    elif isinstance(value, list):
        for child in value:
            refs.extend(local_refs(child))
    return refs


def main() -> int:
    failures: list[str] = []
    schema = json.loads((LANE / "INSTANCE_V002.schema.json").read_text())
    defs = set(schema.get("$defs", {}))
    refs = local_refs(schema)
    unresolved = sorted(set(refs) - defs)
    if unresolved:
        failures.append(f"unresolved schema refs: {unresolved}")
    if schema.get("additionalProperties") is not False:
        failures.append("top schema is open")
    if set(schema.get("required", [])) != set(schema.get("properties", {})):
        failures.append("top required/properties mismatch")
    v001_matched, v001_total, v001_failures = verify_manifest(WORKSPACE / "LANE_GRA_O_GF_CONTRACT/MANIFEST.sha256")
    verify_matched, verify_total, verify_failures = verify_manifest(WORKSPACE / "LANE_GRA_O_GF_CONTRACT/VERIFY_CODEX/MANIFEST.sha256")
    if v001_failures:
        failures.append(f"V001 modified: {v001_failures}")
    if verify_failures:
        failures.append(f"V001 verifier modified: {verify_failures}")
    source_matches = 0
    for relative, expected in v.REQUIRED_FRAMEWORKS.items():
        if digest(WORKSPACE / relative) == expected:
            source_matches += 1
        else:
            failures.append(f"framework mismatch: {relative}")
    if digest(WORKSPACE / v.PRINCIPAL_DECISION_PATH) != v.PRINCIPAL_DECISION_SHA256:
        failures.append("principal decision mismatch")
    baseline = v.evaluate_instance(fixture.build_fixture(), WORKSPACE)
    if not baseline.get("accepted") or baseline.get("custody", {}).get("disposition") != "QUALIFIED":
        failures.append("synthetic baseline not qualified")
    if baseline.get("actual_platform_present") is not False:
        failures.append("synthetic baseline claims actual platform")
    if any(value != "NO_PROOF_OUTPUT" for value in baseline.get("authoritative_proof_outputs", {}).values()):
        failures.append("synthetic baseline promoted proof")
    result = {
        "check_id": "GRA-O-GF-CONTRACT-V002-STRUCTURAL",
        "overall_result": "PASS" if not failures else "FAIL",
        "schema_draft": schema.get("$schema"),
        "schema_definitions": len(defs),
        "schema_local_refs": len(refs),
        "schema_unresolved_refs": unresolved,
        "v001_manifest": {"matched": v001_matched, "total": v001_total},
        "v001_verifier_manifest": {"matched": verify_matched, "total": verify_total},
        "frozen_frameworks": {"matched": source_matches, "total": len(v.REQUIRED_FRAMEWORKS)},
        "principal_decision_sha256": v.PRINCIPAL_DECISION_SHA256,
        "synthetic_baseline": {
            "accepted": baseline.get("accepted"),
            "custody": baseline.get("custody", {}).get("disposition"),
            "actual_platform_present": baseline.get("actual_platform_present"),
            "proof_outputs": baseline.get("authoritative_proof_outputs"),
        },
        "failures": failures,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
