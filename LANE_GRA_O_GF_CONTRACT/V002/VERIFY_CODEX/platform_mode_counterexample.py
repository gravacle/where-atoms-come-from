#!/usr/bin/env python3
"""Read-only in-memory V002 platform/mode isolation counterexample."""

from __future__ import annotations

import base64
import copy
import hashlib
import json
import sys
from pathlib import Path


LANE = Path(__file__).resolve().parent.parent
WORKSPACE = LANE.parents[1]
sys.path.insert(0, str(LANE))

import synthetic_positive_fixture_v002 as fixture  # noqa: E402
import validator_v002 as v  # noqa: E402


def inline_artifact(aid: str, kind: str, role: str, payload: dict) -> dict:
    raw = v.canonical_bytes(payload)
    return {
        "artifact_id": aid,
        "kind": kind,
        "role": role,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "byte_length": len(raw),
        "storage": "INLINE_JSON",
        "locator": None,
        "payload_b64": base64.b64encode(raw).decode("ascii"),
        "parents": [],
    }


def main() -> int:
    baseline = fixture.build_fixture()
    baseline_result = v.evaluate_instance(copy.deepcopy(baseline), WORKSPACE)

    candidate = copy.deepcopy(baseline)
    candidate["package_mode"] = "SCIENTIFIC"
    platform_id = "PLATFORM.EMPTY-ASSERTION"
    candidate["artifacts"].append(
        inline_artifact(
            platform_id,
            "PLATFORM_INSTANTIATION",
            "DEVELOPMENT",
            {
                "platform_id": "",
                "surface_ids": [],
                "synthetic_only": False,
                "platform_map_artifact_ids": [],
                "freeze_time": "",
            },
        )
    )
    candidate["bindings"]["platform_instantiation_artifact_id"] = platform_id
    result = v.evaluate_instance(candidate, WORKSPACE)

    promoted = {
        claim: value
        for claim, value in result.get("authoritative_proof_outputs", {}).items()
        if value != "NO_PROOF_OUTPUT"
    }
    decisive = (
        baseline_result.get("accepted") is True
        and baseline_result.get("actual_platform_present") is False
        and all(
            value == "NO_PROOF_OUTPUT"
            for value in baseline_result.get("authoritative_proof_outputs", {}).values()
        )
        and result.get("accepted") is True
        and result.get("custody", {}).get("disposition") == "QUALIFIED"
        and result.get("actual_platform_present") is True
        and bool(promoted)
    )
    report = {
        "audit_id": "GRA-O-GF-CONTRACT-V002-PLATFORM-MODE-COUNTEREXAMPLE",
        "contract_disposition": "REFUTED" if decisive else "NO_COUNTEREXAMPLE",
        "mutation": {
            "package_mode": "SYNTHETIC_TEST -> SCIENTIFIC",
            "binding": "platform_instantiation_artifact_id -> PLATFORM.EMPTY-ASSERTION",
            "added_platform_payload": {
                "platform_id": "",
                "surface_ids": [],
                "synthetic_only": False,
                "platform_map_artifact_ids": [],
                "freeze_time": "",
            },
        },
        "baseline": {
            "accepted": baseline_result.get("accepted"),
            "custody": baseline_result.get("custody", {}).get("disposition"),
            "actual_platform_present": baseline_result.get("actual_platform_present"),
            "authoritative_proof_outputs": baseline_result.get("authoritative_proof_outputs"),
        },
        "mutant": {
            "accepted": result.get("accepted"),
            "custody": result.get("custody", {}).get("disposition"),
            "custody_errors": result.get("custody", {}).get("errors"),
            "actual_platform_present": result.get("actual_platform_present"),
            "candidate_milestones": result.get("candidate_milestones"),
            "authoritative_milestones": result.get("authoritative_milestones"),
            "authoritative_proof_outputs": result.get("authoritative_proof_outputs"),
            "promoted_outputs": promoted,
            "product_reproduction": result.get("product_reproduction"),
            "semantic_diagnostics": result.get("semantic_diagnostics"),
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if decisive else 1


if __name__ == "__main__":
    raise SystemExit(main())
