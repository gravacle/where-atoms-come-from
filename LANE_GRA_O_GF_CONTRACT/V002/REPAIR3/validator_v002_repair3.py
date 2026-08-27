#!/usr/bin/env python3
"""Narrow Repair3 overlay: immutable, monotone Repair1 + observation evaluation.

The original instance is passed unchanged to Repair1 exactly once.  Repair2's typed
formation predicate reads that same original instance.  Combining results can only
add diagnostics and remove promotion; no artifact, payload, digest, parent, binding,
or candidate is rewritten or re-evaluated.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any


REPAIR3_DIR = Path(__file__).resolve().parent
REPAIR2_DIR = REPAIR3_DIR.parent / "REPAIR2"
REPAIR1_DIR = REPAIR3_DIR.parent / "REPAIR1"
V002_DIR = REPAIR3_DIR.parent
for path in (REPAIR2_DIR, REPAIR1_DIR, V002_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import validator_v002 as base  # noqa: E402
import validator_v002_repair1 as repair1  # noqa: E402
import validator_v002_repair2 as repair2  # noqa: E402


REPAIR_ID = "GRA-O-GF-CONTRACT-V002-REPAIR3"
REPAIR1_EVALUATE_ORIGINAL = repair1.evaluate_instance
OBSERVATION_VALIDATE_ORIGINAL = repair2.validate_actual_formation_observations


def _monotone_refuse(result: dict[str, Any], diagnostics: list[str]) -> dict[str, Any]:
    combined = copy.deepcopy(result)
    custody = combined.setdefault("custody", {})
    custody["disposition"] = "REFUSE"
    custody["errors"] = sorted(set(custody.get("errors", []) + diagnostics))
    combined["accepted"] = False
    combined["actual_platform_present"] = False
    combined["semantic_diagnostics"] = sorted(set(combined.get("semantic_diagnostics", []) + diagnostics))
    combined["scientific_gates"] = {gate: "UNSCOREABLE" for gate in base.GATES}
    combined["candidate_milestones"] = {claim: "UNSCOREABLE" for claim in base.CLAIMS}
    combined["authoritative_milestones"] = {claim: "UNSCOREABLE" for claim in base.CLAIMS}
    combined["authoritative_proof_outputs"] = {claim: "NO_PROOF_OUTPUT" for claim in base.CLAIMS}
    # Preserve Repair2's public diagnostic shape without invoking its projecting
    # evaluator; both fields report the one read-only check of the original.
    combined["formation_observation_repair2"] = {"valid": False, "diagnostics": diagnostics}
    combined["formation_observation_repair3"] = {"valid": False, "diagnostics": diagnostics}
    return combined


def evaluate_instance(instance: Any, workspace: Path | None = None) -> dict[str, Any]:
    workspace = (workspace or Path.cwd()).resolve()

    # Custody/platform authority: exactly the original bytes and links, exactly once.
    repair1_result = REPAIR1_EVALUATE_ORIGINAL(instance, workspace)
    if not isinstance(instance, dict):
        return repair1_result

    # Observation authority: read-only evaluation of the same original candidate.
    observations_valid, observation_diagnostics = OBSERVATION_VALIDATE_ORIGINAL(instance, workspace)

    # Monotone conjunction: the new layer cannot repair any Repair1 failure.
    if repair1_result.get("accepted") is not True:
        combined = copy.deepcopy(repair1_result)
        combined["actual_platform_present"] = False
        combined["formation_observation_repair2"] = {
            "valid": observations_valid,
            "diagnostics": observation_diagnostics,
        }
        combined["formation_observation_repair3"] = {
            "valid": observations_valid,
            "diagnostics": observation_diagnostics,
            "upstream_repair1_accepted": False,
        }
        return combined
    if not observations_valid:
        return _monotone_refuse(repair1_result, observation_diagnostics)

    combined = copy.deepcopy(repair1_result)
    combined["formation_observation_repair2"] = {"valid": True, "diagnostics": []}
    combined["formation_observation_repair3"] = {
        "valid": True,
        "diagnostics": [],
        "upstream_repair1_accepted": True,
    }
    return combined


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    args = parser.parse_args()
    try:
        candidate = base.load_candidate(args.candidate)
        result = evaluate_instance(candidate, args.workspace)
    except Exception as exc:
        result = {
            "accepted": False,
            "custody": {"disposition": "REFUSE", "errors": [f"REPAIR3_ERROR:{type(exc).__name__}:{exc}"], "missing_artifact_ids": []},
            "actual_platform_present": False,
            "authoritative_proof_outputs": {claim: "NO_PROOF_OUTPUT" for claim in base.CLAIMS},
        }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if result.get("accepted") else 2


if __name__ == "__main__":
    raise SystemExit(main())
