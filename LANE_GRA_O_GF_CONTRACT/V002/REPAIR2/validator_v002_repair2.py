#!/usr/bin/env python3
"""Narrow Repair2 overlay: all-member typed formation observations.

Sealed V002 and Repair1 bytes remain unchanged.  This layer validates every member
of each actual-platform formation observation list before reading any member field.
It adds no geometric, metric, stress-energy, or IR-source requirement.
"""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


REPAIR2_DIR = Path(__file__).resolve().parent
REPAIR1_DIR = REPAIR2_DIR.parent / "REPAIR1"
V002_DIR = REPAIR2_DIR.parent
for path in (REPAIR1_DIR, V002_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import validator_v002 as base  # noqa: E402
import validator_v002_repair1 as repair1  # noqa: E402


REPAIR_ID = "GRA-O-GF-CONTRACT-V002-REPAIR2"
OBSERVATION_KEYS = {
    "predicate_id", "event_id", "surface_id", "stage", "role", "time",
    "value", "unit", "source_artifact_id", "reproducible",
}
LIFECYCLE = (
    ("FORMATION", "FORMATION"),
    ("CLOSURE", "CLOSURE"),
    ("PERSISTENCE", "PERSISTENCE"),
)
LIFECYCLE_STAGES = {stage for stage, _ in LIFECYCLE}
ROLE_BY_STAGE = dict(LIFECYCLE)

REPAIR1_EVALUATE = repair1.evaluate_instance


def _time(value: Any) -> datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def _artifact_payload(aid: str, registry: dict[str, dict[str, Any]], payloads: dict[str, Any]) -> dict[str, Any] | None:
    if registry.get(aid, {}).get("kind") != "GATE_EVIDENCE":
        return None
    payload = payloads.get(aid)
    return payload if isinstance(payload, dict) else None


def validate_actual_formation_observations(instance: dict[str, Any], workspace: Path) -> tuple[bool, list[str]]:
    diagnostics: list[str] = []
    bindings = instance.get("bindings")
    if not isinstance(bindings, dict) or bindings.get("platform_instantiation_artifact_id") is None:
        return True, []
    registry, payloads, artifact_errors, _ = base.validate_artifacts(instance, workspace)
    platform_aid = bindings.get("platform_instantiation_artifact_id")
    platform = payloads.get(platform_aid)
    # Old/incomplete platform schemas remain Repair1's responsibility.
    if not isinstance(platform, dict) or "formation_evidence_artifact_ids" not in platform:
        return True, []
    formation_ids = platform.get("formation_evidence_artifact_ids")
    if not isinstance(formation_ids, list) or not formation_ids:
        return False, ["FORMATION_EVIDENCE_REGISTRY_EMPTY"]

    freeze_time = _time(platform.get("freeze_time"))
    platform_surfaces = set(platform.get("surface_ids") or [])
    for evidence_id in formation_ids:
        evidence = _artifact_payload(evidence_id, registry, payloads)
        if evidence is None:
            if any(evidence_id in error and "ARTIFACT_JSON_INVALID" in error for error in artifact_errors):
                diagnostics.append(f"FORMATION_OBSERVATION_NONFINITE_OR_INVALID_JSON:{evidence_id}")
            else:
                diagnostics.append(f"FORMATION_OBSERVATION_PAYLOAD_INVALID:{evidence_id}")
            continue
        if evidence.get("gate_id") != "FORM.ALLOW0":
            diagnostics.append(f"FORMATION_OBSERVATION_GATE_INVALID:{evidence_id}")
        observations = evidence.get("observations")
        if not isinstance(observations, list) or not observations:
            diagnostics.append(f"FORMATION_OBSERVATION_LIST_EMPTY:{evidence_id}")
            continue

        # Critical Repair2 order: validate every member type before any `.get` call.
        if not all(isinstance(member, dict) for member in observations):
            diagnostics.append(f"FORMATION_OBSERVATION_MEMBER_NOT_OBJECT:{evidence_id}")
            continue
        if not all(set(member) == OBSERVATION_KEYS for member in observations):
            diagnostics.append(f"FORMATION_OBSERVATION_KEYS_INVALID:{evidence_id}")
            continue

        typed: list[dict[str, Any]] = observations
        evidence_sources = set(evidence.get("source_artifact_ids") or [])
        for index, member in enumerate(typed):
            prefix = f"{evidence_id}:{index}"
            for key in ("predicate_id", "event_id", "surface_id", "stage", "role", "time", "unit", "source_artifact_id"):
                if not isinstance(member[key], str) or not member[key].strip():
                    diagnostics.append(f"FORMATION_OBSERVATION_IDENTITY_MISSING:{prefix}:{key}")
            value = member["value"]
            if isinstance(value, bool):
                finite_value = True
            else:
                finite_value = isinstance(value, (int, float)) and math.isfinite(float(value))
            if not finite_value:
                diagnostics.append(f"FORMATION_OBSERVATION_VALUE_NONFINITE_OR_INVALID:{prefix}")
            if value is not True:
                diagnostics.append(f"FORMATION_OBSERVATION_VALUE_NOT_CERTIFIED:{prefix}")
            if member["reproducible"] is not True:
                diagnostics.append(f"FORMATION_OBSERVATION_NOT_REPRODUCIBLE:{prefix}")
            if member["stage"] not in LIFECYCLE_STAGES or ROLE_BY_STAGE.get(member["stage"]) != member["role"]:
                diagnostics.append(f"FORMATION_OBSERVATION_STAGE_ROLE_INVALID:{prefix}")
            if member["surface_id"] not in platform_surfaces:
                diagnostics.append(f"FORMATION_OBSERVATION_SURFACE_UNKNOWN:{prefix}")
            if member["source_artifact_id"] not in evidence_sources:
                diagnostics.append(f"FORMATION_OBSERVATION_SOURCE_JOIN_INVALID:{prefix}")
            observed_time = _time(member["time"])
            if observed_time is None or freeze_time is None or observed_time >= freeze_time:
                diagnostics.append(f"FORMATION_OBSERVATION_TIME_INVALID:{prefix}")

        if {member["stage"] for member in typed} != LIFECYCLE_STAGES:
            diagnostics.append(f"FORMATION_LIFECYCLE_COVERAGE_INVALID:{evidence_id}")
        if len(typed) != len(LIFECYCLE):
            diagnostics.append(f"FORMATION_LIFECYCLE_CARDINALITY_INVALID:{evidence_id}")
        if len({member["predicate_id"] for member in typed}) != len(typed):
            diagnostics.append(f"FORMATION_PREDICATE_ID_DUPLICATE:{evidence_id}")
        if len({member["event_id"] for member in typed}) != 1:
            diagnostics.append(f"FORMATION_LIFECYCLE_CROSS_EVENT:{evidence_id}")
        if len({member["surface_id"] for member in typed}) != 1:
            diagnostics.append(f"FORMATION_LIFECYCLE_CROSS_SURFACE:{evidence_id}")
        if len({member["unit"] for member in typed}) != 1:
            diagnostics.append(f"FORMATION_LIFECYCLE_UNIT_UNSTABLE:{evidence_id}")
        if len({member["source_artifact_id"] for member in typed}) != 1:
            diagnostics.append(f"FORMATION_LIFECYCLE_SOURCE_UNSTABLE:{evidence_id}")
        by_stage = {member["stage"]: _time(member["time"]) for member in typed if member["stage"] in LIFECYCLE_STAGES}
        ordered_times = [by_stage.get(stage) for stage, _ in LIFECYCLE]
        if any(value is None for value in ordered_times) or not all(ordered_times[i] < ordered_times[i + 1] for i in range(len(ordered_times) - 1)):
            diagnostics.append(f"FORMATION_LIFECYCLE_TIME_ORDER_INVALID:{evidence_id}")
    return not diagnostics, sorted(set(diagnostics))


def _project_for_repair1(instance: dict[str, Any]) -> dict[str, Any]:
    projected = copy.deepcopy(instance)
    bindings = projected.get("bindings")
    if not isinstance(bindings, dict):
        return projected
    platform_aid = bindings.get("platform_instantiation_artifact_id")
    if not isinstance(platform_aid, str):
        return projected
    registry = {item.get("artifact_id"): item for item in projected.get("artifacts", []) if isinstance(item, dict)}
    platform_artifact = registry.get(platform_aid)
    if not isinstance(platform_artifact, dict) or platform_artifact.get("storage") != "INLINE_JSON":
        return projected
    try:
        platform = base.strict_json_loads(base64.b64decode(platform_artifact["payload_b64"], validate=True))
    except Exception:
        return projected
    formation_ids = platform.get("formation_evidence_artifact_ids") if isinstance(platform, dict) else None
    if not isinstance(formation_ids, list):
        return projected
    for evidence_id in formation_ids:
        evidence_artifact = registry.get(evidence_id)
        if not isinstance(evidence_artifact, dict) or evidence_artifact.get("storage") != "INLINE_JSON":
            continue
        try:
            evidence = base.strict_json_loads(base64.b64decode(evidence_artifact["payload_b64"], validate=True))
        except Exception:
            continue
        observations = evidence.get("observations") if isinstance(evidence, dict) else None
        if isinstance(observations, list):
            evidence["observations"] = [
                {"predicate_id": member["predicate_id"], "value": member["value"], "reproducible": member["reproducible"]}
                if isinstance(member, dict) and {"predicate_id", "value", "reproducible"}.issubset(member)
                else member
                for member in observations
            ]
            raw = base.canonical_bytes(evidence)
            evidence_artifact["payload_b64"] = base64.b64encode(raw).decode()
            evidence_artifact["sha256"] = hashlib.sha256(raw).hexdigest()
            evidence_artifact["byte_length"] = len(raw)
    registry = {item.get("artifact_id"): item for item in projected.get("artifacts", []) if isinstance(item, dict)}
    for item in projected.get("artifacts", []):
        if not isinstance(item, dict):
            continue
        for parent in item.get("parents", []):
            parent_id = parent.get("artifact_id") if isinstance(parent, dict) else None
            if parent_id in registry:
                parent["sha256"] = registry[parent_id]["sha256"]
    return projected


def _refuse(result: dict[str, Any], diagnostics: list[str]) -> dict[str, Any]:
    result = copy.deepcopy(result)
    custody = result.setdefault("custody", {})
    custody["disposition"] = "REFUSE"
    custody["errors"] = sorted(set(custody.get("errors", []) + diagnostics))
    result["accepted"] = False
    result["actual_platform_present"] = False
    result["semantic_diagnostics"] = sorted(set(result.get("semantic_diagnostics", []) + diagnostics))
    result["scientific_gates"] = {gate: "UNSCOREABLE" for gate in base.GATES}
    result["candidate_milestones"] = {claim: "UNSCOREABLE" for claim in base.CLAIMS}
    result["authoritative_milestones"] = {claim: "UNSCOREABLE" for claim in base.CLAIMS}
    result["authoritative_proof_outputs"] = {claim: "NO_PROOF_OUTPUT" for claim in base.CLAIMS}
    result["formation_observation_repair2"] = {"valid": False, "diagnostics": diagnostics}
    return result


def evaluate_instance(instance: Any, workspace: Path | None = None) -> dict[str, Any]:
    workspace = (workspace or Path.cwd()).resolve()
    if not isinstance(instance, dict):
        return REPAIR1_EVALUATE(instance, workspace)
    valid, diagnostics = validate_actual_formation_observations(instance, workspace)
    projected = _project_for_repair1(instance)
    result = REPAIR1_EVALUATE(projected, workspace)
    if not valid:
        return _refuse(result, diagnostics)
    result["formation_observation_repair2"] = {"valid": True, "diagnostics": []}
    return result


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
            "custody": {"disposition": "REFUSE", "errors": [f"REPAIR2_ERROR:{type(exc).__name__}:{exc}"], "missing_artifact_ids": []},
            "actual_platform_present": False,
            "authoritative_proof_outputs": {claim: "NO_PROOF_OUTPUT" for claim in base.CLAIMS},
        }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if result.get("accepted") else 2


if __name__ == "__main__":
    raise SystemExit(main())
