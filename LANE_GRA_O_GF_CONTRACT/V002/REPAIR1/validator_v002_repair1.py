#!/usr/bin/env python3
"""Narrow Repair1 overlay for V002 PLATFORM_INSTANTIATION semantics.

All sealed V002 bytes remain unchanged.  This module calls the sealed evaluator and
then refuses any platform binding that is not a fully joined, non-synthetic,
pre-access platform packet.  `actual_platform_present` is derived only from that
conjunction; a label or boolean can never establish it.
"""

from __future__ import annotations

import argparse
import base64
import json
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


PARENT = Path(__file__).resolve().parent.parent
if str(PARENT) not in sys.path:
    sys.path.insert(0, str(PARENT))

import validator_v002 as base  # noqa: E402


REPAIR_ID = "GRA-O-GF-CONTRACT-V002-REPAIR1"
OLD_PLATFORM_KEYS = set(base.PAYLOAD_KEYS["PLATFORM_INSTANTIATION"])
PLATFORM_KEYS = {
    "platform_id", "package_id", "surface_ids", "synthetic_only",
    "platform_map_artifact_ids", "formation_evidence_artifact_ids",
    "gamma_process_artifact_ids", "seed_definition_artifact_id",
    "dilation_artifact_id", "transport_artifact_id",
    "source_custody_artifact_ids", "freeze_time",
}
PLATFORM_MAP_KEYS = {
    "platform_id", "package_id", "surface_id", "process_artifact_id",
    "seed_definition_artifact_id", "dilation_artifact_id",
    "transport_artifact_id", "formation_evidence_artifact_ids",
    "source_custody_artifact_ids",
}
SOURCE_CUSTODY_KEYS = {
    "platform_id", "package_id", "source_artifact_ids", "acquisition_ids",
    "specimen_ids", "independent_unit_ids", "provenance", "license_id",
    "synthetic_test",
}

BASE_EVALUATE = base.evaluate_instance


def configure_repair_types() -> None:
    base.PAYLOAD_KEYS["PLATFORM_INSTANTIATION"] = set(PLATFORM_KEYS)
    base.PAYLOAD_KEYS["PLATFORM_MAP"] = set(PLATFORM_MAP_KEYS)
    base.PAYLOAD_KEYS["PLATFORM_SOURCE_CUSTODY"] = set(SOURCE_CUSTODY_KEYS)
    base.KIND_ROLES["PLATFORM_MAP"] = {"DEVELOPMENT"}
    base.KIND_ROLES["PLATFORM_SOURCE_CUSTODY"] = {"DEVELOPMENT", "EXTERNAL_CALIBRATION"}


configure_repair_types()


def _payload(aid: Any, kind: str, registry: dict[str, dict[str, Any]], payloads: dict[str, Any], diagnostics: list[str]) -> dict[str, Any] | None:
    if not isinstance(aid, str) or aid not in registry:
        diagnostics.append(f"PLATFORM_ARTIFACT_UNRESOLVED:{aid}")
        return None
    if registry[aid].get("kind") != kind:
        diagnostics.append(f"PLATFORM_ARTIFACT_KIND_INVALID:{aid}:{kind}")
        return None
    value = payloads.get(aid)
    if not isinstance(value, dict):
        diagnostics.append(f"PLATFORM_ARTIFACT_PAYLOAD_INVALID:{aid}")
        return None
    return value


def _nonempty_unique_strings(value: Any) -> bool:
    return isinstance(value, list) and bool(value) and len(set(value)) == len(value) and all(isinstance(x, str) and bool(x.strip()) for x in value)


def _time(value: Any) -> datetime | None:
    if not isinstance(value, str) or not value:
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def _source_is_nonsynthetic(aid: str, registry: dict[str, dict[str, Any]], payloads: dict[str, Any]) -> bool:
    artifact = registry.get(aid, {})
    payload = payloads.get(aid)
    if artifact.get("kind") != "SOURCE_DATA" or artifact.get("role") not in {"DEVELOPMENT", "VALIDATION", "EXTERNAL_CALIBRATION", "GENERIC"} or not isinstance(payload, dict):
        return False
    try:
        raw = base64.b64decode(payload["content_b64"], validate=True)
    except Exception:
        return False
    label = str(payload.get("source_id", "")).lower()
    return "synthetic" not in label and b"synthetic" not in raw.lower()


def validate_platform_packet(instance: dict[str, Any], workspace: Path, sealed_result: dict[str, Any]) -> tuple[bool, list[str]]:
    diagnostics: list[str] = []
    bindings = instance.get("bindings")
    if not isinstance(bindings, dict):
        return False, ["PLATFORM_BINDINGS_INVALID"]
    platform_aid = bindings.get("platform_instantiation_artifact_id")
    if platform_aid is None:
        return False, ["ACTUAL_PLATFORM_ABSENT"]

    registry, payloads, artifact_errors, _ = base.validate_artifacts(instance, workspace)
    platform = _payload(platform_aid, "PLATFORM_INSTANTIATION", registry, payloads, diagnostics)
    if artifact_errors:
        diagnostics.append("PLATFORM_PACKET_CUSTODY_INVALID")
    if not isinstance(platform, dict) or set(platform) != PLATFORM_KEYS:
        return False, sorted(set(diagnostics + ["PLATFORM_INSTANCE_SCHEMA_INVALID"]))

    platform_id = platform.get("platform_id")
    if not isinstance(platform_id, str) or not re.fullmatch(r"[A-Z0-9][A-Z0-9._:-]{2,127}", platform_id):
        diagnostics.append("PLATFORM_ID_EMPTY_OR_INVALID")
    if platform.get("package_id") != instance.get("package_id"):
        diagnostics.append("PLATFORM_CROSS_PACKAGE_BINDING")
    if instance.get("package_mode") != "SCIENTIFIC" or platform.get("synthetic_only") is not False:
        diagnostics.append("PLATFORM_NOT_ACTUAL_SCIENTIFIC")

    surfaces = platform.get("surface_ids")
    maps = platform.get("platform_map_artifact_ids")
    formation_ids = platform.get("formation_evidence_artifact_ids")
    process_ids = platform.get("gamma_process_artifact_ids")
    custody_ids = platform.get("source_custody_artifact_ids")
    for label, value in (
        ("SURFACES", surfaces), ("MAPS", maps), ("FORMATION", formation_ids),
        ("PROCESSES", process_ids), ("SOURCE_CUSTODY", custody_ids),
    ):
        if not _nonempty_unique_strings(value):
            diagnostics.append(f"PLATFORM_{label}_EMPTY_OR_DUPLICATE")

    expected_process_ids = bindings.get("gamma_process_artifact_ids")
    if set(process_ids or []) != set(expected_process_ids or []):
        diagnostics.append("PLATFORM_PROCESS_BINDING_MISMATCH")
    process_surfaces: dict[str, str] = {}
    for process_id in process_ids or []:
        process = _payload(process_id, "PROCESS_REPRESENTATION", registry, payloads, diagnostics)
        if isinstance(process, dict):
            surface = process.get("surface_id")
            if not isinstance(surface, str) or not surface:
                diagnostics.append(f"PLATFORM_PROCESS_SURFACE_INVALID:{process_id}")
            elif surface in process_surfaces:
                diagnostics.append(f"PLATFORM_SURFACE_DUPLICATE_PROCESS:{surface}")
            else:
                process_surfaces[surface] = process_id
    if set(surfaces or []) != set(process_surfaces):
        diagnostics.append("PLATFORM_SURFACE_UNKNOWN_OR_INCOMPLETE")

    if platform.get("seed_definition_artifact_id") != bindings.get("seed_definition_artifact_id"):
        diagnostics.append("PLATFORM_SEED_BINDING_MISMATCH")
    if platform.get("dilation_artifact_id") != bindings.get("dilation_artifact_id"):
        diagnostics.append("PLATFORM_DILATION_BINDING_MISMATCH")
    if platform.get("transport_artifact_id") != bindings.get("transport_artifact_id"):
        diagnostics.append("PLATFORM_TRANSPORT_BINDING_MISMATCH")

    transport = _payload(platform.get("transport_artifact_id"), "CALIBRATED_TRANSPORT", registry, payloads, diagnostics)
    if isinstance(transport, dict) and transport.get("surface_id") not in set(surfaces or []):
        diagnostics.append("PLATFORM_TRANSPORT_SURFACE_MISMATCH")
    dilation = _payload(platform.get("dilation_artifact_id"), "FULL_DILATION", registry, payloads, diagnostics)
    if isinstance(dilation, dict) and (dilation.get("surface_id") not in set(surfaces or []) or not dilation.get("measurement_artifact_ids") or not dilation.get("exchange_ledger")):
        diagnostics.append("PLATFORM_DILATION_IDENTITY_INCOMPLETE")

    custody_source_ids: set[str] = set()
    for custody_id in custody_ids or []:
        custody = _payload(custody_id, "PLATFORM_SOURCE_CUSTODY", registry, payloads, diagnostics)
        if not isinstance(custody, dict):
            continue
        if custody.get("platform_id") != platform_id or custody.get("package_id") != instance.get("package_id"):
            diagnostics.append(f"PLATFORM_SOURCE_CUSTODY_SCOPE_MISMATCH:{custody_id}")
        if custody.get("synthetic_test") is not False or not all(
            _nonempty_unique_strings(custody.get(key))
            for key in ("source_artifact_ids", "acquisition_ids", "specimen_ids", "independent_unit_ids")
        ) or not custody.get("provenance") or not custody.get("license_id"):
            diagnostics.append(f"PLATFORM_SOURCE_CUSTODY_INCOMPLETE:{custody_id}")
        for source_id in custody.get("source_artifact_ids") or []:
            custody_source_ids.add(source_id)
            if not _source_is_nonsynthetic(source_id, registry, payloads):
                diagnostics.append(f"PLATFORM_SOURCE_NOT_NONSYNTHETIC:{source_id}")
    if not custody_source_ids:
        diagnostics.append("PLATFORM_NONSYNTHETIC_SOURCE_CUSTODY_EMPTY")

    for evidence_id in formation_ids or []:
        evidence = _payload(evidence_id, "GATE_EVIDENCE", registry, payloads, diagnostics)
        if not isinstance(evidence, dict):
            continue
        evidence_sources = set(evidence.get("source_artifact_ids") or [])
        observations = evidence.get("observations")
        if evidence.get("gate_id") != "FORM.ALLOW0" or evidence.get("reproducible") is not True or not isinstance(observations, list) or not observations or any(item.get("value") is not True or item.get("reproducible") is not True for item in observations if isinstance(item, dict)) or not evidence_sources or not evidence_sources.issubset(custody_source_ids):
            diagnostics.append(f"PLATFORM_FORMATION_EVIDENCE_INVALID:{evidence_id}")

    map_surfaces: dict[str, str] = {}
    for map_id in maps or []:
        platform_map = _payload(map_id, "PLATFORM_MAP", registry, payloads, diagnostics)
        if not isinstance(platform_map, dict):
            continue
        surface = platform_map.get("surface_id")
        if surface in map_surfaces:
            diagnostics.append(f"PLATFORM_MAP_SURFACE_DUPLICATE:{surface}")
        else:
            map_surfaces[surface] = map_id
        expected = {
            "platform_id": platform_id,
            "package_id": instance.get("package_id"),
            "process_artifact_id": process_surfaces.get(surface),
            "seed_definition_artifact_id": bindings.get("seed_definition_artifact_id"),
            "dilation_artifact_id": bindings.get("dilation_artifact_id"),
            "transport_artifact_id": bindings.get("transport_artifact_id"),
        }
        for key, value in expected.items():
            if platform_map.get(key) != value:
                diagnostics.append(f"PLATFORM_MAP_JOIN_MISMATCH:{map_id}:{key}")
        if set(platform_map.get("formation_evidence_artifact_ids") or []) != set(formation_ids or []) or set(platform_map.get("source_custody_artifact_ids") or []) != set(custody_ids or []):
            diagnostics.append(f"PLATFORM_MAP_EVIDENCE_JOIN_MISMATCH:{map_id}")
    if set(map_surfaces) != set(surfaces or []):
        diagnostics.append("PLATFORM_MAP_SURFACE_COVERAGE_INCOMPLETE")

    parent_map = {
        aid: {parent.get("artifact_id") for parent in artifact.get("parents", []) if isinstance(parent, dict)}
        for aid, artifact in registry.items()
    }
    bound_roots = set(process_ids or []) | set(formation_ids or []) | {
        platform.get("seed_definition_artifact_id"), platform.get("dilation_artifact_id"),
        platform.get("transport_artifact_id"),
    }
    required_source_leaves: set[str] = set()
    stack = [aid for aid in bound_roots if isinstance(aid, str)]
    visited: set[str] = set()
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        if registry.get(current, {}).get("kind") == "SOURCE_DATA":
            required_source_leaves.add(current)
        stack.extend(parent_map.get(current, ()))
    if not required_source_leaves.issubset(custody_source_ids):
        diagnostics.append("PLATFORM_BOUND_SOURCE_CUSTODY_INCOMPLETE")

    freeze = platform.get("freeze_time")
    freeze_time = _time(freeze)
    access = _payload(bindings.get("access_log_artifact_id"), "ACCESS_LOG", registry, payloads, diagnostics)
    response_time = _time(transport.get("response_access_time")) if isinstance(transport, dict) else None
    if freeze_time is None or response_time is None or freeze_time >= response_time:
        diagnostics.append("PLATFORM_FREEZE_NOT_PRE_RESPONSE")
    if isinstance(access, dict):
        package_freeze = _time(access.get("freeze_time"))
        if freeze_time is None or package_freeze is None or freeze_time > package_freeze:
            diagnostics.append("PLATFORM_FREEZE_AFTER_PACKAGE_FREEZE")
        validation_id = None
        validation = payloads.get(bindings.get("validation_root_artifact_id"))
        if isinstance(validation, dict):
            validation_id = validation.get("dataset_id")
        for event in access.get("events") or []:
            if event.get("dataset_id") == validation_id and event.get("event_type") == "NUMERIC_ACCESS":
                event_time = _time(event.get("time"))
                if freeze_time is None or event_time is None or freeze_time >= event_time:
                    diagnostics.append("PLATFORM_FREEZE_AFTER_VALIDATION_ACCESS")

    required_gate_passes = {
        "FORM.ALLOW0", "GAMMA.PROCESS", "SEED.PHI",
        "SCALE.CALIBRATED_TRANSPORT", "ANCESTRY.FULL_PATH",
    }
    gates = sealed_result.get("scientific_gates", {})
    if sealed_result.get("accepted") is not True or sealed_result.get("custody", {}).get("disposition") != "QUALIFIED":
        diagnostics.append("PLATFORM_PACKAGE_CUSTODY_NOT_QUALIFIED")
    if any(gates.get(gate) != "PASS" for gate in required_gate_passes):
        diagnostics.append("PLATFORM_REQUIRED_SCIENCE_IDENTITY_NOT_PASSING")

    return not diagnostics, sorted(set(diagnostics))


def _refuse_platform(result: dict[str, Any], diagnostics: list[str]) -> dict[str, Any]:
    result = json.loads(json.dumps(result))
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
    result["platform_repair1"] = {"valid": False, "diagnostics": diagnostics}
    return result


def evaluate_instance(instance: Any, workspace: Path | None = None) -> dict[str, Any]:
    configure_repair_types()
    workspace = (workspace or Path.cwd()).resolve()
    result = BASE_EVALUATE(instance, workspace)
    if not isinstance(instance, dict):
        return result
    bindings = instance.get("bindings")
    if not isinstance(bindings, dict) or bindings.get("platform_instantiation_artifact_id") is None:
        result["actual_platform_present"] = False
        result["platform_repair1"] = {"valid": False, "diagnostics": ["ACTUAL_PLATFORM_ABSENT"]}
        return result
    valid, diagnostics = validate_platform_packet(instance, workspace, result)
    if not valid:
        return _refuse_platform(result, diagnostics)
    result["actual_platform_present"] = bool(instance.get("package_mode") == "SCIENTIFIC")
    result["platform_repair1"] = {"valid": True, "diagnostics": []}
    return result


def evaluate_legacy_exact(instance: dict[str, Any], workspace: Path | None = None) -> dict[str, Any]:
    """Reproduce only the sealed V002 semantics for the exact historical witness."""
    current = set(base.PAYLOAD_KEYS["PLATFORM_INSTANTIATION"])
    try:
        base.PAYLOAD_KEYS["PLATFORM_INSTANTIATION"] = set(OLD_PLATFORM_KEYS)
        return BASE_EVALUATE(instance, (workspace or Path.cwd()).resolve())
    finally:
        base.PAYLOAD_KEYS["PLATFORM_INSTANTIATION"] = current


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
            "custody": {"disposition": "REFUSE", "errors": [f"REPAIR1_ERROR:{type(exc).__name__}:{exc}"], "missing_artifact_ids": []},
            "authoritative_proof_outputs": {claim: "NO_PROOF_OUTPUT" for claim in base.CLAIMS},
            "actual_platform_present": False,
        }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if result.get("accepted") else 2


if __name__ == "__main__":
    raise SystemExit(main())
