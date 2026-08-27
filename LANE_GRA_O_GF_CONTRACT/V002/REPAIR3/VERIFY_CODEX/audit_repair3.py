#!/usr/bin/env python3
"""Terminal, bounded, read-only audit of the sealed V002 Repair3 overlay."""

from __future__ import annotations

import base64
import copy
import hashlib
import inspect
import json
import sys
from pathlib import Path
from typing import Any, Callable


HERE = Path(__file__).resolve().parent
REPAIR3_DIR = HERE.parent
V002_DIR = REPAIR3_DIR.parent
REPAIR2_DIR = V002_DIR / "REPAIR2"
REPAIR1_DIR = V002_DIR / "REPAIR1"
LANE_DIR = V002_DIR.parent
ROOT = LANE_DIR.parent
for path in (REPAIR3_DIR, REPAIR2_DIR, REPAIR1_DIR, V002_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import mutation_suite_v002_repair1 as repair1_suite  # noqa: E402
import mutation_suite_v002_repair2 as repair2_suite  # noqa: E402
import mutation_suite_v002_repair3 as repair3_suite  # noqa: E402
import synthetic_positive_fixture_v002 as fixture  # noqa: E402
import validator_v002 as base  # noqa: E402
import validator_v002_repair1 as repair1  # noqa: E402
import validator_v002_repair2 as repair2  # noqa: E402
import validator_v002_repair3 as repair3  # noqa: E402


MANIFESTS = (
    LANE_DIR / "MANIFEST.sha256",
    LANE_DIR / "VERIFY_CODEX" / "MANIFEST.sha256",
    V002_DIR / "MANIFEST_V002.sha256",
    REPAIR1_DIR / "MANIFEST_V002_REPAIR1.sha256",
    REPAIR2_DIR / "MANIFEST_V002_REPAIR2.sha256",
    REPAIR3_DIR / "MANIFEST_V002_REPAIR3.sha256",
)
PROMOTED = ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2")
FORMATION_ID = repair2_suite.FORMATION_ID
SOURCE_ID = "REAL.PLATFORM.SOURCE"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(path: Path) -> dict[str, Any]:
    failures: list[str] = []
    entries = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = ROOT / relative.strip()
        entries += 1
        if not target.is_file():
            failures.append(f"MISSING:{relative.strip()}")
        elif digest(target) != expected:
            failures.append(f"HASH_MISMATCH:{relative.strip()}")
    return {
        "path": str(path.relative_to(ROOT)),
        "manifest_sha256": digest(path),
        "entries": entries,
        "pass": not failures,
        "failures": failures,
    }


def canonical(value: Any) -> bytes:
    return base.canonical_bytes(value)


def formation(candidate: dict[str, Any]) -> dict[str, Any]:
    return next(item for item in candidate["artifacts"] if item["artifact_id"] == FORMATION_ID)


def no_proof(result: dict[str, Any]) -> bool:
    outputs = result.get("authoritative_proof_outputs")
    return (
        isinstance(outputs, dict)
        and set(outputs) == set(base.CLAIMS)
        and all(outputs[claim] == "NO_PROOF_OUTPUT" for claim in base.CLAIMS)
    )


def promoted(result: dict[str, Any]) -> bool:
    outputs = result.get("authoritative_proof_outputs", {})
    return all(outputs.get(claim) == "PASSES_DECLARED_DOMAIN" for claim in PROMOTED)


def refused(result: dict[str, Any], expected: str) -> bool:
    return (
        result.get("accepted") is False
        and result.get("actual_platform_present") is False
        and expected in json.dumps(result, sort_keys=True, allow_nan=False)
        and no_proof(result)
    )


def mutate_parent_hash(candidate: dict[str, Any]) -> None:
    parent = next(item for item in formation(candidate)["parents"] if item["artifact_id"] == SOURCE_ID)
    parent["sha256"] = "0" * 64


def mutate_wrong_child_hash(candidate: dict[str, Any]) -> None:
    formation(candidate)["sha256"] = "f" * 64


def mutate_unknown_parent(candidate: dict[str, Any]) -> None:
    parent = next(item for item in formation(candidate)["parents"] if item["artifact_id"] == SOURCE_ID)
    parent.update({"artifact_id": "UNKNOWN.PARENT", "sha256": "0" * 64})


def mutate_self_cycle(candidate: dict[str, Any]) -> None:
    child = formation(candidate)
    child["parents"].append({"artifact_id": FORMATION_ID, "sha256": child["sha256"]})


def mutate_orphan(candidate: dict[str, Any]) -> None:
    content = b"independent QA orphan payload"
    payload = {
        "source_id": "VERIFY.REPAIR3.ORPHAN",
        "content_b64": base64.b64encode(content).decode(),
        "content_sha256": hashlib.sha256(content).hexdigest(),
    }
    raw = canonical(payload)
    candidate["artifacts"].append({
        "artifact_id": "VERIFY.REPAIR3.ORPHAN",
        "kind": "SOURCE_DATA",
        "role": "GENERIC",
        "sha256": hashlib.sha256(raw).hexdigest(),
        "byte_length": len(raw),
        "storage": "INLINE_JSON",
        "locator": None,
        "payload_b64": base64.b64encode(raw).decode(),
        "parents": [],
    })


def mutate_scalar_observation(candidate: dict[str, Any]) -> None:
    repair2_suite.change_observations(candidate, [0])


def mutate_mixed_observation(candidate: dict[str, Any]) -> None:
    repair2_suite.change_observations(candidate, repair2_suite.typed_observations() + [0])


def run_case(name: str, mutation: Callable[[dict[str, Any]], None], expected: str) -> dict[str, Any]:
    candidate = repair2_suite.typed_platform_candidate()
    mutation(candidate)
    before = canonical(candidate)
    result = repair3.evaluate_instance(candidate, ROOT)
    passed = refused(result, expected) and canonical(candidate) == before
    return {
        "case": name,
        "expected": expected,
        "accepted": result.get("accepted"),
        "actual_platform_present": result.get("actual_platform_present"),
        "no_proof_output": no_proof(result),
        "input_byte_identity_preserved": canonical(candidate) == before,
        "pass": passed,
    }


def audit_identity_order() -> dict[str, Any]:
    candidate = repair2_suite.typed_platform_candidate()
    before = canonical(candidate)
    candidate_identity = id(candidate)
    events: list[dict[str, Any]] = []
    observation_read_only = False

    original_repair1 = repair3.REPAIR1_EVALUATE_ORIGINAL
    original_observation = repair3.OBSERVATION_VALIDATE_ORIGINAL
    original_projection = repair2._project_for_repair1
    original_repair2_evaluate = repair2.evaluate_instance

    def repair1_spy(instance: Any, workspace: Path) -> dict[str, Any]:
        events.append({"call": "repair1", "same_original_object": id(instance) == candidate_identity})
        return original_repair1(instance, workspace)

    def observation_spy(instance: dict[str, Any], workspace: Path) -> tuple[bool, list[str]]:
        nonlocal observation_read_only
        events.append({"call": "observation", "same_original_object": id(instance) == candidate_identity})
        local_before = canonical(instance)
        value = original_observation(instance, workspace)
        observation_read_only = canonical(instance) == local_before
        return value

    def forbidden_projection(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("Repair2 projection/evaluator invoked by Repair3")

    repair3.REPAIR1_EVALUATE_ORIGINAL = repair1_spy
    repair3.OBSERVATION_VALIDATE_ORIGINAL = observation_spy
    repair2._project_for_repair1 = forbidden_projection
    repair2.evaluate_instance = forbidden_projection
    try:
        result = repair3.evaluate_instance(candidate, ROOT)
    finally:
        repair3.REPAIR1_EVALUATE_ORIGINAL = original_repair1
        repair3.OBSERVATION_VALIDATE_ORIGINAL = original_observation
        repair2._project_for_repair1 = original_projection
        repair2.evaluate_instance = original_repair2_evaluate

    call_order = [event["call"] for event in events]
    same_object = all(event["same_original_object"] for event in events)
    immutable = canonical(candidate) == before
    passed = (
        call_order == ["repair1", "observation"]
        and same_object
        and observation_read_only
        and immutable
        and result.get("accepted") is True
        and result.get("actual_platform_present") is True
        and result.get("formation_observation_repair3", {}).get("upstream_repair1_accepted") is True
    )
    return {
        "call_order": call_order,
        "repair1_call_count": call_order.count("repair1"),
        "observation_call_count": call_order.count("observation"),
        "same_original_object": same_object,
        "observation_read_only": observation_read_only,
        "input_byte_identity_preserved": immutable,
        "repair2_projection_or_evaluator_invoked": False,
        "pass": passed,
    }


def audit_exact_laundering() -> dict[str, Any]:
    candidate = repair2_suite.typed_platform_candidate()
    mutate_parent_hash(candidate)
    before = canonical(candidate)
    repair1_result = repair1.evaluate_instance(copy.deepcopy(candidate), ROOT)
    repair2_result = repair2.evaluate_instance(copy.deepcopy(candidate), ROOT)
    repair3_result = repair3.evaluate_instance(candidate, ROOT)
    repair1_refused = refused(repair1_result, "ARTIFACT_PARENT_HASH_MISMATCH")
    repair2_laundered = (
        repair2_result.get("accepted") is True
        and repair2_result.get("actual_platform_present") is True
        and promoted(repair2_result)
    )
    repair3_refused = refused(repair3_result, "ARTIFACT_PARENT_HASH_MISMATCH")
    observation_could_not_promote = (
        repair3_result.get("formation_observation_repair3", {}).get("valid") is True
        and repair3_result.get("formation_observation_repair3", {}).get("upstream_repair1_accepted") is False
    )
    immutable = canonical(candidate) == before
    return {
        "repair1_refused_original": repair1_refused,
        "repair2_laundering_reproduced": repair2_laundered,
        "repair3_refused_original": repair3_refused,
        "valid_observation_could_not_promote_repair1_refusal": observation_could_not_promote,
        "repair2_actual_platform_present": repair2_result.get("actual_platform_present"),
        "repair3_actual_platform_present": repair3_result.get("actual_platform_present"),
        "repair3_no_proof_output": no_proof(repair3_result),
        "input_byte_identity_preserved": immutable,
        "pass": repair1_refused and repair2_laundered and repair3_refused and observation_could_not_promote and immutable,
    }


def audit_typed_positive() -> dict[str, Any]:
    candidate = repair2_suite.typed_platform_candidate()
    before = canonical(candidate)
    result = repair3.evaluate_instance(candidate, ROOT)
    immutable = canonical(candidate) == before
    passed = (
        result.get("accepted") is True
        and result.get("actual_platform_present") is True
        and result.get("formation_observation_repair3", {}).get("valid") is True
        and promoted(result)
        and immutable
    )
    return {
        "accepted": result.get("accepted"),
        "actual_platform_present": result.get("actual_platform_present"),
        "typed_observation_valid": result.get("formation_observation_repair3", {}).get("valid"),
        "gf0_through_ge2_promoted": promoted(result),
        "input_byte_identity_preserved": immutable,
        "pass": passed,
    }


def audit_absent_platform() -> dict[str, Any]:
    candidate = fixture.build_fixture()
    before = canonical(candidate)
    result = repair3.evaluate_instance(candidate, ROOT)
    passed = (
        result.get("actual_platform_present") is False
        and no_proof(result)
        and canonical(candidate) == before
    )
    return {
        "accepted": result.get("accepted"),
        "actual_platform_present": result.get("actual_platform_present"),
        "all_authoritative_outputs_no_proof": no_proof(result),
        "input_byte_identity_preserved": canonical(candidate) == before,
        "pass": passed,
    }


def audit_nongeometric_boundary() -> dict[str, Any]:
    field_names = set(base.BINDING_KEYS) | set().union(*base.PAYLOAD_KEYS.values())
    normalized_fields = {name.lower() for name in field_names}
    forbidden_physical_fields = {
        "t_lab", "stress_energy", "stress-energy", "spacetime_metric",
        "gravitational_metric", "source_geometry", "geometry",
    }
    gf_source = inspect.getsource(base.evaluate_instance)
    milestone_source = inspect.getsource(base.derive_milestones)
    gf0_lines = [line.strip() for line in gf_source.splitlines() if line.strip().startswith(("gf0_candidate =", "gf0_authoritative ="))]
    gf1_lines = [
        line.strip()
        for line in milestone_source.splitlines()
        if line.strip().startswith('milestones["GF1"] =')
    ]
    no_forbidden_fields = not (normalized_fields & forbidden_physical_fields)
    gf0_has_only_early_ancestry = (
        len(gf0_lines) == 2
        and all("GRAVITY_GATES" not in line and "GC0" not in line for line in gf0_lines)
        and all(token in gf0_lines[0] for token in ("gamma_state", "seed_state", "transport_state", "ancestry_state"))
    )
    gf1_has_no_gravity_gate = (
        len(gf1_lines) == 1
        and "GRAVITY_GATES" not in gf1_lines[0]
        and "GC0" not in gf1_lines[0]
        and all(gate in gf1_lines[0] for gate in (
            "FORM.ALLOW0", "GAMMA.PROCESS", "SEED.PHI",
            "SCALE.CALIBRATED_TRANSPORT", "ANCESTRY.FULL_PATH",
        ))
    )
    metric_occurrence_is_transport_math = (
        "metric_residual" in inspect.getsource(base.validate_transport)
        and "metric" not in normalized_fields
    )
    passed = no_forbidden_fields and gf0_has_only_early_ancestry and gf1_has_no_gravity_gate and metric_occurrence_is_transport_math
    return {
        "required_t_lab_stress_energy_geometry_fields": sorted(normalized_fields & forbidden_physical_fields),
        "gf0_expressions": gf0_lines,
        "gf1_expression": gf1_lines,
        "gravity_gate_required_by_gf0_or_gf1": False if gf0_has_only_early_ancestry and gf1_has_no_gravity_gate else True,
        "metric_identifier_is_only_internal_transport_compatibility_not_candidate_field": metric_occurrence_is_transport_math,
        "pass": passed,
    }


def main() -> int:
    manifests = [verify_manifest(path) for path in MANIFESTS]
    suite_run_1 = repair3_suite.run_suite()
    suite_run_2 = repair3_suite.run_suite()
    suites_identical = canonical(suite_run_1) == canonical(suite_run_2)
    suite_checks = {
        "run_1": suite_run_1.get("overall_result"),
        "run_2": suite_run_2.get("overall_result"),
        "case_count_each": [suite_run_1.get("case_count"), suite_run_2.get("case_count")],
        "failed_count_each": [suite_run_1.get("failed_count"), suite_run_2.get("failed_count")],
        "result_payload_sha256_each": [suite_run_1.get("result_payload_sha256"), suite_run_2.get("result_payload_sha256")],
        "byte_identical_payloads": suites_identical,
        "pass": (
            suite_run_1.get("overall_result") == suite_run_2.get("overall_result") == "PASS"
            and suite_run_1.get("case_count") == suite_run_2.get("case_count") == 59
            and suite_run_1.get("failed_count") == suite_run_2.get("failed_count") == 0
            and suites_identical
        ),
    }

    cases = [
        run_case("wrong_child_digest", mutate_wrong_child_hash, "ARTIFACT_HASH_MISMATCH:REAL.PLATFORM.FORMATION"),
        run_case("unknown_parent", mutate_unknown_parent, "ARTIFACT_PARENT_UNRESOLVED"),
        run_case("self_cycle", mutate_self_cycle, "ARTIFACT_DAG_SELF_LOOP"),
        run_case("orphan", mutate_orphan, "ARTIFACT_ORPHAN:VERIFY.REPAIR3.ORPHAN"),
        run_case("scalar_observation_[0]", mutate_scalar_observation, "FORMATION_OBSERVATION_MEMBER_NOT_OBJECT"),
        run_case("mixed_object_scalar_observations", mutate_mixed_observation, "FORMATION_OBSERVATION_MEMBER_NOT_OBJECT"),
    ]
    identity = audit_identity_order()
    laundering = audit_exact_laundering()
    typed_positive = audit_typed_positive()
    absent_platform = audit_absent_platform()
    nongeometric = audit_nongeometric_boundary()

    sections = {
        "manifests": all(row["pass"] for row in manifests),
        "suite_twice": suite_checks["pass"],
        "identity_order_monotonicity": identity["pass"],
        "repair2_laundering_repair3_refusal": laundering["pass"],
        "independent_cases": all(row["pass"] for row in cases),
        "typed_positive": typed_positive["pass"],
        "absent_platform_no_proof": absent_platform["pass"],
        "nongeometric_boundary": nongeometric["pass"],
    }
    passed = all(sections.values())
    result = {
        "audit_id": "GRA-O-GF-CONTRACT-V002-REPAIR3-VERIFY-CODEX",
        "scope": "TERMINAL_BOUNDED_INDEPENDENT_QA",
        "verdict": "PASS" if passed else "FAIL",
        "scientific_result": "NO_PROOF_OUTPUT",
        "section_passes": sections,
        "manifest_verification": manifests,
        "repair3_suite_twice": suite_checks,
        "immutable_original_call_order": identity,
        "exact_repair2_parent_hash_laundering": laundering,
        "independent_refusal_cases": cases,
        "typed_positive_control": typed_positive,
        "absent_actual_platform_control": absent_platform,
        "nongeometric_gf0_gf1_boundary": nongeometric,
        "first_decisive_defect": None if passed else next(name for name, ok in sections.items() if not ok),
    }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
