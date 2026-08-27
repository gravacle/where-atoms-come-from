#!/usr/bin/env python3
"""Executable semantic validator/evaluator for GRA-O GF contract V002.

JSON Schema is deliberately not the authority for cross-object joins.  This module
recomputes artifact custody, process/seed/transport/data/ancestry predicates, the
four-valued gate results, the fixed GF0--UGE conjunctions, and product reproduction.
No input status, outcome, proof result, or product result is trusted.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Iterable


CONTRACT_ID = "GRA-O-GF0-UGE-CONTRACT-V002"
TASK_ID = "GRA-O-GF-CONTRACT"
CLAIMS = ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2", "UGE")
CORE_GATES = (
    "FORM.ALLOW0",
    "GAMMA.PROCESS",
    "SEED.PHI",
    "SCALE.CALIBRATED_TRANSPORT",
    "ANCESTRY.FULL_PATH",
    "VALIDATION.HELD_OUT",
    "CROSS_SURFACE.UNCHANGED",
    "EVIDENCE.JOINED",
)
GRAVITY_GATES = tuple(f"GC{i:02d}" for i in range(1, 15))
FINAL_GATES = ("COVERAGE.CONSTRUCTIVE", "WORLD.REPRODUCTION")
GATES = CORE_GATES + GRAVITY_GATES + FINAL_GATES
ANCESTRY_CHECKS = (
    "IDENTITY", "SUPPORT", "UNITS", "CAUSALITY", "SYMMETRY", "CONSERVATION", "INTERVENTION"
)
PRINCIPAL_DECISION_PATH = "LANE_GRA_S_JOINT_SEED_DECISION/PRINCIPAL_DECISION.md"
PRINCIPAL_DECISION_SHA256 = "962de48aa3ced22887200933d1dd397951a2dd25cc8bf9ea335893fc888ef82d"

FORBIDDEN_DERIVED_KEYS = {
    "status", "outcome", "custody_status", "scientific_outcome", "proof_outcome",
    "product_outcome", "adequate_sensitivity", "adequate_power", "qualified",
    "milestone_outcomes", "proof_outputs", "scientific_weight",
}

TOP_KEYS = {
    "contract_id", "task_id", "package_id", "package_mode", "producer_id", "principal_id",
    "claim_ids", "gate_ids", "horizon_complete_claim", "bindings", "artifacts", "gate_inputs",
    "gravity_applicability", "milestone_domains",
}
BINDING_KEYS = {
    "framework_artifact_ids", "principal_decision_artifact_id", "gamma_process_artifact_ids",
    "gamma_identification_artifact_ids", "seed_definition_artifact_id", "dilation_artifact_id",
    "platform_instantiation_artifact_id", "transport_artifact_id", "development_root_artifact_id",
    "validation_root_artifact_id", "access_log_artifact_id", "ancestry_artifact_id",
    "joined_evidence_artifact_id", "coverage_artifact_id", "public_release_artifact_id",
    "public_dataset_manifest_artifact_id", "executor_key_artifact_id",
    "execution_output_artifact_id", "execution_report_artifact_id",
    "execution_signature_artifact_id",
}
ARTIFACT_KEYS = {
    "artifact_id", "kind", "role", "sha256", "byte_length", "storage", "locator",
    "payload_b64", "parents",
}

KIND_ROLES = {
    "FRAMEWORK": {"GENERIC"},
    "PRINCIPAL_DECISION": {"THEORY"},
    "PROCESS_REPRESENTATION": {"EXTERNAL_CALIBRATION"},
    "PROCESS_IDENTIFICATION": {"EXTERNAL_CALIBRATION"},
    "SEED_DEFINITION": {"THEORY"},
    "FULL_DILATION": {"EXTERNAL_CALIBRATION", "THEORY"},
    "PLATFORM_INSTANTIATION": {"DEVELOPMENT"},
    "CALIBRATED_TRANSPORT": {"DEVELOPMENT"},
    "DATASET_ROOT": {"DEVELOPMENT", "VALIDATION"},
    "ACCESS_LOG": {"GENERIC"},
    "ANCESTRY_GRAPH": {"DEVELOPMENT"},
    "GATE_RULE": {"THEORY"},
    "GATE_EVIDENCE": {"DEVELOPMENT", "VALIDATION"},
    "POWER_CERTIFICATE": {"THEORY"},
    "JOINED_EVIDENCE": {"VALIDATION"},
    "APPLICABILITY_JUSTIFICATION": {"THEORY"},
    "DOMAIN_WITNESS": {"THEORY"},
    "COVERAGE_THEOREM": {"THEORY"},
    "PUBLIC_DATASET_MANIFEST": {"PRODUCT"},
    "PUBLIC_RELEASE": {"PRODUCT"},
    "EXECUTOR_KEY": {"PRODUCT"},
    "EXECUTION_OUTPUT": {"PRODUCT"},
    "EXECUTION_REPORT": {"PRODUCT"},
    "EXECUTION_SIGNATURE": {"PRODUCT"},
    "SOURCE_DATA": {"DEVELOPMENT", "VALIDATION", "EXTERNAL_CALIBRATION", "GENERIC", "PRODUCT"},
}

PAYLOAD_KEYS: dict[str, set[str]] = {
    "FRAMEWORK": {"document_id", "source_path", "source_sha256"},
    "PRINCIPAL_DECISION": {"decision_id", "decision", "source_path", "source_sha256"},
    "PROCESS_REPRESENTATION": {"surface_id", "carrier_class", "representation_type", "equation", "state_space", "native_units", "support", "bath", "clock", "measured_channels"},
    "PROCESS_IDENTIFICATION": {"surface_id", "process_artifact_id", "measurement_artifact_ids", "method", "projection_map_artifact_ids", "independent_process_identification"},
    "SEED_DEFINITION": {"seed_id", "construction", "blocks", "pre_response_inputs", "operator_type", "units", "support", "roles", "conserved_quantity", "derivation_stage", "principal_decision_artifact_id", "dilation_artifact_id"},
    "FULL_DILATION": {"surface_id", "system", "environment", "exchange_ledger", "physical_quantities", "measurement_artifact_ids"},
    "PLATFORM_INSTANTIATION": {"platform_id", "surface_ids", "synthetic_only", "platform_map_artifact_ids", "freeze_time"},
    "CALIBRATED_TRANSPORT": {"transport_id", "surface_id", "scales", "edges", "covariance", "composition", "refinements", "ancestry_component", "freeze_time", "response_access_time"},
    "DATASET_ROOT": {"dataset_id", "dataset_role", "raw_root_sha256", "source_ids", "acquisition_ids", "specimen_ids", "independent_unit_ids", "event_ids", "outcome_content_ids", "source_artifact_ids", "acquired_time"},
    "ACCESS_LOG": {"freeze_time", "events", "signed_by", "signature_artifact_id"},
    "ANCESTRY_GRAPH": {"seed_node_id", "endpoint_node_id", "nodes", "arrows"},
    "GATE_RULE": {"gate_id", "predicate_ids", "decision_region", "freeze_time", "rule_version"},
    "GATE_EVIDENCE": {"gate_id", "observations", "source_artifact_ids", "taxonomy_match", "reproducible"},
    "POWER_CERTIFICATE": {"gate_id", "alpha", "beta", "false_pass_upper", "achieved_power_lower", "monte_carlo_error", "mc_tolerance", "independent_units", "effect_size_core"},
    "JOINED_EVIDENCE": {"formation_artifact_ids", "growth_artifact_ids", "response_artifact_ids", "join_keys"},
    "APPLICABILITY_JUSTIFICATION": {"characteristic_id", "scope", "physical_reason", "evidence_artifact_ids"},
    "DOMAIN_WITNESS": {"claim_id", "member_ids", "quantifier", "coverage_artifact_id"},
    "COVERAGE_THEOREM": {"theorem_id", "surface_classes", "constructive_map", "preservation_obligations", "checker_artifact_ids"},
    "PUBLIC_DATASET_MANIFEST": {"release_id", "dataset_artifact_ids", "dataset_sha256s", "public_locators", "license_ids"},
    "PUBLIC_RELEASE": {"release_id", "contract_id", "package_id", "claim_ids", "dataset_manifest_artifact_id", "dataset_manifest_sha256", "validator_sha256", "expected_output_artifact_id", "expected_output_sha256"},
    "EXECUTOR_KEY": {"executor_id", "n_hex", "e"},
    "EXECUTION_OUTPUT": {"release_id", "result_digest", "claim_results"},
    "EXECUTION_REPORT": {"release_artifact_id", "output_artifact_id", "executor_id", "environment_sha256", "no_private_logic", "mismatch_count"},
    "EXECUTION_SIGNATURE": {"key_artifact_id", "release_artifact_id", "output_artifact_id", "report_artifact_id", "message_b64", "signature_hex"},
    "SOURCE_DATA": {"source_id", "content_b64", "content_sha256"},
}


class DuplicateKeyError(ValueError):
    pass


def _pairs_no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in pairs:
        if key in out:
            raise DuplicateKeyError(f"duplicate JSON key: {key}")
        out[key] = value
    return out


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number: {value}")


def strict_json_loads(raw: str | bytes) -> Any:
    return json.loads(raw, object_pairs_hook=_pairs_no_duplicates, parse_constant=_reject_constant)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8")


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def exact_keys(value: Any, keys: set[str], label: str, errors: list[str]) -> bool:
    if not isinstance(value, dict):
        errors.append(f"{label}:OBJECT_REQUIRED")
        return False
    unknown = sorted(set(value) - keys)
    missing = sorted(keys - set(value))
    if unknown:
        errors.append(f"{label}:UNKNOWN_KEYS:{','.join(unknown)}")
    if missing:
        errors.append(f"{label}:MISSING_KEYS:{','.join(missing)}")
    return not unknown and not missing


def finite_tree(value: Any, path: str, errors: list[str]) -> None:
    if isinstance(value, float) and not math.isfinite(value):
        errors.append(f"NONFINITE_NUMBER:{path}")
    elif isinstance(value, dict):
        for key, child in value.items():
            if key in FORBIDDEN_DERIVED_KEYS:
                errors.append(f"UNTRUSTED_DERIVED_FIELD:{path}.{key}")
            finite_tree(child, f"{path}.{key}", errors)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            finite_tree(child, f"{path}[{index}]", errors)


def decode_artifact(artifact: dict[str, Any], workspace: Path, errors: list[str]) -> tuple[bytes | None, Any | None]:
    aid = artifact.get("artifact_id", "?")
    if artifact.get("storage") == "INLINE_JSON":
        if artifact.get("locator") is not None or not isinstance(artifact.get("payload_b64"), str):
            errors.append(f"ARTIFACT_STORAGE_INVALID:{aid}")
            return None, None
        try:
            raw = base64.b64decode(artifact["payload_b64"], validate=True)
        except Exception:
            errors.append(f"ARTIFACT_BASE64_INVALID:{aid}")
            return None, None
        try:
            payload = strict_json_loads(raw)
        except Exception as exc:
            errors.append(f"ARTIFACT_JSON_INVALID:{aid}:{type(exc).__name__}")
            return raw, None
    elif artifact.get("storage") == "WORKSPACE_FILE":
        if artifact.get("payload_b64") is not None or not isinstance(artifact.get("locator"), str):
            errors.append(f"ARTIFACT_STORAGE_INVALID:{aid}")
            return None, None
        path = Path(artifact["locator"])
        if path.is_absolute() or ".." in path.parts:
            errors.append(f"ARTIFACT_LOCATOR_UNSAFE:{aid}")
            return None, None
        path = workspace / path
        if not path.is_file():
            errors.append(f"ARTIFACT_FILE_MISSING:{aid}")
            return None, None
        raw = path.read_bytes()
        payload = None
    else:
        errors.append(f"ARTIFACT_STORAGE_INVALID:{aid}")
        return None, None
    if artifact.get("sha256") != sha256_bytes(raw):
        errors.append(f"ARTIFACT_HASH_MISMATCH:{aid}")
    if artifact.get("byte_length") != len(raw):
        errors.append(f"ARTIFACT_LENGTH_MISMATCH:{aid}")
    return raw, payload


def referenced_artifact_ids(value: Any) -> set[str]:
    refs: set[str] = set()
    if isinstance(value, dict):
        for key, child in value.items():
            if key.endswith("_artifact_id") and isinstance(child, str):
                refs.add(child)
            elif key.endswith("_artifact_ids") and isinstance(child, list):
                refs.update(x for x in child if isinstance(x, str))
            refs.update(referenced_artifact_ids(child))
    elif isinstance(value, list):
        for child in value:
            refs.update(referenced_artifact_ids(child))
    return refs


def validate_artifacts(instance: dict[str, Any], workspace: Path) -> tuple[dict[str, dict[str, Any]], dict[str, Any], list[str], list[str]]:
    errors: list[str] = []
    missing: list[str] = []
    registry: dict[str, dict[str, Any]] = {}
    payloads: dict[str, Any] = {}
    digests: dict[str, str] = {}
    artifacts = instance.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        return {}, {}, ["ARTIFACT_REGISTRY_EMPTY"], []
    for index, artifact in enumerate(artifacts):
        if not exact_keys(artifact, ARTIFACT_KEYS, f"artifact[{index}]", errors):
            continue
        aid = artifact.get("artifact_id")
        if not isinstance(aid, str) or not re.fullmatch(r"[A-Z0-9][A-Z0-9._:-]{1,159}", aid):
            errors.append(f"ARTIFACT_ID_INVALID:{aid}")
            continue
        if aid in registry:
            errors.append(f"ARTIFACT_ID_DUPLICATE:{aid}")
            continue
        registry[aid] = artifact
        digest = artifact.get("sha256")
        if digest in digests:
            errors.append(f"ARTIFACT_DIGEST_REUSE:{digests[digest]}:{aid}")
        elif isinstance(digest, str):
            digests[digest] = aid
        kind = artifact.get("kind")
        role = artifact.get("role")
        if kind not in KIND_ROLES:
            errors.append(f"ARTIFACT_KIND_UNKNOWN:{aid}:{kind}")
        elif role not in KIND_ROLES[kind]:
            errors.append(f"KIND_ROLE_MISMATCH:{aid}:{kind}:{role}")
        _, payload = decode_artifact(artifact, workspace, errors)
        if artifact.get("storage") == "INLINE_JSON" and payload is not None:
            payloads[aid] = payload
            expected = PAYLOAD_KEYS.get(kind)
            if expected is None:
                errors.append(f"PAYLOAD_KIND_UNVALIDATED:{aid}:{kind}")
            else:
                exact_keys(payload, expected, f"payload:{aid}", errors)
                finite_tree(payload, f"payload:{aid}", errors)
                if kind == "SOURCE_DATA" and isinstance(payload, dict):
                    try:
                        source_bytes = base64.b64decode(payload["content_b64"], validate=True)
                    except Exception:
                        errors.append(f"SOURCE_DATA_BASE64_INVALID:{aid}")
                    else:
                        if payload.get("content_sha256") != sha256_bytes(source_bytes):
                            errors.append(f"SOURCE_DATA_HASH_MISMATCH:{aid}")

    graph: dict[str, set[str]] = {aid: set() for aid in registry}
    for aid, artifact in registry.items():
        parents = artifact.get("parents")
        if not isinstance(parents, list):
            errors.append(f"ARTIFACT_PARENTS_INVALID:{aid}")
            continue
        seen: set[str] = set()
        for parent in parents:
            if not exact_keys(parent, {"artifact_id", "sha256"}, f"parent:{aid}", errors):
                continue
            pid = parent.get("artifact_id")
            if pid in seen:
                errors.append(f"ARTIFACT_PARENT_DUPLICATE:{aid}:{pid}")
            seen.add(pid)
            if pid not in registry:
                missing.append(pid)
                errors.append(f"ARTIFACT_PARENT_UNRESOLVED:{aid}:{pid}")
                continue
            if parent.get("sha256") != registry[pid].get("sha256"):
                errors.append(f"ARTIFACT_PARENT_HASH_MISMATCH:{aid}:{pid}")
            graph[aid].add(pid)
        if aid in graph[aid]:
            errors.append(f"ARTIFACT_DAG_SELF_LOOP:{aid}")
        if aid in payloads:
            refs = referenced_artifact_ids(payloads[aid])
            if refs != seen:
                errors.append(f"ARTIFACT_PARENT_JOIN_MISMATCH:{aid}")

    color: dict[str, int] = {aid: 0 for aid in registry}
    def visit(aid: str) -> None:
        if color[aid] == 1:
            errors.append(f"ARTIFACT_DAG_CYCLE:{aid}")
            return
        if color[aid] == 2:
            return
        color[aid] = 1
        for pid in graph[aid]:
            visit(pid)
        color[aid] = 2
    for aid in registry:
        visit(aid)

    roots = referenced_artifact_ids({
        "bindings": instance.get("bindings", {}),
        "gate_inputs": instance.get("gate_inputs", []),
        "gravity_applicability": instance.get("gravity_applicability", []),
        "milestone_domains": instance.get("milestone_domains", []),
    })
    roots.update(instance.get("bindings", {}).get("framework_artifact_ids", []) if isinstance(instance.get("bindings"), dict) else [])
    reachable: set[str] = set()
    stack = [root for root in roots if root in registry]
    while stack:
        aid = stack.pop()
        if aid in reachable:
            continue
        reachable.add(aid)
        stack.extend(graph.get(aid, ()))
    for root in sorted(roots):
        if root not in registry:
            missing.append(root)
            errors.append(f"BINDING_UNRESOLVED:{root}")
    for aid in sorted(set(registry) - reachable):
        errors.append(f"ARTIFACT_ORPHAN:{aid}")
    return registry, payloads, sorted(set(errors)), sorted(set(missing))


# ---- small deterministic real-matrix kernel ---------------------------------

Matrix = list[list[float]]
TOL = 1e-8


def shape(a: Matrix) -> tuple[int, int]:
    if not isinstance(a, list) or not a or not all(isinstance(row, list) and row for row in a):
        raise ValueError("matrix must be nonempty")
    width = len(a[0])
    if any(len(row) != width for row in a):
        raise ValueError("ragged matrix")
    for row in a:
        for value in row:
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise ValueError("matrix entry must be finite real")
    return len(a), width


def eye(n: int) -> Matrix:
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def zeros(rows: int, cols: int) -> Matrix:
    return [[0.0 for _ in range(cols)] for _ in range(rows)]


def transpose(a: Matrix) -> Matrix:
    rows, cols = shape(a)
    return [[float(a[i][j]) for i in range(rows)] for j in range(cols)]


def matmul(a: Matrix, b: Matrix) -> Matrix:
    ar, ac = shape(a)
    br, bc = shape(b)
    if ac != br:
        raise ValueError(f"matrix shape mismatch {ar}x{ac} @ {br}x{bc}")
    return [[sum(float(a[i][k]) * float(b[k][j]) for k in range(ac)) for j in range(bc)] for i in range(ar)]


def matsub(a: Matrix, b: Matrix) -> Matrix:
    ar, ac = shape(a)
    if shape(b) != (ar, ac):
        raise ValueError("matrix subtraction shape mismatch")
    return [[float(a[i][j]) - float(b[i][j]) for j in range(ac)] for i in range(ar)]


def diag(values: list[float]) -> Matrix:
    return [[float(values[i]) if i == j else 0.0 for j in range(len(values))] for i in range(len(values))]


def frob(a: Matrix) -> float:
    shape(a)
    return math.sqrt(sum(float(x) * float(x) for row in a for x in row))


def close_matrix(a: Matrix, b: Matrix, tolerance: float = TOL) -> bool:
    try:
        return frob(matsub(a, b)) <= tolerance
    except (TypeError, ValueError):
        return False


def matrix_rank(a: Matrix, tolerance: float = TOL) -> int:
    work = [[float(x) for x in row] for row in a]
    rows, cols = shape(work)
    rank = 0
    for col in range(cols):
        pivot = max(range(rank, rows), key=lambda r: abs(work[r][col]), default=rank)
        if rank >= rows or abs(work[pivot][col]) <= tolerance:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][col]
        work[rank] = [x / scale for x in work[rank]]
        for row in range(rows):
            if row == rank:
                continue
            factor = work[row][col]
            work[row] = [work[row][j] - factor * work[rank][j] for j in range(cols)]
        rank += 1
        if rank == rows:
            break
    return rank


def jacobi_eigen_symmetric(a: Matrix, tolerance: float = 1e-12, iterations: int = 200) -> tuple[list[float], Matrix]:
    n, m = shape(a)
    if n != m or not close_matrix(a, transpose(a), 1e-9):
        raise ValueError("matrix is not symmetric")
    work = [[float(x) for x in row] for row in a]
    vectors = eye(n)
    for _ in range(iterations):
        p, q, largest = 0, 0, 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(work[i][j]) > largest:
                    p, q, largest = i, j, abs(work[i][j])
        if largest <= tolerance:
            break
        phi = 0.5 * math.atan2(2.0 * work[p][q], work[q][q] - work[p][p])
        c, s = math.cos(phi), math.sin(phi)
        rotation = eye(n)
        rotation[p][p], rotation[q][q] = c, c
        rotation[p][q], rotation[q][p] = -s, s
        work = matmul(transpose(rotation), matmul(work, rotation))
        vectors = matmul(vectors, rotation)
    values = [work[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: values[i], reverse=True)
    values = [values[i] for i in order]
    vectors = [[vectors[r][i] for i in order] for r in range(n)]
    return values, vectors


def spectral_parts(w: Matrix, threshold: float) -> tuple[Matrix, Matrix, Matrix, Matrix, list[float]]:
    values, vectors = jacobi_eigen_symmetric(w)
    if any(value < -1e-8 for value in values):
        raise ValueError("W not positive semidefinite")
    kept = [value > threshold for value in values]
    p = matmul(vectors, matmul(diag([1.0 if flag else 0.0 for flag in kept]), transpose(vectors)))
    n = matsub(eye(len(values)), p)
    sqrt_w = matmul(vectors, matmul(diag([math.sqrt(max(value, 0.0)) for value in values]), transpose(vectors)))
    pinv_sqrt = matmul(vectors, matmul(diag([1.0 / math.sqrt(value) if flag else 0.0 for value, flag in zip(values, kept)]), transpose(vectors)))
    return p, n, sqrt_w, pinv_sqrt, values


def adequate_power(power: Any) -> bool:
    keys = {"alpha", "beta", "false_pass_upper", "achieved_power_lower", "monte_carlo_error", "mc_tolerance", "independent_units", "effect_size_core"}
    if not isinstance(power, dict) or set(power) != keys:
        return False
    numeric = [power[k] for k in ("alpha", "beta", "false_pass_upper", "achieved_power_lower", "monte_carlo_error", "mc_tolerance", "effect_size_core")]
    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(float(x)) for x in numeric):
        return False
    return (
        0 < power["alpha"] < 1 and 0 < power["beta"] < 1
        and power["false_pass_upper"] <= power["alpha"]
        and power["achieved_power_lower"] >= 1.0 - power["beta"]
        and power["monte_carlo_error"] <= power["mc_tolerance"]
        and isinstance(power["independent_units"], int) and power["independent_units"] >= 2
        and power["effect_size_core"] > 0
    )


SCALE_KEYS = {"scale_id", "order", "physical_scale", "units", "W", "eigenvectors", "eigenvalues", "eigenvalue_intervals", "rank_threshold", "max_condition", "P", "N", "sqrt_W", "pinv_sqrt_W", "reference", "calibration_artifact_id"}
EDGE_KEYS = {"edge_id", "source", "target", "F", "Btilde", "causal_superchannel_artifact_id", "null_leakage_norm", "transport_loss_norm", "transport_rank", "condition_number", "null_margin", "power"}
REFINEMENT_KEYS = {"control_id", "type", "same_event", "frozen_before_response", "new_independent_mode", "transformation_artifact_id", "paired_event_artifact_id", "scales", "edges", "covariance", "power"}
REF_SCALE_KEYS = {"scale_id", "V", "R", "W_refined", "P_refined", "N_refined", "U", "base_stream_artifact_id", "refined_stream_artifact_id"}
REF_EDGE_KEYS = {"edge_id", "F_refined", "Btilde_refined"}


def validate_transport(payload: Any) -> tuple[str, list[str]]:
    """Return PASS/FAIL/UNSCOREABLE/REFUSE and deterministic diagnostics."""
    diagnostics: list[str] = []
    if not isinstance(payload, dict) or set(payload) != PAYLOAD_KEYS["CALIBRATED_TRANSPORT"]:
        return "REFUSE", ["TRANSPORT_PAYLOAD_CLOSED_SCHEMA"]
    if payload.get("ancestry_component") != "BTILDE":
        diagnostics.append("ANCESTRY_USES_UNCALIBRATED_TRANSPORT")
    if not isinstance(payload.get("freeze_time"), str) or not isinstance(payload.get("response_access_time"), str) or payload["freeze_time"] >= payload["response_access_time"]:
        diagnostics.append("TRANSPORT_FREEZE_CHRONOLOGY_VIOLATION")
    scales_raw = payload.get("scales")
    if not isinstance(scales_raw, list) or len(scales_raw) < 3:
        return "UNSCOREABLE", diagnostics + ["MULTISCALE_COVERAGE_INCOMPLETE"]
    scale_map: dict[str, dict[str, Any]] = {}
    derived: dict[str, dict[str, Any]] = {}
    try:
        for scale in scales_raw:
            if not isinstance(scale, dict) or set(scale) != SCALE_KEYS:
                raise ValueError("TRANSPORT_SCALE_CLOSED_SCHEMA")
            sid = scale["scale_id"]
            if sid in scale_map:
                raise ValueError("TRANSPORT_SCALE_ID_DUPLICATE")
            scale_map[sid] = scale
            dim, dim2 = shape(scale["W"])
            if dim != dim2 or len(scale["reference"]) != dim:
                raise ValueError("TRANSPORT_CALIBRATION_DIMENSION")
            if not isinstance(scale["rank_threshold"], (int, float)) or scale["rank_threshold"] <= 0:
                raise ValueError("TRANSPORT_RANK_THRESHOLD_INVALID")
            p, n, sqrt_w, pinv, values = spectral_parts(scale["W"], float(scale["rank_threshold"]))
            if not close_matrix(scale["W"], matmul(scale["eigenvectors"], matmul(diag(scale["eigenvalues"]), transpose(scale["eigenvectors"]))), 1e-7):
                raise ValueError("TRANSPORT_EIGENSYSTEM_MISMATCH")
            if any(abs(a - b) > 1e-7 for a, b in zip(values, scale["eigenvalues"])):
                raise ValueError("TRANSPORT_EIGENVALUES_MISMATCH")
            if not (close_matrix(scale["P"], p) and close_matrix(scale["N"], n) and close_matrix(scale["sqrt_W"], sqrt_w) and close_matrix(scale["pinv_sqrt_W"], pinv)):
                raise ValueError("TRANSPORT_PROJECTOR_RECOMPUTE_MISMATCH")
            intervals = scale["eigenvalue_intervals"]
            if len(intervals) != dim or any(not isinstance(x, list) or len(x) != 2 for x in intervals):
                raise ValueError("TRANSPORT_RANK_INTERVAL_INVALID")
            if any(float(lo) <= scale["rank_threshold"] <= float(hi) for lo, hi in intervals):
                diagnostics.append(f"UNSCOREABLE_RANK:{sid}")
            derived[sid] = {"P": p, "N": n, "sqrt": sqrt_w, "pinv": pinv, "rank": matrix_rank(p)}
        ordered = sorted(scales_raw, key=lambda x: x["order"])
        if [x["order"] for x in ordered] != list(range(len(ordered))) or any(ordered[i]["physical_scale"] >= ordered[i + 1]["physical_scale"] for i in range(len(ordered) - 1)):
            raise ValueError("SCALE_NOT_STRICTLY_ORDERED")
    except (KeyError, TypeError, ValueError) as exc:
        return "REFUSE", diagnostics + [str(exc)]

    edges_raw = payload.get("edges")
    edge_map: dict[str, dict[str, Any]] = {}
    try:
        if not isinstance(edges_raw, list):
            raise ValueError("TRANSPORT_EDGES_INVALID")
        for edge in edges_raw:
            if not isinstance(edge, dict) or set(edge) != EDGE_KEYS:
                raise ValueError("TRANSPORT_EDGE_CLOSED_SCHEMA")
            eid = edge["edge_id"]
            if eid in edge_map:
                raise ValueError("TRANSPORT_EDGE_ID_DUPLICATE")
            edge_map[eid] = edge
            source, target = edge["source"], edge["target"]
            if source not in scale_map or target not in scale_map or scale_map[source]["order"] >= scale_map[target]["order"]:
                raise ValueError("RAW_TRANSPORT_CAUSALITY_INVALID")
            expected = matmul(derived[target]["sqrt"], matmul(derived[target]["P"], matmul(edge["F"], matmul(derived[source]["P"], derived[source]["pinv"]))))
            if not close_matrix(edge["Btilde"], expected, 1e-7):
                raise ValueError("BTILDE_RECOMPUTE_MISMATCH")
            leakage = matmul(derived[target]["P"], matmul(edge["F"], derived[source]["N"]))
            loss = matmul(derived[target]["N"], matmul(edge["F"], derived[source]["P"]))
            if abs(frob(leakage) - edge["null_leakage_norm"]) > 1e-7 or abs(frob(loss) - edge["transport_loss_norm"]) > 1e-7:
                raise ValueError("TRANSPORT_NULL_OUTPUT_MISMATCH")
            rank = matrix_rank(expected)
            if edge["transport_rank"] != rank:
                raise ValueError("TRANSPORT_RANK_OUTPUT_MISMATCH")
            singular_squared, _ = jacobi_eigen_symmetric(matmul(transpose(expected), expected))
            singular = [math.sqrt(max(value, 0.0)) for value in singular_squared if value > TOL * TOL]
            computed_condition = (max(singular) / min(singular)) if singular else math.inf
            if abs(float(edge["condition_number"]) - computed_condition) > 1e-6:
                raise ValueError("TRANSPORT_CONDITION_OUTPUT_MISMATCH")
            if rank == 0:
                diagnostics.append("UNSCOREABLE_COMMON_SUPPORT")
            if not adequate_power(edge["power"]):
                diagnostics.append(f"UNSCOREABLE_UNDERPOWERED:{eid}")
            elif edge["null_leakage_norm"] > edge["null_margin"]:
                diagnostics.append(f"FAIL_NULL_LEAKAGE:{eid}")
            if not isinstance(edge["condition_number"], (int, float)) or not math.isfinite(edge["condition_number"]) or edge["condition_number"] > scale_map[source]["max_condition"]:
                diagnostics.append(f"UNSCOREABLE_CONDITION:{eid}")
        required_edges = {(ordered[0]["scale_id"], ordered[1]["scale_id"]), (ordered[1]["scale_id"], ordered[2]["scale_id"]), (ordered[0]["scale_id"], ordered[2]["scale_id"])}
        if not required_edges.issubset({(x["source"], x["target"]) for x in edges_raw}):
            diagnostics.append("MULTISCALE_COVERAGE_INCOMPLETE")
    except (KeyError, TypeError, ValueError) as exc:
        return "REFUSE", diagnostics + [str(exc)]

    covariance = payload.get("covariance")
    cov_keys = {"method", "uncertain_components", "shared_cross_scale", "recompute_sqrt_each_draw", "recompute_pinv_each_draw", "recompute_rank_each_draw", "rank_boundary_crossed", "common_support_stable"}
    if not isinstance(covariance, dict) or set(covariance) != cov_keys:
        diagnostics.append("UNSCOREABLE_TRANSPORT_COVARIANCE")
    else:
        required_uncertain = {"z", "F", "W", "reference", "rank", "shared_sources"}
        if set(covariance["uncertain_components"]) != required_uncertain or not all(covariance[k] is True for k in ("shared_cross_scale", "recompute_sqrt_each_draw", "recompute_pinv_each_draw", "recompute_rank_each_draw", "common_support_stable")):
            diagnostics.append("UNSCOREABLE_TRANSPORT_COVARIANCE")
        if covariance["rank_boundary_crossed"]:
            diagnostics.append("UNSCOREABLE_RANK")

    composition = payload.get("composition")
    comp_keys = {"source", "middle", "target", "direct_edge", "first_edge", "second_edge", "P_in_common", "P_out_common", "residual_norm", "residual_margin", "power"}
    try:
        if not isinstance(composition, dict) or set(composition) != comp_keys:
            raise ValueError("CALIBRATED_COMPOSITION_CLOSED_SCHEMA")
        direct = edge_map[composition["direct_edge"]]["Btilde"]
        sequential = matmul(edge_map[composition["second_edge"]]["Btilde"], edge_map[composition["first_edge"]]["Btilde"])
        if (
            (edge_map[composition["direct_edge"]]["source"], edge_map[composition["direct_edge"]]["target"]) != (composition["source"], composition["target"])
            or (edge_map[composition["first_edge"]]["source"], edge_map[composition["first_edge"]]["target"]) != (composition["source"], composition["middle"])
            or (edge_map[composition["second_edge"]]["source"], edge_map[composition["second_edge"]]["target"]) != (composition["middle"], composition["target"])
        ):
            raise ValueError("CALIBRATED_COMPOSITION_EDGE_JOIN_INVALID")
        residual = matmul(composition["P_out_common"], matmul(matsub(direct, sequential), composition["P_in_common"]))
        if matrix_rank(composition["P_in_common"]) == 0 or matrix_rank(composition["P_out_common"]) == 0:
            diagnostics.append("UNSCOREABLE_COMMON_SUPPORT")
        if abs(frob(residual) - composition["residual_norm"]) > 1e-7:
            raise ValueError("CALIBRATED_COMPOSITION_OUTPUT_MISMATCH")
        if not adequate_power(composition["power"]):
            diagnostics.append("UNSCOREABLE_COMPOSITION_UNDERPOWERED")
        elif composition["residual_norm"] > composition["residual_margin"]:
            diagnostics.append("FAIL_CALIBRATED_COMPOSITION")
    except (KeyError, TypeError, ValueError) as exc:
        return "REFUSE", diagnostics + [str(exc)]

    refinements = payload.get("refinements")
    if not isinstance(refinements, list) or len(refinements) != 3 or len({x.get("control_id") for x in refinements if isinstance(x, dict)}) != 3 or {x.get("type") for x in refinements if isinstance(x, dict)} != {"RELABEL", "SUBDIVIDE", "DUPLICATE"}:
        diagnostics.append("REFINEMENT_COVERAGE_INCOMPLETE")
    else:
        for control in refinements:
            try:
                if set(control) != REFINEMENT_KEYS:
                    raise ValueError("REFINEMENT_CLOSED_SCHEMA")
                cid = control["control_id"]
                if control["same_event"] is not True or control["frozen_before_response"] is not True:
                    diagnostics.append(f"UNSCOREABLE_REFINEMENT_IDENTITY:{cid}")
                if control["new_independent_mode"] is True:
                    diagnostics.append(f"UNSCOREABLE_NEW_INDEPENDENT_MODE:{cid}")
                control_power_ok = adequate_power(control["power"])
                if not control_power_ok:
                    diagnostics.append(f"UNSCOREABLE_REFINEMENT_UNDERPOWERED:{cid}")
                ref_scales = {x["scale_id"]: x for x in control["scales"] if isinstance(x, dict)}
                if set(ref_scales) != set(scale_map):
                    diagnostics.append(f"REFINEMENT_COVERAGE_INCOMPLETE:{cid}")
                    continue
                ref_derived: dict[str, dict[str, Any]] = {}
                for sid, rs in ref_scales.items():
                    if set(rs) != REF_SCALE_KEYS:
                        raise ValueError("REFINEMENT_SCALE_CLOSED_SCHEMA")
                    base = scale_map[sid]
                    dbase = derived[sid]
                    v, r = rs["V"], rs["R"]
                    vr, vc = shape(v)
                    if vc != shape(base["W"])[0] or matrix_rank(v) != vc:
                        raise ValueError("REFINEMENT_EMBEDDING_INVALID")
                    if not close_matrix(matmul(r, matmul(v, dbase["P"])), dbase["P"]):
                        raise ValueError("REFINEMENT_EMBEDDING_INVALID")
                    pr, nr, sr, pir, _ = spectral_parts(rs["W_refined"], float(base["rank_threshold"]))
                    expected_pr = matmul(v, matmul(dbase["P"], transpose(v)))
                    if matrix_rank(pr) != dbase["rank"] or not close_matrix(rs["P_refined"], pr) or not close_matrix(pr, expected_pr) or not close_matrix(rs["N_refined"], nr):
                        raise ValueError("REFINEMENT_RANK_PROJECTOR_INVALID")
                    if frob(matmul(pr, matmul(v, dbase["N"]))) > TOL or frob(matmul(sr, matmul(v, dbase["N"]))) > TOL:
                        raise ValueError("REFINEMENT_NULL_COMPATIBILITY_INVALID")
                    metric_residual = matsub(matmul(dbase["P"], matmul(transpose(v), matmul(rs["W_refined"], matmul(v, dbase["P"])))), matmul(dbase["P"], matmul(base["W"], dbase["P"])))
                    if frob(metric_residual) > TOL:
                        raise ValueError("REFINEMENT_METRIC_MISMATCH")
                    u_expected = matmul(sr, matmul(v, matmul(dbase["P"], dbase["pinv"])))
                    if not close_matrix(rs["U"], u_expected) or not close_matrix(matmul(transpose(rs["U"]), rs["U"]), dbase["P"]):
                        raise ValueError("REFINEMENT_U_DERIVATION_INVALID")
                    ref_derived[sid] = {"P": pr, "N": nr, "sqrt": sr, "pinv": pir, "U": rs["U"]}
                ref_edges = {x["edge_id"]: x for x in control["edges"] if isinstance(x, dict)}
                if set(ref_edges) != set(edge_map):
                    diagnostics.append(f"REFINEMENT_EDGE_COVERAGE_INCOMPLETE:{cid}")
                    continue
                for eid, re in ref_edges.items():
                    if set(re) != REF_EDGE_KEYS:
                        raise ValueError("REFINEMENT_EDGE_CLOSED_SCHEMA")
                    base_edge = edge_map[eid]
                    source, target = base_edge["source"], base_edge["target"]
                    vin, vout = ref_scales[source]["V"], ref_scales[target]["V"]
                    if not close_matrix(matmul(re["F_refined"], vin), matmul(vout, base_edge["F"])) and control_power_ok:
                        diagnostics.append(f"FAIL_RAW_REFINEMENT_NATURALITY:{cid}:{eid}")
                    b_expected = matmul(ref_derived[target]["sqrt"], matmul(ref_derived[target]["P"], matmul(re["F_refined"], matmul(ref_derived[source]["P"], ref_derived[source]["pinv"]))))
                    if not close_matrix(re["Btilde_refined"], b_expected):
                        raise ValueError("REFINED_BTILDE_RECOMPUTE_MISMATCH")
                    if not close_matrix(matmul(re["Btilde_refined"], ref_derived[source]["U"]), matmul(ref_derived[target]["U"], base_edge["Btilde"])) and control_power_ok:
                        diagnostics.append(f"FAIL_CALIBRATED_REFINEMENT_NATURALITY:{cid}:{eid}")
                rcov = control["covariance"]
                required_rcov = {"paired_events", "includes_cross_covariance", "pushforward_checked", "recompute_per_draw", "rank_boundary_crossed"}
                if not isinstance(rcov, dict) or set(rcov) != required_rcov or not all(rcov[k] is True for k in ("paired_events", "includes_cross_covariance", "pushforward_checked", "recompute_per_draw")):
                    diagnostics.append(f"UNSCOREABLE_REFINEMENT_COVARIANCE:{cid}")
                elif rcov["rank_boundary_crossed"]:
                    diagnostics.append(f"UNSCOREABLE_REFINEMENT_RANK:{cid}")
            except (KeyError, TypeError, ValueError) as exc:
                return "REFUSE", diagnostics + [f"{exc}:{control.get('control_id', '?')}"]

    if any(item.startswith("FAIL_") for item in diagnostics):
        return "FAIL", sorted(set(diagnostics))
    if diagnostics:
        return "UNSCOREABLE" if not any(item.startswith(("ANCESTRY_USES_", "TRANSPORT_FREEZE_")) for item in diagnostics) else "REFUSE", sorted(set(diagnostics))
    return "PASS", []


REQUIRED_FRAMEWORKS = {
    "HANDOFF_2026-08-22.md": "61c678f62e60816083e6b16bcd0d4a6176b579bc9c176aa4fd2b57472eb7eaa1",
    "LANE_GRA_K_DISCOVERY_FRAMEWORK/DERIVATION_AND_PROTOCOL.md": "56cb64594c564f2b24e68f160d541265ad0fb94535a0e4cb8655172803f23720",
    "LANE_GRA_L_CHARACTERISTICS_MATRIX/GRAVITY_CHARACTERISTICS_MATRIX.md": "a784a10cb462e27233b8bd5fa2afcbd8f00e4e383b81302d32f5947009df4308",
    "LANE_GRA_N_POWERED_RULES/CONTENT_ADDRESSED_ARTIFACT_CONTRACT.md": "ad77c16504d5f029e70e730e968e0db134d2e91c0db0d4d75ac43186a4de4022",
    "LANE_GRA_N_POWERED_RULES/PROTOCOL_AND_DECISION_RULES.md": "b46ae996850628a024e716d890e1f33dece28ecf7e05ce851e1f50194364d11a",
    "LANE_GRA_N_OPERATOR_GENEALOGY_SEARCH/PROTOCOL_FREEZE.md": "56993957cf9f4b545eed93a412c0bb94f82f07b9b1f0ac65f218f8c28916041e",
    "LANE_GRA_N_OPERATOR_GENEALOGY_SEARCH/FINDINGS.md": "efe1eba345e1eee0972129d14afb2a1b50143effe1a24f49c3ce947bd30070b4",
    "LANE_GRA_N_OPERATOR_GENEALOGY_SEARCH/DEVELOPMENT_SCALE_FLOW.md": "86106f5bbc1c463750068648585bf11f418e1ac23f5e15fc3da49e73422d4e56",
    "LANE_GRA_N_VALIDATION_PAIRS/PAIR_SELECTION_FREEZE.md": "0d755fadfb3264dec47533be1632ae55bf493fb732e9abae77b98581dece48a4",
    "LANE_GRA_N_VALIDATION_PAIRS/ASI_DEVELOPMENT_PREREGISTRATION.md": "3b9b8226b2aa210be7c098401cfa6f81d4c237995b3c90761b3488064c0f9dbe",
    "LANE_GRA_N_VALIDATION_PAIRS/SEARCH_REPORT.md": "358e20dd7dac41149b8e6f9b6613121629692fac54f7b8f8d87ec96b7c663956",
    "LANE_GRA_N_VALIDATION_PAIRS/BOUNDED_CLOSE_CHECKPOINT.md": "732dbff75f7d0b327b660a309437ffd121a0a3cc6bd3fc81159358253ae4a003",
    "LANE_GRA_O_GAMMA_OPERATOR_DECISION/DECISION_CARD_V004.md": "93774521323b62cbeb7489b0bae4100f0b44e4af8c86297e14b08715e431d17e",
}


def artifact_payload(aid: Any, expected_kind: str, registry: dict[str, dict[str, Any]], payloads: dict[str, Any], errors: list[str]) -> Any | None:
    if not isinstance(aid, str) or aid not in registry:
        errors.append(f"REQUIRED_ARTIFACT_UNRESOLVED:{aid}")
        return None
    if registry[aid].get("kind") != expected_kind:
        errors.append(f"ARTIFACT_BINDING_KIND_INVALID:{aid}:{expected_kind}")
        return None
    if aid not in payloads:
        errors.append(f"ARTIFACT_INLINE_PAYLOAD_REQUIRED:{aid}")
        return None
    return payloads[aid]


def validate_framework_bindings(bindings: dict[str, Any], registry: dict[str, dict[str, Any]], errors: list[str]) -> None:
    aids = bindings.get("framework_artifact_ids")
    if not isinstance(aids, list) or len(aids) != len(REQUIRED_FRAMEWORKS) or len(set(aids)) != len(aids):
        errors.append("FRAMEWORK_REGISTRY_INCOMPLETE")
        return
    observed: dict[str, str] = {}
    for aid in aids:
        artifact = registry.get(aid)
        if not artifact or artifact.get("kind") != "FRAMEWORK" or artifact.get("storage") != "WORKSPACE_FILE":
            errors.append(f"FRAMEWORK_BINDING_INVALID:{aid}")
            continue
        observed[str(artifact.get("locator"))] = str(artifact.get("sha256"))
    if observed != REQUIRED_FRAMEWORKS:
        errors.append("FRAMEWORK_IDENTITY_SET_MISMATCH")
    decision_id = bindings.get("principal_decision_artifact_id")
    decision = registry.get(decision_id)
    if not decision or decision.get("kind") != "PRINCIPAL_DECISION" or decision.get("storage") != "WORKSPACE_FILE":
        errors.append("SEED_AUTHORIZATION_MISSING")
    elif decision.get("locator") != PRINCIPAL_DECISION_PATH or decision.get("sha256") != PRINCIPAL_DECISION_SHA256:
        errors.append("SEED_AUTHORIZATION_IDENTITY_INVALID")


def validate_gamma_seed(bindings: dict[str, Any], registry: dict[str, dict[str, Any]], payloads: dict[str, Any]) -> tuple[str, str, list[str], list[dict[str, Any]]]:
    errors: list[str] = []
    process_ids = bindings.get("gamma_process_artifact_ids")
    identification_ids = bindings.get("gamma_identification_artifact_ids")
    processes: list[dict[str, Any]] = []
    if not isinstance(process_ids, list) or not process_ids:
        errors.append("GAMMA_PROCESS_UNRESOLVED")
    else:
        for aid in process_ids:
            payload = artifact_payload(aid, "PROCESS_REPRESENTATION", registry, payloads, errors)
            if isinstance(payload, dict):
                processes.append(payload)
                if payload.get("representation_type") not in {"HAMILTONIAN", "LIOUVILLIAN", "STOCHASTIC_GENERATOR", "PROCESS_TENSOR", "NON_MARKOVIAN_TRANSITION_LAW"}:
                    errors.append(f"GAMMA_PROCESS_KIND_INVALID:{aid}")
                for key in ("equation", "state_space", "support", "bath", "clock", "carrier_class"):
                    if not isinstance(payload.get(key), str) or not payload[key].strip():
                        errors.append(f"GAMMA_PROCESS_FIELD_MISSING:{aid}:{key}")
                for key in ("native_units", "measured_channels"):
                    if not isinstance(payload.get(key), list) or not payload[key]:
                        errors.append(f"GAMMA_PROCESS_FIELD_MISSING:{aid}:{key}")
    if not isinstance(identification_ids, list) or len(identification_ids) != len(process_ids or []):
        errors.append("GAMMA_IDENTIFICATION_MISSING")
    else:
        identified: set[str] = set()
        for aid in identification_ids:
            payload = artifact_payload(aid, "PROCESS_IDENTIFICATION", registry, payloads, errors)
            if isinstance(payload, dict):
                process_id = payload.get("process_artifact_id")
                identified.add(process_id)
                if process_id not in process_ids or payload.get("independent_process_identification") is not True or not payload.get("measurement_artifact_ids"):
                    errors.append(f"GAMMA_IDENTIFICATION_INVALID:{aid}")
        if identified != set(process_ids or []):
            errors.append("GAMMA_IDENTIFICATION_JOIN_MISMATCH")

    seed_id = bindings.get("seed_definition_artifact_id")
    seed = artifact_payload(seed_id, "SEED_DEFINITION", registry, payloads, errors)
    seed_errors: list[str] = []
    if isinstance(seed, dict):
        if seed.get("seed_id") != "JOINT_SEED" or seed.get("principal_decision_artifact_id") != bindings.get("principal_decision_artifact_id"):
            seed_errors.append("SEED_AUTHORIZATION_SCOPE_MISMATCH")
        blocks = seed.get("blocks")
        if not isinstance(blocks, list) or {x.get("family") for x in blocks if isinstance(x, dict)} != {"B", "C", "D"}:
            seed_errors.append("JOINT_SEED_BLOCKS_INCOMPLETE")
        else:
            for block in blocks:
                if set(block) != {"family", "operator_ids", "physical_type", "units", "nonredundancy_witness"} or not block["operator_ids"] or not block["units"] or not block["nonredundancy_witness"]:
                    seed_errors.append(f"JOINT_SEED_BLOCK_INVALID:{block.get('family', '?')}")
        if seed.get("derivation_stage") != "PRE_RESPONSE" or not seed.get("pre_response_inputs"):
            seed_errors.append("SEED_RESPONSE_DERIVED")
        roles = set(seed.get("roles") or [])
        physical = bool(roles & {"CURRENT", "STRESS", "FLUX", "WORK", "HEAT", "MOMENTUM", "ENTROPY_PRODUCTION"}) or seed.get("conserved_quantity") != "NONE"
        dilation_id = bindings.get("dilation_artifact_id")
        if physical or any(x.get("family") == "D" for x in blocks or []):
            dilation = artifact_payload(dilation_id, "FULL_DILATION", registry, payloads, seed_errors)
            if not isinstance(dilation, dict) or seed.get("dilation_artifact_id") != dilation_id or not dilation.get("exchange_ledger") or not dilation.get("measurement_artifact_ids"):
                seed_errors.append("PHYSICAL_CURRENT_DILATION_REQUIRED")
        elif seed.get("dilation_artifact_id") is not None:
            seed_errors.append("DILATION_BINDING_INCONSISTENT")
    errors.extend(seed_errors)
    gamma_state = "PASS" if not any(x.startswith("GAMMA_") for x in errors) else "UNSCOREABLE"
    seed_state = "PASS" if not seed_errors and isinstance(seed, dict) else "UNSCOREABLE"
    return gamma_state, seed_state, errors, processes


def validate_data_topology(bindings: dict[str, Any], registry: dict[str, dict[str, Any]], payloads: dict[str, Any]) -> tuple[str, list[str]]:
    diagnostics: list[str] = []
    development = artifact_payload(bindings.get("development_root_artifact_id"), "DATASET_ROOT", registry, payloads, diagnostics)
    validation = artifact_payload(bindings.get("validation_root_artifact_id"), "DATASET_ROOT", registry, payloads, diagnostics)
    access = artifact_payload(bindings.get("access_log_artifact_id"), "ACCESS_LOG", registry, payloads, diagnostics)
    if not all(isinstance(x, dict) for x in (development, validation, access)):
        return "UNSCOREABLE", diagnostics + ["DATA_TOPOLOGY_INCOMPLETE"]
    if development["dataset_role"] != "DEVELOPMENT" or validation["dataset_role"] != "VALIDATION":
        diagnostics.append("DATASET_ROLE_INVALID")
    axes = {
        "RAW_ROOT": ({development["raw_root_sha256"]}, {validation["raw_root_sha256"]}),
        "SOURCE": (set(development["source_ids"]), set(validation["source_ids"])),
        "ACQUISITION": (set(development["acquisition_ids"]), set(validation["acquisition_ids"])),
        "SPECIMEN": (set(development["specimen_ids"]), set(validation["specimen_ids"])),
        "INDEPENDENT_UNIT": (set(development["independent_unit_ids"]), set(validation["independent_unit_ids"])),
        "EVENT": (set(development["event_ids"]), set(validation["event_ids"])),
        "OUTCOME_CONTENT": (set(development["outcome_content_ids"]), set(validation["outcome_content_ids"])),
        "CONTENT_ANCESTRY": (set(development["source_artifact_ids"]), set(validation["source_artifact_ids"])),
    }
    for axis, (left, right) in axes.items():
        if not left or not right:
            diagnostics.append(f"DV_AXIS_EMPTY:{axis}")
        if left & right:
            diagnostics.append(f"DV_OVERLAP_{axis}")
    parent_map = {
        aid: {parent.get("artifact_id") for parent in artifact.get("parents", []) if isinstance(parent, dict)}
        for aid, artifact in registry.items()
    }
    def ancestors(aid: str) -> set[str]:
        found: set[str] = set()
        stack = list(parent_map.get(aid, ()))
        while stack:
            current = stack.pop()
            if current in found:
                continue
            found.add(current)
            stack.extend(parent_map.get(current, ()))
        return found
    d_sources = set(development["source_artifact_ids"])
    v_sources = set(validation["source_artifact_ids"])
    if any(d in ancestors(vs) or vs in ancestors(d) for d in d_sources for vs in v_sources):
        diagnostics.append("DV_CONTENT_ANCESTRY_OVERLAP")
    for aid in d_sources:
        if registry.get(aid, {}).get("role") != "DEVELOPMENT":
            diagnostics.append(f"DATASET_SOURCE_ROLE_INVALID:DEVELOPMENT:{aid}")
    for aid in v_sources:
        if registry.get(aid, {}).get("role") != "VALIDATION":
            diagnostics.append(f"DATASET_SOURCE_ROLE_INVALID:VALIDATION:{aid}")
    if access.get("signed_by") == "" or access.get("signature_artifact_id") != bindings.get("principal_decision_artifact_id"):
        diagnostics.append("ACCESS_LOG_AUTHENTICITY_INVALID")
    freeze = access.get("freeze_time")
    events = access.get("events")
    if not isinstance(events, list) or not events:
        diagnostics.append("ACCESS_LOG_MISSING")
    else:
        validation_numeric = False
        for event in events:
            if not isinstance(event, dict) or set(event) != {"dataset_id", "event_type", "time"}:
                diagnostics.append("ACCESS_LOG_EVENT_INVALID")
                continue
            if event["dataset_id"] == validation["dataset_id"] and event["event_type"] == "NUMERIC_ACCESS":
                validation_numeric = True
                if not isinstance(freeze, str) or event["time"] <= freeze:
                    diagnostics.append("VALIDATION_ACCESSED_BEFORE_FREEZE")
        if not validation_numeric:
            diagnostics.append("VALIDATION_NUMERIC_ACCESS_UNRECORDED")
    return ("PASS" if not diagnostics else "UNSCOREABLE"), sorted(set(diagnostics))


def validate_ancestry(payload: Any, transport_id: str) -> tuple[str, list[str]]:
    diagnostics: list[str] = []
    if not isinstance(payload, dict) or set(payload) != PAYLOAD_KEYS["ANCESTRY_GRAPH"]:
        return "REFUSE", ["ANCESTRY_CLOSED_SCHEMA"]
    nodes = payload.get("nodes")
    arrows = payload.get("arrows")
    if not isinstance(nodes, list) or not nodes or not isinstance(arrows, list) or not arrows:
        return "UNSCOREABLE", ["ANCESTRY_GRAPH_EMPTY"]
    node_map: dict[str, dict[str, Any]] = {}
    for node in nodes:
        if not isinstance(node, dict) or set(node) != {"node_id", "surface_id", "scale_id", "scale_order", "operator_type", "transport_artifact_id"}:
            return "REFUSE", ["ANCESTRY_NODE_CLOSED_SCHEMA"]
        if node["node_id"] in node_map:
            return "REFUSE", ["ANCESTRY_NODE_ID_DUPLICATE"]
        if node["transport_artifact_id"] != transport_id:
            return "REFUSE", ["ANCESTRY_TRANSPORT_JOIN_INVALID"]
        node_map[node["node_id"]] = node
    seed, endpoint = payload.get("seed_node_id"), payload.get("endpoint_node_id")
    if seed not in node_map or endpoint not in node_map:
        return "REFUSE", ["ANCESTRY_ENDPOINT_UNRESOLVED"]
    adjacency: dict[str, set[str]] = {node: set() for node in node_map}
    arrow_ids: set[str] = set()
    for arrow in arrows:
        if not isinstance(arrow, dict) or set(arrow) != {"arrow_id", "source", "target", "transport_artifact_id", "transport_component", "checks"}:
            return "REFUSE", ["ANCESTRY_ARROW_CLOSED_SCHEMA"]
        if arrow["arrow_id"] in arrow_ids:
            return "REFUSE", ["ANCESTRY_ARROW_ID_DUPLICATE"]
        arrow_ids.add(arrow["arrow_id"])
        source, target = arrow["source"], arrow["target"]
        if source not in node_map or target not in node_map:
            return "REFUSE", ["ANCESTRY_ARROW_ENDPOINT_UNRESOLVED"]
        if source == target:
            return "REFUSE", ["ANCESTRY_SELF_LOOP"]
        if node_map[source]["scale_order"] >= node_map[target]["scale_order"]:
            return "REFUSE", ["ANCESTRY_SCALE_ORDER_INVALID"]
        if arrow["transport_artifact_id"] != transport_id or arrow["transport_component"] != "BTILDE":
            return "REFUSE", ["ANCESTRY_USES_UNCALIBRATED_TRANSPORT"]
        checks = arrow["checks"]
        if not isinstance(checks, list) or {x.get("check_id") for x in checks if isinstance(x, dict)} != set(ANCESTRY_CHECKS):
            return "REFUSE", ["ANCESTRY_CHECK_REGISTRY_INVALID"]
        for check in checks:
            if set(check) != {"check_id", "observed", "threshold", "direction", "power"}:
                return "REFUSE", ["ANCESTRY_CHECK_CLOSED_SCHEMA"]
            if not adequate_power(check["power"]):
                diagnostics.append(f"UNSCOREABLE_ANCESTRY_UNDERPOWERED:{arrow['arrow_id']}:{check['check_id']}")
            direction = check["direction"]
            passed = check["observed"] >= check["threshold"] if direction == "GE" else check["observed"] <= check["threshold"] if direction == "LE" else False
            if not passed and adequate_power(check["power"]):
                diagnostics.append(f"FAIL_ANCESTRY_{check['check_id']}:{arrow['arrow_id']}")
        adjacency[source].add(target)
    color = {node: 0 for node in node_map}
    cycle = False
    def visit(node: str) -> None:
        nonlocal cycle
        if color[node] == 1:
            cycle = True
            return
        if color[node] == 2:
            return
        color[node] = 1
        for child in adjacency[node]:
            visit(child)
        color[node] = 2
    for node in node_map:
        visit(node)
    if cycle:
        return "REFUSE", ["ANCESTRY_CYCLE"]
    reachable = {seed}
    stack = [seed]
    while stack:
        current = stack.pop()
        for child in adjacency[current]:
            if child not in reachable:
                reachable.add(child)
                stack.append(child)
    if endpoint not in reachable or reachable != set(node_map):
        return "REFUSE", ["ANCESTRY_DISCONNECTED"]
    if any(item.startswith("FAIL_") for item in diagnostics):
        return "FAIL", sorted(set(diagnostics))
    if diagnostics:
        return "UNSCOREABLE", sorted(set(diagnostics))
    return "PASS", []


PREREQUISITES: dict[str, tuple[str, ...]] = {
    "FORM.ALLOW0": (),
    "GAMMA.PROCESS": (),
    "SEED.PHI": ("GAMMA.PROCESS",),
    "SCALE.CALIBRATED_TRANSPORT": ("SEED.PHI",),
    "ANCESTRY.FULL_PATH": ("SCALE.CALIBRATED_TRANSPORT",),
    "VALIDATION.HELD_OUT": ("ANCESTRY.FULL_PATH",),
    "CROSS_SURFACE.UNCHANGED": ("VALIDATION.HELD_OUT",),
    "EVIDENCE.JOINED": ("ANCESTRY.FULL_PATH",),
    **{gate: ("EVIDENCE.JOINED",) for gate in GRAVITY_GATES},
    "COVERAGE.CONSTRUCTIVE": ("CROSS_SURFACE.UNCHANGED",),
    "WORLD.REPRODUCTION": ("VALIDATION.HELD_OUT",),
}


def validate_gate_registry(instance: dict[str, Any], registry: dict[str, dict[str, Any]], payloads: dict[str, Any], semantic_overrides: dict[str, str], applicability: dict[str, str]) -> tuple[dict[str, str], list[str]]:
    diagnostics: list[str] = []
    inputs = instance.get("gate_inputs")
    if not isinstance(inputs, list):
        return {gate: "UNSCOREABLE" for gate in GATES}, ["GATE_INPUTS_INVALID"]
    gate_map: dict[str, dict[str, Any]] = {}
    for index, item in enumerate(inputs):
        if not isinstance(item, dict) or set(item) != {"gate_id", "rule_artifact_id", "evidence_artifact_id", "power_artifact_id", "prerequisite_gate_ids"}:
            diagnostics.append(f"GATE_INPUT_CLOSED_SCHEMA:{index}")
            continue
        gate = item.get("gate_id")
        if gate not in GATES:
            diagnostics.append(f"UNKNOWN_GATE_ID:{gate}")
            continue
        if gate in gate_map:
            diagnostics.append(f"DUPLICATE_GATE_ID:{gate}")
            continue
        gate_map[gate] = item
    if set(gate_map) != set(GATES):
        diagnostics.append("REQUIRED_GATE_MISSING")
    outcomes: dict[str, str] = {}
    for gate in GATES:
        item = gate_map.get(gate)
        if item is None:
            outcomes[gate] = "UNSCOREABLE"
            continue
        if tuple(item["prerequisite_gate_ids"]) != PREREQUISITES[gate]:
            diagnostics.append(f"GATE_PREREQUISITE_REGISTRY_INVALID:{gate}")
        rule = artifact_payload(item["rule_artifact_id"], "GATE_RULE", registry, payloads, diagnostics)
        evidence = artifact_payload(item["evidence_artifact_id"], "GATE_EVIDENCE", registry, payloads, diagnostics)
        power = artifact_payload(item["power_artifact_id"], "POWER_CERTIFICATE", registry, payloads, diagnostics)
        outcome = "UNSCOREABLE"
        if all(isinstance(x, dict) for x in (rule, evidence, power)):
            if rule["gate_id"] != gate or evidence["gate_id"] != gate or power["gate_id"] != gate:
                diagnostics.append(f"GATE_ARTIFACT_JOIN_MISMATCH:{gate}")
            elif not isinstance(rule["predicate_ids"], list) or not rule["predicate_ids"] or len(set(rule["predicate_ids"])) != len(rule["predicate_ids"]):
                diagnostics.append(f"GATE_RULE_EMPTY:{gate}")
            elif not isinstance(evidence["source_artifact_ids"], list) or not evidence["source_artifact_ids"]:
                diagnostics.append(f"GATE_EVIDENCE_EMPTY:{gate}")
            elif not adequate_power({key: power[key] for key in power if key != "gate_id"}):
                outcome = "UNSCOREABLE"
                diagnostics.append(f"UNSCOREABLE_UNDERPOWERED:{gate}")
            elif evidence["reproducible"] is not True:
                diagnostics.append(f"GATE_EVIDENCE_NOT_REPRODUCIBLE:{gate}")
            elif evidence["taxonomy_match"] == "OUTSIDE":
                outcome = "UNCLASSIFIED"
            elif evidence["taxonomy_match"] != "KNOWN":
                diagnostics.append(f"GATE_TAXONOMY_INVALID:{gate}")
            else:
                observations = evidence["observations"]
                if not isinstance(observations, list) or not observations:
                    diagnostics.append(f"GATE_EVIDENCE_EMPTY:{gate}")
                elif any(not isinstance(x, dict) or set(x) != {"predicate_id", "value", "reproducible"} for x in observations):
                    diagnostics.append(f"GATE_OBSERVATION_CLOSED_SCHEMA:{gate}")
                elif {x["predicate_id"] for x in observations} != set(rule["predicate_ids"]) or any(x["reproducible"] is not True for x in observations):
                    diagnostics.append(f"GATE_EVIDENCE_RULE_JOIN_INVALID:{gate}")
                else:
                    outcome = "PASS" if all(x["value"] is True for x in observations) else "FAIL"
        if any(outcomes.get(prereq) != "PASS" for prereq in PREREQUISITES[gate]):
            outcome = "UNSCOREABLE"
        if gate in applicability and applicability[gate] != "APPLICABLE":
            outcome = "UNSCOREABLE"
        override = semantic_overrides.get(gate)
        if override == "REFUSE":
            diagnostics.append(f"SEMANTIC_GATE_REFUSE:{gate}")
            outcome = "UNSCOREABLE"
        elif override in {"PASS", "FAIL", "UNCLASSIFIED", "UNSCOREABLE"}:
            if override != "PASS":
                outcome = override
            elif outcome == "PASS":
                outcome = "PASS"
        outcomes[gate] = outcome
    return outcomes, sorted(set(diagnostics))


def validate_applicability(instance: dict[str, Any], registry: dict[str, dict[str, Any]], payloads: dict[str, Any]) -> tuple[dict[str, str], list[str]]:
    diagnostics: list[str] = []
    rows = instance.get("gravity_applicability")
    mapping: dict[str, str] = {}
    if not isinstance(rows, list):
        return mapping, ["APPLICABILITY_REGISTRY_INVALID"]
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"characteristic_id", "applicability", "justification_artifact_id"}:
            diagnostics.append("APPLICABILITY_ROW_CLOSED_SCHEMA")
            continue
        gate = row["characteristic_id"]
        if gate not in GRAVITY_GATES or gate in mapping:
            diagnostics.append(f"APPLICABILITY_ID_INVALID:{gate}")
            continue
        mapping[gate] = row["applicability"]
        if gate in GRAVITY_GATES[:13] and row["applicability"] != "APPLICABLE":
            diagnostics.append(f"MANDATORY_GRAVITY_GATE_NOT_APPLICABLE:{gate}")
        if gate == "GC14" and instance.get("horizon_complete_claim") and row["applicability"] != "APPLICABLE":
            diagnostics.append("MANDATORY_HORIZON_GATE_NOT_APPLICABLE:GC14")
        if row["applicability"] == "APPLICABLE" and row["justification_artifact_id"] is not None:
            diagnostics.append(f"APPLICABILITY_JUSTIFICATION_UNEXPECTED:{gate}")
        if row["applicability"] == "NOT_APPLICABLE":
            justification = artifact_payload(row["justification_artifact_id"], "APPLICABILITY_JUSTIFICATION", registry, payloads, diagnostics)
            if not isinstance(justification, dict) or justification.get("characteristic_id") != gate or not justification.get("physical_reason") or not justification.get("evidence_artifact_ids"):
                diagnostics.append(f"APPLICABILITY_JUSTIFICATION_INVALID:{gate}")
    if set(mapping) != set(GRAVITY_GATES):
        diagnostics.append("APPLICABILITY_REGISTRY_INCOMPLETE")
    return mapping, sorted(set(diagnostics))


def validate_domains(instance: dict[str, Any], bindings: dict[str, Any], registry: dict[str, dict[str, Any]], payloads: dict[str, Any]) -> tuple[dict[str, list[str]], list[str]]:
    diagnostics: list[str] = []
    mapping: dict[str, list[str]] = {}
    rows = instance.get("milestone_domains")
    if not isinstance(rows, list):
        return mapping, ["DOMAIN_REGISTRY_INVALID"]
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"claim_id", "domain_artifact_id"}:
            diagnostics.append("DOMAIN_ROW_CLOSED_SCHEMA")
            continue
        claim = row["claim_id"]
        if claim not in CLAIMS or claim in mapping:
            diagnostics.append(f"DOMAIN_CLAIM_ID_INVALID:{claim}")
            continue
        payload = artifact_payload(row["domain_artifact_id"], "DOMAIN_WITNESS", registry, payloads, diagnostics)
        if not isinstance(payload, dict) or payload.get("claim_id") != claim:
            diagnostics.append(f"DOMAIN_BINDING_INVALID:{claim}")
            continue
        members = payload.get("member_ids")
        if not isinstance(members, list) or not members:
            diagnostics.append(f"DOMAIN_EMPTY:{claim}")
            mapping[claim] = []
        elif len(set(members)) != len(members):
            diagnostics.append(f"DOMAIN_MEMBER_DUPLICATE:{claim}")
            mapping[claim] = members
        elif payload.get("quantifier") != "FOR_ALL_DECLARED_MEMBERS" or payload.get("coverage_artifact_id") != bindings.get("coverage_artifact_id"):
            diagnostics.append(f"DOMAIN_WITNESS_INVALID:{claim}")
            mapping[claim] = members
        else:
            mapping[claim] = members
    if set(mapping) != set(CLAIMS):
        diagnostics.append("DOMAIN_REGISTRY_INCOMPLETE")
    return mapping, sorted(set(diagnostics))


def verify_rsa_pkcs1_v15_sha256(n_hex: str, exponent: Any, message: bytes, signature_hex: str, minimum_bits: int) -> bool:
    try:
        n = int(n_hex, 16)
        e = int(exponent)
        signature = int(signature_hex, 16)
    except (TypeError, ValueError):
        return False
    if n.bit_length() < minimum_bits or e < 3 or e % 2 == 0 or signature >= n:
        return False
    length = (n.bit_length() + 7) // 8
    encoded = pow(signature, e, n).to_bytes(length, "big")
    digest_info = bytes.fromhex("3031300d060960864801650304020105000420") + hashlib.sha256(message).digest()
    padding_length = length - len(digest_info) - 3
    expected = b"\x00\x01" + b"\xff" * padding_length + b"\x00" + digest_info
    return padding_length >= 8 and encoded == expected


def validate_product(instance: dict[str, Any], bindings: dict[str, Any], registry: dict[str, dict[str, Any]], payloads: dict[str, Any]) -> tuple[str, list[str]]:
    diagnostics: list[str] = []
    kinds = {
        "public_dataset_manifest_artifact_id": "PUBLIC_DATASET_MANIFEST",
        "public_release_artifact_id": "PUBLIC_RELEASE",
        "executor_key_artifact_id": "EXECUTOR_KEY",
        "execution_output_artifact_id": "EXECUTION_OUTPUT",
        "execution_report_artifact_id": "EXECUTION_REPORT",
        "execution_signature_artifact_id": "EXECUTION_SIGNATURE",
    }
    objects: dict[str, Any] = {}
    for binding, kind in kinds.items():
        objects[binding] = artifact_payload(bindings.get(binding), kind, registry, payloads, diagnostics)
    if not all(isinstance(value, dict) for value in objects.values()):
        return "UNSCOREABLE", sorted(set(diagnostics + ["PRODUCT_BINDING_UNRESOLVED"]))
    manifest = objects["public_dataset_manifest_artifact_id"]
    release = objects["public_release_artifact_id"]
    key = objects["executor_key_artifact_id"]
    output = objects["execution_output_artifact_id"]
    report = objects["execution_report_artifact_id"]
    signature = objects["execution_signature_artifact_id"]
    if key["executor_id"] in {instance.get("producer_id"), instance.get("principal_id")}:
        diagnostics.append("EXECUTOR_NOT_INDEPENDENT")
    if release["contract_id"] != CONTRACT_ID or release["package_id"] != instance.get("package_id") or set(release["claim_ids"]) != set(CLAIMS):
        diagnostics.append("PUBLIC_RELEASE_SCOPE_MISMATCH")
    if release["dataset_manifest_artifact_id"] != bindings["public_dataset_manifest_artifact_id"]:
        diagnostics.append("PUBLIC_RELEASE_DATASET_JOIN_MISMATCH")
    if release["dataset_manifest_sha256"] != registry[bindings["public_dataset_manifest_artifact_id"]]["sha256"]:
        diagnostics.append("PUBLIC_RELEASE_DATASET_HASH_MISMATCH")
    if release["expected_output_artifact_id"] != bindings["execution_output_artifact_id"] or release["expected_output_sha256"] != registry[bindings["execution_output_artifact_id"]]["sha256"]:
        diagnostics.append("RELEASE_EXECUTION_MISMATCH")
    validator_hash = sha256_bytes(Path(__file__).read_bytes())
    if release["validator_sha256"] != validator_hash:
        diagnostics.append("VALIDATOR_IDENTITY_MISMATCH")
    if output["release_id"] != release["release_id"] or report["release_artifact_id"] != bindings["public_release_artifact_id"] or report["output_artifact_id"] != bindings["execution_output_artifact_id"]:
        diagnostics.append("PRODUCT_REPORT_JOIN_MISMATCH")
    if report["executor_id"] != key["executor_id"] or report["no_private_logic"] is not True or report["mismatch_count"] != 0:
        diagnostics.append("PRODUCT_EXECUTION_REPORT_FAILED")
    expected_datasets = {bindings["development_root_artifact_id"], bindings["validation_root_artifact_id"]}
    expected_dataset_hashes = {aid: registry[aid]["sha256"] for aid in expected_datasets if aid in registry}
    if set(manifest["dataset_artifact_ids"]) != expected_datasets or manifest["dataset_sha256s"] != expected_dataset_hashes or not manifest["public_locators"] or not manifest["license_ids"]:
        diagnostics.append("PUBLIC_DATASET_MANIFEST_INCOMPLETE")
    try:
        message = base64.b64decode(signature["message_b64"], validate=True)
        decoded = strict_json_loads(message)
    except Exception:
        message, decoded = b"", None
        diagnostics.append("EXECUTION_SIGNATURE_MESSAGE_INVALID")
    expected_message = {
        "contract_id": CONTRACT_ID,
        "package_id": instance.get("package_id"),
        "release_sha256": registry[bindings["public_release_artifact_id"]]["sha256"],
        "output_sha256": registry[bindings["execution_output_artifact_id"]]["sha256"],
        "report_sha256": registry[bindings["execution_report_artifact_id"]]["sha256"],
    }
    if decoded != expected_message:
        diagnostics.append("EXECUTION_SIGNATURE_SCOPE_MISMATCH")
    if signature["key_artifact_id"] != bindings["executor_key_artifact_id"] or signature["release_artifact_id"] != bindings["public_release_artifact_id"] or signature["output_artifact_id"] != bindings["execution_output_artifact_id"] or signature["report_artifact_id"] != bindings["execution_report_artifact_id"]:
        diagnostics.append("EXECUTION_SIGNATURE_BINDING_MISMATCH")
    minimum_bits = 2048 if instance.get("package_mode") == "SCIENTIFIC" else 512
    if not verify_rsa_pkcs1_v15_sha256(key["n_hex"], key["e"], message, signature["signature_hex"], minimum_bits):
        diagnostics.append("EXECUTION_SIGNATURE_INVALID")
    return ("PASS" if not diagnostics else "UNSCOREABLE"), sorted(set(diagnostics))


def combine_outcomes(values: Iterable[str]) -> str:
    values = list(values)
    if values and all(value == "PASS" for value in values):
        return "PASS"
    if "FAIL" in values:
        return "FAIL"
    if "UNCLASSIFIED" in values:
        return "UNCLASSIFIED"
    return "UNSCOREABLE"


def derive_milestones(gates: dict[str, str], gf0_state: str, pair_state: str, carrier_count: int, joined_state: str, coverage_state: str, product_state: str, horizon_complete: bool) -> dict[str, str]:
    milestones: dict[str, str] = {"GF0": gf0_state}
    milestones["GF1"] = combine_outcomes([milestones["GF0"], gates["FORM.ALLOW0"], gates["GAMMA.PROCESS"], gates["SEED.PHI"], gates["SCALE.CALIBRATED_TRANSPORT"], gates["ANCESTRY.FULL_PATH"]])
    milestones["GF2"] = combine_outcomes([milestones["GF1"], pair_state, gates["VALIDATION.HELD_OUT"]])
    milestones["GF3"] = combine_outcomes([milestones["GF2"], "PASS" if carrier_count >= 2 else "UNSCOREABLE", gates["CROSS_SURFACE.UNCHANGED"]])
    milestones["GE1"] = combine_outcomes([milestones["GF1"], joined_state, gates["EVIDENCE.JOINED"]])
    gravity = list(GRAVITY_GATES if horizon_complete else GRAVITY_GATES[:13])
    milestones["GE2"] = combine_outcomes([milestones["GE1"]] + [gates[gate] for gate in gravity])
    milestones["UGE"] = combine_outcomes([milestones["GF3"], milestones["GE2"], coverage_state, gates["COVERAGE.CONSTRUCTIVE"], gates["WORLD.REPRODUCTION"], product_state])
    return milestones


def proof_mapping(outcome: str) -> str:
    return {"PASS": "PASSES_DECLARED_DOMAIN", "FAIL": "REFUTED", "UNCLASSIFIED": "NO_PROOF_OUTPUT", "UNSCOREABLE": "NO_PROOF_OUTPUT"}[outcome]


def evaluate_instance(instance: Any, workspace: Path | None = None) -> dict[str, Any]:
    workspace = (workspace or Path.cwd()).resolve()
    structural: list[str] = []
    if not exact_keys(instance, TOP_KEYS, "instance", structural):
        return {"accepted": False, "custody": {"disposition": "REFUSE", "errors": sorted(set(structural)), "missing_artifact_ids": []}, "scientific_gates": {}, "candidate_milestones": {}, "authoritative_proof_outputs": {}, "product_reproduction": "UNSCOREABLE", "scientific_weight_of_product": "NONE"}
    finite_tree(instance, "instance", structural)
    if instance.get("contract_id") != CONTRACT_ID:
        structural.append("CONTRACT_ID_INVALID")
    if instance.get("task_id") != TASK_ID:
        structural.append("TASK_ID_INVALID")
    if instance.get("package_mode") not in {"SCIENTIFIC", "SYNTHETIC_TEST"}:
        structural.append("PACKAGE_MODE_INVALID")
    if not isinstance(instance.get("package_id"), str) or not re.fullmatch(r"[A-Z0-9][A-Z0-9._-]{2,127}", instance["package_id"]):
        structural.append("PACKAGE_ID_INVALID")
    if not isinstance(instance.get("producer_id"), str) or not instance["producer_id"] or not isinstance(instance.get("principal_id"), str) or not instance["principal_id"]:
        structural.append("PARTY_ID_INVALID")
    if set(instance.get("claim_ids", [])) != set(CLAIMS) or len(instance.get("claim_ids", [])) != len(CLAIMS):
        structural.append("CLAIM_REGISTRY_INVALID")
    if set(instance.get("gate_ids", [])) != set(GATES) or len(instance.get("gate_ids", [])) != len(GATES):
        structural.append("GATE_REGISTRY_INVALID")
    bindings = instance.get("bindings")
    exact_keys(bindings, BINDING_KEYS, "bindings", structural)
    registry, payloads, artifact_errors, missing = validate_artifacts(instance, workspace)
    structural.extend(artifact_errors)
    if isinstance(bindings, dict):
        validate_framework_bindings(bindings, registry, structural)

    applicability, applicability_diagnostics = validate_applicability(instance, registry, payloads)
    domains, domain_diagnostics = validate_domains(instance, bindings if isinstance(bindings, dict) else {}, registry, payloads)
    structural.extend(applicability_diagnostics)
    structural.extend(domain_diagnostics)

    gamma_state, seed_state, gamma_seed_diagnostics, processes = validate_gamma_seed(bindings, registry, payloads) if isinstance(bindings, dict) else ("UNSCOREABLE", "UNSCOREABLE", ["BINDINGS_INVALID"], [])
    structural.extend(gamma_seed_diagnostics)
    transport_payload = payloads.get(bindings.get("transport_artifact_id")) if isinstance(bindings, dict) else None
    if not isinstance(bindings, dict) or registry.get(bindings.get("transport_artifact_id"), {}).get("kind") != "CALIBRATED_TRANSPORT":
        transport_state, transport_diagnostics = "UNSCOREABLE", ["TRANSPORT_BINDING_INVALID"]
    else:
        transport_state, transport_diagnostics = validate_transport(transport_payload)
    if transport_state == "REFUSE":
        structural.extend(transport_diagnostics)
    ancestry_payload = payloads.get(bindings.get("ancestry_artifact_id")) if isinstance(bindings, dict) else None
    ancestry_state, ancestry_diagnostics = validate_ancestry(ancestry_payload, bindings.get("transport_artifact_id")) if isinstance(bindings, dict) else ("UNSCOREABLE", ["ANCESTRY_BINDING_INVALID"])
    if ancestry_state == "REFUSE":
        structural.extend(ancestry_diagnostics)
    pair_state, pair_diagnostics = validate_data_topology(bindings, registry, payloads) if isinstance(bindings, dict) else ("UNSCOREABLE", ["DATA_TOPOLOGY_INVALID"])

    joined_payload = payloads.get(bindings.get("joined_evidence_artifact_id")) if isinstance(bindings, dict) else None
    joined_state = "PASS"
    joined_diagnostics: list[str] = []
    if not isinstance(joined_payload, dict) or registry.get(bindings.get("joined_evidence_artifact_id"), {}).get("kind") != "JOINED_EVIDENCE" or any(not joined_payload.get(key) for key in ("formation_artifact_ids", "growth_artifact_ids", "response_artifact_ids", "join_keys")):
        joined_state = "UNSCOREABLE"
        joined_diagnostics.append("JOINED_EVIDENCE_MISSING_OR_EMPTY")
    coverage_payload = payloads.get(bindings.get("coverage_artifact_id")) if isinstance(bindings, dict) else None
    coverage_state = "PASS"
    coverage_diagnostics: list[str] = []
    if not isinstance(coverage_payload, dict) or registry.get(bindings.get("coverage_artifact_id"), {}).get("kind") != "COVERAGE_THEOREM" or not coverage_payload.get("surface_classes") or not coverage_payload.get("constructive_map") or not coverage_payload.get("preservation_obligations") or not coverage_payload.get("checker_artifact_ids"):
        coverage_state = "UNSCOREABLE"
        coverage_diagnostics.append("CONSTRUCTIVE_COVERAGE_INCOMPLETE")
    product_state, product_diagnostics = validate_product(instance, bindings, registry, payloads) if isinstance(bindings, dict) else ("UNSCOREABLE", ["PRODUCT_BINDINGS_INVALID"])

    platform_id = bindings.get("platform_instantiation_artifact_id") if isinstance(bindings, dict) else None
    platform = payloads.get(platform_id)
    actual_platform = bool(platform_id and registry.get(platform_id, {}).get("kind") == "PLATFORM_INSTANTIATION" and isinstance(platform, dict) and platform.get("synthetic_only") is False and instance.get("package_mode") == "SCIENTIFIC")
    platform_candidate = bool(platform_id and registry.get(platform_id, {}).get("kind") == "PLATFORM_INSTANTIATION" and isinstance(platform, dict))
    carrier_count = len({p.get("carrier_class") for p in processes if isinstance(p, dict) and p.get("carrier_class")})

    semantic_overrides = {
        "GAMMA.PROCESS": gamma_state,
        "SEED.PHI": seed_state,
        "SCALE.CALIBRATED_TRANSPORT": transport_state,
        "ANCESTRY.FULL_PATH": ancestry_state,
        "VALIDATION.HELD_OUT": pair_state,
        "CROSS_SURFACE.UNCHANGED": "PASS" if carrier_count >= 2 else "UNSCOREABLE",
        "EVIDENCE.JOINED": joined_state,
        "COVERAGE.CONSTRUCTIVE": coverage_state,
        "WORLD.REPRODUCTION": pair_state,
    }
    gates, gate_diagnostics = validate_gate_registry(instance, registry, payloads, semantic_overrides, applicability)
    structural.extend(item for item in gate_diagnostics if item.startswith(("UNKNOWN_", "DUPLICATE_", "REQUIRED_GATE_", "GATE_INPUT_CLOSED_", "GATE_PREREQUISITE_REGISTRY_", "GATE_ARTIFACT_JOIN_")))

    structural = sorted(set(structural))
    custody = "REFUSE" if structural else "UNSCOREABLE" if missing else "QUALIFIED"
    if custody != "QUALIFIED":
        gates = {gate: "UNSCOREABLE" for gate in GATES}
    gf0_candidate = "PASS" if custody == "QUALIFIED" and platform_candidate and gamma_state == seed_state == transport_state == ancestry_state == "PASS" else "UNSCOREABLE"
    # The authoritative GF0 additionally requires a real, non-synthetic platform.
    gf0_authoritative = "PASS" if gf0_candidate == "PASS" and actual_platform else "UNSCOREABLE"
    candidate_milestones = derive_milestones(gates, gf0_candidate, pair_state, carrier_count, joined_state, coverage_state, product_state, bool(instance.get("horizon_complete_claim")))
    authoritative_milestones = derive_milestones(gates, gf0_authoritative, pair_state, carrier_count, joined_state, coverage_state, product_state, bool(instance.get("horizon_complete_claim")))
    if instance.get("package_mode") != "SCIENTIFIC":
        proof_outputs = {claim: "NO_PROOF_OUTPUT" for claim in CLAIMS}
    else:
        proof_outputs = {claim: proof_mapping(authoritative_milestones[claim]) for claim in CLAIMS}
    diagnostics = sorted(set(
        transport_diagnostics + ancestry_diagnostics + pair_diagnostics + joined_diagnostics
        + coverage_diagnostics + product_diagnostics + gate_diagnostics
    ))
    return {
        "accepted": custody != "REFUSE",
        "custody": {"disposition": custody, "errors": structural, "missing_artifact_ids": sorted(set(missing))},
        "semantic_diagnostics": diagnostics,
        "scientific_gates": gates,
        "candidate_milestones": candidate_milestones,
        "authoritative_milestones": authoritative_milestones,
        "authoritative_proof_outputs": proof_outputs,
        "product_reproduction": product_state,
        "scientific_weight_of_product": "NONE",
        "principal_seed_decision": "SELECT_JOINT_SEED",
        "actual_platform_present": actual_platform,
    }


def load_candidate(path: Path) -> Any:
    return strict_json_loads(path.read_bytes())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    args = parser.parse_args()
    try:
        candidate = load_candidate(args.candidate)
        result = evaluate_instance(candidate, args.workspace)
    except Exception as exc:
        result = {
            "accepted": False,
            "custody": {"disposition": "REFUSE", "errors": [f"PARSE_OR_EVALUATION_ERROR:{type(exc).__name__}:{exc}"], "missing_artifact_ids": []},
            "scientific_gates": {}, "candidate_milestones": {}, "authoritative_milestones": {},
            "authoritative_proof_outputs": {}, "product_reproduction": "UNSCOREABLE",
            "scientific_weight_of_product": "NONE",
        }
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if result.get("accepted") else 2


if __name__ == "__main__":
    raise SystemExit(main())
