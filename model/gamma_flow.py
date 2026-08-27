"""Origin-neutral public URM gamma-flow evaluator.

This layer binds the sealed GF V002 Repair3 evaluator and the principal's
non-geometric clarification.  It adds typed scale, response, ancestry-stage,
optional metric, gravity-characteristic, and parallel metrology joins without
adding a microscopic metric or stress-energy premise.  It never edits or projects
the supplied GF candidate.  The Repair3 evaluator sees that original candidate
first; this overlay can only remove authority.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
V002 = ROOT / "LANE_GRA_O_GF_CONTRACT" / "V002"
for _path in (V002, V002 / "REPAIR1", V002 / "REPAIR2", V002 / "REPAIR3"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

# These are fixed local imports of the sealed contract, not plugin or candidate code.
import validator_v002 as _gf_base  # noqa: E402
import validator_v002_repair3 as _gf_repair3  # noqa: E402


SCHEMA = "WAC_GAMMA_FLOW_V001"
CERTIFICATE_SCHEMA = "WAC_GAMMA_FLOW_CERTIFICATE_V001"
REPAIR3_MANIFEST_PATH = (
    "LANE_GRA_O_GF_CONTRACT/V002/REPAIR3/MANIFEST_V002_REPAIR3.sha256"
)
REPAIR3_MANIFEST_SHA256 = (
    "a0d5eacb33641b7ae32f1dfce149a3cda9a10d4839d76304236b6ec7827b9911"
)
UPSTREAM_MANIFESTS = (
    (
        "LANE_GRA_O_GF_CONTRACT/V002/MANIFEST_V002.sha256",
        "274b19ce8c640f41b37421c961d6e06668501e6020064894fb6ebd6327d565c2",
    ),
    (
        "LANE_GRA_O_GF_CONTRACT/V002/REPAIR1/MANIFEST_V002_REPAIR1.sha256",
        "89a77c7394c69b3c5233a6fa91646ebf406feb3445be39dd3a1fceabcf8f3c2a",
    ),
    (
        "LANE_GRA_O_GF_CONTRACT/V002/REPAIR2/MANIFEST_V002_REPAIR2.sha256",
        "02f4fcebd15680dfc0bd46b14a7a9be90368e17dc0f7fb98654598a14d88730b",
    ),
    (REPAIR3_MANIFEST_PATH, REPAIR3_MANIFEST_SHA256),
)
PRINCIPAL_DECISION_PATH = (
    "LANE_GRA_S_JOINT_SEED_DECISION/PRINCIPAL_DECISION.md"
)
PRINCIPAL_DECISION_SHA256 = (
    "962de48aa3ced22887200933d1dd397951a2dd25cc8bf9ea335893fc888ef82d"
)
PRINCIPAL_CLARIFICATION_PATH = (
    "LANE_GRA_S_JOINT_SEED_DECISION/PRINCIPAL_CLARIFICATION_001.md"
)
PRINCIPAL_CLARIFICATION_SHA256 = (
    "2268762250f69c1ee8297ecd22cc3e67d490ccfb44062bf94ec18db9f768b8cb"
)
CHARACTERISTIC_REGISTRY_PATH = "model/gamma_flow_characteristics_v001.json"
CHARACTERISTIC_REGISTRY_SHA256 = (
    "36577b79ed22a984e2f85be46571573c59217133be3305d8b91a04fb1b8ad498"
)

DISCOVERY_STATES = {"PASS", "FAIL", "UNCLASSIFIED", "UNSCOREABLE"}
CLAIMS = ("GF0", "GF1", "GF2", "GF3", "GE1", "GE2", "UGE")
SCALE_KINDS = {
    "SPATIAL",
    "TEMPORAL",
    "CORRELATION",
    "GRAPH_HIERARCHY",
    "POPULATION",
    "MIXED",
    "OTHER_CALIBRATED",
}
EMERGENCE_STAGES = {
    "SEED",
    "COLLECTIVE",
    "RELATIONAL",
    "METRIC_CANDIDATE",
    "IR_GRAVITY",
}
STAGE_ORDER = {
    "SEED": 0,
    "COLLECTIVE": 1,
    "RELATIONAL": 2,
    "METRIC_CANDIDATE": 3,
    "IR_GRAVITY": 4,
}
REPRESENTATION_CLASSES = {
    "NON_GEOMETRIC",
    "GEOMETRIC",
    "MIXED",
    "UNCLASSIFIED",
}
METRIC_PROBE_FAMILIES = {"CLOCK", "MATTER", "LIGHT", "INDEPENDENT_PROBE"}
METROLOGY_ROLES = {
    "KNOWN_PHYSICS_CONTROL",
    "IR_SOURCE_NORMALIZATION_AND_FALSIFICATION",
}
MATERIAL_METROLOGY_SCOPE = (
    "FULL_MATERIAL_SOURCE_METROLOGY_CONTROL_AND_IR_NORMALIZATION"
)
MATERIAL_METROLOGY_UNIT_SYSTEM = "SI"
MATERIAL_COVERAGE_SI_UNITS = {
    "MASS_ENERGY": {"J", "kg"},
    "MOMENTUM_CURRENT": {"N", "kg*m/s"},
    "STRESS_PRESSURE": {"Pa"},
    "MATERIAL_COMPOSITION": {"mol/m^3", "kg/m^3"},
    "GEOMETRY_SUPPORT": {"m"},
    "WRITER_BATH_PROBE_EXCHANGE": {"J", "W"},
    "UNCERTAINTY_COVARIANCE": {"J^2", "Pa^2"},
}
MATERIAL_COVERAGE_CATEGORIES = tuple(MATERIAL_COVERAGE_SI_UNITS)
MATERIAL_COVERAGE_OPERATOR_IDS = {
    category: f"material.control.{category.lower()}"
    for category in MATERIAL_COVERAGE_CATEGORIES
}
MATERIAL_COVERAGE_ROW_KEYS = {
    "operator_id",
    "evidence_status",
    "quantity_type",
    "si_unit",
    "source_artifact_ids",
}
MATERIAL_COVERAGE_STATUSES = {"MEASURED", "BOUNDED"}

ROOT_KEYS = {
    "schema",
    "flow_id",
    "package_mode",
    "theory_bindings",
    "contract_candidate",
    "extension_artifacts",
    "scale_semantics",
    "ancestry_stages",
    "response_kernels",
    "metric_reconstructions",
    "metrology_lane",
}
BINDING_NAMES = {
    "principal_decision",
    "principal_clarification",
    "repair3_manifest",
    "gravity_characteristic_registry",
}
REF_KEYS = {"path", "sha256"}
EXTENSION_ARTIFACT_KEYS = {"artifact_id", "kind", "sha256", "payload"}
SCALE_ROW_KEYS = {"scale_id", "scale_kind", "scale_definition_artifact_id"}
STAGE_ROW_KEYS = {
    "node_id",
    "emergence_stage",
    "representation_class",
    "operator_ids",
}
RESPONSE_ROW_KEYS = {
    "surface_id",
    "scale_id",
    "input_operator_ids",
    "output_probe_ids",
    "input_policy",
    "response_kernel_artifact_id",
    "detector_map_artifact_id",
    "material_metrology_artifact_ids",
    "known_force_control_artifact_ids",
    "freeze_time",
    "response_access_time",
}
METRIC_ROW_KEYS = {
    "node_id",
    "probe_families",
    "metric_artifact_id",
    "connection_artifact_id",
    "common_cone_artifact_id",
    "powered_alternatives_artifact_id",
}
METROLOGY_LANE_KEYS = {
    "roles",
    "carrier_artifact_ids",
    "writer_artifact_ids",
    "bath_artifact_ids",
    "support_artifact_ids",
    "probe_artifact_ids",
    "byproduct_artifact_ids",
    "ordinary_force_artifact_ids",
    "apparatus_geometry_artifact_ids",
    "clock_artifact_ids",
    "read_backaction_artifact_ids",
    "uncertainty_artifact_ids",
    "material_metrology_artifact_ids",
    "known_force_control_artifact_ids",
}
METROLOGY_CATEGORIES = tuple(
    key for key in METROLOGY_LANE_KEYS if key not in {"roles"}
)

PAYLOAD_KEYS = {
    "SCALE_DEFINITION": {
        "scale_id",
        "scale_kind",
        "units",
        "support_identity",
        "aggregation_rule",
        "monotone_order_justification",
    },
    "RESPONSE_KERNEL": {
        "surface_id",
        "scale_id",
        "input_operator_ids",
        "output_probe_ids",
        "kernel_representation",
        "source_artifact_ids",
    },
    "DETECTOR_MAP": {
        "surface_id",
        "output_probe_ids",
        "detector_channels",
        "calibration_artifact_ids",
    },
    "MATERIAL_METROLOGY": {
        "operator_ids",
        "scope",
        "quantity_types",
        "unit_system",
        "source_artifact_ids",
        "coverage_map",
    },
    "KNOWN_FORCE_CONTROL": {
        "force_family",
        "null_operator_ids",
        "source_artifact_ids",
    },
    "METROLOGY_CUSTODY": {"custody_role", "source_artifact_ids"},
    "METRIC_RECONSTRUCTION": {
        "probe_families",
        "relational_observables",
        "source_artifact_ids",
    },
    "CONNECTION_RECONSTRUCTION": {
        "probe_families",
        "transport_observables",
        "source_artifact_ids",
    },
    "COMMON_CONE_TEST": {
        "probe_families",
        "propagation_observables",
        "source_artifact_ids",
    },
    "POWERED_ALTERNATIVES": {
        "alternative_ids",
        "rule_artifact_ids",
        "power_artifact_ids",
    },
}
FORBIDDEN_DERIVED_KEYS = {
    "status",
    "outcome",
    "proof",
    "proof_output",
    "scientific_verdict",
    "authorized",
    "passed",
    "failed",
}
SHA256 = re.compile(r"^[0-9a-f]{64}$")
IDENTIFIER = re.compile(r"^[A-Z0-9][A-Z0-9._:-]{1,159}$")
UTC = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


class GammaFlowRefusal(ValueError):
    """The supplied gamma-flow envelope violates the closed public contract."""


def _refuse(message: str) -> None:
    raise GammaFlowRefusal("GAMMA FLOW REFUSES: " + message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_bytes(value: Any) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        _refuse(f"value is not finite canonical JSON: {exc}")


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _refuse(f"duplicate JSON member name {key!r}")
        result[key] = value
    return result


def _parse_constant(value: str) -> None:
    _refuse(f"nonfinite JSON constant {value!r}")


def _closed(value: Any, keys: set[str], field: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        _refuse(f"{field} must be an object")
    if set(value) != keys:
        _refuse(
            f"{field} key closure failed; "
            f"missing={sorted(keys - set(value))}, extra={sorted(set(value) - keys)}"
        )
    return value


def _nonempty(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        _refuse(f"{field} must be a nonempty string")
    return value


def _identifier(value: Any, field: str) -> str:
    value = _nonempty(value, field)
    if not IDENTIFIER.fullmatch(value):
        _refuse(f"{field} is not a closed identifier")
    return value


def _digest(value: Any, field: str) -> str:
    if not isinstance(value, str) or not SHA256.fullmatch(value):
        _refuse(f"{field} must be a lowercase SHA-256")
    return value


def _unique_strings(value: Any, field: str, *, allow_empty: bool = False) -> list[str]:
    if not isinstance(value, list) or (not value and not allow_empty):
        _refuse(f"{field} must be a unique string list")
    if any(not isinstance(item, str) or not item.strip() for item in value):
        _refuse(f"{field} contains an empty or non-string member")
    if len(value) != len(set(value)):
        _refuse(f"{field} contains duplicates")
    return value


def _utc(value: Any, field: str) -> datetime:
    if not isinstance(value, str) or not UTC.fullmatch(value):
        _refuse(f"{field} must be second-resolution UTC")
    try:
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=timezone.utc
        )
    except ValueError:
        _refuse(f"{field} is not a real UTC time")


def _forbid_derived_keys(value: Any, field: str) -> None:
    if isinstance(value, dict):
        bad = set(value) & FORBIDDEN_DERIVED_KEYS
        if bad:
            _refuse(f"{field} supplies forbidden derived keys {sorted(bad)}")
        for key, child in value.items():
            _forbid_derived_keys(child, f"{field}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _forbid_derived_keys(child, f"{field}[{index}]")
    elif isinstance(value, float) and not math.isfinite(value):
        _refuse(f"{field} contains a nonfinite number")


def _combine(values: list[str]) -> str:
    if values and all(value == "PASS" for value in values):
        return "PASS"
    if "FAIL" in values:
        return "FAIL"
    if "UNCLASSIFIED" in values:
        return "UNCLASSIFIED"
    return "UNSCOREABLE"


def _verify_sha256_manifest(relative: str, expected_manifest_hash: str) -> None:
    manifest_path = ROOT / relative
    if (
        not manifest_path.is_file()
        or manifest_path.is_symlink()
        or sha256_file(manifest_path) != expected_manifest_hash
    ):
        _refuse(f"sealed manifest {relative} is absent or hash-mismatched")
    try:
        lines = manifest_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        _refuse(f"sealed manifest {relative} is unreadable: {exc}")
    if not lines:
        _refuse(f"sealed manifest {relative} is empty")
    for index, line in enumerate(lines):
        match = re.fullmatch(r"([0-9a-f]{64})  ([^\s].*)", line)
        if not match:
            _refuse(f"sealed manifest {relative} has malformed line {index + 1}")
        expected, item_relative = match.groups()
        item = ROOT / item_relative
        try:
            item.resolve(strict=True).relative_to(ROOT.resolve(strict=True))
        except (OSError, ValueError):
            _refuse(f"sealed manifest {relative} item escapes the workspace")
        if not item.is_file() or item.is_symlink() or sha256_file(item) != expected:
            _refuse(f"sealed manifest {relative} item {item_relative} is hash-mismatched")


def _verify_import_identity() -> None:
    expected = {
        _gf_base: V002 / "validator_v002.py",
        _gf_repair3: V002 / "REPAIR3" / "validator_v002_repair3.py",
        _gf_repair3.repair1: V002 / "REPAIR1" / "validator_v002_repair1.py",
        _gf_repair3.repair2: V002 / "REPAIR2" / "validator_v002_repair2.py",
    }
    for module, path in expected.items():
        try:
            imported = Path(module.__file__).resolve(strict=True)
            registered = path.resolve(strict=True)
        except (AttributeError, OSError) as exc:
            _refuse(f"sealed evaluator import identity is unavailable: {exc}")
        if imported != registered:
            _refuse(f"sealed evaluator import was shadowed: {module.__name__}")


def _theory_bindings(value: Any) -> None:
    bindings = _closed(value, BINDING_NAMES, "theory_bindings")
    expected = {
        "principal_decision": (PRINCIPAL_DECISION_PATH, PRINCIPAL_DECISION_SHA256),
        "principal_clarification": (
            PRINCIPAL_CLARIFICATION_PATH,
            PRINCIPAL_CLARIFICATION_SHA256,
        ),
        "repair3_manifest": (REPAIR3_MANIFEST_PATH, REPAIR3_MANIFEST_SHA256),
        "gravity_characteristic_registry": (
            CHARACTERISTIC_REGISTRY_PATH,
            CHARACTERISTIC_REGISTRY_SHA256,
        ),
    }
    for name, (expected_path, expected_hash) in expected.items():
        ref = _closed(bindings[name], REF_KEYS, f"theory_bindings.{name}")
        if ref["path"] != expected_path or ref["sha256"] != expected_hash:
            _refuse(f"theory_bindings.{name} does not bind the registered parent")
        path = ROOT / expected_path
        if not path.is_file() or path.is_symlink() or sha256_file(path) != expected_hash:
            _refuse(f"registered parent {expected_path} is absent or hash-mismatched")
    for relative, expected_hash in UPSTREAM_MANIFESTS:
        _verify_sha256_manifest(relative, expected_hash)
    _verify_import_identity()


def _characteristic_registry() -> dict[str, Any]:
    path = ROOT / CHARACTERISTIC_REGISTRY_PATH
    try:
        registry = json.loads(
            path.read_bytes(),
            object_pairs_hook=_unique_json_object,
            parse_constant=_parse_constant,
        )
    except GammaFlowRefusal:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        _refuse(f"gravity characteristic registry is unreadable: {exc}")
    registry = _closed(
        registry,
        {"registry_id", "source_clarification_sha256", "characteristics"},
        "gravity characteristic registry",
    )
    if registry["registry_id"] != "WAC_GRAVITY_CHARACTERISTIC_REGISTRY_V001":
        _refuse("gravity characteristic registry ID is not registered")
    if registry["source_clarification_sha256"] != PRINCIPAL_CLARIFICATION_SHA256:
        _refuse("gravity characteristic registry loses its clarification parent")
    rows = registry["characteristics"]
    if not isinstance(rows, list) or len(rows) != 14:
        _refuse("gravity characteristic registry must contain exactly fourteen rows")
    row_keys = {
        "characteristic_id",
        "characteristic_name",
        "allowed_micro_behavior",
        "intermediate_discriminators",
        "ir_endpoint_requirement",
        "terminal_scope",
    }
    ids: list[str] = []
    for index, raw in enumerate(rows):
        row = _closed(raw, row_keys, f"gravity characteristic registry row {index}")
        ids.append(_nonempty(row["characteristic_id"], f"characteristic[{index}].id"))
        for key in row_keys - {"characteristic_id", "terminal_scope"}:
            _nonempty(row[key], f"characteristic[{index}].{key}")
        expected_scope = "IR_OPTIONAL_STRONGER" if index == 13 else "IR_REQUIRED"
        if row["terminal_scope"] != expected_scope:
            _refuse(f"characteristic {row['characteristic_id']} has wrong terminal scope")
    if ids != [f"GC{index:02d}" for index in range(1, 15)]:
        _refuse("gravity characteristic IDs are not exactly GC01 through GC14")
    return registry


def _extension_registry(value: Any) -> tuple[dict[str, dict[str, Any]], list[str]]:
    if not isinstance(value, list) or not value:
        _refuse("extension_artifacts must be a nonempty list")
    registry: dict[str, dict[str, Any]] = {}
    digests: set[str] = set()
    diagnostics: list[str] = []
    for index, raw in enumerate(value):
        field = f"extension_artifacts[{index}]"
        item = _closed(raw, EXTENSION_ARTIFACT_KEYS, field)
        artifact_id = _identifier(item["artifact_id"], f"{field}.artifact_id")
        if artifact_id in registry:
            _refuse(f"duplicate extension artifact ID {artifact_id}")
        kind = _nonempty(item["kind"], f"{field}.kind")
        if kind not in PAYLOAD_KEYS:
            _refuse(f"{field}.kind is not registered")
        payload = _closed(item["payload"], PAYLOAD_KEYS[kind], f"{field}.payload")
        _forbid_derived_keys(payload, f"{field}.payload")
        for key, member in payload.items():
            member_field = f"{field}.payload.{key}"
            if key.endswith("_ids") or key in {
                "input_operator_ids",
                "output_probe_ids",
                "operator_ids",
                "null_operator_ids",
                "quantity_types",
                "detector_channels",
                "probe_families",
                "relational_observables",
                "transport_observables",
                "propagation_observables",
            }:
                _unique_strings(member, member_field)
            elif key == "coverage_map":
                if not isinstance(member, dict):
                    _refuse(f"{member_field} must be an object")
            elif key not in {"source_artifact_ids"}:
                _nonempty(member, member_field)
        if kind == "SCALE_DEFINITION" and payload["scale_kind"] not in SCALE_KINDS:
            _refuse(f"{field}.payload.scale_kind is not registered")
        expected = _digest(item["sha256"], f"{field}.sha256")
        observed = sha256_bytes(canonical_bytes(payload))
        if expected != observed:
            _refuse(f"{field} content hash mismatch")
        if expected in digests:
            _refuse(f"{field} reuses an extension payload digest")
        digests.add(expected)
        registry[artifact_id] = item
    return registry, diagnostics


def _source_ids_resolve(
    extension: dict[str, dict[str, Any]], candidate_ids: set[str], errors: list[str]
) -> None:
    all_ids = set(extension) | candidate_ids
    for artifact_id, item in extension.items():
        payload = item["payload"]
        for key in ("source_artifact_ids", "calibration_artifact_ids", "rule_artifact_ids", "power_artifact_ids"):
            if key not in payload:
                continue
            ids = _unique_strings(payload[key], f"extension artifact {artifact_id}.{key}")
            if any(value not in all_ids for value in ids):
                errors.append(f"EXTENSION_SOURCE_UNRESOLVED:{artifact_id}:{key}")


def _validate_metrology_lane(
    value: Any, extension: dict[str, dict[str, Any]], errors: list[str]
) -> dict[str, Any]:
    lane = _closed(value, METROLOGY_LANE_KEYS, "metrology_lane")
    roles = set(_unique_strings(lane["roles"], "metrology_lane.roles"))
    if roles != METROLOGY_ROLES:
        errors.append("METROLOGY_ROLES_INCOMPLETE")
    seen: set[str] = set()
    custody_kind = {
        "material_metrology_artifact_ids": "MATERIAL_METROLOGY",
        "known_force_control_artifact_ids": "KNOWN_FORCE_CONTROL",
    }
    for category in METROLOGY_CATEGORIES:
        ids = _unique_strings(lane[category], f"metrology_lane.{category}")
        for artifact_id in ids:
            if artifact_id in seen:
                errors.append(f"METROLOGY_ARTIFACT_ROLE_REUSE:{artifact_id}")
            seen.add(artifact_id)
            item = extension.get(artifact_id)
            expected_kind = custody_kind.get(category, "METROLOGY_CUSTODY")
            if not item or item["kind"] != expected_kind:
                errors.append(f"METROLOGY_ARTIFACT_KIND_MISMATCH:{category}:{artifact_id}")
                continue
            if expected_kind == "METROLOGY_CUSTODY" and item["payload"].get("custody_role") != category:
                errors.append(f"METROLOGY_CUSTODY_ROLE_MISMATCH:{category}:{artifact_id}")
    return lane


def _validate_material_metrology(
    lane: dict[str, Any],
    extension: dict[str, dict[str, Any]],
    candidate_registry: dict[str, Any],
    payloads: dict[str, Any],
    bindings: dict[str, Any],
    seed_operators: set[str],
    errors: list[str],
) -> None:
    """Require complete material/source controls without making them a GF0 premise."""
    dilation_id = bindings.get("dilation_artifact_id")
    dilation_item = candidate_registry.get(dilation_id)
    dilation = payloads.get(dilation_id)
    measurement_ids = (
        dilation.get("measurement_artifact_ids")
        if isinstance(dilation, dict)
        else None
    )
    if (
        not isinstance(dilation_id, str)
        or not isinstance(dilation_item, dict)
        or dilation_item.get("kind") != "FULL_DILATION"
        or not isinstance(measurement_ids, list)
        or not measurement_ids
        or any(
            not isinstance(item, str) or item not in candidate_registry
            for item in measurement_ids
        )
    ):
        errors.append("MATERIAL_METROLOGY_BOUND_FULL_DILATION_UNAVAILABLE")
        required_source_ids: set[str] = set()
    else:
        required_source_ids = {dilation_id, *measurement_ids}

    candidate_ids = set(candidate_registry)
    expected_quantity_types = list(MATERIAL_COVERAGE_CATEGORIES)
    for artifact_id in lane["material_metrology_artifact_ids"]:
        item = extension.get(artifact_id)
        if not item or item.get("kind") != "MATERIAL_METROLOGY":
            continue
        payload = item["payload"]
        if payload["scope"] != MATERIAL_METROLOGY_SCOPE:
            errors.append(f"MATERIAL_METROLOGY_SCOPE_UNREGISTERED:{artifact_id}")
        if payload["unit_system"] != MATERIAL_METROLOGY_UNIT_SYSTEM:
            errors.append(f"MATERIAL_METROLOGY_UNIT_SYSTEM_UNREGISTERED:{artifact_id}")

        coverage = _closed(
            payload["coverage_map"],
            set(MATERIAL_COVERAGE_CATEGORIES),
            f"extension artifact {artifact_id}.coverage_map",
        )
        coverage_operator_ids: list[str] = []
        coverage_quantity_types: list[str] = []
        coverage_source_union: set[str] = set()
        for category in MATERIAL_COVERAGE_CATEGORIES:
            field = f"extension artifact {artifact_id}.coverage_map.{category}"
            row = _closed(coverage[category], MATERIAL_COVERAGE_ROW_KEYS, field)
            operator_id = _nonempty(row["operator_id"], f"{field}.operator_id")
            if (
                operator_id != MATERIAL_COVERAGE_OPERATOR_IDS[category]
                or operator_id in seed_operators
            ):
                errors.append(
                    f"MATERIAL_METROLOGY_OPERATOR_NOT_NONSEED_CONTROL:"
                    f"{artifact_id}:{category}"
                )
            if row["evidence_status"] not in MATERIAL_COVERAGE_STATUSES:
                _refuse(f"{field}.evidence_status is not MEASURED or BOUNDED")
            if row["quantity_type"] != category:
                errors.append(
                    f"MATERIAL_METROLOGY_QUANTITY_SEMANTICS_MISMATCH:"
                    f"{artifact_id}:{category}"
                )
            if row["si_unit"] not in MATERIAL_COVERAGE_SI_UNITS[category]:
                errors.append(
                    f"MATERIAL_METROLOGY_SI_UNIT_INVALID:{artifact_id}:{category}"
                )
            source_ids = set(
                _unique_strings(row["source_artifact_ids"], f"{field}.source_artifact_ids")
            )
            if not source_ids <= candidate_ids:
                errors.append(
                    f"MATERIAL_METROLOGY_SOURCE_OUTSIDE_CANDIDATE_CUSTODY:"
                    f"{artifact_id}:{category}"
                )
            if required_source_ids and not required_source_ids <= source_ids:
                errors.append(
                    f"MATERIAL_METROLOGY_DILATION_JOIN_MISSING:"
                    f"{artifact_id}:{category}"
                )
            coverage_operator_ids.append(operator_id)
            coverage_quantity_types.append(row["quantity_type"])
            coverage_source_union.update(source_ids)

        if len(coverage_operator_ids) != len(set(coverage_operator_ids)):
            errors.append(f"MATERIAL_METROLOGY_OPERATOR_REUSE:{artifact_id}")
        if payload["operator_ids"] != coverage_operator_ids:
            errors.append(f"MATERIAL_METROLOGY_OPERATOR_COVERAGE_MISMATCH:{artifact_id}")
        if (
            payload["quantity_types"] != coverage_quantity_types
            or payload["quantity_types"] != expected_quantity_types
        ):
            errors.append(f"MATERIAL_METROLOGY_QUANTITY_COVERAGE_MISMATCH:{artifact_id}")
        top_source_ids = set(payload["source_artifact_ids"])
        if top_source_ids != coverage_source_union:
            errors.append(f"MATERIAL_METROLOGY_SOURCE_COVERAGE_MISMATCH:{artifact_id}")
        if not top_source_ids <= candidate_ids:
            errors.append(
                f"MATERIAL_METROLOGY_SOURCE_OUTSIDE_CANDIDATE_CUSTODY:{artifact_id}:TOP"
            )
        if required_source_ids and not required_source_ids <= top_source_ids:
            errors.append(f"MATERIAL_METROLOGY_DILATION_JOIN_MISSING:{artifact_id}:TOP")


def _candidate_payloads(candidate: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    try:
        registry, payloads, errors, missing = _gf_base.validate_artifacts(candidate, ROOT)
    except Exception as exc:
        return {}, {}, [f"CANDIDATE_PAYLOAD_READ_ERROR:{type(exc).__name__}:{exc}"]
    diagnostics = list(errors) + [f"CANDIDATE_ARTIFACT_MISSING:{item}" for item in missing]
    return registry, payloads, diagnostics


def _seed_operators(candidate: dict[str, Any], payloads: dict[str, Any], errors: list[str]) -> set[str]:
    bindings = candidate.get("bindings", {}) if isinstance(candidate, dict) else {}
    seed = payloads.get(bindings.get("seed_definition_artifact_id"), {})
    blocks = seed.get("blocks") if isinstance(seed, dict) else None
    operators: set[str] = set()
    families: set[str] = set()
    if not isinstance(blocks, list):
        errors.append("GAMMA_FLOW_SEED_BLOCKS_MISSING")
        return operators
    for block in blocks:
        if not isinstance(block, dict):
            errors.append("GAMMA_FLOW_SEED_BLOCK_INVALID")
            continue
        families.add(str(block.get("family")))
        operators.update(block.get("operator_ids") or [])
        if block.get("family") == "D" and "probability" in str(block.get("physical_type", "")).lower():
            errors.append("PROBABILITY_CURRENT_IS_NOT_PHYSICAL_D")
    if families != {"B", "C", "D"} or not operators:
        errors.append("GAMMA_FLOW_JOINT_SEED_INCOMPLETE")
    dilation_id = bindings.get("dilation_artifact_id")
    dilation = payloads.get(dilation_id)
    if not isinstance(dilation, dict) or not dilation.get("exchange_ledger") or not dilation.get("measurement_artifact_ids"):
        errors.append("PHYSICAL_D_EXCHANGE_LEDGER_MISSING")
    return operators


def _validate_scale_semantics(
    value: Any,
    transport: Any,
    extension: dict[str, dict[str, Any]],
    errors: list[str],
) -> set[str]:
    if not isinstance(value, list) or not value:
        _refuse("scale_semantics must be a nonempty list")
    candidate_scales = transport.get("scales") if isinstance(transport, dict) else None
    if not isinstance(candidate_scales, list):
        errors.append("TRANSPORT_SCALES_UNAVAILABLE")
        candidate_scales = []
    expected_ids = [row.get("scale_id") for row in candidate_scales if isinstance(row, dict)]
    seen: list[str] = []
    for index, raw in enumerate(value):
        field = f"scale_semantics[{index}]"
        row = _closed(raw, SCALE_ROW_KEYS, field)
        scale_id = _identifier(row["scale_id"], f"{field}.scale_id")
        if row["scale_kind"] not in SCALE_KINDS:
            _refuse(f"{field}.scale_kind is not registered")
        definition_id = _identifier(
            row["scale_definition_artifact_id"],
            f"{field}.scale_definition_artifact_id",
        )
        item = extension.get(definition_id)
        if not item or item["kind"] != "SCALE_DEFINITION":
            errors.append(f"SCALE_DEFINITION_UNRESOLVED:{scale_id}")
        elif (
            item["payload"]["scale_id"] != scale_id
            or item["payload"]["scale_kind"] != row["scale_kind"]
        ):
            errors.append(f"SCALE_DEFINITION_JOIN_MISMATCH:{scale_id}")
        seen.append(scale_id)
    if seen != expected_ids or len(seen) != len(set(seen)):
        errors.append("SCALE_SEMANTICS_NOT_EXACTLY_TRANSPORT_SCALES")
    return set(seen)


def _validate_response_kernels(
    value: Any,
    scale_ids: set[str],
    transport: Any,
    extension: dict[str, dict[str, Any]],
    lane: dict[str, Any],
    seed_operators: set[str],
    errors: list[str],
) -> tuple[set[str], set[str]]:
    if not isinstance(value, list) or not value:
        _refuse("response_kernels must be a nonempty list")
    seen_scales: set[str] = set()
    output_probes: set[str] = set()
    transport_surface = transport.get("surface_id") if isinstance(transport, dict) else None
    transport_freeze = transport.get("freeze_time") if isinstance(transport, dict) else None
    transport_access = transport.get("response_access_time") if isinstance(transport, dict) else None
    allowed_material = set(lane["material_metrology_artifact_ids"])
    allowed_controls = set(lane["known_force_control_artifact_ids"])
    for index, raw in enumerate(value):
        field = f"response_kernels[{index}]"
        row = _closed(raw, RESPONSE_ROW_KEYS, field)
        surface_id = _nonempty(row["surface_id"], f"{field}.surface_id")
        scale_id = _identifier(row["scale_id"], f"{field}.scale_id")
        if scale_id in seen_scales or scale_id not in scale_ids:
            errors.append(f"RESPONSE_SCALE_DUPLICATE_OR_UNKNOWN:{scale_id}")
        seen_scales.add(scale_id)
        if surface_id != transport_surface:
            errors.append(f"RESPONSE_SURFACE_JOIN_MISMATCH:{scale_id}")
        if row["input_policy"] != "JOINT_FROZEN_BASIS":
            errors.append(f"EXCLUSIVE_SOURCE_POLICY_FORBIDDEN:{scale_id}")
        input_ids = set(_unique_strings(row["input_operator_ids"], f"{field}.input_operator_ids"))
        probe_ids = set(_unique_strings(row["output_probe_ids"], f"{field}.output_probe_ids"))
        material_ids = set(
            _unique_strings(
                row["material_metrology_artifact_ids"],
                f"{field}.material_metrology_artifact_ids",
            )
        )
        control_ids = set(
            _unique_strings(
                row["known_force_control_artifact_ids"],
                f"{field}.known_force_control_artifact_ids",
            )
        )
        if material_ids != allowed_material or control_ids != allowed_controls:
            errors.append(f"RESPONSE_METROLOGY_OUTSIDE_FROZEN_LANE:{scale_id}")
        material_operators = {
            operator
            for artifact_id in material_ids
            for operator in extension.get(artifact_id, {}).get("payload", {}).get("operator_ids", [])
        }
        control_operators = {
            operator
            for artifact_id in control_ids
            for operator in extension.get(artifact_id, {}).get("payload", {}).get("null_operator_ids", [])
        }
        expected_basis = seed_operators | material_operators | control_operators
        if (
            input_ids != expected_basis
            or not input_ids & seed_operators
            or not material_operators
            or not control_operators
        ):
            errors.append(f"RESPONSE_JOINT_BASIS_INCOMPLETE:{scale_id}")
        kernel_id = _identifier(
            row["response_kernel_artifact_id"], f"{field}.response_kernel_artifact_id"
        )
        detector_id = _identifier(
            row["detector_map_artifact_id"], f"{field}.detector_map_artifact_id"
        )
        kernel = extension.get(kernel_id)
        detector = extension.get(detector_id)
        if not kernel or kernel["kind"] != "RESPONSE_KERNEL":
            errors.append(f"RESPONSE_KERNEL_UNRESOLVED:{scale_id}")
        else:
            payload = kernel["payload"]
            if (
                payload["surface_id"] != surface_id
                or payload["scale_id"] != scale_id
                or set(payload["input_operator_ids"]) != input_ids
                or set(payload["output_probe_ids"]) != probe_ids
            ):
                errors.append(f"RESPONSE_KERNEL_JOIN_MISMATCH:{scale_id}")
        if not detector or detector["kind"] != "DETECTOR_MAP":
            errors.append(f"DETECTOR_MAP_UNRESOLVED:{scale_id}")
        elif (
            detector["payload"]["surface_id"] != surface_id
            or set(detector["payload"]["output_probe_ids"]) != probe_ids
        ):
            errors.append(f"DETECTOR_MAP_JOIN_MISMATCH:{scale_id}")
        freeze = _utc(row["freeze_time"], f"{field}.freeze_time")
        access = _utc(row["response_access_time"], f"{field}.response_access_time")
        if not freeze < access:
            errors.append(f"RESPONSE_NOT_FROZEN_BEFORE_ACCESS:{scale_id}")
        if row["freeze_time"] != transport_freeze or row["response_access_time"] != transport_access:
            errors.append(f"RESPONSE_TRANSPORT_CHRONOLOGY_MISMATCH:{scale_id}")
        output_probes.update(probe_ids)
    if seen_scales != scale_ids:
        errors.append("RESPONSE_KERNELS_NOT_EXACTLY_TRANSPORT_SCALES")
    return output_probes, seen_scales


def _validate_metric_rows(
    value: Any,
    stage_by_node: dict[str, str],
    extension: dict[str, dict[str, Any]],
    errors: list[str],
) -> set[str]:
    if not isinstance(value, list):
        _refuse("metric_reconstructions must be a list")
    metric_nodes: set[str] = set()
    kind_by_key = {
        "metric_artifact_id": "METRIC_RECONSTRUCTION",
        "connection_artifact_id": "CONNECTION_RECONSTRUCTION",
        "common_cone_artifact_id": "COMMON_CONE_TEST",
        "powered_alternatives_artifact_id": "POWERED_ALTERNATIVES",
    }
    for index, raw in enumerate(value):
        field = f"metric_reconstructions[{index}]"
        row = _closed(raw, METRIC_ROW_KEYS, field)
        node_id = _identifier(row["node_id"], f"{field}.node_id")
        if node_id in metric_nodes:
            errors.append(f"METRIC_NODE_DUPLICATE:{node_id}")
        metric_nodes.add(node_id)
        probes = set(_unique_strings(row["probe_families"], f"{field}.probe_families"))
        if probes != METRIC_PROBE_FAMILIES:
            errors.append(f"METRIC_MULTI_PROBE_BASIS_INCOMPLETE:{node_id}")
        if stage_by_node.get(node_id) not in {"METRIC_CANDIDATE", "IR_GRAVITY"}:
            errors.append(f"METRIC_RECONSTRUCTION_WRONG_STAGE:{node_id}")
        for key, kind in kind_by_key.items():
            artifact_id = _identifier(row[key], f"{field}.{key}")
            item = extension.get(artifact_id)
            if not item or item["kind"] != kind:
                errors.append(f"METRIC_ARTIFACT_UNRESOLVED:{node_id}:{key}")
            elif "probe_families" in item["payload"] and set(item["payload"]["probe_families"]) != probes:
                errors.append(f"METRIC_ARTIFACT_PROBE_MISMATCH:{node_id}:{key}")
    required_nodes = {
        node_id
        for node_id, stage in stage_by_node.items()
        if stage in {"METRIC_CANDIDATE", "IR_GRAVITY"}
    }
    if metric_nodes != required_nodes:
        errors.append("METRIC_RECONSTRUCTIONS_NOT_EXACTLY_METRIC_STAGE_NODES")
    return metric_nodes


def _validate_ancestry_stages(
    value: Any,
    ancestry: Any,
    seed_operators: set[str],
    output_probes: set[str],
    errors: list[str],
) -> dict[str, str]:
    if not isinstance(value, list) or not value:
        _refuse("ancestry_stages must be a nonempty list")
    nodes = ancestry.get("nodes") if isinstance(ancestry, dict) else None
    arrows = ancestry.get("arrows") if isinstance(ancestry, dict) else None
    if not isinstance(nodes, list) or not isinstance(arrows, list):
        errors.append("ANCESTRY_GRAPH_UNAVAILABLE_FOR_STAGE_JOIN")
        nodes, arrows = [], []
    expected_ids = [node.get("node_id") for node in nodes if isinstance(node, dict)]
    stage_by_node: dict[str, str] = {}
    available_operators = seed_operators | output_probes
    for index, raw in enumerate(value):
        field = f"ancestry_stages[{index}]"
        row = _closed(raw, STAGE_ROW_KEYS, field)
        node_id = _identifier(row["node_id"], f"{field}.node_id")
        if node_id in stage_by_node:
            errors.append(f"ANCESTRY_STAGE_DUPLICATE:{node_id}")
        stage = row["emergence_stage"]
        if stage not in EMERGENCE_STAGES:
            _refuse(f"{field}.emergence_stage is not registered")
        if row["representation_class"] not in REPRESENTATION_CLASSES:
            _refuse(f"{field}.representation_class is not registered")
        operators = set(_unique_strings(row["operator_ids"], f"{field}.operator_ids"))
        if not operators <= available_operators:
            errors.append(f"ANCESTRY_STAGE_OPERATOR_UNRESOLVED:{node_id}")
        stage_by_node[node_id] = stage
    if list(stage_by_node) != expected_ids or len(stage_by_node) != len(expected_ids):
        errors.append("ANCESTRY_STAGES_NOT_EXACTLY_GRAPH_NODES")
    seed_node = ancestry.get("seed_node_id") if isinstance(ancestry, dict) else None
    if stage_by_node.get(seed_node) != "SEED":
        errors.append("ANCESTRY_SEED_STAGE_MISMATCH")
    seed_row = next(
        (row for row in value if isinstance(row, dict) and row.get("node_id") == seed_node),
        None,
    )
    if not seed_row or set(seed_row.get("operator_ids", [])) != seed_operators:
        errors.append("ANCESTRY_SEED_OPERATOR_JOIN_MISMATCH")
    for arrow in arrows:
        if not isinstance(arrow, dict):
            continue
        source_stage = stage_by_node.get(arrow.get("source"))
        target_stage = stage_by_node.get(arrow.get("target"))
        if source_stage is None or target_stage is None:
            errors.append("ANCESTRY_STAGE_ARROW_UNRESOLVED")
        elif STAGE_ORDER[target_stage] < STAGE_ORDER[source_stage]:
            errors.append(f"ANCESTRY_STAGE_REGRESSION:{arrow.get('arrow_id', '?')}")
    return stage_by_node


def _stage_results(
    repair: dict[str, Any], stage_by_node: dict[str, str], metric_nodes: set[str]
) -> tuple[dict[str, str], dict[str, list[str]]]:
    milestones = repair.get("authoritative_milestones", {})
    gates = repair.get("scientific_gates", {})
    actual = repair.get("actual_platform_present") is True
    stages = set(stage_by_node.values())
    early = milestones.get("GF1", "UNSCOREABLE")
    collective = (
        milestones.get("GE1", "UNSCOREABLE")
        if "COLLECTIVE" in stages
        else "UNCLASSIFIED"
    )
    if not ({"METRIC_CANDIDATE", "IR_GRAVITY"} & stages):
        metric = "UNCLASSIFIED"
    elif not actual:
        metric = "UNSCOREABLE"
    else:
        metric = gates.get("GC08", "UNSCOREABLE")
    endpoint = (
        milestones.get("GE2", "UNSCOREABLE")
        if "IR_GRAVITY" in stages
        else "UNCLASSIFIED"
    )
    results = {
        "SEED_OPERATOR_ANCESTRY": early,
        "GENERIC_COLLECTIVE_RESPONSE": collective,
        "OPTIONAL_METRIC_RECONSTRUCTION": metric,
        "CLASSICAL_GRAVITY_ENDPOINT": endpoint,
    }
    for key, value in results.items():
        if value not in DISCOVERY_STATES:
            results[key] = "UNSCOREABLE"
    reasons = {
        "SEED_OPERATOR_ANCESTRY": [f"REPAIR3_GF1:{results['SEED_OPERATOR_ANCESTRY']}"],
        "GENERIC_COLLECTIVE_RESPONSE": [
            "COLLECTIVE_STAGE_PRESENT" if "COLLECTIVE" in stages else "COLLECTIVE_STAGE_ABSENT",
            f"REPAIR3_GE1:{milestones.get('GE1', 'UNSCOREABLE')}",
        ],
        "OPTIONAL_METRIC_RECONSTRUCTION": [
            "NO_METRIC_STAGE_DECLARED_EARLY_METRIC_NOT_REQUIRED"
            if not metric_nodes
            else f"MULTI_PROBE_METRIC_NODES:{len(metric_nodes)}"
        ],
        "CLASSICAL_GRAVITY_ENDPOINT": [
            "IR_GRAVITY_STAGE_ABSENT"
            if "IR_GRAVITY" not in stages
            else f"REPAIR3_GE2:{milestones.get('GE2', 'UNSCOREABLE')}"
        ],
    }
    return results, reasons


def _monotone_proof_outputs(
    repair: dict[str, Any], stages: dict[str, str], extension_errors: list[str]
) -> dict[str, str]:
    upstream = repair.get("authoritative_proof_outputs", {})
    outputs = {claim: upstream.get(claim, "NO_PROOF_OUTPUT") for claim in CLAIMS}
    if extension_errors:
        return {claim: "NO_PROOF_OUTPUT" for claim in CLAIMS}
    if stages["SEED_OPERATOR_ANCESTRY"] != "PASS":
        for claim in ("GF1", "GF2", "GF3", "GE1", "GE2", "UGE"):
            if outputs[claim] == "PASSES_DECLARED_DOMAIN":
                outputs[claim] = "NO_PROOF_OUTPUT"
    if stages["GENERIC_COLLECTIVE_RESPONSE"] != "PASS":
        for claim in ("GE1", "GE2", "UGE"):
            if outputs[claim] == "PASSES_DECLARED_DOMAIN":
                outputs[claim] = "NO_PROOF_OUTPUT"
    if stages["CLASSICAL_GRAVITY_ENDPOINT"] != "PASS":
        for claim in ("GE2", "UGE"):
            if outputs[claim] == "PASSES_DECLARED_DOMAIN":
                outputs[claim] = "NO_PROOF_OUTPUT"
    return outputs


@dataclass(frozen=True)
class GammaFlow:
    manifest_path: Path
    manifest_sha256: str
    manifest: dict[str, Any]
    repair3_result: dict[str, Any]
    extension_errors: tuple[str, ...]
    internal_discovery_states: tuple[tuple[str, str], ...]
    stage_reasons: tuple[tuple[str, tuple[str, ...]], ...]
    proof_outputs: tuple[tuple[str, str], ...]
    candidate_identity_preserved: bool

    def certificate(self) -> dict[str, Any]:
        states = dict(self.internal_discovery_states)
        outputs = dict(self.proof_outputs)
        reasons = {key: list(value) for key, value in self.stage_reasons}
        actual = self.repair3_result.get("actual_platform_present") is True
        accepted = (
            self.repair3_result.get("accepted") is True
            and not self.extension_errors
            and self.candidate_identity_preserved
        )
        return {
            "schema": CERTIFICATE_SCHEMA,
            "flow_id": self.manifest["flow_id"],
            "manifest_sha256": self.manifest_sha256,
            "gamma_flow_source_sha256": sha256_file(Path(__file__).resolve()),
            "repair3_manifest_sha256": REPAIR3_MANIFEST_SHA256,
            "principal_decision_sha256": PRINCIPAL_DECISION_SHA256,
            "principal_clarification_sha256": PRINCIPAL_CLARIFICATION_SHA256,
            "gravity_characteristic_registry_sha256": CHARACTERISTIC_REGISTRY_SHA256,
            "accepted": accepted,
            "contract_candidate_identity_preserved": self.candidate_identity_preserved,
            "repair3_custody": self.repair3_result.get("custody", {}).get(
                "disposition", "REFUSE"
            ),
            "extension_errors": list(self.extension_errors),
            "actual_platform_present": actual,
            "internal_discovery_states": states,
            "stage_reasons": reasons,
            "authoritative_proof_outputs": outputs,
            "scientific_result": (
                "NO_PROOF_OUTPUT"
                if all(value == "NO_PROOF_OUTPUT" for value in outputs.values())
                else "SEE_AUTHORITATIVE_PROOF_OUTPUTS"
            ),
            "metric_required_at_GF0_or_GF1": False,
            "geometry_required_at_GF0_or_GF1": False,
            "stress_energy_required_at_GF0_or_GF1": False,
            "t_lab_exclusive_microscopic_source_authorized": False,
            "full_material_metrology_role": (
                "KNOWN_PHYSICS_CONTROL_AND_IR_SOURCE_NORMALIZATION_FALSIFICATION"
            ),
            "scientific_validation_authorized": bool(
                accepted and actual and outputs["GF2"] == "PASSES_DECLARED_DOMAIN"
            ),
            "record_formation_claim_authorized": bool(
                accepted and actual and outputs["GF1"] == "PASSES_DECLARED_DOMAIN"
            ),
            "scoped_gamma_flow_claim_authorized": bool(
                accepted and actual and outputs["GE1"] == "PASSES_DECLARED_DOMAIN"
            ),
            "gravity_claim_authorized": bool(
                accepted and actual and outputs["GE2"] == "PASSES_DECLARED_DOMAIN"
            ),
            "universal_claim_authorized": bool(
                accepted and actual and outputs["UGE"] == "PASSES_DECLARED_DOMAIN"
            ),
            "program_completion_authorized": bool(
                accepted and actual and outputs["UGE"] == "PASSES_DECLARED_DOMAIN"
            ),
            "product_reproduction_scientific_weight": "NONE",
            "synthetic_fixture_scientific_weight": "ZERO",
        }


def load_gamma_flow(manifest_path: str | Path) -> GammaFlow:
    """Load and evaluate one closed origin-neutral gamma-flow envelope."""
    manifest_path = Path(manifest_path)
    if manifest_path.is_symlink() or manifest_path.parent.is_symlink():
        _refuse("submitted manifest or its directory must not be a symlink")
    try:
        raw = manifest_path.read_bytes()
        manifest = json.loads(
            raw,
            object_pairs_hook=_unique_json_object,
            parse_constant=_parse_constant,
        )
    except GammaFlowRefusal:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        _refuse(f"manifest is not readable JSON: {exc}")
    manifest = _closed(manifest, ROOT_KEYS, "manifest")
    if manifest["schema"] != SCHEMA:
        _refuse(f"unsupported schema {manifest['schema']!r}")
    _identifier(manifest["flow_id"], "flow_id")
    if manifest["package_mode"] not in {"SCIENTIFIC", "SYNTHETIC_TEST"}:
        _refuse("package_mode is not registered")
    _theory_bindings(manifest["theory_bindings"])
    _characteristic_registry()
    extension, _ = _extension_registry(manifest["extension_artifacts"])

    candidate = manifest["contract_candidate"]
    candidate_before = canonical_bytes(candidate)
    repair3_result = _gf_repair3.evaluate_instance(candidate, ROOT)
    candidate_identity_preserved = canonical_bytes(candidate) == candidate_before
    errors: list[str] = []
    if not candidate_identity_preserved:
        errors.append("REPAIR3_OR_EXTENSION_MUTATED_CONTRACT_CANDIDATE")
    if not isinstance(candidate, dict):
        errors.append("CONTRACT_CANDIDATE_NOT_OBJECT")
        candidate = {}
    if candidate.get("package_mode") != manifest["package_mode"]:
        errors.append("PACKAGE_MODE_JOIN_MISMATCH")
    candidate_registry, payloads, payload_errors = _candidate_payloads(candidate)
    errors.extend(payload_errors)
    _source_ids_resolve(extension, set(candidate_registry), errors)
    lane = _validate_metrology_lane(manifest["metrology_lane"], extension, errors)

    bindings = candidate.get("bindings", {}) if isinstance(candidate, dict) else {}
    transport = payloads.get(bindings.get("transport_artifact_id"))
    ancestry = payloads.get(bindings.get("ancestry_artifact_id"))
    seed_operators = _seed_operators(candidate, payloads, errors)
    _validate_material_metrology(
        lane,
        extension,
        candidate_registry,
        payloads,
        bindings,
        seed_operators,
        errors,
    )
    scale_ids = _validate_scale_semantics(
        manifest["scale_semantics"], transport, extension, errors
    )
    output_probes, _ = _validate_response_kernels(
        manifest["response_kernels"],
        scale_ids,
        transport,
        extension,
        lane,
        seed_operators,
        errors,
    )
    stage_by_node = _validate_ancestry_stages(
        manifest["ancestry_stages"],
        ancestry,
        seed_operators,
        output_probes,
        errors,
    )
    metric_nodes = _validate_metric_rows(
        manifest["metric_reconstructions"], stage_by_node, extension, errors
    )

    ir_nodes = {node for node, stage in stage_by_node.items() if stage == "IR_GRAVITY"}
    required_gc = list(_gf_base.GRAVITY_GATES)
    if not candidate.get("horizon_complete_claim"):
        required_gc = required_gc[:13]
    if ir_nodes and (
        repair3_result.get("actual_platform_present") is not True
        or repair3_result.get("authoritative_milestones", {}).get("GE2") != "PASS"
        or any(
            repair3_result.get("scientific_gates", {}).get(gate) != "PASS"
            for gate in required_gc
        )
    ):
        errors.append("IR_GRAVITY_LABEL_WITHOUT_CHARACTERISTIC_CONJUNCTION")

    errors = sorted(set(errors))
    states, reasons = _stage_results(repair3_result, stage_by_node, metric_nodes)
    if errors:
        states = {key: "UNSCOREABLE" for key in states}
        for key in reasons:
            reasons[key].append("GAMMA_FLOW_EXTENSION_INVALID")
    outputs = _monotone_proof_outputs(repair3_result, states, errors)
    try:
        resolved = manifest_path.resolve(strict=True)
    except OSError:
        _refuse("submitted manifest cannot be resolved")
    return GammaFlow(
        manifest_path=resolved,
        manifest_sha256=sha256_bytes(raw),
        manifest=manifest,
        repair3_result=repair3_result,
        extension_errors=tuple(errors),
        internal_discovery_states=tuple(sorted(states.items())),
        stage_reasons=tuple(
            sorted((key, tuple(value)) for key, value in reasons.items())
        ),
        proof_outputs=tuple((claim, outputs[claim]) for claim in CLAIMS),
        candidate_identity_preserved=candidate_identity_preserved,
    )


def certificate_json(flow: GammaFlow) -> str:
    return json.dumps(flow.certificate(), indent=2, sort_keys=True, allow_nan=False) + "\n"
