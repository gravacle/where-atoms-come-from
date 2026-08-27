"""Generic, closed public-URM contract for physical formation datasets.

The contract validates content custody, physical identities, temporal coverage,
frozen rules, disjoint validation datasets, SI metadata, and declarative generic
predicates.  It contains no surface-specific physics or numerical acceptance
threshold.  DEVELOPMENT data may freeze a prediction but can never validate it.

This is an explicit dataset-validation extension.  It does not modify or loosen the
sealed FORMATION_PROTOCOL_V001 experiment contract.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any
import unicodedata


SCHEMA = "WAC_FORMATION_BUNDLE_V002"
CERTIFICATE_SCHEMA = "WAC_FORMATION_BUNDLE_CERTIFICATE_V002"
EXECUTION_SCHEMA = "WAC_FORMATION_EXECUTION_V002"
EXECUTION_CERTIFICATE_SCHEMA = "WAC_FORMATION_EXECUTION_CERTIFICATE_V002"
BASE_CONTRACT_VERSION = "FORMATION_PROTOCOL_V001"
VALIDATION_EXTENSION_VERSION = "FORMATION_DATASET_VALIDATION_V002"
BASE_SCHEMA_SHA256 = "044d447a5007653c1245a2e3ceea725aad60fd977502721a3c06bff1c9be48ce"

ROOT_KEYS = {
    "schema",
    "bundle_id",
    "protocol_id",
    "base_contract_version",
    "validation_extension_version",
    "claim_scope",
    "frozen_objects",
    "datasets",
    "generic_predicates",
    "product_reproduction",
    "signatures",
}
CLAIM_SCOPE_KEYS = {
    "claim_id",
    "surface_scope",
    "formation_mode",
    "provenance",
    "mission_time",
    "requested_result",
    "universal_claim_requested",
    "gravity_claim_requested",
}
FROZEN_OBJECT_NAMES = {
    "protocol",
    "law",
    "prediction",
    "outcome_map",
    "declarative_mapping",
    "analysis",
    "urm_release",
    "urm_source",
    "urm_environment",
}
FROZEN_REF_KEYS = {
    "object_id",
    "path",
    "sha256",
    "media_type",
    "frozen_at_utc",
    "immutable",
}
FILE_REF_KEYS = {"object_id", "path", "sha256", "media_type"}
RELEASE_KEYS = {
    "schema",
    "release_id",
    "formation_input_source_sha256",
    "environment_sha256",
    "public",
    "immutable",
}

DATASET_KEYS = {
    "dataset_id",
    "dataset_role",
    "evidence_class",
    "access_mode",
    "source_id",
    "source_uri",
    "acquisition_id",
    "acquisition_cohort_id",
    "raw_root_sha256",
    "raw_root_manifest",
    "provenance",
    "freeze",
    "prior_use",
    "independent_from_dataset_ids",
    "no_reuse_attestation",
    "lots",
    "physical_units",
    "events",
    "expected_objects",
    "raw_objects",
    "derived_objects",
    "missing_objects",
    "deviations",
    "custody_signatures",
}
PROVENANCE_KEYS = {
    "institution",
    "site_id",
    "operator_ids",
    "apparatus_ids",
    "acquired_start_utc",
    "acquired_end_utc",
    "acquisition_record",
    "physical_origin_documented",
}
FREEZE_KEYS = {
    "frozen_before_numerical_access",
    "numerical_access_first_utc",
    "access_log",
    "source_outcomes_known_before_freeze",
    "response_labels_blinded_until_lock",
}
NO_REUSE_KEYS = {
    "no_resplit",
    "no_rerun",
    "no_reanalysis",
    "no_copied_data",
    "no_overlapping_physical_units",
    "no_derived_data_relabel",
}
LOT_KEYS = {"lot_id", "cohort_id", "description"}
PHYSICAL_UNIT_KEYS = {
    "physical_unit_id",
    "unit_kind",
    "lot_id",
    "parent_specimen_id",
    "identity_sha256",
}
EVENT_KEYS = {
    "event_id",
    "physical_unit_id",
    "lot_id",
    "block_id",
    "event_value",
    "route",
    "allocation_concealed",
    "allocation_commitment_sha256",
    "event_generated_at_utc",
    "event_generated_monotonic_ns",
    "before",
    "formation",
    "c_ext_off",
    "g_surface_close",
    "common_hold",
    "reads",
    "raw_object_ids",
    "validity_status",
}
BEFORE_KEYS = {
    "sealed_at_utc",
    "sealed_monotonic_ns",
    "event_absent_when_sealed",
    "raw_object_ids",
}
FORMATION_KEYS = {
    "started_at_utc",
    "started_monotonic_ns",
    "ended_at_utc",
    "ended_monotonic_ns",
    "physical_target_coupling_measured",
    "raw_object_ids",
}
BOUNDARY_KEYS = {
    "observed_at_utc",
    "observed_monotonic_ns",
    "physically_measured",
    "raw_object_ids",
}
G_KEYS = {
    "status",
    "observed_at_utc",
    "observed_monotonic_ns",
    "separate_from_c_ext_off",
    "raw_object_ids",
    "predicate_ids",
}
HOLD_KEYS = {
    "started_at_utc",
    "started_monotonic_ns",
    "ended_at_utc",
    "ended_monotonic_ns",
    "writer_off_throughout",
    "raw_object_ids",
}
READ_KEYS = {
    "read_id",
    "observed_at_utc",
    "observed_monotonic_ns",
    "query_id",
    "mission_read",
    "raw_object_ids",
}
RAW_OBJECT_KEYS = {
    "object_id",
    "path",
    "sha256",
    "size_bytes",
    "media_type",
    "stage",
    "role",
    "event_id",
    "physical_unit_id",
    "instrument_id",
    "clock_id",
    "created_at_utc",
    "monotonic_start_ns",
    "monotonic_end_ns",
    "channels",
    "vendor_original",
    "complete",
}
CHANNEL_KEYS = {
    "name",
    "quantity",
    "unit",
    "uncertainty_value",
    "uncertainty_unit",
    "calibration_id",
}
DERIVED_OBJECT_KEYS = {
    "object_id",
    "path",
    "sha256",
    "size_bytes",
    "media_type",
    "stage",
    "role",
    "event_id",
    "physical_unit_id",
    "parent_object_ids",
    "transformation_object_id",
    "transformation_sha256",
    "environment_sha256",
}
EXPECTED_OBJECT_KEYS = {
    "expected_object_id",
    "event_id",
    "physical_unit_id",
    "role",
}
MISSING_KEYS = {
    "expected_object_id",
    "event_id",
    "physical_unit_id",
    "role",
    "reason",
    "critical",
    "disposition",
}
DEVIATION_KEYS = {
    "deviation_id",
    "event_id",
    "description",
    "critical",
    "allocation_blind",
    "disposition",
}
SIGNATURE_KEYS = {
    "signature_id",
    "role",
    "signer_id",
    "signed_sha256",
    "signed_at_utc",
    "algorithm",
    "signature",
}
PREDICATE_KEYS = {
    "predicate_id",
    "dataset_id",
    "predicate_role",
    "stage",
    "metric_id",
    "acceptance",
    "source_object_ids",
    "rule_object_id",
    "rule_sha256",
    "frozen_before_access",
}
OBSERVED_KEYS = {"value", "unit", "standard_uncertainty"}
ACCEPTANCE_KEYS = {"operator", "lower", "upper", "unit"}
EXECUTION_KEYS = {
    "schema",
    "input_manifest_sha256",
    "release_object_id",
    "release_sha256",
    "source_sha256",
    "environment_sha256",
    "analysis_object_id",
    "analysis_sha256",
    "measurements",
    "invocation",
    "computed_j",
    "evaluated_falsifier",
    "result",
    "stdout",
    "stderr",
}
MEASUREMENT_KEYS = {
    "predicate_id",
    "dataset_id",
    "observed",
    "source_object_ids",
    "measurement_object_id",
    "measurement_sha256",
}
PRODUCT_KEYS = {
    "status",
    "external_physicist_ids",
    "public_urm_release_sha256",
    "used_only_public_instructions",
    "no_private_help",
    "different_real_world_datasets",
    "any_and_every_surface_coverage",
    "evidence_object_ids",
}

DATASET_ROLES = {
    "DEVELOPMENT",
    "VALIDATION",
}
VALIDATION_ROLES = {"VALIDATION"}
EVIDENCE_CLASSES = {"REAL_WORLD_ACTUAL", "SYNTHETIC_TEST_ONLY"}
ACCESS_MODES = {"PROSPECTIVE_NEW_ACQUISITION", "RESPONSE_BLIND_ARCHIVAL_HOLDOUT"}
UNIT_KINDS = {"DEVICE", "SPECIMEN", "SURFACE_INSTANCE"}
ROUTES = {"WRITE", "SHAM"}
VALIDITY = {"VALID", "UNSCOREABLE"}
STAGES = {
    "IDENTITY",
    "CALIBRATION",
    "BEFORE",
    "FORMATION_OR_SHAM",
    "C_EXT_OFF",
    "G_SURFACE_CLOSE",
    "COMMON_HOLD",
    "READ",
    "BLIND_LOCK",
}
RAW_ROLES = {
    "PHYSICAL_IDENTITY",
    "UNCERTAINTY_CALIBRATION",
    "BEFORE_STATE",
    "FORMATION_INPUT",
    "ROUTING_OR_CONTROL",
    "C_EXT_OFF_ASSAY",
    "G_SURFACE_CLOSE_ASSAY",
    "COMMON_HOLD_ENVIRONMENT",
    "READ_OBSERVATION",
    "READ_BACKACTION_ASSAY",
}
REQUIRED_PREDICATE_ROLES = {
    "BEFORE_ABSENCE",
    "WRITE_COUPLING",
    "SHAM_NO_COUPLING",
    "C_EXT_OFF",
    "G_SURFACE_CLOSE",
    "COMMON_HOLD_EQUIVALENCE",
    "EVENT_LEVEL_STATE",
    "MISSION_EVENT_INFORMATION",
    "READ_BACKACTION",
    "K_PREDICTION",
    "M_PREDICTION",
    "ADEQUATE_SENSITIVITY",
}
PREDICATE_OPERATORS = {
    "GE_LOWER_BOUND",
    "LE_UPPER_BOUND",
    "WITHIN_CLOSED_INTERVAL",
}
SI_UNITS = {
    "1",
    "bit",
    "s",
    "ns",
    "m",
    "m^2",
    "m^3",
    "kg",
    "K",
    "A",
    "V",
    "C",
    "J",
    "W",
    "ohm",
    "S",
    "T",
    "Pa",
    "Hz",
    "mol",
    "m/s",
    "K/s",
    "A/s",
    "V/s",
    "ohm/s",
    "J/s",
}
PROTOCOL_ROOT_KEYS = {
    "protocol_id",
    "contract_version",
    "registration",
    "site",
    "scope",
    "unified_dynamics",
    "sequence",
    "causal_design",
    "before",
    "formation_interaction",
    "sham",
    "common_hold",
    "read",
    "analysis",
    "uncertainty_system",
    "transition_assays",
    "data_custody",
    "urm_execution",
    "independent_reproduction",
}
MAPPING_FIELDS = {
    "S", "T", "X", "P", "D_causal", "A", "I_f", "C_psi",
    "C_EXT_OFF", "G_t", "H", "tau_m", "Q", "Z", "J", "U",
    "K_theta", "M_phi", "G_X", "F", "V_URM",
}
PRE_EXECUTION_MAPPING_FIELDS = (MAPPING_FIELDS - {"J", "F", "V_URM"}) | {"F_policy"}
DATASET_MAPPING_FIELDS = PRE_EXECUTION_MAPPING_FIELDS - {"K_theta", "M_phi", "F_policy"}
MAPPING_ROOT_KEYS = {
    "schema",
    "mapping_id",
    "status",
    "protocol_id",
    "mapping_contains_executable_science",
    "surface_specific_source_modification_authorized",
    "adapter_id",
    "adapter_source_sha256",
    "formation_tuple_mapping",
}
MAPPING_ROW_KEYS = {"contract_field", "sources"}
MAPPING_SOURCE_KEYS = {
    "dataset_id",
    "event_id",
    "object_id",
    "path",
    "selector_kind",
    "selector",
    "value_type",
    "quantity",
    "unit",
    "calibration_id",
    "clock_id",
    "coordinate_transform",
}
MAPPING_SELECTOR_KINDS = {"WHOLE_OBJECT", "JSON_POINTER", "CSV_COLUMN", "HDF5_DATASET"}
MAPPING_FORBIDDEN_TOKENS = (
    "<event_uuid>",
    "READ_<time>",
    "*",
    "entries[]",
    "python",
    "javascript",
    "notebook",
    "plugin",
    "import",
    "eval(",
    "exec(",
    "lambda",
)

SHA256 = re.compile(r"^[0-9a-f]{64}$")
UTC = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


class FormationRefusal(ValueError):
    """The supplied formation bundle violates the closed public contract."""


def _refuse(message: str) -> None:
    raise FormationRefusal("FORMATION BUNDLE REFUSES: " + message)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            _refuse(f"duplicate JSON member name {key!r}")
        result[key] = value
    return result


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


def _digest(value: Any, field: str) -> str:
    if not isinstance(value, str) or not SHA256.fullmatch(value):
        _refuse(f"{field} must be a lowercase SHA-256")
    return value


def _boolean(value: Any, field: str) -> bool:
    if not isinstance(value, bool):
        _refuse(f"{field} must be boolean")
    return value


def _finite(value: Any, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        _refuse(f"{field} must be numeric")
    parsed = float(value)
    if not math.isfinite(parsed):
        _refuse(f"{field} must be finite")
    return parsed


def _integer(value: Any, field: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        _refuse(f"{field} must be an integer >= {minimum}")
    return value


def _utc(value: Any, field: str) -> datetime:
    if not isinstance(value, str) or not UTC.fullmatch(value):
        _refuse(f"{field} must be second-resolution UTC")
    try:
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    except ValueError:
        _refuse(f"{field} is not a real UTC date")


def _unique_strings(value: Any, field: str, *, allow_empty: bool = False) -> list[str]:
    if not isinstance(value, list) or (not value and not allow_empty):
        _refuse(f"{field} must be {'a' if not allow_empty else 'an'} unique string list")
    if any(not isinstance(item, str) or not item.strip() for item in value):
        _refuse(f"{field} contains an empty or non-string value")
    if len(value) != len(set(value)):
        _refuse(f"{field} contains duplicates")
    return value


def _safe_relative(root: Path, value: Any, field: str) -> Path:
    raw = _nonempty(value, field)
    if unicodedata.normalize("NFC", raw) != raw:
        _refuse(f"{field} is not Unicode NFC")
    if "\\" in raw:
        _refuse(f"{field} must use forward-slash relative syntax")
    relative = Path(raw)
    if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
        _refuse(f"{field} must be a safe relative path")
    lexical = root / relative
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            _refuse(f"{field} traverses a symlink")
    if not lexical.is_file():
        _refuse(f"{field} does not name a regular file")
    try:
        resolved = lexical.resolve(strict=True)
        resolved.relative_to(root.resolve(strict=True))
    except (OSError, ValueError):
        _refuse(f"{field} escapes the formation bundle")
    return resolved


def _read_json(path: Path, field: str) -> tuple[bytes, Any]:
    try:
        payload = path.read_bytes()
        value = json.loads(payload, object_pairs_hook=_unique_json_object)
    except FormationRefusal:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        _refuse(f"{field} is not readable JSON: {exc}")
    return payload, value


def _file_identity(path: Path) -> tuple[int, int]:
    stat = path.stat()
    return stat.st_dev, stat.st_ino


def _quantity(value: Any, field: str) -> tuple[float, str, float]:
    obj = _closed(value, {"value", "unit", "standard_uncertainty"}, field)
    number = _finite(obj["value"], f"{field}.value")
    unit = _nonempty(obj["unit"], f"{field}.unit")
    if unit not in SI_UNITS:
        _refuse(f"{field}.unit is not in the closed SI unit registry")
    uncertainty = _finite(obj["standard_uncertainty"], f"{field}.standard_uncertainty")
    if uncertainty < 0:
        _refuse(f"{field}.standard_uncertainty must be >= 0")
    return number, unit, uncertainty


def _schema_type_matches(value: Any, expected: str) -> bool:
    if expected == "object":
        return isinstance(value, dict)
    if expected == "array":
        return isinstance(value, list)
    if expected == "string":
        return isinstance(value, str)
    if expected == "boolean":
        return isinstance(value, bool)
    if expected == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if expected == "null":
        return value is None
    _refuse(f"trusted V001 schema uses unsupported type {expected!r}")


def _schema_resolve(root_schema: dict[str, Any], reference: str) -> dict[str, Any]:
    if not reference.startswith("#/"):
        _refuse("trusted V001 schema contains a nonlocal reference")
    node: Any = root_schema
    for token in reference[2:].split("/"):
        token = token.replace("~1", "/").replace("~0", "~")
        if not isinstance(node, dict) or token not in node:
            _refuse(f"trusted V001 schema has unresolved reference {reference!r}")
        node = node[token]
    if not isinstance(node, dict):
        _refuse(f"trusted V001 schema reference {reference!r} is not an object")
    return node


def _schema_matches(
    value: Any,
    schema: dict[str, Any],
    root_schema: dict[str, Any],
    field: str,
) -> bool:
    try:
        _schema_validate(value, schema, root_schema, field)
    except FormationRefusal:
        return False
    return True


def _schema_validate(
    value: Any,
    schema: dict[str, Any],
    root_schema: dict[str, Any],
    field: str,
) -> None:
    if "$ref" in schema:
        _schema_validate(value, _schema_resolve(root_schema, schema["$ref"]), root_schema, field)
    if "allOf" in schema:
        for index, branch in enumerate(schema["allOf"]):
            _schema_validate(value, branch, root_schema, f"{field}.allOf[{index}]")
    if "oneOf" in schema:
        matches = sum(
            _schema_matches(value, branch, root_schema, field)
            for branch in schema["oneOf"]
        )
        if matches != 1:
            _refuse(f"{field} matches {matches} branches of a V001 oneOf")
    if "not" in schema and _schema_matches(value, schema["not"], root_schema, field):
        _refuse(f"{field} violates a V001 not condition")
    if "if" in schema and _schema_matches(value, schema["if"], root_schema, field):
        if "then" in schema:
            _schema_validate(value, schema["then"], root_schema, field)
    if "const" in schema and value != schema["const"]:
        _refuse(f"{field} differs from its V001 const")
    if "enum" in schema and value not in schema["enum"]:
        _refuse(f"{field} is absent from its V001 enum")
    if "type" in schema:
        expected_types = schema["type"] if isinstance(schema["type"], list) else [schema["type"]]
        if not any(_schema_type_matches(value, expected) for expected in expected_types):
            _refuse(f"{field} has the wrong V001 type")
    if isinstance(value, dict):
        required = schema.get("required", [])
        if any(key not in value for key in required):
            _refuse(f"{field} omits a V001 required key")
        properties = schema.get("properties", {})
        for key, child in properties.items():
            if key in value:
                _schema_validate(value[key], child, root_schema, f"{field}.{key}")
        if schema.get("additionalProperties") is False:
            extras = set(value) - set(properties)
            if extras:
                _refuse(f"{field} has V001-forbidden keys {sorted(extras)}")
    if isinstance(value, list):
        if len(value) < schema.get("minItems", 0):
            _refuse(f"{field} has too few V001 items")
        if schema.get("uniqueItems"):
            normalized = [json.dumps(item, sort_keys=True, separators=(",", ":")) for item in value]
            if len(normalized) != len(set(normalized)):
                _refuse(f"{field} violates V001 uniqueItems")
        if isinstance(schema.get("items"), dict):
            for index, item in enumerate(value):
                _schema_validate(item, schema["items"], root_schema, f"{field}[{index}]")
    if isinstance(value, str):
        if len(value) < schema.get("minLength", 0):
            _refuse(f"{field} is shorter than V001 minLength")
        if "pattern" in schema and re.search(schema["pattern"], value) is None:
            _refuse(f"{field} does not match its V001 pattern")
        if schema.get("format") == "date-time":
            _utc(value, field)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if not math.isfinite(float(value)):
            _refuse(f"{field} is nonfinite")
        if "exclusiveMinimum" in schema and not value > schema["exclusiveMinimum"]:
            _refuse(f"{field} violates V001 exclusiveMinimum")


@dataclass(frozen=True)
class PredicateScore:
    predicate_id: str
    dataset_id: str
    predicate_role: str
    passed: bool | None


@dataclass(frozen=True)
class FormationInput:
    manifest_path: Path
    manifest_sha256: str
    manifest: dict[str, Any]
    dataset_input_complete: tuple[tuple[str, bool], ...]
    dataset_object_ids: tuple[tuple[str, tuple[str, ...]], ...]
    referenced_hashes: tuple[tuple[str, str], ...]

    def certificate(self) -> dict[str, Any]:
        roles: dict[str, list[str]] = {role: [] for role in sorted(DATASET_ROLES)}
        evidence: dict[str, str] = {}
        for dataset in self.manifest["datasets"]:
            roles[dataset["dataset_role"]].append(dataset["dataset_id"])
            evidence[dataset["dataset_id"]] = dataset["evidence_class"]
        input_complete = dict(self.dataset_input_complete)
        development_ids = roles["DEVELOPMENT"]
        validation_ids = roles["VALIDATION"]
        pair_ids = development_ids + validation_ids
        pair_present = bool(development_ids and validation_ids)
        real_world_pair = pair_present and all(
            evidence[dataset_id] == "REAL_WORLD_ACTUAL" for dataset_id in pair_ids
        )
        all_pair_complete = pair_present and all(
            input_complete.get(dataset_id, False) for dataset_id in pair_ids
        )
        eligible = bool(real_world_pair and all_pair_complete)
        return {
            "schema": CERTIFICATE_SCHEMA,
            "input_contract_sha256": sha256_file(Path(__file__).resolve()),
            "manifest_sha256": self.manifest_sha256,
            "bundle_id": self.manifest["bundle_id"],
            "protocol_id": self.manifest["protocol_id"],
            "base_contract_version": BASE_CONTRACT_VERSION,
            "validation_extension_version": VALIDATION_EXTENSION_VERSION,
            "base_protocol_semantics_modified": False,
            "dataset_roles": roles,
            "dataset_evidence_class": evidence,
            "dataset_input_complete": input_complete,
            "dataset_independence": "EXACT_ID_HASH_PATH_INODE_AND_CONTENT_DISJOINTNESS_VERIFIED",
            "development_data_used_for_validation": False,
            "validation_dataset_eligible": eligible,
            "generic_predicates_scored": False,
            "scientific_verdict": "NONE_NOT_SCORED",
            "scientific_validation_authorized": False,
            "record_formation_claim_authorized": False,
            "scoped_formation_result_authorized": False,
            "universal_claim_authorized": False,
            "gravity_claim_authorized": False,
            "program_completion_authorized": False,
            "product_reproduction": {
                "status": "NOT_DEMONSTRATED",
                "separate_from_dataset_validation": True,
                "attested": False,
            },
            "referenced_hashes": [
                {"path": path, "sha256": digest}
                for path, digest in self.referenced_hashes
            ],
        }


@dataclass(frozen=True)
class FormationExecution:
    formation_input: FormationInput
    execution_path: Path
    execution_sha256: str
    predicate_scores: tuple[PredicateScore, ...]
    referenced_hashes: tuple[tuple[str, str], ...]

    def certificate(self) -> dict[str, Any]:
        certificate = self.formation_input.certificate()
        validation_ids = {
            dataset["dataset_id"]
            for dataset in self.formation_input.manifest["datasets"]
            if dataset["dataset_role"] == "VALIDATION"
        }
        validation_scores = [
            score for score in self.predicate_scores if score.dataset_id in validation_ids
        ]
        if any(score.passed is None for score in validation_scores):
            generic_status = "UNSCOREABLE"
        elif validation_scores and all(score.passed for score in validation_scores):
            generic_status = "ALL_FROZEN_GENERIC_PREDICATES_PASS"
        else:
            generic_status = "ONE_OR_MORE_FROZEN_GENERIC_PREDICATES_FAIL"
        certificate.update({
            "schema": EXECUTION_CERTIFICATE_SCHEMA,
            "execution_sha256": self.execution_sha256,
            "generic_predicates_scored": True,
            "generic_predicate_status": generic_status,
            "predicate_scores": [
                {
                    "predicate_id": score.predicate_id,
                    "dataset_id": score.dataset_id,
                    "predicate_role": score.predicate_role,
                    "status": (
                        "PASS" if score.passed is True
                        else "FAIL" if score.passed is False
                        else "UNSCOREABLE"
                    ),
                }
                for score in self.predicate_scores
            ],
            "scientific_verdict": "NONE_NOT_SCORED",
            "scientific_validation_authorized": False,
            "record_formation_claim_authorized": False,
            "scoped_formation_result_authorized": False,
            "universal_claim_authorized": False,
            "gravity_claim_authorized": False,
            "program_completion_authorized": False,
        })
        certificate["referenced_hashes"] = certificate["referenced_hashes"] + [
            {"path": path, "sha256": digest}
            for path, digest in self.referenced_hashes
        ]
        return certificate


def _score_predicate(spec: dict[str, Any], observed: Any, field: str) -> bool:
    value, unit, uncertainty = _quantity(observed, f"{field}.observed")
    acceptance = _closed(spec["acceptance"], ACCEPTANCE_KEYS, f"{field}.acceptance")
    operator = acceptance["operator"]
    if operator not in PREDICATE_OPERATORS:
        _refuse(f"{field}.acceptance.operator is not registered")
    threshold_unit = _nonempty(acceptance["unit"], f"{field}.acceptance.unit")
    if threshold_unit != unit:
        _refuse(f"{field} observed and acceptance units differ")
    lower = acceptance["lower"]
    upper = acceptance["upper"]
    if operator == "GE_LOWER_BOUND":
        if lower is None or upper is not None:
            _refuse(f"{field} GE_LOWER_BOUND requires lower and null upper")
        return value - uncertainty >= _finite(lower, f"{field}.acceptance.lower")
    if operator == "LE_UPPER_BOUND":
        if lower is not None or upper is None:
            _refuse(f"{field} LE_UPPER_BOUND requires null lower and upper")
        return value + uncertainty <= _finite(upper, f"{field}.acceptance.upper")
    if lower is None or upper is None:
        _refuse(f"{field} WITHIN_CLOSED_INTERVAL requires lower and upper")
    low = _finite(lower, f"{field}.acceptance.lower")
    high = _finite(upper, f"{field}.acceptance.upper")
    if low > high:
        _refuse(f"{field} acceptance interval is reversed")
    return value - uncertainty >= low and value + uncertainty <= high


def load_formation_input(manifest_path: str | Path) -> FormationInput:
    """Load and structurally validate one V002 formation input bundle."""
    manifest_path = Path(manifest_path)
    if manifest_path.is_symlink():
        _refuse("submitted manifest must not be a symlink")
    if manifest_path.parent.is_symlink():
        _refuse("submitted manifest bundle directory must not be a symlink")
    try:
        manifest_bytes = manifest_path.read_bytes()
        manifest = json.loads(manifest_bytes, object_pairs_hook=_unique_json_object)
    except FormationRefusal:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        _refuse(f"manifest is not readable JSON: {exc}")
    manifest = _closed(manifest, ROOT_KEYS, "manifest")
    if manifest["schema"] != SCHEMA:
        _refuse(f"unsupported schema {manifest['schema']!r}")
    _nonempty(manifest["bundle_id"], "bundle_id")
    _nonempty(manifest["protocol_id"], "protocol_id")
    if manifest["base_contract_version"] != BASE_CONTRACT_VERSION:
        _refuse("base_contract_version must preserve FORMATION_PROTOCOL_V001")
    if manifest["validation_extension_version"] != VALIDATION_EXTENSION_VERSION:
        _refuse("validation_extension_version is not registered")
    try:
        resolved_manifest = manifest_path.resolve(strict=True)
        root = manifest_path.absolute().parent.resolve(strict=True)
        resolved_manifest.relative_to(root)
    except (OSError, ValueError):
        _refuse("submitted manifest escapes its lexical bundle root")

    claim = _closed(manifest["claim_scope"], CLAIM_SCOPE_KEYS, "claim_scope")
    for key in ("claim_id", "surface_scope", "formation_mode", "provenance"):
        _nonempty(claim[key], f"claim_scope.{key}")
    if claim["requested_result"] != "SCOPED_FORMATION_DATASET_VALIDATION":
        _refuse("claim_scope.requested_result is not registered")
    if _boolean(claim["universal_claim_requested"], "claim_scope.universal_claim_requested"):
        _refuse("this scoped contract cannot request a universal claim")
    if _boolean(claim["gravity_claim_requested"], "claim_scope.gravity_claim_requested"):
        _refuse("this formation contract cannot request a gravity claim")
    mission_value, mission_unit, mission_uncertainty = _quantity(
        claim["mission_time"], "claim_scope.mission_time"
    )
    if mission_unit != "s" or mission_value <= 0:
        _refuse("claim_scope.mission_time must be positive seconds")

    seen_paths: dict[Path, str] = {}
    seen_identities: dict[tuple[int, int], str] = {
        _file_identity(resolved_manifest): "manifest"
    }
    seen_file_object_ids: dict[str, str] = {}
    referenced_hashes: list[tuple[str, str]] = []

    def register_file(ref: Any, field: str, *, frozen: bool = False) -> tuple[dict[str, Any], Path]:
        keys = FROZEN_REF_KEYS if frozen else FILE_REF_KEYS
        item = _closed(ref, keys, field)
        object_id = _nonempty(item["object_id"], f"{field}.object_id")
        if object_id in seen_file_object_ids:
            _refuse(f"{field}.object_id duplicates {seen_file_object_ids[object_id]}")
        seen_file_object_ids[object_id] = field
        _nonempty(item["media_type"], f"{field}.media_type")
        expected = _digest(item["sha256"], f"{field}.sha256")
        path = _safe_relative(root, item["path"], f"{field}.path")
        if path.stat().st_size == 0:
            _refuse(f"{field}.path is a zero-byte placeholder")
        if path in seen_paths:
            _refuse(f"{field} duplicates resolved path used by {seen_paths[path]}")
        identity = _file_identity(path)
        if identity in seen_identities:
            _refuse(f"{field} hard-link aliases {seen_identities[identity]}")
        if sha256_file(path) != expected:
            _refuse(f"{field} hash mismatch")
        seen_paths[path] = field
        seen_identities[identity] = field
        referenced_hashes.append((str(path.relative_to(root)), expected))
        if frozen:
            _utc(item["frozen_at_utc"], f"{field}.frozen_at_utc")
            if _boolean(item["immutable"], f"{field}.immutable") is not True:
                _refuse(f"{field} must be immutable")
        return item, path

    frozen = _closed(manifest["frozen_objects"], FROZEN_OBJECT_NAMES, "frozen_objects")
    frozen_items: dict[str, dict[str, Any]] = {}
    frozen_paths: dict[str, Path] = {}
    frozen_ids: set[str] = set()
    for name in sorted(FROZEN_OBJECT_NAMES):
        item, path = register_file(frozen[name], f"frozen_objects.{name}", frozen=True)
        if item["object_id"] in frozen_ids:
            _refuse("frozen object IDs must be unique")
        frozen_ids.add(item["object_id"])
        frozen_items[name] = item
        frozen_paths[name] = path

    public_source_hash = sha256_file(Path(__file__).resolve())
    if frozen_items["urm_source"]["sha256"] != public_source_hash:
        _refuse("frozen URM source is not this installed generic formation-input source")
    _, environment = _read_json(frozen_paths["urm_environment"], "frozen URM environment")
    if not isinstance(environment, dict):
        _refuse("frozen URM environment must be a JSON object")
    _, release = _read_json(frozen_paths["urm_release"], "frozen URM release")
    release = _closed(release, RELEASE_KEYS, "frozen URM release root")
    if release["schema"] != "WAC_FORMATION_URM_RELEASE_V002":
        _refuse("frozen URM release schema is not registered")
    _nonempty(release["release_id"], "frozen URM release.release_id")
    if _digest(
        release["formation_input_source_sha256"],
        "frozen URM release.formation_input_source_sha256",
    ) != public_source_hash:
        _refuse("frozen URM release source hash differs from the installed public source")
    if _digest(
        release["environment_sha256"], "frozen URM release.environment_sha256"
    ) != frozen_items["urm_environment"]["sha256"]:
        _refuse("frozen URM release environment hash mismatch")
    if _boolean(release["public"], "frozen URM release.public") is not True:
        _refuse("frozen URM release is not public")
    if _boolean(release["immutable"], "frozen URM release.immutable") is not True:
        _refuse("frozen URM release is not immutable")

    _, protocol = _read_json(frozen_paths["protocol"], "frozen protocol")
    protocol = _closed(protocol, PROTOCOL_ROOT_KEYS, "frozen protocol root")
    base_schema_path = (
        Path(__file__).resolve().parent.parent
        / "LANE_T53_D_FORMATION_GATE"
        / "FORMATION_PROTOCOL_V001.schema.json"
    )
    if not base_schema_path.is_file() or sha256_file(base_schema_path) != BASE_SCHEMA_SHA256:
        _refuse("sealed FORMATION_PROTOCOL_V001 schema is absent or hash-mismatched")
    _, base_schema = _read_json(base_schema_path, "sealed FORMATION_PROTOCOL_V001 schema")
    if not isinstance(base_schema, dict):
        _refuse("sealed FORMATION_PROTOCOL_V001 schema root is not an object")
    _schema_validate(protocol, base_schema, base_schema, "frozen protocol")
    if protocol["protocol_id"] != manifest["protocol_id"]:
        _refuse("frozen protocol_id differs from the bundle")
    if protocol["contract_version"] != BASE_CONTRACT_VERSION:
        _refuse("frozen protocol does not preserve FORMATION_PROTOCOL_V001")
    _, mapping = _read_json(frozen_paths["declarative_mapping"], "declarative mapping")
    mapping = _closed(mapping, MAPPING_ROOT_KEYS, "declarative mapping root")
    if mapping["schema"] != "WAC_FORMATION_NORMALIZED_MAPPING_V002":
        _refuse("declarative mapping is not the normalized V002 machine contract")
    _nonempty(mapping["mapping_id"], "declarative mapping.mapping_id")
    if mapping["status"] != "NORMALIZED_MACHINE_RESOLVABLE":
        _refuse("declarative mapping remains unresolved")
    if mapping["protocol_id"] != manifest["protocol_id"]:
        _refuse("declarative mapping protocol_id differs from the bundle")
    if mapping["mapping_contains_executable_science"] is not False:
        _refuse("declarative mapping must explicitly contain no executable science")
    if mapping["surface_specific_source_modification_authorized"] is not False:
        _refuse("declarative mapping must forbid surface-specific source modification")
    _nonempty(mapping["adapter_id"], "declarative mapping.adapter_id")
    if _digest(
        mapping["adapter_source_sha256"], "declarative mapping.adapter_source_sha256"
    ) != frozen_items["urm_source"]["sha256"]:
        _refuse("declarative mapping adapter source is not the frozen generic URM source")
    rows = mapping["formation_tuple_mapping"]
    if not isinstance(rows, list):
        _refuse("declarative mapping lacks formation_tuple_mapping")
    mapped_fields: set[str] = set()
    mapping_sources: list[tuple[str, dict[str, Any]]] = []
    for index, raw_row in enumerate(rows):
        row_field = f"declarative mapping row {index}"
        row = _closed(raw_row, MAPPING_ROW_KEYS, row_field)
        contract_field = _nonempty(row["contract_field"], f"{row_field}.contract_field")
        if contract_field in mapped_fields:
            _refuse(f"declarative mapping duplicates field {contract_field}")
        mapped_fields.add(contract_field)
        sources = row["sources"]
        if not isinstance(sources, list) or not sources:
            _refuse(f"{row_field}.sources must be nonempty")
        for source_index, raw_source in enumerate(sources):
            source_field = f"{row_field}.sources[{source_index}]"
            source = _closed(raw_source, MAPPING_SOURCE_KEYS, source_field)
            for nullable_id in ("dataset_id", "event_id", "calibration_id", "clock_id"):
                if source[nullable_id] is not None:
                    _nonempty(source[nullable_id], f"{source_field}.{nullable_id}")
            object_id = _nonempty(source["object_id"], f"{source_field}.object_id")
            path_text = _nonempty(source["path"], f"{source_field}.path")
            if unicodedata.normalize("NFC", path_text) != path_text or "\\" in path_text:
                _refuse(f"{source_field}.path is not a normalized relative path")
            path_value = Path(path_text)
            if path_value.is_absolute() or any(part in {"", ".", ".."} for part in path_value.parts):
                _refuse(f"{source_field}.path is not a safe relative path")
            if source["selector_kind"] not in MAPPING_SELECTOR_KINDS:
                _refuse(f"{source_field}.selector_kind is not registered")
            selector = _nonempty(source["selector"], f"{source_field}.selector")
            quantity = _nonempty(source["quantity"], f"{source_field}.quantity")
            value_type = source["value_type"]
            if value_type not in {"IDENTITY", "CATEGORICAL", "BOOLEAN", "INTEGER", "FLOAT", "TIMESTAMP"}:
                _refuse(f"{source_field}.value_type is not registered")
            transform = source["coordinate_transform"]
            if transform not in {"IDENTITY", "REGISTERED_LINEAR_SI", "REGISTERED_ENUMERATION"}:
                _refuse(f"{source_field}.coordinate_transform is not registered")
            text_fields = " ".join((path_text, selector, quantity, object_id, transform)).lower()
            if any(token.lower() in text_fields for token in MAPPING_FORBIDDEN_TOKENS):
                _refuse(f"{source_field} contains executable or unresolved mapping syntax")
            if value_type in {"INTEGER", "FLOAT", "TIMESTAMP"}:
                unit = _nonempty(source["unit"], f"{source_field}.unit")
                if unit not in SI_UNITS:
                    _refuse(f"{source_field}.unit is not in the closed SI registry")
                _nonempty(source["calibration_id"], f"{source_field}.calibration_id")
                _nonempty(source["clock_id"], f"{source_field}.clock_id")
            else:
                if source["unit"] is not None or source["calibration_id"] is not None or source["clock_id"] is not None:
                    _refuse(f"{source_field} nonnumeric mapping must use null unit/calibration/clock")
            mapping_sources.append((contract_field, source))
    if mapped_fields != PRE_EXECUTION_MAPPING_FIELDS:
        _refuse("declarative mapping does not cover exactly the pre-execution formation tuple")

    datasets = manifest["datasets"]
    if not isinstance(datasets, list) or not datasets:
        _refuse("datasets must be a nonempty list")
    dataset_ids: set[str] = set()
    dataset_by_id: dict[str, dict[str, Any]] = {}
    dataset_times: dict[str, datetime] = {}
    dataset_objects: dict[str, set[str]] = {}
    dataset_scoreable: dict[str, bool] = {}
    g_predicate_refs: list[tuple[str, str, tuple[str, ...]]] = []

    disjoint_values: dict[str, dict[str, str]] = {
        "source_id": {},
        "source_uri": {},
        "acquisition_id": {},
        "acquisition_cohort_id": {},
        "acquisition_record_sha256": {},
        "access_log_sha256": {},
        "raw_root_sha256": {},
        "lot_id": {},
        "lot_cohort_id": {},
        "parent_specimen_id": {},
        "physical_unit_id": {},
        "physical_identity_sha256": {},
        "event_id": {},
        "expected_object_id": {},
        "object_id": {},
        "raw_object_id": {},
        "derived_object_id": {},
        "data_content_sha256": {},
    }

    def disjoint(kind: str, value: str, dataset_id: str) -> None:
        previous = disjoint_values[kind].get(value)
        if previous is not None and previous != dataset_id:
            _refuse(f"dataset overlap: {kind} {value!r} reused by {previous} and {dataset_id}")
        disjoint_values[kind][value] = dataset_id

    for d_index, raw_dataset in enumerate(datasets):
        field = f"datasets[{d_index}]"
        dataset = _closed(raw_dataset, DATASET_KEYS, field)
        dataset_id = _nonempty(dataset["dataset_id"], f"{field}.dataset_id")
        if dataset_id in dataset_ids:
            _refuse(f"duplicate dataset_id {dataset_id}")
        dataset_ids.add(dataset_id)
        dataset_by_id[dataset_id] = dataset
        role = dataset["dataset_role"]
        if role not in DATASET_ROLES:
            _refuse(f"{field}.dataset_role is not registered")
        if dataset["evidence_class"] not in EVIDENCE_CLASSES:
            _refuse(f"{field}.evidence_class is not registered")
        if dataset["access_mode"] not in ACCESS_MODES:
            _refuse(f"{field}.access_mode is not registered")
        source_id = _nonempty(dataset["source_id"], f"{field}.source_id")
        source_uri = _nonempty(dataset["source_uri"], f"{field}.source_uri")
        if not source_uri.startswith(("https://", "doi:", "urn:")):
            _refuse(f"{field}.source_uri must use https, doi, or urn")
        acquisition_id = _nonempty(dataset["acquisition_id"], f"{field}.acquisition_id")
        cohort_id = _nonempty(dataset["acquisition_cohort_id"], f"{field}.acquisition_cohort_id")
        root_hash = _digest(dataset["raw_root_sha256"], f"{field}.raw_root_sha256")
        for kind, value in (
            ("source_id", source_id),
            ("source_uri", source_uri),
            ("acquisition_id", acquisition_id),
            ("acquisition_cohort_id", cohort_id),
            ("raw_root_sha256", root_hash),
        ):
            disjoint(kind, value, dataset_id)
        root_ref, root_path = register_file(dataset["raw_root_manifest"], f"{field}.raw_root_manifest")
        if root_ref["sha256"] != root_hash:
            _refuse(f"{field}.raw_root_sha256 differs from raw_root_manifest hash")

        provenance = _closed(dataset["provenance"], PROVENANCE_KEYS, f"{field}.provenance")
        for key in ("institution", "site_id"):
            _nonempty(provenance[key], f"{field}.provenance.{key}")
        _unique_strings(provenance["operator_ids"], f"{field}.provenance.operator_ids")
        _unique_strings(provenance["apparatus_ids"], f"{field}.provenance.apparatus_ids")
        start = _utc(provenance["acquired_start_utc"], f"{field}.provenance.acquired_start_utc")
        end = _utc(provenance["acquired_end_utc"], f"{field}.provenance.acquired_end_utc")
        if end < start:
            _refuse(f"{field}.provenance acquisition interval is reversed")
        acquisition_record, _ = register_file(
            provenance["acquisition_record"], f"{field}.provenance.acquisition_record"
        )
        disjoint("acquisition_record_sha256", acquisition_record["sha256"], dataset_id)
        _boolean(provenance["physical_origin_documented"], f"{field}.provenance.physical_origin_documented")

        freeze = _closed(dataset["freeze"], FREEZE_KEYS, f"{field}.freeze")
        access_time = _utc(
            freeze["numerical_access_first_utc"],
            f"{field}.freeze.numerical_access_first_utc",
        )
        dataset_times[dataset_id] = access_time
        access_log, _ = register_file(freeze["access_log"], f"{field}.freeze.access_log")
        disjoint("access_log_sha256", access_log["sha256"], dataset_id)
        frozen_before = _boolean(
            freeze["frozen_before_numerical_access"],
            f"{field}.freeze.frozen_before_numerical_access",
        )
        known_before = _boolean(
            freeze["source_outcomes_known_before_freeze"],
            f"{field}.freeze.source_outcomes_known_before_freeze",
        )
        response_blinded = _boolean(
            freeze["response_labels_blinded_until_lock"],
            f"{field}.freeze.response_labels_blinded_until_lock",
        )
        if role in VALIDATION_ROLES:
            if dataset["evidence_class"] != "REAL_WORLD_ACTUAL":
                _refuse(f"{field} validation role requires REAL_WORLD_ACTUAL evidence")
            if not frozen_before or known_before or not response_blinded:
                _refuse(f"{field} validation data were not frozen before numerical access")
            if dataset["prior_use"] != "UNUSED_BEFORE_VALIDATION":
                _refuse(f"{field} validation dataset was previously used")
            if provenance["physical_origin_documented"] is not True:
                _refuse(f"{field} validation dataset lacks documented physical origin")
            for name, item in frozen_items.items():
                if _utc(item["frozen_at_utc"], f"frozen_objects.{name}.frozen_at_utc") >= access_time:
                    _refuse(f"{field} opened before frozen object {name}")
        else:
            if dataset["prior_use"] != "DEVELOPMENT_ONLY":
                _refuse(f"{field} DEVELOPMENT dataset prior_use must be DEVELOPMENT_ONLY")

        no_reuse = _closed(
            dataset["no_reuse_attestation"], NO_REUSE_KEYS, f"{field}.no_reuse_attestation"
        )
        for key, value in no_reuse.items():
            if _boolean(value, f"{field}.no_reuse_attestation.{key}") is not True:
                _refuse(f"{field} does not attest {key}")

        lots = dataset["lots"]
        if not isinstance(lots, list) or not lots:
            _refuse(f"{field}.lots must be nonempty")
        lot_ids: set[str] = set()
        for index, raw_lot in enumerate(lots):
            lot = _closed(raw_lot, LOT_KEYS, f"{field}.lots[{index}]")
            lot_id = _nonempty(lot["lot_id"], f"{field}.lots[{index}].lot_id")
            if lot_id in lot_ids:
                _refuse(f"{field}.lots duplicates {lot_id}")
            lot_ids.add(lot_id)
            disjoint("lot_id", lot_id, dataset_id)
            lot_cohort = _nonempty(lot["cohort_id"], f"{field}.lots[{index}].cohort_id")
            disjoint("lot_cohort_id", lot_cohort, dataset_id)
            _nonempty(lot["description"], f"{field}.lots[{index}].description")

        physical_units = dataset["physical_units"]
        if not isinstance(physical_units, list) or not physical_units:
            _refuse(f"{field}.physical_units must be nonempty")
        physical_unit_ids: set[str] = set()
        for index, raw_unit in enumerate(physical_units):
            unit = _closed(raw_unit, PHYSICAL_UNIT_KEYS, f"{field}.physical_units[{index}]")
            unit_id = _nonempty(
                unit["physical_unit_id"], f"{field}.physical_units[{index}].physical_unit_id"
            )
            if unit_id in physical_unit_ids:
                _refuse(f"{field}.physical_units duplicates {unit_id}")
            physical_unit_ids.add(unit_id)
            disjoint("physical_unit_id", unit_id, dataset_id)
            if unit["unit_kind"] not in UNIT_KINDS:
                _refuse(f"{field}.physical_units[{index}].unit_kind is not registered")
            if unit["lot_id"] not in lot_ids:
                _refuse(f"{field}.physical_units[{index}] references an unknown lot")
            if unit["parent_specimen_id"] is not None:
                parent_specimen = _nonempty(
                    unit["parent_specimen_id"],
                    f"{field}.physical_units[{index}].parent_specimen_id",
                )
                disjoint("parent_specimen_id", parent_specimen, dataset_id)
            identity_hash = _digest(
                unit["identity_sha256"], f"{field}.physical_units[{index}].identity_sha256"
            )
            disjoint("physical_identity_sha256", identity_hash, dataset_id)

        raw_objects = dataset["raw_objects"]
        if not isinstance(raw_objects, list) or not raw_objects:
            _refuse(f"{field}.raw_objects must be nonempty")
        raw_by_id: dict[str, dict[str, Any]] = {}
        object_ids: set[str] = set()
        raw_root_rows: list[dict[str, str]] = []
        roles_present: set[str] = set()
        all_raw_complete = True
        for index, raw_item in enumerate(raw_objects):
            obj_field = f"{field}.raw_objects[{index}]"
            item = _closed(raw_item, RAW_OBJECT_KEYS, obj_field)
            object_id = _nonempty(item["object_id"], f"{obj_field}.object_id")
            if object_id in object_ids:
                _refuse(f"{field} duplicates object_id {object_id}")
            object_ids.add(object_id)
            raw_by_id[object_id] = item
            disjoint("object_id", object_id, dataset_id)
            disjoint("raw_object_id", object_id, dataset_id)
            expected_hash = _digest(item["sha256"], f"{obj_field}.sha256")
            disjoint("data_content_sha256", expected_hash, dataset_id)
            item_ref = {
                "object_id": object_id,
                "path": item["path"],
                "sha256": expected_hash,
                "media_type": item["media_type"],
            }
            _, path = register_file(item_ref, obj_field)
            if path.stat().st_size != _integer(item["size_bytes"], f"{obj_field}.size_bytes", 1):
                _refuse(f"{obj_field}.size_bytes mismatch")
            if item["stage"] not in STAGES:
                _refuse(f"{obj_field}.stage is not registered")
            if item["role"] not in RAW_ROLES:
                _refuse(f"{obj_field}.role is not registered")
            roles_present.add(item["role"])
            if item["event_id"] is not None:
                _nonempty(item["event_id"], f"{obj_field}.event_id")
            if item["physical_unit_id"] is not None and item["physical_unit_id"] not in physical_unit_ids:
                _refuse(f"{obj_field} references an unknown physical unit")
            _nonempty(item["instrument_id"], f"{obj_field}.instrument_id")
            _nonempty(item["clock_id"], f"{obj_field}.clock_id")
            _utc(item["created_at_utc"], f"{obj_field}.created_at_utc")
            mono_start = _integer(item["monotonic_start_ns"], f"{obj_field}.monotonic_start_ns")
            mono_end = _integer(item["monotonic_end_ns"], f"{obj_field}.monotonic_end_ns")
            if mono_end < mono_start:
                _refuse(f"{obj_field} monotonic interval is reversed")
            channels = item["channels"]
            if not isinstance(channels, list) or not channels:
                _refuse(f"{obj_field}.channels must be nonempty")
            channel_names: set[str] = set()
            for c_index, raw_channel in enumerate(channels):
                channel_field = f"{obj_field}.channels[{c_index}]"
                channel = _closed(raw_channel, CHANNEL_KEYS, channel_field)
                name = _nonempty(channel["name"], f"{channel_field}.name")
                if name in channel_names:
                    _refuse(f"{obj_field} duplicates channel {name}")
                channel_names.add(name)
                _nonempty(channel["quantity"], f"{channel_field}.quantity")
                unit = _nonempty(channel["unit"], f"{channel_field}.unit")
                uncertainty_unit = _nonempty(
                    channel["uncertainty_unit"], f"{channel_field}.uncertainty_unit"
                )
                if unit not in SI_UNITS or uncertainty_unit != unit:
                    _refuse(f"{channel_field} must use matching registered SI units")
                uncertainty = _finite(
                    channel["uncertainty_value"], f"{channel_field}.uncertainty_value"
                )
                if uncertainty < 0:
                    _refuse(f"{channel_field}.uncertainty_value must be >= 0")
                _nonempty(channel["calibration_id"], f"{channel_field}.calibration_id")
            _boolean(item["vendor_original"], f"{obj_field}.vendor_original")
            all_raw_complete &= _boolean(item["complete"], f"{obj_field}.complete")
            raw_root_rows.append({"object_id": object_id, "sha256": expected_hash})

        derived_objects = dataset["derived_objects"]
        if not isinstance(derived_objects, list):
            _refuse(f"{field}.derived_objects must be a list")
        derived_by_id: dict[str, dict[str, Any]] = {}
        derived_root_rows: list[dict[str, str]] = []
        for index, raw_item in enumerate(derived_objects):
            obj_field = f"{field}.derived_objects[{index}]"
            item = _closed(raw_item, DERIVED_OBJECT_KEYS, obj_field)
            object_id = _nonempty(item["object_id"], f"{obj_field}.object_id")
            if object_id in object_ids:
                _refuse(f"{field} duplicates object_id {object_id}")
            object_ids.add(object_id)
            derived_by_id[object_id] = item
            disjoint("object_id", object_id, dataset_id)
            disjoint("derived_object_id", object_id, dataset_id)
            expected_hash = _digest(item["sha256"], f"{obj_field}.sha256")
            disjoint("data_content_sha256", expected_hash, dataset_id)
            item_ref = {
                "object_id": object_id,
                "path": item["path"],
                "sha256": expected_hash,
                "media_type": item["media_type"],
            }
            _, path = register_file(item_ref, obj_field)
            if path.stat().st_size != _integer(item["size_bytes"], f"{obj_field}.size_bytes", 1):
                _refuse(f"{obj_field}.size_bytes mismatch")
            if item["stage"] not in STAGES:
                _refuse(f"{obj_field}.stage is not registered")
            _nonempty(item["role"], f"{obj_field}.role")
            if item["event_id"] is not None:
                _nonempty(item["event_id"], f"{obj_field}.event_id")
            if item["physical_unit_id"] is not None and item["physical_unit_id"] not in physical_unit_ids:
                _refuse(f"{obj_field} references an unknown physical unit")
            parents = _unique_strings(item["parent_object_ids"], f"{obj_field}.parent_object_ids")
            if any(parent not in raw_by_id for parent in parents):
                _refuse(f"{obj_field} parent crosses a dataset or is not raw")
            if item["transformation_object_id"] != frozen_items["analysis"]["object_id"]:
                _refuse(f"{obj_field} transformation is not the frozen analysis object")
            if _digest(item["transformation_sha256"], f"{obj_field}.transformation_sha256") != frozen_items["analysis"]["sha256"]:
                _refuse(f"{obj_field} transformation hash mismatch")
            _digest(item["environment_sha256"], f"{obj_field}.environment_sha256")
            derived_root_rows.append({"object_id": object_id, "sha256": expected_hash})

        _, root_inventory = _read_json(root_path, f"{field}.raw_root_manifest")
        root_inventory = _closed(
            root_inventory,
            {"schema", "dataset_id", "raw_objects", "derived_objects"},
            f"{field}.raw_root_manifest root",
        )
        if root_inventory["schema"] != "WAC_FORMATION_RAW_ROOT_V002":
            _refuse(f"{field}.raw_root_manifest schema is not registered")
        if root_inventory["dataset_id"] != dataset_id:
            _refuse(f"{field}.raw_root_manifest dataset_id mismatch")
        if root_inventory["raw_objects"] != sorted(raw_root_rows, key=lambda row: row["object_id"]):
            _refuse(f"{field}.raw_root_manifest raw inventory mismatch")
        if root_inventory["derived_objects"] != sorted(derived_root_rows, key=lambda row: row["object_id"]):
            _refuse(f"{field}.raw_root_manifest derived inventory mismatch")

        events = dataset["events"]
        if not isinstance(events, list) or not events:
            _refuse(f"{field}.events must be nonempty")
        event_ids: set[str] = set()
        event_routes: set[str] = set()
        event_values: set[str] = set()
        valid_event_count = 0
        repeated_read_present = False
        for index, raw_event in enumerate(events):
            event_field = f"{field}.events[{index}]"
            event = _closed(raw_event, EVENT_KEYS, event_field)
            event_id = _nonempty(event["event_id"], f"{event_field}.event_id")
            if event_id in event_ids:
                _refuse(f"{field}.events duplicates {event_id}")
            event_ids.add(event_id)
            disjoint("event_id", event_id, dataset_id)
            if event["physical_unit_id"] not in physical_unit_ids:
                _refuse(f"{event_field} references an unknown physical unit")
            if event["lot_id"] not in lot_ids:
                _refuse(f"{event_field} references an unknown lot")
            _nonempty(event["block_id"], f"{event_field}.block_id")
            event_value = _nonempty(event["event_value"], f"{event_field}.event_value")
            event_values.add(event_value)
            if event["route"] not in ROUTES:
                _refuse(f"{event_field}.route is not WRITE or SHAM")
            event_routes.add(event["route"])
            if _boolean(event["allocation_concealed"], f"{event_field}.allocation_concealed") is not True:
                _refuse(f"{event_field} allocation was not concealed")
            _digest(event["allocation_commitment_sha256"], f"{event_field}.allocation_commitment_sha256")
            generated_utc = _utc(event["event_generated_at_utc"], f"{event_field}.event_generated_at_utc")
            generated_ns = _integer(
                event["event_generated_monotonic_ns"],
                f"{event_field}.event_generated_monotonic_ns",
            )
            before = _closed(event["before"], BEFORE_KEYS, f"{event_field}.before")
            before_utc = _utc(before["sealed_at_utc"], f"{event_field}.before.sealed_at_utc")
            before_ns = _integer(before["sealed_monotonic_ns"], f"{event_field}.before.sealed_monotonic_ns")
            if _boolean(before["event_absent_when_sealed"], f"{event_field}.before.event_absent_when_sealed") is not True:
                _refuse(f"{event_field} BEFORE did not precede event existence")
            if not (before_utc < generated_utc and before_ns < generated_ns):
                _refuse(f"{event_field} BEFORE does not precede event generation")
            before_ids = _unique_strings(before["raw_object_ids"], f"{event_field}.before.raw_object_ids")
            formation = _closed(event["formation"], FORMATION_KEYS, f"{event_field}.formation")
            formation_start_utc = _utc(formation["started_at_utc"], f"{event_field}.formation.started_at_utc")
            formation_start_ns = _integer(formation["started_monotonic_ns"], f"{event_field}.formation.started_monotonic_ns")
            formation_end_utc = _utc(formation["ended_at_utc"], f"{event_field}.formation.ended_at_utc")
            formation_end_ns = _integer(formation["ended_monotonic_ns"], f"{event_field}.formation.ended_monotonic_ns")
            if not (generated_utc <= formation_start_utc <= formation_end_utc):
                _refuse(f"{event_field} formation UTC order is invalid")
            if not (generated_ns <= formation_start_ns <= formation_end_ns):
                _refuse(f"{event_field} formation monotonic order is invalid")
            if _boolean(formation["physical_target_coupling_measured"], f"{event_field}.formation.physical_target_coupling_measured") is not True:
                _refuse(f"{event_field} lacks measured target/control coupling")
            formation_ids = _unique_strings(formation["raw_object_ids"], f"{event_field}.formation.raw_object_ids")
            c_ext = _closed(event["c_ext_off"], BOUNDARY_KEYS, f"{event_field}.c_ext_off")
            c_utc = _utc(c_ext["observed_at_utc"], f"{event_field}.c_ext_off.observed_at_utc")
            c_ns = _integer(c_ext["observed_monotonic_ns"], f"{event_field}.c_ext_off.observed_monotonic_ns")
            if _boolean(c_ext["physically_measured"], f"{event_field}.c_ext_off.physically_measured") is not True:
                _refuse(f"{event_field} C_EXT_OFF is nominal rather than measured")
            c_ids = _unique_strings(c_ext["raw_object_ids"], f"{event_field}.c_ext_off.raw_object_ids")
            if not (formation_end_utc <= c_utc and formation_end_ns <= c_ns):
                _refuse(f"{event_field} C_EXT_OFF precedes formation end")
            g = _closed(event["g_surface_close"], G_KEYS, f"{event_field}.g_surface_close")
            if g["status"] not in {"CERTIFIED", "FAILED"}:
                _refuse(f"{event_field}.g_surface_close.status is not registered")
            if _boolean(g["separate_from_c_ext_off"], f"{event_field}.g_surface_close.separate_from_c_ext_off") is not True:
                _refuse(f"{event_field} conflates C_EXT_OFF and G_SURFACE_CLOSE")
            g_ids = _unique_strings(g["raw_object_ids"], f"{event_field}.g_surface_close.raw_object_ids")
            if set(c_ids) & set(g_ids):
                _refuse(f"{event_field} reuses one raw object for C_EXT_OFF and G_SURFACE_CLOSE")
            g_predicate_ids = _unique_strings(
                g["predicate_ids"], f"{event_field}.g_surface_close.predicate_ids"
            )
            g_predicate_refs.append((dataset_id, event_id, tuple(g_predicate_ids)))
            all_event_ids = _unique_strings(event["raw_object_ids"], f"{event_field}.raw_object_ids")
            stage_ids = set(before_ids + formation_ids + c_ids + g_ids)
            if not stage_ids <= set(all_event_ids):
                _refuse(f"{event_field} stage object is absent from event raw_object_ids")
            if any(object_id not in raw_by_id for object_id in all_event_ids):
                _refuse(f"{event_field} references an unknown or cross-dataset raw object")
            for object_id in all_event_ids:
                owner = raw_by_id[object_id]["event_id"]
                if owner != event_id:
                    _refuse(f"{event_field} raw object {object_id} has wrong event owner")
            for object_id in before_ids:
                if raw_by_id[object_id]["role"] != "BEFORE_STATE":
                    _refuse(f"{event_field} BEFORE object has wrong role")
            expected_formation_role = "FORMATION_INPUT" if event["route"] == "WRITE" else "ROUTING_OR_CONTROL"
            if not any(raw_by_id[object_id]["role"] == expected_formation_role for object_id in formation_ids):
                _refuse(f"{event_field} formation/control role is missing")
            if not all(raw_by_id[object_id]["role"] == "C_EXT_OFF_ASSAY" for object_id in c_ids):
                _refuse(f"{event_field} C_EXT_OFF object has wrong role")
            if not all(raw_by_id[object_id]["role"] == "G_SURFACE_CLOSE_ASSAY" for object_id in g_ids):
                _refuse(f"{event_field} G object has wrong role")

            reads = event["reads"]
            if not isinstance(reads, list):
                _refuse(f"{event_field}.reads must be a list")
            read_ids: set[str] = set()
            read_times: list[int] = []
            read_utcs: list[datetime] = []
            mission_reads = 0
            for r_index, raw_read in enumerate(reads):
                read_field = f"{event_field}.reads[{r_index}]"
                read = _closed(raw_read, READ_KEYS, read_field)
                read_id = _nonempty(read["read_id"], f"{read_field}.read_id")
                if read_id in read_ids:
                    _refuse(f"{event_field}.reads duplicates {read_id}")
                read_ids.add(read_id)
                read_utc = _utc(read["observed_at_utc"], f"{read_field}.observed_at_utc")
                read_ns = _integer(read["observed_monotonic_ns"], f"{read_field}.observed_monotonic_ns")
                if read_times and (read_ns <= read_times[-1] or read_utc <= read_utcs[-1]):
                    _refuse(f"{event_field}.reads are not strictly ordered")
                read_times.append(read_ns)
                read_utcs.append(read_utc)
                _nonempty(read["query_id"], f"{read_field}.query_id")
                mission_flag = _boolean(read["mission_read"], f"{read_field}.mission_read")
                mission_reads += int(mission_flag)
                read_object_ids = _unique_strings(read["raw_object_ids"], f"{read_field}.raw_object_ids")
                if not set(read_object_ids) <= set(all_event_ids):
                    _refuse(f"{read_field} object is absent from event raw_object_ids")
                if not all(raw_by_id[object_id]["role"] in {"READ_OBSERVATION", "READ_BACKACTION_ASSAY"} for object_id in read_object_ids):
                    _refuse(f"{read_field} contains a non-read object")
            if len(reads) >= 2:
                repeated_read_present = True

            validity = event["validity_status"]
            if validity not in VALIDITY:
                _refuse(f"{event_field}.validity_status is not registered")
            if g["status"] == "FAILED" and validity != "UNSCOREABLE":
                _refuse(f"{event_field} failed G must make the event UNSCOREABLE")
            if g["status"] == "CERTIFIED" and (
                g["observed_at_utc"] is None or g["observed_monotonic_ns"] is None
            ):
                _refuse(f"{event_field} certified G lacks time")
            hold = event["common_hold"]
            if g["status"] == "CERTIFIED" and validity == "VALID":
                g_utc = _utc(g["observed_at_utc"], f"{event_field}.g_surface_close.observed_at_utc")
                g_ns = _integer(g["observed_monotonic_ns"], f"{event_field}.g_surface_close.observed_monotonic_ns")
                if not (c_utc <= g_utc and c_ns <= g_ns):
                    _refuse(f"{event_field} G precedes C_EXT_OFF")
                hold = _closed(hold, HOLD_KEYS, f"{event_field}.common_hold")
                hold_start_utc = _utc(hold["started_at_utc"], f"{event_field}.common_hold.started_at_utc")
                hold_start_ns = _integer(hold["started_monotonic_ns"], f"{event_field}.common_hold.started_monotonic_ns")
                hold_end_utc = _utc(hold["ended_at_utc"], f"{event_field}.common_hold.ended_at_utc")
                hold_end_ns = _integer(hold["ended_monotonic_ns"], f"{event_field}.common_hold.ended_monotonic_ns")
                if hold_start_utc != g_utc or hold_start_ns != g_ns:
                    _refuse(f"{event_field} hold does not start at G")
                if hold_end_utc < hold_start_utc or hold_end_ns < hold_start_ns:
                    _refuse(f"{event_field} hold interval is reversed")
                if _boolean(hold["writer_off_throughout"], f"{event_field}.common_hold.writer_off_throughout") is not True:
                    _refuse(f"{event_field} common hold is not writer-off")
                hold_ids = _unique_strings(hold["raw_object_ids"], f"{event_field}.common_hold.raw_object_ids")
                if not set(hold_ids) <= set(all_event_ids):
                    _refuse(f"{event_field} hold object is absent from event raw_object_ids")
                if not all(raw_by_id[object_id]["role"] == "COMMON_HOLD_ENVIRONMENT" for object_id in hold_ids):
                    _refuse(f"{event_field} hold object has wrong role")
                if not reads or mission_reads != 1:
                    _refuse(f"{event_field} valid event requires exactly one mission read")
                if any(time <= hold_end_ns for time in read_times) or any(
                    read_utc <= hold_end_utc for read_utc in read_utcs
                ):
                    _refuse(f"{event_field} READ does not follow the common hold")
                mission_read = next(read for read in reads if read["mission_read"])
                delta_s = (mission_read["observed_monotonic_ns"] - g_ns) / 1e9
                mission_utc = _utc(
                    mission_read["observed_at_utc"],
                    f"{event_field}.mission_read.observed_at_utc",
                )
                utc_delta_s = (mission_utc - g_utc).total_seconds()
                if (
                    abs(delta_s - mission_value) > mission_uncertainty
                    or abs(utc_delta_s - mission_value) > mission_uncertainty
                ):
                    _refuse(f"{event_field} mission read is outside frozen mission-time uncertainty")
                valid_event_count += 1
            else:
                if g["status"] == "FAILED":
                    if g["observed_at_utc"] is not None or g["observed_monotonic_ns"] is not None:
                        _refuse(f"{event_field} failed G must not declare a closure time")
                if hold is not None or reads:
                    _refuse(f"{event_field} unclosed/unscoreable event cannot declare hold or READ")

        for raw_id, item in raw_by_id.items():
            if item["event_id"] is not None and item["event_id"] not in event_ids:
                _refuse(f"{field}.raw_objects object {raw_id} references an unknown event")
        if role in VALIDATION_ROLES:
            if event_routes != ROUTES:
                _refuse(f"{field} validation dataset requires both WRITE and SHAM")
            if len(event_values) < 2:
                _refuse(f"{field} validation dataset requires at least two event values")
            if valid_event_count == 0:
                _refuse(f"{field} validation dataset has no valid closed event")
            if not repeated_read_present:
                _refuse(f"{field} validation dataset lacks repeated READ")
            required_roles = {
                "BEFORE_STATE", "FORMATION_INPUT", "ROUTING_OR_CONTROL",
                "C_EXT_OFF_ASSAY", "G_SURFACE_CLOSE_ASSAY",
                "COMMON_HOLD_ENVIRONMENT", "READ_OBSERVATION",
                "READ_BACKACTION_ASSAY", "UNCERTAINTY_CALIBRATION",
            }
            if not required_roles <= roles_present:
                _refuse(f"{field} validation raw roles are incomplete")

        expected_objects = dataset["expected_objects"]
        if not isinstance(expected_objects, list) or not expected_objects:
            _refuse(f"{field}.expected_objects must be nonempty")
        expected_by_id: dict[str, dict[str, Any]] = {}
        for index, raw_expected in enumerate(expected_objects):
            expected_field = f"{field}.expected_objects[{index}]"
            expected = _closed(raw_expected, EXPECTED_OBJECT_KEYS, expected_field)
            expected_id = _nonempty(
                expected["expected_object_id"], f"{expected_field}.expected_object_id"
            )
            if expected_id in expected_by_id:
                _refuse(f"{field}.expected_objects duplicates {expected_id}")
            expected_by_id[expected_id] = expected
            disjoint("expected_object_id", expected_id, dataset_id)
            if expected["event_id"] is not None and expected["event_id"] not in event_ids:
                _refuse(f"{expected_field} references an unknown event")
            if (
                expected["physical_unit_id"] is not None
                and expected["physical_unit_id"] not in physical_unit_ids
            ):
                _refuse(f"{expected_field} references an unknown physical unit")
            _nonempty(expected["role"], f"{expected_field}.role")
            if expected_id in raw_by_id:
                present = raw_by_id[expected_id]
            elif expected_id in derived_by_id:
                present = derived_by_id[expected_id]
            else:
                present = None
            if present is not None and any(
                present[key] != expected[key]
                for key in ("event_id", "physical_unit_id", "role")
            ):
                _refuse(f"{expected_field} metadata differs from the present object")

        missing = dataset["missing_objects"]
        if not isinstance(missing, list):
            _refuse(f"{field}.missing_objects must be a list")
        critical_missing = False
        missing_ids: set[str] = set()
        for index, raw_missing in enumerate(missing):
            missing_field = f"{field}.missing_objects[{index}]"
            item = _closed(raw_missing, MISSING_KEYS, missing_field)
            expected_id = _nonempty(
                item["expected_object_id"], f"{missing_field}.expected_object_id"
            )
            if expected_id in missing_ids:
                _refuse(f"{field}.missing_objects duplicates {expected_id}")
            missing_ids.add(expected_id)
            if expected_id not in expected_by_id:
                _refuse(f"{missing_field} is absent from expected_objects")
            if expected_id in object_ids:
                _refuse(f"{missing_field} marks a present object missing")
            if item["event_id"] is not None and item["event_id"] not in event_ids:
                _refuse(f"{missing_field} references an unknown event")
            if item["physical_unit_id"] is not None and item["physical_unit_id"] not in physical_unit_ids:
                _refuse(f"{missing_field} references an unknown physical unit")
            _nonempty(item["role"], f"{missing_field}.role")
            expected = expected_by_id[expected_id]
            if any(
                item[key] != expected[key]
                for key in ("event_id", "physical_unit_id", "role")
            ):
                _refuse(f"{missing_field} metadata differs from expected_objects")
            _nonempty(item["reason"], f"{missing_field}.reason")
            critical_missing |= _boolean(item["critical"], f"{missing_field}.critical")
            if item["disposition"] not in {"RETAIN", "UNSCOREABLE_EVENT", "UNSCOREABLE_DATASET"}:
                _refuse(f"{missing_field}.disposition is not registered")
        if set(expected_by_id) != object_ids | missing_ids:
            _refuse(f"{field} expected/present/missing object inventory is not exact")

        deviations = dataset["deviations"]
        if not isinstance(deviations, list):
            _refuse(f"{field}.deviations must be a list")
        critical_deviation = False
        deviation_ids: set[str] = set()
        for index, raw_deviation in enumerate(deviations):
            deviation_field = f"{field}.deviations[{index}]"
            item = _closed(raw_deviation, DEVIATION_KEYS, deviation_field)
            deviation_id = _nonempty(item["deviation_id"], f"{deviation_field}.deviation_id")
            if deviation_id in deviation_ids:
                _refuse(f"{field}.deviations duplicates {deviation_id}")
            deviation_ids.add(deviation_id)
            if item["event_id"] is not None and item["event_id"] not in event_ids:
                _refuse(f"{deviation_field} references an unknown event")
            _nonempty(item["description"], f"{deviation_field}.description")
            critical_deviation |= _boolean(item["critical"], f"{deviation_field}.critical")
            _boolean(item["allocation_blind"], f"{deviation_field}.allocation_blind")
            if item["disposition"] not in {"RETAIN", "UNSCOREABLE_EVENT", "UNSCOREABLE_DATASET"}:
                _refuse(f"{deviation_field}.disposition is not registered")

        signatures = dataset["custody_signatures"]
        if not isinstance(signatures, list) or not signatures:
            _refuse(f"{field}.custody_signatures must be nonempty")
        signature_ids: set[str] = set()
        signature_roles: set[str] = set()
        for index, raw_signature in enumerate(signatures):
            signature_field = f"{field}.custody_signatures[{index}]"
            item = _closed(raw_signature, SIGNATURE_KEYS, signature_field)
            signature_id = _nonempty(item["signature_id"], f"{signature_field}.signature_id")
            if signature_id in signature_ids:
                _refuse(f"{field}.custody_signatures duplicates {signature_id}")
            signature_ids.add(signature_id)
            signature_roles.add(_nonempty(item["role"], f"{signature_field}.role"))
            _nonempty(item["signer_id"], f"{signature_field}.signer_id")
            _digest(item["signed_sha256"], f"{signature_field}.signed_sha256")
            _utc(item["signed_at_utc"], f"{signature_field}.signed_at_utc")
            _nonempty(item["algorithm"], f"{signature_field}.algorithm")
            _nonempty(item["signature"], f"{signature_field}.signature")
        if role in VALIDATION_ROLES and not {"DATA_CUSTODIAN", "BLIND_LOCK"} <= signature_roles:
            _refuse(f"{field} validation custody signatures are incomplete")

        dataset_objects[dataset_id] = object_ids
        dataset_scoreable[dataset_id] = bool(
            valid_event_count > 0
            and all_raw_complete
            and not critical_missing
            and not critical_deviation
        )

    all_dataset_ids = set(dataset_ids)
    for dataset_id, dataset in dataset_by_id.items():
        declared = set(_unique_strings(
            dataset["independent_from_dataset_ids"],
            f"dataset {dataset_id}.independent_from_dataset_ids",
            allow_empty=True,
        ))
        if declared != all_dataset_ids - {dataset_id}:
            _refuse(f"dataset {dataset_id} independence declaration is incomplete")

    frozen_by_id = {item["object_id"]: item for item in frozen_items.values()}
    dataset_object_records: dict[str, dict[str, dict[str, Any]]] = {}
    for dataset in datasets:
        dataset_object_records[dataset["dataset_id"]] = {
            item["object_id"]: item
            for item in dataset["raw_objects"] + dataset["derived_objects"]
        }
    mapping_coverage: dict[str, set[str]] = {
        dataset_id: set() for dataset_id in dataset_ids
    }
    seen_mapping_sources: set[tuple[str, str | None, str | None, str, str]] = set()
    for contract_field, source in mapping_sources:
        dataset_id = source["dataset_id"]
        event_id = source["event_id"]
        object_id = source["object_id"]
        selector = source["selector"]
        source_key = (contract_field, dataset_id, event_id, object_id, selector)
        if source_key in seen_mapping_sources:
            _refuse("declarative mapping duplicates a normalized source binding")
        seen_mapping_sources.add(source_key)
        if dataset_id is None:
            if event_id is not None or object_id not in frozen_by_id:
                _refuse("declarative mapping frozen source has invalid identity")
            if contract_field not in {"K_theta", "M_phi", "F_policy"}:
                _refuse(f"declarative mapping field {contract_field} requires dataset evidence")
            object_record = frozen_by_id[object_id]
        else:
            if dataset_id not in dataset_object_records:
                _refuse("declarative mapping references an unknown dataset")
            object_record = dataset_object_records[dataset_id].get(object_id)
            if object_record is None:
                _refuse("declarative mapping references an unknown or cross-dataset object")
            if object_record["event_id"] != event_id:
                _refuse("declarative mapping event/object join is not exact")
            mapping_coverage[dataset_id].add(contract_field)
        if source["path"] != object_record["path"]:
            _refuse("declarative mapping path differs from the content-addressed object")
        selector_kind = source["selector_kind"]
        media_type = object_record["media_type"].lower()
        if selector_kind == "WHOLE_OBJECT" and selector != "/":
            _refuse("WHOLE_OBJECT mapping selector must be exactly /")
        if selector_kind == "JSON_POINTER" and (
            not selector.startswith("/") or "json" not in media_type
        ):
            _refuse("JSON_POINTER mapping is not resolvable against a JSON object")
        if selector_kind == "CSV_COLUMN" and (
            not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.:-]*", selector)
            or "csv" not in media_type
        ):
            _refuse("CSV_COLUMN mapping is not a closed column selector")
        if selector_kind == "HDF5_DATASET" and (
            not selector.startswith("/")
            or not any(token in media_type for token in ("hdf", "h5"))
        ):
            _refuse("HDF5_DATASET mapping is not resolvable against an HDF5 object")
        if source["value_type"] in {"INTEGER", "FLOAT", "TIMESTAMP"}:
            if dataset_id is None:
                _refuse("numeric frozen mappings require an explicit dataset calibration")
            if object_id in {
                raw["object_id"] for raw in dataset_by_id[dataset_id]["raw_objects"]
            }:
                calibration_ids = {
                    channel["calibration_id"] for channel in object_record["channels"]
                }
                clock_ids = {object_record["clock_id"]}
            else:
                raw_by_object_id = {
                    raw["object_id"]: raw for raw in dataset_by_id[dataset_id]["raw_objects"]
                }
                parents = [raw_by_object_id[parent] for parent in object_record["parent_object_ids"]]
                calibration_ids = {
                    channel["calibration_id"]
                    for parent in parents
                    for channel in parent["channels"]
                }
                clock_ids = {parent["clock_id"] for parent in parents}
            if source["calibration_id"] not in calibration_ids:
                _refuse("declarative mapping calibration reference is unresolved")
            if source["clock_id"] not in clock_ids:
                _refuse("declarative mapping clock reference is unresolved")

    for dataset_id, dataset in dataset_by_id.items():
        if dataset["dataset_role"] in VALIDATION_ROLES:
            missing_fields = DATASET_MAPPING_FIELDS - mapping_coverage[dataset_id]
            if missing_fields:
                _refuse(
                    f"declarative mapping omits validation dataset {dataset_id} fields "
                    f"{sorted(missing_fields)}"
                )

    predicates = manifest["generic_predicates"]
    if not isinstance(predicates, list):
        _refuse("generic_predicates must be a list")
    predicate_ids: set[str] = set()
    predicates_by_id: dict[str, dict[str, Any]] = {}
    predicate_roles_by_dataset: dict[str, list[str]] = {dataset_id: [] for dataset_id in dataset_ids}
    for index, raw_predicate in enumerate(predicates):
        field = f"generic_predicates[{index}]"
        item = _closed(raw_predicate, PREDICATE_KEYS, field)
        predicate_id = _nonempty(item["predicate_id"], f"{field}.predicate_id")
        if predicate_id in predicate_ids:
            _refuse(f"duplicate predicate_id {predicate_id}")
        predicate_ids.add(predicate_id)
        predicates_by_id[predicate_id] = item
        dataset_id = item["dataset_id"]
        if dataset_id not in dataset_ids:
            _refuse(f"{field} references an unknown dataset")
        predicate_role = item["predicate_role"]
        if predicate_role not in REQUIRED_PREDICATE_ROLES:
            _refuse(f"{field}.predicate_role is not registered")
        predicate_roles_by_dataset[dataset_id].append(predicate_role)
        if item["stage"] not in STAGES:
            _refuse(f"{field}.stage is not registered")
        _nonempty(item["metric_id"], f"{field}.metric_id")
        source_ids = _unique_strings(item["source_object_ids"], f"{field}.source_object_ids")
        if any(object_id not in dataset_objects[dataset_id] for object_id in source_ids):
            _refuse(f"{field} source object crosses datasets or is unknown")
        matching_names = [
            name for name, frozen_item in frozen_items.items()
            if frozen_item["object_id"] == item["rule_object_id"]
        ]
        if len(matching_names) != 1:
            _refuse(f"{field}.rule_object_id does not identify one frozen object")
        rule_name = matching_names[0]
        if _digest(item["rule_sha256"], f"{field}.rule_sha256") != frozen_items[rule_name]["sha256"]:
            _refuse(f"{field}.rule_sha256 mismatch")
        frozen_before_access = _boolean(item["frozen_before_access"], f"{field}.frozen_before_access")
        if dataset_by_id[dataset_id]["dataset_role"] in VALIDATION_ROLES and not frozen_before_access:
            _refuse(f"{field} validation predicate was not frozen before access")
        acceptance = _closed(item["acceptance"], ACCEPTANCE_KEYS, f"{field}.acceptance")
        operator = acceptance["operator"]
        if operator not in PREDICATE_OPERATORS:
            _refuse(f"{field}.acceptance.operator is not registered")
        unit = _nonempty(acceptance["unit"], f"{field}.acceptance.unit")
        if unit not in SI_UNITS:
            _refuse(f"{field}.acceptance.unit is not in the closed SI unit registry")
        lower = acceptance["lower"]
        upper = acceptance["upper"]
        if operator == "GE_LOWER_BOUND":
            if lower is None or upper is not None:
                _refuse(f"{field} GE_LOWER_BOUND requires lower and null upper")
            _finite(lower, f"{field}.acceptance.lower")
        elif operator == "LE_UPPER_BOUND":
            if lower is not None or upper is None:
                _refuse(f"{field} LE_UPPER_BOUND requires null lower and upper")
            _finite(upper, f"{field}.acceptance.upper")
        else:
            if lower is None or upper is None:
                _refuse(f"{field} WITHIN_CLOSED_INTERVAL requires lower and upper")
            low = _finite(lower, f"{field}.acceptance.lower")
            high = _finite(upper, f"{field}.acceptance.upper")
            if low > high:
                _refuse(f"{field} acceptance interval is reversed")

    for dataset_id, dataset in dataset_by_id.items():
        if dataset["dataset_role"] in VALIDATION_ROLES:
            roles = predicate_roles_by_dataset[dataset_id]
            if len(roles) != len(set(roles)) or set(roles) != REQUIRED_PREDICATE_ROLES:
                _refuse(f"dataset {dataset_id} must have exactly one fully specified predicate per required role")

    for dataset_id, event_id, g_ids in g_predicate_refs:
        for predicate_id in g_ids:
            predicate = predicates_by_id.get(predicate_id)
            if predicate is None:
                _refuse(f"event {event_id} references an unknown G predicate")
            if (
                predicate["dataset_id"] != dataset_id
                or predicate["predicate_role"] != "G_SURFACE_CLOSE"
            ):
                _refuse(f"event {event_id} G predicate crosses dataset or role")

    product = _closed(manifest["product_reproduction"], PRODUCT_KEYS, "product_reproduction")
    if product["status"] != "NOT_DEMONSTRATED":
        _refuse("formation input cannot assert the separate product-reproduction gate")
    physicists = _unique_strings(
        product["external_physicist_ids"],
        "product_reproduction.external_physicist_ids",
        allow_empty=True,
    )
    release_hash = product["public_urm_release_sha256"]
    for key in (
        "used_only_public_instructions",
        "no_private_help",
        "different_real_world_datasets",
        "any_and_every_surface_coverage",
    ):
        _boolean(product[key], f"product_reproduction.{key}")
    evidence_ids = _unique_strings(
        product["evidence_object_ids"],
        "product_reproduction.evidence_object_ids",
        allow_empty=True,
    )
    if release_hash is not None or physicists or evidence_ids:
        _refuse("undemonstrated product reproduction cannot name release, physicists, or evidence")
    if any(product[key] for key in (
        "used_only_public_instructions",
        "no_private_help",
        "different_real_world_datasets",
        "any_and_every_surface_coverage",
    )):
        _refuse("undemonstrated product reproduction fields must remain false")

    signatures = manifest["signatures"]
    if not isinstance(signatures, list) or not signatures:
        _refuse("signatures must be nonempty")
    signature_ids: set[str] = set()
    for index, raw_signature in enumerate(signatures):
        field = f"signatures[{index}]"
        item = _closed(raw_signature, SIGNATURE_KEYS, field)
        signature_id = _nonempty(item["signature_id"], f"{field}.signature_id")
        if signature_id in signature_ids:
            _refuse(f"duplicate root signature_id {signature_id}")
        signature_ids.add(signature_id)
        _nonempty(item["role"], f"{field}.role")
        _nonempty(item["signer_id"], f"{field}.signer_id")
        _digest(item["signed_sha256"], f"{field}.signed_sha256")
        _utc(item["signed_at_utc"], f"{field}.signed_at_utc")
        _nonempty(item["algorithm"], f"{field}.algorithm")
        _nonempty(item["signature"], f"{field}.signature")

    return FormationInput(
        manifest_path=resolved_manifest,
        manifest_sha256=sha256_bytes(manifest_bytes),
        manifest=manifest,
        dataset_input_complete=tuple(sorted(dataset_scoreable.items())),
        dataset_object_ids=tuple(
            sorted(
                (dataset_id, tuple(sorted(object_ids)))
                for dataset_id, object_ids in dataset_objects.items()
            )
        ),
        referenced_hashes=tuple(sorted(referenced_hashes)),
    )


def attach_formation_execution(
    formation_input: FormationInput,
    execution_path: str | Path,
) -> FormationExecution:
    """Attach numerical outputs to a previously validated, immutable V002 input."""
    execution_path = Path(execution_path)
    if execution_path.is_symlink() or execution_path.parent.is_symlink():
        _refuse("submitted execution or its bundle directory must not be a symlink")
    root = formation_input.manifest_path.parent.resolve(strict=True)
    try:
        resolved_execution = execution_path.resolve(strict=True)
        resolved_execution.relative_to(root)
    except (OSError, ValueError):
        _refuse("submitted execution escapes the formation bundle")
    execution_bytes, execution = _read_json(resolved_execution, "formation execution")
    execution = _closed(execution, EXECUTION_KEYS, "formation execution")
    if execution["schema"] != EXECUTION_SCHEMA:
        _refuse("formation execution schema is not registered")
    if _digest(
        execution["input_manifest_sha256"],
        "formation execution.input_manifest_sha256",
    ) != formation_input.manifest_sha256:
        _refuse("formation execution does not bind the exact input manifest")

    frozen = formation_input.manifest["frozen_objects"]
    expected_bindings = {
        "release_object_id": frozen["urm_release"]["object_id"],
        "release_sha256": frozen["urm_release"]["sha256"],
        "source_sha256": frozen["urm_source"]["sha256"],
        "environment_sha256": frozen["urm_environment"]["sha256"],
        "analysis_object_id": frozen["analysis"]["object_id"],
        "analysis_sha256": frozen["analysis"]["sha256"],
    }
    for key, expected in expected_bindings.items():
        value = execution[key]
        if key.endswith("sha256"):
            _digest(value, f"formation execution.{key}")
        else:
            _nonempty(value, f"formation execution.{key}")
        if value != expected:
            _refuse(f"formation execution.{key} does not match the frozen input")

    used_paths: set[Path] = {formation_input.manifest_path, resolved_execution}
    used_identities: set[tuple[int, int]] = {
        _file_identity(formation_input.manifest_path),
        _file_identity(resolved_execution),
    }
    for relative, _ in formation_input.referenced_hashes:
        path = (root / relative).resolve(strict=True)
        used_paths.add(path)
        used_identities.add(_file_identity(path))
    execution_hashes: list[tuple[str, str]] = []
    input_object_ids = {
        item["object_id"] for item in formation_input.manifest["frozen_objects"].values()
    }
    for dataset in formation_input.manifest["datasets"]:
        input_object_ids.update(item["object_id"] for item in dataset["raw_objects"])
        input_object_ids.update(item["object_id"] for item in dataset["derived_objects"])
        input_object_ids.add(dataset["raw_root_manifest"]["object_id"])
        input_object_ids.add(dataset["provenance"]["acquisition_record"]["object_id"])
        input_object_ids.add(dataset["freeze"]["access_log"]["object_id"])
    execution_object_ids: set[str] = set()
    for name in (
        "invocation",
        "computed_j",
        "evaluated_falsifier",
        "result",
        "stdout",
        "stderr",
    ):
        field = f"formation execution.{name}"
        ref = _closed(execution[name], FILE_REF_KEYS, field)
        object_id = _nonempty(ref["object_id"], f"{field}.object_id")
        if object_id in execution_object_ids or object_id in input_object_ids:
            _refuse("formation execution file object IDs must be unique and new")
        execution_object_ids.add(object_id)
        _nonempty(ref["media_type"], f"{field}.media_type")
        expected_hash = _digest(ref["sha256"], f"{field}.sha256")
        path = _safe_relative(root, ref["path"], f"{field}.path")
        if path.stat().st_size == 0:
            _refuse(f"{field}.path is a zero-byte placeholder")
        if path in used_paths:
            _refuse(f"{field} duplicates a formation-input or execution path")
        identity = _file_identity(path)
        if identity in used_identities:
            _refuse(f"{field} aliases a formation-input or execution inode")
        if sha256_file(path) != expected_hash:
            _refuse(f"{field} hash mismatch")
        used_paths.add(path)
        used_identities.add(identity)
        execution_hashes.append((str(path.relative_to(root)), expected_hash))

    dataset_objects = dict(formation_input.dataset_object_ids)
    object_hashes: dict[str, str] = {}
    object_dataset: dict[str, str] = {}
    for dataset in formation_input.manifest["datasets"]:
        dataset_id = dataset["dataset_id"]
        for obj in dataset["raw_objects"] + dataset["derived_objects"]:
            object_hashes[obj["object_id"]] = obj["sha256"]
            object_dataset[obj["object_id"]] = dataset_id

    predicate_specs = {
        item["predicate_id"]: item
        for item in formation_input.manifest["generic_predicates"]
    }
    measurements = execution["measurements"]
    if not isinstance(measurements, list):
        _refuse("formation execution.measurements must be a list")
    seen_predicates: set[str] = set()
    scores: list[PredicateScore] = []
    for index, raw_measurement in enumerate(measurements):
        field = f"formation execution.measurements[{index}]"
        measurement = _closed(raw_measurement, MEASUREMENT_KEYS, field)
        predicate_id = _nonempty(measurement["predicate_id"], f"{field}.predicate_id")
        if predicate_id in seen_predicates:
            _refuse(f"formation execution duplicates predicate {predicate_id}")
        seen_predicates.add(predicate_id)
        spec = predicate_specs.get(predicate_id)
        if spec is None:
            _refuse(f"{field} references an unknown frozen predicate")
        dataset_id = measurement["dataset_id"]
        if dataset_id != spec["dataset_id"]:
            _refuse(f"{field}.dataset_id differs from the frozen predicate")
        source_ids = _unique_strings(measurement["source_object_ids"], f"{field}.source_object_ids")
        if source_ids != spec["source_object_ids"]:
            _refuse(f"{field}.source_object_ids differ from the frozen predicate")
        measurement_object_id = _nonempty(
            measurement["measurement_object_id"], f"{field}.measurement_object_id"
        )
        if measurement_object_id not in dataset_objects[dataset_id]:
            _refuse(f"{field}.measurement_object_id crosses datasets or is unknown")
        if measurement_object_id not in source_ids:
            _refuse(f"{field}.measurement_object_id is not a frozen predicate source")
        if object_dataset[measurement_object_id] != dataset_id:
            _refuse(f"{field}.measurement_object_id has the wrong dataset owner")
        measurement_hash = _digest(
            measurement["measurement_sha256"], f"{field}.measurement_sha256"
        )
        if measurement_hash != object_hashes[measurement_object_id]:
            _refuse(f"{field}.measurement_sha256 does not bind the named object")
        passed = _score_predicate(spec, measurement["observed"], field)
        scores.append(PredicateScore(
            predicate_id=predicate_id,
            dataset_id=dataset_id,
            predicate_role=spec["predicate_role"],
            passed=passed,
        ))
    if seen_predicates != set(predicate_specs):
        missing = sorted(set(predicate_specs) - seen_predicates)
        _refuse(f"formation execution omits frozen predicates {missing}")

    return FormationExecution(
        formation_input=formation_input,
        execution_path=resolved_execution,
        execution_sha256=sha256_bytes(execution_bytes),
        predicate_scores=tuple(scores),
        referenced_hashes=tuple(sorted(execution_hashes)),
    )


def assess_validation_pair(
    formation_input: FormationInput | str | Path,
) -> dict[str, Any]:
    """Return structural D<->V eligibility without promoting a science result."""
    if not isinstance(formation_input, FormationInput):
        formation_input = load_formation_input(formation_input)
    certificate = formation_input.certificate()
    return {
        "schema": "WAC_FORMATION_VALIDATION_PAIR_ASSESSMENT_V002",
        "manifest_sha256": formation_input.manifest_sha256,
        "dataset_roles": certificate["dataset_roles"],
        "dataset_independence": certificate["dataset_independence"],
        "validation_dataset_eligible": certificate["validation_dataset_eligible"],
        "scientific_verdict": "NONE_NOT_SCORED",
        "scientific_validation_authorized": False,
        "record_formation_claim_authorized": False,
        "universal_claim_authorized": False,
        "gravity_claim_authorized": False,
        "program_completion_authorized": False,
        "product_reproduction_attested": False,
    }


def certificate_json(observation: FormationInput | FormationExecution) -> str:
    return json.dumps(observation.certificate(), indent=2, sort_keys=True) + "\n"
