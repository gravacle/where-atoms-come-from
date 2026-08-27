"""Missing-data-aware, proof-first execution frontier for the URM.

This layer coordinates work across incomplete public datasets.  It never converts
missing data into a failed physical prediction, never converts a paper summary into
event-level evidence, and never emits an authoritative proof or scientific-readiness
decision.  Its strongest completion output is a nonauthoritative input assertion that
still requires separate authoritative admission and scoring.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any
import unicodedata


SCHEMA = "WAC_PROOF_FRONTIER_V001"
CERTIFICATE_SCHEMA = "WAC_PROOF_FRONTIER_CERTIFICATE_V001"

ROOT_KEYS = {
    "schema",
    "frontier_id",
    "profile_id",
    "policy",
    "external_actions",
    "proofs",
    "datasets",
    "objects",
    "obligations",
    "theory_tracks",
}
POLICY_KEYS = {
    "goal",
    "public_data_only",
    "validation_preserved",
    "theory_fallback_rule",
    "policy_object_ids",
}
EXTERNAL_ACTION_KEYS = {
    "authorized",
    "executed",
    "third_party_contact",
    "private_data",
    "new_acquisition",
}
PROOF_KEYS = {
    "proof_id",
    "target",
    "required_obligation_ids",
    "authoritative_output",
}
DATASET_KEYS = {
    "dataset_id",
    "role",
    "availability",
    "evidence_class",
    "public",
    "opened",
    "source_id",
    "acquisition_id",
    "physical_unit_namespace",
    "acquisition_ancestry",
    "prior_use",
    "independent_from_dataset_ids",
    "frozen_before_response_access",
    "role_locked_before_response_access",
    "upstream_status",
    "qualification_object_ids",
    "measured_channel_ids",
    "missing_channel_ids",
    "object_ids",
}
OBJECT_KEYS = {
    "object_id",
    "path",
    "sha256",
    "media_type",
    "evidence_class",
    "dataset_id",
    "stage",
    "role",
}
OBLIGATION_KEYS = {
    "obligation_id",
    "proof_id",
    "label",
    "required",
    "dependencies",
    "status",
    "dataset_ids",
    "evidence_object_ids",
    "derivation",
    "missing_channel_ids",
    "theory_track_ids",
}
DERIVATION_KEYS = {
    "mode",
    "transform_object_id",
    "frozen_before_response_access",
    "identifiable_on_declared_support",
    "uncertainty_propagated",
    "uses_tested_response_as_control",
}
THEORY_KEYS = {
    "theory_id",
    "label",
    "mathematics_closed",
    "adversarial_checks_pass",
    "predictions",
    "independent_event_level_domain_ids",
    "held_out_validation_result_object_ids",
    "remaining_blockers",
}
PREDICTION_KEYS = {
    "prediction_id",
    "status",
    "dataset_ids",
    "evidence_object_ids",
}
BLOCKER_KEYS = {
    "blocker_id",
    "blocker_class",
    "description",
    "affected_obligation_ids",
    "missing_channel_ids",
}

PROOF_IDS = {"URF", "GE"}
AUTHORITATIVE_OUTPUT = "NO_PROOF_OUTPUT"
GOAL = "PROOF"
THEORY_FALLBACK = "SUBSTANTIALLY_PROVEN_ONLY_UNAVAILABLE_DATA"
PROFILES = {"PUBLIC_TWO_PROOF_PROFILE_V001", "SYNTHETIC_TEST_PROFILE_V001"}
CANONICAL_REQUIRED_OBLIGATIONS = {
    "URF": {
        "URF_EVENT_SPINE",
        "URF_BEFORE",
        "URF_CAUSAL_CONTRAST",
        "URF_COUPLING_OFF",
        "URF_G",
        "URF_HOLD",
        "URF_READ",
        "URF_D_CAUSAL",
        "URF_DEVELOPMENT",
        "URF_VALIDATION",
        "URF_UNIVERSALITY",
        "URF_REPRODUCTION",
        "URF_TARGET",
    },
    "GE": {
        "GE_JOIN",
        "GE_D",
        "GE_ANCESTRY",
        "GE_REMOTE_PROBE",
        "GE_ETA_PRED",
        "GE_COVARIANT_RESPONSE",
        "GE_CONTROLS",
        "GE_VALIDATION",
        "GE_TARGET",
    },
}

DATASET_ROLES = {
    "CALIBRATION",
    "DEVELOPMENT",
    "VALIDATION",
    "REPRODUCTION",
    "THEORY_ONLY",
}
AVAILABILITIES = {
    "AVAILABLE_PUBLIC",
    "RESERVED_UNOPENED_PUBLIC",
    "UNAVAILABLE_PUBLIC_DATA",
}
EVIDENCE_CLASSES = {
    "REAL_WORLD_EVENT_LEVEL",
    "PAPER_LEVEL_SUMMARY",
    "SYNTHETIC_TEST_ONLY",
    "THEORY_ONLY",
}
UPSTREAM_STATUSES = {
    "QUALIFIED_EVENT_LEVEL",
    "PARTIAL_OR_UNSCOREABLE",
    "PAPER_SUMMARY_ONLY",
    "SYNTHETIC_ONLY",
    "THEORY_ONLY",
}
PRIOR_USES = {"NONE", "CALIBRATION", "MODEL_SELECTION", "DEVELOPMENT"}
OBJECT_STAGES = {
    "POLICY",
    "METADATA",
    "RAW_EVENT",
    "CALIBRATION",
    "DERIVATION",
    "PAPER_SUMMARY",
    "THEORY",
    "SYNTHETIC",
    "VALIDATION_RESPONSE",
    "RESULT",
}
OBJECT_ROLES = {
    "PUBLIC_DATA_POLICY",
    "POLICY_AMENDMENT",
    "EVENT_LEVEL_QUALIFICATION_CERTIFICATE",
    "RAW_EVENT_DATA",
    "CALIBRATION_DATA",
    "DERIVED_EVENT_DATA",
    "PAPER_SUMMARY_DATA",
    "SYNTHETIC_TEST_FIXTURE",
    "THEORY_DERIVATION",
    "THEORY_TEST_RESULT",
    "VALIDATION_TEST_RESULT",
    "TRANSFORM_SPECIFICATION",
    "METADATA",
}
OBLIGATION_STATUSES = {
    "PASS_DIRECT_PUBLIC",
    "PASS_DERIVED_IDENTIFIABLE",
    "FAIL_EMPIRICAL",
    "UNSCOREABLE_MISSING_DATA",
    "BLOCKED_BY_PARENT",
    "NOT_RUN_AVAILABLE",
    "THEOREM_ONLY",
}
DERIVATION_MODES = {"NONE", "DIRECT", "FROZEN_IDENTIFIABLE_TRANSFORM"}
PREDICTION_STATUSES = {"PASS", "FAIL", "NOT_RUN"}
BLOCKER_CLASSES = {
    "UNAVAILABLE_PUBLIC_DATA",
    "THEORY_CHOICE",
    "IDENTIFIABILITY",
    "MODEL_FAILURE",
    "UNEXECUTED_AVAILABLE_TEST",
}

SHA256 = re.compile(r"^[0-9a-f]{64}$")


class ProofFrontierRefusal(ValueError):
    """The supplied frontier violates the closed public-data proof policy."""


def _refuse(message: str) -> None:
    raise ProofFrontierRefusal("PROOF FRONTIER REFUSES: " + message)


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


def _reject_constant(value: str) -> None:
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


def _optional_nonempty(value: Any, field: str) -> str | None:
    if value is None:
        return None
    return _nonempty(value, field)


def _boolean(value: Any, field: str) -> bool:
    if not isinstance(value, bool):
        _refuse(f"{field} must be boolean")
    return value


def _digest(value: Any, field: str) -> str:
    if not isinstance(value, str) or not SHA256.fullmatch(value):
        _refuse(f"{field} must be a lowercase SHA-256")
    return value


def _unique_strings(value: Any, field: str, *, allow_empty: bool = True) -> list[str]:
    if not isinstance(value, list) or (not allow_empty and not value):
        _refuse(f"{field} must be a unique string list")
    if any(not isinstance(item, str) or not item.strip() for item in value):
        _refuse(f"{field} contains an empty or non-string value")
    if len(value) != len(set(value)):
        _refuse(f"{field} contains duplicates")
    return value


def _enum(value: Any, allowed: set[str], field: str) -> str:
    text = _nonempty(value, field)
    if text not in allowed:
        _refuse(f"{field} is outside the closed registry")
    return text


def _safe_relative(root: Path, value: Any, field: str) -> Path:
    raw = _nonempty(value, field)
    if unicodedata.normalize("NFC", raw) != raw or "\\" in raw:
        _refuse(f"{field} is not a canonical forward-slash path")
    relative = Path(raw)
    if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
        _refuse(f"{field} must be a safe relative path")
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            _refuse(f"{field} traverses a symlink")
    if not current.is_file():
        _refuse(f"{field} does not name a regular file")
    try:
        resolved = current.resolve(strict=True)
        resolved.relative_to(root.resolve(strict=True))
    except (OSError, ValueError):
        _refuse(f"{field} escapes the frontier root")
    return resolved


def _check_finite(value: Any, field: str = "root") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        _refuse(f"{field} contains a nonfinite number")
    if isinstance(value, list):
        for index, item in enumerate(value):
            _check_finite(item, f"{field}[{index}]")
    elif isinstance(value, dict):
        for key, item in value.items():
            _check_finite(item, f"{field}.{key}")


@dataclass(frozen=True)
class ProofFrontier:
    manifest_path: Path
    manifest_sha256: str
    manifest: dict[str, Any]
    proof_states: dict[str, str]
    theory_states: dict[str, str]
    missing_data_blockers: tuple[dict[str, Any], ...]
    available_execution_frontier: tuple[str, ...]

    def certificate(self) -> dict[str, Any]:
        paper_tests = []
        for theory in self.manifest["theory_tracks"]:
            for prediction in theory["predictions"]:
                classes = {
                    _object_map(self.manifest)[object_id]["evidence_class"]
                    for object_id in prediction["evidence_object_ids"]
                }
                if "PAPER_LEVEL_SUMMARY" in classes:
                    paper_tests.append(
                        {
                            "theory_id": theory["theory_id"],
                            "prediction_id": prediction["prediction_id"],
                            "status": prediction["status"],
                            "proof_eligible": False,
                        }
                    )
        return {
            "schema": CERTIFICATE_SCHEMA,
            "frontier_id": self.manifest["frontier_id"],
            "profile_id": self.manifest["profile_id"],
            "manifest_sha256": self.manifest_sha256,
            "contract_sha256": sha256_file(Path(__file__).resolve()),
            "goal": GOAL,
            "public_data_only": True,
            "proof_states": dict(self.proof_states),
            "theory_states": dict(self.theory_states),
            "missing_data_blockers": [dict(row) for row in self.missing_data_blockers],
            "available_execution_frontier": list(self.available_execution_frontier),
            "paper_level_theory_tests": paper_tests,
            "paper_level_proof_credit": 0,
            "frontier_scientific_weight": "ZERO",
            "scientific_readiness_authorized": False,
            "authoritative_scorer_required": True,
            "obligation_statuses_are_input_assertions": True,
            "theory_states_authoritative": False,
            "authoritative_proof_outputs": {
                proof["proof_id"]: AUTHORITATIVE_OUTPUT
                for proof in self.manifest["proofs"]
            },
            "proofs_earned": [],
            "external_actions_authorized": False,
            "external_actions_executed": False,
        }


def _object_map(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["object_id"]: row for row in manifest["objects"]}


def _proof_state(required: list[dict[str, Any]]) -> str:
    statuses = {row["status"] for row in required}
    if "FAIL_EMPIRICAL" in statuses:
        return "NONAUTHORITATIVE_REPORTED_FAILURE_REQUIRES_AUTHORITATIVE_SCORER"
    if statuses and statuses <= {"PASS_DIRECT_PUBLIC", "PASS_DERIVED_IDENTIFIABLE"}:
        return "NONAUTHORITATIVE_INPUT_COMPLETE_REQUIRES_AUTHORITATIVE_SCORER"
    if "UNSCOREABLE_MISSING_DATA" in statuses:
        return "BLOCKED_MISSING_DATA"
    return "IN_PROGRESS"


def _theory_state(
    theory: dict[str, Any],
    datasets: dict[str, dict[str, Any]],
    objects: dict[str, dict[str, Any]],
) -> str:
    predictions = theory["predictions"]
    statuses = {row["status"] for row in predictions}
    if "FAIL" in statuses:
        return "REPORTED_REFUTATION_REQUIRES_SCIENTIFIC_AUDIT"
    if not theory["mathematics_closed"] or not theory["adversarial_checks_pass"]:
        return "INCOMPLETE_THEORY"
    if not predictions:
        return "INCOMPLETE_THEORY"

    unexecuted_available = False
    for prediction in predictions:
        if prediction["status"] != "NOT_RUN":
            continue
        if not prediction["dataset_ids"]:
            unexecuted_available = True
        elif any(
            datasets[dataset_id]["availability"] != "UNAVAILABLE_PUBLIC_DATA"
            for dataset_id in prediction["dataset_ids"]
        ):
            unexecuted_available = True

    blockers = {row["blocker_class"] for row in theory["remaining_blockers"]}
    domains = theory["independent_event_level_domain_ids"]
    ancestries = {datasets[dataset_id]["acquisition_ancestry"] for dataset_id in domains}
    heldout = theory["held_out_validation_result_object_ids"]
    heldout_valid = bool(heldout) and all(
        objects[object_id]["role"] == "VALIDATION_TEST_RESULT"
        and datasets[objects[object_id]["dataset_id"]]["role"] == "VALIDATION"
        for object_id in heldout
    )

    substantial = (
        not unexecuted_available
        and blockers <= {"UNAVAILABLE_PUBLIC_DATA"}
        and len(domains) >= 2
        and len(ancestries) >= 2
        and heldout_valid
    )
    return (
        "CANDIDATE_SUBSTANTIAL_SUPPORT_REQUIRES_SCIENTIFIC_AUDIT"
        if substantial
        else "PARTIAL_SUPPORT_ONLY"
    )


def load_proof_frontier(manifest_path: str | Path) -> ProofFrontier:
    path = Path(manifest_path)
    if path.is_symlink() or not path.is_file():
        _refuse("manifest must be a non-symlink regular file")
    root = path.absolute().parent.resolve(strict=True)
    try:
        path.resolve(strict=True).relative_to(root)
        manifest_bytes = path.read_bytes()
        manifest = json.loads(
            manifest_bytes,
            object_pairs_hook=_unique_json_object,
            parse_constant=_reject_constant,
        )
    except ProofFrontierRefusal:
        raise
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        _refuse(f"manifest is not readable strict JSON: {exc}")
    _check_finite(manifest)
    _closed(manifest, ROOT_KEYS, "root")
    if manifest["schema"] != SCHEMA:
        _refuse("root.schema mismatch")
    frontier_id = _nonempty(manifest["frontier_id"], "root.frontier_id")
    profile_id = _enum(manifest["profile_id"], PROFILES, "root.profile_id")
    if profile_id == "SYNTHETIC_TEST_PROFILE_V001" and not frontier_id.startswith("SYNTHETIC_"):
        _refuse("synthetic profile is restricted to an explicit SYNTHETIC_ frontier")

    policy = _closed(manifest["policy"], POLICY_KEYS, "root.policy")
    if policy["goal"] != GOAL:
        _refuse("policy.goal must remain PROOF")
    if not _boolean(policy["public_data_only"], "policy.public_data_only"):
        _refuse("policy.public_data_only must be true")
    if not _boolean(policy["validation_preserved"], "policy.validation_preserved"):
        _refuse("policy.validation_preserved must be true")
    if policy["theory_fallback_rule"] != THEORY_FALLBACK:
        _refuse("policy.theory_fallback_rule mismatch")
    policy_object_ids = _unique_strings(
        policy["policy_object_ids"], "policy.policy_object_ids", allow_empty=False
    )

    actions = _closed(manifest["external_actions"], EXTERNAL_ACTION_KEYS, "external_actions")
    for key, value in actions.items():
        if _boolean(value, f"external_actions.{key}"):
            _refuse(f"external_actions.{key} must remain false")

    if not isinstance(manifest["objects"], list) or not manifest["objects"]:
        _refuse("objects must be a nonempty list")
    objects: dict[str, dict[str, Any]] = {}
    paths: dict[Path, str] = {}
    for index, item in enumerate(manifest["objects"]):
        field = f"objects[{index}]"
        obj = _closed(item, OBJECT_KEYS, field)
        object_id = _nonempty(obj["object_id"], f"{field}.object_id")
        if object_id in objects:
            _refuse(f"duplicate object_id {object_id}")
        evidence_class = _enum(obj["evidence_class"], EVIDENCE_CLASSES, f"{field}.evidence_class")
        dataset_id = _optional_nonempty(obj["dataset_id"], f"{field}.dataset_id")
        stage = _enum(obj["stage"], OBJECT_STAGES, f"{field}.stage")
        role = _enum(obj["role"], OBJECT_ROLES, f"{field}.role")
        file_path = _safe_relative(root, obj["path"], f"{field}.path")
        if file_path in paths:
            _refuse(f"{field}.path aliases object {paths[file_path]}")
        expected = _digest(obj["sha256"], f"{field}.sha256")
        if sha256_file(file_path) != expected:
            _refuse(f"{field}.sha256 mismatch")
        _nonempty(obj["media_type"], f"{field}.media_type")
        if evidence_class == "PAPER_LEVEL_SUMMARY" and stage not in {"PAPER_SUMMARY", "RESULT"}:
            _refuse(f"{field} paper evidence must remain paper-summary/result stage")
        if role == "PAPER_SUMMARY_DATA" and evidence_class != "PAPER_LEVEL_SUMMARY":
            _refuse(f"{field} paper role/evidence mismatch")
        if role in {"PUBLIC_DATA_POLICY", "POLICY_AMENDMENT"}:
            if dataset_id is not None or stage != "POLICY" or evidence_class != "THEORY_ONLY":
                _refuse(f"{field} policy object typing mismatch")
        elif dataset_id is None:
            _refuse(f"{field}.dataset_id is required for non-policy objects")
        objects[object_id] = obj
        paths[file_path] = object_id

    for object_id in policy_object_ids:
        if object_id not in objects or objects[object_id]["role"] not in {
            "PUBLIC_DATA_POLICY",
            "POLICY_AMENDMENT",
        }:
            _refuse("policy_object_ids must bind policy-role objects")
    if {objects[object_id]["role"] for object_id in policy_object_ids} != {
        "PUBLIC_DATA_POLICY",
        "POLICY_AMENDMENT",
    }:
        _refuse("both public-data policy and evidence-tier amendment are required")

    if not isinstance(manifest["datasets"], list) or not manifest["datasets"]:
        _refuse("datasets must be a nonempty list")
    datasets: dict[str, dict[str, Any]] = {}
    ancestry_owner: dict[str, str] = {}
    namespace_owner: dict[str, str] = {}
    for index, item in enumerate(manifest["datasets"]):
        field = f"datasets[{index}]"
        dataset = _closed(item, DATASET_KEYS, field)
        dataset_id = _nonempty(dataset["dataset_id"], f"{field}.dataset_id")
        if dataset_id in datasets:
            _refuse(f"duplicate dataset_id {dataset_id}")
        role = _enum(dataset["role"], DATASET_ROLES, f"{field}.role")
        availability = _enum(dataset["availability"], AVAILABILITIES, f"{field}.availability")
        evidence_class = _enum(dataset["evidence_class"], EVIDENCE_CLASSES, f"{field}.evidence_class")
        if not _boolean(dataset["public"], f"{field}.public"):
            _refuse(f"{field}.public must be true")
        opened = _boolean(dataset["opened"], f"{field}.opened")
        if opened and availability != "AVAILABLE_PUBLIC":
            _refuse(f"{field} cannot be opened when public bytes are unavailable/reserved")
        _nonempty(dataset["source_id"], f"{field}.source_id")
        _nonempty(dataset["acquisition_id"], f"{field}.acquisition_id")
        namespace = _nonempty(dataset["physical_unit_namespace"], f"{field}.physical_unit_namespace")
        ancestry = _nonempty(dataset["acquisition_ancestry"], f"{field}.acquisition_ancestry")
        prior_use = _enum(dataset["prior_use"], PRIOR_USES, f"{field}.prior_use")
        independent = _unique_strings(
            dataset["independent_from_dataset_ids"],
            f"{field}.independent_from_dataset_ids",
        )
        frozen = _boolean(
            dataset["frozen_before_response_access"],
            f"{field}.frozen_before_response_access",
        )
        role_locked = _boolean(
            dataset["role_locked_before_response_access"],
            f"{field}.role_locked_before_response_access",
        )
        upstream = _enum(dataset["upstream_status"], UPSTREAM_STATUSES, f"{field}.upstream_status")
        qualification_ids = _unique_strings(
            dataset["qualification_object_ids"], f"{field}.qualification_object_ids"
        )
        measured = _unique_strings(dataset["measured_channel_ids"], f"{field}.measured_channel_ids")
        missing = _unique_strings(dataset["missing_channel_ids"], f"{field}.missing_channel_ids")
        if set(measured) & set(missing):
            _refuse(f"{field} measures and misses the same channel")
        object_ids = _unique_strings(dataset["object_ids"], f"{field}.object_ids", allow_empty=False)
        for object_id in qualification_ids + object_ids:
            if object_id not in objects:
                _refuse(f"{field} references unknown object {object_id}")
            if objects[object_id]["dataset_id"] != dataset_id:
                _refuse(f"{field} references an object owned by another dataset")
            if objects[object_id]["evidence_class"] != evidence_class:
                _refuse(f"{field} object evidence class mismatch")
        if not set(qualification_ids) <= set(object_ids):
            _refuse(f"{field}.qualification_object_ids must be in object_ids")
        if upstream == "QUALIFIED_EVENT_LEVEL":
            if evidence_class != "REAL_WORLD_EVENT_LEVEL" or not qualification_ids:
                _refuse(f"{field} qualified event-level typing is incomplete")
            if any(objects[x]["role"] != "EVENT_LEVEL_QUALIFICATION_CERTIFICATE" for x in qualification_ids):
                _refuse(f"{field} qualification objects have wrong role")
        if evidence_class == "PAPER_LEVEL_SUMMARY" and upstream != "PAPER_SUMMARY_ONLY":
            _refuse(f"{field} paper evidence cannot claim event qualification")
        if evidence_class == "SYNTHETIC_TEST_ONLY" and upstream != "SYNTHETIC_ONLY":
            _refuse(f"{field} synthetic evidence typing mismatch")
        if evidence_class == "THEORY_ONLY" and upstream != "THEORY_ONLY":
            _refuse(f"{field} theory evidence typing mismatch")
        if role == "THEORY_ONLY" and evidence_class not in {"THEORY_ONLY", "SYNTHETIC_TEST_ONLY"}:
            _refuse(f"{field} theory-only role has empirical evidence class")
        if role == "VALIDATION":
            if evidence_class != "REAL_WORLD_EVENT_LEVEL":
                _refuse(f"{field} validation must be real-world event-level")
            if prior_use != "NONE" or not frozen or not role_locked or not independent:
                _refuse(f"{field} validation custody is not preserved")
        if role == "VALIDATION" or evidence_class == "REAL_WORLD_EVENT_LEVEL":
            if ancestry in ancestry_owner:
                _refuse(f"{field} duplicates acquisition ancestry of {ancestry_owner[ancestry]}")
            if namespace in namespace_owner:
                _refuse(f"{field} duplicates physical-unit namespace of {namespace_owner[namespace]}")
            ancestry_owner[ancestry] = dataset_id
            namespace_owner[namespace] = dataset_id
        datasets[dataset_id] = dataset

    for dataset_id, dataset in datasets.items():
        for other in dataset["independent_from_dataset_ids"]:
            if other not in datasets or other == dataset_id:
                _refuse(f"dataset {dataset_id} has invalid independence reference")
        if dataset["role"] == "VALIDATION":
            nonvalidation = {key for key, row in datasets.items() if row["role"] != "VALIDATION"}
            if not set(dataset["independent_from_dataset_ids"]) & nonvalidation:
                _refuse(f"validation dataset {dataset_id} lacks a development/calibration comparator")
    for object_id, obj in objects.items():
        if obj["dataset_id"] is None:
            continue
        if obj["dataset_id"] not in datasets:
            _refuse(f"object {object_id} names an unknown dataset")
        if object_id not in datasets[obj["dataset_id"]]["object_ids"]:
            _refuse(f"object {object_id} is not closed into its dataset object list")

    if not isinstance(manifest["proofs"], list) or not manifest["proofs"]:
        _refuse("proofs must be a nonempty list")
    proofs: dict[str, dict[str, Any]] = {}
    for index, item in enumerate(manifest["proofs"]):
        field = f"proofs[{index}]"
        proof = _closed(item, PROOF_KEYS, field)
        proof_id = _enum(proof["proof_id"], PROOF_IDS, f"{field}.proof_id")
        if proof_id in proofs:
            _refuse(f"duplicate proof_id {proof_id}")
        _nonempty(proof["target"], f"{field}.target")
        _unique_strings(
            proof["required_obligation_ids"],
            f"{field}.required_obligation_ids",
            allow_empty=False,
        )
        if proof["authoritative_output"] != AUTHORITATIVE_OUTPUT:
            _refuse(f"{field}.authoritative_output must remain NO_PROOF_OUTPUT")
        proofs[proof_id] = proof
    if set(proofs) != PROOF_IDS:
        _refuse("both URF and GE proof targets are required")

    if not isinstance(manifest["obligations"], list) or not manifest["obligations"]:
        _refuse("obligations must be a nonempty list")
    obligations: dict[str, dict[str, Any]] = {}
    for index, item in enumerate(manifest["obligations"]):
        field = f"obligations[{index}]"
        obligation = _closed(item, OBLIGATION_KEYS, field)
        obligation_id = _nonempty(obligation["obligation_id"], f"{field}.obligation_id")
        if obligation_id in obligations:
            _refuse(f"duplicate obligation_id {obligation_id}")
        proof_id = _enum(obligation["proof_id"], PROOF_IDS, f"{field}.proof_id")
        _nonempty(obligation["label"], f"{field}.label")
        required = _boolean(obligation["required"], f"{field}.required")
        dependencies = _unique_strings(obligation["dependencies"], f"{field}.dependencies")
        status = _enum(obligation["status"], OBLIGATION_STATUSES, f"{field}.status")
        dataset_ids = _unique_strings(obligation["dataset_ids"], f"{field}.dataset_ids")
        evidence_ids = _unique_strings(
            obligation["evidence_object_ids"], f"{field}.evidence_object_ids"
        )
        missing_ids = _unique_strings(
            obligation["missing_channel_ids"], f"{field}.missing_channel_ids"
        )
        _unique_strings(obligation["theory_track_ids"], f"{field}.theory_track_ids")
        for dataset_id in dataset_ids:
            if dataset_id not in datasets:
                _refuse(f"{field} references unknown dataset {dataset_id}")
        for object_id in evidence_ids:
            if object_id not in objects:
                _refuse(f"{field} references unknown evidence object {object_id}")
            owner = objects[object_id]["dataset_id"]
            if owner not in dataset_ids:
                _refuse(f"{field} evidence object owner is not in dataset_ids")
        derivation = _closed(obligation["derivation"], DERIVATION_KEYS, f"{field}.derivation")
        mode = _enum(derivation["mode"], DERIVATION_MODES, f"{field}.derivation.mode")
        transform_id = _optional_nonempty(
            derivation["transform_object_id"], f"{field}.derivation.transform_object_id"
        )
        derivation_flags = {
            key: _boolean(derivation[key], f"{field}.derivation.{key}")
            for key in {
                "frozen_before_response_access",
                "identifiable_on_declared_support",
                "uncertainty_propagated",
                "uses_tested_response_as_control",
            }
        }
        if mode == "FROZEN_IDENTIFIABLE_TRANSFORM":
            if transform_id not in objects or objects[transform_id]["role"] != "TRANSFORM_SPECIFICATION":
                _refuse(f"{field} derived pass lacks a transform specification")
            if transform_id not in evidence_ids:
                _refuse(f"{field} transform specification is not evidence-bound")
        elif transform_id is not None:
            _refuse(f"{field} non-derived obligation cannot name a transform")

        if status in {"PASS_DIRECT_PUBLIC", "PASS_DERIVED_IDENTIFIABLE", "FAIL_EMPIRICAL"}:
            if not dataset_ids or not evidence_ids or missing_ids:
                _refuse(f"{field} passing obligation has incomplete evidence or missing channels")
            for dataset_id in dataset_ids:
                dataset = datasets[dataset_id]
                if (
                    dataset["evidence_class"] != "REAL_WORLD_EVENT_LEVEL"
                    or dataset["availability"] != "AVAILABLE_PUBLIC"
                    or not dataset["opened"]
                    or dataset["upstream_status"] != "QUALIFIED_EVENT_LEVEL"
                    or dataset["role"] == "THEORY_ONLY"
                ):
                    _refuse(f"{field} pass uses non-proof-grade data")
            if any(objects[object_id]["evidence_class"] != "REAL_WORLD_EVENT_LEVEL" for object_id in evidence_ids):
                _refuse(f"{field} pass uses non-event-level evidence")
            empirical_roles = {
                "RAW_EVENT_DATA",
                "CALIBRATION_DATA",
                "DERIVED_EVENT_DATA",
                "VALIDATION_TEST_RESULT",
            }
            if not any(objects[object_id]["role"] in empirical_roles for object_id in evidence_ids):
                _refuse(f"{field} pass lacks an empirical event object")
            if status == "PASS_DIRECT_PUBLIC" and mode != "DIRECT":
                _refuse(f"{field} direct pass derivation mismatch")
            if status == "PASS_DERIVED_IDENTIFIABLE" and mode != "FROZEN_IDENTIFIABLE_TRANSFORM":
                _refuse(f"{field} derived pass mode mismatch")
            if status in {"PASS_DERIVED_IDENTIFIABLE", "FAIL_EMPIRICAL"} and mode == "FROZEN_IDENTIFIABLE_TRANSFORM":
                if not all(
                    derivation_flags[key]
                    for key in {
                        "frozen_before_response_access",
                        "identifiable_on_declared_support",
                        "uncertainty_propagated",
                    }
                ) or derivation_flags["uses_tested_response_as_control"]:
                    _refuse(f"{field} derivation is not lawful")
            if status == "FAIL_EMPIRICAL" and mode not in {
                "DIRECT",
                "FROZEN_IDENTIFIABLE_TRANSFORM",
            }:
                _refuse(f"{field} empirical failure lacks a lawful direct/derived instrument")
        elif status == "UNSCOREABLE_MISSING_DATA":
            if not missing_ids or evidence_ids:
                _refuse(f"{field} missing-data status must name missing channels and no scored evidence")
            if mode != "NONE":
                _refuse(f"{field} missing-data status cannot claim a derivation")
        elif status == "NOT_RUN_AVAILABLE":
            if missing_ids or not dataset_ids:
                _refuse(f"{field} available test must bind datasets and no missing channel")
            if any(
                datasets[x]["availability"] != "AVAILABLE_PUBLIC"
                or not datasets[x]["opened"]
                or datasets[x]["role"] == "VALIDATION"
                for x in dataset_ids
            ):
                _refuse(f"{field} is not an opened, non-validation executable dataset")
        elif status == "THEOREM_ONLY":
            if required or any(
                datasets[x]["evidence_class"] == "REAL_WORLD_EVENT_LEVEL" for x in dataset_ids
            ):
                _refuse(f"{field} theorem-only cannot be required empirical proof")
        elif status == "BLOCKED_BY_PARENT" and not dependencies:
            _refuse(f"{field} blocked status lacks a parent")
        obligations[obligation_id] = obligation

    # Close the obligation DAG and prevent parent laundering.
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(obligation_id: str) -> None:
        if obligation_id in visiting:
            _refuse("obligation dependency graph contains a cycle")
        if obligation_id in visited:
            return
        visiting.add(obligation_id)
        row = obligations[obligation_id]
        for parent_id in row["dependencies"]:
            if parent_id not in obligations:
                _refuse(f"obligation {obligation_id} references unknown parent {parent_id}")
            parent = obligations[parent_id]
            if parent["proof_id"] != row["proof_id"]:
                _refuse(f"obligation {obligation_id} crosses proof DAGs")
            visit(parent_id)
        visiting.remove(obligation_id)
        visited.add(obligation_id)

    for obligation_id in obligations:
        visit(obligation_id)
    pass_statuses = {"PASS_DIRECT_PUBLIC", "PASS_DERIVED_IDENTIFIABLE"}
    empirical_terminal_statuses = pass_statuses | {"FAIL_EMPIRICAL"}
    for obligation_id, row in obligations.items():
        parent_statuses = {obligations[x]["status"] for x in row["dependencies"]}
        if row["status"] in empirical_terminal_statuses and not parent_statuses <= pass_statuses:
            _refuse(f"obligation {obligation_id} scores an unpassed parent")
        if row["status"] == "NOT_RUN_AVAILABLE" and not parent_statuses <= pass_statuses:
            _refuse(f"obligation {obligation_id} advertises work before its parents pass")
        if row["status"] == "BLOCKED_BY_PARENT" and parent_statuses <= pass_statuses:
            _refuse(f"obligation {obligation_id} claims a nonexistent parent block")

    for proof_id, proof in proofs.items():
        required_ids = proof["required_obligation_ids"]
        for obligation_id in required_ids:
            if obligation_id not in obligations:
                _refuse(f"proof {proof_id} references unknown obligation {obligation_id}")
            row = obligations[obligation_id]
            if row["proof_id"] != proof_id or not row["required"]:
                _refuse(f"proof {proof_id} required-obligation typing mismatch")
        actual_required = {
            key for key, row in obligations.items() if row["proof_id"] == proof_id and row["required"]
        }
        if set(required_ids) != actual_required:
            _refuse(f"proof {proof_id} omits or invents required obligations")
        if (
            profile_id == "PUBLIC_TWO_PROOF_PROFILE_V001"
            and set(required_ids) != CANONICAL_REQUIRED_OBLIGATIONS[proof_id]
        ):
            _refuse(f"proof {proof_id} does not implement the canonical public profile")

    if not isinstance(manifest["theory_tracks"], list):
        _refuse("theory_tracks must be a list")
    theories: dict[str, dict[str, Any]] = {}
    for index, item in enumerate(manifest["theory_tracks"]):
        field = f"theory_tracks[{index}]"
        theory = _closed(item, THEORY_KEYS, field)
        theory_id = _nonempty(theory["theory_id"], f"{field}.theory_id")
        if theory_id in theories:
            _refuse(f"duplicate theory_id {theory_id}")
        _nonempty(theory["label"], f"{field}.label")
        _boolean(theory["mathematics_closed"], f"{field}.mathematics_closed")
        _boolean(theory["adversarial_checks_pass"], f"{field}.adversarial_checks_pass")
        if not isinstance(theory["predictions"], list):
            _refuse(f"{field}.predictions must be a list")
        prediction_ids: set[str] = set()
        for pindex, item in enumerate(theory["predictions"]):
            pfield = f"{field}.predictions[{pindex}]"
            prediction = _closed(item, PREDICTION_KEYS, pfield)
            prediction_id = _nonempty(prediction["prediction_id"], f"{pfield}.prediction_id")
            if prediction_id in prediction_ids:
                _refuse(f"{field} duplicates prediction {prediction_id}")
            prediction_ids.add(prediction_id)
            status = _enum(prediction["status"], PREDICTION_STATUSES, f"{pfield}.status")
            dataset_ids = _unique_strings(prediction["dataset_ids"], f"{pfield}.dataset_ids")
            evidence_ids = _unique_strings(
                prediction["evidence_object_ids"], f"{pfield}.evidence_object_ids"
            )
            for dataset_id in dataset_ids:
                if dataset_id not in datasets:
                    _refuse(f"{pfield} references unknown dataset {dataset_id}")
            for object_id in evidence_ids:
                if object_id not in objects or objects[object_id]["dataset_id"] not in dataset_ids:
                    _refuse(f"{pfield} evidence ownership mismatch")
            if status in {"PASS", "FAIL"} and (not dataset_ids or not evidence_ids):
                _refuse(f"{pfield} executed result lacks evidence")
            if status == "NOT_RUN" and evidence_ids:
                _refuse(f"{pfield} unexecuted prediction cannot carry result evidence")
        domains = _unique_strings(
            theory["independent_event_level_domain_ids"],
            f"{field}.independent_event_level_domain_ids",
        )
        for dataset_id in domains:
            if dataset_id not in datasets:
                _refuse(f"{field} references unknown independent domain {dataset_id}")
            dataset = datasets[dataset_id]
            if (
                dataset["evidence_class"] != "REAL_WORLD_EVENT_LEVEL"
                or dataset["upstream_status"] != "QUALIFIED_EVENT_LEVEL"
                or dataset["availability"] != "AVAILABLE_PUBLIC"
                or not dataset["opened"]
            ):
                _refuse(f"{field} independent domain is not qualified event-level data")
            if not any(
                prediction["status"] == "PASS"
                and dataset_id in prediction["dataset_ids"]
                for prediction in theory["predictions"]
            ):
                _refuse(f"{field} independent domain lacks a passed prediction")
        heldout = _unique_strings(
            theory["held_out_validation_result_object_ids"],
            f"{field}.held_out_validation_result_object_ids",
        )
        for object_id in heldout:
            if object_id not in objects:
                _refuse(f"{field} references unknown held-out result")
            obj = objects[object_id]
            if (
                obj["role"] != "VALIDATION_TEST_RESULT"
                or obj["evidence_class"] != "REAL_WORLD_EVENT_LEVEL"
                or datasets[obj["dataset_id"]]["role"] != "VALIDATION"
                or not datasets[obj["dataset_id"]]["opened"]
            ):
                _refuse(f"{field} held-out result is not qualified validation")
            if not any(
                prediction["status"] == "PASS"
                and object_id in prediction["evidence_object_ids"]
                for prediction in theory["predictions"]
            ):
                _refuse(f"{field} held-out result is not bound to a passed prediction")
        if not isinstance(theory["remaining_blockers"], list):
            _refuse(f"{field}.remaining_blockers must be a list")
        blocker_ids: set[str] = set()
        for bindex, item in enumerate(theory["remaining_blockers"]):
            bfield = f"{field}.remaining_blockers[{bindex}]"
            blocker = _closed(item, BLOCKER_KEYS, bfield)
            blocker_id = _nonempty(blocker["blocker_id"], f"{bfield}.blocker_id")
            if blocker_id in blocker_ids:
                _refuse(f"{field} duplicates blocker {blocker_id}")
            blocker_ids.add(blocker_id)
            _enum(blocker["blocker_class"], BLOCKER_CLASSES, f"{bfield}.blocker_class")
            _nonempty(blocker["description"], f"{bfield}.description")
            affected = _unique_strings(
                blocker["affected_obligation_ids"], f"{bfield}.affected_obligation_ids"
            )
            missing = _unique_strings(blocker["missing_channel_ids"], f"{bfield}.missing_channel_ids")
            if blocker["blocker_class"] == "UNAVAILABLE_PUBLIC_DATA" and not missing:
                _refuse(f"{bfield} unavailable-data blocker lacks missing channels")
            for obligation_id in affected:
                if obligation_id not in obligations:
                    _refuse(f"{bfield} references unknown obligation")
        theories[theory_id] = theory

    for obligation_id, row in obligations.items():
        for theory_id in row["theory_track_ids"]:
            if theory_id not in theories:
                _refuse(f"obligation {obligation_id} references unknown theory track")

    proof_states = {
        proof_id: _proof_state([obligations[x] for x in proof["required_obligation_ids"]])
        for proof_id, proof in proofs.items()
    }
    theory_states = {
        theory_id: _theory_state(theory, datasets, objects)
        for theory_id, theory in theories.items()
    }

    missing_rows = tuple(
        {
            "obligation_id": row["obligation_id"],
            "proof_id": row["proof_id"],
            "missing_channel_ids": list(row["missing_channel_ids"]),
        }
        for row in manifest["obligations"]
        if row["status"] == "UNSCOREABLE_MISSING_DATA"
    )
    available_frontier = tuple(
        row["obligation_id"]
        for row in manifest["obligations"]
        if row["status"] == "NOT_RUN_AVAILABLE"
        and all(
            datasets[x]["availability"] == "AVAILABLE_PUBLIC"
            and datasets[x]["opened"]
            and datasets[x]["role"] != "VALIDATION"
            for x in row["dataset_ids"]
        )
    )

    return ProofFrontier(
        manifest_path=path.resolve(strict=True),
        manifest_sha256=sha256_bytes(manifest_bytes),
        manifest=manifest,
        proof_states=proof_states,
        theory_states=theory_states,
        missing_data_blockers=missing_rows,
        available_execution_frontier=available_frontier,
    )


def certificate_json(frontier: ProofFrontier) -> str:
    return json.dumps(frontier.certificate(), indent=2, sort_keys=True) + "\n"
