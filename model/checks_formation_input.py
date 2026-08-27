#!/usr/bin/env python3
"""Fixed synthetic checks for the generic two-phase formation-input contract.

Fixtures exercise contract behavior only.  They are never scientific evidence and
the resulting certificates must keep every science authorization false.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
from typing import Any, Callable

import formation_input as contract
from formation_input import (
    FormationRefusal,
    attach_formation_execution,
    load_formation_input,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE_SCHEMA = ROOT / "LANE_T53_D_FORMATION_GATE" / "FORMATION_PROTOCOL_V001.schema.json"


def digest_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def digest(path: Path) -> str:
    return digest_bytes(path.read_bytes())


def write_payload(root: Path, relative: str, payload: bytes | str) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(payload, str):
        path.write_text(payload, encoding="utf-8")
    else:
        path.write_bytes(payload)
    return path


def write_json(root: Path, relative: str, value: Any) -> Path:
    return write_payload(
        root,
        relative,
        json.dumps(value, indent=2, sort_keys=True, separators=(",", ": ")) + "\n",
    )


def file_ref(root: Path, path: Path, object_id: str, *, frozen: bool = False) -> dict[str, Any]:
    ref: dict[str, Any] = {
        "object_id": object_id,
        "path": str(path.relative_to(root)),
        "sha256": digest(path),
        "media_type": "application/json" if path.suffix == ".json" else "text/plain",
    }
    if frozen:
        ref.update({"frozen_at_utc": "2026-01-01T00:00:00Z", "immutable": True})
    return ref


def _resolve(schema: dict[str, Any], root_schema: dict[str, Any]) -> dict[str, Any]:
    if "$ref" not in schema:
        return copy.deepcopy(schema)
    node: Any = root_schema
    for token in schema["$ref"][2:].split("/"):
        node = node[token.replace("~1", "/").replace("~0", "~")]
    merged = _merge_schema(copy.deepcopy(node), {k: v for k, v in schema.items() if k != "$ref"}, root_schema)
    return merged


def _merge_schema(left: dict[str, Any], right: dict[str, Any], root_schema: dict[str, Any]) -> dict[str, Any]:
    left = _resolve(left, root_schema) if "$ref" in left else copy.deepcopy(left)
    right = _resolve(right, root_schema) if "$ref" in right else copy.deepcopy(right)
    result = copy.deepcopy(left)
    for key, value in right.items():
        if key == "required":
            result[key] = list(dict.fromkeys(result.get(key, []) + value))
        elif key == "properties":
            properties = result.setdefault(key, {})
            for name, child in value.items():
                properties[name] = _merge_schema(properties.get(name, {}), child, root_schema)
        elif key == "allOf":
            result[key] = result.get(key, []) + copy.deepcopy(value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def protocol_example(schema: dict[str, Any], root_schema: dict[str, Any]) -> Any:
    schema = _resolve(schema, root_schema)
    unconditional: list[dict[str, Any]] = []
    conditional: list[dict[str, Any]] = []
    for branch in schema.pop("allOf", []):
        (conditional if "if" in branch else unconditional).append(branch)
    for branch in unconditional:
        schema = _merge_schema(schema, branch, root_schema)
    if "const" in schema:
        return copy.deepcopy(schema["const"])
    if "enum" in schema:
        return copy.deepcopy(schema["enum"][0])
    if "oneOf" in schema:
        return protocol_example(schema["oneOf"][0], root_schema)
    value_type = schema.get("type")
    if value_type == "object" or "properties" in schema:
        result: dict[str, Any] = {}
        properties = schema.get("properties", {})
        for key in schema.get("required", []):
            result[key] = protocol_example(properties[key], root_schema)
        for branch in conditional:
            if contract._schema_matches(result, branch["if"], root_schema, "fixture"):
                then_schema = _resolve(branch.get("then", {}), root_schema)
                then_properties = then_schema.get("properties", {})
                for key in then_schema.get("required", []):
                    result[key] = protocol_example(properties[key], root_schema)
                for key, child in then_properties.items():
                    if key in result:
                        combined = _merge_schema(properties.get(key, {}), child, root_schema)
                        result[key] = protocol_example(combined, root_schema)
        return result
    if value_type == "array":
        count = max(1, int(schema.get("minItems", 0)))
        return [protocol_example(schema["items"], root_schema) for _ in range(count)]
    if value_type == "boolean":
        return False
    if value_type == "integer":
        return int(schema.get("exclusiveMinimum", -1)) + 1
    if value_type == "number":
        return float(schema.get("exclusiveMinimum", -1)) + 1.0
    if value_type == "string" or "minLength" in schema or "pattern" in schema:
        if schema.get("format") == "date-time":
            return "2026-01-01T00:00:00Z"
        if schema.get("pattern") == "^[a-f0-9]{64}$":
            return "1" * 64
        return "fixture"
    raise AssertionError(f"cannot synthesize V001 schema fragment {schema}")


def make_protocol(protocol_id: str) -> dict[str, Any]:
    schema = json.loads(BASE_SCHEMA.read_text(encoding="utf-8"))
    protocol = protocol_example(schema, schema)
    protocol["protocol_id"] = protocol_id
    # The first causal-design enum is randomized; its conditional member is required.
    if protocol["causal_design"]["design"] == "RANDOMIZED_INTERVENTION":
        definition = schema["$defs"]["randomizedIntervention"]
        protocol["causal_design"]["randomized_intervention"] = protocol_example(definition, schema)
    # The first interaction-origin enum is engineered; its conditional members are required.
    if protocol["scope"]["interaction_origin"] == "ENGINEERED_EXOGENOUS":
        formation_properties = schema["$defs"]["formationStage"]["properties"]
        for key in ("external_writer", "writer_target_coupling", "writer_off_boundary"):
            protocol["formation_interaction"][key] = protocol_example(formation_properties[key], schema)
    return protocol


def raw_object(
    root: Path,
    dataset_id: str,
    object_id: str,
    role: str,
    stage: str,
    event_id: str | None,
    physical_unit_id: str | None,
    monotonic_ns: int,
    ordinal: int,
) -> dict[str, Any]:
    path = write_json(
        root,
        f"data/{dataset_id}/{object_id}.json",
        {"dataset": dataset_id, "object": object_id, "ordinal": ordinal, "value": ordinal + 1},
    )
    return {
        "object_id": object_id,
        "path": str(path.relative_to(root)),
        "sha256": digest(path),
        "size_bytes": path.stat().st_size,
        "media_type": "application/json",
        "stage": stage,
        "role": role,
        "event_id": event_id,
        "physical_unit_id": physical_unit_id,
        "instrument_id": f"instrument-{dataset_id}",
        "clock_id": f"clock-{dataset_id}",
        "created_at_utc": "2026-02-01T00:00:00Z",
        "monotonic_start_ns": monotonic_ns,
        "monotonic_end_ns": monotonic_ns + 1,
        "channels": [{
            "name": "generic_coordinate",
            "quantity": "dimensionless registered coordinate",
            "unit": "1",
            "uncertainty_value": 0.01,
            "uncertainty_unit": "1",
            "calibration_id": f"cal-{dataset_id}",
        }],
        "vendor_original": True,
        "complete": True,
    }


def make_dataset(root: Path, dataset_id: str, role: str) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    lot_id = f"lot-{dataset_id}"
    objects: list[dict[str, Any]] = []
    by_label: dict[str, dict[str, Any]] = {}
    identity = raw_object(root, dataset_id, f"{dataset_id}-identity", "PHYSICAL_IDENTITY", "IDENTITY", None, None, 1, 0)
    calibration = raw_object(root, dataset_id, f"{dataset_id}-calibration", "UNCERTAINTY_CALIBRATION", "CALIBRATION", None, None, 2, 1)
    objects.extend((identity, calibration))
    by_label.update(identity=identity, calibration=calibration)
    events: list[dict[str, Any]] = []
    units: list[dict[str, Any]] = []
    for event_index, route in enumerate(("WRITE", "SHAM")):
        event_id = f"event-{dataset_id}-{event_index}"
        unit_id = f"unit-{dataset_id}-{event_index}"
        units.append({
            "physical_unit_id": unit_id,
            "unit_kind": "SURFACE_INSTANCE",
            "lot_id": lot_id,
            "parent_specimen_id": f"specimen-{dataset_id}-{event_index}",
            "identity_sha256": digest_bytes(f"identity-{dataset_id}-{event_index}".encode()),
        })
        base_ns = (event_index + 1) * 100_000_000_000
        definitions = [
            ("before", "BEFORE_STATE", "BEFORE", 0),
            ("formation", "FORMATION_INPUT" if route == "WRITE" else "ROUTING_OR_CONTROL", "FORMATION_OR_SHAM", 2),
            ("c", "C_EXT_OFF_ASSAY", "C_EXT_OFF", 4),
            ("g", "G_SURFACE_CLOSE_ASSAY", "G_SURFACE_CLOSE", 5),
            ("hold", "COMMON_HOLD_ENVIRONMENT", "COMMON_HOLD", 6),
            ("read1", "READ_OBSERVATION", "READ", 7),
            ("back1", "READ_BACKACTION_ASSAY", "READ", 7),
            ("read2", "READ_OBSERVATION", "READ", 15),
            ("back2", "READ_BACKACTION_ASSAY", "READ", 15),
        ]
        event_objects: dict[str, dict[str, Any]] = {}
        for local_index, (label, object_role, stage, offset) in enumerate(definitions):
            item = raw_object(
                root,
                dataset_id,
                f"{event_id}-{label}",
                object_role,
                stage,
                event_id,
                unit_id,
                base_ns + offset * 1_000_000_000,
                10 + event_index * 20 + local_index,
            )
            objects.append(item)
            event_objects[label] = item
            by_label[f"{event_index}:{label}"] = item
        minute = event_index * 60
        def utc(second: int) -> str:
            return f"2026-02-01T00:{minute // 60:02d}:{minute % 60 + second:02d}Z"
        event_object_ids = [item["object_id"] for item in event_objects.values()]
        events.append({
            "event_id": event_id,
            "physical_unit_id": unit_id,
            "lot_id": lot_id,
            "block_id": f"block-{dataset_id}",
            "event_value": "event-a" if event_index == 0 else "event-b",
            "route": route,
            "allocation_concealed": True,
            "allocation_commitment_sha256": digest_bytes(f"allocation-{event_id}".encode()),
            "event_generated_at_utc": utc(1),
            "event_generated_monotonic_ns": base_ns + 1_000_000_000,
            "before": {
                "sealed_at_utc": utc(0),
                "sealed_monotonic_ns": base_ns,
                "event_absent_when_sealed": True,
                "raw_object_ids": [event_objects["before"]["object_id"]],
            },
            "formation": {
                "started_at_utc": utc(2),
                "started_monotonic_ns": base_ns + 2_000_000_000,
                "ended_at_utc": utc(3),
                "ended_monotonic_ns": base_ns + 3_000_000_000,
                "physical_target_coupling_measured": True,
                "raw_object_ids": [event_objects["formation"]["object_id"]],
            },
            "c_ext_off": {
                "observed_at_utc": utc(4),
                "observed_monotonic_ns": base_ns + 4_000_000_000,
                "physically_measured": True,
                "raw_object_ids": [event_objects["c"]["object_id"]],
            },
            "g_surface_close": {
                "status": "CERTIFIED",
                "observed_at_utc": utc(5),
                "observed_monotonic_ns": base_ns + 5_000_000_000,
                "separate_from_c_ext_off": True,
                "raw_object_ids": [event_objects["g"]["object_id"]],
                "predicate_ids": [f"pred-{dataset_id}-G_SURFACE_CLOSE"],
            },
            "common_hold": {
                "started_at_utc": utc(5),
                "started_monotonic_ns": base_ns + 5_000_000_000,
                "ended_at_utc": utc(6),
                "ended_monotonic_ns": base_ns + 6_000_000_000,
                "writer_off_throughout": True,
                "raw_object_ids": [event_objects["hold"]["object_id"]],
            },
            "reads": [
                {
                    "read_id": f"read-{event_id}-1",
                    "observed_at_utc": utc(7),
                    "observed_monotonic_ns": base_ns + 7_000_000_000,
                    "query_id": "generic-query-1",
                    "mission_read": False,
                    "raw_object_ids": [event_objects["read1"]["object_id"], event_objects["back1"]["object_id"]],
                },
                {
                    "read_id": f"read-{event_id}-2",
                    "observed_at_utc": utc(15),
                    "observed_monotonic_ns": base_ns + 15_000_000_000,
                    "query_id": "generic-query-2",
                    "mission_read": True,
                    "raw_object_ids": [event_objects["read2"]["object_id"], event_objects["back2"]["object_id"]],
                },
            ],
            "raw_object_ids": event_object_ids,
            "validity_status": "VALID",
        })
    root_inventory = {
        "schema": "WAC_FORMATION_RAW_ROOT_V002",
        "dataset_id": dataset_id,
        "raw_objects": sorted(
            ({"object_id": item["object_id"], "sha256": item["sha256"]} for item in objects),
            key=lambda row: row["object_id"],
        ),
        "derived_objects": [],
    }
    root_path = write_json(root, f"custody/{dataset_id}-raw-root.json", root_inventory)
    acquisition_path = write_json(root, f"custody/{dataset_id}-acquisition.json", {"dataset_id": dataset_id, "physical": True})
    access_path = write_json(root, f"custody/{dataset_id}-access.json", {"dataset_id": dataset_id, "first_access": "2026-02-02T00:00:00Z"})
    expected = [{
        "expected_object_id": item["object_id"],
        "event_id": item["event_id"],
        "physical_unit_id": item["physical_unit_id"],
        "role": item["role"],
    } for item in objects]
    validation = role == "VALIDATION"
    signatures = [{
        "signature_id": f"signature-{dataset_id}-custody",
        "role": "DATA_CUSTODIAN",
        "signer_id": f"custodian-{dataset_id}",
        "signed_sha256": digest(root_path),
        "signed_at_utc": "2026-02-02T00:00:01Z",
        "algorithm": "TEST_FIXTURE_SIGNATURE",
        "signature": f"signature-{dataset_id}",
    }]
    if validation:
        signatures.append({
            "signature_id": f"signature-{dataset_id}-blind",
            "role": "BLIND_LOCK",
            "signer_id": f"blind-{dataset_id}",
            "signed_sha256": digest(access_path),
            "signed_at_utc": "2026-02-02T00:00:02Z",
            "algorithm": "TEST_FIXTURE_SIGNATURE",
            "signature": f"blind-signature-{dataset_id}",
        })
    dataset = {
        "dataset_id": dataset_id,
        "dataset_role": role,
        "evidence_class": "REAL_WORLD_ACTUAL",
        "access_mode": "RESPONSE_BLIND_ARCHIVAL_HOLDOUT",
        "source_id": f"source-{dataset_id}",
        "source_uri": f"urn:wac:fixture:{dataset_id}",
        "acquisition_id": f"acquisition-{dataset_id}",
        "acquisition_cohort_id": f"cohort-{dataset_id}",
        "raw_root_sha256": digest(root_path),
        "raw_root_manifest": file_ref(root, root_path, f"object-{dataset_id}-raw-root"),
        "provenance": {
            "institution": f"institution-{dataset_id}",
            "site_id": f"site-{dataset_id}",
            "operator_ids": [f"operator-{dataset_id}"],
            "apparatus_ids": [f"apparatus-{dataset_id}"],
            "acquired_start_utc": "2025-01-01T00:00:00Z",
            "acquired_end_utc": "2025-01-02T00:00:00Z",
            "acquisition_record": file_ref(root, acquisition_path, f"object-{dataset_id}-acquisition"),
            "physical_origin_documented": True,
        },
        "freeze": {
            "frozen_before_numerical_access": validation,
            "numerical_access_first_utc": "2026-02-02T00:00:00Z",
            "access_log": file_ref(root, access_path, f"object-{dataset_id}-access"),
            "source_outcomes_known_before_freeze": not validation,
            "response_labels_blinded_until_lock": validation,
        },
        "prior_use": "UNUSED_BEFORE_VALIDATION" if validation else "DEVELOPMENT_ONLY",
        "independent_from_dataset_ids": [],
        "no_reuse_attestation": {
            "no_resplit": True,
            "no_rerun": True,
            "no_reanalysis": True,
            "no_copied_data": True,
            "no_overlapping_physical_units": True,
            "no_derived_data_relabel": True,
        },
        "lots": [{"lot_id": lot_id, "cohort_id": f"lot-cohort-{dataset_id}", "description": "identified actual fixture cohort"}],
        "physical_units": units,
        "events": events,
        "expected_objects": expected,
        "raw_objects": objects,
        "derived_objects": [],
        "missing_objects": [],
        "deviations": [],
        "custody_signatures": signatures,
    }
    return dataset, by_label


def make_predicates(dataset: dict[str, Any], labels: dict[str, dict[str, Any]], analysis_ref: dict[str, Any]) -> list[dict[str, Any]]:
    role_sources = {
        "BEFORE_ABSENCE": labels["0:before"],
        "WRITE_COUPLING": labels["0:formation"],
        "SHAM_NO_COUPLING": labels["1:formation"],
        "C_EXT_OFF": labels["0:c"],
        "G_SURFACE_CLOSE": labels["0:g"],
        "COMMON_HOLD_EQUIVALENCE": labels["0:hold"],
        "EVENT_LEVEL_STATE": labels["0:read1"],
        "MISSION_EVENT_INFORMATION": labels["0:read2"],
        "READ_BACKACTION": labels["0:back2"],
        "K_PREDICTION": labels["0:read2"],
        "M_PREDICTION": labels["0:read2"],
        "ADEQUATE_SENSITIVITY": labels["calibration"],
    }
    if dataset["dataset_role"] == "DEVELOPMENT":
        role_sources = {"G_SURFACE_CLOSE": labels["0:g"]}
    predicates: list[dict[str, Any]] = []
    for role, source in role_sources.items():
        predicates.append({
            "predicate_id": f"pred-{dataset['dataset_id']}-{role}",
            "dataset_id": dataset["dataset_id"],
            "predicate_role": role,
            "stage": source["stage"],
            "metric_id": f"metric-{role.lower()}",
            "acceptance": {"operator": "GE_LOWER_BOUND", "lower": 0.0, "upper": None, "unit": "1"},
            "source_object_ids": [source["object_id"]],
            "rule_object_id": analysis_ref["object_id"],
            "rule_sha256": analysis_ref["sha256"],
            "frozen_before_access": True,
        })
    return predicates


def make_mapping(
    root: Path,
    protocol_id: str,
    source_sha256: str,
    datasets: list[dict[str, Any]],
    labels: dict[str, dict[str, dict[str, Any]]],
    frozen_refs: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    field_labels = {
        "S": "identity", "T": "0:before", "X": "0:formation", "P": "identity",
        "D_causal": "0:formation", "A": "1:formation", "I_f": "0:formation",
        "C_psi": "0:formation", "C_EXT_OFF": "0:c", "G_t": "0:g",
        "H": "0:hold", "tau_m": "0:read2", "Q": "0:read1", "Z": "0:read2",
        "U": "calibration", "G_X": "0:formation",
    }
    numeric_fields = {"I_f", "C_psi", "C_EXT_OFF", "G_t", "H", "tau_m", "Q", "Z", "U"}
    rows: list[dict[str, Any]] = []
    for field in sorted(contract.PRE_EXECUTION_MAPPING_FIELDS):
        sources: list[dict[str, Any]] = []
        if field in {"K_theta", "M_phi", "F_policy"}:
            frozen_name = {"K_theta": "prediction", "M_phi": "law", "F_policy": "outcome_map"}[field]
            ref = frozen_refs[frozen_name]
            sources.append({
                "dataset_id": None, "event_id": None, "object_id": ref["object_id"],
                "path": ref["path"], "selector_kind": "WHOLE_OBJECT", "selector": "/",
                "value_type": "IDENTITY", "quantity": f"frozen {field} identity", "unit": None,
                "calibration_id": None, "clock_id": None, "coordinate_transform": "IDENTITY",
            })
        else:
            for dataset in datasets:
                if dataset["dataset_role"] != "VALIDATION":
                    continue
                source = labels[dataset["dataset_id"]][field_labels[field]]
                numeric = field in numeric_fields
                sources.append({
                    "dataset_id": dataset["dataset_id"],
                    "event_id": source["event_id"],
                    "object_id": source["object_id"],
                    "path": source["path"],
                    "selector_kind": "WHOLE_OBJECT",
                    "selector": "/",
                    "value_type": "FLOAT" if numeric else "IDENTITY",
                    "quantity": f"registered generic {field} coordinate",
                    "unit": "1" if numeric else None,
                    "calibration_id": f"cal-{dataset['dataset_id']}" if numeric else None,
                    "clock_id": f"clock-{dataset['dataset_id']}" if numeric else None,
                    "coordinate_transform": "IDENTITY",
                })
        rows.append({"contract_field": field, "sources": sources})
    return {
        "schema": "WAC_FORMATION_NORMALIZED_MAPPING_V002",
        "mapping_id": "generic-normalized-fixture-v002",
        "status": "NORMALIZED_MACHINE_RESOLVABLE",
        "protocol_id": protocol_id,
        "mapping_contains_executable_science": False,
        "surface_specific_source_modification_authorized": False,
        "adapter_id": "INSTALLED_GENERIC_FORMATION_INPUT_V002",
        "adapter_source_sha256": source_sha256,
        "formation_tuple_mapping": rows,
    }


def write_bundle(root: Path, *, surface_shape: str = "distributed mechanical marks") -> tuple[Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    protocol_id = "generic-formation-fixture-v002"
    development, development_labels = make_dataset(root, "development-data", "DEVELOPMENT")
    validation, validation_labels = make_dataset(root, "validation-data", "VALIDATION")
    datasets = [development, validation]
    for dataset in datasets:
        dataset["independent_from_dataset_ids"] = [
            other["dataset_id"] for other in datasets if other is not dataset
        ]

    protocol_path = write_json(root, "frozen/protocol.json", make_protocol(protocol_id))
    law_path = write_json(root, "frozen/law.json", {"schema": "GENERIC_LAW_FIXTURE", "law": "frozen"})
    prediction_path = write_json(root, "frozen/prediction.json", {"schema": "GENERIC_PREDICTION_FIXTURE", "status": "frozen"})
    outcome_path = write_json(root, "frozen/outcome-map.json", {"schema": "GENERIC_OUTCOME_POLICY_FIXTURE", "status": "frozen"})
    analysis_path = write_json(root, "frozen/analysis.json", {"schema": "GENERIC_CLOSED_COMPARATOR_FIXTURE", "operators": sorted(contract.PREDICATE_OPERATORS)})
    environment_path = write_json(root, "frozen/environment.json", {"schema": "GENERIC_ENVIRONMENT_FIXTURE", "runtime": "python"})
    source_path = write_payload(root, "frozen/formation_input.py", Path(contract.__file__).read_bytes())
    release_path = write_json(root, "frozen/release.json", {
        "schema": "WAC_FORMATION_URM_RELEASE_V002",
        "release_id": "generic-public-test-release",
        "formation_input_source_sha256": digest(source_path),
        "environment_sha256": digest(environment_path),
        "public": True,
        "immutable": True,
    })
    frozen_refs = {
        "protocol": file_ref(root, protocol_path, "frozen-protocol", frozen=True),
        "law": file_ref(root, law_path, "frozen-law", frozen=True),
        "prediction": file_ref(root, prediction_path, "frozen-prediction", frozen=True),
        "outcome_map": file_ref(root, outcome_path, "frozen-outcome-map", frozen=True),
        "analysis": file_ref(root, analysis_path, "frozen-analysis", frozen=True),
        "urm_release": file_ref(root, release_path, "frozen-urm-release", frozen=True),
        "urm_source": file_ref(root, source_path, "frozen-urm-source", frozen=True),
        "urm_environment": file_ref(root, environment_path, "frozen-urm-environment", frozen=True),
    }
    labels = {"development-data": development_labels, "validation-data": validation_labels}
    mapping = make_mapping(root, protocol_id, digest(source_path), datasets, labels, frozen_refs)
    mapping_path = write_json(root, "frozen/normalized-mapping.json", mapping)
    frozen_refs["declarative_mapping"] = file_ref(root, mapping_path, "frozen-mapping", frozen=True)

    predicates = make_predicates(development, development_labels, frozen_refs["analysis"])
    predicates += make_predicates(validation, validation_labels, frozen_refs["analysis"])
    manifest = {
        "schema": "WAC_FORMATION_BUNDLE_V002",
        "bundle_id": "generic-two-dataset-fixture",
        "protocol_id": protocol_id,
        "base_contract_version": "FORMATION_PROTOCOL_V001",
        "validation_extension_version": "FORMATION_DATASET_VALIDATION_V002",
        "claim_scope": {
            "claim_id": "fixture-structural-eligibility-only",
            "surface_scope": surface_shape,
            "formation_mode": "engineered or autonomous generic interaction",
            "provenance": "identified event-to-target graph",
            "mission_time": {"value": 10.0, "unit": "s", "standard_uncertainty": 0.0},
            "requested_result": "SCOPED_FORMATION_DATASET_VALIDATION",
            "universal_claim_requested": False,
            "gravity_claim_requested": False,
        },
        "frozen_objects": frozen_refs,
        "datasets": datasets,
        "generic_predicates": predicates,
        "product_reproduction": {
            "status": "NOT_DEMONSTRATED",
            "external_physicist_ids": [],
            "public_urm_release_sha256": None,
            "used_only_public_instructions": False,
            "no_private_help": False,
            "different_real_world_datasets": False,
            "any_and_every_surface_coverage": False,
            "evidence_object_ids": [],
        },
        "signatures": [{
            "signature_id": "root-freeze-signature",
            "role": "PROTOCOL_CUSTODIAN",
            "signer_id": "fixture-custodian",
            "signed_sha256": digest(protocol_path),
            "signed_at_utc": "2026-01-01T00:00:01Z",
            "algorithm": "TEST_FIXTURE_SIGNATURE",
            "signature": "root-fixture-signature",
        }],
    }
    manifest_path = write_json(root, "formation-input.json", manifest)

    execution_files: dict[str, dict[str, Any]] = {}
    for name in ("invocation", "computed_j", "evaluated_falsifier", "result", "stdout", "stderr"):
        path = write_json(root, f"execution/{name}.json", {"role": name, "fixture": True})
        execution_files[name] = file_ref(root, path, f"execution-{name}")
    object_lookup = {
        item["object_id"]: item
        for dataset in datasets
        for item in dataset["raw_objects"] + dataset["derived_objects"]
    }
    measurements = []
    for predicate in predicates:
        measurement_id = predicate["source_object_ids"][0]
        measurements.append({
            "predicate_id": predicate["predicate_id"],
            "dataset_id": predicate["dataset_id"],
            "observed": {"value": 1.0, "unit": "1", "standard_uncertainty": 0.01},
            "source_object_ids": predicate["source_object_ids"],
            "measurement_object_id": measurement_id,
            "measurement_sha256": object_lookup[measurement_id]["sha256"],
        })
    execution = {
        "schema": "WAC_FORMATION_EXECUTION_V002",
        "input_manifest_sha256": digest(manifest_path),
        "release_object_id": frozen_refs["urm_release"]["object_id"],
        "release_sha256": frozen_refs["urm_release"]["sha256"],
        "source_sha256": frozen_refs["urm_source"]["sha256"],
        "environment_sha256": frozen_refs["urm_environment"]["sha256"],
        "analysis_object_id": frozen_refs["analysis"]["object_id"],
        "analysis_sha256": frozen_refs["analysis"]["sha256"],
        "measurements": measurements,
        **execution_files,
    }
    execution_path = write_json(root, "formation-execution.json", execution)
    return manifest_path, execution_path


def mutate_manifest(path: Path, mutation: Callable[[dict[str, Any]], None]) -> None:
    value = json.loads(path.read_text(encoding="utf-8"))
    mutation(value)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def refuses(path: Path, text: str) -> None:
    try:
        load_formation_input(path)
    except FormationRefusal as exc:
        assert text in str(exc), (text, str(exc))
    else:
        raise AssertionError(f"expected refusal containing {text!r}")


def mutate_execution(path: Path, mutation: Callable[[dict[str, Any]], None]) -> None:
    value = json.loads(path.read_text(encoding="utf-8"))
    mutation(value)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def execution_refuses(manifest_path: Path, execution_path: Path, text: str) -> None:
    formation = load_formation_input(manifest_path)
    try:
        attach_formation_execution(formation, execution_path)
    except FormationRefusal as exc:
        assert text in str(exc), (text, str(exc))
    else:
        raise AssertionError(f"expected execution refusal containing {text!r}")


def mutate_frozen_json(
    manifest_path: Path,
    frozen_name: str,
    mutation: Callable[[dict[str, Any]], None],
) -> None:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    ref = manifest["frozen_objects"][frozen_name]
    path = manifest_path.parent / ref["path"]
    value = json.loads(path.read_text(encoding="utf-8"))
    mutation(value)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    ref["sha256"] = digest(path)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    checks = 0
    with tempfile.TemporaryDirectory(prefix="wac-formation-input-") as directory:
        base = Path(directory)

        # 8 schema and custody checks.
        manifest_path, execution_path = write_bundle(base / "s01-positive")
        formation = load_formation_input(manifest_path)
        certificate = formation.certificate()
        assert certificate["validation_dataset_eligible"] is True
        assert certificate["scientific_verdict"] == "NONE_NOT_SCORED"
        assert not any(
            certificate[key]
            for key in (
                "scientific_validation_authorized",
                "record_formation_claim_authorized",
                "scoped_formation_result_authorized",
                "universal_claim_authorized",
                "gravity_claim_authorized",
                "program_completion_authorized",
            )
        )
        execution = attach_formation_execution(formation, execution_path)
        assert execution.certificate()["generic_predicate_status"] == "ALL_FROZEN_GENERIC_PREDICATES_PASS"
        assert execution.certificate()["scientific_validation_authorized"] is False
        assert contract.certificate_json(formation) == contract.certificate_json(load_formation_input(manifest_path))
        checks += 1

        manifest_path, _ = write_bundle(base / "s02-root-closure")
        mutate_manifest(manifest_path, lambda value: value.__setitem__("authoritative_result", True))
        refuses(manifest_path, "key closure")
        checks += 1

        manifest_path, _ = write_bundle(base / "s03-duplicate-json")
        text = manifest_path.read_text(encoding="utf-8")
        manifest_path.write_text(text.replace("{\n", "{\n  \"schema\": \"DUPLICATE\",\n", 1), encoding="utf-8")
        refuses(manifest_path, "duplicate JSON member")
        checks += 1

        manifest_path, _ = write_bundle(base / "s04-traversal")
        mutate_manifest(
            manifest_path,
            lambda value: value["frozen_objects"]["law"].__setitem__("path", "../escape.json"),
        )
        refuses(manifest_path, "safe relative path")
        checks += 1

        manifest_path, _ = write_bundle(base / "s05-byte-mutation")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        prediction = manifest_path.parent / manifest["frozen_objects"]["prediction"]["path"]
        prediction.write_text(prediction.read_text(encoding="utf-8") + "mutated\n", encoding="utf-8")
        refuses(manifest_path, "hash mismatch")
        checks += 1

        manifest_path, _ = write_bundle(base / "s06-size")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][0]["raw_objects"][0].__setitem__(
                "size_bytes", value["datasets"][0]["raw_objects"][0]["size_bytes"] + 1
            ),
        )
        refuses(manifest_path, "size_bytes mismatch")
        checks += 1

        manifest_path, _ = write_bundle(base / "s07-inventory")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][0]["expected_objects"].pop(),
        )
        refuses(manifest_path, "expected/present/missing object inventory is not exact")
        checks += 1

        manifest_path, _ = write_bundle(base / "s08-role-file-alias")
        def alias_frozen(value: dict[str, Any]) -> None:
            protocol_ref = value["frozen_objects"]["protocol"]
            law_ref = value["frozen_objects"]["law"]
            for key in ("path", "sha256", "media_type"):
                law_ref[key] = protocol_ref[key]
        mutate_manifest(manifest_path, alias_frozen)
        refuses(manifest_path, "duplicates resolved path")
        checks += 1

        # 8 physical tuple, event, C-off and G chronology checks.
        manifest_path, _ = write_bundle(base / "t09-actual")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1].__setitem__("evidence_class", "SYNTHETIC_TEST_ONLY"),
        )
        refuses(manifest_path, "validation role requires REAL_WORLD_ACTUAL")
        checks += 1

        manifest_path, _ = write_bundle(base / "t10-target-join")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1]["events"][0].__setitem__("physical_unit_id", "unknown-unit"),
        )
        refuses(manifest_path, "unknown physical unit")
        checks += 1

        manifest_path, _ = write_bundle(base / "t11-before-order")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1]["events"][0]["before"].__setitem__(
                "sealed_monotonic_ns", value["datasets"][1]["events"][0]["event_generated_monotonic_ns"]
            ),
        )
        refuses(manifest_path, "BEFORE does not precede")
        checks += 1

        manifest_path, _ = write_bundle(base / "t12-boundary-alias")
        def alias_boundaries(value: dict[str, Any]) -> None:
            event = value["datasets"][1]["events"][0]
            event["g_surface_close"]["raw_object_ids"] = event["c_ext_off"]["raw_object_ids"][:]
        mutate_manifest(manifest_path, alias_boundaries)
        refuses(manifest_path, "reuses one raw object")
        checks += 1

        manifest_path, _ = write_bundle(base / "t13-g-before-c")
        def reverse_g(value: dict[str, Any]) -> None:
            event = value["datasets"][1]["events"][0]
            event["g_surface_close"]["observed_at_utc"] = "2026-02-01T00:00:03Z"
            event["g_surface_close"]["observed_monotonic_ns"] = 103_000_000_000
            event["common_hold"]["started_at_utc"] = "2026-02-01T00:00:03Z"
            event["common_hold"]["started_monotonic_ns"] = 103_000_000_000
        mutate_manifest(manifest_path, reverse_g)
        refuses(manifest_path, "G precedes C_EXT_OFF")
        checks += 1

        manifest_path, _ = write_bundle(base / "t14-hold-off")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1]["events"][0]["common_hold"].__setitem__(
                "writer_off_throughout", False
            ),
        )
        refuses(manifest_path, "common hold is not writer-off")
        checks += 1

        manifest_path, _ = write_bundle(base / "t15-mission-origin")
        def move_mission(value: dict[str, Any]) -> None:
            read = value["datasets"][1]["events"][0]["reads"][1]
            read["observed_at_utc"] = "2026-02-01T00:00:14Z"
            read["observed_monotonic_ns"] = 114_000_000_000
        mutate_manifest(manifest_path, move_mission)
        refuses(manifest_path, "mission read is outside")
        checks += 1

        manifest_path, _ = write_bundle(base / "t16-arms")
        def remove_sham(value: dict[str, Any]) -> None:
            dataset = value["datasets"][1]
            event = dataset["events"][1]
            event["route"] = "WRITE"
            formation_id = event["formation"]["raw_object_ids"][0]
            next(item for item in dataset["raw_objects"] if item["object_id"] == formation_id)["role"] = "FORMATION_INPUT"
            next(item for item in dataset["expected_objects"] if item["expected_object_id"] == formation_id)["role"] = "FORMATION_INPUT"
        mutate_manifest(manifest_path, remove_sham)
        refuses(manifest_path, "requires both WRITE and SHAM")
        checks += 1

        # 6 pre-access-freeze and D<->V dataset-disjointness checks.
        manifest_path, _ = write_bundle(base / "d17-freeze")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1]["freeze"].__setitem__("frozen_before_numerical_access", False),
        )
        refuses(manifest_path, "not frozen before numerical access")
        checks += 1

        for name, key in (("source", "source_id"), ("acquisition", "acquisition_id")):
            manifest_path, _ = write_bundle(base / f"d18-{name}")
            mutate_manifest(
                manifest_path,
                lambda value, key=key: value["datasets"][1].__setitem__(key, value["datasets"][0][key]),
            )
            refuses(manifest_path, "dataset overlap")
        checks += 1

        manifest_path, _ = write_bundle(base / "d19-raw-root")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1].__setitem__(
                "raw_root_sha256", value["datasets"][0]["raw_root_sha256"]
            ),
        )
        refuses(manifest_path, "dataset overlap")
        checks += 1

        manifest_path, _ = write_bundle(base / "d20-specimen")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1]["physical_units"][0].__setitem__(
                "physical_unit_id", value["datasets"][0]["physical_units"][0]["physical_unit_id"]
            ),
        )
        refuses(manifest_path, "dataset overlap")
        checks += 1

        manifest_path, _ = write_bundle(base / "d21-event")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1]["events"][0].__setitem__(
                "event_id", value["datasets"][0]["events"][0]["event_id"]
            ),
        )
        refuses(manifest_path, "dataset overlap")
        checks += 1

        manifest_path, _ = write_bundle(base / "d22-content")
        mutate_manifest(
            manifest_path,
            lambda value: value["datasets"][1]["raw_objects"][0].__setitem__(
                "sha256", value["datasets"][0]["raw_objects"][0]["sha256"]
            ),
        )
        refuses(manifest_path, "dataset overlap")
        checks += 1

        # 5 K/M/F-policy and two-phase I/O identity checks.
        manifest_path, _ = write_bundle(base / "k23-rule-hash")
        mutate_manifest(
            manifest_path,
            lambda value: value["generic_predicates"][0].__setitem__("rule_sha256", "0" * 64),
        )
        refuses(manifest_path, "rule_sha256 mismatch")
        checks += 1

        manifest_path, _ = write_bundle(base / "k24-unresolved-mapping")
        mutate_frozen_json(
            manifest_path,
            "declarative_mapping",
            lambda value: value["formation_tuple_mapping"][0]["sources"][0].__setitem__(
                "selector", "<event_uuid>"
            ),
        )
        refuses(manifest_path, "executable or unresolved mapping syntax")
        checks += 1

        manifest_path, _ = write_bundle(base / "k25-circular-input")
        mutate_manifest(
            manifest_path,
            lambda value: value.__setitem__("record_formation_claim_authorized", True),
        )
        refuses(manifest_path, "key closure")
        checks += 1

        manifest_path, execution_path = write_bundle(base / "k26-execution-input")
        mutate_execution(execution_path, lambda value: value.__setitem__("input_manifest_sha256", "0" * 64))
        execution_refuses(manifest_path, execution_path, "exact input manifest")
        checks += 1

        for name, mutation, expected in (
            (
                "source",
                lambda value: value.__setitem__("source_sha256", "0" * 64),
                "does not match the frozen input",
            ),
            (
                "measurement",
                lambda value: value["measurements"][0].__setitem__("measurement_sha256", "0" * 64),
                "does not bind the named object",
            ),
        ):
            manifest_path, execution_path = write_bundle(base / f"k27-{name}")
            mutate_execution(execution_path, mutation)
            execution_refuses(manifest_path, execution_path, expected)
        checks += 1

        # 5 no-promotion, genericity, failed-G and sensitivity controls.
        manifest_path, execution_path = write_bundle(base / "n28-generic-fail")
        def fail_validation_predicate(value: dict[str, Any]) -> None:
            measurement = next(
                item for item in value["measurements"]
                if item["dataset_id"] == "validation-data"
            )
            measurement["observed"]["value"] = -1.0
        mutate_execution(execution_path, fail_validation_predicate)
        result = attach_formation_execution(load_formation_input(manifest_path), execution_path).certificate()
        assert result["generic_predicate_status"] == "ONE_OR_MORE_FROZEN_GENERIC_PREDICATES_FAIL"
        assert result["scientific_validation_authorized"] is False
        checks += 1

        manifest_path, execution_path = write_bundle(base / "n29-sensitivity")
        execution_value = json.loads(execution_path.read_text(encoding="utf-8"))
        sensitivity = next(
            item for item in execution_value["measurements"]
            if item["predicate_id"].endswith("ADEQUATE_SENSITIVITY")
        )
        sensitivity["observed"]["value"] = -1.0
        execution_path.write_text(json.dumps(execution_value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        result = attach_formation_execution(load_formation_input(manifest_path), execution_path).certificate()
        assert result["generic_predicate_status"] == "ONE_OR_MORE_FROZEN_GENERIC_PREDICATES_FAIL"
        assert result["record_formation_claim_authorized"] is False
        checks += 1

        manifest_path, _ = write_bundle(base / "n30-failed-g")
        def fail_one_g(value: dict[str, Any]) -> None:
            event = value["datasets"][1]["events"][0]
            event["g_surface_close"]["status"] = "FAILED"
            event["g_surface_close"]["observed_at_utc"] = None
            event["g_surface_close"]["observed_monotonic_ns"] = None
            event["common_hold"] = None
            event["reads"] = []
            event["validity_status"] = "UNSCOREABLE"
        mutate_manifest(manifest_path, fail_one_g)
        formation = load_formation_input(manifest_path)
        failed = formation.manifest["datasets"][1]["events"][0]
        assert failed["g_surface_close"]["status"] == "FAILED"
        assert failed["validity_status"] == "UNSCOREABLE"
        assert formation.certificate()["scientific_verdict"] == "NONE_NOT_SCORED"
        checks += 1

        first_manifest, _ = write_bundle(base / "n31-shape-one", surface_shape="distributed mechanical marks")
        second_manifest, _ = write_bundle(base / "n31-shape-two", surface_shape="localized biochemical conformations")
        assert load_formation_input(first_manifest).certificate()["validation_dataset_eligible"]
        assert load_formation_input(second_manifest).certificate()["validation_dataset_eligible"]
        source_text = Path(contract.__file__).read_text(encoding="utf-8").lower()
        for forbidden in ("pt/tio2", "memristor", "redox", "0.1 v", "4.0 v", "86400 s"):
            assert forbidden not in source_text
        checks += 1

        manifest_path, execution_path = write_bundle(base / "n32-product-separate")
        certificate = attach_formation_execution(
            load_formation_input(manifest_path), execution_path
        ).certificate()
        assert certificate["validation_dataset_eligible"] is True
        assert certificate["product_reproduction"] == {
            "status": "NOT_DEMONSTRATED",
            "separate_from_dataset_validation": True,
            "attested": False,
        }
        assert certificate["program_completion_authorized"] is False
        checks += 1

    assert checks == 32, checks
    print(f"FORMATION_INPUT_CORE_CHECKS: {checks}/32 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
