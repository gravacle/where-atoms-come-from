#!/usr/bin/env python3
"""Forty-two fixed checks for the origin-neutral gamma-flow core.

All fixtures are synthetic contract tests.  They are not scientific evidence, and
the positive fixture must retain NO_PROOF_OUTPUT for GF0 through UGE.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
import tempfile
from typing import Any

import gamma_flow as gf
from gamma_flow import GammaFlowRefusal, load_gamma_flow

import mutation_suite_v002 as sealed_suite
import mutation_suite_v002_repair2 as repair2_suite
import synthetic_positive_fixture_v002 as sealed_fixture


def extension_artifact(artifact_id: str, kind: str, payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "artifact_id": artifact_id,
        "kind": kind,
        "sha256": gf.sha256_bytes(gf.canonical_bytes(payload)),
        "payload": payload,
    }


def rehash(item: dict[str, Any]) -> None:
    item["sha256"] = gf.sha256_bytes(gf.canonical_bytes(item["payload"]))


def artifact_by_id(manifest: dict[str, Any], artifact_id: str) -> dict[str, Any]:
    return next(
        item for item in manifest["extension_artifacts"]
        if item["artifact_id"] == artifact_id
    )


def build_manifest() -> dict[str, Any]:
    candidate = sealed_fixture.build_fixture()
    registry, payloads, errors, missing = gf._gf_base.validate_artifacts(candidate, gf.ROOT)
    assert not errors and not missing
    bindings = candidate["bindings"]
    transport_id = bindings["transport_artifact_id"]
    transport = payloads[transport_id]
    ancestry = payloads[bindings["ancestry_artifact_id"]]
    seed = payloads[bindings["seed_definition_artifact_id"]]
    seed_operators = sorted(
        operator for block in seed["blocks"] for operator in block["operator_ids"]
    )

    extension: list[dict[str, Any]] = []
    scale_semantics: list[dict[str, Any]] = []
    response_kernels: list[dict[str, Any]] = []
    material_id = "GF.METROLOGY.MATERIAL"
    control_id = "GF.CONTROL.KNOWN_FORCE"
    material_operators = [
        gf.MATERIAL_COVERAGE_OPERATOR_IDS[category]
        for category in gf.MATERIAL_COVERAGE_CATEGORIES
    ]
    material_sources = [
        bindings["dilation_artifact_id"],
        *payloads[bindings["dilation_artifact_id"]]["measurement_artifact_ids"],
    ]
    material_units = {
        "MASS_ENERGY": "J",
        "MOMENTUM_CURRENT": "N",
        "STRESS_PRESSURE": "Pa",
        "MATERIAL_COMPOSITION": "mol/m^3",
        "GEOMETRY_SUPPORT": "m",
        "WRITER_BATH_PROBE_EXCHANGE": "W",
        "UNCERTAINTY_COVARIANCE": "J^2",
    }
    coverage_map = {
        category: {
            "operator_id": operator_id,
            "evidence_status": (
                "BOUNDED" if category == "UNCERTAINTY_COVARIANCE" else "MEASURED"
            ),
            "quantity_type": category,
            "si_unit": material_units[category],
            "source_artifact_ids": list(material_sources),
        }
        for category, operator_id in zip(
            gf.MATERIAL_COVERAGE_CATEGORIES, material_operators
        )
    }
    control_operator = "known.force.null"
    extension.append(extension_artifact(material_id, "MATERIAL_METROLOGY", {
        "operator_ids": list(material_operators),
        "scope": gf.MATERIAL_METROLOGY_SCOPE,
        "quantity_types": list(gf.MATERIAL_COVERAGE_CATEGORIES),
        "unit_system": gf.MATERIAL_METROLOGY_UNIT_SYSTEM,
        "source_artifact_ids": list(material_sources),
        "coverage_map": coverage_map,
    }))
    extension.append(extension_artifact(control_id, "KNOWN_FORCE_CONTROL", {
        "force_family": "registered ordinary-force null basis",
        "null_operator_ids": [control_operator],
        "source_artifact_ids": [transport_id],
    }))

    custody_categories = [
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
    ]
    lane: dict[str, Any] = {
        "roles": sorted(gf.METROLOGY_ROLES),
        "material_metrology_artifact_ids": [material_id],
        "known_force_control_artifact_ids": [control_id],
    }
    for index, category in enumerate(custody_categories):
        artifact_id = f"GF.CUSTODY.{index:02d}"
        lane[category] = [artifact_id]
        extension.append(extension_artifact(artifact_id, "METROLOGY_CUSTODY", {
            "custody_role": category,
            "source_artifact_ids": [transport_id],
        }))

    complete_basis = sorted(seed_operators + material_operators + [control_operator])
    output_by_scale: dict[str, str] = {}
    for index, scale in enumerate(transport["scales"]):
        scale_id = scale["scale_id"]
        definition_id = f"GF.SCALE.{scale_id}.DEFINITION"
        probe_id = f"probe.collective.{scale_id}"
        output_by_scale[scale_id] = probe_id
        kernel_id = f"GF.RESPONSE.{scale_id}"
        detector_id = f"GF.DETECTOR.{scale_id}"
        extension.append(extension_artifact(definition_id, "SCALE_DEFINITION", {
            "scale_id": scale_id,
            "scale_kind": "POPULATION",
            "units": "generation",
            "support_identity": f"blocked support at {scale_id}",
            "aggregation_rule": f"prospectively frozen population block {index + 1}",
            "monotone_order_justification": "larger block order contains no response outcome",
        }))
        extension.append(extension_artifact(kernel_id, "RESPONSE_KERNEL", {
            "surface_id": transport["surface_id"],
            "scale_id": scale_id,
            "input_operator_ids": list(complete_basis),
            "output_probe_ids": [probe_id],
            "kernel_representation": "generic calibrated linear-response container",
            "source_artifact_ids": [transport_id],
        }))
        extension.append(extension_artifact(detector_id, "DETECTOR_MAP", {
            "surface_id": transport["surface_id"],
            "output_probe_ids": [probe_id],
            "detector_channels": [f"synthetic-channel-{index}"],
            "calibration_artifact_ids": [scale["calibration_artifact_id"]],
        }))
        scale_semantics.append({
            "scale_id": scale_id,
            "scale_kind": "POPULATION",
            "scale_definition_artifact_id": definition_id,
        })
        response_kernels.append({
            "surface_id": transport["surface_id"],
            "scale_id": scale_id,
            "input_operator_ids": list(complete_basis),
            "output_probe_ids": [probe_id],
            "input_policy": "JOINT_FROZEN_BASIS",
            "response_kernel_artifact_id": kernel_id,
            "detector_map_artifact_id": detector_id,
            "material_metrology_artifact_ids": [material_id],
            "known_force_control_artifact_ids": [control_id],
            "freeze_time": transport["freeze_time"],
            "response_access_time": transport["response_access_time"],
        })

    stage_names = ("SEED", "COLLECTIVE", "RELATIONAL")
    ancestry_stages: list[dict[str, Any]] = []
    for index, node in enumerate(ancestry["nodes"]):
        ancestry_stages.append({
            "node_id": node["node_id"],
            "emergence_stage": stage_names[index],
            "representation_class": (
                "NON_GEOMETRIC" if index == 0 else "UNCLASSIFIED"
            ),
            "operator_ids": (
                seed_operators if index == 0 else [output_by_scale[node["scale_id"]]]
            ),
        })

    return {
        "schema": gf.SCHEMA,
        "flow_id": "GAMMA.FLOW.SYNTHETIC.V001",
        "package_mode": "SYNTHETIC_TEST",
        "theory_bindings": {
            "principal_decision": {
                "path": gf.PRINCIPAL_DECISION_PATH,
                "sha256": gf.PRINCIPAL_DECISION_SHA256,
            },
            "principal_clarification": {
                "path": gf.PRINCIPAL_CLARIFICATION_PATH,
                "sha256": gf.PRINCIPAL_CLARIFICATION_SHA256,
            },
            "repair3_manifest": {
                "path": gf.REPAIR3_MANIFEST_PATH,
                "sha256": gf.REPAIR3_MANIFEST_SHA256,
            },
            "gravity_characteristic_registry": {
                "path": gf.CHARACTERISTIC_REGISTRY_PATH,
                "sha256": gf.CHARACTERISTIC_REGISTRY_SHA256,
            },
        },
        "contract_candidate": candidate,
        "extension_artifacts": extension,
        "scale_semantics": scale_semantics,
        "ancestry_stages": ancestry_stages,
        "response_kernels": response_kernels,
        "metric_reconstructions": [],
        "metrology_lane": lane,
    }


def write_manifest(root: Path, manifest: dict[str, Any], name: str) -> Path:
    path = root / name
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return path


def evaluate(root: Path, manifest: dict[str, Any], name: str) -> Any:
    return load_gamma_flow(write_manifest(root, manifest, name))


def expect_refusal(root: Path, manifest: dict[str, Any], name: str, text: str) -> None:
    try:
        evaluate(root, manifest, name)
    except GammaFlowRefusal as exc:
        assert text in str(exc), str(exc)
    else:
        raise AssertionError(f"gamma flow accepted mutation {name}")


def add_metric_candidate(manifest: dict[str, Any]) -> None:
    manifest["ancestry_stages"][-1]["emergence_stage"] = "METRIC_CANDIDATE"
    probes = sorted(gf.METRIC_PROBE_FAMILIES)
    transport_id = manifest["contract_candidate"]["bindings"]["transport_artifact_id"]
    definitions = [
        (
            "GF.METRIC.MULTIPROBE",
            "METRIC_RECONSTRUCTION",
            {
                "probe_families": probes,
                "relational_observables": ["clock", "matter", "light", "independent probe"],
                "source_artifact_ids": [transport_id],
            },
        ),
        (
            "GF.CONNECTION.MULTIPROBE",
            "CONNECTION_RECONSTRUCTION",
            {
                "probe_families": probes,
                "transport_observables": ["parallel transport", "trajectory"],
                "source_artifact_ids": [transport_id],
            },
        ),
        (
            "GF.CONE.MULTIPROBE",
            "COMMON_CONE_TEST",
            {
                "probe_families": probes,
                "propagation_observables": ["time of flight", "dispersion"],
                "source_artifact_ids": [transport_id],
            },
        ),
        (
            "GF.METRIC.ALTERNATIVES",
            "POWERED_ALTERNATIVES",
            {
                "alternative_ids": ["acoustic-fit", "carrier-geometry", "graph-distance"],
                "rule_artifact_ids": [transport_id],
                "power_artifact_ids": [transport_id],
            },
        ),
    ]
    manifest["extension_artifacts"].extend(
        extension_artifact(artifact_id, kind, payload)
        for artifact_id, kind, payload in definitions
    )
    manifest["metric_reconstructions"] = [{
        "node_id": manifest["ancestry_stages"][-1]["node_id"],
        "probe_families": probes,
        "metric_artifact_id": "GF.METRIC.MULTIPROBE",
        "connection_artifact_id": "GF.CONNECTION.MULTIPROBE",
        "common_cone_artifact_id": "GF.CONE.MULTIPROBE",
        "powered_alternatives_artifact_id": "GF.METRIC.ALTERNATIVES",
    }]


def material_artifact(manifest: dict[str, Any]) -> dict[str, Any]:
    return next(
        artifact
        for artifact in manifest["extension_artifacts"]
        if artifact["kind"] == "MATERIAL_METROLOGY"
    )


def replace_material_basis(
    manifest: dict[str, Any], old_operators: list[str], new_operators: list[str]
) -> None:
    old = set(old_operators)
    replacement = set(new_operators)
    for row in manifest["response_kernels"]:
        row["input_operator_ids"] = sorted(
            (set(row["input_operator_ids"]) - old) | replacement
        )
        kernel = artifact_by_id(manifest, row["response_kernel_artifact_id"])
        kernel["payload"]["input_operator_ids"] = list(row["input_operator_ids"])
        rehash(kernel)


def actual_ir_manifest() -> dict[str, Any]:
    """Lift the sealed Repair2 typed platform into the gamma overlay."""
    manifest = build_manifest()
    manifest["contract_candidate"] = repair2_suite.typed_platform_candidate()
    manifest["package_mode"] = "SCIENTIFIC"
    for row in manifest["response_kernels"]:
        row["surface_id"] = "SURFACE-0"
    for artifact in manifest["extension_artifacts"]:
        if artifact["kind"] in {"RESPONSE_KERNEL", "DETECTOR_MAP"}:
            artifact["payload"]["surface_id"] = "SURFACE-0"
            rehash(artifact)
    add_metric_candidate(manifest)
    manifest["ancestry_stages"][-1]["emergence_stage"] = "IR_GRAVITY"
    return manifest


def main() -> int:
    checks = 0
    with tempfile.TemporaryDirectory(prefix="wac-gamma-flow-") as directory:
        root = Path(directory)
        positive_manifest = build_manifest()
        positive = evaluate(root, positive_manifest, "positive.json")
        certificate = positive.certificate()

        assert certificate["accepted"] is True
        checks += 1
        assert certificate["repair3_custody"] == "QUALIFIED"
        checks += 1
        assert certificate["actual_platform_present"] is False
        checks += 1
        assert certificate["contract_candidate_identity_preserved"] is True
        checks += 1
        assert set(certificate["authoritative_proof_outputs"].values()) == {"NO_PROOF_OUTPUT"}
        checks += 1
        assert not any(
            certificate[key] for key in (
                "scientific_validation_authorized",
                "record_formation_claim_authorized",
                "scoped_gamma_flow_claim_authorized",
                "gravity_claim_authorized",
                "universal_claim_authorized",
                "program_completion_authorized",
            )
        )
        checks += 1
        assert certificate["metric_required_at_GF0_or_GF1"] is False
        assert certificate["geometry_required_at_GF0_or_GF1"] is False
        assert certificate["stress_energy_required_at_GF0_or_GF1"] is False
        assert certificate["t_lab_exclusive_microscopic_source_authorized"] is False
        checks += 1
        states = certificate["internal_discovery_states"]
        assert states == {
            "SEED_OPERATOR_ANCESTRY": "UNSCOREABLE",
            "GENERIC_COLLECTIVE_RESPONSE": "UNSCOREABLE",
            "OPTIONAL_METRIC_RECONSTRUCTION": "UNCLASSIFIED",
            "CLASSICAL_GRAVITY_ENDPOINT": "UNCLASSIFIED",
        }
        checks += 1
        assert certificate["stage_reasons"]["OPTIONAL_METRIC_RECONSTRUCTION"] == [
            "NO_METRIC_STAGE_DECLARED_EARLY_METRIC_NOT_REQUIRED"
        ]
        checks += 1
        assert {row["scale_kind"] for row in positive_manifest["scale_semantics"]} == {"POPULATION"}
        assert not positive.extension_errors
        checks += 1
        assert set(positive_manifest["metrology_lane"]["roles"]) == gf.METROLOGY_ROLES
        positive_material = material_artifact(positive_manifest)["payload"]
        assert positive_material["scope"] == gf.MATERIAL_METROLOGY_SCOPE
        assert positive_material["unit_system"] == gf.MATERIAL_METROLOGY_UNIT_SYSTEM
        assert list(positive_material["coverage_map"]) == list(
            gf.MATERIAL_COVERAGE_CATEGORIES
        )
        assert positive_material["quantity_types"] == list(
            gf.MATERIAL_COVERAGE_CATEGORIES
        )
        assert {
            row["evidence_status"]
            for row in positive_material["coverage_map"].values()
        } == {"MEASURED", "BOUNDED"}
        checks += 1
        assert certificate["gravity_characteristic_registry_sha256"] == gf.CHARACTERISTIC_REGISTRY_SHA256
        checks += 1
        assert certificate["repair3_manifest_sha256"] == gf.REPAIR3_MANIFEST_SHA256
        assert certificate["principal_clarification_sha256"] == gf.PRINCIPAL_CLARIFICATION_SHA256
        checks += 1

        temporal = build_manifest()
        for row in temporal["scale_semantics"]:
            row["scale_kind"] = "TEMPORAL"
            item = artifact_by_id(temporal, row["scale_definition_artifact_id"])
            item["payload"]["scale_kind"] = "TEMPORAL"
            item["payload"]["units"] = "s"
            rehash(item)
        assert evaluate(root, temporal, "temporal.json").certificate()["accepted"] is True
        checks += 1

        unclassified_seed = build_manifest()
        unclassified_seed["ancestry_stages"][0]["representation_class"] = "UNCLASSIFIED"
        unclassified_result = evaluate(root, unclassified_seed, "unclassified-seed.json")
        assert unclassified_result.certificate()["accepted"] is True
        assert not any("METRIC" in error for error in unclassified_result.extension_errors)
        checks += 1

        metric_like = build_manifest()
        add_metric_candidate(metric_like)
        metric_result = evaluate(root, metric_like, "metric-like.json").certificate()
        assert metric_result["accepted"] is True
        assert metric_result["internal_discovery_states"]["OPTIONAL_METRIC_RECONSTRUCTION"] == "UNSCOREABLE"
        assert metric_result["internal_discovery_states"]["CLASSICAL_GRAVITY_ENDPOINT"] == "UNCLASSIFIED"
        assert metric_result["authoritative_proof_outputs"]["GE2"] == "NO_PROOF_OUTPUT"
        checks += 1

        t_only = build_manifest()
        t_only["response_kernels"][0]["input_policy"] = "T_LAB_ONLY"
        t_result = evaluate(root, t_only, "t-only.json").certificate()
        assert t_result["accepted"] is False
        assert any("EXCLUSIVE_SOURCE_POLICY_FORBIDDEN" in item for item in t_result["extension_errors"])
        checks += 1

        probability = build_manifest()
        seed_id = probability["contract_candidate"]["bindings"]["seed_definition_artifact_id"]
        seed = sealed_suite.payload(probability["contract_candidate"], seed_id)
        next(block for block in seed["blocks"] if block["family"] == "D")["physical_type"] = "probability-current"
        sealed_suite.replace_payload(probability["contract_candidate"], seed_id, seed)
        p_result = evaluate(root, probability, "probability-current.json").certificate()
        assert p_result["accepted"] is False
        assert "PROBABILITY_CURRENT_IS_NOT_PHYSICAL_D" in p_result["extension_errors"]
        checks += 1

        incomplete_metrology = build_manifest()
        incomplete_metrology["metrology_lane"]["known_force_control_artifact_ids"] = []
        expect_refusal(root, incomplete_metrology, "incomplete-metrology.json", "unique string list")
        checks += 1

        ir_label = build_manifest()
        ir_label["ancestry_stages"][-1]["emergence_stage"] = "IR_GRAVITY"
        ir_result = evaluate(root, ir_label, "ir-label.json").certificate()
        assert ir_result["accepted"] is False
        assert "IR_GRAVITY_LABEL_WITHOUT_CHARACTERISTIC_CONJUNCTION" in ir_result["extension_errors"]
        checks += 1

        bad = build_manifest()
        bad["theory_bindings"]["principal_clarification"]["sha256"] = "0" * 64
        expect_refusal(root, bad, "bad-clarification.json", "principal_clarification")
        checks += 1

        bad = build_manifest()
        bad["theory_bindings"]["repair3_manifest"]["sha256"] = "0" * 64
        expect_refusal(root, bad, "bad-repair3.json", "repair3_manifest")
        checks += 1

        bad = build_manifest()
        bad["theory_bindings"]["principal_decision"]["sha256"] = "0" * 64
        expect_refusal(root, bad, "bad-decision.json", "principal_decision")
        checks += 1

        bad = build_manifest()
        bad["theory_bindings"]["gravity_characteristic_registry"]["sha256"] = "0" * 64
        expect_refusal(root, bad, "bad-characteristics.json", "gravity_characteristic_registry")
        checks += 1

        bad = build_manifest()
        bad["injected_outcome"] = "PASS"
        expect_refusal(root, bad, "extra-root.json", "key closure")
        checks += 1

        duplicate_path = root / "duplicate.json"
        duplicate_path.write_text(
            '{"schema":"WAC_GAMMA_FLOW_V001","schema":"WAC_GAMMA_FLOW_V001"}',
            encoding="utf-8",
        )
        try:
            load_gamma_flow(duplicate_path)
        except GammaFlowRefusal as exc:
            assert "duplicate JSON member" in str(exc)
        else:
            raise AssertionError("duplicate JSON member escaped")
        checks += 1

        nonfinite_path = root / "nonfinite.json"
        nonfinite_path.write_text('{"schema":NaN}', encoding="utf-8")
        try:
            load_gamma_flow(nonfinite_path)
        except GammaFlowRefusal as exc:
            assert "nonfinite JSON constant" in str(exc)
        else:
            raise AssertionError("nonfinite JSON escaped")
        checks += 1

        bad = build_manifest()
        bad["extension_artifacts"][0]["sha256"] = "0" * 64
        expect_refusal(root, bad, "extension-hash.json", "content hash mismatch")
        checks += 1

        bad = build_manifest()
        bad["extension_artifacts"].append(copy.deepcopy(bad["extension_artifacts"][0]))
        expect_refusal(root, bad, "duplicate-extension-id.json", "duplicate extension artifact ID")
        checks += 1

        bad = build_manifest()
        duplicate = copy.deepcopy(bad["extension_artifacts"][0])
        duplicate["artifact_id"] = "GF.DUPLICATE.DIGEST"
        bad["extension_artifacts"].append(duplicate)
        expect_refusal(root, bad, "duplicate-extension-digest.json", "reuses an extension payload digest")
        checks += 1

        unresolved = build_manifest()
        custody = next(item for item in unresolved["extension_artifacts"] if item["kind"] == "METROLOGY_CUSTODY")
        custody["payload"]["source_artifact_ids"] = ["UNKNOWN.ARTIFACT"]
        rehash(custody)
        unresolved_result = evaluate(root, unresolved, "unresolved-source.json").certificate()
        assert unresolved_result["accepted"] is False
        assert any("EXTENSION_SOURCE_UNRESOLVED" in item for item in unresolved_result["extension_errors"])
        checks += 1

        missing_scale = build_manifest()
        missing_scale["scale_semantics"].pop()
        missing_result = evaluate(root, missing_scale, "missing-scale.json").certificate()
        assert missing_result["accepted"] is False
        assert "SCALE_SEMANTICS_NOT_EXACTLY_TRANSPORT_SCALES" in missing_result["extension_errors"]
        checks += 1

        bad = build_manifest()
        bad["scale_semantics"][0]["scale_kind"] = "METRIC_LENGTH_ONLY"
        expect_refusal(root, bad, "unknown-scale-kind.json", "scale_kind is not registered")
        checks += 1

        incomplete_basis = build_manifest()
        row = incomplete_basis["response_kernels"][0]
        row["input_operator_ids"].remove("B.transition")
        kernel = artifact_by_id(incomplete_basis, row["response_kernel_artifact_id"])
        kernel["payload"]["input_operator_ids"].remove("B.transition")
        rehash(kernel)
        basis_result = evaluate(root, incomplete_basis, "incomplete-basis.json").certificate()
        assert basis_result["accepted"] is False
        assert any("RESPONSE_JOINT_BASIS_INCOMPLETE" in item for item in basis_result["extension_errors"])
        checks += 1

        late = build_manifest()
        late["response_kernels"][0]["response_access_time"] = late["response_kernels"][0]["freeze_time"]
        late_result = evaluate(root, late, "late-freeze.json").certificate()
        assert late_result["accepted"] is False
        assert any("RESPONSE_NOT_FROZEN_BEFORE_ACCESS" in item for item in late_result["extension_errors"])
        checks += 1

        regressive = build_manifest()
        regressive["ancestry_stages"][-1]["emergence_stage"] = "SEED"
        regression_result = evaluate(root, regressive, "stage-regression.json").certificate()
        assert regression_result["accepted"] is False
        assert any("ANCESTRY_STAGE_REGRESSION" in item for item in regression_result["extension_errors"])
        checks += 1

        qa_record_only = actual_ir_manifest()
        qa_baseline = evaluate(root, qa_record_only, "qa-actual-baseline.json").certificate()
        assert qa_baseline["accepted"] is True
        assert qa_baseline["actual_platform_present"] is True
        assert qa_baseline["authoritative_proof_outputs"]["GE2"] == "PASSES_DECLARED_DOMAIN"
        assert qa_baseline["gravity_claim_authorized"] is True
        qa_material = material_artifact(qa_record_only)
        old_material_operators = list(qa_material["payload"]["operator_ids"])
        qa_material["payload"].update({
            "operator_ids": ["record.count.only"],
            "scope": "record variables only",
            "quantity_types": ["record_count"],
            "unit_system": "dimensionless",
        })
        rehash(qa_material)
        replace_material_basis(
            qa_record_only, old_material_operators, ["record.count.only"]
        )
        qa_result = evaluate(root, qa_record_only, "qa-record-only.json").certificate()
        assert qa_result["accepted"] is False
        assert qa_result["gravity_claim_authorized"] is False
        assert set(qa_result["authoritative_proof_outputs"].values()) == {
            "NO_PROOF_OUTPUT"
        }
        assert any(
            "MATERIAL_METROLOGY_SCOPE_UNREGISTERED" in error
            for error in qa_result["extension_errors"]
        )
        checks += 1

        missing_category = build_manifest()
        missing_material = material_artifact(missing_category)
        del missing_material["payload"]["coverage_map"]["UNCERTAINTY_COVARIANCE"]
        rehash(missing_material)
        expect_refusal(
            root,
            missing_category,
            "material-category-missing.json",
            "coverage_map key closure failed",
        )
        checks += 1

        dimensionless = build_manifest()
        dimensionless_material = material_artifact(dimensionless)
        dimensionless_material["payload"]["coverage_map"]["MASS_ENERGY"][
            "si_unit"
        ] = "dimensionless"
        rehash(dimensionless_material)
        dimensionless_result = evaluate(
            root, dimensionless, "material-dimensionless.json"
        ).certificate()
        assert dimensionless_result["accepted"] is False
        assert any(
            "MATERIAL_METROLOGY_SI_UNIT_INVALID" in error
            for error in dimensionless_result["extension_errors"]
        )
        checks += 1

        substituted_source = build_manifest()
        substituted_material = material_artifact(substituted_source)
        substituted_material["payload"]["coverage_map"]["MASS_ENERGY"][
            "source_artifact_ids"
        ].append(substituted_material["artifact_id"])
        substituted_material["payload"]["source_artifact_ids"].append(
            substituted_material["artifact_id"]
        )
        rehash(substituted_material)
        substituted_result = evaluate(
            root, substituted_source, "material-extension-self-source.json"
        ).certificate()
        assert substituted_result["accepted"] is False
        assert any(
            "MATERIAL_METROLOGY_SOURCE_OUTSIDE_CANDIDATE_CUSTODY" in error
            for error in substituted_result["extension_errors"]
        )
        checks += 1

        operator_mismatch = build_manifest()
        mismatch_material = material_artifact(operator_mismatch)
        old_material_operators = list(mismatch_material["payload"]["operator_ids"])
        mismatch_material["payload"]["operator_ids"][0] = (
            "material.control.forged_mass_energy"
        )
        rehash(mismatch_material)
        replace_material_basis(
            operator_mismatch,
            old_material_operators,
            mismatch_material["payload"]["operator_ids"],
        )
        mismatch_result = evaluate(
            root, operator_mismatch, "material-operator-coverage-mismatch.json"
        ).certificate()
        assert mismatch_result["accepted"] is False
        assert any(
            "MATERIAL_METROLOGY_OPERATOR_COVERAGE_MISMATCH" in error
            for error in mismatch_result["extension_errors"]
        )
        checks += 1

        absent_dilation = build_manifest()
        absent_material = material_artifact(absent_dilation)
        dilation_id = absent_dilation["contract_candidate"]["bindings"][
            "dilation_artifact_id"
        ]
        absent_material["payload"]["source_artifact_ids"].remove(dilation_id)
        for coverage_row in absent_material["payload"]["coverage_map"].values():
            coverage_row["source_artifact_ids"].remove(dilation_id)
        rehash(absent_material)
        absent_result = evaluate(
            root, absent_dilation, "material-dilation-join-absent.json"
        ).certificate()
        assert absent_result["accepted"] is False
        assert any(
            "MATERIAL_METROLOGY_DILATION_JOIN_MISSING" in error
            for error in absent_result["extension_errors"]
        )
        checks += 1

    assert checks == 42, checks
    print(f"GAMMA_FLOW_CORE_CHECKS: {checks}/42 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
