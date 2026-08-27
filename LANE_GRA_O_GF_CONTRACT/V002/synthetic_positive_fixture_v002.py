#!/usr/bin/env python3
"""Deterministic synthetic positive fixture builder for validator V002.

The fixture exercises every mechanical join, including product reproduction and V004
transport/refinement.  It is explicitly synthetic and contains no actual platform;
therefore it can never produce authoritative GF0--UGE proof output.
"""

from __future__ import annotations

import base64
import copy
import hashlib
import json
from pathlib import Path
from typing import Any

import validator_v002 as v


WORKSPACE = Path(__file__).resolve().parents[2]
RSA_N = int("ae87b1f48f74ad7c09e1ee2649a7a5fc5fe524f6784c77f81556f5fefbbcbddec4c991c7276285e3e67a3d4c6feb88fc5a86475f812ab0b6ebe3e531ae9a6aef", 16)
RSA_D = int("9b34ef2b667bc35fc0461c0c0e8a395123525ab998a7a0348d7da50bd980b273b40aa998ef9530ba1ad1ea36fc9655ceee1f6d62fc171118055fa85d1a698c61", 16)
RSA_E = 65537


def power() -> dict[str, Any]:
    return {
        "alpha": 0.05, "beta": 0.2, "false_pass_upper": 0.01,
        "achieved_power_lower": 0.9, "monte_carlo_error": 0.001,
        "mc_tolerance": 0.01, "independent_units": 8, "effect_size_core": 1.0,
    }


def rsa_sign(message: bytes) -> str:
    length = (RSA_N.bit_length() + 7) // 8
    digest_info = bytes.fromhex("3031300d060960864801650304020105000420") + hashlib.sha256(message).digest()
    encoded = b"\x00\x01" + b"\xff" * (length - len(digest_info) - 3) + b"\x00" + digest_info
    return pow(int.from_bytes(encoded, "big"), RSA_D, RSA_N).to_bytes(length, "big").hex()


class ArtifactBuilder:
    def __init__(self) -> None:
        self.artifacts: list[dict[str, Any]] = []
        self.by_id: dict[str, dict[str, Any]] = {}

    def file(self, aid: str, kind: str, role: str, relative_path: str) -> str:
        raw = (WORKSPACE / relative_path).read_bytes()
        artifact = {
            "artifact_id": aid, "kind": kind, "role": role,
            "sha256": hashlib.sha256(raw).hexdigest(), "byte_length": len(raw),
            "storage": "WORKSPACE_FILE", "locator": relative_path, "payload_b64": None,
            "parents": [],
        }
        self.artifacts.append(artifact)
        self.by_id[aid] = artifact
        return aid

    def inline(self, aid: str, kind: str, role: str, payload: dict[str, Any]) -> str:
        refs = sorted(v.referenced_artifact_ids(payload))
        missing = [ref for ref in refs if ref not in self.by_id]
        if missing:
            raise ValueError(f"parents must exist before {aid}: {missing}")
        raw = v.canonical_bytes(payload)
        artifact = {
            "artifact_id": aid, "kind": kind, "role": role,
            "sha256": hashlib.sha256(raw).hexdigest(), "byte_length": len(raw),
            "storage": "INLINE_JSON", "locator": None,
            "payload_b64": base64.b64encode(raw).decode("ascii"),
            "parents": [{"artifact_id": ref, "sha256": self.by_id[ref]["sha256"]} for ref in refs],
        }
        self.artifacts.append(artifact)
        self.by_id[aid] = artifact
        return aid


def source(builder: ArtifactBuilder, aid: str, role: str, text: str) -> str:
    raw = text.encode()
    return builder.inline(aid, "SOURCE_DATA", role, {
        "source_id": aid, "content_b64": base64.b64encode(raw).decode(),
        "content_sha256": hashlib.sha256(raw).hexdigest(),
    })


def scale(scale_id: str, order: int, calibration_artifact_id: str) -> dict[str, Any]:
    ident = v.eye(2)
    zero = v.zeros(2, 2)
    return {
        "scale_id": scale_id, "order": order, "physical_scale": float(order + 1), "units": "calibrated_operator_norm",
        "W": ident, "eigenvectors": ident, "eigenvalues": [1.0, 1.0],
        "eigenvalue_intervals": [[0.9, 1.1], [0.9, 1.1]], "rank_threshold": 0.1,
        "max_condition": 10.0, "P": ident, "N": zero, "sqrt_W": ident,
        "pinv_sqrt_W": ident, "reference": [0.0, 0.0],
        "calibration_artifact_id": calibration_artifact_id,
    }


def edge(edge_id: str, source_id: str, target_id: str, causal_artifact_id: str) -> dict[str, Any]:
    ident = v.eye(2)
    return {
        "edge_id": edge_id, "source": source_id, "target": target_id,
        "F": ident, "Btilde": ident, "causal_superchannel_artifact_id": causal_artifact_id,
        "null_leakage_norm": 0.0, "transport_loss_norm": 0.0, "transport_rank": 2,
        "condition_number": 1.0, "null_margin": 0.01, "power": power(),
    }


def refinement(builder: ArtifactBuilder, control_type: str, index: int, scale_ids: list[str], edge_ids: list[str], base_streams: dict[str, str]) -> dict[str, Any]:
    root2 = 2.0 ** -0.5
    if control_type == "RELABEL":
        embedding = [[0.0, 1.0], [1.0, 0.0]]
    elif control_type == "SUBDIVIDE":
        embedding = [[root2, 0.0], [root2, 0.0], [0.0, 1.0]]
    else:
        embedding = [[root2, 0.0], [0.0, 1.0], [root2, 0.0]]
    aggregation = v.transpose(embedding)
    projector = v.matmul(embedding, v.transpose(embedding))
    null = v.matsub(v.eye(len(embedding)), projector)
    transform = source(builder, f"REF.{index}.TRANSFORM", "GENERIC", f"frozen {control_type} transform")
    paired = source(builder, f"REF.{index}.PAIRED", "GENERIC", f"same-event root {control_type}")
    scales: list[dict[str, Any]] = []
    for scale_id in scale_ids:
        refined_stream = source(builder, f"REF.{index}.{scale_id}.STREAM", "EXTERNAL_CALIBRATION", f"independently rerun {control_type} {scale_id}")
        scales.append({
            "scale_id": scale_id, "V": embedding, "R": aggregation,
            "W_refined": projector, "P_refined": projector, "N_refined": null,
            "U": embedding, "base_stream_artifact_id": base_streams[scale_id],
            "refined_stream_artifact_id": refined_stream,
        })
    edges = [{"edge_id": edge_id, "F_refined": projector, "Btilde_refined": projector} for edge_id in edge_ids]
    return {
        "control_id": f"REF-{index}-{control_type}", "type": control_type,
        "same_event": True, "frozen_before_response": True, "new_independent_mode": False,
        "transformation_artifact_id": transform, "paired_event_artifact_id": paired,
        "scales": scales, "edges": edges,
        "covariance": {"paired_events": True, "includes_cross_covariance": True, "pushforward_checked": True, "recompute_per_draw": True, "rank_boundary_crossed": False},
        "power": power(),
    }


def build_fixture() -> dict[str, Any]:
    b = ArtifactBuilder()
    framework_ids = []
    for index, path in enumerate(v.REQUIRED_FRAMEWORKS):
        framework_ids.append(b.file(f"FRAMEWORK.{index:02d}", "FRAMEWORK", "GENERIC", path))
    decision_id = b.file("PRINCIPAL.JOINT_SEED", "PRINCIPAL_DECISION", "THEORY", v.PRINCIPAL_DECISION_PATH)

    measurements = [source(b, f"MEASURE.{i}", "EXTERNAL_CALIBRATION", f"measured process stream {i}") for i in range(2)]
    process_ids: list[str] = []
    identification_ids: list[str] = []
    for index, carrier in enumerate(("MAGNETIC_LATTICE", "ATOMIC_CLOCK_NETWORK")):
        pid = b.inline(f"GAMMA.{index}", "PROCESS_REPRESENTATION", "EXTERNAL_CALIBRATION", {
            "surface_id": f"SURFACE-{index}", "carrier_class": carrier,
            "representation_type": "LIOUVILLIAN" if index == 0 else "STOCHASTIC_GENERATOR",
            "equation": "d rho / dt = Gamma_rec[rho]", "state_space": "finite measured process space",
            "native_units": ["s^-1"], "support": "measured record support", "bath": "measured bath channels",
            "clock": "traceable monotonic clock", "measured_channels": ["formation", "closure", "transition"],
        })
        process_ids.append(pid)
        identification_ids.append(b.inline(f"GAMMA.ID.{index}", "PROCESS_IDENTIFICATION", "EXTERNAL_CALIBRATION", {
            "surface_id": f"SURFACE-{index}", "process_artifact_id": pid,
            "measurement_artifact_ids": [measurements[index]], "method": "prospective process tomography",
            "projection_map_artifact_ids": [], "independent_process_identification": True,
        }))

    dilation_measure = source(b, "DILATION.MEASURE", "EXTERNAL_CALIBRATION", "synthetic full exchange stream")
    dilation_id = b.inline("DILATION.FULL", "FULL_DILATION", "EXTERNAL_CALIBRATION", {
        "surface_id": "SYNTHETIC-SURFACE", "system": "record writer", "environment": "bath plus support",
        "exchange_ledger": ["energy", "momentum", "work", "heat"],
        "physical_quantities": ["CURRENT", "STRESS", "FLUX"],
        "measurement_artifact_ids": [dilation_measure],
    })
    blocks = [
        {"family": "B", "operator_ids": ["B.transition"], "physical_type": "transition-support", "units": ["s^-1"], "nonredundancy_witness": "rank-B"},
        {"family": "C", "operator_ids": ["C.closure"], "physical_type": "closure-contrast", "units": ["1"], "nonredundancy_witness": "rank-C"},
        {"family": "D", "operator_ids": ["D.current"], "physical_type": "measured-current", "units": ["J s^-1"], "nonredundancy_witness": "rank-D"},
    ]
    seed_id = b.inline("SEED.JOINT", "SEED_DEFINITION", "THEORY", {
        "seed_id": "JOINT_SEED", "construction": "finite typed B union C union D dictionary",
        "blocks": blocks, "pre_response_inputs": ["Gamma_rec", "C", "C_EXT_OFF", "G", "K", "M", "R*", "B", "G_X"],
        "operator_type": "TYPED_JOINT_OPERATOR_SUBSPACE", "units": ["typed-per-block"],
        "support": "record-sufficient pre-response support", "roles": ["CURRENT", "FLUX"],
        "conserved_quantity": "ENERGY_MOMENTUM", "derivation_stage": "PRE_RESPONSE",
        "principal_decision_artifact_id": decision_id, "dilation_artifact_id": dilation_id,
    })

    scale_ids = ["L0", "L1", "L2"]
    edge_ids = ["E10", "E21", "E20"]
    calibration_ids = {sid: source(b, f"CAL.{sid}", "EXTERNAL_CALIBRATION", f"calibration {sid}") for sid in scale_ids}
    base_streams = {sid: source(b, f"BASE.{sid}.STREAM", "EXTERNAL_CALIBRATION", f"base stream {sid}") for sid in scale_ids}
    causal_ids = {eid: source(b, f"CAUSAL.{eid}", "EXTERNAL_CALIBRATION", f"causal superchannel {eid}") for eid in edge_ids}
    refinements = [refinement(b, kind, index, scale_ids, edge_ids, base_streams) for index, kind in enumerate(("RELABEL", "SUBDIVIDE", "DUPLICATE"), 1)]
    transport_id = b.inline("TRANSPORT.V004", "CALIBRATED_TRANSPORT", "DEVELOPMENT", {
        "transport_id": "SYNTHETIC-W-F-BTILDE-V004", "surface_id": "SYNTHETIC-SURFACE",
        "scales": [scale(sid, index, calibration_ids[sid]) for index, sid in enumerate(scale_ids)],
        "edges": [edge("E10", "L0", "L1", causal_ids["E10"]), edge("E21", "L1", "L2", causal_ids["E21"]), edge("E20", "L0", "L2", causal_ids["E20"])],
        "covariance": {"method": "JOINT_MONTE_CARLO", "uncertain_components": ["z", "F", "W", "reference", "rank", "shared_sources"], "shared_cross_scale": True, "recompute_sqrt_each_draw": True, "recompute_pinv_each_draw": True, "recompute_rank_each_draw": True, "rank_boundary_crossed": False, "common_support_stable": True},
        "composition": {"source": "L0", "middle": "L1", "target": "L2", "direct_edge": "E20", "first_edge": "E10", "second_edge": "E21", "P_in_common": v.eye(2), "P_out_common": v.eye(2), "residual_norm": 0.0, "residual_margin": 0.01, "power": power()},
        "refinements": refinements, "ancestry_component": "BTILDE",
        "freeze_time": "2026-08-22T10:00:00Z", "response_access_time": "2026-08-22T12:00:00Z",
    })

    d_source = source(b, "DATA.D.BYTES", "DEVELOPMENT", "development synthetic bytes")
    v_source = source(b, "DATA.V.BYTES", "VALIDATION", "validation synthetic bytes")
    d_root = b.inline("DATA.DEVELOPMENT", "DATASET_ROOT", "DEVELOPMENT", {
        "dataset_id": "D-SYNTH", "dataset_role": "DEVELOPMENT", "raw_root_sha256": hashlib.sha256(b"d-root").hexdigest(),
        "source_ids": ["D-source"], "acquisition_ids": ["D-acq"], "specimen_ids": ["D-specimen"],
        "independent_unit_ids": ["D-unit"], "event_ids": ["D-event"], "outcome_content_ids": ["D-outcome"],
        "source_artifact_ids": [d_source], "acquired_time": "2026-08-21T08:00:00Z",
    })
    v_root = b.inline("DATA.VALIDATION", "DATASET_ROOT", "VALIDATION", {
        "dataset_id": "V-SYNTH", "dataset_role": "VALIDATION", "raw_root_sha256": hashlib.sha256(b"v-root").hexdigest(),
        "source_ids": ["V-source"], "acquisition_ids": ["V-acq"], "specimen_ids": ["V-specimen"],
        "independent_unit_ids": ["V-unit"], "event_ids": ["V-event"], "outcome_content_ids": ["V-outcome"],
        "source_artifact_ids": [v_source], "acquired_time": "2026-08-21T09:00:00Z",
    })
    access_id = b.inline("DATA.ACCESS_LOG", "ACCESS_LOG", "GENERIC", {
        "freeze_time": "2026-08-22T10:00:00Z",
        "events": [
            {"dataset_id": "D-SYNTH", "event_type": "NUMERIC_ACCESS", "time": "2026-08-22T09:00:00Z"},
            {"dataset_id": "V-SYNTH", "event_type": "METADATA_ACCESS", "time": "2026-08-22T09:30:00Z"},
            {"dataset_id": "V-SYNTH", "event_type": "NUMERIC_ACCESS", "time": "2026-08-22T11:00:00Z"},
        ],
        "signed_by": "Brian", "signature_artifact_id": decision_id,
    })

    checks = [{"check_id": check, "observed": 1.0, "threshold": 0.9, "direction": "GE", "power": power()} for check in v.ANCESTRY_CHECKS]
    ancestry_id = b.inline("ANCESTRY.GRAPH", "ANCESTRY_GRAPH", "DEVELOPMENT", {
        "seed_node_id": "N0", "endpoint_node_id": "N2",
        "nodes": [
            {"node_id": f"N{i}", "surface_id": "SYNTHETIC-SURFACE", "scale_id": f"L{i}", "scale_order": i, "operator_type": "JOINT_SEED" if i == 0 else "BLOCKED_DESCENDANT", "transport_artifact_id": transport_id}
            for i in range(3)
        ],
        "arrows": [
            {"arrow_id": "A10", "source": "N0", "target": "N1", "transport_artifact_id": transport_id, "transport_component": "BTILDE", "checks": copy.deepcopy(checks)},
            {"arrow_id": "A21", "source": "N1", "target": "N2", "transport_artifact_id": transport_id, "transport_component": "BTILDE", "checks": copy.deepcopy(checks)},
        ],
    })
    joined_id = b.inline("EVIDENCE.JOINED", "JOINED_EVIDENCE", "VALIDATION", {
        "formation_artifact_ids": [d_source], "growth_artifact_ids": [d_root],
        "response_artifact_ids": [v_source], "join_keys": ["surface", "event", "scale", "clock"],
    })
    checker = source(b, "COVERAGE.CHECKER", "GENERIC", "synthetic constructive checker")
    coverage_id = b.inline("COVERAGE.THEOREM", "COVERAGE_THEOREM", "THEORY", {
        "theorem_id": "SYNTHETIC-CONSTRUCTIVE-COVERAGE", "surface_classes": ["MAGNETIC_LATTICE", "ATOMIC_CLOCK_NETWORK"],
        "constructive_map": "finite typed map from bona-fide surface descriptor to package",
        "preservation_obligations": ["Gamma", "seed", "transport", "falsifier"], "checker_artifact_ids": [checker],
    })

    gate_inputs: list[dict[str, Any]] = []
    for index, gate in enumerate(v.GATES):
        predicate = f"PREDICATE.{index:02d}"
        rule = b.inline(f"RULE.{index:02d}", "GATE_RULE", "THEORY", {
            "gate_id": gate, "predicate_ids": [predicate], "decision_region": "prospective finite equivalence region",
            "freeze_time": "2026-08-22T10:00:00Z", "rule_version": "V002",
        })
        evidence = b.inline(f"EVIDENCE.{index:02d}", "GATE_EVIDENCE", "VALIDATION", {
            "gate_id": gate, "observations": [{"predicate_id": predicate, "value": True, "reproducible": True}],
            "source_artifact_ids": [v_source], "taxonomy_match": "KNOWN", "reproducible": True,
        })
        pcert = b.inline(f"POWER.{index:02d}", "POWER_CERTIFICATE", "THEORY", {"gate_id": gate, **power()})
        gate_inputs.append({"gate_id": gate, "rule_artifact_id": rule, "evidence_artifact_id": evidence, "power_artifact_id": pcert, "prerequisite_gate_ids": list(v.PREREQUISITES[gate])})

    domain_rows: list[dict[str, Any]] = []
    for claim in v.CLAIMS:
        domain_id = b.inline(f"DOMAIN.{claim}", "DOMAIN_WITNESS", "THEORY", {
            "claim_id": claim, "member_ids": ["SYNTHETIC-SURFACE-0", "SYNTHETIC-SURFACE-1"],
            "quantifier": "FOR_ALL_DECLARED_MEMBERS", "coverage_artifact_id": coverage_id,
        })
        domain_rows.append({"claim_id": claim, "domain_artifact_id": domain_id})

    public_manifest = b.inline("PRODUCT.DATA_MANIFEST", "PUBLIC_DATASET_MANIFEST", "PRODUCT", {
        "release_id": "SYNTHETIC-RELEASE", "dataset_artifact_ids": [d_root, v_root],
        "dataset_sha256s": {d_root: b.by_id[d_root]["sha256"], v_root: b.by_id[v_root]["sha256"]},
        "public_locators": ["https://invalid.example/synthetic-d", "https://invalid.example/synthetic-v"],
        "license_ids": ["SYNTHETIC-ONLY"],
    })
    execution_output = b.inline("PRODUCT.OUTPUT", "EXECUTION_OUTPUT", "PRODUCT", {
        "release_id": "SYNTHETIC-RELEASE", "result_digest": hashlib.sha256(b"no-proof-output").hexdigest(),
        "claim_results": [[claim, "NO_PROOF_OUTPUT"] for claim in v.CLAIMS],
    })
    release = b.inline("PRODUCT.RELEASE", "PUBLIC_RELEASE", "PRODUCT", {
        "release_id": "SYNTHETIC-RELEASE", "contract_id": v.CONTRACT_ID, "package_id": "SYNTHETIC.POSITIVE.V002",
        "claim_ids": list(v.CLAIMS), "dataset_manifest_artifact_id": public_manifest,
        "dataset_manifest_sha256": b.by_id[public_manifest]["sha256"],
        "validator_sha256": hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest(),
        "expected_output_artifact_id": execution_output, "expected_output_sha256": b.by_id[execution_output]["sha256"],
    })
    key_id = b.inline("PRODUCT.EXECUTOR_KEY", "EXECUTOR_KEY", "PRODUCT", {
        "executor_id": "INDEPENDENT-SYNTHETIC-EXECUTOR", "n_hex": f"{RSA_N:x}", "e": RSA_E,
    })
    report = b.inline("PRODUCT.REPORT", "EXECUTION_REPORT", "PRODUCT", {
        "release_artifact_id": release, "output_artifact_id": execution_output,
        "executor_id": "INDEPENDENT-SYNTHETIC-EXECUTOR", "environment_sha256": hashlib.sha256(b"synthetic-env").hexdigest(),
        "no_private_logic": True, "mismatch_count": 0,
    })
    message = v.canonical_bytes({
        "contract_id": v.CONTRACT_ID, "package_id": "SYNTHETIC.POSITIVE.V002",
        "release_sha256": b.by_id[release]["sha256"], "output_sha256": b.by_id[execution_output]["sha256"],
        "report_sha256": b.by_id[report]["sha256"],
    })
    signature = b.inline("PRODUCT.SIGNATURE", "EXECUTION_SIGNATURE", "PRODUCT", {
        "key_artifact_id": key_id, "release_artifact_id": release, "output_artifact_id": execution_output,
        "report_artifact_id": report, "message_b64": base64.b64encode(message).decode(), "signature_hex": rsa_sign(message),
    })

    bindings = {
        "framework_artifact_ids": framework_ids, "principal_decision_artifact_id": decision_id,
        "gamma_process_artifact_ids": process_ids, "gamma_identification_artifact_ids": identification_ids,
        "seed_definition_artifact_id": seed_id, "dilation_artifact_id": dilation_id,
        "platform_instantiation_artifact_id": None, "transport_artifact_id": transport_id,
        "development_root_artifact_id": d_root, "validation_root_artifact_id": v_root,
        "access_log_artifact_id": access_id, "ancestry_artifact_id": ancestry_id,
        "joined_evidence_artifact_id": joined_id, "coverage_artifact_id": coverage_id,
        "public_release_artifact_id": release, "public_dataset_manifest_artifact_id": public_manifest,
        "executor_key_artifact_id": key_id, "execution_output_artifact_id": execution_output,
        "execution_report_artifact_id": report, "execution_signature_artifact_id": signature,
    }
    return {
        "contract_id": v.CONTRACT_ID, "task_id": v.TASK_ID, "package_id": "SYNTHETIC.POSITIVE.V002",
        "package_mode": "SYNTHETIC_TEST", "producer_id": "V002-FIXTURE-BUILDER", "principal_id": "Brian",
        "claim_ids": list(v.CLAIMS), "gate_ids": list(v.GATES), "horizon_complete_claim": False,
        "bindings": bindings, "artifacts": b.artifacts, "gate_inputs": gate_inputs,
        "gravity_applicability": [{"characteristic_id": gate, "applicability": "APPLICABLE", "justification_artifact_id": None} for gate in v.GRAVITY_GATES],
        "milestone_domains": domain_rows,
    }


def main() -> None:
    print(json.dumps(build_fixture(), indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
