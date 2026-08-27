#!/usr/bin/env python3
"""Fixed V002 regression suite: positive, V001 witness/bypasses, V004 guards."""

from __future__ import annotations

import argparse
import base64
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any, Callable

import synthetic_positive_fixture_v002 as fixture
import validator_v002 as v


WORKSPACE = Path(__file__).resolve().parents[2]


def artifact(candidate: dict[str, Any], aid: str) -> dict[str, Any]:
    return next(item for item in candidate["artifacts"] if item["artifact_id"] == aid)


def payload(candidate: dict[str, Any], aid: str) -> dict[str, Any]:
    return v.strict_json_loads(base64.b64decode(artifact(candidate, aid)["payload_b64"]))


def replace_payload(candidate: dict[str, Any], aid: str, new_payload: dict[str, Any]) -> None:
    item = artifact(candidate, aid)
    raw = v.canonical_bytes(new_payload)
    item["payload_b64"] = base64.b64encode(raw).decode()
    item["sha256"] = hashlib.sha256(raw).hexdigest()
    item["byte_length"] = len(raw)
    registry = {x["artifact_id"]: x for x in candidate["artifacts"]}
    refs = sorted(v.referenced_artifact_ids(new_payload))
    item["parents"] = [{"artifact_id": ref, "sha256": registry.get(ref, {"sha256": "0" * 64})["sha256"]} for ref in refs]
    for child in candidate["artifacts"]:
        for parent in child["parents"]:
            if parent["artifact_id"] in registry:
                parent["sha256"] = registry[parent["artifact_id"]]["sha256"]


def change_payload(candidate: dict[str, Any], aid: str, change: Callable[[dict[str, Any]], None]) -> None:
    value = payload(candidate, aid)
    change(value)
    replace_payload(candidate, aid, value)


def result_text(result: dict[str, Any]) -> str:
    return json.dumps(result, sort_keys=True)


def proof_suppressed(result: dict[str, Any]) -> bool:
    outputs = result.get("authoritative_proof_outputs", {})
    return not outputs or all(value == "NO_PROOF_OUTPUT" for value in outputs.values())


def mutate_v001_01_empty_registry(c: dict[str, Any]) -> None:
    c["artifacts"] = []


def mutate_v001_02_injected_custody(c: dict[str, Any]) -> None:
    c["custody_status"] = "QUALIFIED"


def mutate_v001_03_seed_self_authorized(c: dict[str, Any]) -> None:
    change_payload(c, "SEED.JOINT", lambda x: x.update(principal_decision_artifact_id="FRAMEWORK.00"))


def mutate_v001_04_current_without_dilation(c: dict[str, Any]) -> None:
    c["bindings"]["dilation_artifact_id"] = None
    change_payload(c, "SEED.JOINT", lambda x: x.update(dilation_artifact_id=None))


def mutate_v001_05_gamma_projection(c: dict[str, Any]) -> None:
    change_payload(c, "GAMMA.0", lambda x: x.update(representation_type="SCALAR_PROJECTION"))


def mutate_v001_06_dv_overlap(c: dict[str, Any]) -> None:
    d = payload(c, "DATA.DEVELOPMENT")
    change_payload(c, "DATA.VALIDATION", lambda x: x.update(source_ids=d["source_ids"]))


def mutate_v001_07_validation_preaccess(c: dict[str, Any]) -> None:
    def change(x: dict[str, Any]) -> None:
        next(event for event in x["events"] if event["dataset_id"] == "V-SYNTH" and event["event_type"] == "NUMERIC_ACCESS")["time"] = "2026-08-22T09:00:00Z"
    change_payload(c, "DATA.ACCESS_LOG", change)


def mutate_v001_08_ancestry_cycle(c: dict[str, Any]) -> None:
    def change(x: dict[str, Any]) -> None:
        x["arrows"][0]["target"] = x["arrows"][0]["source"]
    change_payload(c, "ANCESTRY.GRAPH", change)


def mutate_v001_09_failed_check_promoted(c: dict[str, Any]) -> None:
    def change(x: dict[str, Any]) -> None:
        next(check for check in x["arrows"][0]["checks"] if check["check_id"] == "IDENTITY")["observed"] = 0.0
    change_payload(c, "ANCESTRY.GRAPH", change)


def mutate_v001_10_missing_power(c: dict[str, Any]) -> None:
    c["gate_inputs"][0]["power_artifact_id"] = "POWER.MISSING"


def mutate_v001_11_missing_joined_evidence(c: dict[str, Any]) -> None:
    change_payload(c, "EVIDENCE.JOINED", lambda x: x.update(response_artifact_ids=[]))


def mutate_v001_12_all_gravity_na(c: dict[str, Any]) -> None:
    for row in c["gravity_applicability"]:
        row["applicability"] = "NOT_APPLICABLE"
        row["justification_artifact_id"] = None


def mutate_v001_13_empty_ge2_domain(c: dict[str, Any]) -> None:
    change_payload(c, "DOMAIN.GE2", lambda x: x.update(member_ids=[]))


def mutate_v001_14_gf2_overlap_promotion(c: dict[str, Any]) -> None:
    mutate_v001_06_dv_overlap(c)
    c["proof_outputs"] = {"GF2": "PASSES_DECLARED_DOMAIN"}


def mutate_v001_15_uge_without_evidence(c: dict[str, Any]) -> None:
    change_payload(c, "COVERAGE.THEOREM", lambda x: x.update(checker_artifact_ids=[]))


def mutate_v001_16_legacy_gate(c: dict[str, Any]) -> None:
    c["gate_ids"][0] = "LEGACY.V003.P7"


def mutate_v001_17_product_self_attestation(c: dict[str, Any]) -> None:
    c["product_outcome"] = "PASS"


def mutate_v001_18_digest_reuse(c: dict[str, Any]) -> None:
    left = artifact(c, "GAMMA.0")
    right = artifact(c, "GAMMA.1")
    right["sha256"] = left["sha256"]
    right["byte_length"] = left["byte_length"]
    right["payload_b64"] = left["payload_b64"]


def mutate_v001_19_proof_injection(c: dict[str, Any]) -> None:
    c["authoritative_proof_outputs"] = {claim: "PASSES_DECLARED_DOMAIN" for claim in v.CLAIMS}


def transport_change(c: dict[str, Any], change: Callable[[dict[str, Any]], None]) -> None:
    change_payload(c, "TRANSPORT.V004", change)


def tr_w_not_psd(x: dict[str, Any]) -> None:
    x["scales"][0]["W"][1][1] = -1.0


def tr_rank_ambiguous(x: dict[str, Any]) -> None:
    x["scales"][0]["eigenvalue_intervals"][0] = [0.05, 1.1]


def tr_btilde_mismatch(x: dict[str, Any]) -> None:
    x["edges"][0]["Btilde"][0][0] = 0.5


def tr_causal_reverse(x: dict[str, Any]) -> None:
    x["edges"][0]["source"], x["edges"][0]["target"] = x["edges"][0]["target"], x["edges"][0]["source"]


def tr_raw_ancestry(x: dict[str, Any]) -> None:
    x["ancestry_component"] = "F"


def tr_covariance_missing(x: dict[str, Any]) -> None:
    x["covariance"]["shared_cross_scale"] = False


def tr_composition_failure(x: dict[str, Any]) -> None:
    factor = 1.02
    direct = next(edge for edge in x["edges"] if edge["edge_id"] == "E20")
    direct["F"] = [[factor, 0.0], [0.0, factor]]
    direct["Btilde"] = [[factor, 0.0], [0.0, factor]]
    x["composition"]["residual_norm"] = (2.0 * (factor - 1.0) ** 2) ** 0.5


def tr_empty_support(x: dict[str, Any]) -> None:
    x["composition"]["P_in_common"] = v.zeros(2, 2)


def tr_edge_underpowered(x: dict[str, Any]) -> None:
    x["edges"][0]["power"]["achieved_power_lower"] = 0.5


def tr_missing_refinement(x: dict[str, Any]) -> None:
    x["refinements"] = [control for control in x["refinements"] if control["type"] != "DUPLICATE"]


def tr_subdivision_alias(x: dict[str, Any]) -> None:
    next(control for control in x["refinements"] if control["type"] == "SUBDIVIDE")["type"] = "SUBDIVISION"


def tr_ref_identity(x: dict[str, Any]) -> None:
    x["refinements"][0]["same_event"] = False


def tr_ref_embedding(x: dict[str, Any]) -> None:
    control = next(item for item in x["refinements"] if item["type"] == "SUBDIVIDE")
    control["scales"][0]["R"] = v.zeros(2, 3)


def tr_ref_artificial_rank(x: dict[str, Any]) -> None:
    control = next(item for item in x["refinements"] if item["type"] == "SUBDIVIDE")
    control["scales"][0]["W_refined"] = v.eye(3)
    control["scales"][0]["P_refined"] = v.eye(3)
    control["scales"][0]["N_refined"] = v.zeros(3, 3)


def tr_ref_raw_naturality(x: dict[str, Any]) -> None:
    control = next(item for item in x["refinements"] if item["type"] == "DUPLICATE")
    control["edges"][0]["F_refined"][0][0] *= 0.5


def tr_ref_calibrated_naturality(x: dict[str, Any]) -> None:
    control = next(item for item in x["refinements"] if item["type"] == "RELABEL")
    control["scales"][0]["U"] = v.eye(2)


def tr_ref_covariance(x: dict[str, Any]) -> None:
    x["refinements"][0]["covariance"]["includes_cross_covariance"] = False


def tr_ref_new_mode(x: dict[str, Any]) -> None:
    x["refinements"][0]["new_independent_mode"] = True


def tr_ref_underpowered(x: dict[str, Any]) -> None:
    x["refinements"][0]["power"]["achieved_power_lower"] = 0.2


def tr_freeze_chronology(x: dict[str, Any]) -> None:
    x["response_access_time"] = "2026-08-22T09:00:00Z"


Case = tuple[str, Callable[[dict[str, Any]], None], str, str]
V001_CASES: list[Case] = [
    ("V001-01-empty-registry", mutate_v001_01_empty_registry, "ARTIFACT_REGISTRY_EMPTY", "REFUSE"),
    ("V001-02-injected-custody", mutate_v001_02_injected_custody, "UNKNOWN_KEYS:custody_status", "REFUSE"),
    ("V001-03-seed-self-authorization", mutate_v001_03_seed_self_authorized, "SEED_AUTHORIZATION_SCOPE_MISMATCH", "REFUSE"),
    ("V001-04-current-without-dilation", mutate_v001_04_current_without_dilation, "PHYSICAL_CURRENT_DILATION_REQUIRED", "REFUSE"),
    ("V001-05-gamma-projection", mutate_v001_05_gamma_projection, "GAMMA_PROCESS_KIND_INVALID", "REFUSE"),
    ("V001-06-dv-overlap", mutate_v001_06_dv_overlap, "DV_OVERLAP_SOURCE", "NO_GF2"),
    ("V001-07-validation-preaccess", mutate_v001_07_validation_preaccess, "VALIDATION_ACCESSED_BEFORE_FREEZE", "NO_GF2"),
    ("V001-08-ancestry-cycle", mutate_v001_08_ancestry_cycle, "ANCESTRY_SELF_LOOP", "REFUSE"),
    ("V001-09-failed-check-promoted", mutate_v001_09_failed_check_promoted, "FAIL_ANCESTRY_IDENTITY", "NO_GF1"),
    ("V001-10-missing-power", mutate_v001_10_missing_power, "BINDING_UNRESOLVED:POWER.MISSING", "REFUSE"),
    ("V001-11-missing-joined-evidence", mutate_v001_11_missing_joined_evidence, "JOINED_EVIDENCE_MISSING_OR_EMPTY", "NO_GE1"),
    ("V001-12-all-gravity-na", mutate_v001_12_all_gravity_na, "MANDATORY_GRAVITY_GATE_NOT_APPLICABLE", "REFUSE"),
    ("V001-13-empty-domain", mutate_v001_13_empty_ge2_domain, "DOMAIN_EMPTY:GE2", "REFUSE"),
    ("V001-14-gf2-overlap-promotion", mutate_v001_14_gf2_overlap_promotion, "UNKNOWN_KEYS:proof_outputs", "REFUSE"),
    ("V001-15-uge-without-evidence", mutate_v001_15_uge_without_evidence, "CONSTRUCTIVE_COVERAGE_INCOMPLETE", "NO_UGE"),
    ("V001-16-legacy-gate", mutate_v001_16_legacy_gate, "GATE_REGISTRY_INVALID", "REFUSE"),
    ("V001-17-product-self-attestation", mutate_v001_17_product_self_attestation, "UNKNOWN_KEYS:product_outcome", "REFUSE"),
    ("V001-18-digest-reuse", mutate_v001_18_digest_reuse, "ARTIFACT_DIGEST_REUSE", "REFUSE"),
    ("V001-19-proof-injection", mutate_v001_19_proof_injection, "UNKNOWN_KEYS:authoritative_proof_outputs", "REFUSE"),
]

V004_CASES: list[Case] = [
    ("V004-01-W-not-PSD", lambda c: transport_change(c, tr_w_not_psd), "W not positive semidefinite", "NO_SCALE"),
    ("V004-02-rank-ambiguous", lambda c: transport_change(c, tr_rank_ambiguous), "UNSCOREABLE_RANK", "NO_SCALE"),
    ("V004-03-Btilde-mismatch", lambda c: transport_change(c, tr_btilde_mismatch), "BTILDE_RECOMPUTE_MISMATCH", "NO_SCALE"),
    ("V004-04-causal-reverse", lambda c: transport_change(c, tr_causal_reverse), "RAW_TRANSPORT_CAUSALITY_INVALID", "NO_SCALE"),
    ("V004-05-raw-ancestry", lambda c: transport_change(c, tr_raw_ancestry), "ANCESTRY_USES_UNCALIBRATED_TRANSPORT", "NO_SCALE"),
    ("V004-06-covariance-missing", lambda c: transport_change(c, tr_covariance_missing), "UNSCOREABLE_TRANSPORT_COVARIANCE", "NO_SCALE"),
    ("V004-07-composition-failure", lambda c: transport_change(c, tr_composition_failure), "FAIL_CALIBRATED_COMPOSITION", "NO_SCALE"),
    ("V004-08-empty-common-support", lambda c: transport_change(c, tr_empty_support), "UNSCOREABLE_COMMON_SUPPORT", "NO_SCALE"),
    ("V004-09-edge-underpowered", lambda c: transport_change(c, tr_edge_underpowered), "UNSCOREABLE_UNDERPOWERED:E10", "NO_SCALE"),
    ("V004-10-refinement-missing", lambda c: transport_change(c, tr_missing_refinement), "REFINEMENT_COVERAGE_INCOMPLETE", "NO_SCALE"),
    ("V004-11-subdivision-alias", lambda c: transport_change(c, tr_subdivision_alias), "REFINEMENT_COVERAGE_INCOMPLETE", "NO_SCALE"),
    ("V004-12-refinement-identity", lambda c: transport_change(c, tr_ref_identity), "UNSCOREABLE_REFINEMENT_IDENTITY", "NO_SCALE"),
    ("V004-13-refinement-embedding", lambda c: transport_change(c, tr_ref_embedding), "REFINEMENT_EMBEDDING_INVALID", "NO_SCALE"),
    ("V004-14-refinement-artificial-rank", lambda c: transport_change(c, tr_ref_artificial_rank), "REFINEMENT_RANK_PROJECTOR_INVALID", "NO_SCALE"),
    ("V004-15-raw-naturality", lambda c: transport_change(c, tr_ref_raw_naturality), "FAIL_RAW_REFINEMENT_NATURALITY", "NO_SCALE"),
    ("V004-16-calibrated-naturality", lambda c: transport_change(c, tr_ref_calibrated_naturality), "REFINEMENT_U_DERIVATION_INVALID", "NO_SCALE"),
    ("V004-17-refinement-covariance", lambda c: transport_change(c, tr_ref_covariance), "UNSCOREABLE_REFINEMENT_COVARIANCE", "NO_SCALE"),
    ("V004-18-new-independent-mode", lambda c: transport_change(c, tr_ref_new_mode), "UNSCOREABLE_NEW_INDEPENDENT_MODE", "NO_SCALE"),
    ("V004-19-refinement-underpowered", lambda c: transport_change(c, tr_ref_underpowered), "UNSCOREABLE_REFINEMENT_UNDERPOWERED", "NO_SCALE"),
    ("V004-20-freeze-chronology", lambda c: transport_change(c, tr_freeze_chronology), "TRANSPORT_FREEZE_CHRONOLOGY_VIOLATION", "NO_SCALE"),
]


def expected_mode(result: dict[str, Any], mode: str) -> bool:
    if mode == "REFUSE":
        return result.get("accepted") is False
    if mode == "NO_GF1":
        return result.get("candidate_milestones", {}).get("GF1") != "PASS"
    if mode == "NO_GF2":
        return result.get("candidate_milestones", {}).get("GF2") != "PASS"
    if mode == "NO_GE1":
        return result.get("candidate_milestones", {}).get("GE1") != "PASS"
    if mode == "NO_UGE":
        return result.get("candidate_milestones", {}).get("UGE") != "PASS"
    if mode == "NO_SCALE":
        return result.get("scientific_gates", {}).get("SCALE.CALIBRATED_TRANSPORT") != "PASS"
    raise ValueError(mode)


def original_witness() -> dict[str, Any]:
    audit_path = WORKSPACE / "LANE_GRA_O_GF_CONTRACT/VERIFY_CODEX/audit_contract.py"
    spec = importlib.util.spec_from_file_location("v001_audit", audit_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load V001 audit")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    schema = json.loads((WORKSPACE / "LANE_GRA_O_GF_CONTRACT/CONTRACT.schema.json").read_text())
    return module.contradictory_witness(schema, module.SchemaEvaluator(schema))


def run_suite() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    baseline = fixture.build_fixture()
    baseline_result = v.evaluate_instance(copy.deepcopy(baseline), WORKSPACE)
    baseline_ok = (
        baseline_result.get("accepted") is True
        and baseline_result.get("custody", {}).get("disposition") == "QUALIFIED"
        and all(value == "PASS" for value in baseline_result.get("scientific_gates", {}).values())
        and baseline_result.get("product_reproduction") == "PASS"
        and baseline_result.get("scientific_weight_of_product") == "NONE"
        and baseline_result.get("actual_platform_present") is False
        and proof_suppressed(baseline_result)
    )
    rows.append({"case": "POSITIVE-synthetic-complete-no-platform", "pass": baseline_ok})

    witness = original_witness()
    witness_hash = hashlib.sha256(v.canonical_bytes(witness)).hexdigest()
    witness_result = v.evaluate_instance(witness, WORKSPACE)
    witness_ok = witness_hash == "d55aaab13a908afc9425746ddd18ae85924fb5ece338b215c574c6fe38320855" and witness_result.get("accepted") is False
    rows.append({"case": "V001-original-contradictory-witness", "pass": witness_ok, "witness_sha256": witness_hash})

    for name, mutation, expected, mode in V001_CASES + V004_CASES:
        candidate = copy.deepcopy(baseline)
        mutation(candidate)
        result = v.evaluate_instance(candidate, WORKSPACE)
        text = result_text(result)
        ok = expected in text and expected_mode(result, mode) and proof_suppressed(result)
        rows.append({"case": name, "pass": ok, "expected": expected, "mode": mode})

    duplicate_key_ok = False
    nonfinite_ok = False
    try:
        v.strict_json_loads('{"x":1,"x":2}')
    except v.DuplicateKeyError:
        duplicate_key_ok = True
    try:
        v.strict_json_loads('{"x":NaN}')
    except ValueError:
        nonfinite_ok = True
    rows.append({"case": "PARSER-duplicate-key", "pass": duplicate_key_ok})
    rows.append({"case": "PARSER-nonfinite", "pass": nonfinite_ok})

    failed = [row["case"] for row in rows if not row["pass"]]
    result = {
        "suite_id": "GRA-O-GF-CONTRACT-V002-FIXED-MUTATIONS",
        "overall_result": "PASS" if not failed else "FAIL",
        "case_count": len(rows),
        "positive_controls": 1,
        "original_witness_cases": 1,
        "v001_bypass_families": len(V001_CASES),
        "v004_transport_refinement_guards": len(V004_CASES),
        "parser_guards": 2,
        "failed_count": len(failed),
        "failed_cases": failed,
        "cases": rows,
    }
    result["result_payload_sha256"] = hashlib.sha256(v.canonical_bytes(result)).hexdigest()
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    result = run_suite()
    print(json.dumps(result, sort_keys=True, indent=None if args.compact else 2))
    return 0 if result["overall_result"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
