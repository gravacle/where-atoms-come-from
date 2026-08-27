#!/usr/bin/env python3
"""Fixed adversarial checks for the missing-data-aware proof frontier."""

from __future__ import annotations

import copy
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Any, Callable

from proof_frontier import ProofFrontierRefusal, load_proof_frontier


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_payload(root: Path, relative: str, content: str) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def object_row(
    root: Path,
    object_id: str,
    relative: str,
    content: str,
    evidence_class: str,
    dataset_id: str | None,
    stage: str,
    role: str,
) -> dict[str, Any]:
    path = write_payload(root, relative, content)
    return {
        "object_id": object_id,
        "path": relative,
        "sha256": digest(path),
        "media_type": "application/json" if relative.endswith(".json") else "text/plain",
        "evidence_class": evidence_class,
        "dataset_id": dataset_id,
        "stage": stage,
        "role": role,
    }


def dataset_row(
    dataset_id: str,
    role: str,
    availability: str,
    evidence_class: str,
    opened: bool,
    ancestry: str,
    namespace: str,
    upstream_status: str,
    object_ids: list[str],
    qualification_object_ids: list[str],
    *,
    prior_use: str = "NONE",
    independent: list[str] | None = None,
    measured: list[str] | None = None,
    missing: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "dataset_id": dataset_id,
        "role": role,
        "availability": availability,
        "evidence_class": evidence_class,
        "public": True,
        "opened": opened,
        "source_id": f"PUBLIC:{dataset_id}",
        "acquisition_id": f"ACQ:{dataset_id}",
        "physical_unit_namespace": namespace,
        "acquisition_ancestry": ancestry,
        "prior_use": prior_use,
        "independent_from_dataset_ids": list(independent or []),
        "frozen_before_response_access": True,
        "role_locked_before_response_access": True,
        "upstream_status": upstream_status,
        "qualification_object_ids": qualification_object_ids,
        "measured_channel_ids": list(measured or []),
        "missing_channel_ids": list(missing or []),
        "object_ids": object_ids,
    }


def derivation(mode: str = "NONE", transform: str | None = None) -> dict[str, Any]:
    lawful = mode == "FROZEN_IDENTIFIABLE_TRANSFORM"
    return {
        "mode": mode,
        "transform_object_id": transform,
        "frozen_before_response_access": lawful,
        "identifiable_on_declared_support": lawful,
        "uncertainty_propagated": lawful,
        "uses_tested_response_as_control": False,
    }


def obligation(
    obligation_id: str,
    proof_id: str,
    status: str,
    *,
    required: bool = True,
    dependencies: list[str] | None = None,
    datasets: list[str] | None = None,
    evidence: list[str] | None = None,
    missing: list[str] | None = None,
    derivation_row: dict[str, Any] | None = None,
    theories: list[str] | None = None,
) -> dict[str, Any]:
    return {
        "obligation_id": obligation_id,
        "proof_id": proof_id,
        "label": obligation_id.replace("_", " "),
        "required": required,
        "dependencies": list(dependencies or []),
        "status": status,
        "dataset_ids": list(datasets or []),
        "evidence_object_ids": list(evidence or []),
        "derivation": derivation_row or derivation(),
        "missing_channel_ids": list(missing or []),
        "theory_track_ids": list(theories or []),
    }


def build_manifest(root: Path) -> dict[str, Any]:
    objects = [
        object_row(root, "policy", "policy/POLICY.md", "public-data-only\n", "THEORY_ONLY", None, "POLICY", "PUBLIC_DATA_POLICY"),
        object_row(root, "amendment", "policy/A001.json", "{}\n", "THEORY_ONLY", None, "POLICY", "POLICY_AMENDMENT"),
        object_row(root, "dev1_qual", "dev1/qual.json", "qualified\n", "REAL_WORLD_EVENT_LEVEL", "dev1", "RESULT", "EVENT_LEVEL_QUALIFICATION_CERTIFICATE"),
        object_row(root, "dev1_raw", "dev1/raw.bin", "raw-event-1\n", "REAL_WORLD_EVENT_LEVEL", "dev1", "RAW_EVENT", "RAW_EVENT_DATA"),
        object_row(root, "dev1_result", "dev1/result.json", "prediction-pass\n", "REAL_WORLD_EVENT_LEVEL", "dev1", "RESULT", "THEORY_TEST_RESULT"),
        object_row(root, "dev1_transform", "dev1/transform.json", "frozen-transform\n", "REAL_WORLD_EVENT_LEVEL", "dev1", "DERIVATION", "TRANSFORM_SPECIFICATION"),
        object_row(root, "dev2_qual", "dev2/qual.json", "qualified\n", "REAL_WORLD_EVENT_LEVEL", "dev2", "RESULT", "EVENT_LEVEL_QUALIFICATION_CERTIFICATE"),
        object_row(root, "dev2_raw", "dev2/raw.bin", "raw-event-2\n", "REAL_WORLD_EVENT_LEVEL", "dev2", "RAW_EVENT", "RAW_EVENT_DATA"),
        object_row(root, "val_qual", "val/qual.json", "qualified\n", "REAL_WORLD_EVENT_LEVEL", "val", "RESULT", "EVENT_LEVEL_QUALIFICATION_CERTIFICATE"),
        object_row(root, "val_raw", "val/raw.bin", "heldout-event\n", "REAL_WORLD_EVENT_LEVEL", "val", "VALIDATION_RESPONSE", "RAW_EVENT_DATA"),
        object_row(root, "val_result", "val/result.json", "heldout-pass\n", "REAL_WORLD_EVENT_LEVEL", "val", "RESULT", "VALIDATION_TEST_RESULT"),
        object_row(root, "missing_meta", "missing/metadata.json", "metadata-only\n", "REAL_WORLD_EVENT_LEVEL", "missing", "METADATA", "METADATA"),
        object_row(root, "paper_summary", "paper/summary.json", "reported interval\n", "PAPER_LEVEL_SUMMARY", "paper", "PAPER_SUMMARY", "PAPER_SUMMARY_DATA"),
        object_row(root, "theory_math", "theory/derivation.md", "theorem\n", "THEORY_ONLY", "theory", "THEORY", "THEORY_DERIVATION"),
    ]
    datasets = [
        dataset_row(
            "dev1", "DEVELOPMENT", "AVAILABLE_PUBLIC", "REAL_WORLD_EVENT_LEVEL", True,
            "ANCESTRY-DEV1", "UNITS-DEV1", "QUALIFIED_EVENT_LEVEL",
            ["dev1_qual", "dev1_raw", "dev1_result", "dev1_transform"], ["dev1_qual"],
            prior_use="MODEL_SELECTION", measured=["BEFORE", "FORMATION", "D_CAUSAL"],
        ),
        dataset_row(
            "dev2", "REPRODUCTION", "AVAILABLE_PUBLIC", "REAL_WORLD_EVENT_LEVEL", True,
            "ANCESTRY-DEV2", "UNITS-DEV2", "QUALIFIED_EVENT_LEVEL",
            ["dev2_qual", "dev2_raw"], ["dev2_qual"], independent=["dev1"],
            measured=["BEFORE", "FORMATION", "D_CAUSAL"],
        ),
        dataset_row(
            "val", "VALIDATION", "AVAILABLE_PUBLIC", "REAL_WORLD_EVENT_LEVEL", True,
            "ANCESTRY-VAL", "UNITS-VAL", "QUALIFIED_EVENT_LEVEL",
            ["val_qual", "val_raw", "val_result"], ["val_qual"],
            independent=["dev1", "dev2"], measured=["VALIDATION_RESPONSE"],
        ),
        dataset_row(
            "missing", "DEVELOPMENT", "UNAVAILABLE_PUBLIC_DATA", "REAL_WORLD_EVENT_LEVEL", False,
            "ANCESTRY-MISSING", "UNITS-MISSING", "PARTIAL_OR_UNSCOREABLE",
            ["missing_meta"], [], missing=["PHYSICAL_G", "REMOTE_PROBE"],
        ),
        dataset_row(
            "paper", "DEVELOPMENT", "AVAILABLE_PUBLIC", "PAPER_LEVEL_SUMMARY", True,
            "ANCESTRY-PAPER", "UNITS-PAPER", "PAPER_SUMMARY_ONLY",
            ["paper_summary"], [], prior_use="MODEL_SELECTION", measured=["REPORTED_INTERVAL"],
        ),
        dataset_row(
            "theory", "THEORY_ONLY", "AVAILABLE_PUBLIC", "THEORY_ONLY", True,
            "ANCESTRY-THEORY", "UNITS-THEORY", "THEORY_ONLY", ["theory_math"], [],
        ),
    ]
    obligations = [
        obligation("u_before", "URF", "PASS_DIRECT_PUBLIC", datasets=["dev1"], evidence=["dev1_raw"], derivation_row=derivation("DIRECT"), theories=["T_RECORD"]),
        obligation("u_closure", "URF", "UNSCOREABLE_MISSING_DATA", dependencies=["u_before"], datasets=["missing"], missing=["PHYSICAL_G"]),
        obligation("u_hold", "URF", "BLOCKED_BY_PARENT", dependencies=["u_closure"]),
        obligation("u_cal", "URF", "NOT_RUN_AVAILABLE", datasets=["dev2"]),
        obligation("u_theorem", "URF", "THEOREM_ONLY", required=False, datasets=["theory"], evidence=["theory_math"], theories=["T_RECORD"]),
        obligation("g_transport", "GE", "PASS_DIRECT_PUBLIC", datasets=["dev1"], evidence=["dev1_raw"], derivation_row=derivation("DIRECT"), theories=["T_RECORD"]),
        obligation("g_probe", "GE", "UNSCOREABLE_MISSING_DATA", dependencies=["g_transport"], datasets=["missing"], missing=["REMOTE_PROBE"]),
        obligation("g_endpoint", "GE", "BLOCKED_BY_PARENT", dependencies=["g_probe"]),
        obligation("g_cal", "GE", "NOT_RUN_AVAILABLE", datasets=["dev2"]),
    ]
    return {
        "schema": "WAC_PROOF_FRONTIER_V001",
        "frontier_id": "SYNTHETIC_FRONTIER_V001",
        "profile_id": "SYNTHETIC_TEST_PROFILE_V001",
        "policy": {
            "goal": "PROOF",
            "public_data_only": True,
            "validation_preserved": True,
            "theory_fallback_rule": "SUBSTANTIALLY_PROVEN_ONLY_UNAVAILABLE_DATA",
            "policy_object_ids": ["policy", "amendment"],
        },
        "external_actions": {
            "authorized": False,
            "executed": False,
            "third_party_contact": False,
            "private_data": False,
            "new_acquisition": False,
        },
        "proofs": [
            {
                "proof_id": "URF",
                "target": "UNIVERSAL_RECORD_FORMATION_PROOF",
                "required_obligation_ids": ["u_before", "u_closure", "u_hold", "u_cal"],
                "authoritative_output": "NO_PROOF_OUTPUT",
            },
            {
                "proof_id": "GE",
                "target": "GRAVITY_EMERGENCE_PROOF",
                "required_obligation_ids": ["g_transport", "g_probe", "g_endpoint", "g_cal"],
                "authoritative_output": "NO_PROOF_OUTPUT",
            },
        ],
        "datasets": datasets,
        "objects": objects,
        "obligations": obligations,
        "theory_tracks": [
            {
                "theory_id": "T_RECORD",
                "label": "record-process theory",
                "mathematics_closed": True,
                "adversarial_checks_pass": True,
                "predictions": [
                    {
                        "prediction_id": "paper_prediction",
                        "status": "PASS",
                        "dataset_ids": ["paper"],
                        "evidence_object_ids": ["paper_summary"],
                    },
                    {
                        "prediction_id": "event_prediction",
                        "status": "PASS",
                        "dataset_ids": ["dev1"],
                        "evidence_object_ids": ["dev1_result"],
                    },
                    {
                        "prediction_id": "unavailable_prediction",
                        "status": "NOT_RUN",
                        "dataset_ids": ["missing"],
                        "evidence_object_ids": [],
                    },
                ],
                "independent_event_level_domain_ids": ["dev1"],
                "held_out_validation_result_object_ids": [],
                "remaining_blockers": [
                    {
                        "blocker_id": "B_DATA",
                        "blocker_class": "UNAVAILABLE_PUBLIC_DATA",
                        "description": "physical closure and remote response are not public",
                        "affected_obligation_ids": ["u_closure", "g_probe"],
                        "missing_channel_ids": ["PHYSICAL_G", "REMOTE_PROBE"],
                    }
                ],
            }
        ],
    }


def write_manifest(root: Path, manifest: dict[str, Any], *, allow_nan: bool = False) -> Path:
    path = root / "frontier.json"
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=allow_nan) + "\n",
        encoding="utf-8",
    )
    return path


def run_case(mutator: Callable[[Path, dict[str, Any]], None] | None = None):
    with tempfile.TemporaryDirectory(prefix="wac-proof-frontier-") as directory:
        root = Path(directory) / "bundle"
        root.mkdir()
        manifest = build_manifest(root)
        if mutator:
            mutator(root, manifest)
        return load_proof_frontier(write_manifest(root, manifest))


def refuses(mutator: Callable[[Path, dict[str, Any]], None], *, allow_nan: bool = False) -> bool:
    with tempfile.TemporaryDirectory(prefix="wac-proof-frontier-") as directory:
        root = Path(directory) / "bundle"
        root.mkdir()
        manifest = build_manifest(root)
        mutator(root, manifest)
        path = write_manifest(root, manifest, allow_nan=allow_nan)
        try:
            load_proof_frontier(path)
        except ProofFrontierRefusal:
            return True
        return False


def duplicate_json_refuses() -> bool:
    with tempfile.TemporaryDirectory(prefix="wac-proof-frontier-") as directory:
        root = Path(directory) / "bundle"
        root.mkdir()
        manifest = build_manifest(root)
        payload = json.dumps(manifest, sort_keys=True)
        payload = payload[:-1] + ',"schema":"WAC_PROOF_FRONTIER_V001"}'
        path = root / "frontier.json"
        path.write_text(payload, encoding="utf-8")
        try:
            load_proof_frontier(path)
        except ProofFrontierRefusal:
            return True
        return False


def obligation_row(manifest: dict[str, Any], obligation_id: str) -> dict[str, Any]:
    return next(row for row in manifest["obligations"] if row["obligation_id"] == obligation_id)


def dataset(manifest: dict[str, Any], dataset_id: str) -> dict[str, Any]:
    return next(row for row in manifest["datasets"] if row["dataset_id"] == dataset_id)


def obj(manifest: dict[str, Any], object_id: str) -> dict[str, Any]:
    return next(row for row in manifest["objects"] if row["object_id"] == object_id)


def theory(manifest: dict[str, Any]) -> dict[str, Any]:
    return manifest["theory_tracks"][0]


def substantial(root: Path, manifest: dict[str, Any]) -> None:
    del root
    track = theory(manifest)
    track["independent_event_level_domain_ids"] = ["dev1", "dev2"]
    track["held_out_validation_result_object_ids"] = ["val_result"]
    track["predictions"].extend(
        [
            {
                "prediction_id": "event_prediction_2",
                "status": "PASS",
                "dataset_ids": ["dev2"],
                "evidence_object_ids": ["dev2_raw"],
            },
            {
                "prediction_id": "heldout_validation_prediction",
                "status": "PASS",
                "dataset_ids": ["val"],
                "evidence_object_ids": ["val_result"],
            },
        ]
    )


def all_proof_obligations_pass(root: Path, manifest: dict[str, Any]) -> None:
    del root
    for row in manifest["obligations"]:
        if not row["required"]:
            continue
        row["status"] = "PASS_DIRECT_PUBLIC"
        row["dataset_ids"] = ["dev1"]
        row["evidence_object_ids"] = ["dev1_raw"]
        row["missing_channel_ids"] = []
        row["derivation"] = derivation("DIRECT")


def main() -> int:
    checks = 0

    def check(condition: bool, label: str) -> None:
        nonlocal checks
        if not condition:
            raise AssertionError(label)
        checks += 1

    baseline = run_case()
    certificate = baseline.certificate()
    check(baseline.manifest["frontier_id"] == "SYNTHETIC_FRONTIER_V001", "baseline load")
    check(set(certificate["authoritative_proof_outputs"].values()) == {"NO_PROOF_OUTPUT"}, "no proof output")
    check(baseline.proof_states == {"URF": "BLOCKED_MISSING_DATA", "GE": "BLOCKED_MISSING_DATA"}, "blocked proof states")
    check(baseline.theory_states == {"T_RECORD": "PARTIAL_SUPPORT_ONLY"}, "partial theory state")
    check(certificate["paper_level_theory_tests"] == [{"theory_id": "T_RECORD", "prediction_id": "paper_prediction", "status": "PASS", "proof_eligible": False}], "paper result classified")
    check(certificate["paper_level_proof_credit"] == 0, "paper has zero proof credit")
    check({row["obligation_id"] for row in certificate["missing_data_blockers"]} == {"u_closure", "g_probe"}, "missing blockers surfaced")
    check(set(certificate["available_execution_frontier"]) == {"u_cal", "g_cal"}, "available work continues")
    check(run_case().certificate()["proof_states"] == certificate["proof_states"], "deterministic reload")

    check(refuses(lambda r, m: m.__setitem__("schema", "WRONG")), "schema closure")
    check(refuses(lambda r, m: m.__setitem__("extra", True)), "extra root key")
    check(refuses(lambda r, m: m.pop("proofs")), "missing root key")
    check(duplicate_json_refuses(), "duplicate JSON")
    check(refuses(lambda r, m: m.__setitem__("nonfinite", math.nan), allow_nan=True), "nonfinite JSON")
    check(refuses(lambda r, m: m["external_actions"].__setitem__("authorized", True)), "authorization forbidden")
    check(refuses(lambda r, m: m["external_actions"].__setitem__("executed", True)), "execution forbidden")
    check(refuses(lambda r, m: m["external_actions"].__setitem__("third_party_contact", True)), "contact forbidden")
    check(refuses(lambda r, m: m["policy"].__setitem__("goal", "THEORY")), "proof goal immutable")
    check(refuses(lambda r, m: m["policy"].__setitem__("public_data_only", False)), "public-only immutable")
    check(refuses(lambda r, m: m["policy"].__setitem__("validation_preserved", False)), "validation preserved")
    check(refuses(lambda r, m: m["policy"].__setitem__("theory_fallback_rule", "EASY_THEORY")), "fallback rule immutable")
    check(refuses(lambda r, m: m["policy"].__setitem__("policy_object_ids", ["policy"])), "amendment required")
    check(refuses(lambda r, m: obj(m, "dev1_raw").__setitem__("sha256", "0" * 64)), "object hash")

    def symlink_object(root: Path, manifest: dict[str, Any]) -> None:
        row = obj(manifest, "dev1_raw")
        path = root / row["path"]
        target = path.with_name("target.bin")
        path.rename(target)
        os.symlink(target.name, path)
        row["sha256"] = digest(target)

    check(refuses(symlink_object), "object symlink")

    def escape_path(root: Path, manifest: dict[str, Any]) -> None:
        outside = root.parent / "outside.txt"
        outside.write_text("outside\n", encoding="utf-8")
        row = obj(manifest, "dev1_raw")
        row["path"] = "../outside.txt"
        row["sha256"] = digest(outside)

    check(refuses(escape_path), "path escape")
    check(refuses(lambda r, m: m["objects"].append(copy.deepcopy(m["objects"][0]))), "duplicate object id")

    def duplicate_path(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obj(manifest, "dev2_raw")
        source = obj(manifest, "dev1_raw")
        row["path"], row["sha256"] = source["path"], source["sha256"]

    check(refuses(duplicate_path), "duplicate path")
    check(refuses(lambda r, m: obj(m, "theory_math").__setitem__("dataset_id", "ghost")), "unknown object dataset")

    def unlisted_object(root: Path, manifest: dict[str, Any]) -> None:
        del root
        dataset(manifest, "dev1")["object_ids"].remove("dev1_result")

    check(refuses(unlisted_object), "unlisted object")
    check(refuses(lambda r, m: obj(m, "paper_summary").__setitem__("stage", "RAW_EVENT")), "paper stage")
    check(refuses(lambda r, m: obj(m, "paper_summary").__setitem__("evidence_class", "REAL_WORLD_EVENT_LEVEL")), "paper role mismatch")
    check(refuses(lambda r, m: dataset(m, "dev1").__setitem__("public", False)), "nonpublic dataset")
    check(refuses(lambda r, m: dataset(m, "missing").__setitem__("opened", True)), "open unavailable")
    check(refuses(lambda r, m: dataset(m, "paper").__setitem__("upstream_status", "QUALIFIED_EVENT_LEVEL")), "paper qualification laundering")
    check(refuses(lambda r, m: dataset(m, "val").__setitem__("prior_use", "MODEL_SELECTION")), "validation prior use")
    check(refuses(lambda r, m: dataset(m, "val").__setitem__("independent_from_dataset_ids", [])), "validation independence")
    check(refuses(lambda r, m: dataset(m, "dev2").__setitem__("acquisition_ancestry", "ANCESTRY-DEV1")), "duplicate ancestry")

    def paper_pass(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_before")
        row["dataset_ids"] = ["paper"]
        row["evidence_object_ids"] = ["paper_summary"]

    check(refuses(paper_pass), "paper cannot pass proof")
    check(refuses(lambda r, m: dataset(m, "dev1").__setitem__("upstream_status", "PARTIAL_OR_UNSCOREABLE")), "unqualified event cannot pass")

    def qualification_only(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_before")["evidence_object_ids"] = ["dev1_qual"]

    check(refuses(qualification_only), "certificate alone cannot pass")

    def derived_no_transform(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_before")
        row["status"] = "PASS_DERIVED_IDENTIFIABLE"
        row["derivation"] = derivation("FROZEN_IDENTIFIABLE_TRANSFORM")

    check(refuses(derived_no_transform), "derived transform required")

    def derived_not_frozen(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_before")
        row["status"] = "PASS_DERIVED_IDENTIFIABLE"
        row["evidence_object_ids"] = ["dev1_raw", "dev1_transform"]
        row["derivation"] = derivation("FROZEN_IDENTIFIABLE_TRANSFORM", "dev1_transform")
        row["derivation"]["frozen_before_response_access"] = False

    check(refuses(derived_not_frozen), "derivation freeze")

    def derived_response_control(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_before")
        row["status"] = "PASS_DERIVED_IDENTIFIABLE"
        row["evidence_object_ids"] = ["dev1_raw", "dev1_transform"]
        row["derivation"] = derivation("FROZEN_IDENTIFIABLE_TRANSFORM", "dev1_transform")
        row["derivation"]["uses_tested_response_as_control"] = True

    check(refuses(derived_response_control), "response-derived control")
    check(refuses(lambda r, m: obligation_row(m, "u_closure").__setitem__("missing_channel_ids", [])), "missing channel required")

    def unavailable_notrun(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_cal")["dataset_ids"] = ["missing"]

    check(refuses(unavailable_notrun), "not-run must be available")

    def reserved_validation_notrun(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_cal")["dataset_ids"] = ["val"]

    check(refuses(reserved_validation_notrun), "validation is never directly runnable")

    def reserved_public_notrun(root: Path, manifest: dict[str, Any]) -> None:
        del root
        dataset(manifest, "dev2")["availability"] = "RESERVED_UNOPENED_PUBLIC"
        dataset(manifest, "dev2")["opened"] = False

    check(refuses(reserved_public_notrun), "reserved public data are not runnable")

    def unopened_notrun(root: Path, manifest: dict[str, Any]) -> None:
        del root
        dataset(manifest, "dev2")["opened"] = False

    check(refuses(unopened_notrun), "unopened data are not runnable")

    def premature_notrun(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_cal")["dependencies"] = ["u_closure"]

    check(refuses(premature_notrun), "runnable work requires passed parents")
    check(refuses(lambda r, m: obligation_row(m, "u_hold").__setitem__("dependencies", ["ghost"])), "unknown dependency")

    def cycle(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_before")["dependencies"] = ["u_hold"]

    check(refuses(cycle), "dependency cycle")

    def pass_unpassed_parent(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_hold")
        row["status"] = "PASS_DIRECT_PUBLIC"
        row["dataset_ids"] = ["dev1"]
        row["evidence_object_ids"] = ["dev1_raw"]
        row["derivation"] = derivation("DIRECT")

    check(refuses(pass_unpassed_parent), "parent cannot be laundered")
    check(refuses(lambda r, m: m["proofs"][0]["required_obligation_ids"].remove("u_cal")), "proof cannot omit required")

    def theorem_required(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_theorem")
        row["required"] = True
        manifest["proofs"][0]["required_obligation_ids"].append("u_theorem")

    check(refuses(theorem_required), "theorem cannot replace empirical proof")

    def result_without_evidence(root: Path, manifest: dict[str, Any]) -> None:
        del root
        theory(manifest)["predictions"][0]["evidence_object_ids"] = []

    check(refuses(result_without_evidence), "theory result evidence")

    def fail_without_evidence(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_before")
        row["status"] = "FAIL_EMPIRICAL"
        row["evidence_object_ids"] = []

    check(refuses(fail_without_evidence), "empirical failure needs evidence")

    def paper_failure(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_before")
        row["status"] = "FAIL_EMPIRICAL"
        row["dataset_ids"] = ["paper"]
        row["evidence_object_ids"] = ["paper_summary"]

    check(refuses(paper_failure), "paper cannot refute proof")

    def valid_empirical_failure(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_before")["status"] = "FAIL_EMPIRICAL"

    check(run_case(valid_empirical_failure).proof_states["URF"] == "NONAUTHORITATIVE_REPORTED_FAILURE_REQUIRES_AUTHORITATIVE_SCORER", "qualified empirical refutation remains nonauthoritative")

    def failure_with_unpassed_parent(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_hold")
        row["status"] = "FAIL_EMPIRICAL"
        row["dataset_ids"] = ["dev1"]
        row["evidence_object_ids"] = ["dev1_raw"]
        row["derivation"] = derivation("DIRECT")

    check(refuses(failure_with_unpassed_parent), "empirical failure needs passed parents")

    def paper_domain(root: Path, manifest: dict[str, Any]) -> None:
        del root
        theory(manifest)["independent_event_level_domain_ids"].append("paper")

    check(refuses(paper_domain), "paper cannot count as event domain")

    def paper_validation(root: Path, manifest: dict[str, Any]) -> None:
        del root
        theory(manifest)["held_out_validation_result_object_ids"] = ["paper_summary"]

    check(refuses(paper_validation), "paper cannot count as validation")
    check(run_case(substantial).theory_states["T_RECORD"] == "CANDIDATE_SUBSTANTIAL_SUPPORT_REQUIRES_SCIENTIFIC_AUDIT", "substantial theory candidate remains nonauthoritative")

    def disconnected_domain(root: Path, manifest: dict[str, Any]) -> None:
        substantial(root, manifest)
        theory(manifest)["predictions"] = [
            row for row in theory(manifest)["predictions"]
            if row["prediction_id"] != "event_prediction_2"
        ]

    check(refuses(disconnected_domain), "independent domain must have passed prediction")

    def unlinked_heldout(root: Path, manifest: dict[str, Any]) -> None:
        substantial(root, manifest)
        theory(manifest)["predictions"] = [
            row for row in theory(manifest)["predictions"]
            if row["prediction_id"] != "heldout_validation_prediction"
        ]

    check(refuses(unlinked_heldout), "heldout object must have passed prediction")

    def unrun_available(root: Path, manifest: dict[str, Any]) -> None:
        substantial(root, manifest)
        prediction = theory(manifest)["predictions"][2]
        prediction["dataset_ids"] = ["dev2"]

    check(run_case(unrun_available).theory_states["T_RECORD"] == "PARTIAL_SUPPORT_ONLY", "available falsifier cannot be skipped")

    def failed_prediction(root: Path, manifest: dict[str, Any]) -> None:
        del root
        theory(manifest)["predictions"][0]["status"] = "FAIL"

    check(run_case(failed_prediction).theory_states["T_RECORD"] == "REPORTED_REFUTATION_REQUIRES_SCIENTIFIC_AUDIT", "theory refutation requires audit")

    def incomplete_math(root: Path, manifest: dict[str, Any]) -> None:
        del root
        theory(manifest)["mathematics_closed"] = False

    check(run_case(incomplete_math).theory_states["T_RECORD"] == "INCOMPLETE_THEORY", "math closure")

    def theory_choice(root: Path, manifest: dict[str, Any]) -> None:
        substantial(root, manifest)
        theory(manifest)["remaining_blockers"][0]["blocker_class"] = "THEORY_CHOICE"

    check(run_case(theory_choice).theory_states["T_RECORD"] == "PARTIAL_SUPPORT_ONLY", "theory choice blocks fallback")
    readiness = run_case(all_proof_obligations_pass)
    check(set(readiness.proof_states.values()) == {"NONAUTHORITATIVE_INPUT_COMPLETE_REQUIRES_AUTHORITATIVE_SCORER"}, "input completeness is nonauthoritative")
    check(set(readiness.certificate()["authoritative_proof_outputs"].values()) == {"NO_PROOF_OUTPUT"}, "input completeness is not proof")
    check(readiness.certificate()["scientific_readiness_authorized"] is False, "frontier cannot authorize readiness")

    check(refuses(lambda r, m: m.__setitem__("profile_id", "PUBLIC_TWO_PROOF_PROFILE_V001")), "public profile requires canonical obligations")
    check(refuses(lambda r, m: (m.__setitem__("frontier_id", "NOT_SYNTHETIC"))), "synthetic profile requires explicit synthetic frontier")

    def unknown_theory(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_before")["theory_track_ids"] = ["ghost"]

    check(refuses(unknown_theory), "unknown theory binding")

    def nonexistent_parent_block(root: Path, manifest: dict[str, Any]) -> None:
        del root
        row = obligation_row(manifest, "u_hold")
        obligation_row(manifest, "u_closure")["status"] = "PASS_DIRECT_PUBLIC"
        obligation_row(manifest, "u_closure")["dataset_ids"] = ["dev1"]
        obligation_row(manifest, "u_closure")["evidence_object_ids"] = ["dev1_raw"]
        obligation_row(manifest, "u_closure")["missing_channel_ids"] = []
        obligation_row(manifest, "u_closure")["derivation"] = derivation("DIRECT")
        row["status"] = "BLOCKED_BY_PARENT"

    check(refuses(nonexistent_parent_block), "false parent block")

    def proof_mismatch(root: Path, manifest: dict[str, Any]) -> None:
        del root
        obligation_row(manifest, "u_cal")["proof_id"] = "GE"

    check(refuses(proof_mismatch), "proof ownership")

    if checks != 76:
        raise AssertionError(f"check count changed: {checks}/76")
    print(f"PROOF_FRONTIER_CORE_CHECKS: {checks}/76 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
