#!/usr/bin/env python3
"""Author-side exact verifier for the GL6FJ V003 m001 launch binding."""

from __future__ import annotations

import copy
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Callable

import admission_contract_v003 as admission

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "m001_full_cu_c64_projection_v001.cpp"
POST = HERE / "reconstruct_m001_full_response_v002.py"
CONTROLLER = HERE / "run_gl6fj_v003_pilot.py"
RESULT_VERIFIER = HERE / "result_packet_verifier_v003.py"
PLAN = HERE / "MEASUREMENT_PLAN.json"
HISTORY = HERE / "EXPECTED_HISTORY.json"
DEPENDENCIES = HERE / "DEPENDENCIES.sha256"
DEPENDENCY_CONTRACT = HERE / "DEPENDENCY_CONTRACT.json"
RESULT_NAME = "RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
SOURCE_SHA256 = "859e64d2081c1a5dc7b537616ad9014b0878dcffbf25e4bdadfb198c0436c55e"
RECONSTRUCTOR_SHA256 = "5b12083b7a37d2ff22ddbfe8a0acfe207245ff895fb0f2a5093343a02a539b27"
RECONSTRUCTOR_SUPPORT_SHA256 = (
    "92d55a6296bb3288e9d6824008ac4a065a8b6a3ce10a177202de5a9c687b0b72")
RAW_CEILING = (
    "SELECTED_COMPONENT_L4_H6_M001__FULL_REAL_RESPONSE_"
    "RECONSTRUCTIBLE_ONLY_AFTER_ALL_300_SAME_STATE_RAYS__NO_"
    "RANGE_QUOTIENT_OR_WARD_CLAIM"
)
CHECKS: list[str] = []
NEGATIVE_FIXTURES = 0


class Failure(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise Failure(message)
    CHECKS.append(message)


def negative_fixture(message: str) -> None:
    global NEGATIVE_FIXTURES
    NEGATIVE_FIXTURES += 1
    check(True, message)


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def strict_json(path: Path) -> Any:
    try:
        return admission.strict_json(path)
    except admission.Failure as error:
        raise Failure(str(error)) from error


def parse_dependencies() -> dict[str, str]:
    try:
        typed = admission.validate_dependency_contract(
            DEPENDENCIES, DEPENDENCY_CONTRACT, ROOT, HERE)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(len(typed["entries"]) == 76, "exact 76 dependency rows")
    check(typed["role_counts"] == {
        "runtime_input": 4, "audited_reduction_gate": 6,
        "audited_acoustic_gate": 7, "audited_orbit_gate": 7,
        "retained_evidence": 52,
    }, "exact typed dependency role census")
    check(typed["upstream_chain_current"] is True and
          typed["launch_authorized"] is False,
          "audited upstream chain current; own audit still locks launch")
    fg = strict_json(
        ROOT / "AUDIT_G_GL6FG_SIGNED_CHARACTER_ACOUSTIC_REDUCTION_GATE_V003"
        / "INDEPENDENT_RESULT.json")
    check(type(fg) is dict, "GL6FG V003 audit object")
    admission.deep_exact(
        fg.get("schema"), "GL6FG_V003_INDEPENDENT_HOSTILE_AUDIT_RESULT_V001",
        "GL6FG V003 audit schema")
    check(admission.exact_int(fg.get("material_defect_count"),
                              "GL6FG material defect count", 0, 0) == 0 and
          admission.exact_int(fg.get("physics_defect_count"),
                              "GL6FG physics defect count", 0, 0) == 0 and
          admission.exact_string(fg.get("disposition"),
                                 "GL6FG disposition").startswith("PASS__"),
          "GL6FG V003 independently audited and load bearing")
    fh = strict_json(
        ROOT / "AUDIT_G_GL6FH_PARENT_DERIVED_ACOUSTIC_BUNDLE_GATE_V004"
        / "INDEPENDENT_RESULT.json")
    check(type(fh) is dict and
          fh.get("schema") == "GL6FH_V004_INDEPENDENT_HOSTILE_AUDIT_EXECUTION_V001" and
          fh.get("target") ==
          "DEVELOPMENT_G_GL6FH_PARENT_DERIVED_ACOUSTIC_BUNDLE_GATE_V004" and
          fh.get("target_modified_by_auditor") is False and
          fh.get("disposition") == "PASS" and
          fh.get("material_admission_defects") == 0 and
          fh.get("physics_defects") == 0,
          "GL6FH V004 independently audited and load bearing")
    fi = strict_json(
        ROOT / "AUDIT_G_GL6FI_MINIMAL_OFF_NYQUIST_WARD_ORBIT_GATE_V003"
        / "INDEPENDENT_RESULT.json")
    check(type(fi) is dict and
          fi.get("schema") == "GL6FI_V003_INDEPENDENT_HOSTILE_AUDIT_RESULT_V001" and
          fi.get("target") ==
          "DEVELOPMENT_G_GL6FI_MINIMAL_OFF_NYQUIST_WARD_ORBIT_GATE_V003" and
          fi.get("target_modified_by_auditor") is False and
          fi.get("disposition") == "PASS" and
          fi.get("material_defect_count") == 0 and
          fi.get("physics_defect_count") == 0,
          "GL6FI V003 independently audited and load bearing")
    return {path: row["sha256"] for path, row in typed["entries"].items()}


def load_postprocessor():
    specification = importlib.util.spec_from_file_location(
        "gl6fj_postprocessor_author_test", POST)
    check(specification is not None and specification.loader is not None,
          "postprocessor import specification")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def verify_design_and_surfaces() -> None:
    plan = strict_json(PLAN)
    expected_plan = {
        "schema": "GL6FJ_M001_FULL_CU_C64_RESPONSE_PLAN_V003",
        "status": "LAUNCH_BINDING__UPSTREAM_CHAIN_AUDITED__TARGET_AUDIT_REQUIRED",
        "dependency_state": {
            "gl6fg_v003": "INDEPENDENTLY_AUDITED_LOAD_BEARING",
            "gl6fh_v004": "INDEPENDENTLY_AUDITED_LOAD_BEARING",
            "gl6fi_v003": "INDEPENDENTLY_AUDITED_LOAD_BEARING",
            "upstream_chain_current": True,
            "target_independent_audit": "OPEN",
            "launch_authorized": False,
        },
        "physics_domain": {
            "L": 4, "volume": 64, "character": [0, 0, 1],
            "character_label": "m001", "self_conjugate": False,
            "real_source_channels": 24, "symmetric_entries": 300,
            "canonical_scalar_runs": 300, "redundant_scalar_runs": 0,
            "total_scalar_runs": 300,
        },
        "schedule": {
            "population": 1024, "burn": 2560, "windows": 6,
            "window_length": 256, "checkpoints": [64, 128, 256],
            "stay_units": 4096, "writer_units": 64,
            "random_seed": 96511001, "population_or_seed_ladder": False,
        },
        "symmetry_controls": {
            "signed_minus_lift": [0, 0, -1],
            "signed_conjugate_profiles_checked": True,
            "translation_forced_anomalous_expectation_zero": True,
            "finite_run_anomalous_zero_imposed": False,
            "integer_lift_centering_phases_retained": True,
            "real_source_chart_unaliased": True,
        },
        "owners": [
            "stay_contact", "writer_contact", "stay_spectral",
            "writer_spectral", "stay_writer_interference", "disconnected",
        ],
        "output": {
            "result_packet": RESULT_NAME,
            "raw_files": 300, "result_payload_files": 312,
            "write_once_staging": True, "atomic_no_replace_promotion": True,
        },
        "claim_ceiling": {
            "full_m001_response_only": True,
            "m123_composition_performed": False,
            "range_inverse_performed": False,
            "physical_ward_rank_claimed": False,
            "gravity_claimed": False,
        },
    }
    try:
        admission.deep_exact(plan, expected_plan, "exact measurement plan")
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(True, "exact measurement plan")

    source = SOURCE.read_text()
    check(digest(SOURCE) == SOURCE_SHA256,
          "m001 C++ engine bytes preserved from V002")
    for token in (
        "constexpr gl6er::Character character() { return {{0, 0, 1}}; }",
        "GL6FJ_SCALAR_PROJECTION_RAW_V001", "canonical_full_response",
        "translation_forced_anomalous_expectation_zero",
        "finite_run_anomalous_zero_imposed", "real_source_chart_unaliased",
        "GL6FJ_V001_SUBORDINATE_300_RAY_RAW_NOT_PRODUCTION",
        "stay_writer_interference", "print_scalar_checkpoint",
    ):
        check(token in source, "source required surface: " + token)
    for forbidden in ("priority-mixed", "m123", "m002", "graviton"):
        check(forbidden not in source, "source excludes " + forbidden)

    post = POST.read_text()
    check(digest(POST) == RECONSTRUCTOR_SHA256,
          "m001 reconstruction bytes preserved from V002")
    check(digest(HERE / "admission_contract_v002.py") ==
          RECONSTRUCTOR_SUPPORT_SHA256,
          "byte-preserved reconstructor support module retained")
    for token in (
        "exact 300 canonical raw census", "all 300 projections share exact histories",
        "lawful_range_inverse", "lawful_optical_schur", "normal_anomalous",
        "real_ao_blocks", '"certified_ward_rank": None',
        "finite_run_anomalous_zero_imposed",
    ):
        check(token in post, "postprocessor required surface: " + token)
    for forbidden in ("m123", "m002", "mixed_directions", "612"):
        check(forbidden not in post, "postprocessor excludes " + forbidden)

    controller = CONTROLLER.read_text()
    for token in (
        "AUTHORIZE_GL6FJ_V003_SINGLE_300_RUN_C64_M001_PILOT",
        "VERIFY_GL6FJ_V003_STAGING_BEFORE_ATOMIC_PROMOTION",
        "os.O_EXCL", "acquire_output_lock", "rename_no_replace",
        "target/audit/dependencies unchanged during production",
        "final target/audit/dependency reauthentication",
        "closed in-memory 300-raw census", "redundant_probe_files",
        "population_or_seed_ladder_run", "UPSTREAM_DEPENDENCY_CHAIN_CURRENT = True",
        "GL6FJ upstream dependency chain must be current",
    ):
        check(token in controller, "controller custody surface: " + token)
    check(controller.index("check(args.authorization == TOKEN") <
          controller.index("acquire_output_lock(bindings)"),
          "invalid launch authorization precedes output lock")

    result_verifier = RESULT_VERIFIER.read_text()
    for token in (
        "exact 312-file result payload census",
        "exact 300 raw-file name and manifest census",
        "full response reconstruction replay", "raw semantic digest",
        'observed_response["certified_ward_rank"] is None',
        "gl6fj_m001_runner_frozen.bin", "admission_contract_v003.py",
        "admission_contract_v002.py", "reconstructor_support_sha256",
        "FROZEN_DEPENDENCY_CONTRACT.json",
    ):
        check(token in result_verifier,
              "result verifier required surface: " + token)


def fixture(module: Any) -> dict[int, dict[str, Any]]:
    history_payload = strict_json(HISTORY)
    checkpoints = history_payload["checkpoints"]
    np = __import__("numpy")
    rng = np.random.default_rng(96513001)
    stay = rng.normal(size=(24, 24)); stay = (stay + stay.T) / 20.0
    writer = rng.normal(size=(24, 24)); writer = (writer + writer.T) / 25.0
    stay_s = rng.normal(size=(24, 24)); stay_s = (stay_s + stay_s.T) / 30.0
    writer_s = rng.normal(size=(24, 24)); writer_s = (writer_s + writer_s.T) / 35.0
    interference = rng.normal(size=(24, 24)); interference = (interference + interference.T) / 40.0
    stay_first = rng.normal(size=24) / 10.0
    writer_first = rng.normal(size=24) / 12.0
    first = stay_first + writer_first
    base_owners = (stay, writer, stay_s, writer_s, interference)
    rows: dict[int, dict[str, Any]] = {}
    for ray in range(300):
        direction = module.canonical_direction(ray)
        left, right = module.ray_pair(ray)
        raw_checkpoints = []
        for history_row in checkpoints:
            factor = (1.0 + history_row["window"] / 10.0 +
                      history_row["path_steps"] / 1000.0)
            first_values = [float(direction @ stay_first),
                            float(direction @ writer_first),
                            float(direction @ first)]
            owner_values = [float(direction @ (factor * matrix) @ direction)
                            for matrix in base_owners[:5]]
            owner_values.append(-first_values[2] ** 2)
            total = sum(owner_values)
            raw_connected = (history_row["lambda_geometric"] * total /
                             history_row["path_steps"])
            raw_checkpoints.append({
                **history_row,
                "stay_score_mean": first_values[0],
                "writer_score_mean": first_values[1],
                "first_mean": first_values[2],
                "first_score_rate": first_values[2] / history_row["path_steps"],
                "stay_contact": owner_values[0], "writer_contact": owner_values[1],
                "stay_spectral": owner_values[2], "writer_spectral": owner_values[3],
                "stay_writer_interference": owner_values[4],
                "disconnected": owner_values[5], "owner_total": total,
                "log_hessian_rate": total / history_row["path_steps"],
                "raw_connected_perron_hessian": raw_connected,
                "density_response": raw_connected / module.DENOMINATOR,
            })
        rows[ray] = {
            "schema": "GL6FJ_SCALAR_PROJECTION_RAW_V001",
            "mode": "subordinate_raw", "L": 4, "volume": 64,
            "component": "GL6CC_SELECTED_DENSE_H6_COMPONENT",
            "action_order": "through_h6", "momentum_label": "m001",
            "character": [0, 0, 1], "channels": 24, "symmetric_entries": 300,
            "design_kind": "canonical_full_response", "ray_index": ray,
            "ray_left": left, "ray_right": right,
            "ray_kind": "unit" if left == right else "normalized_pair_sum",
            "population": 1024, "burn": 2560, "window_length": 256,
            "windows": 6, "stay_units": 4096, "writer_units": 64,
            "stay_weight": 64, "random_seed": 96511001,
            "support_cache_dense_entries": 4008960,
            "burn_stays": history_payload["burn_stays"],
            "burn_writers": history_payload["burn_writers"],
            "burn_digest": history_payload["burn_digest"],
            "checkpoints": raw_checkpoints,
            "support_cache_hits": 1, "support_cache_misses": 0,
            "writer_cache_hits": 1, "writer_cache_misses": 0,
            "sampled_cold_residual": 0.0, "production_execution": False,
            "physical_ward_rank_claimed": False,
            "translation_forced_anomalous_expectation_zero": True,
            "finite_run_anomalous_zero_imposed": False,
            "real_source_chart_dimension": 24,
            "real_source_chart_unaliased": True, "ceiling": RAW_CEILING,
        }
    return rows


def reject_mutation(rows: dict[int, dict[str, Any]],
                    mutate: Callable[[dict[int, dict[str, Any]]], None],
                    module: Any, label: str) -> None:
    candidate = copy.deepcopy(rows)
    mutate(candidate)
    try:
        module.validate_common_raw(candidate)
    except module.Failure:
        negative_fixture(label)
        return
    raise Failure("negative fixture accepted: " + label)


def semantic_tests() -> None:
    module = load_postprocessor()
    result = module.selftest()
    check(result["result"] == "PASS" and result["ray_count"] == 300,
          "postprocessor analytic selftest")
    rows = fixture(module)
    module.validate_common_raw(rows)
    check(True, "synthetic complete 300-row semantic fixture accepted")
    reconstructed = module.reconstruct(rows)
    check(len(reconstructed["checkpoints"]) == 18 and
          reconstructed["certified_ward_rank"] is None and
          reconstructed["physical_ward_null_claimed"] is False and
          reconstructed["design"]["finite_run_anomalous_zero_imposed"] is False,
          "synthetic response reconstructs without Ward/anomaly imposition")

    reject_mutation(rows, lambda r: r[0].__setitem__("extra", 1), module,
                    "extra raw field rejected")
    reject_mutation(rows, lambda r: r[0].pop("ceiling"), module,
                    "missing raw field rejected")
    reject_mutation(rows, lambda r: r[0].__setitem__("ray_index", True), module,
                    "Boolean ray index rejected")
    reject_mutation(rows, lambda r: r[0].__setitem__("character", [0, 0, 2]), module,
                    "wrong character rejected")
    reject_mutation(rows, lambda r: r[0].__setitem__("real_source_chart_unaliased", False), module,
                    "aliased chart declaration rejected")
    reject_mutation(rows, lambda r: r[0].__setitem__("finite_run_anomalous_zero_imposed", True), module,
                    "forced finite anomaly rejected")
    reject_mutation(rows, lambda r: r[0]["checkpoints"][0].__setitem__("stay_contact", 999.0), module,
                    "owner-once corruption rejected")
    reject_mutation(rows, lambda r: r[0]["checkpoints"][0].__setitem__("decision_digest", 0), module,
                    "frozen history drift rejected")
    reject_mutation(rows, lambda r: r[0].__setitem__("support_cache_hits", -1), module,
                    "negative cache counter rejected")
    reject_mutation(rows, lambda r: r[1].__setitem__("support_cache_hits", 2), module,
                    "cross-ray same-history cache mismatch rejected")
    reject_mutation(rows, lambda r: r[0].__setitem__(
        "support_cache_hits", admission.U64_MAX + 1), module,
                    "uint64 cache overflow rejected")
    reject_mutation(rows, lambda r: r[0]["checkpoints"][0].__setitem__(
        "stay_contact", True), module,
                    "Boolean physical number rejected")
    reject_mutation(rows, lambda r: r[0]["checkpoints"][0].__setitem__(
        "stay_contact", math.inf), module,
                    "nonfinite nested physical number rejected")
    reject_mutation(rows, lambda r: r[0]["checkpoints"][0].__setitem__(
        "decision_digest", admission.U64_MAX + 1), module,
                    "uint64 nested history overflow rejected")
    missing = dict(rows); missing.pop(299)
    try:
        module.validate_common_raw(missing)
    except module.Failure:
        negative_fixture("missing ray rejected")
    else:
        raise Failure("missing ray accepted")


def expect_admission_failure(call: Callable[[], Any], label: str) -> None:
    try:
        call()
    except (admission.Failure, Failure, FileNotFoundError, OSError,
            ValueError):
        negative_fixture(label)
        return
    raise Failure("admission negative fixture accepted: " + label)


def write_packet(directory: Path) -> None:
    directory.mkdir()
    nested = directory / "payload"
    nested.mkdir()
    data = nested / "row.txt"
    data.write_text("bounded fixture\n")
    manifest = directory / "MANIFEST.sha256"
    manifest.write_text(f"{digest(data)}  payload/row.txt\n")
    (directory / "SEAL.sha256").write_text(
        f"{digest(manifest)}  MANIFEST.sha256\n")


def admission_negative_tests() -> None:
    for payload, label in (
            (b'{"x":1,"x":2}', "duplicate JSON key rejected recursively"),
            (b'{"x":NaN}', "NaN JSON constant rejected recursively"),
            (b'{"x":1e9999}', "floating overflow rejected recursively"),
            (b'{"x":18446744073709551616}',
             "positive integer overflow rejected recursively"),
            (b'{"x":-9223372036854775809}',
             "negative integer overflow rejected recursively")):
        expect_admission_failure(
            lambda payload=payload: admission.strict_json_bytes(payload), label)
    expect_admission_failure(
        lambda: admission.deep_exact({"x": [True]}, {"x": [1]},
                                     "nested scalar"),
        "recursive Boolean/integer coercion rejected")
    expect_admission_failure(
        lambda: admission.deep_exact({"x": [(1,)]}, {"x": [(1,)]},
                                     "nested scalar"),
        "recursive non-JSON container rejected")
    expect_admission_failure(
        lambda: admission.deep_exact(
            {"x": [admission.U64_MAX + 1]},
            {"x": [admission.U64_MAX + 1]}, "nested scalar"),
        "recursive equal integer overflow rejected")

    contract = strict_json(DEPENDENCY_CONTRACT)
    admission.validate_dependency_contract(
        DEPENDENCIES, DEPENDENCY_CONTRACT, ROOT, HERE)
    check(True, "typed dependency contract baseline accepted")
    with tempfile.TemporaryDirectory(prefix="gl6fj_v003_admission_") as name:
        temporary = Path(name)
        for index, mutate in enumerate((
                lambda row: row.__setitem__("entry_count", True),
                lambda row: row.__setitem__("extra", 1),
                lambda row: row["roles"]["audited_reduction_gate"].__setitem__(
                    "promotion_state", "AUTHENTICATED_INPUT"),
                lambda row: (
                    row["roles"]["runtime_input"]["paths"].__setitem__(
                        0, row["roles"]["audited_orbit_gate"]["paths"][0]),
                    row["roles"]["audited_orbit_gate"]["paths"].__setitem__(
                        0, contract["roles"]["runtime_input"]["paths"][0])),
                lambda row: row["roles"]["runtime_input"]["paths"].append(
                    row["roles"]["audited_orbit_gate"]["paths"][0]))):
            candidate = copy.deepcopy(contract)
            mutate(candidate)
            path = temporary / f"dependency_contract_{index}.json"
            path.write_text(json.dumps(candidate, sort_keys=True,
                                       allow_nan=False))
            expect_admission_failure(
                lambda path=path: admission.validate_dependency_contract(
                    DEPENDENCIES, path, ROOT, HERE),
                "typed dependency contract mutation rejected")
        changed_ledger = temporary / "DEPENDENCIES.sha256"
        changed_ledger.write_bytes(DEPENDENCIES.read_bytes() + b"\n")
        expect_admission_failure(
            lambda: admission.validate_dependency_contract(
                changed_ledger, DEPENDENCY_CONTRACT, ROOT, HERE),
            "dependency ledger substitution rejected by exact digest")
        coupled = copy.deepcopy(contract)
        coupled["ledger_sha256"] = digest(changed_ledger)
        coupled_path = temporary / "coupled_dependency_contract.json"
        coupled_path.write_text(json.dumps(coupled, sort_keys=True,
                                           allow_nan=False))
        expect_admission_failure(
            lambda: admission.validate_dependency_contract(
                changed_ledger, coupled_path, ROOT, HERE),
            "coupled dependency ledger and contract substitution rejected")

        repository = temporary / "dependency_repository"
        repository.mkdir()
        packet = repository / "packet"
        packet.mkdir()
        real_directory = repository / "real_dependency"
        real_directory.mkdir()
        (real_directory / "row.bin").write_bytes(b"bounded")
        os.symlink("real_dependency", repository / "aliased_dependency")
        expect_admission_failure(
            lambda: admission.ordinary_dependency_path(
                repository, packet, "../aliased_dependency/row.bin"),
            "intermediate dependency-directory symlink rejected")

        valid = temporary / "valid_packet"
        write_packet(valid)
        admission.verify_closed_packet(valid)
        check(True, "closed-tree baseline accepted")
        for label, mutate in (
                ("nested root-control name rejected",
                 lambda row: (row / "payload" / "MANIFEST.sha256").write_text("x")),
                ("packet symlink rejected",
                 lambda row: os.symlink("payload/row.txt", row / "linked")),
                ("packet FIFO rejected",
                 lambda row: os.mkfifo(row / "fifo")),
                ("empty packet directory rejected",
                 lambda row: (row / "empty").mkdir())):
            candidate = temporary / ("bad_" + str(NEGATIVE_FIXTURES))
            shutil.copytree(valid, candidate)
            mutate(candidate)
            expect_admission_failure(
                lambda candidate=candidate: admission.verify_closed_packet(candidate),
                label)

    manifest_hash = "a" * 64
    seal_hash = "b" * 64
    audit = {
        "schema": admission.AUDIT_SCHEMA,
        "target": admission.AUDIT_TARGET,
        "target_modified_by_auditor": False,
        "disposition": admission.AUDIT_DISPOSITION,
        "independent_checks": 1,
        "target_manifest_sha256": manifest_hash,
        "target_seal_file_sha256": seal_hash,
        "material_defects": [], "physics_defects": [],
        "production_run": False,
        "promotion": {
            "target_admission_pass": True,
            "production_result_promoted": False,
            "physical_ward_rank_claimed": False,
            "gravity_claimed": False,
        },
        "claim_ceiling": admission.AUDIT_CEILING,
    }
    admission.validate_audit_result(audit, manifest_hash, seal_hash)
    check(True, "closed audit-result baseline accepted")
    for label, mutate in (
            ("audit extra field rejected", lambda row: row.__setitem__("extra", 1)),
            ("audit promotion extra field rejected",
             lambda row: row["promotion"].__setitem__("extra", False)),
            ("audit integer Boolean rejected",
             lambda row: row.__setitem__("independent_checks", True)),
            ("audit defect type substitution rejected",
             lambda row: row.__setitem__("material_defects", {})),
            ("audit admission Boolean type substitution rejected",
             lambda row: row["promotion"].__setitem__(
                 "target_admission_pass", 1)),
            ("audit production promotion rejected",
             lambda row: row["promotion"].__setitem__(
                 "production_result_promoted", True)),
            ("audit gravity promotion rejected",
             lambda row: row["promotion"].__setitem__(
                 "gravity_claimed", True))):
        candidate = copy.deepcopy(audit)
        mutate(candidate)
        expect_admission_failure(
            lambda candidate=candidate: admission.validate_audit_result(
                candidate, manifest_hash, seal_hash), label)


def build_and_smoke() -> None:
    modes = {
        "ordinary": ["-O0"],
        "optimized": ["-O3", "-DNDEBUG"],
        "ubsan": ["-O1", "-fsanitize=undefined", "-fno-sanitize-recover=undefined"],
    }
    histories = []
    with tempfile.TemporaryDirectory(prefix="gl6fj_author_verify_") as name:
        temporary = Path(name)
        for label, flags in modes.items():
            executable = temporary / ("runner_" + label)
            command = ["clang++", "-std=c++20", *flags, "-Wall", "-Wextra",
                       "-pedantic", str(SOURCE), "-o", str(executable)]
            build = subprocess.run(command, cwd=ROOT, text=True,
                                   capture_output=True, timeout=600)
            check(build.returncode == 0, label + " build")
            for selftest in ("catalog-selftest", "character-selftest"):
                process = subprocess.run([str(executable), selftest], cwd=ROOT,
                                         text=True, capture_output=True, timeout=120)
                payload = (admission.strict_json_bytes(
                    process.stdout.encode("utf-8"))
                    if process.returncode == 0 else {})
                check(process.returncode == 0 and type(payload) is dict and
                      payload.get("result") == "PASS",
                      label + " " + selftest)
            for ray in (0, 1, 299):
                process = subprocess.run([str(executable), "smoke-projection", str(ray)],
                                         cwd=ROOT, text=True, capture_output=True,
                                         timeout=300)
                check(process.returncode == 0, label + f" smoke ray {ray}")
                raw = admission.strict_json_bytes(
                    process.stdout.encode("utf-8"))
                check(raw["schema"] == "GL6FJ_SCALAR_PROJECTION_RAW_V001" and
                      raw["momentum_label"] == "m001" and
                      raw["character"] == [0, 0, 1] and raw["ray_index"] == ray and
                      raw["production_execution"] is False and
                      raw["finite_run_anomalous_zero_imposed"] is False,
                      label + f" smoke semantics ray {ray}")
                for checkpoint in raw["checkpoints"]:
                    owner_sum = sum(checkpoint[key] for key in (
                        "stay_contact", "writer_contact", "stay_spectral",
                        "writer_spectral", "stay_writer_interference", "disconnected"))
                    check(abs(owner_sum - checkpoint["owner_total"]) <=
                          3e-11 * max(abs(checkpoint["owner_total"]), 1.0),
                          label + " smoke owner-once identity")
                histories.append((label, ray, raw["burn_digest"], tuple(
                    row["decision_digest"] for row in raw["checkpoints"])))
    for label in modes:
        values = {(burn, checkpoints) for mode, _, burn, checkpoints in histories
                  if mode == label}
        check(len(values) == 1, label + " source-zero smoke histories identical")
    cross = {(burn, checkpoints) for _, _, burn, checkpoints in histories}
    check(len(cross) == 1, "ordinary optimized UBSan histories identical")


def controller_negative_tests() -> None:
    result = ROOT / RESULT_NAME
    staging = ROOT / (RESULT_NAME + ".STAGING")
    lock = ROOT / ("." + RESULT_NAME + ".LOCK")
    check(not result.exists() and not staging.exists() and not lock.exists(),
          "GL6FJ output/lock paths absent before controller tests")
    bad = subprocess.run([
        sys.executable, "-B", str(CONTROLLER), "launch", "--workers", "1",
        "--authorization", "INVALID",
    ], cwd=ROOT, text=True, capture_output=True, timeout=60)
    check(bad.returncode != 0 and "exact production authorization token" in bad.stderr,
          "invalid authorization fails before launch")
    check(not result.exists() and not staging.exists() and not lock.exists(),
          "invalid authorization leaves no output or stale lock")
    preflight = subprocess.run([sys.executable, "-B", str(CONTROLLER), "preflight"],
                               cwd=ROOT, text=True, capture_output=True, timeout=60)
    preflight_payload = (admission.strict_json_bytes(
        preflight.stdout.encode("utf-8"))
        if preflight.returncode == 0 else {})
    check(preflight.returncode == 0 and
          type(preflight_payload) is dict and
          preflight_payload.get("production_run") is False,
          "preflight never launches production")
    check(not result.exists() and not staging.exists() and not lock.exists(),
          "preflight leaves no output or lock")


def main() -> int:
    parse_dependencies()
    verify_design_and_surfaces()
    semantic_tests()
    admission_negative_tests()
    build_and_smoke()
    controller_negative_tests()
    check(os.access(RESULT_VERIFIER, os.R_OK), "result verifier readable")
    print(json.dumps({
        "schema": "GL6FJ_AUTHOR_VERIFICATION_V003",
        "checks": len(CHECKS),
        "result": f"PASS__{len(CHECKS)}/{len(CHECKS)}",
        "builds": ["ordinary", "optimized", "ubsan"],
        "negative_fixtures": NEGATIVE_FIXTURES,
        "production_run": False,
        "physical_ward_rank_claimed": False,
        "dependency_chain_current": True,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (Failure, admission.Failure, FileNotFoundError, OSError, ValueError,
            subprocess.SubprocessError) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
