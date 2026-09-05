#!/usr/bin/env python3
"""Standalone fail-closed V003 verifier copied into a GL6FJ result packet."""

from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath
import sys
from typing import Any

import admission_contract_v003 as admission

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET_NAME = "DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
AUDIT_NAME = "AUDIT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
RESULT_NAME = "RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
TARGET = ROOT / TARGET_NAME
AUDIT = ROOT / AUDIT_NAME
AUDIT_RESULT = AUDIT / "INDEPENDENT_RESULT.json"
AUDIT_SCHEMA = admission.AUDIT_SCHEMA
AUDIT_DISPOSITION = admission.AUDIT_DISPOSITION
LOCAL_RECONSTRUCTOR = HERE / "reconstruct_m001_full_response_v002.py"
LOCAL_HISTORY = HERE / "EXPECTED_HISTORY.json"
LOCAL_DEPENDENCIES = HERE / "FROZEN_DEPENDENCIES.sha256"
LOCAL_DEPENDENCY_CONTRACT = HERE / "FROZEN_DEPENDENCY_CONTRACT.json"
LOCAL_ADMISSION_CONTRACT = HERE / "admission_contract_v003.py"
LOCAL_RECONSTRUCTOR_SUPPORT = HERE / "admission_contract_v002.py"
LOCAL_EXECUTABLE = HERE / "gl6fj_m001_runner_frozen.bin"
STAGING_TOKEN = "VERIFY_GL6FJ_V003_STAGING_BEFORE_ATOMIC_PROMOTION"
HEX = frozenset("0123456789abcdef")
CHECKS: list[str] = []

PIN_KEYS = {
    "target_manifest_sha256", "target_seal_file_sha256",
    "audit_manifest_sha256", "audit_seal_file_sha256",
    "audit_result_sha256", "source_sha256", "reconstructor_sha256",
    "expected_history_sha256", "result_verifier_sha256",
    "controller_sha256", "measurement_plan_sha256", "dependencies_sha256",
    "dependency_contract_sha256", "admission_contract_sha256",
    "reconstructor_support_sha256",
}
CUSTODY_KEYS = {
    "schema", *PIN_KEYS, "compiler_path", "compiler_sha256",
    "compiler_version", "build_command_without_temporary_output",
    "executable_sha256", "workers", "started_unix_ns",
    "completed_unix_ns", "elapsed_seconds", "raw_files",
    "canonical_files", "redundant_probe_files", "character",
    "raw_semantic_sha256", "raw_manifest_sha256",
    "raw_seal_file_sha256", "postrun_authorization_reauthenticated",
    "output_lock_retained_through_promotion", "promotion_primitive",
    "production_run", "population_or_seed_ladder_run",
    "translation_forced_anomalous_expectation_zero",
    "finite_run_anomalous_zero_imposed",
    "physical_ward_rank_claimed",
}


class Failure(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise Failure(message)
    CHECKS.append(message)


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise Failure("duplicate JSON key: " + key)
        result[key] = value
    return result


def strict_json_bytes(payload: bytes) -> Any:
    try:
        return admission.strict_json_bytes(payload)
    except admission.Failure as error:
        raise Failure(str(error)) from error


def strict_json(path: Path) -> Any:
    return strict_json_bytes(path.read_bytes())


def _safe_packet_name(name: str) -> bool:
    pure = PurePosixPath(name)
    return (name != "" and not pure.is_absolute() and
            all(part not in {"", ".", ".."} for part in pure.parts))


def hashes(path: Path, *, packet_names: bool = True) -> dict[str, str]:
    try:
        return admission.parse_hashes(path, packet_names=packet_names)
    except admission.Failure as error:
        raise Failure(str(error)) from error


def verify_closed_packet(directory: Path) -> tuple[str, str, dict[str, str]]:
    try:
        return admission.verify_closed_packet(directory)
    except admission.Failure as error:
        raise Failure(str(error)) from error


def verify_external_dependencies(ledger: Path) -> str:
    try:
        admission.validate_dependency_contract(
            ledger, LOCAL_DEPENDENCY_CONTRACT, ROOT, HERE)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    return digest(ledger)


def load_reconstructor():
    specification = importlib.util.spec_from_file_location(
        "gl6fj_m001_frozen_result_reconstructor", LOCAL_RECONSTRUCTOR)
    check(specification is not None and specification.loader is not None,
          "local reconstructor import specification")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def verify_invocation() -> None:
    if HERE.name == RESULT_NAME + ".STAGING":
        check(sys.argv == [sys.argv[0], STAGING_TOKEN],
              "exact staging verification token and arity")
    else:
        check(HERE.name == RESULT_NAME, "exact final result directory name")
        check(sys.argv == [sys.argv[0]], "exact final verifier arity")


def expected_raw_names() -> set[str]:
    return {f"RAW/m001_ray_{ray:03d}.json" for ray in range(300)}


def expected_result_payload_names() -> set[str]:
    return expected_raw_names() | {
        "FULL_RESPONSE_RECONSTRUCTION.json",
        "RAW_MANIFEST.sha256",
        "RAW_SEAL.sha256",
        "LAUNCH_CUSTODY.json",
        "reconstruct_m001_full_response_v002.py",
        "EXPECTED_HISTORY.json",
        "FROZEN_DEPENDENCIES.sha256",
        "FROZEN_DEPENDENCY_CONTRACT.json",
        "admission_contract_v003.py",
        "admission_contract_v002.py",
        "gl6fj_m001_runner_frozen.bin",
        "verify_result_packet.py",
    }


def verify_audit(target_manifest_hash: str,
                 target_seal_hash: str) -> tuple[str, str, str]:
    audit_manifest_hash, audit_seal_hash, audit_manifest = \
        verify_closed_packet(AUDIT)
    payload = strict_json(AUDIT_RESULT)
    try:
        admission.validate_audit_result(
            payload, target_manifest_hash, target_seal_hash)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    audit_result_hash = digest(AUDIT_RESULT)
    check(audit_manifest.get("INDEPENDENT_RESULT.json") == audit_result_hash,
          "audit manifest binds audit result")
    return audit_manifest_hash, audit_seal_hash, audit_result_hash


def main() -> int:
    verify_invocation()
    manifest_hash, seal_hash, manifest = verify_closed_packet(HERE)
    check(set(manifest) == expected_result_payload_names(),
          "exact 312-file result payload census")

    raw_manifest = hashes(HERE / "RAW_MANIFEST.sha256")
    check(set(raw_manifest) == expected_raw_names(),
          "exact 300 raw-file name and manifest census")
    check(hashes(HERE / "RAW_SEAL.sha256") == {
        "RAW_MANIFEST.sha256": digest(HERE / "RAW_MANIFEST.sha256")},
        "raw seal binds raw manifest")
    check(set(raw_manifest) == {name for name in manifest
                                if name.startswith("RAW/")},
          "raw manifest equals recursive result census")
    for name, expected in raw_manifest.items():
        check(digest(HERE / name) == expected, "raw payload hash: " + name)

    custody = strict_json(HERE / "LAUNCH_CUSTODY.json")
    check(type(custody) is dict and set(custody) == CUSTODY_KEYS,
          "exact launch-custody schema keys")
    check(custody["schema"] == "GL6FJ_PRODUCTION_LAUNCH_CUSTODY_V003",
          "launch-custody schema")
    for key in PIN_KEYS | {"compiler_sha256", "executable_sha256",
                           "raw_semantic_sha256", "raw_manifest_sha256",
                           "raw_seal_file_sha256"}:
        check(type(custody[key]) is str and len(custody[key]) == 64 and
              set(custody[key]) <= HEX, "custody SHA-256 field: " + key)
    for key in ("workers", "started_unix_ns", "completed_unix_ns",
                "raw_files", "canonical_files", "redundant_probe_files"):
        admission.exact_u64(custody[key], "custody uint64: " + key)
    check(type(custody["elapsed_seconds"]) is float and
          math.isfinite(custody["elapsed_seconds"]),
          "custody exact finite floating elapsed time")
    check(custody["workers"] in (1, 2, 3, 4) and
          custody["raw_files"] == 300 and
          custody["canonical_files"] == 300 and
          custody["redundant_probe_files"] == 0 and
          type(custody["character"]) is list and
          len(custody["character"]) == 3 and
          all(type(value) is int for value in custody["character"]) and
          custody["character"] == [0, 0, 1],
          "custody exact launch census")
    check(custody["started_unix_ns"] > 0 and
          custody["completed_unix_ns"] > custody["started_unix_ns"] and
          custody["elapsed_seconds"] > 0,
          "custody ordered positive runtime")
    check(abs(custody["elapsed_seconds"] -
              (custody["completed_unix_ns"] - custody["started_unix_ns"])
              / 1_000_000_000.0) <= 1e-9,
          "custody elapsed-time identity")
    true_fields = (
        "postrun_authorization_reauthenticated",
        "output_lock_retained_through_promotion", "production_run",
        "translation_forced_anomalous_expectation_zero")
    check(all(custody[key] is True for key in true_fields),
          "custody required true controls")
    false_fields = ("population_or_seed_ladder_run",
                    "finite_run_anomalous_zero_imposed",
                    "physical_ward_rank_claimed")
    check(all(custody[key] is False for key in false_fields),
          "custody required false controls")
    check(custody["promotion_primitive"] ==
          "darwin_renamex_np_RENAME_EXCL", "atomic no-replace promotion")
    check(type(custody["compiler_path"]) is str and
          Path(custody["compiler_path"]).is_absolute() and
          type(custody["compiler_version"]) is str and
          custody["compiler_version"] != "",
          "compiler identity strings")
    expected_build = [
        custody["compiler_path"], "-std=c++20", "-O3", "-DNDEBUG",
        "-Wall", "-Wextra", "-pedantic",
        str(TARGET / "m001_full_cu_c64_projection_v001.cpp")]
    admission.deep_exact(custody["build_command_without_temporary_output"],
                         expected_build,
                         "exact compiler command excluding temporary output")

    target_manifest_hash, target_seal_hash, _ = verify_closed_packet(TARGET)
    check(custody["target_manifest_sha256"] == target_manifest_hash and
          custody["target_seal_file_sha256"] == target_seal_hash,
          "custody target packet binding")
    audit_manifest_hash, audit_seal_hash, audit_result_hash = verify_audit(
        target_manifest_hash, target_seal_hash)
    check(custody["audit_manifest_sha256"] == audit_manifest_hash and
          custody["audit_seal_file_sha256"] == audit_seal_hash and
          custody["audit_result_sha256"] == audit_result_hash,
          "custody audit packet binding")

    file_pins = {
        "source_sha256": TARGET / "m001_full_cu_c64_projection_v001.cpp",
        "reconstructor_sha256": TARGET / "reconstruct_m001_full_response_v002.py",
        "expected_history_sha256": TARGET / "EXPECTED_HISTORY.json",
        "result_verifier_sha256": TARGET / "result_packet_verifier_v003.py",
        "controller_sha256": TARGET / "run_gl6fj_v003_pilot.py",
        "measurement_plan_sha256": TARGET / "MEASUREMENT_PLAN.json",
        "dependencies_sha256": TARGET / "DEPENDENCIES.sha256",
        "dependency_contract_sha256": TARGET / "DEPENDENCY_CONTRACT.json",
        "admission_contract_sha256": TARGET / "admission_contract_v003.py",
        "reconstructor_support_sha256": TARGET / "admission_contract_v002.py",
    }
    for key, path in file_pins.items():
        check(custody[key] == digest(path), "custody target file pin: " + key)
    check(custody["result_verifier_sha256"] == digest(Path(__file__)) and
          custody["reconstructor_sha256"] == digest(LOCAL_RECONSTRUCTOR) and
          custody["expected_history_sha256"] == digest(LOCAL_HISTORY) and
          custody["dependencies_sha256"] == digest(LOCAL_DEPENDENCIES) and
          custody["dependency_contract_sha256"] ==
          digest(LOCAL_DEPENDENCY_CONTRACT) and
          custody["admission_contract_sha256"] ==
          digest(LOCAL_ADMISSION_CONTRACT) and
          custody["reconstructor_support_sha256"] ==
          digest(LOCAL_RECONSTRUCTOR_SUPPORT) and
          custody["executable_sha256"] == digest(LOCAL_EXECUTABLE),
          "copied verifier/reconstructor/history/contracts/dependencies/executable pins")
    check((LOCAL_DEPENDENCIES.read_bytes() ==
           (TARGET / "DEPENDENCIES.sha256").read_bytes()),
          "frozen and target dependency ledgers identical")
    check((LOCAL_DEPENDENCY_CONTRACT.read_bytes() ==
           (TARGET / "DEPENDENCY_CONTRACT.json").read_bytes()) and
          (LOCAL_ADMISSION_CONTRACT.read_bytes() ==
           (TARGET / "admission_contract_v003.py").read_bytes()) and
          (LOCAL_RECONSTRUCTOR_SUPPORT.read_bytes() ==
           (TARGET / "admission_contract_v002.py").read_bytes()),
          "frozen and target admission contracts identical")
    check(verify_external_dependencies(LOCAL_DEPENDENCIES) ==
          custody["dependencies_sha256"], "all external dependencies replay")
    # Compiler path/version/hash are historical launch custody.  The exact
    # executed binary is copied into this immutable packet and hash-checked;
    # a later operating-system compiler update must not invalidate the result.

    check(custody["raw_manifest_sha256"] ==
          digest(HERE / "RAW_MANIFEST.sha256") and
          custody["raw_seal_file_sha256"] == digest(HERE / "RAW_SEAL.sha256"),
          "custody raw seal binding")
    raw_rows: dict[int, dict[str, Any]] = {}
    semantic_rows: dict[str, dict[str, Any]] = {}
    for relative in sorted(raw_manifest):
        raw = strict_json(HERE / relative)
        check(type(raw) is dict, "raw JSON object: " + relative)
        semantic_rows[relative] = raw
        check(raw.get("schema") == "GL6FJ_SCALAR_PROJECTION_RAW_V001",
              "exact m001 raw schema")
        ray = raw.get("ray_index")
        check(type(ray) is int and ray not in raw_rows,
              "unique canonical raw key")
        raw_rows[ray] = raw
    semantic_digest = sha256(json.dumps(
        semantic_rows, sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode()).hexdigest()
    check(custody["raw_semantic_sha256"] == semantic_digest,
          "raw semantic digest")

    module = load_reconstructor()
    module.validate_common_raw(raw_rows)
    expected_response = module.reconstruct(raw_rows)
    observed_response = strict_json(HERE / "FULL_RESPONSE_RECONSTRUCTION.json")
    admission.deep_exact(observed_response, expected_response,
                         "full response reconstruction replay")
    check(observed_response["quotient"]["performed"] is False,
          "result does not invent a quotient")
    check(observed_response["certified_ward_rank"] is None and
          observed_response["physical_ward_null_claimed"] is False,
          "scalar response does not claim a physical Ward rank")
    print(json.dumps({
        "schema": "GL6FJ_RESULT_PACKET_VERIFICATION_V003",
        "checks": len(CHECKS),
        "result": f"PASS__{len(CHECKS)}/{len(CHECKS)}",
        "manifest_sha256": manifest_hash,
        "seal_file_sha256": seal_hash,
        "production_reexecuted": False,
        "physical_ward_rank_claimed": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (Failure, admission.Failure, FileNotFoundError, OSError,
            ValueError) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
