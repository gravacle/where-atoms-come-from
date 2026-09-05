#!/usr/bin/env python3
"""Fail-closed replay verifier for the sealed GL6FJ V003 m001 target."""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import admission_contract_v003 as admission

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PAYLOAD = {
    "AUTHOR_RESULT.json",
    "BUILD.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "DEPENDENCY_CONTRACT.json",
    "DESIGN.md",
    "EXACT_RESULT.json",
    "EXPECTED_HISTORY.json",
    "MEASUREMENT_PLAN.json",
    "README.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "VERIFICATION.txt",
    "admission_contract_v003.py",
    "admission_contract_v002.py",
    "m001_full_cu_c64_projection_v001.cpp",
    "reconstruct_m001_full_response_v002.py",
    "result_packet_verifier_v003.py",
    "run_gl6fj_v003_pilot.py",
    "verify_gl6fj_v003.py",
    "verify_packet.py",
}
CHECKS = 0


class Failure(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise Failure(message)


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def hashes(path: Path) -> dict[str, str]:
    try:
        return admission.parse_hashes(path)
    except admission.Failure as error:
        raise Failure(str(error)) from error


def run_author(optimized: bool) -> dict:
    command = [sys.executable] + (["-O"] if optimized else []) + [
        "-B", str(HERE / "verify_gl6fj_v003.py")]
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONPYCACHEPREFIX"] = "/private/tmp/gl6fj_packet_pyc"
    process = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                             env=environment)
    check(process.returncode == 0,
          "author replay exits zero: " + process.stderr[-2000:])
    try:
        payload = admission.strict_json_bytes(process.stdout.encode("utf-8"))
    except (json.JSONDecodeError, admission.Failure) as error:
        raise Failure("author replay JSON: " + str(error)) from error
    check(type(payload) is dict and set(payload) == {
        "schema", "checks", "result", "builds", "negative_fixtures",
        "production_run", "physical_ward_rank_claimed",
        "dependency_chain_current",
    }, "author replay exact closed schema")
    check(payload["schema"] == "GL6FJ_AUTHOR_VERIFICATION_V003" and
          type(payload["checks"]) is int and payload["checks"] > 0 and
          payload["result"] ==
          f"PASS__{payload['checks']}/{payload['checks']}" and
          payload["builds"] == ["ordinary", "optimized", "ubsan"] and
          payload["negative_fixtures"] == 42 and
          payload["production_run"] is False and
          payload["physical_ward_rank_claimed"] is False and
          payload["dependency_chain_current"] is True,
          "author replay exact typed result")
    return payload


def main() -> int:
    try:
        _manifest_hash, _seal_hash, manifest = \
            admission.verify_closed_packet(HERE)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(set(manifest) == PAYLOAD, "manifest exact payload census")
    try:
        typed = admission.validate_dependency_contract(
            HERE / "DEPENDENCIES.sha256", HERE / "DEPENDENCY_CONTRACT.json",
            ROOT, HERE)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(len(typed["entries"]) == 76 and typed["role_counts"] == {
        "runtime_input": 4, "audited_reduction_gate": 6,
        "audited_acoustic_gate": 7, "audited_orbit_gate": 7,
        "retained_evidence": 52,
    }, "exact typed 76-file dependency census")
    check(typed["upstream_chain_current"] is True and
          typed["launch_authorized"] is False,
          "audited upstream chain current; own audit hard locks launch")

    ordinary = run_author(False)
    optimized = run_author(True)
    try:
        admission.deep_exact(ordinary, optimized,
                             "ordinary/optimized author replay identity")
    except admission.Failure as error:
        raise Failure(str(error)) from error

    author = admission.strict_json(HERE / "AUTHOR_RESULT.json")
    expected_author = {
        "schema": "GL6FJ_PRELAUNCH_AUTHOR_RESULT_V003",
        "status": (
            "PASS__AUDITED_UPSTREAM_CHAIN__M001_LAUNCH_BINDING_"
            "TARGET__OWN_AUDIT_OPEN"),
        "checks": ordinary["checks"],
        "canonical_scalar_runs": 300,
        "redundant_scalar_runs": 0,
        "total_scalar_runs": 300,
        "negative_fixtures": 42,
        "production_run": False,
        "physical_ward_rank_claimed": False,
        "independent_audit": "OPEN",
        "ceiling": (
            "PROSPECTIVE_TARGET_ONLY__NO_PRODUCTION_RESPONSE_RANGE_WARD_"
            "EINSTEIN_GRAVITY_C_R_OR_G"),
    }
    try:
        admission.deep_exact(author, expected_author,
                             "exact author-result tree")
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(True, "author-result recursive exact schema/types/values")

    exact = admission.strict_json(HERE / "EXACT_RESULT.json")
    expected_exact = {
        "schema": "GL6FJ_PROSPECTIVE_EXACT_RESULT_V003",
        "status": (
            "UPSTREAM_CHAIN_AUDITED__AUTHOR_VERIFIED__TARGET_AUDIT_OPEN__"
            "LAUNCH_HARD_LOCKED"),
        "character": [0, 0, 1],
        "self_conjugate": False,
        "real_source_chart_dimension": 24,
        "real_source_chart_unaliased": True,
        "translation_forced_anomalous_expectation_zero": True,
        "finite_run_anomalous_zero_imposed": False,
        "six_owner_classes_retained": True,
        "same_state_history_required": True,
        "total_scalar_runs": 300,
        "dependency_chain_current": True,
        "independent_audit": "OPEN",
        "production_run": False,
        "response_measured": False,
        "range_inverse_performed": False,
        "ward_spectrum_computed": False,
        "physical_ward_rank_claimed": False,
        "gravity_claimed": False,
        "admission_repairs": [
            "RECURSIVE_STRICT_SCHEMA_TYPE_DEEP_EQUALITY_NONFINITE_AND_OVERFLOW",
            "TYPED_EXACT_DEPENDENCY_ROLE_PATH_HASH_CENSUS",
            "UINT64_CACHE_DOMAINS_AND_CROSS_RAY_SAME_HISTORY_CACHE_EQUALITY",
            "CLOSED_PACKET_TREE_ROOT_CONTROLS_SYMLINK_SPECIAL_AND_EMPTY_DIRECTORY_REJECTION",
            "EXACT_CLOSED_AUDIT_RESULT_SCHEMA_AND_PROMOTION_TYPING",
            "AUDITED_GL6FH_V004_AND_GL6FI_V003_TARGET_AUDIT_REBIND",
        ],
        "next_required_action": (
            "COMPLETE_DISTINCT_GL6FJ_V003_HOSTILE_AUDIT_BEFORE_ANY_"
            "AUTHORIZED_PRODUCTION_LAUNCH"),
    }
    try:
        admission.deep_exact(exact, expected_exact,
                             "exact prelaunch result tree")
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(True, "prelaunch exact-result recursive schema/types/values")

    result_name = "RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
    check(not (ROOT / result_name).exists() and
          not (ROOT / (result_name + ".STAGING")).exists() and
          not (ROOT / ("." + result_name + ".LOCK")).exists(),
          "packet verification does not launch or lock production")
    print(json.dumps({
        "schema": "GL6FJ_FROZEN_PACKET_VERIFICATION_V003",
        "checks": CHECKS,
        "result": f"PASS__{CHECKS}/{CHECKS}",
        "author_checks": ordinary["checks"],
        "production_run": False,
        "physical_ward_rank_claimed": False,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (Failure, admission.Failure) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
