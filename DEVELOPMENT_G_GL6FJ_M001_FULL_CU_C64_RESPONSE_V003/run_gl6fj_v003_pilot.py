#!/usr/bin/env python3
"""Strict write-once controller for the GL6FJ V003 m001 pilot.

Exactly 300 canonical projections reconstruct Sym(24) at L4 m001. Every
completed raw payload is streamed to an exclusively created staging packet.
No result is promoted unless the frozen target, every dependency, and a
distinct hostile audit authenticate before and after flight.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import ctypes
import errno
from hashlib import sha256
import importlib.util
import json
import math
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterator

import admission_contract_v003 as admission

sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = HERE / "m001_full_cu_c64_projection_v001.cpp"
RECONSTRUCTOR = HERE / "reconstruct_m001_full_response_v002.py"
EXPECTED_HISTORY = HERE / "EXPECTED_HISTORY.json"
RESULT_VERIFIER = HERE / "result_packet_verifier_v003.py"
CONTROLLER = Path(__file__).resolve()
PLAN = HERE / "MEASUREMENT_PLAN.json"
DEPENDENCIES = HERE / "DEPENDENCIES.sha256"
DEPENDENCY_CONTRACT = HERE / "DEPENDENCY_CONTRACT.json"
ADMISSION_CONTRACT = HERE / "admission_contract_v003.py"
RECONSTRUCTOR_SUPPORT = HERE / "admission_contract_v002.py"
AUDIT = ROOT / "AUDIT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
AUDIT_RESULT = AUDIT / "INDEPENDENT_RESULT.json"
AUDIT_SCHEMA = admission.AUDIT_SCHEMA
AUDIT_DISPOSITION = admission.AUDIT_DISPOSITION
RESULT_NAME = "RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
RESULT = ROOT / RESULT_NAME
STAGING = ROOT / (RESULT_NAME + ".STAGING")
OUTPUT_LOCK = ROOT / ("." + RESULT_NAME + ".LOCK")
TOKEN = "AUTHORIZE_GL6FJ_V003_SINGLE_300_RUN_C64_M001_PILOT"
STAGING_VERIFICATION_TOKEN = "VERIFY_GL6FJ_V003_STAGING_BEFORE_ATOMIC_PROMOTION"
RAW_TOKEN = "GL6FJ_V001_SUBORDINATE_300_RAY_RAW_NOT_PRODUCTION"
UPSTREAM_DEPENDENCY_CHAIN_CURRENT = True
HEX = frozenset("0123456789abcdef")


class Failure(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise Failure(message)


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


def strict_json_file(path: Path) -> Any:
    return strict_json_bytes(path.read_bytes())


def _safe_packet_name(name: str) -> bool:
    pure = PurePosixPath(name)
    return (name != "" and not pure.is_absolute() and
            all(part not in {"", ".", ".."} for part in pure.parts))


def parse_hashes(path: Path, *, packet_names: bool = True) -> dict[str, str]:
    try:
        return admission.parse_hashes(path, packet_names=packet_names)
    except admission.Failure as error:
        raise Failure(str(error)) from error


def verify_closed_packet(directory: Path) -> tuple[str, str]:
    try:
        manifest_hash, seal_hash, _manifest = \
            admission.verify_closed_packet(directory)
        return manifest_hash, seal_hash
    except admission.Failure as error:
        raise Failure(str(error)) from error


def verify_external_dependencies() -> str:
    try:
        typed = admission.validate_dependency_contract(
            DEPENDENCIES, DEPENDENCY_CONTRACT, ROOT, HERE)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(typed["upstream_chain_current"] is True and
          typed["launch_authorized"] is False,
          "upstream chain current but dependency contract alone cannot launch")
    return digest(DEPENDENCIES)


def code_pins() -> dict[str, str]:
    return {
        "source_sha256": digest(SOURCE),
        "reconstructor_sha256": digest(RECONSTRUCTOR),
        "expected_history_sha256": digest(EXPECTED_HISTORY),
        "result_verifier_sha256": digest(RESULT_VERIFIER),
        "controller_sha256": digest(CONTROLLER),
        "measurement_plan_sha256": digest(PLAN),
        "dependencies_sha256": digest(DEPENDENCIES),
        "dependency_contract_sha256": digest(DEPENDENCY_CONTRACT),
        "admission_contract_sha256": digest(ADMISSION_CONTRACT),
        "reconstructor_support_sha256": digest(RECONSTRUCTOR_SUPPORT),
    }


def verify_authorization() -> dict[str, str]:
    check(UPSTREAM_DEPENDENCY_CHAIN_CURRENT,
          "GL6FJ upstream dependency chain must be current")
    target_manifest_hash, target_seal_hash = verify_closed_packet(HERE)
    dependency_hash = verify_external_dependencies()
    audit_manifest_hash, audit_seal_hash = verify_closed_packet(AUDIT)
    audit = strict_json_file(AUDIT_RESULT)
    try:
        admission.validate_audit_result(
            audit, target_manifest_hash, target_seal_hash)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check(parse_hashes(AUDIT / "MANIFEST.sha256").get(
        "INDEPENDENT_RESULT.json") == digest(AUDIT_RESULT),
        "audit manifest binds authorization result")
    pins = code_pins()
    check(pins["dependencies_sha256"] == dependency_hash,
          "dependency ledger pin is stable")
    return {
        "target_manifest_sha256": target_manifest_hash,
        "target_seal_file_sha256": target_seal_hash,
        "audit_manifest_sha256": audit_manifest_hash,
        "audit_seal_file_sha256": audit_seal_hash,
        "audit_result_sha256": digest(AUDIT_RESULT),
        **pins,
    }


def load_reconstructor():
    specification = importlib.util.spec_from_file_location(
        "gl6fj_m001_reconstruct_controller", RECONSTRUCTOR)
    check(specification is not None and specification.loader is not None,
          "reconstructor import specification")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def exclusive_write(path: Path, payload: bytes, mode: int = 0o444) -> None:
    try:
        descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, mode)
    except FileExistsError as error:
        raise Failure("exclusive output already exists: " + str(path)) from error
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
    except Exception:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise


def canonical_json(payload: Any) -> bytes:
    return (json.dumps(payload, sort_keys=True, indent=2,
                       allow_nan=False) + "\n").encode()


def run_job(executable: Path, ray: int) -> tuple[str, dict[str, Any], bytes]:
    check(type(ray) is int and 0 <= ray < 300, "canonical job domain")
    name = f"RAW/m001_ray_{ray:03d}.json"
    command = [str(executable), "pilot-projection-raw", str(ray), RAW_TOKEN]
    try:
        process = subprocess.run(command, cwd=ROOT, text=False,
                                 capture_output=True, timeout=1200)
    except subprocess.TimeoutExpired as error:
        raise Failure("subordinate timeout: " + name) from error
    check(process.returncode == 0,
          f"subordinate {name} failed: " +
          process.stderr[-1000:].decode("utf-8", errors="replace"))
    payload = strict_json_bytes(process.stdout)
    check(type(payload) is dict, "subordinate JSON object: " + name)
    encoded = canonical_json(payload)
    return name, payload, encoded


def stream_jobs(executable: Path,
                jobs: list[int],
                workers: int) -> Iterator[tuple[str, dict[str, Any], bytes]]:
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(run_job, executable, job): job
                   for job in jobs}
        try:
            for future in as_completed(futures):
                yield future.result()
        except Exception:
            for pending in futures:
                pending.cancel()
            raise


def write_hash_ledger(path: Path, rows: dict[str, bytes]) -> None:
    payload = "".join(
        f"{sha256(rows[name]).hexdigest()}  {name}\n" for name in sorted(rows))
    exclusive_write(path, payload.encode())


def compiler_identity() -> tuple[Path, str, str]:
    requested = os.environ.get("CXX", "clang++")
    check(requested != "" and not any(character.isspace()
                                       for character in requested),
          "CXX must name one executable without arguments")
    found = shutil.which(requested)
    check(found is not None, "C++ compiler exists")
    compiler = Path(found).resolve()
    check(compiler.is_file(), "ordinary C++ compiler file")
    version = subprocess.run([str(compiler), "--version"], cwd=ROOT,
                             text=True, capture_output=True, timeout=30)
    check(version.returncode == 0, "compiler version query")
    version_text = version.stdout.strip()
    check(version_text != "", "nonempty compiler version")
    return compiler, digest(compiler), version_text


def acquire_output_lock(bindings: dict[str, str]) -> None:
    check(OUTPUT_LOCK.parent.resolve() == ROOT.resolve(),
          "output lock has exact repository parent")
    payload = {
        "schema": "GL6FJ_EXCLUSIVE_OUTPUT_LOCK_V003",
        "pid": os.getpid(),
        "target_manifest_sha256": bindings["target_manifest_sha256"],
        "audit_result_sha256": bindings["audit_result_sha256"],
        "result_packet": RESULT_NAME,
    }
    exclusive_write(OUTPUT_LOCK, canonical_json(payload), 0o444)


def rename_no_replace(source: Path, destination: Path) -> None:
    """Atomically rename a directory without replacing an existing target."""
    check(sys.platform == "darwin",
          "atomic no-replace promotion requires Darwin renamex_np")
    libc = ctypes.CDLL(None, use_errno=True)
    function = libc.renamex_np
    function.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
    function.restype = ctypes.c_int
    rename_excl = 0x00000004
    result = function(os.fsencode(source), os.fsencode(destination), rename_excl)
    if result != 0:
        number = ctypes.get_errno()
        if number == errno.EEXIST:
            raise Failure("final result appeared before no-replace promotion")
        raise Failure("renamex_np failed: " + os.strerror(number))


def result_manifest_bytes(directory: Path) -> bytes:
    try:
        files, _directories = admission._walk_closed_tree(directory)
    except admission.Failure as error:
        raise Failure(str(error)) from error
    check("MANIFEST.sha256" not in files and "SEAL.sha256" not in files,
          "staging root controls absent before sealing")
    payload_paths = sorted(directory / name for name in files)
    return "".join(
        f"{digest(path)}  {path.relative_to(directory)}\n"
        for path in payload_paths).encode()


def run_production(args: argparse.Namespace,
                   bindings_before: dict[str, str]) -> str:
    check(RESULT.parent.resolve() == ROOT.resolve() and
          STAGING.parent.resolve() == ROOT.resolve(), "exact result parents")
    check(RESULT.name == RESULT_NAME and
          STAGING.name == RESULT_NAME + ".STAGING", "exact result names")
    check(not RESULT.exists() and not STAGING.exists(),
          "write-once result and staging paths absent under lock")
    STAGING.mkdir(mode=0o755)
    (STAGING / "RAW").mkdir(mode=0o755)
    module = load_reconstructor()
    raw_bytes: dict[str, bytes] = {}
    semantic_rows: dict[str, dict[str, Any]] = {}
    raw_rows: dict[int, dict[str, Any]] = {}
    compiler, compiler_hash, compiler_version = compiler_identity()
    build_command = [str(compiler), "-std=c++20", "-O3", "-DNDEBUG",
                     "-Wall", "-Wextra", "-pedantic", str(SOURCE)]
    started_ns = time.time_ns()

    with tempfile.TemporaryDirectory(prefix="gl6fj_m001_production_") as temporary:
        executable = Path(temporary) / "gl6fj_m001_runner"
        complete_build_command = build_command + ["-o", str(executable)]
        build = subprocess.run(complete_build_command, cwd=ROOT, text=True,
                               capture_output=True, timeout=600)
        check(build.returncode == 0,
              "production build: " + build.stderr[-1000:])
        executable_hash = digest(executable)
        executable_bytes = executable.read_bytes()

        canonical_jobs = list(range(300))
        for completed, (name, payload, encoded) in enumerate(
                stream_jobs(executable, canonical_jobs, args.workers), start=1):
            exclusive_write(STAGING / name, encoded)
            raw_bytes[name] = encoded
            semantic_rows[name] = payload
            ray = payload.get("ray_index")
            check(type(ray) is int and ray not in raw_rows,
                  "unique typed canonical result key")
            raw_rows[ray] = payload
            if completed % 20 == 0 or completed == 300:
                print(f"GL6FJ m001 canonical {completed}/300", file=sys.stderr,
                      flush=True)

        check(len(raw_bytes) == 300 and set(raw_rows) == set(range(300)),
              "closed in-memory 300-raw census")
        module.validate_common_raw(raw_rows)
        bindings_after = verify_authorization()
        try:
            admission.deep_exact(bindings_after, bindings_before,
                                 "target/audit/dependencies unchanged during production")
        except admission.Failure as error:
            raise Failure(str(error)) from error
        reconstruction = module.reconstruct(raw_rows)
        exclusive_write(STAGING / "FULL_RESPONSE_RECONSTRUCTION.json",
                        canonical_json(reconstruction))
        write_hash_ledger(STAGING / "RAW_MANIFEST.sha256", raw_bytes)
        exclusive_write(
            STAGING / "RAW_SEAL.sha256",
            (f"{digest(STAGING / 'RAW_MANIFEST.sha256')}  "
             "RAW_MANIFEST.sha256\n").encode())
        semantic_digest = sha256(json.dumps(
            semantic_rows, sort_keys=True, separators=(",", ":"),
            allow_nan=False).encode()).hexdigest()
        exclusive_write(STAGING / "reconstruct_m001_full_response_v002.py",
                        RECONSTRUCTOR.read_bytes())
        exclusive_write(STAGING / "EXPECTED_HISTORY.json",
                        EXPECTED_HISTORY.read_bytes())
        exclusive_write(STAGING / "FROZEN_DEPENDENCIES.sha256",
                        DEPENDENCIES.read_bytes())
        exclusive_write(STAGING / "FROZEN_DEPENDENCY_CONTRACT.json",
                        DEPENDENCY_CONTRACT.read_bytes())
        exclusive_write(STAGING / "admission_contract_v003.py",
                        ADMISSION_CONTRACT.read_bytes())
        exclusive_write(STAGING / "admission_contract_v002.py",
                        RECONSTRUCTOR_SUPPORT.read_bytes())
        exclusive_write(STAGING / "gl6fj_m001_runner_frozen.bin",
                        executable_bytes)

    completed_ns = time.time_ns()
    custody = {
        "schema": "GL6FJ_PRODUCTION_LAUNCH_CUSTODY_V003",
        **bindings_before,
        "compiler_path": str(compiler),
        "compiler_sha256": compiler_hash,
        "compiler_version": compiler_version,
        "build_command_without_temporary_output": build_command,
        "executable_sha256": executable_hash,
        "workers": args.workers,
        "started_unix_ns": started_ns,
        "completed_unix_ns": completed_ns,
        "elapsed_seconds": (completed_ns - started_ns) / 1_000_000_000.0,
        "character": [0, 0, 1],
        "raw_files": 300,
        "canonical_files": 300,
        "redundant_probe_files": 0,
        "raw_semantic_sha256": semantic_digest,
        "raw_manifest_sha256": digest(STAGING / "RAW_MANIFEST.sha256"),
        "raw_seal_file_sha256": digest(STAGING / "RAW_SEAL.sha256"),
        "postrun_authorization_reauthenticated": True,
        "output_lock_retained_through_promotion": True,
        "promotion_primitive": "darwin_renamex_np_RENAME_EXCL",
        "production_run": True,
        "population_or_seed_ladder_run": False,
        "translation_forced_anomalous_expectation_zero": True,
        "finite_run_anomalous_zero_imposed": False,
        "physical_ward_rank_claimed": False,
    }
    check(math.isfinite(custody["elapsed_seconds"]) and
          custody["elapsed_seconds"] > 0, "finite positive elapsed time")
    exclusive_write(STAGING / "LAUNCH_CUSTODY.json", canonical_json(custody))
    exclusive_write(STAGING / "verify_result_packet.py",
                    RESULT_VERIFIER.read_bytes())
    exclusive_write(STAGING / "MANIFEST.sha256",
                    result_manifest_bytes(STAGING))
    exclusive_write(
        STAGING / "SEAL.sha256",
        (f"{digest(STAGING / 'MANIFEST.sha256')}  MANIFEST.sha256\n").encode())

    verification = subprocess.run(
        [sys.executable, "-B", str(STAGING / "verify_result_packet.py"),
         STAGING_VERIFICATION_TOKEN],
        cwd=ROOT, text=True, capture_output=True, timeout=1800)
    check(verification.returncode == 0,
          "staging result verification: " + verification.stderr[-2000:])
    check(OUTPUT_LOCK.is_file(), "exclusive output lock retained")
    check(not RESULT.exists(), "final result remains absent before promotion")
    try:
        admission.deep_exact(
            verify_authorization(), bindings_before,
            "final target/audit/dependency reauthentication")
    except admission.Failure as error:
        raise Failure(str(error)) from error
    rename_no_replace(STAGING, RESULT)
    final_verification = subprocess.run(
        [sys.executable, "-B", str(RESULT / "verify_result_packet.py")],
        cwd=ROOT, text=True, capture_output=True, timeout=1800)
    check(final_verification.returncode == 0,
          "final result verification: " + final_verification.stderr[-2000:])
    return final_verification.stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("preflight", "launch"))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--authorization", default="")
    args = parser.parse_args()
    check(type(args.workers) is int and 1 <= args.workers <= 4,
          "workers in fixed nonphysics range 1..4")
    if args.mode == "preflight":
        try:
            bindings = verify_authorization()
            available = (not RESULT.exists() and not STAGING.exists() and
                         not OUTPUT_LOCK.exists())
        except (Failure, FileNotFoundError, OSError) as error:
            print(json.dumps({
                "schema": "GL6FJ_PRODUCTION_PREFLIGHT_V003",
                "authorized": False,
                "reason": str(error),
                "production_run": False,
            }, sort_keys=True))
            return 0
        print(json.dumps({
            "schema": "GL6FJ_PRODUCTION_PREFLIGHT_V003",
            "authorized": available,
            "output_available": available,
            "bindings": bindings,
            "production_run": False,
        }, sort_keys=True))
        return 0

    check(args.authorization == TOKEN, "exact production authorization token")
    bindings = verify_authorization()
    acquired = False
    try:
        acquire_output_lock(bindings)
        acquired = True
        output = run_production(args, bindings)
        print(output)
        return 0
    finally:
        if acquired and OUTPUT_LOCK.exists():
            OUTPUT_LOCK.unlink()


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (Failure, admission.Failure, OSError,
            subprocess.SubprocessError) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
