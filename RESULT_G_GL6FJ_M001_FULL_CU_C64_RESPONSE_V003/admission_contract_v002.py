#!/usr/bin/env python3
"""Shared fail-closed admission primitives for GL6FJ V002.

This module changes no physics.  It closes JSON, dependency, and filesystem
admission paths used before or after the unchanged 300-ray flight.
"""

from __future__ import annotations

from hashlib import sha256
import json
import math
import os
from pathlib import Path, PurePosixPath
import stat
from typing import Any, Optional


HEX = frozenset("0123456789abcdef")
U64_MAX = (1 << 64) - 1
I64_MIN = -(1 << 63)
Failure = RuntimeError
DEPENDENCY_LEDGER_SHA256 = (
    "7ad9e4bfac738b124c536e637ac808504ac75479b5246d3b092e929176cb7c29")
DEPENDENCY_ENTRY_COUNT = 70


def fail(message: str) -> None:
    raise Failure(message)


def digest(path: Path) -> str:
    value = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail("duplicate JSON key: " + key)
        result[key] = value
    return result


def _bounded_json_integer(text: str) -> int:
    value = int(text, 10)
    if not I64_MIN <= value <= U64_MAX:
        fail("JSON integer outside admitted 64-bit domain")
    return value


def _finite_json_float(text: str) -> float:
    value = float(text)
    if not math.isfinite(value):
        fail("nonfinite or overflowing JSON float")
    return value


def strict_json_bytes(payload: bytes) -> Any:
    def reject_constant(value: str) -> None:
        fail("nonfinite JSON constant: " + value)

    try:
        value = json.loads(
            payload.decode("utf-8"), object_pairs_hook=_unique_object,
            parse_constant=reject_constant, parse_int=_bounded_json_integer,
            parse_float=_finite_json_float)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise Failure("invalid JSON") from error
    validate_finite_tree(value)
    return value


def strict_json(path: Path) -> Any:
    return strict_json_bytes(Path(path).read_bytes())


def validate_finite_tree(value: Any) -> None:
    stack = [value]
    while stack:
        child = stack.pop()
        if type(child) is dict:
            if not all(type(key) is str for key in child):
                fail("JSON object key is not a string")
            stack.extend(child.values())
        elif type(child) is list:
            stack.extend(child)
        elif type(child) is float:
            if not math.isfinite(child):
                fail("nonfinite JSON number")
        elif type(child) is int:
            if not I64_MIN <= child <= U64_MAX:
                fail("JSON integer outside admitted 64-bit domain")
        elif child is not None and type(child) not in (str, bool):
            fail("non-JSON scalar type")


def deep_exact(actual: Any, expected: Any, label: str = "deep exact value") -> None:
    if type(actual) is not type(expected):
        fail(label + ": type mismatch")
    if type(expected) is dict:
        if set(actual) != set(expected):
            fail(label + ": object-key mismatch")
        for key in expected:
            deep_exact(actual[key], expected[key], label + "/" + key)
    elif type(expected) is list:
        if len(actual) != len(expected):
            fail(label + ": list-length mismatch")
        for index, (left, right) in enumerate(zip(actual, expected)):
            deep_exact(left, right, label + f"/{index}")
    elif type(expected) is float:
        if not math.isfinite(actual) or actual != expected:
            fail(label + ": float mismatch")
    elif type(expected) is int:
        if (not I64_MIN <= actual <= U64_MAX or
                not I64_MIN <= expected <= U64_MAX or actual != expected):
            fail(label + ": bounded integer mismatch")
    elif type(expected) in (str, bool, type(None)):
        if actual != expected:
            fail(label + ": value mismatch")
    else:
        fail(label + ": unsupported non-JSON type")


def exact_object(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict or set(value) != keys:
        fail(label + ": exact object schema")
    return value


def exact_bool(value: Any, label: str) -> bool:
    if type(value) is not bool:
        fail(label + ": exact Boolean")
    return value


def exact_string(value: Any, label: str, *, nonempty: bool = True) -> str:
    if type(value) is not str or (nonempty and value == ""):
        fail(label + ": exact string")
    return value


def exact_int(value: Any, label: str, lower: int = I64_MIN,
              upper: int = U64_MAX) -> int:
    if type(value) is not int or not lower <= value <= upper:
        fail(label + ": exact bounded integer")
    return value


def exact_u64(value: Any, label: str) -> int:
    return exact_int(value, label, 0, U64_MAX)


def exact_number(value: Any, label: str, *, lower: Optional[float] = None,
                 upper: Optional[float] = None) -> float:
    if type(value) not in (int, float):
        fail(label + ": exact finite number")
    answer = float(value)
    if not math.isfinite(answer):
        fail(label + ": finite number")
    if lower is not None and answer < lower:
        fail(label + ": number below domain")
    if upper is not None and answer > upper:
        fail(label + ": number above domain")
    return answer


def safe_relative_name(name: str) -> bool:
    if type(name) is not str or name == "" or "\\" in name:
        return False
    pure = PurePosixPath(name)
    return (not pure.is_absolute() and
            all(part not in {"", ".", ".."} for part in pure.parts))


def safe_dependency_name(name: str) -> bool:
    if type(name) is not str or name == "" or "\\" in name:
        return False
    pure = PurePosixPath(name)
    parts = pure.parts
    return (not pure.is_absolute() and len(parts) >= 2 and parts[0] == ".."
            and all(part not in {"", ".", ".."} for part in parts[1:]))


def ordinary_dependency_path(repository_root: Path, ledger_base: Path,
                             relative: str) -> Path:
    """Resolve one exact ``../...`` dependency without hiding symlink aliases."""
    if not safe_dependency_name(relative):
        fail("unsafe dependency path")
    root = Path(repository_root)
    base = Path(ledger_base)
    root_mode = os.lstat(root).st_mode
    base_mode = os.lstat(base).st_mode
    if (not stat.S_ISDIR(root_mode) or stat.S_ISLNK(root_mode) or
            not stat.S_ISDIR(base_mode) or stat.S_ISLNK(base_mode)):
        fail("ordinary unlinked dependency root/base required")
    if Path(os.path.abspath(base.parent)) != Path(os.path.abspath(root)):
        fail("dependency ledger base must be a direct repository child")
    parts = PurePosixPath(relative).parts[1:]
    cursor = root
    for part in parts[:-1]:
        cursor = cursor / part
        mode = os.lstat(cursor).st_mode
        if not stat.S_ISDIR(mode) or stat.S_ISLNK(mode):
            fail("dependency intermediate component is not an ordinary directory")
    target = cursor / parts[-1]
    mode = os.lstat(target).st_mode
    if not stat.S_ISREG(mode) or stat.S_ISLNK(mode):
        fail("dependency is not an ordinary unlinked file: " + relative)
    try:
        target.resolve().relative_to(root.resolve())
    except ValueError as error:
        raise Failure("dependency escapes repository") from error
    return target


def parse_hashes(path: Path, *, packet_names: bool = True) -> dict[str, str]:
    try:
        text = Path(path).read_bytes().decode("utf-8")
    except UnicodeDecodeError as error:
        raise Failure("hash ledger is not UTF-8") from error
    rows: dict[str, str] = {}
    lines = text.splitlines()
    if not lines:
        fail("empty hash ledger")
    for line in lines:
        if not line.strip():
            fail("blank hash ledger row")
        try:
            value, name = line.split(None, 1)
        except ValueError as error:
            raise Failure("malformed hash row") from error
        name = name.strip()
        if name.startswith("*"):
            fail("star-prefixed hash paths are not admitted")
        if len(value) != 64 or set(value) > HEX:
            fail("lowercase SHA-256 syntax")
        if name in rows:
            fail("duplicate hash path")
        if packet_names and not safe_relative_name(name):
            fail("unsafe relative packet path")
        rows[name] = value
    return rows


def _walk_closed_tree(directory: Path) -> tuple[set[str], set[str]]:
    root_mode = os.lstat(directory).st_mode
    if not stat.S_ISDIR(root_mode) or stat.S_ISLNK(root_mode):
        fail("ordinary unlinked packet directory required")
    files: set[str] = set()
    directories: set[str] = set()

    def visit(current: Path) -> None:
        entries = list(os.scandir(current))
        if current != directory and not entries:
            fail("empty packet directory: " + str(current.relative_to(directory)))
        for entry in entries:
            path = Path(entry.path)
            relative = str(path.relative_to(directory))
            mode = entry.stat(follow_symlinks=False).st_mode
            if stat.S_ISLNK(mode):
                fail("packet symlink: " + relative)
            if stat.S_ISDIR(mode):
                directories.add(relative)
                visit(path)
            elif stat.S_ISREG(mode):
                if (path.name in {"MANIFEST.sha256", "SEAL.sha256"}
                        and path.parent != directory):
                    fail("nested packet control: " + relative)
                files.add(relative)
            else:
                fail("packet special object: " + relative)

    visit(directory)
    return files, directories


def verify_closed_packet(directory: Path) -> tuple[str, str, dict[str, str]]:
    directory = Path(directory)
    files, directories = _walk_closed_tree(directory)
    controls = {"MANIFEST.sha256", "SEAL.sha256"}
    if not controls <= files:
        fail("packet root controls missing")
    manifest_path = directory / "MANIFEST.sha256"
    seal_path = directory / "SEAL.sha256"
    manifest = parse_hashes(manifest_path)
    seal = parse_hashes(seal_path)
    if seal != {"MANIFEST.sha256": digest(manifest_path)}:
        fail("seal does not bind singleton root manifest")
    payload_files = files - controls
    if set(manifest) != payload_files:
        fail("closed packet file census")
    implied_directories: set[str] = set()
    for name in manifest:
        parent = PurePosixPath(name).parent
        while str(parent) != ".":
            implied_directories.add(str(parent))
            parent = parent.parent
    if directories != implied_directories:
        fail("closed packet directory census")
    for name, expected in manifest.items():
        path = directory / name
        if digest(path) != expected:
            fail("packet payload hash: " + name)
    return digest(manifest_path), digest(seal_path), manifest


CONTRACT_KEYS = {
    "schema", "status", "ledger_sha256", "entry_count", "roles",
    "launch_gate",
}
ROLE_KEYS = {
    "runtime_input", "audited_reduction_gate", "pending_acoustic_gate",
    "pending_orbit_gate", "retained_evidence",
}
ROLE_STATES = {
    "runtime_input": "AUTHENTICATED_INPUT",
    "audited_reduction_gate": "AUDITED_LOAD_BEARING",
    "pending_acoustic_gate": "PENDING_NOT_LAUNCH_BEARING",
    "pending_orbit_gate": "PENDING_NOT_LAUNCH_BEARING",
    "retained_evidence": "AUTHENTICATED_INPUT",
}
ROLE_PATHS = {
    "runtime_input": [
        "../DEVELOPMENT_G_GL6ER_COMPILED_FULL_CU_H6_PATH_SCORE_ENGINE_V001/LOCAL_OWNER_TABLES.bin",
        "../DEVELOPMENT_G_GL6ER_COMPILED_FULL_CU_H6_PATH_SCORE_ENGINE_V002/compiled_full_cu_h6_engine_v002.cpp",
        "../DEVELOPMENT_G_GL6FA_DIRECTIONALLY_CONTRACTED_CU_STAY_ENGINE_V002/v001_engine_with_renamed_entrypoint.inc",
        "../DEVELOPMENT_G_GL6FE_M123_M002_FULL_CU_C64_RESPONSE_V001/joint_full_cu_c64_projection_v001.cpp",
    ],
    "audited_reduction_gate": [
        "../DEVELOPMENT_G_GL6FG_SIGNED_CHARACTER_ACOUSTIC_REDUCTION_GATE_V003/EXACT_RESULT.json",
        "../DEVELOPMENT_G_GL6FG_SIGNED_CHARACTER_ACOUSTIC_REDUCTION_GATE_V003/MANIFEST.sha256",
        "../DEVELOPMENT_G_GL6FG_SIGNED_CHARACTER_ACOUSTIC_REDUCTION_GATE_V003/SEAL.sha256",
        "../AUDIT_G_GL6FG_SIGNED_CHARACTER_ACOUSTIC_REDUCTION_GATE_V003/INDEPENDENT_RESULT.json",
        "../AUDIT_G_GL6FG_SIGNED_CHARACTER_ACOUSTIC_REDUCTION_GATE_V003/MANIFEST.sha256",
        "../AUDIT_G_GL6FG_SIGNED_CHARACTER_ACOUSTIC_REDUCTION_GATE_V003/SEAL.sha256",
    ],
    "pending_acoustic_gate": [
        "../DEVELOPMENT_G_GL6FH_PARENT_DERIVED_ACOUSTIC_BUNDLE_GATE_V003/EXACT_RESULT.json",
        "../DEVELOPMENT_G_GL6FH_PARENT_DERIVED_ACOUSTIC_BUNDLE_GATE_V003/THEOREM.md",
        "../DEVELOPMENT_G_GL6FH_PARENT_DERIVED_ACOUSTIC_BUNDLE_GATE_V003/MANIFEST.sha256",
        "../DEVELOPMENT_G_GL6FH_PARENT_DERIVED_ACOUSTIC_BUNDLE_GATE_V003/SEAL.sha256",
    ],
    "pending_orbit_gate": [
        "../DEVELOPMENT_G_GL6FI_MINIMAL_OFF_NYQUIST_WARD_ORBIT_GATE_V001/EXACT_RESULT.json",
        "../DEVELOPMENT_G_GL6FI_MINIMAL_OFF_NYQUIST_WARD_ORBIT_GATE_V001/THEOREM.md",
        "../DEVELOPMENT_G_GL6FI_MINIMAL_OFF_NYQUIST_WARD_ORBIT_GATE_V001/MANIFEST.sha256",
        "../DEVELOPMENT_G_GL6FI_MINIMAL_OFF_NYQUIST_WARD_ORBIT_GATE_V001/SEAL.sha256",
    ],
}

AUDIT_RESULT_KEYS = {
    "schema", "target", "target_modified_by_auditor", "disposition",
    "independent_checks", "target_manifest_sha256",
    "target_seal_file_sha256", "material_defects", "physics_defects",
    "production_run", "promotion", "claim_ceiling",
}
AUDIT_PROMOTION_KEYS = {
    "target_admission_pass", "production_result_promoted",
    "physical_ward_rank_claimed", "gravity_claimed",
}
AUDIT_SCHEMA = "GL6FJ_INDEPENDENT_HOSTILE_AUDIT_RESULT_V002"
AUDIT_TARGET = "DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V002"
AUDIT_DISPOSITION = (
    "PASS__GL6FJ_V002_M001_FULL_RESPONSE_PRELAUNCH_TARGET__"
    "FIVE_ADMISSION_DEFECTS_CLOSED__NO_PRODUCTION"
)
AUDIT_CEILING = (
    "PRELAUNCH_TARGET_AUDIT_ONLY__NO_PRODUCTION_RESPONSE_RANGE_PHYSICAL_"
    "WARD_EINSTEIN_GRAVITY_C_R_OR_G"
)


def validate_audit_result(payload: Any, target_manifest_hash: str,
                          target_seal_hash: str) -> None:
    row = exact_object(payload, AUDIT_RESULT_KEYS, "audit result")
    deep_exact(row["schema"], AUDIT_SCHEMA, "audit schema")
    deep_exact(row["target"], AUDIT_TARGET, "audit target")
    if exact_bool(row["target_modified_by_auditor"],
                  "audit target-modified field"):
        fail("auditor modified target")
    deep_exact(row["disposition"], AUDIT_DISPOSITION, "audit disposition")
    exact_int(row["independent_checks"], "audit independent checks", 1,
              U64_MAX)
    for key, expected in (
            ("target_manifest_sha256", target_manifest_hash),
            ("target_seal_file_sha256", target_seal_hash)):
        value = exact_string(row[key], "audit " + key)
        if len(value) != 64 or set(value) > HEX or value != expected:
            fail("audit target hash binding: " + key)
    deep_exact(row["material_defects"], [], "audit material defects")
    deep_exact(row["physics_defects"], [], "audit physics defects")
    if exact_bool(row["production_run"], "audit production run"):
        fail("prelaunch audit ran production")
    promotion = exact_object(row["promotion"], AUDIT_PROMOTION_KEYS,
                             "audit promotion")
    if not exact_bool(promotion["target_admission_pass"],
                      "audit target admission pass"):
        fail("audit did not pass target admission")
    for key in ("production_result_promoted", "physical_ward_rank_claimed",
                "gravity_claimed"):
        if exact_bool(promotion[key], "audit promotion " + key):
            fail("forbidden audit promotion: " + key)
    deep_exact(row["claim_ceiling"], AUDIT_CEILING, "audit claim ceiling")


def validate_dependency_contract(ledger_path: Path, contract_path: Path,
                                 repository_root: Path,
                                 ledger_base: Path) -> dict[str, Any]:
    contract = exact_object(strict_json(contract_path), CONTRACT_KEYS,
                            "dependency contract")
    deep_exact(contract["schema"], "GL6FJ_TYPED_DEPENDENCY_CONTRACT_V002",
               "dependency contract schema")
    deep_exact(contract["status"],
               "GL6FG_V003_AUDITED__GL6FH_V003_AND_GL6FI_V001_PENDING__LAUNCH_LOCKED",
               "dependency contract status")
    expected_digest = exact_string(contract["ledger_sha256"],
                                   "dependency ledger digest")
    if len(expected_digest) != 64 or set(expected_digest) > HEX:
        fail("dependency ledger digest syntax")
    deep_exact(expected_digest, DEPENDENCY_LEDGER_SHA256,
               "exact dependency ledger digest")
    if digest(ledger_path) != expected_digest:
        fail("exact dependency ledger path/hash census")
    ledger = parse_hashes(ledger_path, packet_names=False)
    exact_int(contract["entry_count"], "dependency entry count", 0, U64_MAX)
    if (contract["entry_count"] != DEPENDENCY_ENTRY_COUNT or
            len(ledger) != DEPENDENCY_ENTRY_COUNT):
        fail("dependency entry census")

    roles = exact_object(contract["roles"], ROLE_KEYS, "dependency roles")
    listed: dict[str, str] = {}
    promotion: dict[str, str] = {}
    role_counts: dict[str, int] = {}
    for role in ROLE_KEYS - {"retained_evidence"}:
        row = exact_object(roles[role], {"promotion_state", "paths"},
                           "dependency role " + role)
        state = exact_string(row["promotion_state"], role + " state")
        deep_exact(state, ROLE_STATES[role], role + " exact promotion state")
        if type(row["paths"]) is not list or not row["paths"]:
            fail(role + " exact nonempty path list")
        deep_exact(row["paths"], ROLE_PATHS[role], role + " exact path census")
        for path in row["paths"]:
            exact_string(path, role + " path")
            if path in listed:
                fail("dependency path assigned twice")
            listed[path] = role
            promotion[path] = state
        role_counts[role] = len(row["paths"])
    retained = exact_object(
        roles["retained_evidence"],
        {"promotion_state", "path_rule", "expected_count"},
        "retained-evidence role")
    deep_exact(retained["path_rule"],
               "EVERY_LEDGER_PATH_NOT_LISTED_IN_ANOTHER_ROLE",
               "retained path rule")
    retained_state = exact_string(retained["promotion_state"],
                                  "retained evidence state")
    deep_exact(retained_state, ROLE_STATES["retained_evidence"],
               "retained evidence exact promotion state")
    exact_int(retained["expected_count"], "retained expected count", 0,
              U64_MAX)
    if retained["expected_count"] != 52:
        fail("exact retained-evidence count")
    if not set(listed) <= set(ledger):
        fail("typed dependency path absent from ledger")
    remainder = set(ledger) - set(listed)
    if len(remainder) != retained["expected_count"]:
        fail("retained-evidence dependency census")
    for path in remainder:
        listed[path] = "retained_evidence"
        promotion[path] = retained_state
    role_counts["retained_evidence"] = len(remainder)
    if set(listed) != set(ledger):
        fail("every dependency has exactly one typed role")

    gate = exact_object(contract["launch_gate"], {
        "audited_roles", "pending_roles",
        "all_pending_must_be_rebound_and_independently_audited",
        "launch_authorized",
    }, "dependency launch gate")
    deep_exact(gate["audited_roles"], ["audited_reduction_gate"],
               "audited dependency roles")
    deep_exact(gate["pending_roles"],
               ["pending_acoustic_gate", "pending_orbit_gate"],
               "pending dependency roles")
    if not exact_bool(gate["all_pending_must_be_rebound_and_independently_audited"],
                      "pending audit requirement"):
        fail("pending dependencies must require rebind and audit")
    if exact_bool(gate["launch_authorized"], "dependency launch"):
        fail("dependency contract must remain launch locked")

    for relative, expected in ledger.items():
        target = ordinary_dependency_path(
            repository_root, ledger_base, relative)
        if digest(target) != expected:
            fail("dependency hash: " + relative)

    return {
        "entries": {
            path: {
                "sha256": ledger[path], "role": listed[path],
                "promotion_state": promotion[path],
            } for path in sorted(ledger)
        },
        "role_counts": role_counts,
        "ledger_sha256": expected_digest,
        "launch_authorized": False,
    }
