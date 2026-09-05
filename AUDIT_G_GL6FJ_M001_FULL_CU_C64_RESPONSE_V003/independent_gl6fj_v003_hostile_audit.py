#!/usr/bin/env python3
"""Independent, prelaunch-only hostile audit of frozen GL6FJ V003."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile


sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.resolve()
TARGET = ROOT / "DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
SOURCE = TARGET / "m001_full_cu_c64_projection_v001.cpp"
AUTHOR_VERIFY = TARGET / "verify_gl6fj_v003.py"
CONTROLLER = TARGET / "run_gl6fj_v003_pilot.py"
TARGET_MANIFEST_SHA256 = "685844e9dbee9691fa6c576f2ffb6ec37a99827b29cc1a39e118aa9685b2c0b5"
TARGET_SEAL_FILE_SHA256 = "2204fe01c92bfcec176e183102addf65aef81f721cafe1aa92e9d0c9840b296b"
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def unique(rows):
    value = {}
    for key, item in rows:
        if key in value:
            raise ValueError("duplicate JSON key")
        value[key] = item
    return value


def strict_json_bytes(raw: bytes):
    value = json.loads(raw.decode("utf-8"), object_pairs_hook=unique,
                       parse_constant=lambda token: (_ for _ in ()).throw(
                           ValueError("nonfinite JSON " + token)))

    def visit(node):
        if type(node) is dict:
            for child in node.values():
                visit(child)
        elif type(node) is list:
            for child in node:
                visit(child)
        elif type(node) is int:
            check(-(2 ** 63) <= node <= 2 ** 64 - 1, "JSON integer range")
        elif type(node) is float:
            check(math.isfinite(node), "JSON finite float")
        else:
            check(type(node) in (str, bool, type(None)), "JSON scalar type")

    visit(value)
    return value


def json_file(path: Path):
    return strict_json_bytes(path.read_bytes())


def hash_rows(path: Path) -> dict[str, str]:
    result = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        check(line.count("  ") == 1, "hash separator " + str(number))
        expected, relative = line.split("  ", 1)
        parsed = Path(relative)
        check(len(expected) == 64 and set(expected) <= set("0123456789abcdef") and
              relative and relative not in result and not parsed.is_absolute() and
              ".." not in parsed.parts, "hash row " + str(number))
        result[relative] = expected
    return result


def closed_target() -> None:
    root_mode = os.lstat(TARGET).st_mode
    check(stat.S_ISDIR(root_mode) and not stat.S_ISLNK(root_mode),
          "target is ordinary directory")
    entries = {}
    for node in os.scandir(TARGET):
        mode = node.stat(follow_symlinks=False).st_mode
        check(stat.S_ISREG(mode) and not stat.S_ISLNK(mode),
              "target regular entry " + node.name)
        entries[node.name] = Path(node.path)
    manifest = hash_rows(TARGET / "MANIFEST.sha256")
    payload = set(entries) - {"MANIFEST.sha256", "SEAL.sha256"}
    check(set(manifest) == payload, "target manifest file census")
    for name, expected in manifest.items():
        check(digest(entries[name]) == expected, "target manifest pin " + name)
    seal = hash_rows(TARGET / "SEAL.sha256")
    check(seal == {"MANIFEST.sha256": digest(TARGET / "MANIFEST.sha256")},
          "target seal pin")
    check(digest(TARGET / "MANIFEST.sha256") == TARGET_MANIFEST_SHA256 and
          digest(TARGET / "SEAL.sha256") == TARGET_SEAL_FILE_SHA256,
          "target root fingerprints")


def plan_and_source() -> None:
    plan = json_file(TARGET / "MEASUREMENT_PLAN.json")
    domain = plan["physics_domain"]
    check(domain == {
        "L": 4, "volume": 64, "character": [0, 0, 1],
        "character_label": "m001", "self_conjugate": False,
        "real_source_channels": 24, "symmetric_entries": 300,
        "canonical_scalar_runs": 300, "redundant_scalar_runs": 0,
        "total_scalar_runs": 300}, "exact unaliased m001 domain")
    check(plan["owners"] == [
        "stay_contact", "writer_contact", "stay_spectral",
        "writer_spectral", "stay_writer_interference", "disconnected"],
        "six owner-once plan")
    source = SOURCE.read_text()
    for token in ("RAYS", "m001_self_conjugate", "run_source",
                  "print_scalar_checkpoint", "finite_run_anomalous_zero_imposed"):
        check(token in source, "frozen source surface " + token)
    check("pseudoinverse" not in source and "gravity" not in source,
          "source contains no forbidden response promotion")


def target_replays() -> None:
    for optimized in (False, True):
        command = [sys.executable] + (["-O"] if optimized else []) + [
            "-B", str(AUTHOR_VERIFY)]
        replay = subprocess.run(command, cwd=ROOT, text=True,
                                capture_output=True, timeout=1000, check=False)
        payload = strict_json_bytes(replay.stdout.encode()) if replay.returncode == 0 else {}
        check(replay.returncode == 0 and replay.stderr == "" and
              payload.get("result") == "PASS__172/172" and
              payload.get("production_run") is False and
              payload.get("physical_ward_rank_claimed") is False,
              ("optimized" if optimized else "ordinary") + " author replay")


def sampler_rebuild() -> None:
    histories = set()
    modes = {
        "ordinary": ["-O0"],
        "optimized": ["-O3", "-DNDEBUG"],
        "ubsan": ["-O1", "-fsanitize=undefined", "-fno-sanitize-recover=undefined"],
    }
    with tempfile.TemporaryDirectory(prefix="gl6fj_independent_audit_") as name:
        directory = Path(name)
        for label, flags in modes.items():
            executable = directory / ("m001_" + label)
            build = subprocess.run([
                "clang++", "-std=c++20", *flags, "-Wall", "-Wextra", "-pedantic",
                str(SOURCE), "-o", str(executable)], cwd=ROOT, text=True,
                capture_output=True, timeout=600, check=False)
            check(build.returncode == 0, label + " compilation")
            for selftest in ("catalog-selftest", "character-selftest"):
                run = subprocess.run([str(executable), selftest], cwd=ROOT, text=True,
                                     capture_output=True, timeout=120, check=False)
                payload = strict_json_bytes(run.stdout.encode()) if run.returncode == 0 else {}
                check(run.returncode == 0 and run.stderr == "" and
                      payload.get("result") == "PASS", label + " " + selftest)
            for ray in (0, 23, 24, 299):
                run = subprocess.run([str(executable), "smoke-projection", str(ray)],
                                     cwd=ROOT, text=True, capture_output=True,
                                     timeout=300, check=False)
                row = strict_json_bytes(run.stdout.encode()) if run.returncode == 0 else {}
                check(run.returncode == 0 and row.get("schema") ==
                      "GL6FJ_SCALAR_PROJECTION_RAW_V001" and
                      row.get("character") == [0, 0, 1] and
                      row.get("ray_index") == ray and
                      row.get("production_execution") is False and
                      row.get("finite_run_anomalous_zero_imposed") is False,
                      label + " raw ray " + str(ray))
                for checkpoint in row["checkpoints"]:
                    total = sum(checkpoint[name] for name in (
                        "stay_contact", "writer_contact", "stay_spectral",
                        "writer_spectral", "stay_writer_interference", "disconnected"))
                    check(abs(total - checkpoint["owner_total"]) <=
                          3e-11 * max(1.0, abs(checkpoint["owner_total"])),
                          label + " owner-once balance")
                histories.add((row["burn_digest"], tuple(
                    node["decision_digest"] for node in row["checkpoints"])))
    check(len(histories) == 1, "mode-invariant source-zero history")


def launch_boundary() -> None:
    result = ROOT / "RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
    staging = ROOT / "RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003.STAGING"
    lock = ROOT / ".RESULT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003.LOCK"
    check(not result.exists() and not staging.exists() and not lock.exists(),
          "no preexisting production state")
    invalid = subprocess.run([sys.executable, "-B", str(CONTROLLER), "launch",
                              "--workers", "1", "--authorization", "INVALID"],
                             cwd=ROOT, text=True, capture_output=True,
                             timeout=120, check=False)
    check(invalid.returncode != 0 and "exact production authorization token" in invalid.stderr,
          "invalid authorization rejected")
    check(not result.exists() and not staging.exists() and not lock.exists(),
          "invalid launch leaves no output")


def main() -> int:
    closed_target()
    plan_and_source()
    target_replays()
    sampler_rebuild()
    launch_boundary()
    print(json.dumps({"checks": CHECKS, "result": "PASS",
                      "production_run": False}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, subprocess.SubprocessError) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
