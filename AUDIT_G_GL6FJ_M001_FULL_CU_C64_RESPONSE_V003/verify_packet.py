#!/usr/bin/env python3
"""Verify the closed GL6FJ V003 independent prelaunch-audit packet."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "DEVELOPMENT_G_GL6FJ_M001_FULL_CU_C64_RESPONSE_V003"
PAYLOAD = {
    "AUDIT_REPORT.md", "INDEPENDENT_RESULT.json", "README.md", "TARGET.sha256",
    "VERIFICATION.txt", "independent_gl6fj_v003_hostile_audit.py", "verify_packet.py",
}
CHECKS = 0


def check(value: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(label)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        check(len(expected) == 64 and relative not in out, "well-formed hash row")
        out[relative] = expected
    return out


def closed() -> None:
    seen = set()
    for entry in os.scandir(HERE):
        mode = entry.stat(follow_symlinks=False).st_mode
        check(stat.S_ISREG(mode) and not stat.S_ISLNK(mode), "ordinary audit entry")
        seen.add(entry.name)
    check(seen == PAYLOAD | {"MANIFEST.sha256", "SEAL.sha256"},
          "exact audit file census")
    manifest = rows(HERE / "MANIFEST.sha256")
    check(set(manifest) == PAYLOAD, "exact audit payload census")
    for name, expected in manifest.items():
        check(digest(HERE / name) == expected, "audit manifest hash " + name)
    check(rows(HERE / "SEAL.sha256") ==
          {"MANIFEST.sha256": digest(HERE / "MANIFEST.sha256")},
          "audit seal")


def admission() -> None:
    specification = importlib.util.spec_from_file_location(
        "gl6fj_admission", TARGET / "admission_contract_v003.py")
    check(specification is not None and specification.loader is not None,
          "admission import")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    target_manifest, target_seal, _ = module.verify_closed_packet(TARGET)
    result = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text())
    module.validate_audit_result(result, target_manifest, target_seal)
    check(result["independent_checks"] == 1627, "frozen independent check count")


def replay() -> None:
    run = subprocess.run([sys.executable, "-B", str(
        HERE / "independent_gl6fj_v003_hostile_audit.py")], cwd=ROOT,
        text=True, capture_output=True, timeout=2400, check=False)
    payload = json.loads(run.stdout) if run.returncode == 0 else {}
    check(run.returncode == 0 and run.stderr == "" and payload == {
        "checks": 1627, "production_run": False, "result": "PASS"},
        "independent replay")


if __name__ == "__main__":
    try:
        closed()
        admission()
        replay()
        print("GL6FJ_V003_AUDIT_PACKET_PASS " + str(CHECKS) + "/" + str(CHECKS))
    except (AssertionError, OSError, ValueError, subprocess.SubprocessError) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
