#!/usr/bin/env python3
"""Fail-closed custody gate for the repaired GL6AV author packet."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def require(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(f"FAIL__{label}")
    checks += 1


required = {
    "README.md",
    "RESULT.md",
    "THEOREM.md",
    "DEPENDENCIES.md",
    "SELF_AUDIT.md",
    "VERIFICATION.txt",
    "verify_record_conditioned_metric_bridge.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
}

actual = {path.name for path in HERE.iterdir() if path.is_file()}
require(actual == required, "EXACT_TEN_FILE_INVENTORY")
require(not any(path.is_dir() for path in HERE.iterdir()), "NO_PACKET_SUBDIRECTORIES")


manifest_rows: dict[str, str] = {}
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    require(len(expected) == 64, f"MANIFEST_HASH_LENGTH__{relative}")
    require(all(ch in "0123456789abcdef" for ch in expected),
            f"MANIFEST_HASH_HEX__{relative}")
    require(relative not in manifest_rows, f"MANIFEST_UNIQUE__{relative}")
    target = ROOT / relative
    require(target.is_file(), f"MANIFEST_TARGET_EXISTS__{relative}")
    require(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
            f"MANIFEST_HASH_MATCH__{relative}")
    manifest_rows[relative] = expected

expected_rows = {
    f"{HERE.name}/{name}" for name in required
    if name not in {"MANIFEST.sha256", "SEAL.sha256"}
}
require(set(manifest_rows) == expected_rows, "MANIFEST_EXACT_EIGHT_FILE_COVERAGE")

seal_rows = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
             if line.strip()]
require(len(seal_rows) == 1, "ONE_SEAL_ROW")
seal_hash, seal_target = seal_rows[0].split(maxsplit=1)
require(seal_target == f"{HERE.name}/MANIFEST.sha256", "SEAL_TARGETS_MANIFEST")
require(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
        == seal_hash, "SEAL_HASH_MATCH")

env = os.environ.copy()
env["PYTHONPYCACHEPREFIX"] = "/private/tmp/gl6av_packet_pycache"
exact = HERE / "verify_record_conditioned_metric_bridge.py"
for optimized in (False, True):
    command = [sys.executable]
    if optimized:
        command.append("-O")
    command.append(str(exact))
    completed = subprocess.run(command, cwd=ROOT, env=env,
                               capture_output=True, text=True, check=False)
    require(completed.returncode == 0,
            f"EXACT_REPLAY_EXIT_{'OPT' if optimized else 'NORMAL'}")
    require("PASS__GL6AV_EXACT_REPLAY_COMPLETE__67/67" in completed.stdout,
            f"EXACT_REPLAY_COUNT_{'OPT' if optimized else 'NORMAL'}")

theorem = (HERE / "THEOREM.md").read_text()
result = (HERE / "RESULT.md").read_text()
self_audit = (HERE / "SELF_AUDIT.md").read_text()
for phrase in (
    "authenticated finite binary word",
    "not claimed to be the complete conditioned",
    "not a claim that a global Hamiltonian or global unitary",
    "formal homogeneous coefficient",
    "same invariant state",
    "one nondegenerate quadratic Lorentzian characteristic",
    "single Lorentzian quadratic principal symbol/conformal-",
    "not yet a retained",
    "preparation+phase+single-cone continuum+clock+constitutive+update",
    "gravity; or `G`",
):
    require(phrase in theorem, f"THEOREM_CEILING__{hashlib.sha256(phrase.encode()).hexdigest()[:12]}")
require("formal smooth coupling chart" in result, "RESULT_FORMAL_CHART_TYPING")
require("Is that smooth rank-four chart" in self_audit,
        "SELF_AUDIT_BINARY_LOG_CHART_ATTACK")

print(f"PASS__GL6AV_AUTHOR_PACKET__{checks}/{checks}")
