#!/usr/bin/env python3
"""Verify custody, replay, and seal for the FW independent hostile audit."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent

checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_manifest(path: Path):
    entries = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        expected, relative = line.split("  ", 1)
        entries.append((expected, relative))
    return entries


core = parse_manifest(AUDIT / "CORE_CUSTODY.sha256")
check(len(core) == 17, "core custody lists nine FW files and eight FO/FV dependencies")
check(len({relative for _, relative in core}) == 17,
      "core custody paths are unique")
for expected, relative in core:
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"core custody regular file: {relative}")
    check(digest(path) == expected, f"core custody hash: {relative}")

result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text(encoding="utf-8"))
check(result["disposition"] == "PASS_AFTER_PRE_FREEZE_WITNESS_SCOPE_REPAIR",
      "audit disposition records the repaired scope")
check(result["independent_checks"] == {"passed": 144, "total": 144},
      "audit result records 144/144 independent checks")
check(result["builder_checks"] == {"passed": 114, "total": 114},
      "audit result records 114/114 builder checks")
check(result["remaining_material_defects"] == [],
      "no material defect remains in the frozen witness claim")
check(result["exact_hierarchy"] == {
    "fv_family_offshell_nonidentity_rank": 6,
    "fv_witness_component_rank_mod_identity": 5,
    "source_to_commutator_rank": 3,
    "ground_retarded_spectral_rank": 2,
    "first_nonzero_commutator_moment_rank": 2,
}, "machine-readable hierarchy is exact")

report = (AUDIT / "INDEPENDENT_HOSTILE_AUDIT.md").read_text(encoding="utf-8")
for phrase in (
    "PASS_AFTER_PRE_FREEZE_WITNESS_SCOPE_REPAIR",
    "Q_{\\rm FV-WITNESS}",
    "144/144",
    "generated `Q_diag^(2,4,6)`",
    "not for FV's complete fixed-order source",
    "not a Ward identity",
    "gravity emergence",
):
    check(phrase in report, f"audit report retains scope phrase: {phrase}")

environment = dict(os.environ)
environment["PYTHONDONTWRITEBYTECODE"] = "1"
independent = subprocess.run(
    [sys.executable, str(AUDIT / "audit_fw_witness_response.py")],
    cwd=LANE, env=environment, text=True, capture_output=True, check=False)
check(independent.returncode == 0,
      "independent hostile replay exits successfully")
check("SUMMARY 144/144 independent hostile-audit checks passed" in independent.stdout,
      "independent hostile replay reproduces 144/144")
check("FV_family=6 component_mod_I=5 adH=3 ground_retarded=2 M1=2" in
      independent.stdout, "independent replay reproduces the full hierarchy")
check("PASS for FV-WITNESS only" in independent.stdout,
      "independent replay retains the witness-only ceiling")

builder = subprocess.run(
    [sys.executable, str(LANE / "verify_projected_response.py")],
    cwd=LANE, env=environment, text=True, capture_output=True, check=False)
check(builder.returncode == 0, "frozen FW builder replay exits successfully")
check("SUMMARY 114/114 projected-response checks passed" in builder.stdout,
      "frozen FW builder replay reproduces 114/114")
check("CEILING FV-WITNESS only" in builder.stdout,
      "frozen FW builder replay retains its repaired ceiling")

manifest_entries = parse_manifest(AUDIT / "AUDIT_MANIFEST.sha256")
expected_audit_files = {
    "AUDIT_RESULT.json",
    "CORE_CUSTODY.sha256",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "VERIFICATION.txt",
    "audit_fw_witness_response.py",
    "verify_independent_audit.py",
}
check({relative for _, relative in manifest_entries} == expected_audit_files,
      "audit manifest lists exactly the six frozen audit artifacts")
for expected, relative in manifest_entries:
    path = AUDIT / relative
    check(path.is_file() and not path.is_symlink(),
          f"audit manifest regular file: {relative}")
    check(digest(path) == expected, f"audit manifest hash: {relative}")

seal = parse_manifest(AUDIT / "AUDIT_SEAL.sha256")
check(seal == [(digest(AUDIT / "AUDIT_MANIFEST.sha256"),
                "AUDIT_MANIFEST.sha256")],
      "audit seal binds the final audit manifest")

allowed = {9, 10, 13}
for relative in expected_audit_files:
    payload = (AUDIT / relative).read_bytes()
    check(not any(byte < 32 and byte not in allowed for byte in payload),
          f"audit artifact has no forbidden control byte: {relative}")

print(f"SUMMARY {checks}/{checks} independent-audit custody checks passed")
print("FINAL PASS_AFTER_PRE_FREEZE_WITNESS_SCOPE_REPAIR")
