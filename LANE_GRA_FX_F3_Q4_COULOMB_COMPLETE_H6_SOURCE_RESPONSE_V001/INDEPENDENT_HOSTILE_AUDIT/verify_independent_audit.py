#!/usr/bin/env python3
"""Verify custody, independent replay, builder replay, and FX audit seal."""

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
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        expected, relative = line.split("  ", 1)
        rows.append((expected, relative))
    return rows


core = parse_manifest(AUDIT / "CORE_CUSTODY.sha256")
check(len(core) == 24 and len({relative for _, relative in core}) == 24,
      "core custody lists nine FX files and fifteen dependencies exactly once")
for expected, relative in core:
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"core custody regular file: {relative}")
    check(digest(path) == expected, f"core custody hash: {relative}")

result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text(encoding="utf-8"))
check(result["disposition"] == "PASS", "machine-readable disposition is PASS")
check(result["independent_checks"] == {"passed": 139, "total": 139},
      "machine-readable independent count is 139/139")
check(result["builder_checks"] == {"passed": 109, "total": 109},
      "machine-readable builder count is 109/109")
check(result["remaining_material_defects"] == [],
      "machine-readable audit has no remaining material defect")
check(result["source_off_coefficients"] ==
      {"a2": "-60", "a4": "-35", "a6": "-893/9"},
      "machine-readable folded coefficients are exact")
check(result["exact_hierarchy"] == {
    "generic_operator_rank_mod_identity": 5,
    "generic_source_to_commutator_rank": 3,
    "generic_ground_retarded_spectral_rank": 2,
    "generic_first_commutator_moment_rank": 2,
    "finite_root_operator_rank_mod_identity": 4,
    "finite_root_source_to_commutator_rank": 2,
    "finite_root_ground_retarded_spectral_rank": 2,
    "finite_root_first_commutator_moment_rank": 2,
}, "machine-readable generic/root response hierarchies are exact")
check(result["finite_polynomial_root"]["physical_threshold"] is False,
      "machine-readable root is not promoted to a threshold")

report = (AUDIT / "INDEPENDENT_HOSTILE_AUDIT.md").read_text(encoding="utf-8")
for phrase in ("**Disposition:** `PASS`", "139/139", "109/109",
               "count-state dynamic program", "formal-series",
               "all eighteen", "rho\\longmapsto\\rho_E",
               "not a physical threshold", "nonzero-momentum", "Ward",
               "thermodynamic", "RGRL-B", "gravity", "Newton"):
    check(phrase in report, f"audit report retains: {phrase}")

environment = dict(os.environ)
environment["PYTHONDONTWRITEBYTECODE"] = "1"
independent = subprocess.run(
    [sys.executable, str(AUDIT / "audit_complete_h6_source.py")],
    cwd=LANE, env=environment, text=True, capture_output=True, check=False)
check(independent.returncode == 0,
      "independent hostile replay exits successfully")
check("SUMMARY 139/139 independent hostile-audit checks passed" in
      independent.stdout, "independent hostile replay reproduces 139/139")
check("EXACT a2=-60 a4=-35 a6=-893/9; all 18 derivative rows replayed" in
      independent.stdout, "independent replay reproduces all folded rows")
check("HIERARCHY generic=5,3,2,2 finite_root=4,2,2,2" in
      independent.stdout, "independent replay reproduces both hierarchies")
check("VERDICT PASS; finite homogeneous selected FO component through H6 only"
      in independent.stdout, "independent replay retains its finite ceiling")

builder = subprocess.run(
    [sys.executable, str(LANE / "derive_complete_h6_source.py")],
    cwd=LANE, env=environment, text=True, capture_output=True, check=False)
check(builder.returncode == 0, "frozen FX builder replay exits successfully")
check("SUMMARY 109/109 complete-H6 source-response checks passed" in
      builder.stdout, "frozen FX builder replay reproduces 109/109")
check("GENERIC_RANKS operator_mod_identity=5 adH=3 retarded=2 M1=2" in
      builder.stdout, "frozen FX builder reproduces the generic hierarchy")
check("TRUNCATED_ROOT x=0.5398271903 ranks=4,2,2,2; NOT_A_THRESHOLD" in
      builder.stdout, "frozen FX builder retains the root ceiling")

manifest = parse_manifest(AUDIT / "AUDIT_MANIFEST.sha256")
expected_files = {
    "AUDIT_RESULT.json", "CORE_CUSTODY.sha256",
    "INDEPENDENT_HOSTILE_AUDIT.md", "VERIFICATION.txt",
    "audit_complete_h6_source.py", "verify_independent_audit.py",
}
check({relative for _, relative in manifest} == expected_files,
      "audit manifest lists exactly six frozen audit payloads")
for expected, relative in manifest:
    path = AUDIT / relative
    check(path.is_file() and not path.is_symlink(),
          f"audit manifest regular file: {relative}")
    check(digest(path) == expected, f"audit manifest hash: {relative}")

seal = parse_manifest(AUDIT / "AUDIT_SEAL.sha256")
check(seal == [(digest(AUDIT / "AUDIT_MANIFEST.sha256"),
                "AUDIT_MANIFEST.sha256")],
      "audit seal binds exactly the audit manifest")

allowed = {9, 10, 13}
for relative in expected_files:
    payload = (AUDIT / relative).read_bytes()
    check(not any(byte < 32 and byte not in allowed for byte in payload),
          f"audit payload has no forbidden control byte: {relative}")

print(f"SUMMARY {checks}/{checks} independent-audit custody checks passed")
print("FINAL PASS")
