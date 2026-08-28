#!/usr/bin/env python3
"""Sealed-packet replay for the independent GF V005 design-contract audit."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
TARGET = AUDIT.parent
ROOT = TARGET.parent
checks = 0


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path):
    for raw in path.read_text().splitlines():
        if raw.strip():
            digest, rel = raw.split("  ", 1)
            yield digest, rel


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


manifest = list(rows(AUDIT / "AUDIT_MANIFEST.sha256"))
check(len(manifest) == 7, "audit manifest freezes seven V005 re-audit artifacts")
for expected, rel in manifest:
    path = AUDIT / rel
    check(path.is_file(), f"audit artifact exists: {rel}")
    check(sha(path) == expected, f"audit artifact hash replays: {rel}")

seal = list(rows(AUDIT / "AUDIT_SEAL.sha256"))
check(seal == [(sha(AUDIT / "AUDIT_MANIFEST.sha256"), "AUDIT_MANIFEST.sha256"), (sha(AUDIT / "VERIFICATION.txt"), "VERIFICATION.txt")], "audit seal owns manifest and verification summary")

custody = list(rows(AUDIT / "TARGET_CUSTODY.sha256"))
check(len(custody) == 11, "target custody retains eleven V005 author files")
for expected, rel in custody:
    check(sha(TARGET / rel) == expected, f"target custody replays: {rel}")

prior = list(rows(AUDIT / "PRIOR_REJECTIONS_CUSTODY.sha256"))
check(len(prior) == 9, "prior custody retains all nine V004-rejection artifacts")
for expected, rel in prior:
    check(sha(TARGET / rel) == expected, f"prior rejection custody replays: {rel}")

run = subprocess.run(["python3", str(AUDIT / "independent_reaudit_gf_v005.py")], cwd=ROOT, check=False, capture_output=True, text=True)
check(run.returncode == 0, "independent V005 re-audit exits successfully")
check("SUMMARY 270/270" in run.stdout, "independent V005 re-audit reproduces 270/270")
check("PASS -- NARROW_DESIGN_CONTRACT_SEAL" in run.stdout, "independent replay reproduces narrow PASS")
check("disjoint and total" in run.stdout, "independent replay retains complete classifier result")
check("diag(L,1) upper divergence FAIL" in run.stdout, "independent replay retains V004 regression")
check("design contract only" in run.stdout, "independent replay retains strict ceiling")

report = (AUDIT / "INDEPENDENT_HOSTILE_REAUDIT_V005.md").read_text()
result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text())
check("Exhaustive amplitude decision" in report and "V004 upper-bound regression -- PASS" in report, "report retains exhaustive classifier and upper-bound regression")
check("V001--V004" in report, "report preserves full prior rejection history")
check("no positive G2" in report and "GE was neither audited nor repinned" in report, "report retains design-only and no-GE ceilings")
check(result["verdict"] == "PASS_NARROW_DESIGN_CONTRACT_SEAL", "machine verdict is narrow design PASS")
check(result["ge_status"] == "NOT_AUDITED_NOT_REPINNED", "machine verdict does not promote GE")
check(result["ceiling"].startswith("design contract only"), "machine result retains strict physics ceiling")

print(f"SUMMARY {checks}/{checks} sealed GF V005 re-audit packet checks passed")
print("VERDICT PASS -- NARROW DESIGN-CONTRACT SEAL; PHYSICS UNEXECUTED")
