#!/usr/bin/env python3
"""Sealed-packet replay for the independent GF V004 hostile re-audit."""

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
check(len(manifest) == 7, "audit manifest freezes seven V004 re-audit artifacts")
for expected, rel in manifest:
    path = AUDIT / rel
    check(path.is_file(), f"audit artifact exists: {rel}")
    check(sha(path) == expected, f"audit artifact hash replays: {rel}")

seal = list(rows(AUDIT / "AUDIT_SEAL.sha256"))
check(seal == [(sha(AUDIT / "AUDIT_MANIFEST.sha256"), "AUDIT_MANIFEST.sha256"), (sha(AUDIT / "VERIFICATION.txt"), "VERIFICATION.txt")], "audit seal owns manifest and verification summary")

custody = list(rows(AUDIT / "TARGET_CUSTODY.sha256"))
check(len(custody) == 11, "target custody retains eleven V004 author files")
for expected, rel in custody:
    check(sha(TARGET / rel) == expected, f"target custody replays: {rel}")

prior = list(rows(AUDIT / "PRIOR_REJECTIONS_CUSTODY.sha256"))
check(len(prior) == 9, "prior custody retains all nine V003-rejection artifacts")
for expected, rel in prior:
    check(sha(TARGET / rel) == expected, f"prior rejection custody replays: {rel}")

run = subprocess.run(["python3", str(AUDIT / "independent_reaudit_gf_v004.py")], cwd=ROOT, check=False, capture_output=True, text=True)
check(run.returncode == 0, "independent V004 physics re-audit exits successfully")
check("SUMMARY 208/208" in run.stdout, "independent V004 physics re-audit reproduces 208/208")
check("VERDICT REJECT -- REPAIR_REQUIRED" in run.stdout, "independent replay reproduces the rejection")
check("finite-sigma_max" in run.stdout, "independent replay retains the material upper-bound defect")
check("diag(L,1)" in run.stdout, "independent replay retains the exact divergent-vertex counterexample")

report = (AUDIT / "INDEPENDENT_HOSTILE_REAUDIT_V004.md").read_text()
result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text())
check("Exact counterexample" in report and "false,false,false" in report, "report retains the exact uncovered branch vector")
check("Requested four-case partition -- PASS" in report, "report records the repaired V003 bounded-domain partition")
check("V001--V003" in report, "report preserves the full earlier rejection history")
check(result["verdict"] == "REJECT_REPAIR_REQUIRED", "machine verdict is rejection")
check(result["ge_repin"].startswith("FORBIDDEN"), "machine verdict forbids GE repin")
check("upper-bound divergence" in result["material_defect"], "machine result records the upper-bound defect")

print(f"SUMMARY {checks}/{checks} sealed GF V004 re-audit packet checks passed")
print("VERDICT REJECT -- REPAIR_REQUIRED; GE REPIN FORBIDDEN")
