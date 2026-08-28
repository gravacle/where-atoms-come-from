#!/usr/bin/env python3
"""Sealed-packet replay for the independent GF V003 hostile re-audit."""

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
check(len(manifest) == 7, "audit manifest freezes seven V003 re-audit artifacts")
for expected, rel in manifest:
    path = AUDIT / rel
    check(path.is_file(), f"audit artifact exists: {rel}")
    check(sha(path) == expected, f"audit artifact hash replays: {rel}")

seal = list(rows(AUDIT / "AUDIT_SEAL.sha256"))
check(
    seal
    == [
        (sha(AUDIT / "AUDIT_MANIFEST.sha256"), "AUDIT_MANIFEST.sha256"),
        (sha(AUDIT / "VERIFICATION.txt"), "VERIFICATION.txt"),
    ],
    "audit seal owns the manifest and frozen verification summary",
)

custody = list(rows(AUDIT / "TARGET_CUSTODY.sha256"))
check(len(custody) == 11, "target custody retains eleven V003 author files")
for expected, rel in custody:
    check(sha(TARGET / rel) == expected, f"target custody replays: {rel}")

prior = list(rows(AUDIT / "PRIOR_REJECTIONS_CUSTODY.sha256"))
check(len(prior) == 9, "prior custody retains all nine V002-rejection artifacts")
for expected, rel in prior:
    check(sha(TARGET / rel) == expected, f"prior rejection custody replays: {rel}")

run = subprocess.run(
    ["python3", str(AUDIT / "independent_reaudit_gf_v003.py")],
    cwd=ROOT,
    check=False,
    capture_output=True,
    text=True,
)
check(run.returncode == 0, "independent V003 physics re-audit exits successfully")
check("SUMMARY 190/190" in run.stdout, "independent V003 physics re-audit reproduces 190/190")
check("VERDICT REJECT -- REPAIR_REQUIRED" in run.stdout, "independent replay reproduces the rejection")
check("disjoint but not total" in run.stdout, "independent replay retains the material totality defect")
check("N_raw,L=I_2" in run.stdout, "independent replay retains the exact no-D counterexample")

report = (AUDIT / "INDEPENDENT_HOSTILE_REAUDIT_V003.md").read_text()
result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text())
check("Exact counterexample" in report and "false,false,false" in report, "report retains the exact uncovered branch vector")
check("V002 overlap separation -- PASS" in report, "report records the repaired V002 regression")
check("V001 and V002" in report, "report preserves both earlier rejection histories")
check(result["verdict"] == "REJECT_REPAIR_REQUIRED", "machine verdict is rejection")
check(result["ge_repin"].startswith("FORBIDDEN"), "machine verdict forbids GE repin")
check("not total" in result["material_defect"], "machine result records the totality defect")

print(f"SUMMARY {checks}/{checks} sealed GF V003 re-audit packet checks passed")
print("VERDICT REJECT -- REPAIR_REQUIRED; GE REPIN FORBIDDEN")
