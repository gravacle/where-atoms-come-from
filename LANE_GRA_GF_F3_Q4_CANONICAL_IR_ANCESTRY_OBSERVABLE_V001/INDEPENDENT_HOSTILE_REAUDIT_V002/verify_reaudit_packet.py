#!/usr/bin/env python3
"""Sealed-packet replay for the independent GF V002 hostile re-audit."""

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
check(len(manifest) == 7, "audit manifest freezes seven V002 re-audit artifacts")
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
check(len(custody) == 11, "target custody retains eleven repaired-author files")
for expected, rel in custody:
    check(sha(TARGET / rel) == expected, f"target custody replays: {rel}")

run = subprocess.run(
    ["python3", str(AUDIT / "independent_reaudit_gf_v002.py")],
    cwd=ROOT,
    check=False,
    capture_output=True,
    text=True,
)
check(run.returncode == 0, "independent V002 physics re-audit exits successfully")
check("SUMMARY 153/153" in run.stdout, "independent V002 physics re-audit reproduces 153/153")
check("VERDICT REJECT -- REPAIR_REQUIRED" in run.stdout, "independent replay reproduces the rejection")
check("PASS and FAIL simultaneously" in run.stdout, "independent replay retains the material overlap")

report = (AUDIT / "INDEPENDENT_HOSTILE_REAUDIT_V002.md").read_text()
result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text())
prior = (AUDIT / "PRIOR_V001_REJECTION.md").read_text()
check("Exact counterexample" in report and "Omega_L" in report, "report retains the exact counterexample")
check("outer positive-g2 conjunction still rejects" in report.lower(), "report distinguishes final G2 from amplitude-label ambiguity")
check(result["verdict"] == "REJECT_REPAIR_REQUIRED", "machine verdict is rejection")
check(result["ge_repin"].startswith("FORBIDDEN"), "machine verdict forbids GE repin")
check("V001 hostile audit rejected" in prior, "packet preserves the prior V001 rejection")

print(f"SUMMARY {checks}/{checks} sealed GF V002 re-audit packet checks passed")
print("VERDICT REJECT -- REPAIR_REQUIRED; GE REPIN FORBIDDEN")
