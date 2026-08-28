#!/usr/bin/env python3
"""Seal and replay verifier for the independent HUST completion audit."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent

AUDIT_PAYLOAD = {
    "CORE_CUSTODY.sha256",
    "LIVE_REQUERY_NORMALIZED.json",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "VERIFICATION_TRANSCRIPT.md",
    "audit_public_completion_search.py",
    "verify_independent_audit.py",
}

checks = 0


def check(statement: bool, label: str) -> None:
    global checks
    if not statement:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entries(path: Path) -> dict[str, str]:
    answer: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line:
            value, relative = line.split("  ", 1)
            answer[relative] = value
    return answer


def main() -> None:
    core = entries(AUDIT / "CORE_CUSTODY.sha256")
    check(len(core) == 9, "nine repaired builder files are pinned")
    for relative, expected in core.items():
        check(digest(LANE / relative) == expected, f"core custody {relative}")

    manifest = entries(AUDIT / "AUDIT_MANIFEST.sha256")
    check(set(manifest) == AUDIT_PAYLOAD, "audit manifest exact membership")
    for relative, expected in manifest.items():
        check(digest(AUDIT / relative) == expected, f"audit payload {relative}")

    seal = entries(AUDIT / "AUDIT_SEAL.sha256")
    check(seal == {
        "AUDIT_MANIFEST.sha256": digest(AUDIT / "AUDIT_MANIFEST.sha256"),
        "VERIFICATION_TRANSCRIPT.md": digest(AUDIT / "VERIFICATION_TRANSCRIPT.md"),
    }, "audit seal binds manifest and transcript")

    for relative in sorted(AUDIT_PAYLOAD | {
            "AUDIT_MANIFEST.sha256", "AUDIT_SEAL.sha256"}):
        raw = (AUDIT / relative).read_bytes()
        check(all(byte in (9, 10) or byte >= 32 for byte in raw)
              and 127 not in raw, f"audit byte hygiene {relative}")
        raw.decode("utf-8")

    for relative in (
            "audit_public_completion_search.py", "verify_independent_audit.py"):
        compile((AUDIT / relative).read_text(encoding="utf-8"),
                relative, "exec")
        check(True, f"Python compile {relative}")

    receipt = json.loads(
        (AUDIT / "LIVE_REQUERY_NORMALIZED.json").read_text(encoding="utf-8"))
    check(receipt["schema"] ==
          "WAC_HUST_PUBLIC_COMPLETION_INDEPENDENT_REQUERY_V001",
          "normalized requery receipt schema")
    check(receipt["nature"]["associated_object_count"] == 7,
          "normalized receipt exact Nature count")
    check(receipt["audit_observation"]["new_qualifying_completion_root_observed"]
          is False, "normalized receipt admits no new root")

    report = (AUDIT / "INDEPENDENT_HOSTILE_AUDIT.md").read_text(
        encoding="utf-8")
    report_flat = " ".join(report.split())
    for phrase in (
        "PASS_AFTER_COMPLETENESS_LEAD_TYPING_AND_CUSTODY_REPAIR",
        "reproducibility/completeness overstatement, repaired",
        "second dissertation lead overtyped, repaired",
        "seven-object custody claim was not pinned, repaired",
        "103/103",
        "48/48",
        "not executable completeness certificates",
        "No numerical value of `G`",
        "one unverified title-only record",
    ):
        check(phrase in report_flat, f"report boundary {phrase}")

    independent = subprocess.run(
        [sys.executable, str(AUDIT / "audit_public_completion_search.py")],
        cwd=LANE.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    check(independent.returncode == 0, "independent replay exits zero")
    check("SUMMARY 103/103 independent checks passed" in independent.stdout,
          "independent replay count")
    check("NATURE_OBJECTS 7" in independent.stdout,
          "independent replay Nature inventory")
    check("DATACITE_RECORDS 3" in independent.stdout,
          "independent replay DataCite count")
    check("ZENODO_RECORDS 2" in independent.stdout,
          "independent replay Zenodo count")
    check("FIGSHARE_MATCHES 0" in independent.stdout,
          "independent replay Figshare result")
    check("PASS_AFTER_COMPLETENESS_LEAD_TYPING_AND_CUSTODY_REPAIR"
          in independent.stdout, "independent replay verdict")

    builder = subprocess.run(
        [sys.executable, str(LANE / "verify_hust_public_completion_search.py")],
        cwd=LANE.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    check(builder.returncode == 0, "repaired builder replay exits zero")
    check("SUMMARY 48/48 checks passed" in builder.stdout,
          "repaired builder replay count")

    print(f"SUMMARY {checks}/{checks} independent-audit seal checks passed")


if __name__ == "__main__":
    main()
