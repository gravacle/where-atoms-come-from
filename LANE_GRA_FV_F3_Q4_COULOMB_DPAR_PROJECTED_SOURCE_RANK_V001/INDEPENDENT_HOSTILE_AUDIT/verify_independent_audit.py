#!/usr/bin/env python3
"""Seal and replay verifier for the independent FV hostile audit."""

from hashlib import sha256
from pathlib import Path
import subprocess
import sys


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent

AUDIT_PAYLOAD = {
    "CORE_CUSTODY.sha256",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "VERIFICATION_TRANSCRIPT.md",
    "audit_projected_source_rank.py",
    "verify_independent_audit.py",
}

checks = 0


def check(statement, label):
    global checks
    if not statement:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def entries(path):
    answer = {}
    for line in path.read_text().splitlines():
        if line:
            value, relative = line.split("  ", 1)
            answer[relative] = value
    return answer


def main():
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

    for relative in sorted(AUDIT_PAYLOAD | {"AUDIT_MANIFEST.sha256", "AUDIT_SEAL.sha256"}):
        raw = (AUDIT / relative).read_bytes()
        check(all(byte in (9, 10) or byte >= 32 for byte in raw) and 127 not in raw,
              f"audit byte hygiene {relative}")
        raw.decode("utf-8")

    for relative in ("audit_projected_source_rank.py", "verify_independent_audit.py"):
        compile((AUDIT / relative).read_text(), relative, "exec")
        check(True, f"Python compile {relative}")

    report = (AUDIT / "INDEPENDENT_HOSTILE_AUDIT.md").read_text()
    report_flat = " ".join(report.split())
    for phrase in (
        "PASS_AFTER_FV_PURE_PREMISE_AND_BYTE_HYGIENE_REPAIR",
        "material premise-inheritance defect, repaired",
        "literal carriage-return byte",
        "83/83",
        "-\\frac{4678629417}{256}",
        "operator evaluations, not a rank of six coefficient vectors",
        "not a retarded/CTP response theorem",
        "Failure of any antecedent requires a fresh complete-source calculation",
    ):
        check(phrase in report_flat, f"report boundary {phrase}")

    independent = subprocess.run(
        [sys.executable, str(AUDIT / "audit_projected_source_rank.py")],
        cwd=LANE.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    check(independent.returncode == 0, "independent scientific replay exits zero")
    check("SUMMARY 83/83" in independent.stdout, "independent replay count")
    check("WITNESS_DETERMINANT -4678629417/256" in independent.stdout,
          "independent replay determinant")
    check("no CTP/Ward/pole/gravity/G" in independent.stdout,
          "independent replay retains ceiling")

    builder = subprocess.run(
        [sys.executable, str(LANE / "verify_projected_source_rank.py")],
        cwd=LANE.parent,
        text=True,
        capture_output=True,
        check=False,
    )
    check(builder.returncode == 0, "repaired builder replay exits zero")
    check("SUMMARY 89/89" in builder.stdout, "repaired builder replay count")
    check("PROJECTED_OPERATOR_RANK 6" in builder.stdout,
          "repaired builder replay rank")

    print(f"SUMMARY {checks}/{checks} independent-audit seal checks passed")


if __name__ == "__main__":
    main()
