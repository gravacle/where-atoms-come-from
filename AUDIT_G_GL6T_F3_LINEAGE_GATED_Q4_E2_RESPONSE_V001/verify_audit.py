#!/usr/bin/env python3
"""Verify frozen custody, independent algebra, and scope for GL6T audit."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6T_F3_LINEAGE_GATED_Q4_E2_RESPONSE_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


for ledger in ("AUDITED_TARGETS.sha256", "DEPENDENCIES.sha256"):
    for line in (HERE / ledger).read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = ROOT / relative
        check(path.is_file() and not path.is_symlink(), f"regular custody {relative}")
        check(digest(path) == expected, f"digest custody {relative}")

run = subprocess.run(
    [sys.executable, "-B", str(HERE / "independent_gl6t_replay.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS__INDEPENDENT_GL6T_REPLAY__" in run.stdout,
      "independent exact replay marker")

frozen = (TARGET / "VERIFICATION.txt").read_text()
check("PASS__GL6T_LINEAGE_GATED_Q4_E2_RESPONSE__70/70" in frozen,
      "frozen author replay marker")
check("PASS__GL6T_PACKET__62/62" in frozen, "frozen author packet marker")

theorem = (TARGET / "THEOREM.md").read_text()
for phrase in (
    "Specialize FPSS to `N=0`",
    "twelve explicit\nnonedges",
    "not a success-filtered unconditional experiment",
    "omit both the optional FPSS saturation pulse `U_KX`",
    "for every `a` and every fixed vector of the other\nthree event alternatives",
    "No such product is asserted after mixing different\ninstrument branches",
    "D^{\\rm KEEP}(\\tau)",
    "D_{E_2}&=-8hx(1-z^2)",
    "D^{\\rm BREAK}(\\tau)=0",
    "multi-link **pair-coordinate response**",
    "this theorem does not prove an inter-link interaction, collective phase,",
    "unqualified occupied register in the same state would reproduce (T11)",
    "not called positive Ricci stiffness",
):
    check(phrase in theorem, f"frozen theorem clause {phrase}")

audit = (HERE / "AUDIT.md").read_text()
for phrase in (
    "PASS_AT_EXACT_FINITE_BRANCHWISE_PAIR_RESPONSE_CAPACITY_SCOPE",
    "not a success-filtered unconditional law",
    "without consuming the active-link blank",
    "does not assert that the active pair observable is a record",
    "not an inter-link interaction, collective phase, feedback law, or stiffness",
    "The six pair sources remain formal and unnormalized",
    "**Audit verdict: PASS.**",
):
    check(phrase in audit, f"audit ceiling clause {phrase}")

print(run.stdout, end="")
print(f"PASS__GL6T_HOSTILE_AUDIT__{checks}/{checks}")
