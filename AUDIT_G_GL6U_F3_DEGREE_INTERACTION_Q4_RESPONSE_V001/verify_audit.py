#!/usr/bin/env python3
"""Verify frozen custody, independent algebra, and scope for GL6U audit."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6U_F3_DEGREE_INTERACTION_Q4_RESPONSE_V001"
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
    [sys.executable, "-B", str(HERE / "independent_gl6u_replay.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS__INDEPENDENT_GL6U_REPLAY__261/261" in run.stdout,
      "independent exact replay marker")

frozen = (TARGET / "VERIFICATION.txt").read_text()
check("PASS__GL6U_DEGREE_INTERACTION_Q4_RESPONSE__72/72" in frozen,
      "frozen author replay marker")
check("PASS__GL6U_PACKET__50/50" in frozen, "frozen author packet marker")

theorem = (TARGET / "THEOREM.md").read_text()
for phrase in (
    "This packet restores (U02) on the same parent.  No Hamiltonian term is added.",
    "The unique active FPSS parent has degree `d`; its four children have degrees",
    "the other three parent sites and all twelve raw-domain nonedges",
    "\\delta_r=r\\Delta+U_dr(r+1-4d_\\star)",
    "D^{\\rm BREAK}(\\tau)=0",
    "D^{\\rm KEEP}=-8hxI_6-4hyA_L",
    "D_{E_2}=-8h(x-y)",
    "\\mathfrak C_{XZZ}:=y-xz^2",
    "-{16\\over3}h^3U_ds^4",
    "not called a full connected cumulant",
    "An unqualified active `K` word in\nthe same physical state would reproduce it",
    "not yet a collective thermodynamic phase",
):
    check(phrase in theorem, f"frozen theorem clause {phrase}")

audit = (HERE / "AUDIT.md").read_text()
for phrase in (
    "PASS_AT_EXACT_FINITE_INTERACTION_OWNED_FACTORIZATION_DEFECT_AND_FULL_PAIR_RESPONSE_SCOPE",
    "one active parent, three guard parents",
    "No child, guard, or\nnonedge energy has been silently discarded",
    "checks all thirty-six double\ncommutators",
    "same operator law",
    "interaction-owned inter-link factorization defect",
    "An unqualified active `K` word in the same physical\nstate would reproduce the response",
    "not a normalized physical\nsource/read CTP Hessian",
    "**Audit verdict: PASS.**",
):
    check(phrase in audit, f"audit clause {phrase}")

print(run.stdout, end="")
print(f"PASS__GL6U_HOSTILE_AUDIT__{checks}/{checks}")
