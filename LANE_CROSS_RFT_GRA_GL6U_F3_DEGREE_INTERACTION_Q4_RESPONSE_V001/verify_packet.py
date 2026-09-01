#!/usr/bin/env python3
"""Structural and custody checks for the frozen GL6U author packet."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
checks: list[str] = []


def check(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    checks.append(label)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


required = [
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "INTERACTING_RESPONSE_LEDGER.json",
    "verify_degree_interaction_q4_response.py",
    "verify_packet.py",
    "DEPENDENCIES.sha256",
    "VERIFICATION.txt",
]
for name in required:
    check(f"required {name}", (HERE / name).is_file())

theorem = (HERE / "THEOREM.md").read_text()
result = (HERE / "RESULT.md").read_text()
self_audit = (HERE / "SELF_AUDIT.md").read_text()
ledger = json.loads((HERE / "INTERACTING_RESPONSE_LEDGER.json").read_text())

for phrase in (
    "delta_r=r\\Delta+U_dr(r+1-4d_\\star)",
    "D^{\\rm BREAK}(\\tau)=0",
    "D^{\\rm KEEP}=-8hxI_6-4hyA_L",
    "D_{E_2}=-8h(x-y)",
    "\\mathfrak C_{XZZ}:=y-xz^2",
    "-{16\\over3}h^3U_ds^4",
    "2\\delta_1-\\delta_2=-2U_d",
    "not called a full connected cumulant",
    "small punctured positive-prewait interval",
    "curvature template",
    "Ricci/Einstein is earned only if",
):
    check(f"theorem phrase {phrase}", phrase in theorem)

check("result interaction defect", "factorization defect" in result)
check("result operator neutral", "no Ricci form is assumed" in result)
check("self audit inherited interaction", "inherited BS06/FJ degree term" in self_audit)
check("native degree term", "degree" in ledger["restored_native_term"])
check("E2 ledger", ledger["s4_eigenvalues"]["E2"] == "-8*h*(x-y)")
check("no G", "no_G" in ledger["ceilings"])

for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, rel = line.split("  ", 1)
    target = REPO / rel
    check(f"dependency exists {rel}", target.is_file())
    check(f"dependency hash {rel}", sha256(target) == expected)

print(f"PASS__GL6U_PACKET__{len(checks)}/{len(checks)}")
