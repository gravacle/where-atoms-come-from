#!/usr/bin/env python3
"""Structural and custody checks for the frozen GL6T author packet."""

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
    "LINEAGE_Q4_RESPONSE_LEDGER.json",
    "verify_lineage_gated_q4_e2_response.py",
    "verify_packet.py",
    "DEPENDENCIES.sha256",
    "VERIFICATION.txt",
]
for name in required:
    check(f"required {name}", (HERE / name).is_file())

theorem = (HERE / "THEOREM.md").read_text()
result = (HERE / "RESULT.md").read_text()
self_audit = (HERE / "SELF_AUDIT.md").read_text()
ledger = json.loads((HERE / "LINEAGE_Q4_RESPONSE_LEDGER.json").read_text())

for phrase in (
    "Specialize FPSS to `N=0`",
    "twelve explicit",
    "not a success-filtered unconditional experiment",
    "omit both the optional FPSS saturation pulse `U_KX`",
    "H_{\\rm full}=H_\\star",
    "Q^{K_a}_{F_a,u_{-a},{\\rm KEEP}}",
    "D^{\\rm KEEP}(\\tau)",
    "-8hx\\,I_6-4hxz^2A_L",
    "D_{E_2}&=-8hx(1-z^2)",
    "D^{\\rm BREAK}(\\tau)=0",
    "does not prove an inter-link interaction",
    "Nor is this a provenance-only force",
    "not called positive Ricci stiffness",
):
    check(f"theorem phrase {phrase}", phrase in theorem)

check("result factorized ceiling", "four link evolutions factorize" in result)
check("result provenance ceiling", "unqualified occupied `K`" in result)
check("self audit per-link record", "other three event" in self_audit)
check("four keep records", ledger["matched_routes"]["KEEP"]["K"] == 4)
check("four break quarantine", ledger["matched_routes"]["BREAK"]["G"] == 4)
check("N0 parent", "FPSS_N0" in ledger["selected_parent"])
check("complete instrument", "no_success_filter" in ledger["instrument_scope"])
check("factorized ceiling", "no_interlink_interaction" in ledger["dynamical_scope"])
check("E2 sector", ledger["s4_eigenvalues"]["E2"] == "-8*h*x*(1-z^2)")
check("no G ceiling", "no_G" in ledger["ceilings"])

for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, rel = line.split("  ", 1)
    target = REPO / rel
    check(f"dependency exists {rel}", target.is_file())
    check(f"dependency hash {rel}", sha256(target) == expected)

print(f"PASS__GL6T_PACKET__{len(checks)}/{len(checks)}")
