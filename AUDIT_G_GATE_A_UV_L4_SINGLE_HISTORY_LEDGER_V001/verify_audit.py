#!/usr/bin/env python3
"""Custody and anti-circularity checks for the L=4 ledger witness audit."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
TARGET = REPO / "DEVELOPMENT_G_GATE_A_UV_L4_SINGLE_HISTORY_LEDGER_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


target_theorem = TARGET / "THEOREM.md"
target_verifier = TARGET / "verify_single_history_ledger.py"
raw_audit = REPO / "AUDIT_G_GL6CY_RAW_SOURCE_COMPOSITION_V001/INDEPENDENT_RESULT.json"
gl6q_audit = REPO / "AUDIT_G_GL6Q_F3_LIFECYCLE_TO_INCIDENCE_RESPONSE_V001/AUDIT.md"
result_file = HERE / "INDEPENDENT_RESULT.json"
for path in (target_theorem, target_verifier, raw_audit, gl6q_audit, result_file):
    check(path.is_file(), f"required custody input {path.name}")

# Target hash pins are filled in before promotion.
EXPECTED = {
    "THEOREM.md": "113ca9798fe60a4afe7bada091d675ebb71608cab30f53b22bbc8ae59d10a06b",
    "verify_single_history_ledger.py": "02f405e695936ad6b20417dee4023897d1c8ab51157a11f52849d4aa259aaef0",
    "raw_source_audit": "df9547fb3ef36030026df47af6321e7a4082ed8e24f8b80d22666f74b575ca25",
    "gl6q_hostile_audit": "59dca3613f4f78c4bab00c239081e01184220c3975c02508a7910a2aa4ecd4d6",
}
actual = {
    "THEOREM.md": digest(target_theorem),
    "verify_single_history_ledger.py": digest(target_verifier),
    "raw_source_audit": digest(raw_audit),
    "gl6q_hostile_audit": digest(gl6q_audit),
}
for key, expected in EXPECTED.items():
    check(len(expected) == 64, f"frozen hash supplied for {key}")
    check(actual[key] == expected, f"frozen target hash for {key}")

theorem = target_theorem.read_text()
verifier = target_verifier.read_text()
raw = raw_audit.read_text()
raw_result = json.loads(raw)
result = json.loads(result_file.read_text())

for phrase in (
    "j_R(t;\\epsilon)={\\hbar r_0\\epsilon\\over\\tau}",
    "\\epsilon_\\star={\\pi\\over4r_0}",
    "W_R\n &:=\\int_0^\\tau",
    "not a parameter-free\nderivation",
    "does **not** close the plan-level Gate A",
):
    check(phrase in theorem, f"target scope and source map: {phrase}")

check("w_source = delta_q" not in verifier,
      "target does not circularly assign source write from retained charge")
check("w_source = (Fraction(1) - cos_two_phi_final) / 2" in verifier,
      "target computes source write from its own integral")
check("boundary_edge_fluxes = [Fraction(0) for _ in range(6)]" in verifier,
      "target declares all six boundary currents")
check(raw_result["disposition"] == "PASS" and
      raw_result["independent_checks"] == 983 and
      "H4=(32128/27)S; H6diag+cycles=(7212448/6075)S"
      in raw_result["raw_result"]["k0"],
      "audited raw source components are present")
check("PASS_AT_EXACT_CONDITIONAL_FINITE_MISSION_RESPONSE_SCOPE" in gl6q_audit.read_text(),
      "GL6Q hostile audit retains conditional mission scope")
check(result["disposition"] == "PASS_AT_CONDITIONAL_SINGLE_HISTORY_WITNESS_SCOPE" and
      result["ledger"]["balance"] == "0" and
      set(result["boundary_edges"].values()) == {"0"},
      "recorded independent result preserves exact scope and six-edge balance")

print("TARGET_HASHES", actual)
print(f"PASS__GATE_A_UV_L4_LEDGER_HOSTILE_CUSTODY__{checks}/{checks}")
