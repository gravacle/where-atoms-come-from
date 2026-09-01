#!/usr/bin/env python3
"""Verify frozen GL6AH custody, exact replays, and hostile scope."""

import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AH_BULK_BOUNDARY_IRREP_SEPARATION_V001"
THEOREM_SHA = "79b04596c6df950a86bfe25fb02f8cf7822d65f6eb6101615deba7d7f88b58eb"
MANIFEST_SHA = "f83a7524592be535a8e21e8fc94fdbcc0ba0f6c6d874779bf904769f3c115367"
SEAL_SHA = "5a94269aca3092c4be83c96738cf977e37369150ea7097a8aed6f000a0c61ef8"
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for ledger_name in ("AUDITED_TARGETS.sha256", "DEPENDENCIES.sha256"):
    for line in (HERE / ledger_name).read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = ROOT / relative
        check(path.is_file() and not path.is_symlink(), "regular " + relative)
        check(digest(path) == expected, "digest " + relative)

check(digest(TARGET / "THEOREM.md") == THEOREM_SHA, "theorem pin")
check(digest(TARGET / "MANIFEST.sha256") == MANIFEST_SHA, "manifest-file pin")
check(digest(TARGET / "SEAL.sha256") == SEAL_SHA, "seal-file pin")
check((HERE / "DEPENDENCIES.sha256").read_bytes()
      == (TARGET / "DEPENDENCIES.sha256").read_bytes(),
      "dependency ledger copied exactly")
author_rows = set((TARGET / "MANIFEST.sha256").read_text().splitlines())
audited_rows = set((HERE / "AUDITED_TARGETS.sha256").read_text().splitlines())
check(author_rows.issubset(audited_rows), "all author manifest members pinned")
check(any(row.endswith("/MANIFEST.sha256") for row in audited_rows),
      "author manifest file pinned")
check(any(row.endswith("/SEAL.sha256") for row in audited_rows),
      "author seal file pinned")
seal_expected, seal_relative = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal_expected == MANIFEST_SHA, "author seal content")
check(ROOT / seal_relative == TARGET / "MANIFEST.sha256", "author seal target")

independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "independent_gl6ah_replay.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS__INDEPENDENT_GL6AH_REPLAY__772/772" in independent.stdout,
      "independent replay marker")
for marker in (
    "DIRECT=FULL6_Q6_MINUS128_U_B_IN;A1_T2_EQUAL;E_NULL",
    "CHAIN=Q12_PLUS63371264_W_BD;Q16_MINUS123422773248_W_BD",
    "SUPPORT=TWO_BRIDGE_RECEIVER_HELPER_CHAIN_OWNS_Q12_Q16",
    "ENDPOINT=ETA_LITERAL_HELPER_COUNT;HOMOGENEOUS_ZERO;N1_SIGNS_RECONCILED",
    "DELTA_SCOPE=EXACT_DISPLAYED_WITNESS_DELTA0;OPEN_NONZERO_NEIGHBORHOOD_ONLY",
    "CEILING=NO_BULK_SHEAR_STATIONARY_MODE_COMMON_CONE_RICCI_GRAVITY_G",
):
    check(marker in independent.stdout, "independent marker " + marker)

author = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_mutable_packet.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
for marker in (
    "PASS__GL6AH_FROZEN_PACKET__131/131",
    "PASS__GL6AH_LOCAL_BULK_BOUNDARY__177/177",
    "PASS__GL6AH_N1_CONNECTOR_SUPPORTS__1293/1293",
    "DIRECT_Q6=-128_H4_UD2_U_B_IN;EXPECTATION=+8/45",
    "FULL_Q12=-63371264_W01;FULL_Q16=+123422773248_W01",
):
    check(marker in author.stdout, "author replay marker " + marker)

ledger = json.loads((TARGET / "EXACT_IRREP_LEDGER.json").read_text())
check(ledger["witness"]["delta"] == 0, "ledger delta witness")
check(ledger["direct_edge"]["raw_order_6"] == "-128 h^4 Ud^2 u_b^in",
      "ledger q6")
check(ledger["helper_chain"]["raw_E_order_12"] == "+63371264 w_bd",
      "ledger q12")
check(ledger["pair_source_chain"]["raw_E_order_16"] == "-123422773248 w_bd",
      "ledger q16")
check(ledger["N1_connector_support_census"]["all_other_supports"] == 0,
      "all other supports zero")
check(ledger["scope"] == {
    "physical_K_not_semantic_REC": True,
    "pair_mean_called_record": False,
    "term_ablation_called_physical_switch": False,
    "fixed_frame_E_only": True,
    "bulk_shear_claimed": False,
    "stationary_mode_claimed": False,
    "common_cone_claimed": False,
    "Ricci_or_gravity_or_G_claimed": False,
}, "ledger scope")
check(Fraction(128, 720) == Fraction(8, 45), "q6 factorial")
check(Fraction(63371264, 479001600) == Fraction(5626, 42525),
      "q12 factorial")
check(Fraction(123422773248, 20922789888000)
      == Fraction(1116019, 189189000), "q16 factorial")

source = (TARGET / "verify_n1_connector_supports.cpp").read_text()
for phrase in (
    "for (int bridges = 0; bridges < 64; ++bridges)",
    "bridge_mobius(matched_q6, support)",
    "bridge_mobius(matched_q12, support)",
    "bridge_mobius(source_pair_q16, support)",
    "inert-source-bridge q12 equality",
    "exactly two nonzero q16 supports",
):
    check(phrase in source, "author connector coverage " + phrase)

theorem = (TARGET / "THEOREM.md").read_text()
for phrase in (
    "author frozen after independent hostile pre-freeze review",
    "-128h^4U_d^2u_b^{\\rm in}",
    "i^6=-1",
    "u_b^{\\rm in}={1\\over2}\\mathbf1+{1\\over2}s_b",
    "E^Tx_q=0\\quad\\hbox{for every order }q",
    "not inferred from the `N=1` cell graph `K4`",
    "+63371264\\,w_{bd}",
    "-123422773248\\,w_{bd}",
    "only nonzero connector supports",
    "eta_{e_c}^{(0)}=\\sum_{d\\ne0,c}w_{0d}=-w_{0c}",
    "open detuning neighborhood",
    "does not extend\nthe exact `-128h^4U_d^2`",
    "not the semantic fact `REC`",
    "No stationary mode, full bulk shear kernel, common\ncone",
    "AUTHOR_FROZEN_POSTFREEZE_AUDIT_REQUIRED",
):
    check(phrase in theorem, "theorem clause " + phrase)

audit = (HERE / "AUDIT.md").read_text()
for phrase in (
    "No author byte was edited",
    "all sixty-four literal subsets",
    "remaining 5-, 9-, and 10-qubit Hamiltonians",
    "equal `A1` and `T2` coefficients",
    "not extrapolated from the `N=1` tetrahedron",
    "neither bridge alone owns either `E` term",
    "boundary ownership only for the displayed q12 helper term",
    "earns only an open neighborhood of\nnonvanishing",
    "not promoted to\nan invariant shear field",
    "**Audit verdict: PASS.**",
):
    check(phrase in audit, "audit clause " + phrase)

print(independent.stdout, end="")
print(author.stdout, end="")
print(f"PASS__GL6AH_POSTFREEZE_HOSTILE_AUDIT__{checks}/{checks}")
