#!/usr/bin/env python3
"""Verify frozen GL6AF custody, independent replay, and hostile scope."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001"
THEOREM_SHA = "6e20f188aa84a369d90eb67d741eb0c8aa0eb89013598b57021fa7464249fc71"
MANIFEST_SHA = "f2f4b416db1caa3a935e641cf145a488728b7281325a0e2a4861d6181e48677f"
SEAL_SHA = "7d440fb501063acd960658bb5318f46483568599e5a272669e48b1992c8a9f86"
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

check(digest(TARGET / "THEOREM.md") == THEOREM_SHA, "frozen theorem pin")
check(digest(TARGET / "MANIFEST.sha256") == MANIFEST_SHA,
      "frozen manifest-file pin")
check(digest(TARGET / "SEAL.sha256") == SEAL_SHA, "frozen seal-file pin")
check((HERE / "DEPENDENCIES.sha256").read_bytes()
      == (TARGET / "DEPENDENCIES.sha256").read_bytes(),
      "audit dependencies equal author dependencies")

author_rows = set((TARGET / "MANIFEST.sha256").read_text().splitlines())
audited_rows = set((HERE / "AUDITED_TARGETS.sha256").read_text().splitlines())
check(author_rows.issubset(audited_rows), "every frozen author byte pinned")
seal_expected, seal_relative = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal_expected == MANIFEST_SHA, "author seal content")
check(ROOT / seal_relative == TARGET / "MANIFEST.sha256", "author seal target")

independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "independent_gl6af_replay.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS__INDEPENDENT_GL6AF_REPLAY__22035/22035" in independent.stdout,
      "independent exact marker")
for marker in (
    "PAULI=ALL_16_PATTERNS_COMPLETE_6X6_DOUBLE_COMMUTATOR_EXACT",
    "RANKS=FORMED_COUNT_0_1_TO0_2_TO1_3_4_TO2",
    "COEFFICIENT=TWO_FORMED_MINUS16_H_X_ONE_MINUS_Z_WTW",
    "MIXING=BROKEN_S4_FIXED_E_RESTRICTION_NOT_INVARIANT_BLOCK",
    "BRANCH=ORTHOGONAL_K_PROJECTORS_H_M_GL6V_SOURCE_PRESERVE_SECTORS",
    "CTP=BRANCHWISE_SLOPE_PLUS8_ESTAR2_H_X_ONE_MINUS_Z_OVER_HBAR2",
    "ANCESTRY=PHYSICAL_K_MEDIATOR_NOT_SEMANTIC_REC_OR_PAIR_RECORD",
    "SCOPE=NO_COLLECTIVE_STIFFNESS_BULK_STRESS_RICCI_GRAVITY_G",
):
    check(marker in independent.stdout, "independent clause " + marker)

author_exact = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_formation_pattern_e2_source.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
author_packet = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS GL6AF exact formation-pattern checks 33/33" in author_exact.stdout,
      "author exact marker")
check("PASS GL6AF packet checks 72/72" in author_packet.stdout,
      "author packet marker")

ledger = json.loads((TARGET / "FORMATION_SOURCE_LEDGER.json").read_text())
check(ledger["rank_by_formed_count"] == {
    "0": 0, "1": 0, "2": 1, "3": 2, "4": 2,
}, "ledger ranks")
check(ledger["two_formed_E_restriction"]
      == "-16 h x (1-z) transpose(w_ab) w_ab", "ledger coefficient")
check(ledger["normalized_ctp_slope"]
      == "8 E_star^2 h x (1-z)/hbar^2 transpose(w_ab) w_ab",
      "ledger normalized slope")
for scope in (
    "source-independent conserved branch projector",
    "fixed program-frame E x E restriction, not invariant E2 block for broken S4",
    "physical K support mediates response; semantic REC is not dynamically read",
    "no collective stiffness, bulk source, stress, Ricci, gravity, or G",
):
    check(scope in ledger["scope"], "ledger scope " + scope)

theorem = (TARGET / "THEOREM.md").read_text()
for phrase in (
    "physical branch restriction of one operator",
    "retained-pattern projector is",
    "support projectors commute with the response Hamiltonian",
    "They also commute with the GL6V source/read dilation",
    "probability of a retained `kappa` branch is independent",
    "restricted response bilinear",
    "rank is\ntherefore basis independent",
    "need not preserve the `E` plane",
    "not call that restriction an invariant `E2` Schur block",
    "-16hx(1-z)",
    "three formed lineages and four formed lineages have\nrank-two `E x E` restrictions",
    "Physical normalized entrance slope",
    "Its sign and normalization are inherited",
    "not the semantic\nfact that a later query certifies `REC`",
    "otherwise supplied `K=1` support\nsector would reproduce",
    "pair-coordinate response, not collective stiffness",
    "NO_COLLECTIVE_STIFFNESS_BULK_STRESS_RICCI_GRAVITY_OR_G_CLAIM",
):
    check(phrase in theorem, "frozen theorem clause " + phrase)

audit = (HERE / "AUDIT.md").read_text()
for phrase in (
    "imports no author algebra",
    "all\nthirty-six pair-response entries",
    "one Hamiltonian, not a branch-labeled coupling law",
    "The replay checks the operator and expectation stages separately",
    "six patterns use the three covector lines",
    "partial branch mixes the fixed `E` plane with its complement",
    "mutually orthogonal, idempotent, and resolve",
    "source cannot change the branch probability",
    "including the sign, factor of one half",
    "not an authentication-sensitive force",
    "not an inter-link interaction, collective feedback, or\nstiffness",
    "**Audit verdict: PASS.**",
):
    check(phrase in audit, "audit clause " + phrase)

print(independent.stdout, end="")
print(author_exact.stdout, end="")
print(author_packet.stdout, end="")
print(f"PASS__GL6AF_HOSTILE_AUDIT__{checks}/{checks}")
