#!/usr/bin/env python3
"""Verify frozen GL6AG custody, exact replays, and hostile scope."""

import hashlib
import json
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AG_N1_FORMATION_NEIGHBOR_PROPAGATION_V001"
THEOREM_SHA = "8551a4dc37183b8ab83ac48a0774dc0b82bd280faa5f9e469272f8c7ef898dea"
MANIFEST_SHA = "f69fd10402b1ec428b95f826ba19ec5fb9635a2079b05c0b49514b3084979b6e"
SEAL_SHA = "12b3e091624a8be4233b342a5303ec09439e140003ea5b770a7c7323e91d55c7"
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
check(author_rows.issubset(audited_rows), "every author manifest member pinned")
for special in ("MANIFEST.sha256", "SEAL.sha256"):
    check(any(row.endswith("/" + special) for row in audited_rows),
          "author " + special + " pinned")
seal_expected, seal_relative = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal_expected == MANIFEST_SHA, "author seal content")
check(ROOT / seal_relative == TARGET / "MANIFEST.sha256", "author seal target")

with tempfile.TemporaryDirectory(prefix="gl6ag-audit-") as temporary:
    temporary = Path(temporary)
    independent_binary = temporary / "independent_gl6ag_spot"
    author_binary = temporary / "verify_gl6ag_exact"
    compile_common = [
        "c++", "-O3", "-std=c++17", "-I/opt/homebrew/include",
        "-L/opt/homebrew/lib",
    ]
    subprocess.run(
        compile_common + [str(HERE / "independent_gl6ag_spot.cpp"),
                          "-lgmpxx", "-lgmp", "-o", str(independent_binary)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    independent = subprocess.run(
        [str(independent_binary)], cwd=ROOT, check=True,
        capture_output=True, text=True,
    )
    subprocess.run(
        compile_common + [str(TARGET / "verify_n1_matched_formation_propagation.cpp"),
                          "-lgmpxx", "-lgmp", "-o", str(author_binary)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    author_exact = subprocess.run(
        [str(author_binary)], cwd=ROOT, check=True,
        capture_output=True, text=True,
    )

check("PASS__INDEPENDENT_GL6AG_SPOT__2378/2378" in independent.stdout,
      "independent exact spot marker")
for marker in (
    "HAMILTONIAN=FULL16_24_WITHIN_6_SHARED_ONLY_SOURCE_K_X_SUPPORT_VARIES",
    "MATCHED=REFERENCE_0000_NO_ABSOLUTE_RECEIVER_INFERENCE",
    "E_TYPE=FIXED_BROKEN_S4_RESTRICTION_W_COLUMNS_IN_R2",
    "Q4=PLUS96_W01_COEFFICIENT_PLUS4",
    "Q12=MINUS63371264_KAPPA_C_W0C_COEFFICIENT_MINUS5626_OVER42525",
    "Q16=MOBIUS_PLUS123422773248_SUPPORTED_RECEIVER",
    "ABLATION=FOUR_CELL_FACTORIZATION_DIAGNOSTIC_NOT_PHYSICAL_SWITCH",
    "BRANCH=PHYSICAL_K_PROJECTORS_NOT_SEMANTIC_REC",
    "SCOPE=NO_BULK_CONE_STRESS_CONSERVATION_RICCI_GRAVITY_G",
):
    check(marker in independent.stdout, "independent clause " + marker)
check("PASS GL6AG exact matched-formation checks 3955/3955"
      in author_exact.stdout, "author full exact marker")

author_packet = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
author_structure = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_structure_and_ledger.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS GL6AG frozen packet checks 62/62" in author_packet.stdout,
      "author packet marker")
check("PASS GL6AG structure/ledger checks 144/144" in author_structure.stdout,
      "author structure marker")

source = (TARGET / "verify_n1_matched_formation_propagation.cpp").read_text()
for phrase in (
    "constexpr int kLinks = 16",
    "constexpr int kDimension = 1 << kLinks",
    "for (int mask = 0; mask < 16; ++mask)",
    "tables[mask][cell][coordinate][12] -",
    "Integer(63371264)",
    "tables[mask][cell][coordinate][16] -",
    "tables[1 << a][cell][coordinate][16] -",
    "Integer(123422773248)",
    "shared_bridges ? 30 : 24",
):
    check(phrase in source, "author exact coverage " + phrase)

ledger = json.loads((TARGET / "EXACT_MATCHED_LEDGER.json").read_text())
check(ledger["parent"] == {
    "N": 1,
    "active_links": 16,
    "within_cell_interactions": 24,
    "shared_child_bridges": 6,
    "total_interactions": 30,
    "witness": {
        "h": "E_0", "U_d": "E_0", "Delta": "6 E_0",
        "d_star": 2, "delta": 0,
    },
}, "ledger common parent")
check(ledger["matched_intervention"]["reference_source_pattern"] == "0000",
      "matched reference")
check(ledger["matched_intervention"]["absolute_neighbor_formula_used"] is False,
      "no absolute receiver formula")
check(ledger["embedding"]["w_ab_type"] ==
      "column vector E_(ab),:^T in R^2", "ledger column typing")
check(ledger["branch_replay"]["all_order_twelve_patterns"] == 16,
      "all sixteen patterns")
check(ledger["matched_receiver"]["signed_raw_order_12"] == -63371264,
      "ledger q12")
check(ledger["matched_receiver"]["order_12_nonzero_iff"] == "kappa_c = 1",
      "conditional q12 support")
check(ledger["remote_pair_mobius"]["signed_raw_order_16"] == 123422773248,
      "ledger q16")
check(ledger["bridge_off"]["authenticated_physical_switch_claimed"] is False,
      "diagnostic not physical switch")
check(Fraction(-63371264, 479001600) == Fraction(-5626, 42525),
      "q12 factorial reduction")
check(Fraction(123422773248, 20922789888000)
      == Fraction(1116019, 189189000), "q16 factorial reduction")

theorem = (TARGET / "THEOREM.md").read_text()
for phrase in (
    "all sixteen active links, twenty-four within-cell parent-clique\ninteractions, and all six shared-child bridges are retained",
    "The reference `widehat H_0` is obtained by setting only the four displayed\nsource transverse supports to zero",
    "Dynamics in (AG08) reads `P^K`, not the semantic terminal predicate\n`REC`",
    "No statement about the absolute receiver mean is used or\nneeded",
    "w_{ab}:=E_{(ab),:}^{T}\\in\\mathbb R^2",
    "only the declared fixed-frame\ntwo-coordinate restriction",
    "full sixteen-\npattern order-twelve census",
    "nonzero exactly when `\\kappa_c=1`",
    "genuine pair Möbius contrast",
    "term-ablation diagnostic proves",
    "not presented\nas an authenticated physical intervention",
    "no length, common cone,\nstress tensor, conservation law, continuum limit, Ricci response, gravity,\nor `G`",
):
    check(phrase in theorem, "frozen theorem clause " + phrase)

audit = (HERE / "AUDIT.md").read_text()
for phrase in (
    "No author byte was edited",
    "It does not read\nthe JSON response ledger",
    "reconstructs the diagonal from four occupied-link\ncounts",
    "only\nthe four source-cell transverse supports vary",
    "dynamic mediator is physical `K`",
    "does not infer\npropagation from an unsubtracted absolute receiver background",
    "correctly typed as\n`w_{ab}=E_{(ab),:}^T` in `R^2`",
    "makes no\ninvariant `E2`-sector claim",
    "order-twelve additivity is real but is not promoted to exact linearity",
    "term-ablation diagnostic",
    "does\nnot claim that any lifecycle controller physically switches those terms",
    "Ricci/Einstein form, gravity, or\nNewton's `G`",
    "tautological all-pattern\nassertion",
    "non-material metadata residue",
    "**Audit verdict: PASS.**",
):
    check(phrase in audit, "audit clause " + phrase)

# The stale labels are intentionally pinned and disclosed, not silently fixed.
check("AUTHOR_DRAFT_HOSTILE_REVIEW_REQUIRED_BEFORE_FREEZE" in theorem,
      "stale theorem label pinned")
check("This mutable author packet" in (TARGET / "README.md").read_text(),
      "stale README label pinned")
check("author frozen after independent hostile pre-freeze review" in theorem,
      "authoritative frozen front matter")

print(independent.stdout, end="")
print(author_exact.stdout, end="")
print(author_packet.stdout, end="")
print(author_structure.stdout, end="")
print(f"PASS__GL6AG_POSTFREEZE_HOSTILE_AUDIT__{checks}/{checks}")
