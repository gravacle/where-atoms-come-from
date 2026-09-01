#!/usr/bin/env python3
"""Verify frozen GL6AI custody, independent replay, and hostile scope."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001"
THEOREM_SHA = "a51e802f6ba148e5f9848e95f41a80073795b24b7eaf87e36c0766b0856aa494"
MANIFEST_SHA = "fc50cad54dca00aab1c30d7c12ef07147df1242f94483f63955185695073f706"
SEAL_SHA = "12daa03d45cd653db24622ae8b3d8166015291534b3295e3c426eb37180fc918"
AH_THEOREM_SHA = "79b04596c6df950a86bfe25fb02f8cf7822d65f6eb6101615deba7d7f88b58eb"
AH_SEAL_SHA = "5a94269aca3092c4be83c96738cf977e37369150ea7097a8aed6f000a0c61ef8"
AH_AUDIT_SHA = "c393b2c98304d3b83dd4b02d5de0cff8f4879001644590365bad44c17a3122f0"
AH_AUDIT_SEAL_SHA = "0915014ca7cf8468225c42d65e467ec3ad8eb5950e229e6c2dcf95dc9c897a26"
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
check(
    (HERE / "DEPENDENCIES.sha256").read_bytes()
    == (TARGET / "DEPENDENCIES.sha256").read_bytes(),
    "dependency ledger copied exactly",
)
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

dependencies = dict(
    (relative, expected)
    for expected, relative in (
        line.split("  ", 1)
        for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    )
)
check(dependencies[
    "LANE_CROSS_RFT_GRA_GL6AH_BULK_BOUNDARY_IRREP_SEPARATION_V001/THEOREM.md"
] == AH_THEOREM_SHA, "GL6AH theorem dependency")
check(dependencies[
    "LANE_CROSS_RFT_GRA_GL6AH_BULK_BOUNDARY_IRREP_SEPARATION_V001/SEAL.sha256"
] == AH_SEAL_SHA, "GL6AH seal dependency")
check(dependencies[
    "AUDIT_G_GL6AH_BULK_BOUNDARY_IRREP_SEPARATION_V001/AUDIT.md"
] == AH_AUDIT_SHA, "GL6AH audit dependency")
check(dependencies[
    "AUDIT_G_GL6AH_BULK_BOUNDARY_IRREP_SEPARATION_V001/SEAL.sha256"
] == AH_AUDIT_SEAL_SHA, "GL6AH audit seal dependency")

independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "independent_gl6ai_replay.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
for marker in (
    "PASS__INDEPENDENT_GL6AI_REPLAY__",
    "SPLIT=DEGREE_SQUARE_ONSITE_PLUS_ONCE_COUNTED_2UD_PAIR_EXACT",
    "TOPOLOGY=LINK_DEGREE_Q_CHILD_PLUS2_LE6;D_LINK_GE_D_CELL_EXACT",
    "INFLUENCE=DIAGONAL_DRESSING_PLUS_OFFDIAGONAL_ROW_2JDEG;LAMBDA_48UD_OVER_HBAR",
    "SOURCE=NORMALIZED_BETA1_MINUS_BETA0;V_S_MINUS_H_X_S;DUHAMEL_2H_OVER_HBAR",
    "TAIL=AI18_COEFFICIENTWISE_EXACT;MARKED_EXPONENTIAL;FINITE_Y_ONLY_CARDINALITY",
    "CEILING=ANALYTIC_QUASILOCAL_ENVELOPE_NOT_EXACT_SPEED_LORENTZ_RICCI_GRAVITY_G",
):
    check(marker in independent.stdout, "independent marker " + marker)

author_packet = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS__GL6AI_PACKET__47/47" in author_packet.stdout, "author packet replay")
author_exact = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_f3_relational_envelope.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS 123777/123777" in author_exact.stdout, "author exact replay")

ledger = json.loads((TARGET / "ENVELOPE_LEDGER.json").read_text())
check(ledger["exact_split"]["pair"] == "sum_{e~f} 2 Ud n_e n_f", "ledger pair")
check(ledger["topology"]["uniform_link_degree"] == 6, "ledger degree")
check(ledger["influence_matrix"]["row_sum_bound"] == "2 J Delta_L", "ledger row")
check(
    ledger["certified_constants"]["lambda_F3"]
    == "4 J Delta_L/hbar = 48 |Ud|/hbar",
    "ledger lambda",
)
check("beta_s=1 minus beta_s=0" in ledger["source_contrast"], "ledger source typing")
check("not exact finite-speed support" in ledger["ceiling"], "ledger analytic ceiling")

theorem = (TARGET / "THEOREM.md").read_text()
for phrase in (
    "author frozen after independent hostile pre-freeze review",
    r"\varepsilon_\star=\Delta+2U_d(1-2d_\star)",
    r"H_{\rm pair}=\sum_{\{e,f\}\in E(L_N)}2U_dn_en_f",
    r"\deg_{L_N}(e)=3+[q(c)-1]=q(c)+2\le6",
    r"J\deg_{L_N}(e),&e=f",
    r"\|\mathsf J\|_\infty\le2J\Delta_L",
    r"={48|U_d|\over\hbar}",
    "No sharper constant is claimed",
    r"V_s=-hX_s,qquad\|V_s\|=h",
    r"beta_s=1,\beta_{-s}",
    r"beta_s=0,\beta_{-s}",
    "not\nsuccess filtering and not a formal deletion switch",
    r"{2h\|B_Y\|\over\hbar}",
    r"d_L(e,f)\ge d_{\mathcal G_N}(m,n)",
    r"d_{\rm cell}(s,Y):=",
    r"T_d(x)\le e^{xe^\mu-\mu d}",
    "At every nonzero time\nthe analytic tail",
    "not a Lorentz cone",
    "Ricci through the\nconditional `GL6L` bridge",
    "AUTHOR_FROZEN_POSTFREEZE_AUDIT_REQUIRED",
):
    check(phrase in theorem, "theorem clause " + phrase)

audit = (HERE / "AUDIT.md").read_text()
for phrase in (
    "No\nauthor byte was edited",
    "one advancing and one non-advancing channel",
    "No\nsupport-cardinality factor is hidden in this rate",
    "normalized `beta_s=1` minus `beta_s=0`",
    "projects to a genuine cell walk",
    "analytic tail is generally nonzero",
    "missing backslash before `qquad`",
    "No Lorentz or common physical cone",
    "**Audit verdict: PASS.**",
):
    check(phrase in audit, "audit clause " + phrase)

print(independent.stdout, end="")
print(author_packet.stdout, end="")
print(author_exact.stdout, end="")
print(f"PASS__GL6AI_POSTFREEZE_HOSTILE_AUDIT__{checks}/{checks}")
