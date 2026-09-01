#!/usr/bin/env python3
"""Verify frozen GL6AA custody, independent replay, and claim scope."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001"
THEOREM_SHA = "faea49e3dcd5f2b4d5b3bab9026432d741192339dd789a939ef2318236848c0e"
MANIFEST_FILE_SHA = "9bebed96f2b864738639571f870ae7a34dc440dea5ab4418bc4e7cc9c2eb2a63"
SEAL_FILE_SHA = "e7d6b64a713fa37170fa70ecd7978a93028fc1594b4e2bd6e8befffad61d3004"
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
        check(path.is_file() and not path.is_symlink(),
              "regular custody " + relative)
        check(digest(path) == expected, "digest custody " + relative)

check(digest(TARGET / "THEOREM.md") == THEOREM_SHA,
      "explicit frozen theorem digest")
check(digest(TARGET / "MANIFEST.sha256") == MANIFEST_FILE_SHA,
      "explicit frozen manifest-file digest")
check(digest(TARGET / "SEAL.sha256") == SEAL_FILE_SHA,
      "explicit frozen seal-file digest")
check((HERE / "DEPENDENCIES.sha256").read_bytes()
      == (TARGET / "DEPENDENCIES.sha256").read_bytes(),
      "audit dependency ledger equals author ledger")

author_rows = set((TARGET / "MANIFEST.sha256").read_text().splitlines())
audited_rows = set((HERE / "AUDITED_TARGETS.sha256").read_text().splitlines())
check(author_rows.issubset(audited_rows),
      "audit target ledger covers every author manifest row")
check(any(row.startswith(MANIFEST_FILE_SHA) and row.endswith(
    "LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001/MANIFEST.sha256")
    for row in audited_rows), "audit target ledger pins author manifest")
check(any(row.startswith(SEAL_FILE_SHA) and row.endswith(
    "LANE_CROSS_RFT_GRA_GL6AA_RECORD_AUTHENTICATED_SHARED_CHILD_ATLAS_V001/SEAL.sha256")
    for row in audited_rows), "audit target ledger pins author seal")

seal_expected, seal_relative = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal_expected == MANIFEST_FILE_SHA, "author seal content digest")
check(ROOT / seal_relative == TARGET / "MANIFEST.sha256", "author seal target")

replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "independent_gl6aa_replay.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS__INDEPENDENT_GL6AA_REPLAY__207446/207446" in replay.stdout,
      "independent replay marker")
for marker in (
    "BOUNDARY=ONE_PORT_SIX_EXTENSIONS_NONTRIVIAL_HOLONOMY",
    "COPY=DIRECT_SUM_BLANK_CODE_INVOLUTION_DISJOINT_TAPS",
    "QUERY=COMPLETE_NONEXCLUSIVE_FLAGS_MATCH_IFF_ALL_FALSE",
    "SHARED_CHILD=EQUAL_ID_IFF_LITERAL_CHILD_ON_MATCH",
    "ATLAS=TRANSLATION_AND_S4_SIX_PAIR_COCYCLES_EXACT",
    "SCOPE=SELECTED_RELATIONAL_INCIDENCE_NOT_SPACE_CONE_RICCI_GRAVITY_G",
):
    check(marker in replay.stdout, "independent replay clause " + marker)

# The frozen author executables are run only after the independent replay.
author_exact = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_record_authenticated_atlas.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
author_packet = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=ROOT, check=True, capture_output=True, text=True,
)
check("PASS GL6AA exact checks 1686208/1686208" in author_exact.stdout,
      "live frozen author replay marker")
check("PASS GL6AA packet checks 57/57" in author_packet.stdout,
      "live frozen author packet marker")

frozen = (TARGET / "VERIFICATION.txt").read_text()
check("PASS GL6AA exact checks 1686208/1686208" in frozen,
      "frozen author exact transcript")
check("PASS GL6AA packet checks 37/37" in frozen,
      "frozen stale author packet transcript is pinned")
check("Independent pre-freeze hostile verdict: CLEAN" in frozen,
      "pre-freeze hostile review marker")

ledger = json.loads((TARGET / "ATLAS_LEDGER.json").read_text())
check(ledger["completion_type"]
      == "selected record-native physical port completion",
      "ledger selected completion typing")
check(ledger["new_bulk_interactions"] == 0, "ledger no bulk interaction")
check("forbidden to read p_N" in ledger["actual_ancestry"],
      "ledger actual ancestry independence")
check(ledger["expected_ancestry"] == "retained FPSS program tuple record",
      "ledger expected ancestry")
check("MATCH iff all false" in ledger["terminal_output"],
      "ledger complete flag output")
for proved in (
    "single-port response does not identify a full S4 transition",
    "equal tapped child IDs iff same literal child on MATCH",
    "translation inverse and all-cycle cocycle",
    "complete four-port and induced six-pair transition cocycle",
    "operational indexing of GL6Z on selected relational atlas",
):
    check(proved in ledger["proved"], "ledger proved clause " + proved)
for open_clause in (
    "autonomous support and address selection",
    "physical length, momentum, and proper time",
    "E2 propagation",
    "common causal cone and refinement",
    "infrared operator, gravity, and G",
):
    check(open_clause in ledger["open"], "ledger open clause " + open_clause)

theorem = (TARGET / "THEOREM.md").read_text()
for phrase in (
    "author-frozen after independent hostile pre-freeze review",
    "giving `3!=6` extensions",
    "response bytes alone",
    "No factor in that source/controller path\nmay read `p_N`",
    "will be a scored `MATCH` result,\nnot the preparation rule",
    "exact direct-sum copy dilation",
    "makes no universal\ncloning claim",
    "conditional finite derivative\nrecords",
    "flags need not be mutually\nexclusive",
    "`K_e=1` **iff**",
    "postprocessing of a complete query, not success filtering",
    "queried injectivity of the independent physical-site source `sigma`",
    "g_{nm}=\\lambda_n^{-1}\\lambda_m",
    "T^{(6)}_{\\ell n}T^{(6)}_{nm}=T^{(6)}_{\\ell m}",
    "response is not\npostselected on `MATCH`",
    "ID and label banks are therefore spectators",
    "dual physical\nwavevector remain calibration data",
    "infrared operator content without a\nRicci ansatz",
):
    check(phrase in theorem, "frozen theorem clause " + phrase)

audit = (HERE / "AUDIT.md").read_text()
for phrase in (
    "PASS__OLD_PAIRWISE_QUERY_NONIDENTIFYING",
    "This is a stale transcript count, not a physics or custody\nfailure",
    "exact information boundary",
    "not used to prepare `sigma`",
    "makes no universal-cloning\nclaim",
    "conditional finite\nderivative record",
    "Dynamic `n_e` on expected active edges remains read but unconstrained",
    "not inferred from equality of response amplitudes",
    "graph adjacency on\nthe selected incidence atlas; it is not an SI length",
    "exhausts all\n`24^3` triples",
    "response is not postselected on MATCH",
    "does not derive autonomous address/support selection",
    "**Audit verdict: PASS.**",
):
    check(phrase in audit, "audit scope clause " + phrase)

print(replay.stdout, end="")
print(author_exact.stdout, end="")
print(author_packet.stdout, end="")
print("PASS__GL6AA_HOSTILE_AUDIT__%d/%d" % (checks, checks))
