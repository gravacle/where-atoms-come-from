#!/usr/bin/env python3
"""Fail-closed custody, exact-replay, and scope verifier for GL6CN."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(label)


def strict_json(path):
    def no_duplicates(pairs):
        answer = {}
        for key, value in pairs:
            if key in answer:
                raise ValueError(f"duplicate JSON key: {key}")
            answer[key] = value
        return answer

    return json.loads(path.read_text(), object_pairs_hook=no_duplicates)


required = (
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "EXACT_LEDGER.json",
    "VERIFICATION.txt",
    "derive_complete_diagonal_h6_t2_source.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)
for name in required:
    check((HERE / name).is_file(), f"required file: {name}")


# Dependency custody.
dependency_paths = set()
for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"dependency hash: {relative}")
    check(relative not in dependency_paths, f"dependency unique: {relative}")
    dependency_paths.add(relative)
check(len(dependency_paths) == 18, "exact dependency count")
check(all(any(label in path for label in ("GL6AO", "GL6CG", "GL6CH"))
          for path in dependency_paths), "only declared lanes imported")
for label in ("GL6AO", "GL6CG", "GL6CH"):
    check(sum(label in path for path in dependency_paths) == 6,
          f"six author/audit custody rows for {label}")
    check(any(label in path and path.endswith("/THEOREM.md")
              for path in dependency_paths), f"{label} theorem pinned")
    check(any(path.startswith("AUDIT_G_") and label in path and
              (path.endswith("/AUDIT.md") or path.endswith("/AUDIT_REPORT.md"))
              for path in dependency_paths), f"{label} audit report pinned")

for relative in sorted(path for path in dependency_paths
                       if path.endswith("/SEAL.sha256")):
    seal_path = ROOT / relative
    seal_lines = [row for row in seal_path.read_text().splitlines() if row.strip()]
    check(len(seal_lines) == 1, f"one upstream seal row: {relative}")
    expected, manifest_relative = seal_lines[0].split(maxsplit=1)
    manifest_path = ROOT / manifest_relative
    if not manifest_path.is_file():
        manifest_path = seal_path.parent / manifest_relative
    check(manifest_path.is_file(), f"upstream seal target exists: {relative}")
    check(hashlib.sha256(manifest_path.read_bytes()).hexdigest() == expected,
          f"upstream seal closes: {relative}")

for audit_relative in (
    "AUDIT_G_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/AUDIT.md",
    "AUDIT_G_GL6CG_DENSE_PARENT_GLOBAL_H4_SOURCE_OPERATOR_V001/AUDIT_REPORT.md",
    "AUDIT_G_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/AUDIT_REPORT.md",
):
    audit_text = (ROOT / audit_relative).read_text()
    check("PASS" in audit_text, f"upstream independent audit passes: {audit_relative}")


# Packet manifest and seal.
manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"manifest target exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"manifest hash: {relative}")
    check(relative not in manifest_paths, f"manifest row unique: {relative}")
    manifest_paths.add(relative)
for name in required:
    if name not in ("MANIFEST.sha256", "SEAL.sha256"):
        check(f"{HERE.name}/{name}" in manifest_paths,
              f"manifest coverage: {name}")
check(len(manifest_paths) == 10, "exact manifest count")

seal_lines = [line for line in (HERE / "SEAL.sha256").read_text().splitlines()
              if line.strip()]
check(len(seal_lines) == 1, "one seal row")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "seal targets manifest")
check(hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
      == expected, "seal hash")


# Frozen exact ledger.
ledger = strict_json(HERE / "EXACT_LEDGER.json")
check(ledger["schema"] == "GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE_V001",
      "ledger schema")
check(ledger["checks"] == 10775, "frozen exact-check count")
check(ledger["kato_formula"] ==
      "K6=T6-bX4+b^2A3-dA2; b=-M/2; d=-7M/24", "Kato formula")
check(ledger["differentiated_formula"] ==
      "K6_T'=T6_T'-bX4_T'+b^2A3_T'-dA2_T'", "differentiated formula")
check(ledger["complete_diagonal_h6_T2_vertex"] == "zero pointwise",
      "main zero result")
check(ledger["all_complete_t_contractions_zero"] is True,
      "all complete literal contractions vanish")
check(ledger["all_individual_kato_terms_t_contractions_zero"] is True,
      "all individual literal contractions vanish")
check(ledger["literal_Q4_case_count"] == 9, "nine literal Q4 cases")
check(ledger["row_sha256"] ==
      "beb0b2a54528b3d3f1cfdbb698670889f39c1c5b1b0670e93e3ee907f78f80a5",
      "frozen literal-row digest")

dual = ledger["dual_number_crosscheck"]
check(dual["resolvent_rule"] ==
      "d[-1/(E+eta m)]/deta at zero = +m/E^2", "resolvent derivative sign")
check(dual["triple_signatures_checked"] == 22,
      "all labelled triple signatures dual checked")
check(dual["A2_prime_single_score"] == "-1/4", "A2 derivative sign")
check(dual["A3_prime_single_score"] == "3/16", "A3 derivative sign")
check(dual["triple_rows_sha256"] ==
      "7713d0a7b35134e381864dd744fa69ed93dc6a79914f301806e0342163263b97",
      "frozen dual triple digest")

pair = ledger["history_derivative_kernels"]["pair"]
check(pair["2"]["T6"] == {"1": "3/16", "2": "3/16", "3": "1/4"},
      "pair-energy-two T6 kernel")
check(pair["4"]["T6"] == {"1": "3/64", "2": "3/64", "3": "1/32"},
      "pair-energy-four T6 kernel")
check(pair["6"]["T6"] == {"1": "1/48", "2": "1/48", "3": "1/108"},
      "pair-energy-six T6 kernel")
check(pair["2"]["X4"] == {"1": "-1/2", "2": "-1/2", "3": "-1/2"},
      "pair-energy-two X4 kernel")
check(pair["4"]["X4"] == {"1": "-7/32", "2": "-7/32", "3": "-3/32"},
      "pair-energy-four X4 kernel")
check(pair["6"]["X4"] == {"1": "-5/36", "2": "-5/36", "3": "-1/27"},
      "pair-energy-six X4 kernel")
triples = ledger["history_derivative_kernels"]["triple_canonical"]
check(len(triples) == 7, "seven canonical triple kernels")
check(triples["((4, 4, 4), 6)"]["7"] == "1/64",
      "matching all-three score coefficient")
check(triples["((2, 2, 4), 2)"]["7"] == "25/64",
      "opposite path all-three score coefficient")
check(triples["((4, 6, 6), 10)"]["7"] == "49/14400",
      "equal path all-three score coefficient")

rooted = ledger["universal_rooted_census"]
check(len(rooted["locked_word_rows"]) == 6, "all six locked words")
check(len({tuple(row["locked_word"]) for row in rooted["locked_word_rows"]}) == 6,
      "locked words unique")
check(all(len(set(row["singleton_coefficients"])) == 1
          for row in rooted["locked_word_rows"]), "singleton port uniformity")
check(all(len(set(row["triple_coefficients"])) == 1
          for row in rooted["locked_word_rows"]), "triple port uniformity")
check(rooted["T6_singleton"] ==
      "(15/128)M^2+(3049/3456)M+8653/4800", "T6 singleton polynomial")
check(rooted["T6_triple"] == "49/576", "T6 triple coefficient")
check(rooted["X4_singleton"] == "-(5/16)M-487/432", "X4 singleton polynomial")
check(rooted["complete_K6_singleton"] ==
      "(1/128)M^2+(283/1152)M+8653/4800", "complete singleton polynomial")
check(rooted["complete_K6_triple"] == "49/576",
      "complete triple coefficient")
check(rooted["pure_T2_conclusion"] ==
      "zero pointwise for all six locked words", "universal rooted conclusion")
for row in ledger["literal_Q4_cases"]:
    check(row["pair_sets_with_local_support"] == 1014,
          f"all supported pairs: {row['case']}")
    check(row["triple_sets_with_local_support"] == 128020,
          f"all supported triples: {row['case']}")
    for key in ("T6_prime", "X4_prime", "A3_prime", "A2_prime",
                "K6_prime_T_contractions"):
        check(row[key] == ["0", "0", "0"], f"{key} vanishes: {row['case']}")


# Fresh exact science replay.
replay = subprocess.run(
    [sys.executable, str(HERE / "derive_complete_diagonal_h6_t2_source.py")],
    cwd=HERE, text=True, capture_output=True, timeout=180, check=False)
check(replay.returncode == 0, "science replay exits zero")
check("PASS__GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE__10775/10775"
      in replay.stdout, "science replay exact pass token")
check("DIAGONAL_H6_T2=ZERO_POINTWISE;COMPLETE_FIRST_T2_THROUGH_H6=GL6CH_WRITER"
      in replay.stdout, "science replay conclusion token")


# Textual claim guards.
theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "R'_S(0)=+{m(S)\\over E(S)^2}",
    "K'_{6,T}=T'_{6,T}-bX'_{4,T}+b^2A'_{3,T}-dA'_{2,T}",
    "all six central locked words are evaluated exactly",
    "\\sum_{a=0}^3w_a=0_6",
    "\\langle s|K'_{6,T}|s\\rangle=0",
    "every spatially nonuniform pure-`T2` first-source profile",
    "complete **first-source** effective operator through sixth order",
    "{105\\over16}{h^6\\over U_d^6}",
    "mixed two-site or nonuniform source-source contacts remain open",
    "No phase, record, metric, Ricci, gravity, graviton, or numerical-`G` claim",
):
    check(token in theorem, f"theorem token: {token}")
check("not a selected-background cancellation" in result,
      "result distinguishes universal proof from samples")
check("does not turn the tensor writer into" in readme,
      "README blocks promotion")
check("This is an author self-audit, not an independent hostile audit."
      in self_audit, "authorship/audit status explicit")
check("PASS__GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE__10775/10775"
      in verification, "verification records physics replay")
check("PASS__GL6CN_PACKET__264/264"
      in verification, "verification records packet replay")

for forbidden in (
    "proves gravity",
    "is gravity",
    "derives the Ricci tensor",
    "calculates G",
    "proves a stationary phase",
    "proves a causal cone",
    "is a graviton",
):
    check(forbidden not in theorem + " " + result + " " + readme,
          f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6CN_PACKET__{checks}/{checks}")
