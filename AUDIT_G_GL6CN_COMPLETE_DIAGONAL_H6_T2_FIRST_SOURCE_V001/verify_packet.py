#!/usr/bin/env python3
"""Fail-closed custody and independent replay verifier for the GL6CN audit."""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET_NAME = "LANE_CROSS_RFT_GRA_GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE_V001"
TARGET = ROOT / TARGET_NAME
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


required = {
    "README.md", "AUDIT_REPORT.md", "INDEPENDENT_RESULT.json",
    "TARGET.sha256", "VERIFICATION.txt", "verify_gl6cn_independent.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required audit file {name}")

target_payload = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EXACT_LEDGER.json",
    "VERIFICATION.txt", "derive_complete_diagonal_h6_t2_source.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
target_lines = [line for line in (HERE / "TARGET.sha256").read_text().splitlines()
                if line.strip()]
check(len(target_lines) == len(target_payload), "twelve frozen target pins")
pinned = set()
for line in target_lines:
    expected, relative = line.split("  ", 1)
    path = Path(relative)
    check(not path.is_absolute() and ".." not in path.parts,
          f"safe target path {relative}")
    check(len(path.parts) == 2 and path.parts[0] == TARGET_NAME,
          f"target-scoped pin {relative}")
    check(path.name not in pinned, f"unique target pin {path.name}")
    pinned.add(path.name)
    check((ROOT / path).is_file(), f"target byte exists {relative}")
    check(digest(ROOT / path) == expected, f"frozen target hash {relative}")
check(pinned == target_payload, "target pins cover payload, manifest, and seal")
target_seal = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
check(target_seal == [digest(TARGET / "MANIFEST.sha256"),
                      f"{TARGET_NAME}/MANIFEST.sha256"],
      "target seal closes frozen target manifest")
check(digest(TARGET / "MANIFEST.sha256") ==
      "79cd181ce77f5f663f2847c4b816fb8cbe32344ef8d2b9bb222ff902a9966c10",
      "audited target manifest identity")

# The independent result is checked before and after replay.  The replay is
# the separately authored implementation only; the target author program is
# deliberately never invoked here.
result = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text(),
                    object_pairs_hook=unique_object)
check(result["verdict"] == "PASS", "independent verdict")
check(result["checks"] == 438506, "independent science check count")
check(result["independence"] == "author derivation was neither imported nor executed",
      "independent implementation custody")
check((result["pair_kernel_classes"], result["triple_canonical_classes"],
       result["triple_labelled_signatures"]) == (3, 7, 22),
      "complete independent history classification")
check(result["dual_number_rule"] ==
      "d[-1/(E+eta m)]/deta=+m/E^2", "resolvent derivative sign")
check(result["fold_signs"] ==
      "K6'=T6'-bX4'+b^2A3'-dA2'; A2'=-m/4; A3'=+3m/16",
      "Kato and fold signs")
check("simple degree-four bipartite girth-at-least-six" in result["rooted_domain"],
      "exact audited graph domain")
check("sole pure-T2 first-source operator through h6" in
      result["integrated_corollary"], "integrated conditional corollary")
check("gravity, and G remain open" in result["ceiling"],
      "promotion ceiling retained")

witnesses = result["independent_parent_witnesses"]
check(len(witnesses) == 24, "twenty-four independent parent/word witnesses")
expected_occupancies = set(combinations(range(4), 2))
seen = {}
for witness in witnesses:
    n = witness["n"]
    m = witness["M"]
    occupied = tuple(witness["occupied_ports"])
    check(n in {17, 19, 23, 29}, f"declared witness size {n}")
    check(m == 4 * n, f"degree-four edge count n={n}")
    check(occupied in expected_occupancies, f"locked port word n={n} {occupied}")
    check((n, occupied) not in seen, f"unique parent/word witness n={n} {occupied}")
    seen[(n, occupied)] = True
    check(witness["pairs_touching_source"] ==
          m * (m - 1) // 2 - (m - 4) * (m - 5) // 2,
          f"complete pair support n={n} {occupied}")
    check(witness["triples_touching_source"] ==
          m * (m - 1) * (m - 2) // 6 -
          (m - 4) * (m - 5) * (m - 6) // 6,
          f"complete triple support n={n} {occupied}")
    check(F(witness["T6_singleton"]) ==
          F(15, 128) * m * m + F(3049, 3456) * m + F(8653, 4800),
          f"T6 singleton polynomial n={n} {occupied}")
    check(F(witness["T6_complement_triple"]) == F(49, 576),
          f"T6 triple coefficient n={n} {occupied}")
    check(F(witness["X4_singleton"]) ==
          -F(5, 16) * m - F(487, 432),
          f"X4 singleton polynomial n={n} {occupied}")
    check(F(witness["K6_singleton"]) ==
          F(1, 128) * m * m + F(283, 1152) * m + F(8653, 4800),
          f"K6 singleton polynomial n={n} {occupied}")
for n in (17, 19, 23, 29):
    check({occupied for size, occupied in seen if size == n} == expected_occupancies,
          f"all six locked words represented at n={n}")

target_ledger = json.loads((TARGET / "EXACT_LEDGER.json").read_text(),
                           object_pairs_hook=unique_object)
check(target_ledger["complete_diagonal_h6_T2_vertex"] == "zero pointwise",
      "target pointwise conclusion")
check(target_ledger["all_complete_t_contractions_zero"] is True,
      "target complete contractions zero")
check(target_ledger["all_individual_kato_terms_t_contractions_zero"] is True,
      "target individual Kato terms zero")
check(target_ledger["dual_number_crosscheck"]["triple_signatures_checked"] == 22,
      "target checks all labelled triple signatures")
check(target_ledger["dual_number_crosscheck"]["A2_prime_single_score"] == "-1/4" and
      target_ledger["dual_number_crosscheck"]["A3_prime_single_score"] == "3/16",
      "target fold derivative signs agree")
rooted = target_ledger["universal_rooted_census"]
check(rooted["T6_singleton"] ==
      "(15/128)M^2+(3049/3456)M+8653/4800", "target T6 census polynomial")
check(rooted["T6_triple"] == "49/576", "target triple coefficient")
check(rooted["X4_singleton"] == "-(5/16)M-487/432",
      "target X4 census polynomial")
check(rooted["complete_K6_singleton"] ==
      "(1/128)M^2+(283/1152)M+8653/4800", "target K6 census polynomial")
check(rooted["pure_T2_conclusion"] == "zero pointwise for all six locked words",
      "target universal local conclusion")
check(target_ledger["integrated_corollary"] ==
      "through h6, the complete pure-T2 first-source effective operator is the GL6CH off-diagonal six-cycle writer",
      "target integrated writer-only conclusion")

report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
for token in (
    "**PASS**, within the exact domain and ceilings stated by the target",
    "neither imported nor executed the target's author derivation",
    "22 distinct labelled signatures",
    "24 independent parent/word witnesses",
    "Girth at least six is load-bearing",
    "port-uniform separately for singleton and complementary triple masks",
    "without cancellation between differently weighted ports",
    "first-source only",
    "an `h8` diagonal first-source term is not excluded",
    "No material defect was found",
):
    check(token in report, f"audit report claim/scope token {token}")
check("without importing or executing the author derivation" in readme,
      "README independence statement")
check("gravity, and `G` remain open" in readme, "README promotion ceiling")

replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_gl6cn_independent.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "independent replay exits zero")
for token in (
    "PASS__GL6CN_INDEPENDENT_HOSTILE_SCIENCE__438506/438506",
    "HISTORY_KERNELS=3_PAIR;7_CANONICAL_TRIPLE;22_LABELLED_TRIPLE",
    "DUAL_SIGNS=RESOLVENT_PLUS;A2_MINUS_1_4;A3_PLUS_3_16;KATO_SIGNS_PASS",
    "ROOTED_CENSUS=EXHAUSTIVE_DEGREE4_GIRTH6;24_INDEPENDENT_PARENT_WORD_WITNESSES",
    "POINTWISE_T2=PORT_UNIFORM_SINGLETON_PLUS_COMPLEMENT_TRIPLE_CANCEL",
    "INTEGRATED_THROUGH_H6=OFFDIAGONAL_GL6CH_WRITER_ONLY_CONDITIONAL_ON_PINNED_INPUTS",
    "SOURCE_SECOND_H8_PHASE_RECORD_BULK_RICCI_GRAVITY_G_OPEN",
):
    check(token in replay.stdout, f"independent replay token {token}")
replayed = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text(),
                      object_pairs_hook=unique_object)
check(replayed == result, "independent result is byte-semantic deterministic")

manifest_lines = [line for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
                  if line.strip()]
manifest_names = set()
for line in manifest_lines:
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"manifest local path {name}")
    check(name not in manifest_names, f"manifest unique path {name}")
    manifest_names.add(name)
    check(digest(HERE / name) == expected, f"manifest hash {name}")
check(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest covers every non-custody audit byte")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "audit seal names and hashes manifest")

print(f"PASS__GL6CN_INDEPENDENT_AUDIT_PACKET__{checks}/{checks}")
