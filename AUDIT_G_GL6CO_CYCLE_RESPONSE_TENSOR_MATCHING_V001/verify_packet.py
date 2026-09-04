#!/usr/bin/env python3
"""Fail-closed verifier for the independent GL6CO hostile audit."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING_V001"
checks = 0


def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def strict_json(path):
    def no_duplicates(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    return json.loads(path.read_text(), object_pairs_hook=no_duplicates)


required = {
    "README.md",
    "AUDIT_REPORT.md",
    "TARGET.sha256",
    "INDEPENDENT_RESULT.json",
    "VERIFICATION.txt",
    "verify_gl6co_independent.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
}
for name in required:
    check((HERE / name).is_file(), f"required file: {name}")


# Pin all twelve repaired author bytes.
target_paths = set()
for line in (HERE / "TARGET.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    path = ROOT / relative
    check(path.is_file(), f"target exists: {relative}")
    check(digest(path) == expected, f"target hash: {relative}")
    check(relative not in target_paths, f"target row unique: {relative}")
    target_paths.add(relative)
check(len(target_paths) == 12, "exact target byte count")
check(digest(TARGET / "MANIFEST.sha256") ==
      "f085a73c4d7590a44ac89117f53bdb00583646153396a16552a84829ebd323b6",
      "repaired target manifest digest")
check(digest(TARGET / "SEAL.sha256") ==
      "1f9847a10993954ba6baca4b5b127f0f9c215c2a172faa0f0a0e533bb3ea7657",
      "repaired target seal-file digest")

target_manifest_names = set()
for line in (TARGET / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, name = line.split(maxsplit=1)
    check((TARGET / name).is_file(), f"target manifest file exists: {name}")
    check(digest(TARGET / name) == expected, f"target manifest hash: {name}")
    check(name not in target_manifest_names, f"target manifest unique: {name}")
    target_manifest_names.add(name)
check(len(target_manifest_names) == 10, "exact target manifest count")
target_seal = [row for row in (TARGET / "SEAL.sha256").read_text().splitlines()
               if row.strip()]
check(len(target_seal) == 1, "one target seal row")
expected, name = target_seal[0].split(maxsplit=1)
check(name == "MANIFEST.sha256", "target seal names local manifest")
check(digest(TARGET / "MANIFEST.sha256") == expected, "target seal closes")

target_verification = (TARGET / "VERIFICATION.txt").read_text()
check("PASS__GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING__217/217"
      in target_verification, "target science pass recorded")
check("PASS__GL6CO_PACKET__180/180" in target_verification,
      "target packet pass recorded")

# Recheck every dependency byte frozen by the target.
dependency_paths = set()
for line in (TARGET / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    path = ROOT / relative
    check(path.is_file(), f"target dependency exists: {relative}")
    check(digest(path) == expected, f"target dependency hash: {relative}")
    check(relative not in dependency_paths, f"target dependency unique: {relative}")
    dependency_paths.add(relative)
check(len(dependency_paths) == 26, "exact target dependency count")


# Independent exact result.
result = strict_json(HERE / "INDEPENDENT_RESULT.json")
check(result["audit"] == "GL6CO independent hostile reconstruction",
      "independent result identity")
check(result["disposition"] == "PASS", "independent mathematical disposition")
check(result["checks"] == 277, "independent mathematical check count")
group = result["group_and_invariants"]
check(group["group_order"] == 24, "full S4 enumerated")
check(group["quadratic_invariant_dimension"] == "5",
      "quadratic invariant dimension")
check(group["constant_invariant_dimension"] == "2",
      "constant invariant dimension")
check(group["basis_labels"] == ["alpha", "eta", "b", "c", "d"],
      "five independent basis labels")

writer = result["writer"]
check(writer["offset_count"] == 12, "twelve centered offsets")
check(writer["offset_norm_squared"] == "11/4", "offset norm")
check(writer["offset_second_moment"] == "11 I", "offset moment")
check(writer["pullback"] == {
          "constant": "8 kappa I",
          "A": "-2 kappa+8 b",
          "B": "-16 kappa+8 c",
          "C": "12 kappa+8 d",
          "alpha_eta_through_k2": "zero",
      }, "independent writer pullback")
check(writer["writer_digest"] ==
      "758c1da77563edf1772ef80443401585a81859b563bbcf13d6cbc9c42eb85ad4",
      "independent writer-jet digest")

tensor = result["tensor_extension_and_reference"]
check(tensor["SO3_extension_iff"] == "B+C=0", "SO3 extension plane")
check(tensor["cycle_condition"] == "c+d=kappa/2", "cycle matching relation")
check(tensor["reference_additional_condition"] == "A=0 iff b=kappa/4",
      "stronger reference condition")
check(tensor["reference_TT"] == "D-O", "reference T-T shape")
check(tensor["reference_rank"] == 3, "full reference rank")
check(any(tensor["reference_matrix"][i][j] != 0
          for i in range(3) for j in range(3, 6)),
      "full reference cross blocks nonzero")

contact = result["contact"]
check(contact["A"] == "(4/3)(1-4p)", "contact A coefficient")
check(contact["B"] == "0", "contact B coefficient")
check(contact["C_and_mismatch"] == "(8/3)(2p-1)",
      "contact C/mismatch")
check(contact["Q4_witness"] == "p=109/128 gives mismatch 15/8",
      "contact witness")
check(contact["energy_hessian_sign"].startswith("-"), "energy Hessian sign")
check(contact["connected_functional_sign"].startswith("+"),
      "connected-functional sign")

normalization = result["common_sublattice_normalization"]
check(normalization["normalized_writer"] == "B_+/sqrt(2)",
      "normalized common writer")
check(normalization["normalized_cycle_block"] ==
      "(mu^2/2) B_+^* K_cycle B_+", "normalized cycle Hessian")
check(normalization["normalized_contact_block"] == "g_ct(8+Re C)",
      "normalized contact Hessian")
check(normalization["normalized_extension_equation"].startswith("(mu^2/2)"),
      "normalized combined extension equation")
check("multiply every displayed contact coefficient by 2" in
      normalization["equivalent_unnormalized_rule"],
      "equivalent unnormalized convention")
check("one stationary state" in
      result["conditional_total"]["same_state_guard"], "same-state guard")
check("1PI Ricci kernel" in
      result["conditional_total"]["one_particle_irreducible_guard"],
      "response-to-1PI guard")

replay = subprocess.run(
    [sys.executable, str(HERE / "verify_gl6co_independent.py")],
    cwd=HERE, text=True, capture_output=True, timeout=60, check=False)
check(replay.returncode == 0, "independent replay exits zero")
for token in (
    "PASS__GL6CO_INDEPENDENT_HOSTILE_REPLAY__278/278",
    "S4_QUADRATIC_INVARIANT_DIMENSION=5;CONSTANT_DIMENSION=2;COMMON_NULL_APPLIED",
    "PULLBACK=A_-2K+8B__B_-16K+8C__C_12K+8D;ALPHA_ETA_ABSENT_AT_K2",
    "SO3_EXTENSION_IFF=B+C_ZERO_IFF_C+D=KAPPA/2;REFERENCE_ADDS_B=KAPPA/4",
    "COMMON_NORMALIZATION=NORMALIZED_CYCLE_MU2_OVER2;CONTACT_G_CT;EQUIVALENT_UNNORMALIZED_CONTACT_TIMES2",
    "DISPOSITION=PASS;MATCHING_ONLY;NO_1PI_RICCI_GRAVITY_OR_G",
):
    check(token in replay.stdout, f"independent replay token: {token}")


# Author-repair and claim guards.
target_theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
for token in (
    "These coordinates are **not** an ordinary three-vector",
    "j_+=(j_P+j_C)/2`, not the orthonormal parent/child coordinate",
    "{\\mu^2\\over2}[-4\\kappa+8(c+d)]",
    "{\\mu^2\\over2}[-2\\kappa+8b]",
    "response-to-1PI operation",
    "same stationary state",
    "This is power counting, not a phase-transition claim",
):
    check(token in target_theorem, f"repaired target token: {token}")
check('"These coordinates are **not** an ordinary three-vector"' in
      (TARGET / "verify_packet.py").read_text(),
      "target verifier enforces explicit tensor negation")

audit = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "**PASS_AFTER_AUTHOR_REPAIR**",
    "missing a relative factor of two",
    "The author repaired the first guard",
    "{\\cal H}^{H6}_{+,T}={\\mu^2\\over2}B_+^*K_{\\rm cyc}B_+",
    "These coordinates are not an ordinary three-vector",
    "passing (A-CO09) in a three-dimensional subblock is not a Ricci proof",
    "No match is asserted",
    "contemporaneous GL6CN result",
    "No source-second contact is double-counted",
):
    check(token in audit, f"audit report token: {token}")
check("Disposition: `PASS_AFTER_AUTHOR_REPAIR`" in readme,
      "README disposition")
check("PASS__GL6CO_INDEPENDENT_HOSTILE_REPLAY__278/278" in verification,
      "verification records independent replay")
check("PASS__AUDIT_G_GL6CO_PACKET__258/258"
      in verification, "verification records audit packet replay")

for forbidden in (
    "proves gravity",
    "is gravity",
    "derives the Ricci tensor",
    "calculates G",
    "is a graviton",
):
    check(forbidden not in audit + " " + readme,
          f"forbidden promotion absent: {forbidden}")


# Audit manifest and seal.
manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    path = ROOT / relative
    check(path.is_file(), f"audit manifest target exists: {relative}")
    check(digest(path) == expected, f"audit manifest hash: {relative}")
    check(relative not in manifest_paths, f"audit manifest unique: {relative}")
    manifest_paths.add(relative)
for name in required - {"MANIFEST.sha256", "SEAL.sha256"}:
    check(f"{HERE.name}/{name}" in manifest_paths,
          f"audit manifest coverage: {name}")
check(len(manifest_paths) == 7, "exact audit manifest count")

seal_lines = [row for row in (HERE / "SEAL.sha256").read_text().splitlines()
              if row.strip()]
check(len(seal_lines) == 1, "one audit seal row")
expected, relative = seal_lines[0].split(maxsplit=1)
check(relative == f"{HERE.name}/MANIFEST.sha256", "audit seal target")
check(digest(HERE / "MANIFEST.sha256") == expected, "audit seal closes")

print(f"PASS__AUDIT_G_GL6CO_PACKET__{checks}/{checks}")
