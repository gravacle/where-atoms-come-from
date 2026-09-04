#!/usr/bin/env python3
"""Fail-closed custody, replay, and scope verifier for GL6CO."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
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
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EXACT_LEDGER.json",
    "VERIFICATION.txt", "derive_cycle_response_tensor_matching.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required packet file {name}")

dependency_lines = [
    line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_lines) == 26, "twenty-six exact dependency bytes pinned")
dependency_names = set()
for line in dependency_lines:
    expected, relative = line.split("  ", 1)
    check(relative not in dependency_names, f"unique dependency {relative}")
    dependency_names.add(relative)
    path = ROOT / relative
    check(path.is_file(), f"dependency exists {relative}")
    check(digest(path) == expected, f"dependency hash {relative}")

for packet in (
    "LANE_CROSS_RFT_GRA_GL6CL_GLOBAL_FOURIER_PAIR_WRITER_V001",
    "LANE_CROSS_RFT_GRA_GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001",
    "AUDIT_G_GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001",
    "LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001",
    "AUDIT_G_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001",
):
    manifest = ROOT / packet / "MANIFEST.sha256"
    seal = (ROOT / packet / "SEAL.sha256").read_text().strip().split("  ", 1)
    check(seal[0] == digest(manifest), f"upstream seal hashes manifest {packet}")
    check(seal[1] in {"MANIFEST.sha256", f"{packet}/MANIFEST.sha256"},
          f"upstream seal names manifest {packet}")

ledger = json.loads((HERE / "EXACT_LEDGER.json").read_text(),
                    object_pairs_hook=unique_object)
classification = ledger["classification"]
check(classification["group_order"] == 24, "full tetrahedral group classified")
check(classification["constant_invariant_dimension"] == "2",
      "two constant invariant matrices")
check(classification["quadratic_invariant_dimension"] == "5",
      "five quadratic invariant matrices")
symbol = ledger["general_centered_cycle_symbol"]
check(symbol["normalization"] ==
      "K_cycle=K_bare with K_bare_cc'=2 Re <0|T_c R T_c'|0>; it contains no writer factor",
      "bare cycle susceptibility normalization")
check("actual common pair-source-to-cycle-amplitude derivative here is mu B_T"
      in symbol["relation_to_GL6CM"], "writer scale applied exactly once")
check("no zero-derivative" in symbol["constant_tensor_guard"] and
      "masslessness" in symbol["constant_tensor_guard"],
      "constant susceptibility is not promoted to a massless endpoint")

pullback = ledger["writer_pullback"]
check(pullback["composed_coefficients"] == {
          "A": "-2 kappa+8 b",
          "B": "-16 kappa+8 c",
          "C": "12 kappa+8 d",
          "cubic_mismatch_B_plus_C": "-4 kappa+8(c+d)",
      }, "exact writer pullback coefficients")
check(pullback["SO3_extension_iff"] ==
      "B+C=0 iff c+d=kappa/2", "one tensor-extension relation")
check(pullback["reference_FP_additional_condition"] ==
      "A=0 iff b=kappa/4, in addition to c+d=kappa/2",
      "stronger reference-shape relation")

contact = ledger["h2_contact_separate_block"]
check(contact["SO3_mismatch_B_plus_C"] == "(8/3)(2p-1)",
      "contact mismatch")
check(contact["contact_alone_SO3_condition"] == "p=1/2",
      "contact-only extension point")
check("same stationary state" in contact["same_state_guard"],
      "contact same-state guard")

normalization = ledger["common_source_normalization"]
check(normalization == {
          "GL6CL_coordinates": "j_P=j_plus+j_minus; j_C=j_plus-j_minus, so j_plus=(j_P+j_C)/2",
          "orthonormal_common_coordinate": "jhat_plus=(j_P+j_C)/sqrt(2)=sqrt(2)j_plus",
          "normalized_writer": "delta a=(mu/sqrt(2))B_T jhat_plus",
          "normalized_cycle_hessian": "Hhat=(mu^2/2)B_T^*K_cycle B_T",
      }, "common-sublattice normalization bridge")
same_state = ledger["conditional_same_state_total"]
check(same_state["SO3_extension_condition"].startswith("(mu^2/2)"),
      "normalized combined extension equation")
check(same_state["reference_FP_additional_condition"].startswith("(mu^2/2)"),
      "normalized combined reference-shape equation")

reference = ledger["full_tensor_reference_kill_test"]
check(reference["rank"] == 3, "full reference rank-three quotient")
check(reference["T_T_block"] == [
          ["4", "-6", "-10"],
          ["-6", "9", "-15"],
          ["-10", "-15", "25"],
      ], "generic-momentum T-T reference block")
check(any(reference["bilinear_matrix"][i][j] != "0"
          for i in range(3) for j in range(3, 6)),
      "reference has nonzero A/E-to-T cross blocks")
check("-1/sqrt(2)" in reference["orthonormal_T_solder"],
      "orthonormal solder normalization explicit")

power = ledger["strong_lock_power_counting"]
check(power["cycle_pullback"] == "mu^2/J=O(h^6/U_d^7)",
      "cycle response scaling")
check(power["contact"] == "h^2/(4U_d^3)=O(h^2/U_d^3)",
      "contact scaling")
check(power["guard"] == "power counting only; no phase transition is claimed",
      "power-counting claim ceiling")
check("no phase/spacetime/metric/Ricci/Einstein/gravity/G claim" in
      ledger["ceiling"], "ledger promotion ceiling")

replay = subprocess.run(
    [sys.executable, "-B",
     str(HERE / "derive_cycle_response_tensor_matching.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "exact replay exits zero")
for token in (
    "PASS__GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING__217/217",
    "CYCLE_SYMBOL=2_CONSTANT_INVARIANTS;5_QUADRATIC_INVARIANTS;STATIONARY_A1_NULL",
    "SO3_TENSOR_EXTENSION=ONE_CONDITION_C_PLUS_D_EQ_KAPPA_OVER_2",
    "REFERENCE_FP_SHAPE=SECOND_CONDITION_B_EQ_KAPPA_OVER_4;NOT_RICCI_PROOF",
    "H2_CONTACT=SEPARATE_SAME_STATE_BLOCK;NO_CANCELLATION_OR_PHASE_ASSUMED",
    "NO_SPACETIME_METRIC_RICCI_EINSTEIN_GRAVITY_G",
):
    check(token in replay.stdout, f"exact replay token {token}")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
for token in (
    "matching is possible and has codimension one",
    "stationary common-amplitude null",
    "normalization used here is deliberately **bare**",
    "prevents double counting the writer",
    "its common coordinate is `j_+=(j_P+j_C)/2`, not the orthonormal parent/child coordinate",
    "{\\mu^2\\over2}[-4\\kappa+8(c+d)]",
    "c+d={\\kappa\\over2}",
    "These coordinates are **not** an ordinary three-vector",
    "rank three",
    "nonzero `A1/E2-T2` cross blocks",
    "same stationary state",
    "cannot generically cancel a nonzero leading contact mismatch",
    "power counting, not a phase-transition claim",
    "does not prove the zero-derivative/background-stationarity",
    "No such promotion is made here",
):
    check(token in theorem, f"theorem claim/scope token {token}")
check("not a Ricci proof" in result and
      "No cancellation or phase is asserted" in result,
      "result preserves Ricci and phase ceilings")
check("executable kill test" in self_audit and
      "does not replace an independent hostile audit" in self_audit,
      "self-audit states kill test and independence ceiling")

manifest_lines = [
    line for line in (HERE / "MANIFEST.sha256").read_text().splitlines()
    if line.strip()
]
manifest_names = set()
for line in manifest_lines:
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"manifest local path {name}")
    check(name not in manifest_names, f"manifest unique path {name}")
    manifest_names.add(name)
    check(digest(HERE / name) == expected, f"manifest hash {name}")
check(manifest_names == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "manifest covers every non-custody packet byte")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "seal names and hashes packet manifest")

print(f"PASS__GL6CO_PACKET__{checks}/{checks}")
