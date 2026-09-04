#!/usr/bin/env python3
"""Fail-closed custody, replay, and scope verifier for GL6CQ."""

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
    "VERIFICATION.txt", "derive_stationary_response_moment_sum_rules.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required packet file {name}")

dependency_lines = [
    line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_lines) == 27, "twenty-seven exact dependency bytes pinned")
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
    "AUDIT_G_GL6CL_GLOBAL_FOURIER_PAIR_WRITER_V001",
    "LANE_CROSS_RFT_GRA_GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001",
    "AUDIT_G_GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001",
    "LANE_CROSS_RFT_GRA_GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING_V001",
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
check(ledger["schema"] == "GL6CQ_STATIONARY_RESPONSE_MOMENT_SUM_RULES_V001",
      "ledger schema")
moments = ledger["moments"]
check(moments["expansion"] ==
      "K(k)=Z-(1/2)k_m k_n M^{mn}+o(k^2)", "moment expansion sign")
check(moments["projected_coefficients"] == {
          "kappa": "Z_T/3",
          "alpha": "-M_AA_trace/6",
          "eta": "-(M_AT,x^yz+M_AT,y^zx+M_AT,z^xy)/3",
          "b": "-M_perp/12",
          "c": "-M_parallel/6+M_perp/12",
          "d": "-M_cross/6",
      }, "all six projected coefficient formulas")
check(moments["contractions"] == {
          "Z_T": "sum_i Q_i^T Z Q_i=3 kappa",
          "M_perp": "sum_{i!=m} M_TT,ii^{mm}=-12 b",
          "M_parallel": "sum_i M_TT,ii^{ii}=-6(b+c)",
          "M_cross": "sum_{i!=j} M_TT,ij^{ij}=-6d",
      }, "contracted moment normalizations")
normalization = ledger["normalization"]
check(normalization["cycle_kernel"] ==
      "bare: no lambda_T^2 and no mu^2", "bare kernel has no writer scale")
check(normalization["CL_mu"] ==
      "(105/8)h^6/U_d^6=2 lambda_T", "mu equals twice lambda_T")
check(normalization["GL6CL_common_pullback"] ==
      "H_T^H6=mu^2 B_T^* K_bare B_T", "mu squared once in GL6CL coordinate")
check(normalization["GL6BV_normalized_common_pullback"] ==
      "Hhat_T^H6=(mu^2/2) B_T^* K_bare B_T",
      "orthonormal GL6BV common-source factor one half")
check(normalization["contact_scale"] ==
      "g_ct=h^2/(4U_d^3), with no mu factor", "contact has no mu")

rules = ledger["observable_sum_rules"]
check(rules["tensor_gradient_extension"] ==
      "(mu^2/2)[-2 Z_T+M_perp-2M_parallel-2M_cross]+4g_ct(2p-1)=0",
      "CO29 moment sum rule")
check(rules["reference_shape"] ==
      "-(mu^2/2)[Z_T+M_perp]+2g_ct(1-4p)=0", "CO30 moment sum rule")
check("same stationary state" in rules["guard"], "same-state guard")
check("zero-momentum" in rules["zero_momentum_guard"] and
      "masslessness" in rules["zero_momentum_guard"],
      "zero-momentum and masslessness guard")
check("nonanalytic" in ledger["convergence"]["critical_limit"],
      "critical nonanalyticity guard")
check("no phase" in ledger["scope"] and "gravity" in ledger["scope"],
      "promotion ceiling")

replay = subprocess.run(
    [sys.executable, "-B",
     str(HERE / "derive_stationary_response_moment_sum_rules.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "exact replay exits zero")
for token in (
    "PASS__GL6CQ_STATIONARY_RESPONSE_MOMENT_SUM_RULES__28/28",
    "MOMENTS=ZEROTH_KAPPA;SECOND_ALPHA_ETA_B_C_D",
    "NORMALIZATION=K_BARE_NO_WRITER;MU_SQUARED_ONCE;ORTHONORMAL_COMMON_HALF;CONTACT_NO_MU",
    "SAME_STATE=K_BARE_MOMENTS_AND_P",
    "SUM_RULES=CO29_CO30_AS_EXACT_REAL_SPACE_MOMENT_IDENTITIES",
    "CRITICAL_LIMIT=SECOND_MOMENT_DIVERGENCE_INVALIDATES_ANALYTIC_K2_FORM",
    "NO_PHASE_1PI_RICCI_GRAVITY_G",
):
    check(token in replay.stdout, f"exact replay token {token}")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
for token in (
    "bare connected",
    "translation-covariant thermodynamic limit",
    "entrywise absolute second-moment convergence",
    "quadratic-gradient",
    "not two factors to multiply",
    "same state `|0>`",
    "double count it",
    "zero-momentum tensor term",
    "do not assert that nature",
    "nonanalytic leading kernel",
):
    check(token in theorem, f"theorem claim/scope token {token}")
check("exact if-and-only-if tests" in result and
      "do not assert that the left sides vanish" in result,
      "result preserves matching ceiling")
check("does not replace an independent hostile audit" in self_audit,
      "self-audit independence ceiling")

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

print(f"PASS__GL6CQ_PACKET__{checks}/{checks}")
