#!/usr/bin/env python3
"""Fail-closed custody and independent replay verifier for the GL6CQ audit."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET_NAME = "LANE_CROSS_RFT_GRA_GL6CQ_STATIONARY_RESPONSE_MOMENT_SUM_RULES_V001"
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
    "TARGET.sha256", "VERIFICATION.txt", "verify_gl6cq_independent.py",
    "verify_packet.py", "MANIFEST.sha256", "SEAL.sha256",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required audit file {name}")

target_payload = {
    "README.md", "THEOREM.md", "RESULT.md", "SELF_AUDIT.md",
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "EXACT_LEDGER.json",
    "VERIFICATION.txt", "derive_stationary_response_moment_sum_rules.py",
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
check(target_seal == [digest(TARGET / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "target seal closes target manifest")
check(digest(TARGET / "MANIFEST.sha256") ==
      "bc0c979698a0738db9c841bfabd718023940fe3a43d397f04332bbe5009fe469",
      "audited target manifest identity")

# Verify every exact upstream byte and close every pinned upstream seal.
dependency_lines = [
    line for line in (TARGET / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_lines) == 27, "twenty-seven upstream dependency pins")
dependency_names = set()
dependency_paths = []
for line in dependency_lines:
    expected, relative = line.split("  ", 1)
    path = Path(relative)
    check(not path.is_absolute() and ".." not in path.parts,
          f"safe dependency path {relative}")
    check(relative not in dependency_names, f"unique dependency pin {relative}")
    dependency_names.add(relative)
    dependency_paths.append(path)
    check((ROOT / path).is_file(), f"dependency exists {relative}")
    check(digest(ROOT / path) == expected, f"dependency hash {relative}")
for path in dependency_paths:
    if path.name != "SEAL.sha256":
        continue
    packet = path.parent
    seal = (ROOT / path).read_text().strip().split("  ", 1)
    check(seal[0] == digest(ROOT / packet / "MANIFEST.sha256"),
          f"upstream seal hashes manifest {packet}")
    check(seal[1] in {"MANIFEST.sha256", f"{packet}/MANIFEST.sha256"},
          f"upstream seal names manifest {packet}")
check(any(str(path) ==
          "LANE_CROSS_RFT_GRA_GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING_V001/MANIFEST.sha256"
          for path in dependency_paths), "repaired GL6CO manifest is pinned")
check(digest(ROOT / "LANE_CROSS_RFT_GRA_GL6CO_CYCLE_RESPONSE_TENSOR_MATCHING_V001" /
             "MANIFEST.sha256") ==
      "f085a73c4d7590a44ac89117f53bdb00583646153396a16552a84829ebd323b6",
      "repaired GL6CO manifest identity")

ledger = json.loads((TARGET / "EXACT_LEDGER.json").read_text(),
                    object_pairs_hook=unique_object)
check(ledger["schema"] == "GL6CQ_STATIONARY_RESPONSE_MOMENT_SUM_RULES_V001",
      "target ledger schema")
check(ledger["moments"]["projected_coefficients"] == {
          "kappa": "Z_T/3",
          "alpha": "-M_AA_trace/6",
          "eta": "-(M_AT,x^yz+M_AT,y^zx+M_AT,z^xy)/3",
          "b": "-M_perp/12",
          "c": "-M_parallel/6+M_perp/12",
          "d": "-M_cross/6",
      }, "target six moment projections")
check(ledger["moments"]["contractions"] == {
          "Z_T": "sum_i Q_i^T Z Q_i=3 kappa",
          "M_perp": "sum_{i!=m} M_TT,ii^{mm}=-12 b",
          "M_parallel": "sum_i M_TT,ii^{ii}=-6(b+c)",
          "M_cross": "sum_{i!=j} M_TT,ij^{ij}=-6d",
      }, "target moment contraction multiplicities")
normalization = ledger["normalization"]
check(normalization["cycle_kernel"] == "bare: no lambda_T^2 and no mu^2",
      "target bare kernel has no writer factor")
check(normalization["GL6CL_common_coordinate"] == "j_plus=(j_P+j_C)/2",
      "target unnormalized GL6CL common coordinate")
check(normalization["orthonormal_common_coordinate"] ==
      "jhat_plus=(j_P+j_C)/sqrt(2)=sqrt(2)j_plus",
      "target orthonormal common coordinate")
check(normalization["GL6CL_common_pullback"] ==
      "H_T^H6=mu^2 B_T^* K_bare B_T", "target unnormalized pullback")
check(normalization["GL6BV_normalized_common_pullback"] ==
      "Hhat_T^H6=(mu^2/2) B_T^* K_bare B_T",
      "target normalized pullback factor one half")
check(normalization["contact_scale"] ==
      "g_ct=h^2/(4U_d^3), with no mu factor", "target contact has no mu")
rules = ledger["observable_sum_rules"]
check(rules["tensor_gradient_extension"] ==
      "(mu^2/2)[-2 Z_T+M_perp-2Mparallel-2Mcross]+4g_ct(2p-1)=0".replace(
          "Mparallel", "M_parallel").replace("Mcross", "M_cross"),
      "target normalized extension observable rule")
check(rules["reference_shape"] ==
      "-(mu^2/2)[Z_T+M_perp]+2g_ct(1-4p)=0",
      "target normalized reference observable rule")
check("same stationary state" in rules["guard"] and
      "1/2 converts" in rules["guard"], "target same-state normalization guard")
check("zero-momentum" in rules["zero_momentum_guard"] and
      "masslessness" in rules["zero_momentum_guard"],
      "target zero-mode promotion guard")

result = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text(),
                    object_pairs_hook=unique_object)
check(result["verdict"] == "PASS", "independent verdict")
check(result["independence"] ==
      "target author derivation was neither imported nor executed",
      "independent science custody")
check(result["algebra_checks_before_frozen_result_custody"] == 527,
      "independent algebra check count")
check(result["coefficient_recovery"] == {
          "kappa": "Z_T/3",
          "alpha": "-sum_m M_AA^{mm}/6",
          "eta": "-(M_AT,x^{yz}+M_AT,y^{zx}+M_AT,z^{xy})/3",
          "b": "-M_perp/12",
          "c": "-M_parallel/6+M_perp/12",
          "d": "-M_cross/6",
          "exact_sample_recoveries": [
              ["7", "11", "13", "17", "19", "23"],
              ["2/3", "-5/7", "11/13", "-17/19", "23/29", "-31/37"],
              ["0", "1", "0", "0", "0", "0"],
              ["5", "0", "0", "0", "0", "0"],
          ],
      }, "independent coefficient reconstruction")
check(result["unrescaled_observable_forms"] == {
          "extension": ["-4/3", "2/3", "-4/3", "-4/3"],
          "reference": ["-2/3", "-2/3", "0", "0"],
      }, "independent unrescaled observable forms")
check(result["observable_rules"] == {
          "extension": "(mu^2/2)[-2Z_T+M_perp-2Mparallel-2Mcross]+4g_ct(2p-1)=0".replace(
              "Mparallel", "M_parallel").replace("Mcross", "M_cross"),
          "reference": "-(mu^2/2)[Z_T+M_perp]+2g_ct(1-4p)=0",
      }, "independent boxed observable rules")
check(result["normalization"]["normalized_cycle_hessian"] ==
      "(mu^2/2)B_T^*K_bare B_T", "independent normalized Hessian")
check(result["normalization"]["contact"] ==
      "g_ct multiplies contact with no mu factor", "independent contact scale")
check("divergence requires a nonanalytic test" in result["analytic_scope"],
      "independent critical analytic boundary")
check("no satisfaction" in result["ceiling"] and "gravity" in result["ceiling"],
      "independent promotion ceiling")

target_theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
target_result = " ".join((TARGET / "RESULT.md").read_text().split())
for token in (
    "common coordinate `j_+=(j_P+j_C)/2`",
    "The half is a coordinate normalization, not another writer factor",
    "No `mu` multiplies this contact",
    "necessary-and-sufficient `T2` **quadratic-gradient**",
    "do not assert that nature or a selected phase",
    "do not prove background stationarity, masslessness, or a gauge null",
    "nonanalytic leading kernel must be measured and matched directly",
):
    check(token in target_theorem, f"target theorem claim/scope token {token}")
check("`(mu^2/2) B_T^* K^bare B_T`" in target_result,
      "target result uses normalized common half")
check("do not assert that the left sides vanish" in target_result,
      "target result preserves matching-test ceiling")

audit_report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
audit_readme = " ".join((HERE / "README.md").read_text().split())
for token in (
    "**PASS**, within the exact analytic, same-state, and sector ceilings",
    "No stale common-source normalization remains",
    "All 144 projected second-moment components",
    "there is no missing factor of two in `eta` or `d`",
    "neither a stale factor-two error nor a hidden second writer factor remains",
    "Both equations are exact if-and-only-if tests",
    "No material defect was found",
):
    check(token in audit_report, f"audit report claim token {token}")
check("neither imported nor executed" in audit_readme,
      "audit README independence statement")
check("gravity, or `G`" in audit_readme, "audit README promotion ceiling")

replay = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_gl6cq_independent.py")],
    cwd=ROOT, capture_output=True, text=True, check=False,
)
check(replay.returncode == 0, "independent replay exits zero")
for token in (
    "PASS__GL6CQ_INDEPENDENT_HOSTILE_SCIENCE__529/529",
    "MOMENT_PROJECTION=RAW_ORIENTATION_TO_A1_T2_ALL_COMPONENTS",
    "COEFFICIENTS=KAPPA_ALPHA_ETA_B_C_D_EXACT",
    "NORMALIZATION=MU_EQ_2LAMBDA;COMMON_GRAM_2;NORMALIZED_HESSIAN_MU2_OVER_2",
    "OBSERVABLE_RULES=EXTENSION_AND_REFERENCE_EXACTLY_REDERIVED",
    "NO_STALE_MU2;CONTACT_HAS_NO_MU",
    "NO_SATISFACTION_PHASE_1PI_RICCI_GRAVITY_G",
):
    check(token in replay.stdout, f"independent replay token {token}")
replayed = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text(),
                      object_pairs_hook=unique_object)
check(replayed == result, "independent result is deterministic")

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

print(f"PASS__GL6CQ_INDEPENDENT_AUDIT_PACKET__{checks}/{checks}")
