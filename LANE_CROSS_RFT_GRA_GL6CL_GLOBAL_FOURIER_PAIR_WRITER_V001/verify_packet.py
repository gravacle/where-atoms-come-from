#!/usr/bin/env python3
"""Fail-closed custody, exact-ledger, and scope verifier for GL6CL."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise AssertionError(label)


required = (
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "EXACT_LEDGER.json",
    "VERIFICATION.txt",
    "derive_global_fourier_pair_writer.py",
    "verify_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
)
for name in required:
    check((HERE / name).is_file(), f"required file: {name}")

dependency_paths = set()
for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"dependency hash: {relative}")
    dependency_paths.add(relative)
check(len(dependency_paths) == 12, "exact dependency count")
check(all(("GL6CH" in path) for path in dependency_paths),
      "only GL6CH target and audit custody imported")
check(any(path.endswith("GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/THEOREM.md")
          for path in dependency_paths), "GL6CH theorem pinned")
check(any(path.endswith("GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/AUDIT_REPORT.md")
          for path in dependency_paths), "GL6CH hostile audit pinned")
check(all("GL6CJ" not in path for path in dependency_paths),
      "unaudited same-author GL6CJ not imported")

manifest_paths = set()
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"manifest target exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"manifest hash: {relative}")
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

ledger = json.loads((HERE / "EXACT_LEDGER.json").read_text())
geometry = ledger["geometry"]
symbol = ledger["exact_symbol"]
common = ledger["common_sector"]
relative_sector = ledger["relative_sector"]
combined = ledger["locked_read_and_common_writer"]
scalar = ledger["uniform_A1"]
check(ledger["lane"] == "GL6CL", "ledger lane")
check(len(geometry["orientations"]) == 4, "four ring orientations")
check(geometry["rho_second_moment"] ==
      [["11", "0", "0"], ["0", "11", "0"], ["0", "0", "11"]],
      "isotropic second offset moment")
check(geometry["rho_fourth_diagonal"] == ["83/4"] * 3,
      "fourth diagonal moment")
check(geometry["rho_fourth_mixed"] == ["19/4"] * 3,
      "fourth mixed moment")
check(symbol["canonical_direct_B_plus"] == "2 sum_p cos(k.rho_dp)e_p",
      "canonical-direct common cosine symbol")
check(symbol["canonical_direct_B_minus"] == "-2i sum_p sin(k.rho_dp)e_p",
      "canonical-direct relative sine symbol")
check(symbol["complete_tensor_writer"] ==
      "delta a_d^T=mu[B_plus P_T j_plus+B_minus P_T j_minus]",
      "complete tensor writer projection")
check(symbol["physical_scale"] == "mu=(105/8)h^6/U_d^6",
      "physical normalization")
check("arbitrary-profile A1/E" in symbol["scope_guard"],
      "arbitrary-profile A1/E completion guarded")
check(common["canonical_direct_zero_rank"] == 4,
      "canonical-direct common zero rank four")
check(common["canonical_direct_zero_null"] == "E",
      "canonical-direct common zero null E")
check(common["tensor_zero_rank"] == 3, "complete tensor zero rank three")
check(common["tensor_zero_null"] == "A1+E", "complete tensor zero null A1+E")
check(common["rigorous_near_zero_condition"] ==
      "rank(B_+(k)P_T)=3 for |k|^4<32/363",
      "rigorous near-zero ball")
check(common["T_normal_invariant_form_orthonormal_basis"] ==
      "8 I-2|k|^2 I+12 k k^T-28 diag(kx^2,ky^2,kz^2)+O(|k|^4)",
      "quadratic cubic-structured T normal")
check("cannot distinguish physical rotational anisotropy" in
      common["rotation_ceiling"], "T-only rotational ceiling")
check(common["det_invariant_form"] ==
      "768-1408|k|^2+1072|k|^4-(416/3)(kx^4+ky^4+kz^4)+O(|k|^6)",
      "common determinant series")
check(common["bz_corner"]["T_rank"] == 1, "exact BZ-corner rank loss")
check(relative_sector["zero_rank"] == 0, "relative zero rank")
check(relative_sector["leading_rank"] ==
      "3 generically; 2 on nonzero Cartesian face diagonals; 0 at k=0",
      "relative leading ranks")
check(len(relative_sector["exact_face_diagonal_dependencies"]) == 6,
      "six exact face-diagonal dependencies")
check(combined["zero_rank"] == 6, "combined common zero rank")
check(combined["zero_determinant"] == "524288", "combined determinant")
check(combined["C_definition"] == "C(k)=(D;B_+(k)P_T)",
      "combined map uses complete tensor writer only")
check(combined["unsoldered_zero_rank"] == 9, "unsoldered zero rank")
check(combined["unsoldered_all_k_rank_ceiling"] == 10,
      "unsoldered all-momentum ceiling")
check(scalar["identity"] == "sum_{a<b} Z_a Z_b=2(n-2)^2-2",
      "uniform pair identity")
check(scalar["denominator_derivative"] == "+(315/4)h^6/U_d^6",
      "uniform source denominator derivative")

theorem = " ".join((HERE / "THEOREM.md").read_text().split())
result = " ".join((HERE / "RESULT.md").read_text().split())
readme = " ".join((HERE / "README.md").read_text().split())
self_audit = " ".join((HERE / "SELF_AUDIT.md").read_text().split())
verification = " ".join((HERE / "VERIFICATION.txt").read_text().split())
for token in (
    "B_d^P(k)=\\sum_{p\\subset\\bar d}e^{-ik\\cdot\\rho_{d,p}}e_p^T",
    "B_d^+(k)=2\\sum_{p\\subset\\bar d}\\cos(k\\cdot\\rho_{d,p})e_p^T",
    "W_T(0)^*W_T(0)=8P_T",
    "C(0)^*C(0)=4P_A+16P_E+8P_T",
    "|k|^4<{32\\over363}",
    "N_T(k)=8I-2|k|^2I+12kk^T",
    "rank exactly nine",
    "common-field/soldering law",
    "canonical-direct bookkeeping",
    "U_d \\to U_d+2q",
    "same-parent storage-energy/future-writer linkage",
):
    check(token in theorem, f"theorem token: {token}")
check("that block alone does not decide physical rotational anisotropy" in readme,
      "README guards T-only rotational diagnosis")
check("A fully unsoldered twelve-direction" in result,
      "result exposes unsoldered ceiling")
check("caught an important possible overstatement" in self_audit,
      "author correction recorded")
check("PASS__GL6CL_GLOBAL_FOURIER_PAIR_WRITER__122/122" in verification,
      "science replay recorded")
check("PASS__GL6CL_PACKET__" in verification,
      "packet replay recorded")

for forbidden in (
    "is spacetime",
    "is a metric",
    "is Ricci",
    "is gravity",
    "proves gravity",
    "derives Einstein",
    "calculates G",
    "proves an autonomous source",
    "proves a stationary response",
):
    check(forbidden not in theorem + " " + result + " " + readme,
          f"forbidden promotion absent: {forbidden}")

print(f"PASS__GL6CL_PACKET__{checks}/{checks}")
