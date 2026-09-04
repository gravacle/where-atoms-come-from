#!/usr/bin/env python3
"""Fail-closed custody and replay verifier for the GL6CM hostile audit."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CM_GLOBAL_H6_STATIONARY_WRITER_RESPONSE_V001"
CHECKS = []


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON member {key!r}")
        answer[key] = value
    return answer


def strict_json(path):
    return json.loads(path.read_text(), object_pairs_hook=unique_object)


def check_hash_list(name):
    rows = []
    for line in (HERE / name).read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        path = ROOT / relative
        check(path.is_file(), f"{name} target exists {relative}")
        check(digest(path) == expected, f"{name} hash {relative}")
        rows.append(relative)
    check(len(rows) == len(set(rows)), f"{name} has unique paths")
    return rows


required = {
    "AUDIT_REPORT.md", "CONTEXT.sha256", "INDEPENDENT_RESULT.json",
    "MANIFEST.sha256", "README.md", "SEAL.sha256", "TARGET.sha256",
    "VERIFICATION.txt", "verify_gl6cm_independent.py", "verify_packet.py",
}
for name in sorted(required):
    check((HERE / name).is_file(), f"required audit file {name}")

target_rows = check_hash_list("TARGET.sha256")
check(len(target_rows) == 12, "all twelve target bytes are pinned")
check({Path(path).name for path in target_rows} ==
      {path.name for path in TARGET.iterdir() if path.is_file()},
      "target pin set is complete")
check(digest(TARGET / "MANIFEST.sha256") ==
      "61232da026192a717b09b903bdca9c72f0abfa7a00a002e301f6bb80873397b9",
      "requested stable target manifest hash")
check(digest(TARGET / "SEAL.sha256") ==
      "f9be064a6eea4630b548eafa4d56ea8c475f9a850fe267451ab6a5995ce37881",
      "stable target seal-file hash")
target_seal = (TARGET / "SEAL.sha256").read_text().strip().split("  ", 1)
check(target_seal == [digest(TARGET / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "target seal closes target manifest")

context_rows = check_hash_list("CONTEXT.sha256")
check(len(context_rows) == 14, "four upstream lanes and three hostile audits are context-pinned")

independent_source = (HERE / "verify_gl6cm_independent.py").read_text()
check("import subprocess" not in independent_source and
      "verify_global_h6_stationary_writer_response.py" not in independent_source,
      "independent science replay neither imports nor invokes target replay")
independent = subprocess.run(
    [sys.executable, "-B", str(HERE / "verify_gl6cm_independent.py")],
    cwd=ROOT,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=True,
)
for token in (
    "PASS finite-component PF scope: 771 exhaustive graphs plus general irreducibility proof",
    "PASS exact spectral Gram factorization and kernel: 625 vectors",
    "PASS exact common-rescaling null and W^T K W pullback: 2187 comparisons",
    "PASS isolated K2 zero and shared-star strict relative response",
    "PASS physical source scaling: 175sqrt(2)/32 h^6/U_d^7",
    "PASS claim typing: writer-only; contact, record authentication, bulk gravity and G open",
    "AUDIT DISPOSITION: PASS",
):
    check(token in independent.stdout, f"independent replay emits {token}")

author = subprocess.run(
    [sys.executable, "-B", str(TARGET / "verify_packet.py")],
    cwd=TARGET,
    text=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    check=True,
)
check("PASS GL6CM packet 104/104" in author.stdout,
      "frozen author packet replays 104/104")

result = strict_json(HERE / "INDEPENDENT_RESULT.json")
check(result["schema"] == "AUDIT_G_GL6CM_INDEPENDENT_V001" and
      result["disposition"] == "PASS", "frozen independent PASS disposition")
pf = result["perron_frobenius_scope"]
check(pf["labeled_connected_simple_graphs_n2_through_n5"] ==
      {"2": 1, "3": 4, "4": 38, "5": 728} and
      pf["total_checked"] == 771,
      "finite-component PF census is exact")
check(pf["constructive_pf_witness"] == "(I+A)^(n-1) is entrywise positive" and
      "unique ground ray" in pf["conclusion"],
      "finite PF scope and conclusion are explicit")

spectral = result["spectral_factorization"]
check(spectral["factorization"] ==
      "K=2 lambda_T^2 S^T diag(Delta_n^-1) S",
      "spectral Gram factorization is exact")
check(spectral["reciprocal"] and spectral["positive_semidefinite"] and
      spectral["kernel"] == "K(y,y)=0 iff S y=0 iff Q B_y|0>=0",
      "reciprocity, positivity, and exact kernel are frozen")
check(spectral["primary_transition_rank"] == 3 and
      spectral["primary_kernel_rank"] == 3 and
      spectral["extra_dark_kernel_example_rank"] == 2,
      "uniform null and possible extra dark modes are distinguished")
check(spectral["integer_vectors_checked"] == 625 and
      "<n|H0|0>=0" in spectral["common_rescaling_derivation"],
      "spectral vector census and parent-derived null are frozen")

pullback = result["writer_pullback"]
check(pullback["formula"] == "K_T=W^T K W" and
      pullback["writer_shape"] == [4, 6] and
      pullback["bilinear_exact_comparisons"] == 2187 and
      pullback["positive_semidefinite"],
      "writer pullback is independently exact")

components = result["finite_components"]
check(components["isolated_K2_spectral_response"] == "0",
      "isolated K2 spectral response vanishes")
check(components["star_spectrum_over_J"] == ["-sqrt(2)", "0", "+sqrt(2)"] and
      components["transition_rows"] ==
      {"middle": ["1/2", "-1/2"], "upper": ["0", "0"]},
      "star spectrum and transition rows are exact")
check(components["cycle_kernel"] ==
      "lambda_T^2 sqrt(2)/(4J) [[1,-1],[-1,1]]" and
      components["minus_energy_hessian"] ==
      "sqrt(2)(w0-w1)^2/(4J)",
      "star difference kernel and curvature are exact")
check(components["uniform_mode"] == "null" and
      components["relative_mode"] == "strict",
      "star null and strict directions are exact")

units = result["physical_scaling_and_units"]
check(units["literal_positive_response"] ==
      "(175sqrt(2)/32)h^6/U_d^7" and
      units["unit_positive_response"] ==
      "(175sqrt(2)/64)h^6/U_d^7" and
      units["spectral_only_inverse"] ==
      "(16sqrt(2)/175)U_d^7/h^6",
      "literal, unit, and inverse coefficients are exact")
check(units["lambda_T"].endswith("[dimensionless]") and
      units["cycle_kernel_entries"] == "inverse energy" and
      "inverse energy" in units["scalar_path"],
      "source and response units are consistent")

for name, finding in result["boundary_findings"].items():
    check(finding.startswith("PASS"), f"boundary finding {name} passes explicitly")

theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
for token in (
    "nontrivial connected component of a finite locked flip graph",
    "explicitly writer-only linear family",
    "possible diagonal order-six first-source vertex",
    "two-writer-vertex spectral contribution",
    "K\\mathbf1=0",
    "K_T^{\\rm spec}=W^TKW",
    "not permission to add one copy",
    "phi=dW/ds=-dE/ds",
    "record authentication",
    "neither a Ricci coefficient nor `G`",
):
    check(token in theorem, f"target theorem contains decisive typing guard: {token}")

report = " ".join((HERE / "AUDIT_REPORT.md").read_text().split())
for token in (
    "**Verdict:** **PASS**",
    "771 labeled connected simple graphs",
    "K=2\\lambda_T^2S^T\\operatorname{diag}(\\Delta_n^{-1})S",
    "K\\mathbf 1=0",
    "W^TKW",
    "K_{\\rm star}",
    "{175\\sqrt2\\over32}{h^6\\over U_d^7}",
    "Writer-only versus full source Hamiltonian",
    "Record authentication",
    "Bulk and gravity promotion",
):
    check(token in report, f"audit report contains decisive result: {token}")

verification = (HERE / "VERIFICATION.txt").read_text()
for token in ("DISPOSITION: PASS", "771", "2187", "isolated K2 zero",
              "175sqrt(2)/32", "record authentication", "gravity", "G open"):
    check(token in verification, f"verification summary contains {token}")

manifest_rows = []
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, name = line.split("  ", 1)
    check(Path(name).parent == Path("."), f"manifest path stays local {name}")
    check(digest(HERE / name) == expected, f"manifest hash {name}")
    manifest_rows.append(name)
check(set(manifest_rows) == required - {"MANIFEST.sha256", "SEAL.sha256"},
      "audit manifest covers every non-custody byte")
check(len(manifest_rows) == len(set(manifest_rows)), "audit manifest paths are unique")
seal = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal == [digest(HERE / "MANIFEST.sha256"), "MANIFEST.sha256"],
      "audit seal closes audit manifest")

print(f"PASS target custody: {len(target_rows)}/{len(target_rows)}")
print(f"PASS context custody: {len(context_rows)}/{len(context_rows)}")
print("PASS independent science: PF, spectral Gram/kernel, null, pullback, star, scaling")
print("PASS hostile boundaries: writer-only, contact, stationarity, authentication, accumulation, promotion")
print(f"PASS audit manifest: {len(manifest_rows)}/{len(manifest_rows)}")
print(f"PASS total checks: {len(CHECKS)}/{len(CHECKS)}")
print("AUDIT DISPOSITION: PASS")
