#!/usr/bin/env python3
"""Fresh hostile replay of CHGC V001 without importing either lane verifier.

This executable treats the theorem, calculator, stored result, and production
verifier as untrusted.  It independently reads the pinned parents and official
Supplement, reconstructs every quotient and local uncertainty component, and
audits claim scope and sealing.
"""

from __future__ import annotations

import ast
from decimal import Decimal as D, getcontext
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys


getcontext().prec = 60
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0

PDF_REL = "../LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/SOURCE/41586_2018_431_MOESM1_ESM.pdf"
PDF_SHA256 = "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb"

CENTRAL_PPM = {
    "AAF-I": D("455.40"),
    "AAF-II": D("455.40"),
    "AAF-III": D("25.74"),
    "TOS-I-F1-first": D("0.47"),
    "TOS-I-F1-repeat": D("0.47"),
    "TOS-I-F2": D("7.13"),
    "TOS-I-F3-first": D("0.32"),
    "TOS-I-F3-repeat": D("0.32"),
    "TOS-II-F4-first": D("0.27"),
    "TOS-II-F4-repeat": D("0.27"),
}

STANDARD_U_PPM = {
    "AAF-I": D("1.95"),
    "AAF-II": D("1.95"),
    "AAF-III": D("0.08"),
    "TOS-I-F1-first": D("0.08"),
    "TOS-I-F1-repeat": D("0.08"),
    "TOS-I-F2": D("1.19"),
    "TOS-I-F3-first": D("0.05"),
    "TOS-I-F3-repeat": D("0.05"),
    "TOS-II-F4-first": D("0.08"),
    "TOS-II-F4-repeat": D("0.08"),
}

TOS_FORWARD_ID = {
    "TOS-I-F1-first": "fiber_1_first",
    "TOS-I-F1-repeat": "fiber_1_repeated",
    "TOS-I-F2": "fiber_2",
    "TOS-I-F3-first": "fiber_3_first",
    "TOS-I-F3-repeat": "fiber_3_repeated",
    "TOS-II-F4-first": "fiber_4_first",
    "TOS-II-F4-repeat": "fiber_4_repeated",
}

# Supplementary Table 1 central component values.  They are used only to
# challenge the signs and rounding of the table's directly reported Delta G/G.
COMPONENTS = {
    "AAF-I": (D("2.401e-5"), D("1.199e-5"), D("2.776e-5"), D("6.313e-9"), "AAF"),
    "AAF-II": (D("2.401e-5"), D("1.199e-5"), D("2.776e-5"), D("6.313e-9"), "AAF"),
    "AAF-III": (D("2.404e-5"), D("21.24e-5"), D("2.776e-5"), D("6.313e-9"), "AAF"),
    "TOS-I-F1-first": (D("2.18e-5"), D("1.2e-5"), D("4.7705e-5"), D("12.2e-9"), "TOS"),
    "TOS-I-F1-repeat": (D("2.18e-5"), D("1.2e-5"), D("4.7705e-5"), D("12.2e-9"), "TOS"),
    "TOS-I-F2": (D("2.18e-5"), D("1.2e-5"), D("4.7706e-5"), D("47.4e-9"), "TOS"),
    "TOS-I-F3-first": (D("2.18e-5"), D("1.2e-5"), D("4.7705e-5"), D("10.1e-9"), "TOS"),
    "TOS-I-F3-repeat": (D("2.18e-5"), D("1.2e-5"), D("4.7705e-5"), D("10.1e-9"), "TOS"),
    "TOS-II-F4-first": (D("2.21e-5"), D("1.4e-5"), D("4.6477e-5"), D("10.6e-9"), "TOS"),
    "TOS-II-F4-repeat": (D("2.21e-5"), D("1.4e-5"), D("4.6477e-5"), D("10.6e-9"), "TOS"),
}

EXPECTED_MANIFEST_MEMBERS = {
    "DEPENDENCIES.sha256",
    "HOSTILE_AUDIT_TRANSCRIPT.txt",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "README.md",
    "RESULT.json",
    "SELF_AUDIT.md",
    "SOURCE_SEMANTICS.md",
    "THEOREM.md",
    "VERIFICATION.txt",
    "calculate_conditional_homogeneous_g.py",
    "hostile_audit_conditional_homogeneous_g.py",
    "verify_conditional_homogeneous_g.py",
}


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print("PASS " + " ".join(label.split()))


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def ledger(path):
    parsed = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            value, name = line.split(maxsplit=1)
            parsed[name.strip()] = value
    return parsed


def dec(value):
    return D(str(value))


def close(left, right, relative=D("2e-15"), absolute=D("1e-30")):
    left, right = dec(left), dec(right)
    return abs(left - right) <= max(absolute, relative * max(abs(left), abs(right)))


# 1. Transitive custody and direct official-source semantics.
dependencies = ledger(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 12, "dependency ledger has exactly twelve pinned files")
check(dependencies.get(PDF_REL) == PDF_SHA256,
      "dependency ledger pins the official Supplement PDF")
for relative, expected in dependencies.items():
    path = (HERE / relative).resolve()
    check(path.is_file() and not path.is_symlink(),
          "dependency is a regular non-symlink file: " + relative)
    check(digest(path) == expected, "dependency digest closes: " + relative)
    check(sha256(path.read_bytes() + b"hostile-tamper").hexdigest() != expected,
          "dependency appended-byte tamper is rejected: " + relative)

pdf_path = (HERE / PDF_REL).resolve()
check(digest(pdf_path) == PDF_SHA256, "official Supplement has the independently pinned hash")
pdftotext = shutil.which("pdftotext")
check(pdftotext is not None, "Poppler pdftotext is available for independent source replay")
pdf_run = subprocess.run([pdftotext, "-layout", str(pdf_path), "-"],
                         check=True, capture_output=True, text=True)
pdf_flat = " ".join(pdf_run.stdout.split())
for token in ("455.40(1.95)", "25.74(8)", "0.47(8)", "7.13(1.19)",
              "0.32(5)", "0.27(8)"):
    check(token in pdf_flat, "official Table 1 contains correction field " + token)
check("corrected synchronously in determining" in pdf_flat,
      "official source says source-gravity nonlinearity enters corrected ToS response")
check("corrected for the air density effect" in pdf_flat,
      "official source says AAF campaign acceleration is air-density corrected")
check("sampling rate of the DAS is 20 kHz" in pdf_flat
      and "sampling interval is 1 s" in pdf_flat,
      "official source distinguishes 20-kHz acquisition from released one-second stream")


# 2. Challenge mechanical equations, signs, and source-owned precision.
for row_id, (im, km, inertia, stiffness, method) in COMPONENTS.items():
    if method == "AAF":
        reconstructed_ppm = (stiffness / km) * (im / inertia) * D("1e6")
    else:
        reconstructed_ppm = im * stiffness * stiffness / (inertia * km * km) * D("1e6")
    check(reconstructed_ppm > 0, row_id + " mechanical correction has positive source sign")
    check(abs(reconstructed_ppm - CENTRAL_PPM[row_id]) <= D("0.005"),
          row_id + " component equation rounds within half a 0.01-ppm display unit")


# 3. AST quarantine and absence of accepted-G inputs.
calculator_path = HERE / "calculate_conditional_homogeneous_g.py"
source = calculator_path.read_text(encoding="utf-8")
tree = ast.parse(source)
functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
for name in ("extract_primary_inputs", "compute_primary", "attach_post_comparators", "calculate"):
    check(name in functions, "calculator exposes phase function " + name)
primary_strings = {
    node.value
    for name in ("extract_primary_inputs", "compute_primary")
    for node in ast.walk(functions[name])
    if isinstance(node, ast.Constant) and isinstance(node.value, str)
}
for key in ("processed_coefficient_kg_m-3",
            "published_G_summary_SI_comparison_only", "recomputed_G_SI"):
    check(key not in primary_strings, "primary AST never selects comparator key " + key)
post_strings = {
    node.value for node in ast.walk(functions["attach_post_comparators"])
    if isinstance(node, ast.Constant) and isinstance(node.value, str)
}
for key in ("processed_coefficient_kg_m-3", "published_G_summary_SI_comparison_only"):
    check(key in post_strings, "post phase alone selects comparator key " + key)
check("accepted_G" not in source and "CODATA_G" not in source,
      "calculator contains no accepted-G or CODATA-G input field")
check(not any("codata" in name.lower() or "accepted" in name.lower()
              for name in dependencies),
      "dependency ledger contains no accepted-G or CODATA payload")

calculate_body = functions["calculate"].body
assignment_names = [node.targets[0].id for node in calculate_body
                    if isinstance(node, ast.Assign)
                    and isinstance(node.targets[0], ast.Name)]
check(assignment_names == ["nominal", "forward", "inputs", "primary"],
      "calculate phase order parses parents, extracts inputs, then freezes primary")
return_node = calculate_body[-1]
check(isinstance(return_node, ast.Return)
      and isinstance(return_node.value, ast.Call)
      and isinstance(return_node.value.func, ast.Name)
      and return_node.value.func.id == "attach_post_comparators",
      "post comparator attachment is the final calculation operation")


# 4. Independent numerical replay from pinned parent objects.
nominal = json.loads((ROOT / "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/RESULT.json").read_text())
forward = json.loads((ROOT / "LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/RESULT.json").read_text())
stored = json.loads((HERE / "RESULT.json").read_text())
check(stored["schema"] == "WAC_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001",
      "stored result has the expected schema")
check(len(stored["AAF"]) == 3 and len(stored["TOS"]) == 7,
      "stored result covers three AAF campaigns and seven ToS rows")

nominal_aaf = {row["id"]: row for row in nominal["AAF"]}
forward_aaf = {row["id"]: row for row in forward["AAF_three_processed_coefficient_forwards"]}
for row in stored["AAF"]:
    row_id = row["id"]
    n = nominal_aaf[row_id]
    f = forward_aaf[row_id]
    response = dec(f["alpha_nrad_s-2"]) * D("1e-9")
    u_response = dec(f["alpha_u_nrad_s-2"]) * D("1e-9")
    kernel = dec(n["nominal_homogeneous_coefficient_kg_m-3"])
    u_kernel = dec(n["public_input_sensitivity"]["rss_standard_uncertainty_kg_m-3"])
    factor = D(1) + CENTRAL_PPM[row_id] * D("1e-6")
    g_value = response * factor / kernel
    check(close(row["primary_conditional_G_SI"], g_value),
          row_id + " AAF conditional quotient independently closes")
    check(dec(row["mechanical_correction_central_ppm"]) == CENTRAL_PPM[row_id]
          and close(row["mechanical_factor_held_at_displayed_value"], factor),
          row_id + " AAF uses direct displayed central correction")
    packet = row["local_uncertainty_diagnostics"]
    u_n = abs(g_value * u_response / response)
    u_k = abs(g_value * u_kernel / kernel)
    u_f = abs(g_value * STANDARD_U_PPM[row_id] * D("1e-6") / factor)
    check(close(packet["response_standard_component_SI"], u_n),
          row_id + " AAF response uncertainty propagation closes")
    check(close(packet["kernel_public_input_RSS_component_SI"], u_k),
          row_id + " AAF kernel uncertainty propagation closes")
    check(close(packet["mechanical_correction_standard_component_SI"], u_f),
          row_id + " AAF magnetic-correction uncertainty propagation closes")
    check(close(packet["zero_covariance_partial_RSS_SI"],
                (u_n * u_n + u_k * u_k + u_f * u_f).sqrt()),
          row_id + " AAF zero-covariance partial RSS closes")
    k_mix = dec(n["core_torque_divided_by_full_I_forbidden_mix_kg_m-3"])
    mix_g = response * factor / k_mix
    check(close(row["normalization_collision"]["conditional_G_if_missing_mass_has_zero_m2_torque_SI"], mix_g),
          row_id + " AAF mixed-normalization diagnostic closes")
    check("not a physical apparatus coefficient" in row["normalization_collision"]["meaning"],
          row_id + " AAF mixed diagnostic is not promoted to a realized map")
    post = row["post_calculation_comparator"]
    k_processed = dec(n["post_calculation_processed_comparator"]["processed_coefficient_kg_m-3"])
    displayed_forward = response * factor / k_processed
    check(close(post["processed_coefficient_displayed_factor_forward_G_SI"], displayed_forward),
          row_id + " AAF post-only displayed-factor forward closes")
    identity = dec(row["primary_conditional_G_SI"]) / displayed_forward - k_processed / kernel
    check(abs(identity) < D("5e-15")
          and abs(dec(post["ratio_identity_residual"])) < D("5e-15"),
          row_id + " AAF post-only kernel-ratio identity closes")
    check(post["upstream_component_reconstructed_processed_forward_G_SI"] == f["recomputed_G_SI"],
          row_id + " upstream extra-digit forward remains comparator-only")

nominal_tos = {row["id"]: row for row in nominal["TOS"]}
forward_tos = {row["id"]: row for row in forward["ToS_seven_processed_coefficient_forwards"]}
for row in stored["TOS"]:
    row_id = row["id"]
    n = nominal_tos[row_id]
    f = forward_tos[TOS_FORWARD_ID[row_id]]
    response = dec(f["mean_delta_omega2_s-2"])
    u_response = dec(f["mean_delta_omega2_standard_u_s-2"])
    kernel = dec(n["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"])
    u_kernel = dec(n["public_input_sensitivity"]["rss_standard_uncertainty_kg_m-3"])
    factor = D(1) + CENTRAL_PPM[row_id] * D("1e-6")
    anchor = response * factor / kernel
    slope = response / kernel
    check(close(row["primary_magnetic_only_G_anchor_SI"], anchor),
          row_id + " ToS magnetic-only anchor independently closes")
    check(close(row["primary_affine_family"]["slope_dG_dc_f_SI"], slope),
          row_id + " ToS affine c_f slope independently closes")
    check(row["primary_affine_family"]["c_f_definition"]
          == "c_f=-Delta_K/(I*Delta_omega2)",
          row_id + " ToS anelastic correction retains the source minus sign")
    check(row["primary_affine_family"]["public_identified_interval"] is None,
          row_id + " ToS family does not invent a deterministic interval")
    packet = row["local_uncertainty_diagnostics_at_c_f_zero"]
    u_n = abs(anchor * u_response / response)
    u_k = abs(anchor * u_kernel / kernel)
    u_f = abs(anchor * STANDARD_U_PPM[row_id] * D("1e-6") / factor)
    check(close(packet["response_standard_component_SI"], u_n),
          row_id + " ToS response uncertainty propagation closes")
    check(close(packet["kernel_public_input_RSS_component_SI"], u_k),
          row_id + " ToS kernel uncertainty propagation closes")
    check(close(packet["mechanical_correction_standard_component_SI"], u_f),
          row_id + " ToS magnetic-correction uncertainty propagation closes")
    check(close(packet["zero_covariance_partial_RSS_SI"],
                (u_n * u_n + u_k * u_k + u_f * u_f).sqrt()),
          row_id + " ToS zero-covariance partial RSS closes")
    k_mix = dec(n["core_curvature_divided_by_full_I_forbidden_mix_kg_m-3"])
    check(close(row["normalization_collision_at_c_f_zero"]["conditional_G_SI"],
                response * factor / k_mix),
          row_id + " ToS mixed-normalization diagnostic closes")
    check("not a physical apparatus coefficient" in row["normalization_collision_at_c_f_zero"]["meaning"],
          row_id + " ToS mixed diagnostic is not promoted to a realized map")
    post = row["post_calculation_comparator"]
    k_processed = dec(n["post_calculation_processed_comparator"]["processed_coefficient_kg_m-3"])
    g_published = dec(f["published_G_summary_SI_comparison_only"])
    c_total = g_published * kernel / response - factor
    c_dynamic = g_published * k_processed / response - factor
    check(close(post["c_total_required_with_homogeneous_kernel_ppm"], c_total * D("1e6")),
          row_id + " ToS homogeneous-kernel post-comparator identity closes")
    check(close(post["c_dynamic_required_with_processed_kernel_ppm"], c_dynamic * D("1e6")),
          row_id + " ToS processed-kernel post-comparator identity closes")


# 5. Campaign binding and non-circularity, without independence overclaim.
aaf_figure = stored["AAF_figure_level_unbound_response"]
check(aaf_figure["conditional_G_withheld"] is True
      and aaf_figure["campaign_kernel_binding"] is None,
      "representative AAF figure stream is not turned into a fourth campaign quotient")
check(aaf_figure["two_hour_source_harmonic_nrad_s-2"]
      == forward["AAF_figure_level_acceleration_stream"]["two_tone_source_amplitude_nrad_s-2"],
      "AAF figure harmonic retains inherited released-file custody")

figure = stored["TOS_figure_level_released_response_diagnostic"]
repeat_id = "TOS-I-F1-repeat"
repeat_n = nominal_tos[repeat_id]
repeat_f = forward_tos[TOS_FORWARD_ID[repeat_id]]
kernel = dec(repeat_n["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"])
factor = D(1) + CENTRAL_PPM[repeat_id] * D("1e-6")
forward_figure = forward["ToS_figure_level_response"]
check(figure["row_binding"] == repeat_id,
      "ToS figure diagnostic preserves inherited repeated-fibre-1 binding")
check(close(figure["A_B_A_magnetic_only_G_anchor_SI"],
            dec(forward_figure["A_B_A_background_subtracted_delta_omega2_s-2"])
            * factor / kernel),
      "ToS A-B-A released-response diagnostic closes")
check(close(figure["quadratic_drift_diagnostic_G_anchor_SI"],
            dec(forward_figure["common_quadratic_background_subtracted_delta_omega2_s-2"])
            * factor / kernel),
      "ToS quadratic released-response diagnostic closes")
check("released_file_level_noncircularity" in stored["numerator_custody"]
      and "not statistical or model independence" in
      stored["numerator_custody"]["released_file_level_noncircularity"],
      "result labels accepted-G nonuse as non-circularity, not independence")
check("corrected response summaries" in stored["numerator_custody"]["campaign_response_ceiling"]
      and "not fully independent measurements" in
      stored["numerator_custody"]["campaign_response_ceiling"],
      "result preserves corrected-summary and shared-model ceiling")


# 6. Exact documentary claim scope.
theorem = (HERE / "THEOREM.md").read_text(encoding="utf-8")
semantics = (HERE / "SOURCE_SEMANTICS.md").read_text(encoding="utf-8")
self_audit = (HERE / "SELF_AUDIT.md").read_text(encoding="utf-8")
for phrase in (
    "computational\nnonuse, not statistical independence",
    "c_f=-\\frac{\\Delta K}{I\\Delta\\omega^2}",
    "(1+\\delta_m+c_f)",
    "not asserted to be a physically admissible second HUST\napparatus map",
    "no deterministic compact interval for\nthis conditional family",
    "not a claim of mathematical unboundedness",
    "not a denial of the authors' published",
    "zero cross covariance",
    "not coverage intervals",
):
    check(phrase in theorem, "theorem preserves hostile-audited scope: " + phrase)
check("Extra digits reconstructed" in semantics
      and "not promoted to source-owned precision" in semantics,
      "source semantics limits correction precision to the displayed table")
check("key-level nonuse, not byte-level" in semantics,
      "source semantics states exact comparator quarantine boundary")
check("Released-file nonuse called statistical independence?** No" in self_audit,
      "self-audit rejects numerator-independence overclaim")
check("Authors' published summaries denied?** No" in self_audit,
      "self-audit acknowledges authors' processed-model summaries")
check("zero_covariance" in json.dumps(stored).lower()
      and "coverage interval" in json.dumps(stored).lower(),
      "machine result does not promote marginal propagation to joint coverage")
check("Combined" not in json.dumps(stored),
      "machine result forms no cross-row combined G with absent covariance")


# 7. Reproduction, byte/TeX hygiene, exact payload, and seal.
run = subprocess.run([sys.executable, "-B", str(calculator_path)], cwd=str(ROOT),
                     check=True, capture_output=True, text=True)
check(json.loads(run.stdout) == stored,
      "untrusted calculator independently replays to the stored result object")

for relative in sorted(EXPECTED_MANIFEST_MEMBERS):
    data = (HERE / relative).read_bytes()
    check(not any(byte < 32 and byte not in (9, 10) for byte in data),
          "payload has no forbidden control bytes: " + relative)
for relative in ("THEOREM.md", "SELF_AUDIT.md", "SOURCE_SEMANTICS.md",
                 "INDEPENDENT_HOSTILE_AUDIT.md"):
    text = (HERE / relative).read_text(encoding="utf-8")
    check("\\rm" not in text, "document has no unsafe TeX rm sequence: " + relative)

manifest = ledger(HERE / "MANIFEST.sha256")
check(set(manifest) == EXPECTED_MANIFEST_MEMBERS,
      "manifest has the exact hostile-audited payload member set")
for relative, expected in manifest.items():
    check(digest(HERE / relative) == expected,
          "manifest member digest closes: " + relative)
seal = ledger(HERE / "LANE_SEAL.sha256")
check(seal == {
          "MANIFEST.sha256": digest(HERE / "MANIFEST.sha256"),
          "VERIFICATION.txt": digest(HERE / "VERIFICATION.txt"),
      }, "lane seal closes exactly over manifest and production transcript bytes")
check(sha256((HERE / "MANIFEST.sha256").read_bytes() + b"tamper").hexdigest()
      != seal["MANIFEST.sha256"],
      "lane-seal appended-byte tamper sentinel rejects mutation")

print("SUMMARY {0}/{0} hostile checks passed".format(checks))
print("DISPOSITION ACCEPT_CONDITIONAL_QUOTIENT_AND_PUBLIC_PACKET_NONIDENTIFICATION__NOT_NEW_G__NOT_GFT_EVIDENCE")
