#!/usr/bin/env python3
"""Independent exact-arithmetic and custody verifier for CHGC V001."""

from __future__ import annotations

import ast
from decimal import Decimal as D, getcontext
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


getcontext().prec = 50
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0

MECHANICAL_CORRECTION_CENTRAL_PPM = {
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

MECHANICAL_CORRECTION_STANDARD_U_PPM = {
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


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print("PASS " + " ".join(label.split()))


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def parse_ledger(path):
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            value, name = line.split(maxsplit=1)
            result[name.strip()] = value
    return result


def dec(value):
    return D(str(value))


def close(left, right, relative=D("2e-15"), absolute=D("1e-30")):
    left, right = dec(left), dec(right)
    return abs(left - right) <= max(absolute, relative * max(abs(left), abs(right)))


# Dependency and transitive source custody.
dependencies = parse_ledger(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 12, "dependency ledger has exactly twelve entries")
for relative, expected in dependencies.items():
    path = (HERE / relative).resolve()
    check(path.is_file() and not path.is_symlink(),
          "dependency is a regular file: " + relative)
    check(digest(path) == expected, "dependency digest matches: " + relative)
    check(sha256(path.read_bytes() + b"tamper").hexdigest() != expected,
          "dependency appended-byte tamper fails: " + relative)

nominal_theorem = (ROOT / "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/THEOREM.md").read_text()
forward_theorem = (ROOT / "LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/THEOREM.md").read_text()
semantics = (HERE / "SOURCE_SEMANTICS.md").read_text()
theorem = (HERE / "THEOREM.md").read_text()
self_audit = (HERE / "SELF_AUDIT.md").read_text()
check("without using a measured value of \\(G\\), a gravity response, or a\npublished processed source coefficient" in nominal_theorem,
      "nominal parent declares response and processed-coefficient independence")
check("processed finite-mass source coefficient" in forward_theorem
      and "processed multipole source-response coefficient" in forward_theorem,
      "forward parent classifies both authors source coefficients as processed")
check("raw-like intermediate observations" in forward_theorem
      and "not\nthe original 0.5-second angle stream" in forward_theorem,
      "forward parent preserves the ToS released-data ceiling")
check("representative figure segment, not a\ncampaign-average" in forward_theorem,
      "forward parent withholds campaign binding from the AAF figure segment")
check("source-mass gravitational nonlinearity is\n   corrected synchronously" in semantics,
      "source semantics records ToS correction entanglement")
check("corrected for the air-density effect" in semantics,
      "source semantics records AAF response processing")
check("AAF-I/II \\(455.40(1.95)\\) ppm" in semantics
      and "ToS fibres" in semantics,
      "source semantics records Supplementary-Table-1 mechanical uncertainties")
check("uses these displayed central corrections" in semantics
      and "Extra digits reconstructed" in semantics,
      "source semantics rejects unowned mechanical-factor digits")


# Code-level comparator quarantine.
calculator_path = HERE / "calculate_conditional_homogeneous_g.py"
calculator_source = calculator_path.read_text(encoding="utf-8")
tree = ast.parse(calculator_source)
functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
check(all(name in functions for name in ("extract_primary_inputs", "compute_primary",
                                         "attach_post_comparators", "calculate")),
      "calculator has explicit primary extraction, primary compute, and post phases")
primary_nodes = (functions["extract_primary_inputs"], functions["compute_primary"])
primary_strings = {node.value for function in primary_nodes
                   for node in ast.walk(function)
                   if isinstance(node, ast.Constant) and isinstance(node.value, str)}
for forbidden in ("processed_coefficient_kg_m-3",
                  "published_G_summary_SI_comparison_only",
                  "recomputed_G_SI"):
    check(forbidden not in primary_strings,
          "primary AST excludes comparator key: " + forbidden)
post_strings = {node.value for node in ast.walk(functions["attach_post_comparators"])
                if isinstance(node, ast.Constant) and isinstance(node.value, str)}
for required in ("processed_coefficient_kg_m-3",
                 "published_G_summary_SI_comparison_only"):
    check(required in post_strings, "post phase owns comparator key: " + required)
check("CODATA" not in calculator_source and "accepted_G" not in calculator_source,
      "calculator imports no accepted or CODATA G field")


# Reproduction transcript: calculator output must equal the sealed result object.
run = subprocess.run([sys.executable, "-B", str(calculator_path)],
                     cwd=str(ROOT), check=True, capture_output=True, text=True)
generated = json.loads(run.stdout)
stored = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
check(generated == stored, "calculator exactly reproduces stored RESULT.json")

nominal = json.loads((ROOT / "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/RESULT.json").read_text())
forward = json.loads((ROOT / "LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/RESULT.json").read_text())
check(len(stored["AAF"]) == 3, "result contains all three AAF campaigns")
check(len(stored["TOS"]) == 7, "result contains all seven ToS rows")


# Independent AAF arithmetic.
nominal_aaf = {row["id"]: row for row in nominal["AAF"]}
forward_aaf = {row["id"]: row
               for row in forward["AAF_three_processed_coefficient_forwards"]}
for row in stored["AAF"]:
    n = nominal_aaf[row["id"]]
    f = forward_aaf[row["id"]]
    kernel = dec(n["nominal_homogeneous_coefficient_kg_m-3"])
    response = dec(f["alpha_nrad_s-2"]) * D("1e-9")
    u_response = dec(f["alpha_u_nrad_s-2"]) * D("1e-9")
    factor = D(1) + MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]] * D("1e-6")
    g_value = response * factor / kernel
    check(close(row["primary_conditional_G_SI"], g_value),
          row["id"] + " conditional AAF quotient is exact")
    check(row["response_alpha_nrad_s-2"] == f["alpha_nrad_s-2"],
          row["id"] + " response numerator has exact upstream custody")
    check(row["nominal_homogeneous_kernel_kg_m-3"]
          == n["nominal_homogeneous_coefficient_kg_m-3"],
          row["id"] + " homogeneous kernel has exact upstream custody")
    check(dec(row["mechanical_correction_central_ppm"])
          == MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]],
          row["id"] + " displayed mechanical-correction central value has table custody")
    check(close(row["mechanical_factor_held_at_displayed_value"], factor),
          row["id"] + " displayed AAF mechanical factor is 1 plus table correction")
    packet = row["local_uncertainty_diagnostics"]
    u_response_component = abs(g_value * u_response / response)
    u_kernel_component = abs(
        g_value * dec(n["public_input_sensitivity"]["rss_standard_uncertainty_kg_m-3"]) / kernel)
    u_mechanical = abs(
        g_value * MECHANICAL_CORRECTION_STANDARD_U_PPM[row["id"]]
        * D("1e-6") / factor)
    rss = (u_response_component**2 + u_kernel_component**2
           + u_mechanical**2).sqrt()
    check(close(packet["response_standard_component_SI"], u_response_component),
          row["id"] + " response uncertainty component is exact")
    check(close(packet["kernel_public_input_RSS_component_SI"], u_kernel_component),
          row["id"] + " kernel uncertainty component is exact")
    check(close(packet["mechanical_correction_standard_component_SI"], u_mechanical),
          row["id"] + " mechanical-correction uncertainty component is exact")
    check(dec(row["mechanical_correction_standard_u_ppm"])
          == MECHANICAL_CORRECTION_STANDARD_U_PPM[row["id"]],
          row["id"] + " mechanical-correction uncertainty has exact table custody")
    check(close(packet["zero_covariance_partial_RSS_SI"], rss),
          row["id"] + " partial zero-covariance RSS is exact")
    k_mix = dec(n["core_torque_divided_by_full_I_forbidden_mix_kg_m-3"])
    collision = response * factor / k_mix
    check(close(row["normalization_collision"]["conditional_G_if_missing_mass_has_zero_m2_torque_SI"], collision),
          row["id"] + " mixed-normalization diagnostic quotient is exact")
    post = row["post_calculation_comparator"]
    k_processed = dec(n["post_calculation_processed_comparator"]["processed_coefficient_kg_m-3"])
    g_processed_displayed = response * factor / k_processed
    check(close(post["processed_coefficient_displayed_factor_forward_G_SI"],
                g_processed_displayed),
          row["id"] + " displayed-factor processed comparator is exact")
    check(abs(dec(row["post_calculation_comparator"]["ratio_identity_residual"])) < D("5e-15"),
          row["id"] + " post comparator ratio identity closes")
    check(row["post_calculation_comparator"]["processed_minus_nominal_kernel_kg_m-3"] < 0,
          row["id"] + " processed-minus-nominal source gap retains negative sign")
    check("corrected for the air-density effect" in row["response_class"],
          row["id"] + " response is not relabelled raw")


# Independent ToS arithmetic and correction-family custody.
tos_map = {
    "TOS-I-F1-first": "fiber_1_first",
    "TOS-I-F1-repeat": "fiber_1_repeated",
    "TOS-I-F2": "fiber_2",
    "TOS-I-F3-first": "fiber_3_first",
    "TOS-I-F3-repeat": "fiber_3_repeated",
    "TOS-II-F4-first": "fiber_4_first",
    "TOS-II-F4-repeat": "fiber_4_repeated",
}
nominal_tos = {row["id"]: row for row in nominal["TOS"]}
forward_tos = {row["id"]: row
               for row in forward["ToS_seven_processed_coefficient_forwards"]}
for row in stored["TOS"]:
    n = nominal_tos[row["id"]]
    f = forward_tos[tos_map[row["id"]]]
    kernel = dec(n["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"])
    response = dec(f["mean_delta_omega2_s-2"])
    u_response = dec(f["mean_delta_omega2_standard_u_s-2"])
    factor = D(1) + MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]] * D("1e-6")
    g_anchor = response * factor / kernel
    slope = response / kernel
    check(close(row["primary_magnetic_only_G_anchor_SI"], g_anchor),
          row["id"] + " ToS magnetic-only anchor is exact")
    check(close(row["primary_affine_family"]["slope_dG_dc_f_SI"], slope),
          row["id"] + " ToS affine correction slope is exact")
    check(row["primary_affine_family"]["public_identified_interval"] is None,
          row["id"] + " ToS interval remains unidentified")
    check(dec(row["mechanical_correction_central_ppm"])
          == MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]],
          row["id"] + " displayed ToS mechanical-correction central value has table custody")
    check(close(D(1) + dec(row["magnetic_factor_minus_one_held_at_displayed_value"]),
                factor),
          row["id"] + " displayed ToS mechanical factor is 1 plus table correction")
    packet = row["local_uncertainty_diagnostics_at_c_f_zero"]
    u_n = abs(g_anchor * u_response / response)
    u_k = abs(g_anchor * dec(n["public_input_sensitivity"]["rss_standard_uncertainty_kg_m-3"]) / kernel)
    u_mechanical = abs(
        g_anchor * MECHANICAL_CORRECTION_STANDARD_U_PPM[row["id"]]
        * D("1e-6") / factor)
    check(close(packet["mechanical_correction_standard_component_SI"],
                u_mechanical),
          row["id"] + " ToS mechanical-correction uncertainty is exact")
    check(dec(row["mechanical_correction_standard_u_ppm"])
          == MECHANICAL_CORRECTION_STANDARD_U_PPM[row["id"]],
          row["id"] + " ToS mechanical uncertainty has exact table custody")
    check(close(packet["zero_covariance_partial_RSS_SI"],
                (u_n**2 + u_k**2 + u_mechanical**2).sqrt()),
          row["id"] + " ToS partial RSS is exact")
    k_mix = dec(n["core_curvature_divided_by_full_I_forbidden_mix_kg_m-3"])
    check(close(row["normalization_collision_at_c_f_zero"]["conditional_G_SI"],
                response * factor / k_mix),
          row["id"] + " ToS mixed-normalization diagnostic is exact")
    post = row["post_calculation_comparator"]
    k_processed = dec(n["post_calculation_processed_comparator"]["processed_coefficient_kg_m-3"])
    g_published = dec(f["published_G_summary_SI_comparison_only"])
    c_total = g_published * kernel / response - factor
    c_dynamic = g_published * k_processed / response - factor
    check(close(post["c_total_required_with_homogeneous_kernel_ppm"],
                c_total * D(1000000)),
          row["id"] + " homogeneous total required correction is exact")
    check(close(post["c_dynamic_required_with_processed_kernel_ppm"],
                c_dynamic * D(1000000)),
          row["id"] + " processed-kernel dynamic bracket is exact")
    check("not a source-model-free raw numerator" in row["response_class"],
          row["id"] + " ToS corrected response ceiling is explicit")


# Figure-level numerator boundaries.
aaf_figure = stored["AAF_figure_level_unbound_response"]
check(aaf_figure["conditional_G_withheld"] is True
      and aaf_figure["campaign_kernel_binding"] is None,
      "AAF figure response is not converted into an unbound fourth quotient")
check(close(aaf_figure["two_hour_source_harmonic_nrad_s-2"],
            forward["AAF_figure_level_acceleration_stream"]["two_tone_source_amplitude_nrad_s-2"]),
      "AAF figure response has exact upstream custody")
figure = stored["TOS_figure_level_released_response_diagnostic"]
repeat_n = nominal_tos["TOS-I-F1-repeat"]
repeat_f = forward_tos["fiber_1_repeated"]
kernel = dec(repeat_n["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"])
factor = D(1) + MECHANICAL_CORRECTION_CENTRAL_PPM["TOS-I-F1-repeat"] * D("1e-6")
forward_figure = forward["ToS_figure_level_response"]
check(close(figure["A_B_A_magnetic_only_G_anchor_SI"],
            dec(forward_figure["A_B_A_background_subtracted_delta_omega2_s-2"]) * factor / kernel),
      "ToS A-B-A figure anchor is exact")
check(close(figure["quadratic_drift_diagnostic_G_anchor_SI"],
            dec(forward_figure["common_quadratic_background_subtracted_delta_omega2_s-2"]) * factor / kernel),
      "ToS quadratic figure anchor is exact")


# Documentary claim ceilings.
required_theorem_phrases = (
    "Setting \\(r_{\\mathrm{norm}}=0\\) is an",
    "None is a source-model-free raw numerator",
    "No fourth \\(G\\)-quotient is formed",
    "There is no independently identified deterministic compact interval",
    "mixed value is not asserted to be a physically admissible second HUST",
    "specifies no bounded admissible domain or covariance law",
    "not a claim of mathematical unboundedness",
    "not a denial of the authors' published",
    "u_{G,f}=|G|\\frac{u_f}{|f|}",
    "inserting the processed\n  coefficient merely evaluates the processed-model forward",
    "not a new or independent measurement of \\(G\\)",
    "Further fitting of the\nprocessed \\(G\\) rows cannot supply those fields",
)
for phrase in required_theorem_phrases:
    check(phrase in theorem, "theorem preserves claim ceiling: " + phrase)
check("Accepted \\(G\\) imported?** No" in self_audit,
      "self-audit rejects accepted-G import")
check("Core numerator divided by full inertia promoted to a result?** No" in self_audit,
      "self-audit keeps mixed-normalization diagnostic non-promotional")
check("Mechanical-factor uncertainty omitted or invented?** No" in self_audit,
      "self-audit owns published magnetic-correction uncertainty")
check("not RGRL or Gravity Formation Theory confirmation" in stored["claim_ceiling"]
      and "authors' processed-model G summaries" in stored["claim_ceiling"],
      "machine result withholds GFT confirmation")

for relative in (
        "README.md", "SELF_AUDIT.md", "SOURCE_SEMANTICS.md", "THEOREM.md",
        "calculate_conditional_homogeneous_g.py",
        "verify_conditional_homogeneous_g.py"):
    data = (HERE / relative).read_bytes()
    check(not any(byte < 32 and byte not in (9, 10) for byte in data),
          "core lane text has no forbidden control bytes: " + relative)
check("\\rm" not in theorem and "\\rm" not in self_audit,
      "theorem and self-audit use no unsafe TeX rm sequence")


# Stable payload manifest, when present.
manifest_path = HERE / "MANIFEST.sha256"
if manifest_path.is_file():
    manifest = parse_ledger(manifest_path)
    for relative, expected in manifest.items():
        check(digest(HERE / relative) == expected,
              "lane manifest digest matches: " + relative)

print("SUMMARY {0}/{0} exact checks passed".format(checks))
print("DISPOSITION AAF_RNORM0_CONDITIONAL_QUOTIENTS__TOS_AFFINE_CORRECTION_FAMILIES__NO_PUBLIC_INDEPENDENT_DETERMINISTIC_G_INTERVAL__NO_NEW_G")
