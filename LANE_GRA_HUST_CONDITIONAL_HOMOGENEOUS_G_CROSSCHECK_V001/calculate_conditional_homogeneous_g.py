#!/usr/bin/env python3
"""Conditional HUST homogeneous-kernel G cross-check.

Primary calculations use only the independently reconstructed nominal kernels,
their public-input sensitivities, measured response summaries, and displayed
mechanical corrections.  Published processed source coefficients and derived
G comparator keys are selected only in attach_post_comparators().
"""

from __future__ import annotations

from copy import deepcopy
from decimal import Decimal as D, getcontext
from hashlib import sha256
import json
from pathlib import Path


getcontext().prec = 50
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

NOMINAL_RESULT = ROOT / "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/RESULT.json"
FORWARD_RESULT = ROOT / "LANE_GRA_HUST_2018_DUAL_METHOD_G_FORWARD_V001/RESULT.json"

PINNED = {
    NOMINAL_RESULT: "1dc31ef45bbd5d8441867b2575e2dec9d5290ee6f1a57f71781f80e0730cd875",
    FORWARD_RESULT: "fb503bf0082d7dd1ff7099485020b2a7b364917ac9c43a414ef07e97bd7a4fbb",
}

# Official Supplementary Table 1, displayed central magnetic-damper correction
# Delta G/G and its quoted one-standard-deviation uncertainty, both in ppm.
# These displayed values, not extra digits reconstructed from component fields,
# are the source-owned correction inputs for this lane.  The official SI PDF is
# directly pinned in DEPENDENCIES.sha256.
MECHANICAL_CORRECTION_CENTRAL_PPM = {
    "AAF-I": "455.40",
    "AAF-II": "455.40",
    "AAF-III": "25.74",
    "TOS-I-F1-first": "0.47",
    "TOS-I-F1-repeat": "0.47",
    "TOS-I-F2": "7.13",
    "TOS-I-F3-first": "0.32",
    "TOS-I-F3-repeat": "0.32",
    "TOS-II-F4-first": "0.27",
    "TOS-II-F4-repeat": "0.27",
}

MECHANICAL_CORRECTION_STANDARD_U_PPM = {
    "AAF-I": "1.95",
    "AAF-II": "1.95",
    "AAF-III": "0.08",
    "TOS-I-F1-first": "0.08",
    "TOS-I-F1-repeat": "0.08",
    "TOS-I-F2": "1.19",
    "TOS-I-F3-first": "0.05",
    "TOS-I-F3-repeat": "0.05",
    "TOS-II-F4-first": "0.08",
    "TOS-II-F4-repeat": "0.08",
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(value):
    return D(str(value))


def as_float(value):
    return float(value)


def relative_ppm(left, right):
    return as_float((left / right - D(1)) * D(1000000))


def validate_pins():
    for path, expected in PINNED.items():
        if digest(path) != expected:
            raise RuntimeError("upstream dependency hash mismatch: " + str(path))


def extract_primary_inputs(nominal, forward):
    """Select only response, nominal-kernel, and displayed-transfer fields."""
    aaf_forward = {row["id"]: row
                   for row in forward["AAF_three_processed_coefficient_forwards"]}
    aaf = []
    for row in nominal["AAF"]:
        observed = aaf_forward[row["id"]]
        aaf.append({
            "id": row["id"],
            "K_hom": row["nominal_homogeneous_coefficient_kg_m-3"],
            "u_K_rss": row["public_input_sensitivity"]["rss_standard_uncertainty_kg_m-3"],
            "K_axis_box": row["public_input_sensitivity"]["linearized_axis_box_half_width_kg_m-3"],
            "K_core_over_full_I": row["core_torque_divided_by_full_I_forbidden_mix_kg_m-3"],
            "I0": row["uniform_cuboid_I0_kg_m2"],
            "I_full": row["supplement_full_apparatus_I_kg_m2"],
            "response_nrad_s2": observed["alpha_nrad_s-2"],
            "u_response_nrad_s2": observed["alpha_u_nrad_s-2"],
            "mechanical_factor": as_float(
                D(1) + dec(MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]])
                * D("1e-6")),
            "mechanical_correction_central_ppm":
                MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]],
            "mechanical_correction_standard_u_ppm":
                MECHANICAL_CORRECTION_STANDARD_U_PPM[row["id"]],
        })

    tos_forward_ids = {
        "TOS-I-F1-first": "fiber_1_first",
        "TOS-I-F1-repeat": "fiber_1_repeated",
        "TOS-I-F2": "fiber_2",
        "TOS-I-F3-first": "fiber_3_first",
        "TOS-I-F3-repeat": "fiber_3_repeated",
        "TOS-II-F4-first": "fiber_4_first",
        "TOS-II-F4-repeat": "fiber_4_repeated",
    }
    tos_forward = {row["id"]: row
                   for row in forward["ToS_seven_processed_coefficient_forwards"]}
    tos = []
    for row in nominal["TOS"]:
        observed = tos_forward[tos_forward_ids[row["id"]]]
        tos.append({
            "id": row["id"],
            "forward_id": observed["id"],
            "K_hom": row["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"],
            "u_K_rss": row["public_input_sensitivity"]["rss_standard_uncertainty_kg_m-3"],
            "K_axis_box": row["public_input_sensitivity"]["linearized_axis_box_half_width_kg_m-3"],
            "K_core_over_full_I": row["core_curvature_divided_by_full_I_forbidden_mix_kg_m-3"],
            "I0": row["uniform_cuboid_I0_kg_m2"],
            "I_full": row["supplement_full_apparatus_I_kg_m2"],
            "response_s2": observed["mean_delta_omega2_s-2"],
            "u_response_s2": observed["mean_delta_omega2_standard_u_s-2"],
            "magnetic_factor_minus_one": as_float(
                dec(MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]])
                * D("1e-6")),
            "mechanical_correction_central_ppm":
                MECHANICAL_CORRECTION_CENTRAL_PPM[row["id"]],
            "mechanical_correction_standard_u_ppm":
                MECHANICAL_CORRECTION_STANDARD_U_PPM[row["id"]],
        })

    figure = forward["ToS_figure_level_response"]
    aaf_figure = forward["AAF_figure_level_acceleration_stream"]
    return {
        "AAF": aaf,
        "TOS": tos,
        "TOS_figure_repeated_F1": {
            "A_B_A_response_s2": figure["A_B_A_background_subtracted_delta_omega2_s-2"],
            "quadratic_response_s2": figure["common_quadratic_background_subtracted_delta_omega2_s-2"],
        },
        "AAF_figure_unbound": {
            "source_amplitude_nrad_s2": aaf_figure["two_tone_source_amplitude_nrad_s-2"],
            "sample_count": aaf_figure["caption_scored_two_hour_sample_count"],
        },
    }


def uncertainty_packet(g_value, response, u_response, kernel, u_kernel,
                       kernel_axis_box, mechanical_factor,
                       u_mechanical_factor):
    response_component = abs(g_value * u_response / response)
    kernel_component = abs(g_value * u_kernel / kernel)
    mechanical_component = abs(
        g_value * u_mechanical_factor / mechanical_factor)
    rss = (response_component * response_component
           + kernel_component * kernel_component
           + mechanical_component * mechanical_component).sqrt()
    return {
        "response_standard_component_SI": as_float(response_component),
        "kernel_public_input_RSS_component_SI": as_float(kernel_component),
        "mechanical_correction_standard_component_SI":
            as_float(mechanical_component),
        "zero_covariance_partial_RSS_SI": as_float(rss),
        "linear_sum_of_reported_standard_components_SI": as_float(
            response_component + kernel_component + mechanical_component),
        "kernel_axis_box_diagnostic_component_SI": as_float(
            abs(g_value * kernel_axis_box / kernel)),
        "ceiling": (
            "First-order propagation of released response standard uncertainty, "
            "nominal-kernel public-input sensitivity, and Supplementary-Table-1 "
            "magnetic-correction standard uncertainty. The RSS assumes zero "
            "cross covariance. No uncertainty is available for the full "
            "normalized-kernel remainder; for ToS the missing signed fibre "
            "correction is also unquantified. None of these numbers is a "
            "coverage interval."
        ),
    }


def compute_primary(inputs):
    """Compute before any processed-coefficient or derived-G key is selected."""
    aaf_results = []
    for row in inputs["AAF"]:
        kernel = dec(row["K_hom"])
        response = dec(row["response_nrad_s2"]) * D("1e-9")
        u_response = dec(row["u_response_nrad_s2"]) * D("1e-9")
        factor = dec(row["mechanical_factor"])
        g_value = response * factor / kernel
        kernel_mix = dec(row["K_core_over_full_I"])
        collision = response * factor / kernel_mix
        aaf_results.append({
            "id": row["id"],
            "primary_conditional_G_SI": as_float(g_value),
            "primary_equation": "G=alpha_t*f_m/(K_hom+r_norm), evaluated at r_norm=0",
            "r_norm_premise": (
                "The real apparatus normalized m=2 source coefficient equals "
                "the pairwise-centred homogeneous coefficient. This is a "
                "conditional model premise, not a public-data result."
            ),
            "response_alpha_nrad_s-2": row["response_nrad_s2"],
            "response_class": (
                "Published campaign-average angular-acceleration response, already "
                "corrected for the air-density effect; not raw encoder data. It is "
                "a response numerator distinct from the processed source coefficient."
            ),
            "nominal_homogeneous_kernel_kg_m-3": row["K_hom"],
            "mechanical_factor_held_at_displayed_value": row["mechanical_factor"],
            "mechanical_correction_central_ppm":
                row["mechanical_correction_central_ppm"],
            "mechanical_correction_standard_u_ppm":
                row["mechanical_correction_standard_u_ppm"],
            "local_uncertainty_diagnostics": uncertainty_packet(
                g_value, response, u_response, kernel, dec(row["u_K_rss"]),
                dec(row["K_axis_box"]), factor,
                dec(row["mechanical_correction_standard_u_ppm"]) * D("1e-6")),
            "normalization_collision": {
                "core_numerator_over_full_I_kernel_kg_m-3": row["K_core_over_full_I"],
                "conditional_G_if_missing_mass_has_zero_m2_torque_SI": as_float(collision),
                "collision_minus_primary_relative_ppm": relative_ppm(collision, g_value),
                "meaning": (
                    "This deliberately incomplete mixed-normalization diagnostic "
                    "keeps the measured full inertia while omitting the unreported "
                    "mass from the m=2 numerator. Its different quotient shows the "
                    "numerical consequence of that mismatch; it is not a physical "
                    "apparatus coefficient or a second realized apparatus map."
                ),
            },
        })

    tos_results = []
    for row in inputs["TOS"]:
        kernel = dec(row["K_hom"])
        response = dec(row["response_s2"])
        u_response = dec(row["u_response_s2"])
        magnetic_factor = D(1) + dec(row["magnetic_factor_minus_one"])
        g_anchor = response * magnetic_factor / kernel
        slope = response / kernel
        kernel_mix = dec(row["K_core_over_full_I"])
        collision = response * magnetic_factor / kernel_mix
        tos_results.append({
            "id": row["id"],
            "forward_id": row["forward_id"],
            "primary_magnetic_only_G_anchor_SI": as_float(g_anchor),
            "primary_affine_family": {
                "equation": "G(c_f)=Delta_omega2/K_hom*(1+delta_m+c_f), at r_norm=0",
                "c_f_definition": "c_f=-Delta_K/(I*Delta_omega2)",
                "intercept_SI": as_float(g_anchor),
                "slope_dG_dc_f_SI": as_float(slope),
                "public_identified_interval": None,
                "reason": (
                    "The public fields pinned by the forward lane do not supply "
                    "the signed row-level Delta_K correction or a deterministic "
                    "compact bound for this independently reconstructed family."
                ),
            },
            "r_norm_premise": (
                "The real apparatus normalized near-minus-far stiffness kernel "
                "equals the homogeneous coefficient; this is conditional."
            ),
            "response_Delta_omega2_s-2": row["response_s2"],
            "response_class": (
                "Published mean stiffness response. It is a processed response "
                "summary, and the Supplement states that the source-gravity "
                "nonlinearity was corrected synchronously in Delta_omega2. It is "
                "therefore not a source-model-free raw numerator."
            ),
            "nominal_homogeneous_kernel_kg_m-3": row["K_hom"],
            "magnetic_factor_minus_one_held_at_displayed_value": row["magnetic_factor_minus_one"],
            "mechanical_correction_central_ppm":
                row["mechanical_correction_central_ppm"],
            "mechanical_correction_standard_u_ppm":
                row["mechanical_correction_standard_u_ppm"],
            "local_uncertainty_diagnostics_at_c_f_zero": uncertainty_packet(
                g_anchor, response, u_response, kernel, dec(row["u_K_rss"]),
                dec(row["K_axis_box"]), magnetic_factor,
                dec(row["mechanical_correction_standard_u_ppm"]) * D("1e-6")),
            "normalization_collision_at_c_f_zero": {
                "core_numerator_over_full_I_kernel_kg_m-3": row["K_core_over_full_I"],
                "conditional_G_SI": as_float(collision),
                "collision_minus_primary_relative_ppm": relative_ppm(collision, g_anchor),
                "meaning": (
                    "This deliberately incomplete mixed-normalization diagnostic "
                    "uses the homogeneous-core numerator with the full inertia. "
                    "It is not a physical apparatus coefficient or a second "
                    "realized apparatus map."
                ),
            },
        })

    repeated = next(row for row in tos_results if row["id"] == "TOS-I-F1-repeat")
    repeated_input = next(row for row in inputs["TOS"]
                          if row["id"] == "TOS-I-F1-repeat")
    kernel = dec(repeated_input["K_hom"])
    factor = D(1) + dec(repeated_input["magnetic_factor_minus_one"])
    aba = dec(inputs["TOS_figure_repeated_F1"]["A_B_A_response_s2"])
    quadratic = dec(inputs["TOS_figure_repeated_F1"]["quadratic_response_s2"])
    figure_result = {
        "row_binding": repeated["id"],
        "A_B_A_magnetic_only_G_anchor_SI": as_float(aba * factor / kernel),
        "quadratic_drift_diagnostic_G_anchor_SI": as_float(quadratic * factor / kernel),
        "table_response_magnetic_only_G_anchor_SI": repeated["primary_magnetic_only_G_anchor_SI"],
        "ceiling": (
            "The A-B-A value is extracted from overlapping figure-level three-day "
            "period summaries and has no independent coverage assignment. The "
            "quadratic value is a drift diagnostic. Both retain c_f and r_norm "
            "unknown and are not extra G measurements."
        ),
    }

    aaf_figure_result = {
        "two_hour_source_harmonic_nrad_s-2": inputs["AAF_figure_unbound"]["source_amplitude_nrad_s2"],
        "sample_count": inputs["AAF_figure_unbound"]["sample_count"],
        "campaign_kernel_binding": None,
        "conditional_G_withheld": True,
        "reason": (
            "This is the strongest released-file-level response numerator, but "
            "the public segment is representative rather than a campaign-average "
            "alpha_t and is not lawfully bound to one AAF kernel row."
        ),
    }

    return {
        "schema": "WAC_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001",
        "status": (
            "AAF_RNORM0_CONDITIONAL_QUOTIENTS_COMPUTED__TOS_RNORM0_CF0_ANCHORS_"
            "AND_AFFINE_CF_FAMILIES_COMPUTED__NO_PUBLIC_INDEPENDENT_DETERMINISTIC_"
            "G_INTERVAL__"
            "MATCHED_FULL_NUMERATOR_AND_TOS_DELTAK_REMAIN_MISSING__NO_NEW_G"
        ),
        "primary_calculation_order": [
            "validate upstream hashes",
            "extract only nominal kernels, public-input sensitivities, response summaries, and displayed mechanical factors",
            "compute AAF r_norm=0 conditional quotients",
            "compute ToS r_norm=0,c_f=0 anchors and affine c_f families",
            "compute mixed-normalization diagnostics",
            "only then open processed coefficient and derived-G comparators",
        ],
        "model": {
            "AAF": "alpha_t=G*(K_hom+r_norm)/f_m",
            "TOS": "Delta_omega2=G*(K_hom+r_norm)/(1+delta_m+c_f)",
            "full_apparatus_identifiability": (
                "r_norm is not numerically identified from unprocessed public "
                "mass fields; c_f is additionally unavailable for ToS."
            ),
        },
        "AAF": aaf_results,
        "TOS": tos_results,
        "TOS_figure_level_released_response_diagnostic": figure_result,
        "AAF_figure_level_unbound_response": aaf_figure_result,
        "numerator_custody": {
            "released_file_level_noncircularity": (
                "The local replay never uses accepted G or a processed source "
                "coefficient to extract the figure-level response harmonics/contrasts. "
                "This is computational nonuse, not statistical or model independence."
            ),
            "campaign_response_ceiling": (
                "Table alpha_t and Delta_omega2 are corrected response summaries, "
                "not raw observations. In particular Delta_omega2 includes a "
                "source-gravity nonlinearity correction. The quotients are therefore "
                "conditional cross-checks, not fully independent measurements."
            ),
        },
        "claim_ceiling": (
            "These are conditional model quotients and affine anchors, not a new "
            "measurement or an independently identified deterministic compact "
            "interval for G from this family, not GC16, and not RGRL or Gravity "
            "Formation Theory confirmation. The authors' processed-model G "
            "summaries remain post-calculation comparators."
        ),
    }


def attach_post_comparators(primary, nominal, forward):
    """Select authors' processed-coefficient and derived-G keys after primary."""
    result = deepcopy(primary)
    aaf_nominal = {row["id"]: row for row in nominal["AAF"]}
    aaf_forward = {row["id"]: row
                   for row in forward["AAF_three_processed_coefficient_forwards"]}
    for row in result["AAF"]:
        source = aaf_nominal[row["id"]]["post_calculation_processed_comparator"]
        observed = aaf_forward[row["id"]]
        g_primary = dec(row["primary_conditional_G_SI"])
        g_published = dec(observed["published_G_summary_SI_comparison_only"])
        k_hom = dec(row["nominal_homogeneous_kernel_kg_m-3"])
        k_processed = dec(source["processed_coefficient_kg_m-3"])
        response = dec(row["response_alpha_nrad_s-2"]) * D("1e-9")
        factor = dec(row["mechanical_factor_held_at_displayed_value"])
        g_processed_displayed = response * factor / k_processed
        row["post_calculation_comparator"] = {
            "processed_kernel_kg_m-3": as_float(k_processed),
            "processed_minus_nominal_kernel_kg_m-3": as_float(k_processed - k_hom),
            "processed_minus_nominal_relative_ppm": relative_ppm(k_processed, k_hom),
            "authors_derived_G_SI": as_float(g_published),
            "primary_minus_authors_relative_ppm": relative_ppm(g_primary, g_published),
            "processed_coefficient_displayed_factor_forward_G_SI":
                as_float(g_processed_displayed),
            "upstream_component_reconstructed_processed_forward_G_SI":
                observed["recomputed_G_SI"],
            "ratio_identity_residual": as_float(
                g_primary / g_processed_displayed - k_processed / k_hom),
            "role": (
                "POST_CALCULATION_ONLY. The processed scalar locates the authors' "
                "model-minus-homogeneous remainder but is not an independent "
                "mass-map input and is not used in the primary quotient. The "
                "local displayed-factor forward is distinct from the upstream "
                "component-reconstructed factor, whose extra digits are retained "
                "only as a comparator."
            ),
        }

    tos_nominal = {row["id"]: row for row in nominal["TOS"]}
    tos_forward = {row["id"]: row
                   for row in forward["ToS_seven_processed_coefficient_forwards"]}
    for row in result["TOS"]:
        source = tos_nominal[row["id"]]["post_calculation_processed_comparator"]
        observed = tos_forward[row["forward_id"]]
        response = dec(row["response_Delta_omega2_s-2"])
        k_hom = dec(row["nominal_homogeneous_kernel_kg_m-3"])
        k_processed = dec(source["processed_coefficient_kg_m-3"])
        magnetic_factor = D(1) + dec(row["magnetic_factor_minus_one_held_at_displayed_value"])
        g_published = dec(observed["published_G_summary_SI_comparison_only"])
        g_anchor = dec(row["primary_magnetic_only_G_anchor_SI"])
        c_total_hom = g_published * k_hom / response - magnetic_factor
        c_dynamic_processed = g_published * k_processed / response - magnetic_factor
        row["post_calculation_comparator"] = {
            "processed_kernel_kg_m-3": as_float(k_processed),
            "processed_minus_nominal_kernel_kg_m-3": as_float(k_processed - k_hom),
            "processed_minus_nominal_relative_ppm": relative_ppm(k_processed, k_hom),
            "authors_derived_G_SI": as_float(g_published),
            "primary_anchor_minus_authors_relative_ppm": relative_ppm(g_anchor, g_published),
            "c_total_required_with_homogeneous_kernel_ppm": as_float(c_total_hom * D(1000000)),
            "c_dynamic_required_with_processed_kernel_ppm": as_float(c_dynamic_processed * D(1000000)),
            "meaning": (
                "The first inferred correction combines the independent-kernel "
                "gap with the missing ToS correction. The second uses the authors' "
                "processed kernel and reproduces the already known unowned dynamic "
                "bracket. Both are post-comparison identities, not inputs."
            ),
        }
    return result


def calculate():
    validate_pins()
    nominal = json.loads(NOMINAL_RESULT.read_text(encoding="utf-8"))
    forward = json.loads(FORWARD_RESULT.read_text(encoding="utf-8"))
    inputs = extract_primary_inputs(nominal, forward)
    primary = compute_primary(inputs)
    return attach_post_comparators(primary, nominal, forward)


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2, sort_keys=True))
