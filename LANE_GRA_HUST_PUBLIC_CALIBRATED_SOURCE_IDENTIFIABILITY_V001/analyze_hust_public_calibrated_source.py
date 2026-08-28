#!/usr/bin/env python3
"""Build the HUST public calibrated partial source model without accepted G.

The primary calculation combines the independently reconstructed homogeneous
Newtonian kernels with only the central physical corrections publicly disclosed
in Li et al. (2018).  Authors' processed kernels and derived G values are attached
only after every primary quotient has been fixed, and only as comparators.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIELDS = HERE / "CALIBRATION_FIELDS.json"
SOURCE_CUSTODY = HERE / "SOURCE_CUSTODY.json"
RESULT = HERE / "RESULT.json"

NOMINAL_DIR = ROOT / "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001"
CONDITIONAL_DIR = ROOT / "LANE_GRA_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001"

EXPECTED_HASHES = {
    HERE / "SOURCE/HUST_2018_main_article_public_mirror.pdf": "40756ec0fb8f00c1fde31020b294521a3b220a196bef884a2ea5f3534d77dfaa",
    HERE / "SOURCE/41586_2018_431_MOESM1_ESM.pdf": "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb",
    HERE / "SOURCE/nature_main_table1_error_budget.html": "23436d4be7600a7a9dffa02cc4167a20b6eea032a181e77899bb57bb90aa02e9",
    NOMINAL_DIR / "THEOREM.md": "cd729640537d66c52da0c9209fb94c1a95ff1a9dc7580ac4b6a9e4a7cea8e67c",
    NOMINAL_DIR / "RESULT.json": "1dc31ef45bbd5d8441867b2575e2dec9d5290ee6f1a57f71781f80e0730cd875",
    NOMINAL_DIR / "SOURCE_FIELDS.json": "8540b9162df4363f64f892bbc2a56cfa542bb43767f605844bc9aa94da5571d8",
    NOMINAL_DIR / "reconstruct_hust_nominal_source_kernels.py": "4abac89235677b1757d66539e20fe21d3b8f9a5468796f976f701dc9ab074d59",
    CONDITIONAL_DIR / "THEOREM.md": "3ddcdcc8d4ef9f905c9ff3e07e813efc2848317e0f2cde4141798b9143c0e3a8",
    CONDITIONAL_DIR / "RESULT.json": "84a7d78a3a618345f89aa93639563077ee2de8043d46e2c91763a8604c7751eb",
    CONDITIONAL_DIR / "calculate_conditional_homogeneous_g.py": "d47578841b208ebeaf4f677528fdd03ecc3330349c8b1f597838b0971e2eea51",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_dependencies() -> dict[str, str]:
    observed: dict[str, str] = {}
    for path, expected in EXPECTED_HASHES.items():
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"dependency hash mismatch for {path}: {actual} != {expected}")
        observed[str(path.relative_to(ROOT))] = actual
    return observed


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def ppm_factor(value: float) -> float:
    return 1.0 + value * 1e-6


def rss(values: Iterable[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def outer(vector: list[float]) -> list[list[float]]:
    return [[left * right for right in vector] for left in vector]


def zeros(rows: int, columns: int) -> list[list[float]]:
    return [[0.0 for _ in range(columns)] for _ in range(rows)]


def add_outer(matrix: list[list[float]], vector: list[float]) -> None:
    for row in range(len(vector)):
        for column in range(len(vector)):
            matrix[row][column] += vector[row] * vector[column]


def quadratic(weights: list[float], matrix: list[list[float]]) -> float:
    return sum(
        weights[row] * matrix[row][column] * weights[column]
        for row in range(len(weights))
        for column in range(len(weights))
    )


def normalize_inverse_variance(uncertainties: list[float]) -> list[float]:
    raw = [1.0 / (value * value) for value in uncertainties]
    total = sum(raw)
    return [value / total for value in raw]


def matmul_left_right(
    left: list[list[float]], matrix: list[list[float]]
) -> list[list[float]]:
    # Return left * matrix * left^T.
    rows = len(left)
    inner = len(matrix)
    result = zeros(rows, rows)
    for i in range(rows):
        for j in range(rows):
            result[i][j] = sum(
                left[i][a] * matrix[a][b] * left[j][b]
                for a in range(inner)
                for b in range(inner)
            )
    return result


def aggregate_material_densities(geometry: dict) -> dict:
    pendulums = {}
    for apparatus, body in geometry["pendulums"].items():
        volume = body["length_m"] * body["width_m"] * body["height_m"]
        pendulums[apparatus] = {
            "uniform_core_volume_m3": volume,
            "uniform_core_average_density_kg_m-3": body["mass_kg"] / volume,
            "ceiling": "Average density of the nominal rectangular body, not the measured planar density map.",
        }

    spheres = {}
    for apparatus, bodies in geometry["source_spheres"].items():
        rows = []
        total_mass = 0.0
        total_volume = 0.0
        for body in bodies:
            volume = 4.0 * math.pi * (body["diameter_m"] / 2.0) ** 3 / 3.0
            total_mass += body["mass_kg"]
            total_volume += volume
            rows.append({
                "id": body["id"],
                "sphere_volume_m3": volume,
                "average_density_kg_m-3": body["mass_kg"] / volume,
            })
        spheres[apparatus] = {
            "bodies": rows,
            "aggregate_mass_kg": total_mass,
            "aggregate_volume_m3": total_volume,
            "aggregate_average_density_kg_m-3": total_mass / total_volume,
            "ceiling": "Mass/nominal-sphere-volume density; it does not resolve density inhomogeneity or nonsphericity.",
        }

    aaf_density = spheres["AAF"]["aggregate_average_density_kg_m-3"]
    return {
        "pendulum_nominal_core": pendulums,
        "source_spheres": spheres,
        "AAF_air_check_at_article_approximate_rho_air": {
            "rho_air_kg_m-3": 1.18,
            "rho_sphere_aggregate_kg_m-3": aaf_density,
            "rho_air_over_rho_sphere_ppm": 1.18 / aaf_density * 1e6,
            "article_approximate_density_ratio_ppm": 1.18 / 7965.0 * 1e6,
            "ceiling": "The exact campaign corrections use monitored run-level air density and are 149.90, 147.33 and 148.27 ppm; 1.18 kg/m^3 is explicitly approximate.",
        },
    }


def covariance_model(fields: dict) -> dict:
    inputs = fields["covariance_inputs"]

    aaf_systematics = list(inputs["AAF_systematic_standard_u_ppm"].values())
    aaf_stat = inputs["AAF_statistical_standard_u_ppm"]
    aaf_cov = zeros(3, 3)
    for index, value in enumerate(aaf_stat):
        aaf_cov[index][index] = value * value
    for vector in aaf_systematics:
        add_outer(aaf_cov, vector)
    aaf_implied_totals = [math.sqrt(aaf_cov[i][i]) for i in range(3)]
    aaf_weights = normalize_inverse_variance(inputs["AAF_total_standard_u_ppm"])
    aaf_combined = math.sqrt(quadratic(aaf_weights, aaf_cov))

    tos_systematics = list(inputs["TOS_systematic_standard_u_ppm"].values())
    tos_totals = inputs["TOS_total_standard_u_ppm"]
    tos_stat = inputs["TOS_statistical_standard_u_ppm"]
    tos_cov = zeros(7, 7)
    # Cross-fibre systematic covariance follows the published same-item 100%
    # correlation rule; diagonals are the published rounded total variances.
    for vector in tos_systematics:
        add_outer(tos_cov, vector)
    for index, total in enumerate(tos_totals):
        tos_cov[index][index] = total * total

    groups = [[0, 1], [2], [3, 4], [5, 6]]
    combined_uncertainties = inputs["TOS_combined_fibre_standard_u_ppm"]
    combination = zeros(4, 7)
    inferred_shared_background = {}
    for group_index, indices in enumerate(groups):
        local_uncertainties = [tos_totals[index] for index in indices]
        local_weights = normalize_inverse_variance(local_uncertainties)
        for weight, index in zip(local_weights, indices):
            combination[group_index][index] = weight
        if len(indices) == 2:
            i, j = indices
            wi, wj = local_weights
            target = combined_uncertainties[group_index] ** 2
            covariance = (
                target - wi * wi * tos_totals[i] ** 2 - wj * wj * tos_totals[j] ** 2
            ) / (2.0 * wi * wj)
            tos_cov[i][j] = covariance
            tos_cov[j][i] = covariance
            systematic_covariance = sum(vector[i] * vector[j] for vector in tos_systematics)
            inferred_shared_background[str(group_index)] = {
                "row_indices": indices,
                "total_covariance_ppm2": covariance,
                "systematic_covariance_ppm2": systematic_covariance,
                "shared_background_statistical_covariance_ppm2": covariance - systematic_covariance,
                "shared_background_standard_u_ppm": math.sqrt(max(0.0, covariance - systematic_covariance)),
                "ceiling": "Inferred from rounded row totals and rounded combined-fibre uncertainty; not a released component measurement.",
            }
    tos_fibre_cov = matmul_left_right(combination, tos_cov)
    tos_weights = normalize_inverse_variance(combined_uncertainties)
    tos_combined = math.sqrt(quadratic(tos_weights, tos_fibre_cov))

    return {
        "AAF": {
            "order": inputs["AAF_order"],
            "covariance_ppm2": aaf_cov,
            "implied_row_standard_u_ppm": aaf_implied_totals,
            "published_row_standard_u_ppm": inputs["AAF_total_standard_u_ppm"],
            "inverse_published_variance_weights": aaf_weights,
            "implied_combined_standard_u_ppm": aaf_combined,
            "published_combined_standard_u_ppm": inputs["AAF_reported_combined_standard_u_ppm"],
        },
        "TOS": {
            "run_order": inputs["TOS_order"],
            "run_covariance_ppm2": tos_cov,
            "fibre_combination_matrix": combination,
            "inferred_same_fibre_shared_background": inferred_shared_background,
            "fibre_order": ["F1", "F2", "F3", "F4"],
            "fibre_covariance_ppm2": tos_fibre_cov,
            "inverse_published_variance_weights": tos_weights,
            "implied_combined_standard_u_ppm": tos_combined,
            "published_combined_standard_u_ppm": 11.64,
        },
        "ceiling": fields["correlation_semantics"]["ceiling"],
    }


def build_result() -> dict:
    dependency_hashes = verify_dependencies()
    fields = load_json(FIELDS)
    load_json(SOURCE_CUSTODY)  # syntax and presence check
    nominal = load_json(NOMINAL_DIR / "RESULT.json")
    geometry = load_json(NOMINAL_DIR / "SOURCE_FIELDS.json")
    conditional = load_json(CONDITIONAL_DIR / "RESULT.json")

    nominal_aaf = {row["id"]: row for row in nominal["AAF"]}
    nominal_tos = {row["id"]: row for row in nominal["TOS"]}
    conditional_aaf = {row["id"]: row for row in conditional["AAF"]}
    conditional_tos = {row["id"]: row for row in conditional["TOS"]}

    # Primary calculations are completed without opening any authors-derived G
    # or processed-kernel comparator field.
    aaf_primary: list[dict] = []
    for campaign_id, calibration in fields["AAF"].items():
        source = nominal_aaf[campaign_id]
        response = conditional_aaf[campaign_id]
        mass_terms = calibration["mass_model_G_corrections_ppm"]
        mass_central = sum(value[0] for value in mass_terms.values())
        mass_u = rss(value[1] for value in mass_terms.values())
        homogeneous_kernel = source["nominal_homogeneous_coefficient_kg_m-3"]
        calibrated_kernel = homogeneous_kernel / ppm_factor(mass_central)

        air = calibration["air_density_G_correction_ppm"][0]
        averaging = calibration["data_average_G_correction_ppm"][0]
        derivative = calibration["numeric_derivative_G_correction_ppm"][0]
        transfer_factor = ppm_factor(air) * ppm_factor(averaging) * ppm_factor(derivative)
        alpha_published = response["response_alpha_nrad_s-2"]
        alpha_deprocessed = alpha_published / transfer_factor
        sensor_kernel = calibrated_kernel / transfer_factor
        mechanical_factor = response["mechanical_factor_held_at_displayed_value"]
        primary_g = alpha_published * 1e-9 * mechanical_factor / calibrated_kernel
        deprocessed_g = alpha_deprocessed * 1e-9 * mechanical_factor / sensor_kernel

        omega = calibration["omega_d_rad_s"]
        sinc_average = math.sin(omega * 0.5) / (omega * 0.5)
        sinc_derivative = math.sin(omega * 10.0) / (omega * 10.0)
        formula_average = (1.0 / sinc_average - 1.0) * 1e6
        formula_derivative = (1.0 / (sinc_derivative * sinc_derivative) - 1.0) * 1e6

        aaf_primary.append({
            "id": campaign_id,
            "nominal_homogeneous_kernel_kg_m-3": homogeneous_kernel,
            "public_mass_model_G_correction_sum_ppm": mass_central,
            "public_mass_model_correction_standard_u_RSS_ppm": mass_u,
            "public_calibrated_partial_kernel_kg_m-3": calibrated_kernel,
            "normalization": "K_partial=K_hom/(1+sum_j c_mass,j), because the published c_mass,j are corrections to G for omitted physical pendulum bodies/maps.",
            "mechanical_damper_factor": mechanical_factor,
            "published_air_corrected_alpha_nrad_s-2": alpha_published,
            "algebraically_deprocessed_alpha_nrad_s-2": alpha_deprocessed,
            "air_and_acquisition_sensor_kernel_kg_m-3": sensor_kernel,
            "transfer_factor_from_displayed_corrections": transfer_factor,
            "acquisition_formula_check_ppm": {
                "data_average_calculated": formula_average,
                "data_average_published": averaging,
                "numeric_derivative_calculated": formula_derivative,
                "numeric_derivative_published": derivative,
            },
            "primary_partial_G_SI": primary_g,
            "deprocessed_convention_same_G_SI": deprocessed_g,
            "identity_relative_residual": deprocessed_g / primary_g - 1.0,
            "identified_family": {
                "equation": "G_i(r_i)=alpha_i*f_m,i/(K_partial,i+r_i)",
                "r_i_units": "kg m^-3",
                "public_packet_owned_compact_domain": None,
                "reason": "The authors' processed coefficient makes one comparator remainder numerically inferable, but the public packet does not independently reconstruct it from the physical maps or specify an independently owned deterministic admissible set.",
            },
            "numerator_ceiling": "The campaign alpha is a published processed response summary. The displayed air, averaging and derivative corrections can be algebraically inverted, but this does not recreate raw encoder samples or an independent campaign fit.",
        })

    tos_primary: list[dict] = []
    for run_id, calibration in fields["TOS"].items():
        source = nominal_tos[run_id]
        response = conditional_tos[run_id]
        mass_terms = calibration["mass_model_G_corrections_ppm"]
        mass_central = sum(value[0] for value in mass_terms.values())
        mass_u = rss(value[1] for value in mass_terms.values())
        homogeneous_kernel = source["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"]
        calibrated_kernel = homogeneous_kernel / ppm_factor(mass_central)

        thermal = calibration["thermoelastic_G_correction_ppm"][0]
        gravitational_nonlinearity = calibration["gravitational_nonlinearity_G_correction_ppm"][0]
        response_transfer_factor = ppm_factor(thermal) * ppm_factor(gravitational_nonlinearity)
        delta_published = response["response_Delta_omega2_s-2"]
        delta_deprocessed = delta_published / response_transfer_factor
        anelastic = calibration["anelastic_G_correction_ppm"][0]
        magnetic = calibration["magnetic_damper_G_correction_ppm"][0]
        dynamic_factor = ppm_factor(anelastic + magnetic)
        primary_g = delta_published * dynamic_factor / calibrated_kernel
        deprocessed_g = delta_deprocessed * response_transfer_factor * dynamic_factor / calibrated_kernel

        tos_primary.append({
            "id": run_id,
            "nominal_homogeneous_kernel_kg_m-3": homogeneous_kernel,
            "public_mass_model_G_correction_sum_ppm": mass_central,
            "public_mass_model_correction_standard_u_RSS_ppm": mass_u,
            "public_calibrated_partial_kernel_kg_m-3": calibrated_kernel,
            "normalization": "K_partial=K_hom/(1+sum_j c_mass,j), with identical numerator/inertia ownership.",
            "published_corrected_Delta_omega2_s-2": delta_published,
            "algebraically_deprocessed_Delta_omega2_s-2": delta_deprocessed,
            "response_reapplication_factor": response_transfer_factor,
            "signed_anelastic_G_correction_ppm": anelastic,
            "magnetic_damper_G_correction_ppm": magnetic,
            "dynamic_factor": dynamic_factor,
            "primary_partial_G_SI": primary_g,
            "deprocessed_convention_same_G_SI": deprocessed_g,
            "identity_relative_residual": deprocessed_g / primary_g - 1.0,
            "identified_family": {
                "equation": "G_i(r_i)=Deltaomega2_i*(1+c_anel,i+c_mag,i)/(K_partial,i+r_i)",
                "r_i_units": "kg m^-3",
                "public_packet_owned_compact_domain": None,
                "reason": "The recovered signed anelastic correction removes c_f as a free central parameter. The authors' processed coefficient makes one comparator remainder numerically inferable, but no independently reconstructed physical-harmonic remainder or deterministic admissible set is public.",
            },
            "numerator_ceiling": "Published Delta omega squared was synchronously corrected for thermal and source-gravity nonlinearity. The disclosed central corrections can be inverted algebraically, but the source-gravity correction was itself calculated with the authors' apparatus model and no campaign raw time series is public.",
        })

    # Comparator attachment occurs only after all primary rows are fixed.
    for row in aaf_primary:
        comparator = conditional_aaf[row["id"]]["post_calculation_comparator"]
        processed_kernel = comparator["processed_kernel_kg_m-3"]
        calibrated_kernel = row["public_calibrated_partial_kernel_kg_m-3"]
        authors_g = comparator["authors_derived_G_SI"]
        original_kernel = row["nominal_homogeneous_kernel_kg_m-3"]
        row["post_calculation_comparator"] = {
            "authors_processed_kernel_kg_m-3": processed_kernel,
            "authors_processed_minus_public_partial_kg_m-3": processed_kernel - calibrated_kernel,
            "authors_processed_minus_public_partial_relative_ppm": (processed_kernel / calibrated_kernel - 1.0) * 1e6,
            "original_processed_minus_homogeneous_relative_ppm": (processed_kernel / original_kernel - 1.0) * 1e6,
            "absolute_remainder_reduction_factor": abs((processed_kernel / original_kernel - 1.0) / (processed_kernel / calibrated_kernel - 1.0)),
            "authors_derived_G_SI": authors_g,
            "primary_minus_authors_relative_ppm": (row["primary_partial_G_SI"] / authors_g - 1.0) * 1e6,
            "role": "QUARANTINED POST-CALCULATION ONLY. Neither value formed K_partial or the primary quotient.",
        }

    for row in tos_primary:
        comparator = conditional_tos[row["id"]]["post_calculation_comparator"]
        processed_kernel = comparator["processed_kernel_kg_m-3"]
        calibrated_kernel = row["public_calibrated_partial_kernel_kg_m-3"]
        authors_g = comparator["authors_derived_G_SI"]
        original_kernel = row["nominal_homogeneous_kernel_kg_m-3"]
        row["post_calculation_comparator"] = {
            "authors_processed_kernel_kg_m-3": processed_kernel,
            "authors_processed_minus_public_partial_kg_m-3": processed_kernel - calibrated_kernel,
            "authors_processed_minus_public_partial_relative_ppm": (processed_kernel / calibrated_kernel - 1.0) * 1e6,
            "original_processed_minus_homogeneous_relative_ppm": (processed_kernel / original_kernel - 1.0) * 1e6,
            "absolute_remainder_reduction_factor": abs((processed_kernel / original_kernel - 1.0) / (processed_kernel / calibrated_kernel - 1.0)),
            "authors_derived_G_SI": authors_g,
            "primary_minus_authors_relative_ppm": (row["primary_partial_G_SI"] / authors_g - 1.0) * 1e6,
            "role": "QUARANTINED POST-CALCULATION ONLY. Neither value formed K_partial or the primary quotient.",
        }

    return {
        "schema": "WAC_HUST_PUBLIC_CALIBRATED_SOURCE_IDENTIFIABILITY_RESULT_V001",
        "date": "2026-08-27",
        "status": "PASS__PUBLIC_CALIBRATED_PARTIAL_SOURCE_MODEL_AND_EXACT_IDENTIFIABILITY_CLASS",
        "dependency_hashes": dependency_hashes,
        "accepted_or_CODATA_G_numeric_inputs": [],
        "primary_calculation_order": "Independent homogeneous kernel -> public central physical mass corrections -> published response and disclosed transfer operators -> conditional quotient -> only then processed-kernel/authors-G comparator.",
        "material_density_inventory": aggregate_material_densities(geometry),
        "AAF": aaf_primary,
        "TOS": tos_primary,
        "published_category_covariance": covariance_model(fields),
        "minimal_independent_remainder_theorem": {
            "point_summary_level": "For each released row i, every missing spatial mass/stress field enters the measured harmonic only through one scalar normalized-kernel remainder r_i. One independently owned r_i per row is therefore sufficient for row-wise point evaluation, and at least one distinguishing scalar is necessary whenever two physically admissible undisclosed maps can have different harmonic kernels.",
            "row_remainders": [row["id"] for row in aaf_primary + tos_primary],
            "ownership_qualification": "Supplementary Tables 2 and 3 already publish the authors' processed kernels, so an authors-model comparator r_i is numerically inferable. What remains missing is an independently owned and reproducible physical-harmonic remainder; the ten row coordinates are not asserted to be ten independent physical degrees of freedom.",
            "physical_fields_sufficient_to_compute_them": [
                "measured pendulum density maps",
                "source-sphere density maps and disclosed run orientations",
                "individual residual 3D CMM coordinates, including the two AAF shear coordinates",
                "coating, clamp, ferrule, mirror, glue, edge-flaw and silica-rod maps",
                "AAF shelf, deformation, rotating-part and compensation-block mass maps"
            ],
            "independent_raw_reanalysis_addition": "A source-model-free reanalysis additionally needs campaign-bound raw response samples, the complete correction ledger, and the raw covariance/design matrix. Algebraic inversion of published corrections is not a replacement for those data.",
            "compact_interval_result": "The public packet owns no independently reconstructed deterministic or coverage-certified compact interval for r_i: quoted one-standard-deviation uncertainties and isolated upper bounds are not a joint deterministic domain. An authors-model conventional display band can be formed under additional assumptions; this is not a claim that physical r_i is mathematically unbounded.",
        },
        "published_bound_fragments": fields["published_bound_fragments"],
        "strict_ceilings": [
            "This is a calibrated public partial Newtonian source model, not a full finite-element reconstruction.",
            "The comparator residuals use authors-processed coefficients only after the primary calculation and are not independent evidence for the partial model.",
            "Published corrected response summaries are not raw numerator data; ToS gravitational-nonlinearity deprocessing remains authors-model-mediated.",
            "The ToS invert/reapply check is exact only under the declared multiplicative correction convention; it does not recover the unique historical raw-processing operator.",
            "Standard uncertainties and the reconstructed category covariance are not deterministic bounds or coverage theorems.",
            "No result here is evidence for RGRL, record lineage, beta_TM, gravity emergence, non-Newtonian gravity, or a common metric.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write RESULT.json")
    args = parser.parse_args()
    result = build_result()
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.write:
        RESULT.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
