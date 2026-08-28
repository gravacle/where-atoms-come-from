#!/usr/bin/env python3
"""Independent hostile audit of the frozen HUST calibrated-source lane.

This executable does not import the builder or its verifier.  It uses an
independent primary-table transcription, parses official Nature Table 1
directly, recomputes every forward and covariance quantity, and separates a
public authors-model comparator from an independently owned source statistic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from html.parser import HTMLParser
from pathlib import Path
from typing import Any


AUDIT_DIR = Path(__file__).resolve().parent
LANE = AUDIT_DIR.parent
ROOT = LANE.parent
TRANSCRIPTION = AUDIT_DIR / "INDEPENDENT_SOURCE_TRANSCRIPTION.json"
AUDIT_RESULT = AUDIT_DIR / "AUDIT_RESULT.json"


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def close(a: float, b: float, tolerance: float = 2e-12) -> bool:
    return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)


def factor(ppm: float) -> float:
    return 1.0 + ppm * 1e-6


def inverse_variance_weights(values: list[float]) -> list[float]:
    raw = [1.0 / (value * value) for value in values]
    return [value / sum(raw) for value in raw]


def zero_matrix(n: int) -> list[list[float]]:
    return [[0.0 for _ in range(n)] for _ in range(n)]


def add_outer(matrix: list[list[float]], vector: list[float]) -> None:
    for i in range(len(vector)):
        for j in range(len(vector)):
            matrix[i][j] += vector[i] * vector[j]


def quadratic(weights: list[float], covariance: list[list[float]]) -> float:
    return sum(weights[i] * covariance[i][j] * weights[j] for i in range(len(weights)) for j in range(len(weights)))


def left_right(transform: list[list[float]], matrix: list[list[float]]) -> list[list[float]]:
    output = zero_matrix(len(transform))
    for i in range(len(transform)):
        for j in range(len(transform)):
            output[i][j] = sum(
                transform[i][a] * matrix[a][b] * transform[j][b]
                for a in range(len(matrix))
                for b in range(len(matrix))
            )
    return output


def positive_definite(matrix: list[list[float]], tolerance: float = 1e-12) -> bool:
    n = len(matrix)
    lower = zero_matrix(n)
    for i in range(n):
        for j in range(i + 1):
            remainder = matrix[i][j] - sum(lower[i][k] * lower[j][k] for k in range(j))
            if i == j:
                if remainder <= tolerance:
                    return False
                lower[i][j] = math.sqrt(remainder)
            else:
                lower[i][j] = remainder / lower[j][j]
    return True


class NatureTableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_table = False
        self.in_cell = False
        self.cell: list[str] = []
        self.row: list[str] = []
        self.rows: list[list[str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "table" and not self.in_table:
            self.in_table = True
        elif self.in_table and tag == "tr":
            self.row = []
        elif self.in_table and tag in ("td", "th"):
            self.in_cell = True
            self.cell = []

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.cell.append(data)

    def handle_endtag(self, tag: str) -> None:
        if self.in_table and tag in ("td", "th"):
            self.row.append(" ".join("".join(self.cell).split()))
            self.in_cell = False
        elif self.in_table and tag == "tr" and self.row:
            self.rows.append(self.row)
        elif self.in_table and tag == "table":
            self.in_table = False


def parse_pair(cell: str) -> tuple[float | None, float | None]:
    if not cell or cell in ("−", "-"):
        return None, None
    cell = cell.replace(",", "")
    if "[" in cell:
        first, second = cell.split("[", 1)
        return float(first.strip()), float(second.rstrip("]").strip())
    value = float(cell)
    return value, value


def official_table_vectors(path: Path) -> dict[str, Any]:
    parser = NatureTableParser()
    parser.feed(path.read_text(encoding="utf-8"))
    if len(parser.rows) != 29 or len(parser.rows[0]) != 8:
        raise RuntimeError("unexpected official Nature Table 1 shape")
    by_name = {row[0]: row[1:] for row in parser.rows[1:]}
    map_names = {
        "Dimensions": "pendulum_dimensions",
        "Attitude": "attitude",
        "Density inhomogeneity": "density_inhomogeneity",
        "Coating layer": "coating",
        "Clamp and ferrule": "clamp_ferrule",
        "Others": "others",
        "Masses": "source_masses",
        "Horizontal distance": "horizontal_distance",
        "Vertical distance": "vertical_distance",
        "Positions, alignment": "positions_alignment",
        "Fibre nonlinearity": "fibre_nonlinearity",
        "Fibre anelasticity": "fibre_anelasticity",
        "Thermal effect": "thermal",
        "Time base": "time_base",
        "Gravitational nonlinearity": "gravitational_nonlinearity",
        "Rotating gravity gradient": "rotating_gravity_gradient",
        "Shelf deformation": "shelf_deformation",
        "Magnetic damper": "magnetic_damper",
        "Air density": "air_density",
        "Magnetic field": "magnetic_field",
        "Electrostatic field": "electrostatic",
        "Angle encoder": "angle_encoder",
        "Residual twist angle": "residual_twist",
    }
    tos: dict[str, list[float]] = {}
    aaf: dict[str, list[float]] = {}
    for source_name, target_name in map_names.items():
        cells = by_name[source_name]
        tos_cells = [parse_pair(cell) for cell in cells[:4]]
        expanded = [
            tos_cells[0][0], tos_cells[0][1], tos_cells[1][0],
            tos_cells[2][0], tos_cells[2][1], tos_cells[3][0], tos_cells[3][1],
        ]
        if any(value is not None for value in expanded):
            tos[target_name] = [float(value) for value in expanded if value is not None]
        aaf_values = [parse_pair(cell)[0] for cell in cells[4:]]
        if any(value is not None for value in aaf_values):
            aaf[target_name] = [float(value) for value in aaf_values if value is not None]
    stat_cells = [parse_pair(cell) for cell in by_name["Statistical error of Δω2 or αt"][:4]]
    total_cells = [parse_pair(cell) for cell in by_name["Total"][:4]]
    tos_stat = [stat_cells[0][0], stat_cells[0][1], stat_cells[1][0], stat_cells[2][0], stat_cells[2][1], stat_cells[3][0], stat_cells[3][1]]
    tos_total = [total_cells[0][0], total_cells[0][1], total_cells[1][0], total_cells[2][0], total_cells[2][1], total_cells[3][0], total_cells[3][1]]
    combined = [parse_pair(cell)[0] for cell in by_name["Combined uncertainty"][:4]]
    return {
        "AAF_systematic": aaf,
        "AAF_statistical": [float(parse_pair(cell)[0]) for cell in by_name["Statistical error of Δω2 or αt"][4:]],
        "AAF_total": [float(parse_pair(cell)[0]) for cell in by_name["Total"][4:]],
        "TOS_systematic": tos,
        "TOS_statistical": [float(value) for value in tos_stat],
        "TOS_total": [float(value) for value in tos_total],
        "TOS_combined_fibre": [float(value) for value in combined],
        "table_row_count": len(parser.rows),
    }


def correction_rows(transcription: dict[str, Any], method: str) -> list[dict[str, float | str]]:
    order = transcription[f"{method}_order"]
    values = transcription[method]
    output = []
    for i, row_id in enumerate(order):
        output.append(
            {
                "id": row_id,
                "coating": values["coating_central_ppm"][i],
                "clamp": values["clamp_central_ppm"][i],
                "ferrule": values["ferrule_central_ppm"][i],
                "others": values["others_central_ppm"][i],
            }
        )
    return output


def covariance_from_official_table(table: dict[str, Any]) -> dict[str, Any]:
    aaf_cov = zero_matrix(3)
    for i, value in enumerate(table["AAF_statistical"]):
        aaf_cov[i][i] = value * value
    for vector in table["AAF_systematic"].values():
        add_outer(aaf_cov, vector)
    aaf_weights = inverse_variance_weights(table["AAF_total"])
    aaf_combined = math.sqrt(quadratic(aaf_weights, aaf_cov))

    tos_cov = zero_matrix(7)
    for vector in table["TOS_systematic"].values():
        add_outer(tos_cov, vector)
    for i, total in enumerate(table["TOS_total"]):
        tos_cov[i][i] = total * total
    groups = [[0, 1], [2], [3, 4], [5, 6]]
    transform = [[0.0 for _ in range(7)] for _ in range(4)]
    inferred: list[dict[str, Any]] = []
    for group_index, indices in enumerate(groups):
        local = [table["TOS_total"][i] for i in indices]
        weights = inverse_variance_weights(local)
        for weight, i in zip(weights, indices):
            transform[group_index][i] = weight
        if len(indices) == 2:
            i, j = indices
            wi, wj = weights
            target_variance = table["TOS_combined_fibre"][group_index] ** 2
            covariance = (
                target_variance - wi * wi * table["TOS_total"][i] ** 2
                - wj * wj * table["TOS_total"][j] ** 2
            ) / (2 * wi * wj)
            systematic = sum(vector[i] * vector[j] for vector in table["TOS_systematic"].values())
            tos_cov[i][j] = tos_cov[j][i] = covariance
            inferred.append(
                {
                    "indices": indices,
                    "total_covariance_ppm2": covariance,
                    "systematic_covariance_ppm2": systematic,
                    "inferred_shared_background_standard_u_ppm": math.sqrt(max(0.0, covariance - systematic)),
                }
            )
    fibre_cov = left_right(transform, tos_cov)
    fibre_weights = inverse_variance_weights(table["TOS_combined_fibre"])
    tos_combined = math.sqrt(quadratic(fibre_weights, fibre_cov))
    return {
        "AAF_covariance_ppm2": aaf_cov,
        "AAF_weights": aaf_weights,
        "AAF_combined_standard_u_ppm": aaf_combined,
        "TOS_run_covariance_ppm2": tos_cov,
        "TOS_fibre_transform": transform,
        "TOS_fibre_covariance_ppm2": fibre_cov,
        "TOS_weights": fibre_weights,
        "TOS_combined_standard_u_ppm": tos_combined,
        "TOS_inferred_shared_background": inferred,
        "AAF_positive_definite": positive_definite(aaf_cov),
        "TOS_run_positive_definite": positive_definite(tos_cov),
        "TOS_fibre_positive_definite": positive_definite(fibre_cov),
    }


def audit_result() -> dict[str, Any]:
    transcription = load(TRANSCRIPTION)
    reacquisition = load(AUDIT_DIR / "SOURCE_REACQUISITION.json")
    fields = load(LANE / "CALIBRATION_FIELDS.json")
    result = load(LANE / "RESULT.json")
    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    analyzer = (LANE / "analyze_hust_public_calibrated_source.py").read_text(encoding="utf-8")
    custody = load(LANE / "SOURCE_CUSTODY.json")
    nominal = load(ROOT / "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/RESULT.json")
    conditional = load(ROOT / "LANE_GRA_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001/RESULT.json")

    local_custody = []
    for source in custody["sources"]:
        path = LANE / source["lane_path"]
        local_custody.append(
            {
                "id": source["id"],
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "expected_sha256": source["sha256"],
                "pass": sha256(path) == source["sha256"],
            }
        )
    dependency_custody = []
    for name, entry in custody["sealed_dependencies"].items():
        path = (LANE / entry["path"]).resolve()
        dependency_custody.append(
            {"id": name, "sha256": sha256(path), "expected_sha256": entry["sha256"], "pass": sha256(path) == entry["sha256"]}
        )

    official_table = official_table_vectors(LANE / "SOURCE/nature_main_table1_error_budget.html")
    covariance = covariance_from_official_table(official_table)

    field_transcription_checks: list[dict[str, Any]] = []
    for method in ("AAF", "TOS"):
        order = transcription[f"{method}_order"]
        source = transcription[method]
        keys = ["coating", "clamp", "ferrule", "others"]
        for index, row_id in enumerate(order):
            stored = fields[method][row_id]["mass_model_G_corrections_ppm"]
            for key in keys:
                expected = [source[f"{key}_central_ppm"][index], source[f"{key}_standard_u_ppm"][index]]
                field_transcription_checks.append(
                    {"field": f"{method}/{row_id}/{key}", "expected": expected, "stored": stored[key], "pass": stored[key] == expected}
                )
    for index, row_id in enumerate(transcription["AAF_order"]):
        for label, field_key in (
            ("magnetic_damper", "mechanical_damper_G_correction_ppm"),
            ("air_density", "air_density_G_correction_ppm"),
            ("data_average", "data_average_G_correction_ppm"),
            ("numeric_derivative", "numeric_derivative_G_correction_ppm"),
        ):
            expected = [transcription["AAF"][f"{label}_central_ppm"][index], transcription["AAF"][f"{label}_standard_u_ppm"][index]]
            stored = fields["AAF"][row_id][field_key]
            field_transcription_checks.append({"field": f"AAF/{row_id}/{label}", "expected": expected, "stored": stored, "pass": stored == expected})
    for index, row_id in enumerate(transcription["TOS_order"]):
        for label, field_key in (
            ("anelastic", "anelastic_G_correction_ppm"),
            ("magnetic_damper", "magnetic_damper_G_correction_ppm"),
            ("thermoelastic", "thermoelastic_G_correction_ppm"),
            ("gravitational_nonlinearity", "gravitational_nonlinearity_G_correction_ppm"),
        ):
            expected = [transcription["TOS"][f"{label}_central_ppm"][index], transcription["TOS"][f"{label}_standard_u_ppm"][index]]
            stored = fields["TOS"][row_id][field_key]
            field_transcription_checks.append({"field": f"TOS/{row_id}/{label}", "expected": expected, "stored": stored, "pass": stored == expected})

    table_field_checks: list[dict[str, Any]] = []
    cinputs = fields["covariance_inputs"]
    for method in ("AAF", "TOS"):
        for key, vector in official_table[f"{method}_systematic"].items():
            stored = cinputs[f"{method}_systematic_standard_u_ppm"].get(key)
            table_field_checks.append({"field": f"{method}/{key}", "official": vector, "stored": stored, "pass": stored == vector})
        for suffix, official_key in (
            ("statistical_standard_u_ppm", f"{method}_statistical"),
            ("total_standard_u_ppm", f"{method}_total"),
        ):
            stored = cinputs[f"{method}_{suffix}"]
            vector = official_table[official_key]
            table_field_checks.append({"field": f"{method}/{suffix}", "official": vector, "stored": stored, "pass": stored == vector})
    table_field_checks.append(
        {
            "field": "TOS/combined_fibre_standard_u_ppm",
            "official": official_table["TOS_combined_fibre"],
            "stored": cinputs["TOS_combined_fibre_standard_u_ppm"],
            "pass": cinputs["TOS_combined_fibre_standard_u_ppm"] == official_table["TOS_combined_fibre"],
        }
    )

    nominal_aaf = {row["id"]: row for row in nominal["AAF"]}
    nominal_tos = {row["id"]: row for row in nominal["TOS"]}
    conditional_aaf = {row["id"]: row for row in conditional["AAF"]}
    conditional_tos = {row["id"]: row for row in conditional["TOS"]}
    reported_aaf = {row["id"]: row for row in result["AAF"]}
    reported_tos = {row["id"]: row for row in result["TOS"]}

    forwards: dict[str, list[dict[str, Any]]] = {"AAF": [], "TOS": []}
    for index, row_id in enumerate(transcription["AAF_order"]):
        components = correction_rows(transcription, "AAF")[index]
        mass_sum = sum(float(components[key]) for key in ("coating", "clamp", "ferrule", "others"))
        khom = nominal_aaf[row_id]["nominal_homogeneous_coefficient_kg_m-3"]
        kpartial = khom / factor(mass_sum)
        alpha = transcription["AAF"]["published_corrected_response"][index]
        magnetic = transcription["AAF"]["magnetic_damper_central_ppm"][index]
        primary_g = alpha * 1e-9 * factor(magnetic) / kpartial
        transfer_values = [
            transcription["AAF"]["air_density_central_ppm"][index],
            transcription["AAF"]["data_average_central_ppm"][index],
            transcription["AAF"]["numeric_derivative_central_ppm"][index],
        ]
        transfer = math.prod(factor(value) for value in transfer_values)
        deprocessed = (alpha / transfer) * 1e-9 * factor(magnetic) / (kpartial / transfer)
        omega = transcription["AAF"]["omega_d_rad_s"][index]
        average_sinc = math.sin(omega * 0.5) / (omega * 0.5)
        derivative_sinc = math.sin(omega * 10.0) / (omega * 10.0)
        average_formula_ppm = (1.0 / average_sinc - 1.0) * 1e6
        derivative_formula_ppm = (1.0 / (derivative_sinc * derivative_sinc) - 1.0) * 1e6
        processed = transcription["AAF"]["published_processed_kernel_kg_m-3"][index]
        gap = (processed / kpartial - 1.0) * 1e6
        reported = reported_aaf[row_id]
        forwards["AAF"].append(
            {
                "id": row_id,
                "mass_correction_sum_ppm": mass_sum,
                "K_hom_kg_m-3": khom,
                "K_partial_kg_m-3": kpartial,
                "primary_partial_G_SI": primary_g,
                "transfer_factor": transfer,
                "data_average_formula_ppm": average_formula_ppm,
                "data_average_published_ppm": transfer_values[1],
                "numeric_derivative_formula_ppm": derivative_formula_ppm,
                "numeric_derivative_published_ppm": transfer_values[2],
                "source_formula_matches": abs(average_formula_ppm - transfer_values[1]) < 0.01 and abs(derivative_formula_ppm - transfer_values[2]) < 0.01,
                "deprocessed_identity_relative_residual": deprocessed / primary_g - 1.0,
                "authors_processed_minus_partial_ppm": gap,
                "authors_processed_remainder_kg_m-3": processed - kpartial,
                "reported_match": close(kpartial, reported["public_calibrated_partial_kernel_kg_m-3"]) and close(primary_g, reported["primary_partial_G_SI"], 2e-13) and close(gap, reported["post_calculation_comparator"]["authors_processed_minus_public_partial_relative_ppm"]),
                "source_response_match": close(alpha, conditional_aaf[row_id]["response_alpha_nrad_s-2"]),
            }
        )

    for index, row_id in enumerate(transcription["TOS_order"]):
        components = correction_rows(transcription, "TOS")[index]
        mass_sum = sum(float(components[key]) for key in ("coating", "clamp", "ferrule", "others"))
        khom = nominal_tos[row_id]["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"]
        kpartial = khom / factor(mass_sum)
        response = transcription["TOS"]["published_corrected_response"][index]
        anelastic = transcription["TOS"]["anelastic_central_ppm"][index]
        magnetic = transcription["TOS"]["magnetic_damper_central_ppm"][index]
        dynamic = factor(anelastic + magnetic)
        primary_g = response * dynamic / kpartial
        thermal = transcription["TOS"]["thermoelastic_central_ppm"][index]
        gravity_nonlin = transcription["TOS"]["gravitational_nonlinearity_central_ppm"][index]
        chosen_reapplication = factor(thermal) * factor(gravity_nonlin)
        additive_reapplication = factor(thermal + gravity_nonlin)
        deprocessed = (response / chosen_reapplication) * chosen_reapplication * dynamic / kpartial
        processed = transcription["TOS"]["published_processed_kernel_kg_m-3"][index]
        gap = (processed / kpartial - 1.0) * 1e6
        reported = reported_tos[row_id]
        forwards["TOS"].append(
            {
                "id": row_id,
                "mass_correction_sum_ppm": mass_sum,
                "K_hom_kg_m-3": khom,
                "K_partial_kg_m-3": kpartial,
                "signed_anelastic_ppm": anelastic,
                "magnetic_ppm": magnetic,
                "dynamic_total_ppm": anelastic + magnetic,
                "primary_partial_G_SI": primary_g,
                "chosen_multiplicative_reapplication_factor": chosen_reapplication,
                "additive_reapplication_factor": additive_reapplication,
                "multiplicative_minus_additive_ppm": (chosen_reapplication / additive_reapplication - 1.0) * 1e6,
                "deprocessed_identity_relative_residual": deprocessed / primary_g - 1.0,
                "authors_processed_minus_partial_ppm": gap,
                "authors_processed_remainder_kg_m-3": processed - kpartial,
                "reported_match": close(kpartial, reported["public_calibrated_partial_kernel_kg_m-3"]) and close(primary_g, reported["primary_partial_G_SI"], 2e-13) and close(gap, reported["post_calculation_comparator"]["authors_processed_minus_public_partial_relative_ppm"]),
                "source_response_match": close(response, conditional_tos[row_id]["response_Delta_omega2_s-2"]),
            }
        )

    reported_cov = result["published_category_covariance"]
    covariance_matches = {
        "AAF_matrix": all(close(covariance["AAF_covariance_ppm2"][i][j], reported_cov["AAF"]["covariance_ppm2"][i][j]) for i in range(3) for j in range(3)),
        "AAF_combined": close(covariance["AAF_combined_standard_u_ppm"], reported_cov["AAF"]["implied_combined_standard_u_ppm"]),
        "TOS_run_matrix": all(close(covariance["TOS_run_covariance_ppm2"][i][j], reported_cov["TOS"]["run_covariance_ppm2"][i][j]) for i in range(7) for j in range(7)),
        "TOS_fibre_matrix": all(close(covariance["TOS_fibre_covariance_ppm2"][i][j], reported_cov["TOS"]["fibre_covariance_ppm2"][i][j]) for i in range(4) for j in range(4)),
        "TOS_combined": close(covariance["TOS_combined_standard_u_ppm"], reported_cov["TOS"]["implied_combined_standard_u_ppm"]),
    }

    remainder_schema = result.get("minimal_independent_remainder_theorem", {})
    strict_ceilings = result.get("strict_ceilings", [])
    core_repair_checks = {
        "M1_theorem_numeric_comparator_inferable": "authors-model comparator \\(r_i\\) is numerically inferable" in theorem,
        "M1_theorem_independent_ownership_target": "independently owned \\(r_i\\) per released row" in theorem,
        "M1_theorem_joint_dimension_ceiling": "does not assert ten independent physical degrees" in theorem,
        "M1_result_schema_renamed": "minimal_independent_remainder_theorem" in result and "minimal_unreported_parameter_theorem" not in result,
        "M1_result_ownership_qualification": "numerically inferable" in remainder_schema.get("ownership_qualification", "") and "independently owned and reproducible" in remainder_schema.get("ownership_qualification", ""),
        "M1_result_joint_dimension_ceiling": "not asserted to be ten independent physical degrees" in remainder_schema.get("ownership_qualification", ""),
        "M1_analyzer_schema_renamed": '"minimal_independent_remainder_theorem"' in analyzer and "minimal_unreported_parameter_theorem" not in analyzer,
        "M1_analyzer_ownership_qualification": "numerically inferable" in analyzer and "independently owned and reproducible physical-harmonic remainder" in analyzer,
        "interval_theorem_authors_model_band": "An authors-model conventional display band can be formed under additional assumptions" in theorem,
        "interval_result_owned_ceiling": "no independently reconstructed deterministic or coverage-certified compact interval" in remainder_schema.get("compact_interval_result", ""),
        "N1_theorem_declared_product_convention": "under the declared multiplicative composition convention" in theorem,
        "N1_theorem_no_unique_raw_operator": "not proof that the convention uniquely reconstructs the historical raw-processing operator" in theorem,
        "N1_result_strict_ceiling": any("declared multiplicative correction convention" in item and "unique historical raw-processing operator" in item for item in strict_ceilings),
        "N1_analyzer_strict_ceiling": "declared multiplicative correction convention" in analyzer and "unique historical raw-processing operator" in analyzer,
    }

    all_pass = (
        all(row["pass"] for row in local_custody)
        and all(row["pass"] for row in dependency_custody)
        and all(obj["sealed_sha256_match"] for obj in reacquisition["objects"])
        and all(row["pass"] for row in field_transcription_checks)
        and all(row["pass"] for row in table_field_checks)
        and all(row["reported_match"] and row["source_response_match"] for method in forwards.values() for row in method)
        and all(row["source_formula_matches"] for row in forwards["AAF"])
        and all(covariance_matches.values())
        and all((row["signed_anelastic_ppm"] < 0 and row["magnetic_ppm"] > 0) for row in forwards["TOS"])
        and result["accepted_or_CODATA_G_numeric_inputs"] == []
        and all(core_repair_checks.values())
    )
    return {
        "audit_schema": "WAC_HUST_PUBLIC_CALIBRATED_SOURCE_INDEPENDENT_HOSTILE_AUDIT_RESULT_V001",
        "date": "2026-08-27",
        "disposition": "PASS__M1_N1_CORE_REPAIRS_INCORPORATED__PUBLICATION_SAFE_WITH_INTRINSIC_CEILINGS" if all_pass else "FAIL",
        "arithmetic_and_custody_all_pass": all_pass,
        "local_source_custody": local_custody,
        "sealed_dependency_custody": dependency_custody,
        "fresh_remote_reacquisition": reacquisition,
        "official_source_transcription_checks": {
            "manual_correction_fields": field_transcription_checks,
            "official_html_error_budget_fields": table_field_checks,
            "all_pass": all(row["pass"] for row in field_transcription_checks + table_field_checks),
        },
        "independent_forwards": forwards,
        "independent_covariance": covariance,
        "covariance_matches_frozen_result": covariance_matches,
        "accepted_or_CODATA_G_audit": {
            "numeric_input_used_in_independent_reconstruction": False,
            "frozen_result_declared_inputs": result["accepted_or_CODATA_G_numeric_inputs"],
            "authors_G_and_processed_kernel_role": "available only as post-calculation comparator; processed kernels were used above solely to calculate comparator gaps after K_partial and primary quotients were fixed",
        },
        "core_repair_audit": {
            "checks": core_repair_checks,
            "all_pass": all(core_repair_checks.values()),
            "M1_disposition": "CLOSED_IN_CORE_LANGUAGE_AND_EXECUTABLE_SCHEMA",
            "N1_disposition": "CLOSED_IN_CORE_LANGUAGE_AND_EXECUTABLE_SCHEMA",
            "publication_disposition": "The repaired builder THEOREM.md, RESULT.json, and analyzer are publication-safe unchanged within their explicit intrinsic ceilings; this audit is independent evidence, not a semantic repair crutch.",
        },
        "identifiability_audit": {
            "fixed_harmonic_forward": "Y_i = G (K_partial_i + r_i)",
            "sufficiency": "Given the released processed response central value Y_i, one independently owned scalar r_i per row is sufficient to evaluate the row quotient.",
            "necessity": "If two admissible unreported physical maps agree on all released independent fields but yield different harmonic kernels, no public-field-only point rule can distinguish them; some statistic that distinguishes their scalar harmonic values is necessary.",
            "critical_ownership_qualification": "The numerical authors-model value of this scalar is already inferable from the public processed kernel minus K_partial. What is missing is an independently owned/reproducible physical-harmonic statistic, not the mere existence of a public number.",
            "joint_dimension_qualification": "Ten row values are sufficient coordinates, but the theorem does not prove ten independent physical degrees of freedom; shared maps may correlate or constrain the row remainders.",
            "raw_reanalysis_qualification": "One r_i is not sufficient for a source-model-free raw-response reanalysis, which additionally needs raw samples, correction/event custody, design matrix, and covariance.",
        },
        "standard_uncertainty_domain_audit": {
            "authors_model_compact_display_bands_possible": True,
            "how": "A published processed coefficient plus/minus a quoted standard uncertainty is a finite conventional authors-model display band, and can be mapped monotonically to a finite quotient band if one accepts that model and a distribution/coverage convention.",
            "independent_deterministic_admissible_domain_owned": False,
            "why_not": "Standard uncertainties are not hard support bounds; the public release gives no joint distribution/covariance for the missing physical maps and no independently reconstructed central r_i. The possible compact display band is therefore neither a deterministic identified set nor an independent source cross-check.",
        },
        "deprocessing_audit": {
            "AAF": "The sinc attenuation corrections are source-owned and independently recomputed; dividing response and kernel by the same product is an exact quotient identity, not recovery of raw encoder data.",
            "TOS": "The source owns the signed correction values and says thermal and nonlinear corrections were applied synchronously, but it does not publish a raw correction operator proving the builder's particular product composition. Its invert/reapply identity is exact only by declared convention and is not raw-data deprocessing. The maximum product-versus-additive difference is reported in the row diagnostics and is not load-bearing to the primary quotient.",
        },
        "normalization_audit": {
            "mapping": "Because the table entries are corrections to G, G_corrected=G_hom(1+c_mass) is represented by K_partial=K_hom/(1+c_mass) at the quotient level. This preserves numerator/inertia ownership better than inserting omitted mass only into I.",
            "ceiling": "The correction scalars remain authors-model-mediated component calibrations transplanted onto the independently integrated homogeneous functional; K_partial is a hybrid calibrated partial model, not a raw-map reconstruction.",
        },
        "digits_audit": "Digits beyond the resolution and uncertainty of the public tables are deterministic computational replay digits, not independently owned physical significant figures.",
        "final_claim_ceiling": [
            "calibrated public partial Newtonian source family only",
            "one independently owned row-harmonic remainder is the minimal point summary only at fixed released response and fixed harmonic",
            "no source-model-free raw numerator, deterministic interval, joint coverage theorem, new G measurement, GR/RGRL/GFT/lineage/common-metric evidence",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = audit_result()
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if AUDIT_RESULT.read_text() != payload:
            raise SystemExit("AUDIT_RESULT.json differs from independent recomputation")
        print("PASS: AUDIT_RESULT.json matches independent recomputation")
        return 0
    if args.write:
        AUDIT_RESULT.write_text(payload)
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
