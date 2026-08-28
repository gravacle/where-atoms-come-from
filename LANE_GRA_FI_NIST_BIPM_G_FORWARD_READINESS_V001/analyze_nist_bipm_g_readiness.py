#!/usr/bin/env python3
"""Reproduce the NIST/BIPM-to-G-forward-model readiness audit.

The sole empirical source is the pinned 2026 NIST/BIPM primary paper.  The
calculation uses published torque, geometry, uncertainty, and correlation
summaries only.  It never imports CODATA, the paper's consensus result, or an
accepted value of G as a generator, calibration, prior, scan center, or check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np
from pypdf import PdfReader


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PDF = ROOT / "LANE_GRA_SPAG_PUBLIC_DATA_SUBSTITUTE_V001/SOURCE/nist_bipm_2026.pdf"
READINESS = HERE / "READINESS.md"
RESULT = HERE / "RESULT.json"

EXPECTED = {
    "pdf": "c79552d62f4d4f4e85cfbbb00f135c1d985b596d9cdcde9bee57cfe4618f33dc",
    "gc_theorem": "cbf0733633ba93756b08dded7486a9be76beb572807693455c761bd36a8f0f5b",
    "gc_protocol": "6ec7d8f0ce9a184d25612107dbfc294dd22d124ebe859831401e9cc0c8e8b819",
    "gc_model": "6c17498d2d65f6420498ac559a97a2c3bbf49e110dd971da34b4c9c9bea2e4e4",
}

DEPENDENCIES = {
    "pdf": PDF,
    "gc_theorem": ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/THEOREM.md",
    "gc_protocol": ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/PROTOCOL.md",
    "gc_model": ROOT / "LANE_GRA_GC_CALIBRATED_FINITE_APPARATUS_G_CROSSCHECK_V001/finite_apparatus_g_model.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def compact(text: str) -> str:
    return re.sub(r"\s+", "", text.replace("\x00", "").replace("\x01", ""))


def check_source() -> tuple[PdfReader, dict[int, str], dict[str, str]]:
    observed = {name: sha256(path) for name, path in DEPENDENCIES.items()}
    if observed != EXPECTED:
        raise RuntimeError(f"dependency hash mismatch: {observed}")
    reader = PdfReader(PDF)
    if len(reader.pages) != 31:
        raise RuntimeError("primary PDF page count changed")
    page_numbers = (6, 7, 10, 11, 15, 18, 23, 25, 26, 27)
    pages = {
        number: compact(reader.pages[number - 1].extract_text() or "")
        for number in page_numbers
    }
    required = {
        6: ("Figure3.LayoutoftheBIPMexperiment", "separatedby37.67◦", "Γmax"),
        7: ("zerosecondderivativeanddoesnotaffectω0", "∆Nmax", "Gf"),
        10: ("Table1.Relativesensitivities", "Rs−5.4−5.4", "Idisk0.110"),
        11: ("Table2.Technicalparameters", "1150.156g", "11191.68g", "5015.60g", "213.9642mm", "0.262pFrad−1"),
        15: ("Table6.Uncertaintybudget", "Γmax=0.485099", "combined16.2"),
        18: ("Table10.Themassintegration", "0.4856410", "Pointmasses&nodisk0.457999−56919"),
        23: ("Table12.Comparisonofmeasuredandcalculatedtorques", "Table13.Measuredperiods", "Idisk=74403gcm2"),
        25: ("Table14.Endpointsaandb", "Servocopper0◦−1.323.714.2", "AllquantitiesinthetableareexpressedrelativetoG"),
        26: (
            "Table15.Torquedifferences",
            "∆N=G16Γmaxmsmt/Rs",
            "Sapphire13.9799±0.000213.9778±0.0003213.96800.485613",
            "Copper0◦31.1979±0.000331.1962±0.0004213.96420.485637",
            "Copper120◦31.1842±0.000331.1828±0.0006213.98170.485467",
            "Copper240◦31.1856±0.000231.1836±0.0003213.98220.485464",
            "Table16.FourindividualGdeterminations",
            "Co.servo6.673642×10−1123",
            "Co.free6.674021×10−1130",
            "Sa.servo6.672637×10−1138",
            "Sa.free6.673636×10−1194",
        ),
        27: (
            "Table17.Relativeuncertaintycontributions",
            "Combined23.230.337.593.9",
            "Table18.Summaryofrelativeuncertaintiesandcorrelations",
            "123.20.420.380.12",
            "230.30.230.23",
            "337.50.25",
            "493.9",
        ),
    }
    for page, tokens in required.items():
        for token in tokens:
            if token not in pages[page]:
                raise RuntimeError(f"custody token absent on PDF page {page}: {token}")
    return reader, pages, observed


def primary_fields() -> dict:
    return {
        "test_mass": {
            "mass_kg": 1.150156,
            "inner_radius_m": 0.002502,
            "outer_radius_m": 0.027481,
            "height_m": 0.054910,
            "mass_circle_radius_m": 0.1200315,
        },
        "source_mass_kg": {"Copper": 11.19168, "Sapphire": 5.01560},
        "disk_moment_of_inertia_kg_m2": 7.45e-3,
        "pendulum_period_s": 121.0,
        "capacitance_pF": {"C12": 5.6, "C13": 22.03, "C23": 23.4},
        "capacitance_gradient_pF_per_rad": {
            "k12": -0.006,
            "k13": 0.262,
            "k23": -0.266,
        },
        "table_1_relative_sensitivity_of_inferred_G": {
            "ms": {"free": -1.0, "servo": -1.0},
            "mt": {"free": -0.1, "servo": -1.0},
            "Rs": {"free": -5.4, "servo": -5.4},
            "Rt": {"free": 2.4, "servo": 4.4},
            "Idisk": {"free": 0.11, "servo": 0.0},
            "phi_t": {"free": 1.0, "servo": -1.0},
        },
        "table_15": [
            {"configuration": "Sapphire", "source_type": "Sapphire", "mode": "free", "torque_nN_m": 13.9799, "type_A_u_nN_m": 0.0002, "Rs_mm": 213.9680, "Gamma": 0.485613},
            {"configuration": "Sapphire", "source_type": "Sapphire", "mode": "servo", "torque_nN_m": 13.9778, "type_A_u_nN_m": 0.0003, "Rs_mm": 213.9680, "Gamma": 0.485613},
            {"configuration": "Copper_0_deg", "source_type": "Copper", "mode": "free", "torque_nN_m": 31.1979, "type_A_u_nN_m": 0.0003, "Rs_mm": 213.9642, "Gamma": 0.485637},
            {"configuration": "Copper_0_deg", "source_type": "Copper", "mode": "servo", "torque_nN_m": 31.1962, "type_A_u_nN_m": 0.0004, "Rs_mm": 213.9642, "Gamma": 0.485637},
            {"configuration": "Copper_120_deg", "source_type": "Copper", "mode": "free", "torque_nN_m": 31.1842, "type_A_u_nN_m": 0.0003, "Rs_mm": 213.9817, "Gamma": 0.485467},
            {"configuration": "Copper_120_deg", "source_type": "Copper", "mode": "servo", "torque_nN_m": 31.1828, "type_A_u_nN_m": 0.0006, "Rs_mm": 213.9817, "Gamma": 0.485467},
            {"configuration": "Copper_240_deg", "source_type": "Copper", "mode": "free", "torque_nN_m": 31.1856, "type_A_u_nN_m": 0.0002, "Rs_mm": 213.9822, "Gamma": 0.485464},
            {"configuration": "Copper_240_deg", "source_type": "Copper", "mode": "servo", "torque_nN_m": 31.1836, "type_A_u_nN_m": 0.0003, "Rs_mm": 213.9822, "Gamma": 0.485464},
        ],
        "four_derived_G_covariance_inputs": {
            "order": ["Copper_servo", "Copper_free", "Sapphire_servo", "Sapphire_free"],
            "table_16_value_SI": [6.673642e-11, 6.674021e-11, 6.672637e-11, 6.673636e-11],
            "table_16_displayed_relative_u_ppm_rounded": [23.0, 30.0, 38.0, 94.0],
            "table_17_18_combined_relative_u_ppm": [23.2, 30.3, 37.5, 93.9],
        },
        "table_18_correlation": [
            [1.0, 0.42, 0.38, 0.12],
            [0.42, 1.0, 0.23, 0.23],
            [0.38, 0.23, 1.0, 0.25],
            [0.12, 0.23, 0.25, 1.0],
        ],
    }


def table_17_categories() -> tuple[list[str], np.ndarray, np.ndarray]:
    names = [
        "Mass integration", "Autocollimator non-linearity", "Type A",
        "Capacitance measurement", "Voltages measurement", "Zero torque",
        "Pressure correction", "SM MCR", "Mass of one SM", "Mass of one TM",
        "MOI disk", "TM MCR", "Anelasticity", "Period measurement",
    ]
    rows = np.array([
        [16.2, 16.2, 16.2, 16.2],
        [10.7, 14.2, 9.3, 85.7],
        [7.6, 4.6, 19.5, 13.9],
        [6.0, 0.0, 6.0, 0.0],
        [6.0, 0.0, 6.0, 0.0],
        [5.4, 5.4, 12.1, 12.1],
        [1.6, 1.6, 21.5, 21.5],
        [1.0, 1.0, 1.0, 1.0],
        [0.2, 0.2, 0.4, 0.4],
        [0.1, 0.0, 0.1, 0.0],
        [0.0, 19.8, 0.0, 19.8],
        [0.0, 2.9, 0.0, 2.9],
        [0.0, 0.7, 0.0, 0.7],
        [0.0, 0.5, 0.0, 0.5],
    ])
    printed = np.array([23.2, 30.3, 37.5, 93.9])
    return names, rows, printed


def analyze() -> dict:
    _, _, observed_hashes = check_source()
    fields = primary_fields()
    mt = fields["test_mass"]["mass_kg"]
    source_masses = fields["source_mass_kg"]

    rows = []
    for item in fields["table_15"]:
        ms = source_masses[item["source_type"]]
        rs = item["Rs_mm"] * 1e-3
        source_column = 16.0 * item["Gamma"] * ms * mt / rs
        torque = item["torque_nN_m"] * 1e-9
        sigma = item["type_A_u_nN_m"] * 1e-9
        ratio = torque / source_column
        ratio_u = sigma / source_column
        rows.append({
            **item,
            "finite_contrast_source_column_kg2_per_m": source_column,
            "jacobian_d_torque_d_G_kg2_per_m": source_column,
            "summary_ratio_G_SI": ratio,
            "type_A_only_ratio_u_SI_denominator_fixed": ratio_u,
            "type_A_only_relative_u": ratio_u / ratio,
        })

    a = np.array([row["finite_contrast_source_column_kg2_per_m"] for row in rows])
    y = np.array([row["torque_nN_m"] for row in rows]) * 1e-9
    sigma_y = np.array([row["type_A_u_nN_m"] for row in rows]) * 1e-9

    # Formal diagonal Type-A diagnostic only; the public paper does not supply
    # the eight-row covariance required to turn this into an estimate.
    info_diag = float(np.sum((a / sigma_y) ** 2))
    g_diag = float(np.sum(a * y / sigma_y**2) / np.sum(a * a / sigma_y**2))
    g_diag_u = float(info_diag ** -0.5)

    configs = [item["configuration"] for item in rows]
    unique_configs = list(dict.fromkeys(configs))
    b_config = np.array([[float(name == target) for target in unique_configs] for name in configs])
    b_method = np.array([[float(item["mode"] == mode) for mode in ("free", "servo")] for item in rows])
    eye = np.eye(len(rows))

    four_g = fields["four_derived_G_covariance_inputs"]
    g4 = np.array(four_g["table_16_value_SI"])
    rel4 = np.array(four_g["table_17_18_combined_relative_u_ppm"]) * 1e-6
    corr = np.array(fields["table_18_correlation"])
    sd4 = g4 * rel4
    cov4 = np.diag(sd4) @ corr @ np.diag(sd4)
    inv4 = np.linalg.inv(cov4)
    ones4 = np.ones(4)
    gls_weights = inv4 @ ones4 / float(ones4 @ inv4 @ ones4)
    gls_value = float(gls_weights @ g4)
    gls_u = math.sqrt(1.0 / float(ones4 @ inv4 @ ones4))

    category_names, category_rows, printed_combined = table_17_categories()
    category_rss = np.sqrt(np.sum(category_rows**2, axis=0))

    result = {
        "schema": "WAC_NIST_BIPM_G_FORWARD_READINESS_V001",
        "status": "PUBLIC_SUMMARY_REDUCED_TORQUE_FORWARD_EXACT__FULL_GC16_REAL_APPARATUS_FIT_NOT_READY",
        "source": {
            "title": "Redetermination of the gravitational constant with the BIPM torsion balance at NIST",
            "paper": "Schlamminger_et_al_Metrologia_63_025012_2026",
            "pdf_pages": 31,
            "sha256": observed_hashes["pdf"],
            "dependency_sha256": observed_hashes,
        },
        "page_table_custody": {
            "geometry_and_torque": {"pdf_page": 6, "printed_page": 5, "locators": ["Figure 3", "equations (1)-(5)"]},
            "free_transfer_and_nominal_zero_source_stiffness": {"pdf_page": 7, "printed_page": 6, "locators": ["equations (6)-(10)", "text below equation (7)"]},
            "mode_sensitivities": {"pdf_page": 10, "printed_page": 9, "locators": ["Table 1", "equation (21)"]},
            "nominal_parameters": {"pdf_page": 11, "printed_page": 10, "locators": ["Table 2"]},
            "mass_integration_sensitivity": {"pdf_page": 15, "printed_page": 14, "locators": ["Table 6", "equations (26)-(33)"]},
            "full_mass_integration_comparators": {"pdf_page": 18, "printed_page": 17, "locators": ["Table 10", "sections 6.7-6.9"]},
            "background_and_inertia_checks": {"pdf_page": 23, "printed_page": 22, "locators": ["Tables 12-13", "sections 8.2-8.4"]},
            "autocollimator_summary": {"pdf_page": 25, "printed_page": 24, "locators": ["Table 14", "equation (62)"]},
            "torque_rows_and_derived_G": {"pdf_page": 26, "printed_page": 25, "locators": ["Tables 15-16"]},
            "uncertainty_and_correlation": {"pdf_page": 27, "printed_page": 26, "locators": ["Tables 17-18", "equation (63)"]},
        },
        "published_fields": fields,
        "reduced_forward": {
            "observation": "published calibrated torque difference Delta_N",
            "formula": "Delta_N_j = G * A_j + r_j; A_j = 16*Gamma_j*ms_j*mt/Rs_j",
            "GC_mapping": "A_j is the finite two-state contrast analogue of GC source column a; it is not silently relabeled as the infinitesimal trajectory derivative",
            "nominal_source_gravitational_stiffness": 0.0,
            "nominal_source_stiffness_basis": "paper states source-mass contribution has zero second derivative at torque extrema",
            "rows": rows,
        },
        "jacobian_and_identifiability": {
            "G_only_source_column_rank": int(np.linalg.matrix_rank(a[:, None])),
            "G_only_identifiable_if_calibration_and_remainders_fixed": bool(np.any(a != 0.0)),
            "G_plus_one_free_remainder_per_row": {
                "design_shape": [8, 9],
                "design_rank": int(np.linalg.matrix_rank(np.column_stack((a, eye)))),
                "G_column_in_nuisance_span": True,
                "G_identifiable": False,
            },
            "G_plus_one_common_remainder_per_configuration": {
                "design_shape": [8, 5],
                "design_rank": int(np.linalg.matrix_rank(np.column_stack((a, b_config)))),
                "nuisance_rank": int(np.linalg.matrix_rank(b_config)),
                "G_column_in_nuisance_span": bool(np.linalg.norm(a - b_config @ np.linalg.lstsq(b_config, a, rcond=None)[0]) < 1e-10),
                "G_identifiable": False,
            },
            "G_plus_two_common_method_offsets": {
                "design_shape": [8, 3],
                "design_rank": int(np.linalg.matrix_rank(np.column_stack((a, b_method)))),
                "nuisance_rank": int(np.linalg.matrix_rank(b_method)),
                "G_identifiable_under_this_extra_restriction_only": True,
            },
            "G_and_free_global_source_scale": {
                "jacobian_rank": 1,
                "parameter_count": 2,
                "identified_combination": "p=G*s",
                "G_separately_identifiable": False,
            },
        },
        "supported_diagnostics_not_real_estimates": {
            "table_15_diagonal_Type_A_only": {
                "explicit_extra_assumptions": [
                    "eight torque summaries independent",
                    "published denominators fixed exactly",
                    "all physical remainders zero",
                    "no calibration covariance",
                ],
                "formal_WLS_G_SI": g_diag,
                "formal_standard_u_SI": g_diag_u,
                "formal_relative_u": g_diag_u / g_diag,
                "status": "ALGEBRA_DIAGNOSTIC_ONLY__NOT_A_G_ESTIMATE",
            },
            "table_17_RSS_check_ppm": {
                "category_order": category_names,
                "recomputed": category_rss.tolist(),
                "printed": printed_combined.tolist(),
                "maximum_rounding_difference_ppm": float(np.max(np.abs(category_rss - printed_combined))),
            },
            "table_18_four_derived_G_covariance": {
                "correlation_eigenvalues": np.linalg.eigvalsh(corr).tolist(),
                "correlation_condition_number": float(np.linalg.cond(corr)),
                "correlation_rank": int(np.linalg.matrix_rank(corr)),
                "covariance_SI": cov4.tolist(),
                "formal_GLS_weights": gls_weights.tolist(),
                "formal_GLS_value_SI": gls_value,
                "formal_GLS_standard_u_SI": gls_u,
                "formal_GLS_relative_u": gls_u / gls_value,
                "status": "DERIVED_G_SUMMARY_COVARIANCE_DIAGNOSTIC_ONLY__EXCLUDES_DARK_UNCERTAINTY_AND_IS_NOT_GC16_RAW_COVARIANCE",
            },
        },
        "fields_blocking_full_real_G_fit": [
            "event-level free-deflection angle/time/frequency data with run and source-position labels",
            "event-level servo voltage, capacitance-gradient, residual-twist, controller, and timing data",
            "complete source/test/disk finite-element mass-coordinate and density files used to compute Gamma and its trajectory derivatives",
            "row-level source/test mass and geometry calibration records plus their joint covariance and global source-scale interval",
            "full mechanical transfer including total inertia, frequency-dependent torsion response, damping, gimbal/support modes, and any auxiliary-mode couplings",
            "raw autocollimator calibration curve/ensemble and row-level readout transfer covariance",
            "eight-row torque observation covariance; Table 18 is only a four-by-four covariance summary of already-derived G values",
            "signed row-level physical remainder columns/bounds for gas thermal torque, background torque, disk model, local gravity, source drive/support reaction, electromagnetic and controller effects",
            "complete conserved apparatus stress/source ledger and finite-dither correction to the nominal k_g=0 extremum",
            "predeclared raw-data nuisance design, held-out rows, likelihood domain, and covariance treatment required by the GC protocol",
        ],
        "accepted_G_or_CODATA_used": False,
        "paper_consensus_or_dark_uncertainty_model_used": False,
        "claim_ceiling": {
            "paper_summary_ratio_reproduced": True,
            "full_GC16_forward_fit_executed": False,
            "independent_real_G_estimate": False,
            "lineage_or_gravity_emergence_test": False,
        },
    }
    return result


def verification_lines(result: dict) -> list[str]:
    rows = result["reduced_forward"]["rows"]
    diag = result["supported_diagnostics_not_real_estimates"]
    ident = result["jacobian_and_identifiability"]
    readiness_text = READINESS.read_text(encoding="utf-8")
    rendered_result = json.dumps(result, indent=2, sort_keys=True) + "\n"
    required_ceiling_phrases = (
        "ALGEBRA_DIAGNOSTIC_ONLY__NOT_A_G_ESTIMATE",
        "Table 18 is only a four-by-",
        "does not provide an independent numerical `G` estimate",
        "No accepted `G`, paper",
    )
    forbidden_promotions = (
        "this packet measures G",
        "this packet confirms Gravity Formation Theory",
        "Table 15 is event-level raw data",
    )
    checks = [
        (result["source"]["sha256"] == EXPECTED["pdf"], "pinned primary PDF hash"),
        (result["source"]["pdf_pages"] == 31, "primary PDF page count"),
        (len(result["page_table_custody"]) == 10, "page/table custody ledger"),
        (len(rows) == 8, "Table 15 eight torque-mode rows"),
        (all(row["finite_contrast_source_column_kg2_per_m"] > 0 for row in rows), "nonzero finite-contrast source columns"),
        (max(abs(row["summary_ratio_G_SI"]) for row in rows) < 1e-9, "summary ratios finite and SI-scaled"),
        (ident["G_only_source_column_rank"] == 1, "G-only source-column rank one"),
        (not ident["G_plus_one_free_remainder_per_row"]["G_identifiable"], "free-row-remainder nonidentifiability"),
        (ident["G_plus_one_common_remainder_per_configuration"]["G_column_in_nuisance_span"], "configuration-remainder exact alias"),
        (ident["G_plus_two_common_method_offsets"]["design_rank"] == 3, "restricted two-method-offset design rank"),
        (ident["G_and_free_global_source_scale"]["jacobian_rank"] == 1, "G/source-scale product rank one"),
        (diag["table_17_RSS_check_ppm"]["maximum_rounding_difference_ppm"] < 0.07, "Table 17 RSS reproduction"),
        (diag["table_18_four_derived_G_covariance"]["correlation_rank"] == 4, "Table 18 correlation full rank"),
        (min(diag["table_18_four_derived_G_covariance"]["correlation_eigenvalues"]) > 0, "Table 18 correlation positive definite"),
        (not result["accepted_G_or_CODATA_used"], "no accepted G or CODATA input"),
        (
            "table_19" not in rendered_result.lower()
            and "6.6743e-11" not in rendered_result,
            "accepted-G/CODATA hierarchy absent from result graph",
        ),
        (not result["paper_consensus_or_dark_uncertainty_model_used"], "paper consensus/dark model excluded"),
        (len(result["fields_blocking_full_real_G_fit"]) == 10, "blocking-field ledger complete"),
        (not result["claim_ceiling"]["full_GC16_forward_fit_executed"], "no full GC16 fit overclaim"),
        (not result["claim_ceiling"]["independent_real_G_estimate"], "no independent G estimate overclaim"),
        (not result["claim_ceiling"]["lineage_or_gravity_emergence_test"], "no lineage/gravity-emergence overclaim"),
        (RESULT.read_text(encoding="utf-8") == rendered_result, "RESULT.json byte-exact reproduction"),
        (all(phrase in readiness_text for phrase in required_ceiling_phrases), "mandatory readiness ceilings present"),
        (not any(phrase.lower() in readiness_text.lower() for phrase in forbidden_promotions), "forbidden promotions absent"),
    ]
    lines = [f"{'PASS' if ok else 'FAIL'}  {label}" for ok, label in checks]
    lines.append(f"SUMMARY {sum(ok for ok, _ in checks)}/{len(checks)} checks passed")
    if not all(ok for ok, _ in checks):
        failed = [label for ok, label in checks if not ok]
        raise RuntimeError(f"verification checks failed: {failed}")
    return lines


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = analyze()
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print("\n".join(verification_lines(result)))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"FAIL  {exc}", file=sys.stderr)
        sys.exit(1)
