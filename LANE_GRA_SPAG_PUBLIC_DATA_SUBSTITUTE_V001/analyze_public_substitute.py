#!/usr/bin/env python3
"""Reproduce the bounded SPAG Lane-A public-data substitute result.

The calculation verifies exact local source bytes, recomputes the public-packet
identifiability ranks, summarizes only native Fuchs noise observables, and
computes an explicitly optimistic NIST/BIPM planning envelope.  It does not
score a lineage effect and accepts no caller-selected data or thresholds.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import NormalDist, median
from typing import Iterable

from pypdf import PdfReader


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULT_PATH = HERE / "RESULT.json"

EXPECTED_SHA256 = {
    "GRAVITY_RGRL_SPAG_PROSPECTIVE_PROTOCOL_V001.md":
        "9495ca2b9edf3ebf1133e077d746e77b78ebda6e0fc061c178c80109506386b9",
    "GRAVITY_SPAG_EIR_CONSISTENCY_AND_EMPIRICAL_PATH_V001.md":
        "1bab3e5901a1566d287a2a7b04563718fbfab321b75ba41e48aa6cd152862ae1",
    "GRAVITY_RGRL_ONSHELL_OFFSHELL_CLARIFICATION_ADOPTION_V001.md":
        "4959f99898b216edc7da3e212ce2e26422287899fcf8f3b41cd34ef5d8bb3ff8",
    "PROGRAM_EXPERIMENT_REGISTER_V001.md":
        "084ed543ddaa5c55fd90e12a76c43688016aedd9218fe8850c34a5b30514c0d5",
    "GRAVITY_EMERGENCE_EXPERIMENT_REGISTER_V001.md":
        "ae7d2e672b3ba59f9c93d160c8562c95541b8e4c00d68fb21bf27d5315c2b58c",
    "LANE_GRA_I_FUCHS_RESPONSE/OUTPUT/normalized_response.csv":
        "b9b86396585cd7a53ce83677a5fab0628693ee475680bdf3f444011c11c2cc68",
    "LANE_GRA_I_FUCHS_RESPONSE/OUTPUT/analysis_summary.json":
        "e0d75137b4c46725f05b1cc87de4c313cc26341ce995c20fb482449aa40f7e83",
    "LANE_CROSS_RFT_MGFT_PAGE_GEILKER_BRANCH_GRAVITY_V001/FROZEN_DATA.json":
        "5cee159f21adeda25ea1df3ffeae95bdcbd9a569826971d0222b40776bb71f0a",
    "LANE_GRA_J_GRAVITY_HOLDOUT/zenodo_10995225_metadata.json":
        "8a7c1820c8af5c3799d2176309a2e1e6c6127ae4cf49158853878c1d90407f98",
    "LANE_GRA_J_GRAVITY_HOLDOUT/panda_2310.01344.pdf":
        "ff80f46b222d8e62bddeab0e5ad473fba5d7a372df89f99256397c46fa895ac9",
}

NIST_PDF_SHA256 = "c79552d62f4d4f4e85cfbbb00f135c1d985b596d9cdcde9bee57cfe4618f33dc"


def sha256(path: Path) -> str:
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"missing, non-file, or symlinked source: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_sources() -> dict[str, str]:
    observed: dict[str, str] = {}
    for relative, expected in EXPECTED_SHA256.items():
        actual = sha256(ROOT / relative)
        if actual != expected:
            raise RuntimeError(f"source hash mismatch: {relative}")
        observed[relative] = actual
    nist_relative = "SOURCE/nist_bipm_2026.pdf"
    actual = sha256(HERE / nist_relative)
    if actual != NIST_PDF_SHA256:
        raise RuntimeError("NIST/BIPM paper hash mismatch")
    observed[nist_relative] = actual
    return observed


def exact_rank(rows: Iterable[Iterable[int]]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    n_rows = len(matrix)
    n_cols = len(matrix[0])
    pivot_row = 0
    for column in range(n_cols):
        pivot = next(
            (row for row in range(pivot_row, n_rows) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(n_rows):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == n_rows:
            break
    return pivot_row


def page_geilker_result() -> dict:
    path = ROOT / (
        "LANE_CROSS_RFT_MGFT_PAGE_GEILKER_BRANCH_GRAVITY_V001/"
        "FROZEN_DATA.json"
    )
    packet = json.loads(path.read_text(encoding="utf-8"))
    x_code = packet["coding"]["decision"]
    m_code = packet["coding"]["mass_configuration"]
    x = [int(x_code[row["decision"]]) for row in packet["rows"]]
    m = [int(m_code[row["mass_configuration"]]) for row in packet["rows"]]
    if x != m or len(x) != 10:
        raise RuntimeError("Page--Geilker deterministic alias changed")

    rank_intercept_x_m = exact_rank([[1, xi, mi] for xi, mi in zip(x, m)])

    # This is deliberately the most generous possible historical recoding:
    # treat the decay decision as if it were T and set the wholly absent dummy
    # factor D to +1.  It still occupies only 2/8 cells and the saturated SPAG
    # design has rank 2/8, so beta_TM is not identifiable.
    spag_proxy_rows = []
    support = set()
    for ti, mi in zip(x, m):
        di = 1
        support.add((mi, ti, di))
        spag_proxy_rows.append(
            [1, mi, ti, di, ti * mi, di * mi, ti * di, ti * di * mi]
        )
    saturated_rank = exact_rank(spag_proxy_rows)
    if rank_intercept_x_m != 2 or saturated_rank != 2 or len(support) != 2:
        raise RuntimeError("Page--Geilker rank result changed")
    return {
        "rows": len(x),
        "observed_decision_mass_rule": "X_EQUALS_M_ON_EVERY_ROW",
        "rank_of_intercept_decision_mass": rank_intercept_x_m,
        "generous_proxy_support_cells_of_8": len(support),
        "generous_proxy_saturated_spag_rank_of_8": saturated_rank,
        "beta_TM_identifiable": False,
        "reason": (
            "decision and mass are deterministic aliases; no independently "
            "randomized target or dummy lineage redistribution exists"
        ),
    }


def fuchs_result() -> dict:
    path = ROOT / "LANE_GRA_I_FUCHS_RESPONSE/OUTPUT/normalized_response.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 24:
        raise RuntimeError("Fuchs trace count changed")
    if any(row["demodulator_physical_mapping"] != "UNKNOWN_IN_DEPOSIT" for row in rows):
        raise RuntimeError("Fuchs channel mapping unexpectedly changed")
    if any(
        name.lower() in {"l_t", "l_d", "target_lineage", "dummy_lineage"}
        for name in rows[0]
    ):
        raise RuntimeError("unexpected lineage column in Fuchs table")

    files = {row["source_filename"]: int(row["source_bytes"]) for row in rows}
    geometries = {(float(row["lateral_cm"]), float(row["height_cm"])) for row in rows}
    by_demod: dict[str, dict] = {}
    for demodulator in sorted({row["demodulator"] for row in rows}):
        selected = [row for row in rows if row["demodulator"] == demodulator]
        control_asd = [
            float(row["welch_control_median_v_per_sqrt_hz"]) for row in selected
        ]
        half_ratio = [
            float(row["half_amplitude_ratio_second_over_first"]) for row in selected
        ]
        by_demod[demodulator] = {
            "trace_count": len(selected),
            "native_residual_control_asd_v_per_sqrt_hz": {
                "minimum": min(control_asd),
                "median": median(control_asd),
                "maximum": max(control_asd),
            },
            "half_amplitude_ratio": {
                "minimum": min(half_ratio),
                "maximum": max(half_ratio),
            },
        }
    spans_h = [float(row["actual_span_s"]) / 3600.0 for row in rows]
    return {
        "files": len(files),
        "source_total_bytes": sum(files.values()),
        "geometries": len(geometries),
        "traces": len(rows),
        "trace_span_hours": {
            "minimum": min(spans_h),
            "median": median(spans_h),
            "maximum": max(spans_h),
        },
        "native_noise_diagnostics_by_unmapped_demodulator": by_demod,
        "lineage_factors_present": False,
        "beta_TM_identifiable": False,
        "force_power_computable": False,
        "reason": (
            "the deposit has no randomized L_T or L_D, no physical demodulator "
            "map, and no complete SI-transfer/covariance/null packet"
        ),
    }


def nist_result() -> dict:
    pdf_path = HERE / "SOURCE/nist_bipm_2026.pdf"
    reader = PdfReader(pdf_path)
    if len(reader.pages) != 31:
        raise RuntimeError("unexpected NIST/BIPM PDF page count")
    page = re.sub(r"\s+", " ", reader.pages[25].extract_text())
    table = {
        "Sapphire": {"free": 13.9799, "u_free": 0.0002, "servo": 13.9778, "u_servo": 0.0003},
        "Copper_0_deg": {"free": 31.1979, "u_free": 0.0003, "servo": 31.1962, "u_servo": 0.0004},
        "Copper_120_deg": {"free": 31.1842, "u_free": 0.0003, "servo": 31.1828, "u_servo": 0.0006},
        "Copper_240_deg": {"free": 31.1856, "u_free": 0.0002, "servo": 31.1836, "u_servo": 0.0003},
    }
    required_rows = (
        r"Sapphire 13\.9799 ±0\.0002 13\.9778 ±0\.0003",
        r"Copper0 ◦ 31\.1979 ±0\.0003 31\.1962 ±0\.0004",
        r"Copper120 ◦ 31\.1842 ±0\.0003 31\.1828 ±0\.0006",
        r"Copper240 ◦ 31\.1856 ±0\.0002 31\.1836 ±0\.0003",
    )
    for pattern in required_rows:
        if not re.search(pattern, page):
            raise RuntimeError(f"NIST Table 15 transcription check failed: {pattern}")

    u_values = [
        row[key]
        for row in table.values()
        for key in ("u_free", "u_servo")
    ]
    u_min = min(u_values)
    u_max = max(u_values)

    # Explicit design conventions for an optimistic single-primary-contrast
    # planning envelope.  This is not the full simultaneous Lane-A rule.
    familywise_alpha = 0.01
    target_power = 0.90
    critical = NormalDist().inv_cdf(1.0 - familywise_alpha / 2.0)
    power_quantile = NormalDist().inv_cdf(target_power)
    factor = critical + power_quantile
    # Table 15 reports a full mass-position torque difference, Delta N. If one
    # independent difference were available in each of the four (T,D) strata,
    # beta_TM=(1/8) sum_{T,D} T DeltaN_TD and u(beta_TM)=u(DeltaN)/4.
    sigma_beta_min = u_min / 4.0
    sigma_beta_max = u_max / 4.0
    detectable_min = factor * sigma_beta_min
    detectable_max = factor * sigma_beta_max
    torque_reference = table["Copper_0_deg"]["free"]

    copper_free = [table[name]["free"] for name in table if name.startswith("Copper")]
    copper_servo = [table[name]["servo"] for name in table if name.startswith("Copper")]
    method_offsets = [
        table[name]["free"] - table[name]["servo"]
        for name in table
        if name.startswith("Copper")
    ]
    return {
        "paper": "Schlamminger_et_al_Metrologia_63_025012_2026",
        "paper_pages": len(reader.pages),
        "table_15_units": "nN m",
        "table_15_type_A_k": 1,
        "table_15": table,
        "public_event_level_lineage_data_present": False,
        "lineage_factors_present": False,
        "beta_TM_identifiable": False,
        "ordinary_torque_reference_nN_m": torque_reference,
        "copper_orientation_spread_nN_m": {
            "free": max(copper_free) - min(copper_free),
            "servo": max(copper_servo) - min(copper_servo),
        },
        "free_minus_servo_offsets_nN_m": method_offsets,
        "optimistic_single_primary_contrast_planning_envelope": {
            "status": "PLANNING_LOWER_BOUND_NOT_A_PROSPECTIVE_DETECTION_LIMIT",
            "assumptions": [
                "four independent mass-position torque differences, one per T-by-D stratum",
                "equal independent Gaussian Type-A uncertainty per stratum difference",
                "NIST Table-15 difference uncertainty transfers unchanged to every new route stratum",
                "no multiplicity penalty beyond one two-sided primary contrast",
                "zero added route, lineage, covariance, drift, or systematic penalty",
            ],
            "two_sided_alpha": familywise_alpha,
            "target_power": target_power,
            "normal_critical_quantile": critical,
            "normal_power_quantile": power_quantile,
            "stratum_mass_difference_standard_uncertainty_range_nN_m": [u_min, u_max],
            "beta_TM_standard_uncertainty_range_nN_m": [sigma_beta_min, sigma_beta_max],
            "minimum_detectable_beta_TM_range_nN_m": [detectable_min, detectable_max],
            "minimum_detectable_fraction_of_31p1979": [
                detectable_min / torque_reference,
                detectable_max / torque_reference,
            ],
        },
    }


def panda_result() -> dict:
    path = ROOT / "LANE_GRA_J_GRAVITY_HOLDOUT/zenodo_10995225_metadata.json"
    packet = json.loads(path.read_text(encoding="utf-8"))
    files = {
        item["key"]: {"bytes": item["size"], "checksum": item["checksum"]}
        for item in packet["files"]
    }
    expected = {
        "Data Extended Fig 1.csv": "md5:1af179a514fae49a032d4d4b93fcf409",
        "Data Fig 3a.csv": "md5:df578555b01f75b3373d42d5336b8190",
        "Data Fig 3b.csv": "md5:0d5d04f094b5a7ff56e8554ab5805d4a",
        "Data Fig 3d.csv": "md5:d77b516441a8f3c12b930fd560172790",
    }
    if {key: value["checksum"] for key, value in files.items()} != expected:
        raise RuntimeError("Panda public inventory changed")
    holdout_dir = ROOT / "LANE_GRA_J_GRAVITY_HOLDOUT"
    response_files_present = any((holdout_dir / name).exists() for name in expected)
    if response_files_present:
        raise RuntimeError("Panda response-bearing holdout file unexpectedly present")
    return {
        "public_inventory": files,
        "response_bearing_files_opened_or_scored_by_this_lane": False,
        "holdout_preserved": True,
        "lineage_factors_present": False,
        "same_parent_torsion_join_present": False,
        "beta_TM_identifiable": False,
        "role": "DOWNSTREAM_ATOM_PROBE_COMPONENT_AND_PRESERVED_HOLDOUT_ONLY",
    }


def decision_rule() -> dict:
    return {
        "normative_lane": "ADOPTED_LANE_A_COMPLETE_SOURCE_MATCHED_DISCOVERY",
        "physical_prediction_under_A04": "BETA_TM_PHYS_EQUALS_ZERO",
        "retrospective_public_rule": {
            "required_support": "RANDOMIZED_M_BY_L_T_BY_L_D_EIGHT_CELL_SAME_PARENT_PACKET",
            "if_support_or_lineage_custody_absent": "PUBLIC_DATA_NO_LANE_A_SCORE",
            "forbidden_repairs": [
                "post_hoc_lineage labels",
                "cross-root pseudo-cells",
                "treating a missing factor as zero",
                "calling ordinary source motion lineage randomization",
            ],
        },
        "prospective_error_rule": {
            "verdict_family_FWER_max": 0.01,
            "primary_power_minimum": 0.90,
            "eta_q_and_all_equivalence_bands": (
                "COMMISSIONED_AND_FROZEN_BEFORE_FIRST_SCORED_RESPONSE"
            ),
            "confidence_set": "SIMULTANEOUS_C_TM",
            "ordinary_error_enlargement": (
                "I_TM_EQUALS_C_TM_MINKOWSKI_PLUS_MINUS_EPSILON_COLL_PLUS_EPSILON_GEOM"
            ),
            "optional_stopping": False,
        },
        "run_A": {
            "failed_upstream_gate": "NO_ANCESTRY_RESULT",
            "zero_outside_I_TM_with_all_controls_passed": (
                "RUN_A_ANCESTRY_CORRELATED_RESIDUAL_CANDIDATE__NO_CONFIRMATION"
            ),
            "I_TM_inside_plus_minus_eta_q": (
                "RUN_A_EXPLORATORY_DECLARED_APPARATUS_BOUND"
            ),
            "overlap": "INCONCLUSIVE",
        },
        "held_out_run_B": {
            "reproduced_nonzero_source_bucket_unresolved": (
                "REPRODUCED_ANCESTRY_CORRELATED_GRAVITY_RESIDUAL__SOURCE_BUCKET_UNRESOLVED"
            ),
            "complete_null": "BOUNDED_NULL__DECLARED_APPARATUS_COLUMN_ONLY",
            "ordinary_channel_owns_response": (
                "COLLATERAL_OWNED__NO_LINEAGE_GRAVITY_RESULT"
            ),
            "source_classification_required_before_named_mechanism_claim": True,
        },
    }


def main() -> None:
    sources = verify_sources()
    result = {
        "schema": "WAC_SPAG_PUBLIC_DATA_SUBSTITUTE_V001",
        "status": (
            "PUBLIC_COMPONENT_FEASIBILITY_AND_OPTIMISTIC_PLANNING_ENVELOPE_ONLY__"
            "NO_PUBLIC_LANE_A_ESTIMAND"
        ),
        "source_sha256": sources,
        "page_geilker": page_geilker_result(),
        "fuchs": fuchs_result(),
        "nist_bipm": nist_result(),
        "panda": panda_result(),
        "lane_A_decision_rule": decision_rule(),
        "cross_packet_theorem": {
            "any_admitted_packet_has_randomized_L_T_and_L_D": False,
            "any_admitted_packet_has_all_eight_same_parent_cells": False,
            "beta_TM_identifiable_from_admitted_packets": False,
            "cross_root_pooling_repairs_identifiability": False,
            "reason": (
                "different apparatuses, outcomes, units, and physical parents do not "
                "supply missing within-parent randomized factors"
            ),
        },
        "exact_claim_ceiling": {
            "public_data_can_test": [
                "ordinary gravity response and branch-following endpoints",
                "raw detector response/noise diagnostics in native units",
                "paper-level torque scale and best-case factorial planning algebra",
                "component feasibility for torsion and atom-probe stages",
            ],
            "public_data_cannot_test": [
                "a causal lineage redistribution effect",
                "a nonzero Lane-A beta_TM",
                "RGRL-C off-shell rank",
                "empirical Gravity Formation Theory confirmation",
                "a common-freefall ancestry response",
                "a numerical lineage source functional or microscopic origin of G",
            ],
            "empirical_lineage_claim": False,
            "gravity_formation_confirmation": False,
        },
    }
    RESULT_PATH.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("SPAG_PUBLIC_DATA_SUBSTITUTE: PASS")


if __name__ == "__main__":
    main()
