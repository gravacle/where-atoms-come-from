#!/usr/bin/env python3
"""Independent narrow verification of the scientific public-substitute result."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import NormalDist

from pypdf import PdfReader


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    require(
        result["status"].endswith("NO_PUBLIC_LANE_A_ESTIMAND"),
        "claim ceiling changed",
    )
    require(not result["exact_claim_ceiling"]["empirical_lineage_claim"], "lineage overclaim")
    require(
        not result["exact_claim_ceiling"]["gravity_formation_confirmation"],
        "GFT overclaim",
    )

    # Independent Page--Geilker support/rank calculation without importing the
    # builder. With T=M and D=1, every saturated row is one of only two vectors.
    page_path = ROOT / (
        "LANE_CROSS_RFT_MGFT_PAGE_GEILKER_BRANCH_GRAVITY_V001/"
        "FROZEN_DATA.json"
    )
    page = json.loads(page_path.read_text(encoding="utf-8"))
    x = [page["coding"]["decision"][row["decision"]] for row in page["rows"]]
    m = [page["coding"]["mass_configuration"][row["mass_configuration"]] for row in page["rows"]]
    require(x == m and set(x) == {-1, 1}, "Page--Geilker alias not reproduced")
    saturated = {
        (1, mi, ti, 1, ti * mi, mi, ti, ti * mi)
        for ti, mi in zip(x, m)
    }
    require(len(saturated) == 2, "Page--Geilker proxy support changed")
    require(result["page_geilker"]["generous_proxy_saturated_spag_rank_of_8"] == 2, "rank changed")
    require(not result["page_geilker"]["beta_TM_identifiable"], "Page beta overclaim")

    # Independent Fuchs inventory and extrema.
    fuchs_path = ROOT / "LANE_GRA_I_FUCHS_RESPONSE/OUTPUT/normalized_response.csv"
    with fuchs_path.open(newline="", encoding="utf-8") as stream:
        fuchs = list(csv.DictReader(stream))
    require(len(fuchs) == 24, "Fuchs trace count changed")
    require(len({row["source_filename"] for row in fuchs}) == 8, "Fuchs file count changed")
    require(
        all(row["demodulator_physical_mapping"] == "UNKNOWN_IN_DEPOSIT" for row in fuchs),
        "Fuchs demodulator map changed",
    )
    all_ratios = [float(row["half_amplitude_ratio_second_over_first"]) for row in fuchs]
    require(math.isclose(min(all_ratios), 0.32177772406854027, rel_tol=0, abs_tol=1e-15), "Fuchs min ratio")
    require(math.isclose(max(all_ratios), 2.7453120596913974, rel_tol=0, abs_tol=1e-15), "Fuchs max ratio")
    require(not result["fuchs"]["force_power_computable"], "Fuchs power overclaim")

    # Independently extract the NIST table and recompute the idealized envelope.
    pdf = HERE / "SOURCE/nist_bipm_2026.pdf"
    require(
        digest(pdf) == "c79552d62f4d4f4e85cfbbb00f135c1d985b596d9cdcde9bee57cfe4618f33dc",
        "NIST PDF custody",
    )
    reader = PdfReader(pdf)
    require(len(reader.pages) == 31, "NIST page count")
    text = re.sub(r"\s+", " ", reader.pages[25].extract_text())
    require("31.1979 ±0.0003 31.1962 ±0.0004" in text, "NIST copper row missing")
    z = NormalDist().inv_cdf(0.995) + NormalDist().inv_cdf(0.90)
    expected = [z * 0.0002 / 4.0, z * 0.0006 / 4.0]
    observed = result["nist_bipm"]["optimistic_single_primary_contrast_planning_envelope"][
        "minimum_detectable_beta_TM_range_nN_m"
    ]
    require(
        all(math.isclose(a, b, rel_tol=0, abs_tol=1e-18) for a, b in zip(expected, observed)),
        "NIST planning envelope changed",
    )
    require(
        result["nist_bipm"]["optimistic_single_primary_contrast_planning_envelope"]["status"]
        == "PLANNING_LOWER_BOUND_NOT_A_PROSPECTIVE_DETECTION_LIMIT",
        "NIST ceiling changed",
    )

    # The Panda response objects remain absent from this holdout directory.
    panda_names = (
        "Data Extended Fig 1.csv",
        "Data Fig 3a.csv",
        "Data Fig 3b.csv",
        "Data Fig 3d.csv",
    )
    holdout = ROOT / "LANE_GRA_J_GRAVITY_HOLDOUT"
    require(not any((holdout / name).exists() for name in panda_names), "Panda holdout opened")
    require(result["panda"]["holdout_preserved"], "Panda holdout flag changed")

    cross = result["cross_packet_theorem"]
    require(not cross["any_admitted_packet_has_randomized_L_T_and_L_D"], "invented factors")
    require(not cross["any_admitted_packet_has_all_eight_same_parent_cells"], "invented support")
    require(not cross["beta_TM_identifiable_from_admitted_packets"], "identifiability overclaim")
    require(not cross["cross_root_pooling_repairs_identifiability"], "illegal pool")

    rule = result["lane_A_decision_rule"]
    require(rule["physical_prediction_under_A04"] == "BETA_TM_PHYS_EQUALS_ZERO", "A04 changed")
    require(
        rule["retrospective_public_rule"]["if_support_or_lineage_custody_absent"]
        == "PUBLIC_DATA_NO_LANE_A_SCORE",
        "retrospective rule changed",
    )
    require(rule["prospective_error_rule"]["optional_stopping"] is False, "optional stopping enabled")

    report = (HERE / "PUBLIC_DATA_SUBSTITUTE.md").read_text(encoding="utf-8")
    retired = (
        "LOCAL_RGRL_C_REDISTRIBUTION_FORCE_COLUMN_PASS",
        "LOCAL_RGRL_C_REDISTRIBUTION_COMMON_FREEFALL_COLUMN_PASS",
    )
    require(not any(label in report for label in retired), "retired SPAG verdict reused")
    require("not a prospective\nSPAG limit" in report, "planning caveat missing")
    require("This is a proof of a data-design ceiling, not evidence that the physical effect\nis zero." in report, "null caveat missing")

    verification = (
        "SPAG_PUBLIC_DATA_SUBSTITUTE_INDEPENDENT_CHECK: PASS\n"
        "Page--Geilker generous proxy support/rank: 2 cells, rank 2/8\n"
        "Fuchs: 8 files, 24 traces, unknown demodulator mapping, no force power\n"
        "NIST/BIPM idealized MDE: 0.00019286904345467502 to "
        "0.000578607130364025 nN m; planning lower bound only\n"
        "Panda response holdout: preserved\n"
        "Empirical lineage/GFT claim: false\n"
    )
    (HERE / "VERIFICATION.txt").write_text(verification, encoding="utf-8")
    print(verification, end="")


if __name__ == "__main__":
    main()
