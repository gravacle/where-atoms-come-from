#!/usr/bin/env python3
"""Verify the bounded SPAG public-data second-pass packet.

This verifier checks frozen documentary custody and claim ceilings. External
source bytes are not redistributed by this lane, so their observed hashes are
custody metadata rather than locally recomputed source seals.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    custody = json.loads((HERE / "SEARCH_CUSTODY.json").read_text(encoding="utf-8"))

    require(result["search"]["query_count"] == len(custody["queries"]) == 28, "query freeze")
    require(len(set(custody["queries"])) == 28, "duplicate query")
    require(all(isinstance(query, str) and query.strip() for query in custody["queries"]), "blank query")
    require(len(custody["candidates"]) == 2, "candidate count")
    require(result["search"]["new_same_parent_lineage_candidates"] == 0, "lineage candidate overclaim")
    search_custody = custody["search_execution_custody"]
    require(search_custody["query_strings_frozen"], "query strings not frozen")
    require(not search_custody["ordered_returned_hit_lists_frozen"], "invented hit-list custody")
    require(not search_custody["ranking_snapshots_frozen"], "invented ranking custody")
    require(not result["search"]["ordered_returned_hit_lists_frozen"], "result hit-list overclaim")
    require(not result["search"]["full_screened_result_set_byte_reproducible"], "search reproducibility overclaim")
    require("NOT_EXHAUSTIVE" in result["search"]["negative_result_scope"], "search scope ceiling")

    spag = result["spag"]
    require(not spag["any_new_candidate_has_native_M_L_T_L_D"], "invented lineage")
    require(not spag["any_new_candidate_has_all_eight_same_parent_cells"], "invented support")
    require(not spag["any_new_candidate_identifies_beta_TM"], "beta overclaim")
    require(not spag["cross_root_pooling_authorized"], "cross-root pseudo-cells")
    require(not spag["physical_null_inferred"], "absence called null")

    by_id = {candidate["id"]: candidate for candidate in custody["candidates"]}
    clock = by_id["C1_ZHENG_CLOCK_NETWORK"]
    require(clock["data_doi"] == "10.5281/zenodo.8184043", "clock root")
    require([f["bytes"] for f in clock["files"]] == [25550, 2876, 958, 369], "clock sizes")
    require(
        [f["remote_md5"] for f in clock["files"]]
        == [
            "a829e2044f0ef2dd450435d7b790e8c7",
            "8edb5d957e81a07b07e9d350afe43e3a",
            "257eab6f29cf3480c536cac73ed3998a",
            "53987b91b94d75d844f865ac0a778e75",
        ],
        "clock remote custody",
    )
    require(
        [f["observed_sha256"] for f in clock["files"]]
        == [
            "5ff66b2229f9c3b57d0fd2fa27e38aa773e9bfb597441c22393f9a65a36fed61",
            "2df463bc8840f1b366d6acd9e361bf58322b844b72f6309641fc66e60735bc30",
            "a4384040d5678893b4df15ea0be2980661023305f5009de876179fab1dea632f",
            "3c14df355b5cf1e6dcf138cf3b3de750f59ad270091e7874c70ad05204fa988d",
        ],
        "clock observed sha256 custody",
    )
    require(clock["files"][0]["shape"] == "721 rows x 5 columns", "clock Fig1 shape")
    require(
        clock["files"][0]["nonblank_numeric_by_ensemble"] == [714, 719, 717, 716, 716],
        "clock Fig1 numeric profile",
    )
    require([f["shape"] for f in clock["files"][1:]] == [
        "12 rows x 17 columns",
        "16 rows x 6 columns",
        "3 rows x 20 columns",
    ], "clock process-data shapes")
    require(not clock["deposited_full_covariance_present"], "clock covariance overclaim")
    require(not clock["beta_TM_identifiable"], "clock beta overclaim")

    force = by_id["C2_YIN_LEVITATED_GRAVITY_DRIVE"]
    require(force["workbook_bytes"] == 15485, "force workbook size")
    require(
        force["workbook_observed_sha256"]
        == "853a2f209d77c2b124d8319f16064944d33d6ba4a73cef7acb8228f9c847fb7b",
        "force workbook custody",
    )
    profile = force["workbook_profile"]["Extended_Figure3d"]
    require(profile == {
        "theoretical_lower_fN": 2.19,
        "theoretical_upper_fN": 2.57,
        "measured_average_fN": 2.33,
        "measured_standard_deviation_fN": 0.33,
    }, "force values")
    require(
        force["workbook_profile"]["Extended_Figure3a"]
        == "91 thermal-noise frequency/PSD rows from 15.310 to 15.400 Hz",
        "force thermal profile",
    )
    require(
        force["workbook_profile"]["Extended_Figure3b"]
        == "91 gravity-drive frequency/PSD rows from 15.310 to 15.400 Hz",
        "force drive profile",
    )
    require(not force["independent_G_crosscheck_possible"], "G overclaim")
    require(not force["beta_TM_identifiable"], "force beta overclaim")

    protected = result["custody"]
    require(not protected["Panda_response_holdout_opened_or_scored"], "Panda holdout")
    require(not protected["Panda_response_filenames_present_in_holdout_directory"], "Panda presence flag")
    holdout = ROOT / "LANE_GRA_J_GRAVITY_HOLDOUT"
    panda = custody["protected_holdout"]
    require(not any((holdout / name).exists() for name in panda["response_filenames"]), "Panda file present")
    metadata = json.loads((holdout / "zenodo_10995225_metadata.json").read_text(encoding="utf-8"))
    observed_panda = {item["key"]: item["checksum"] for item in metadata["files"]}
    require(observed_panda == panda["published_md5"], "Panda inventory changed")
    require(not protected["lineage_labels_manufactured"], "manufactured labels")
    require(not protected["incompatible_roots_pooled"], "pooled roots")
    require(not protected["missing_covariance_replaced_by_independence"], "covariance invention")
    require(not protected["missing_data_called_zero"], "missing-as-zero")
    require(not protected["negative_search_called_physical_null"], "negative-as-null")

    for relative, expected in custody["dependencies"].items():
        path = ROOT / relative
        require(path.is_file() and not path.is_symlink(), f"dependency missing: {relative}")
        require(digest(path) == expected, f"dependency changed: {relative}")

    # The two registers were inputs to the search-time duplicate screen and are
    # then intentionally changed when this accepted result is integrated.  Pin
    # their exact pre-integration bytes as historical context without creating
    # a circular live dependency on their post-integration state.
    contextual = custody["contextual_register_baselines"]
    require(
        contextual["git_commit"] == "a7e54f1cc5295bcfe885415f006521356f683627",
        "register baseline commit",
    )
    require(not contextual["live_hash_required"], "circular live register dependency")
    require(
        contextual["files"]
        == {
            "GRAVITY_EMERGENCE_EXPERIMENT_REGISTER_V001.md": (
                "ae7d2e672b3ba59f9c93d160c8562c95541b8e4c00d68fb21bf27d5315c2b58c"
            ),
            "PROGRAM_EXPERIMENT_REGISTER_V001.md": (
                "084ed543ddaa5c55fd90e12a76c43688016aedd9218fe8850c34a5b30514c0d5"
            ),
        },
        "historical register baselines changed",
    )
    require(
        all((ROOT / relative).is_file() for relative in contextual["files"]),
        "contextual register path absent",
    )

    report = (HERE / "PUBLIC_DATA_SECOND_PASS.md").read_text(encoding="utf-8")
    require("bounded negative result" in report, "bounded-search language missing")
    require("not a physical null" in report, "null ceiling missing")
    require("returned-hit lists" in report, "search reproducibility ceiling missing")
    require("listed repositories were exhaustively enumerated" in report, "exhaustiveness ceiling missing")
    require("No lineage labels were manufactured" in report, "lineage ceiling missing")
    require("Panda response holdout was not opened" in report, "holdout statement missing")

    audit = (HERE / "INDEPENDENT_HOSTILE_AUDIT.md").read_text(encoding="utf-8")
    require("ACCEPT_WITH_EXPLICIT_SEARCH_REPRODUCIBILITY_CEILING" in audit, "independent verdict missing")
    require("all five public files" in audit, "source re-download audit missing")

    reaudit = (HERE / "POST_INTEGRATION_CUSTODY_REAUDIT.md").read_text(encoding="utf-8")
    require("ACCEPT_NONCIRCULAR_CONTEXT_BASELINE_REPAIR" in reaudit, "post-integration verdict missing")
    require("queries, retained roots, and Panda holdout gates are unchanged" in reaudit, "search invariance missing")

    text = (
        "SPAG_PUBLIC_DATA_SECOND_PASS_CHECK: PASS\n"
        "Frozen queries: 28; new component datasets: 2\n"
        "New randomized same-parent lineage datasets: 0\n"
        "Clock root: processed endpoint only; full covariance/lineage absent\n"
        "Force root: fN calibration precedent only; SI covariance/lineage absent\n"
        "Panda response holdout: preserved\n"
        "Search execution: query strings frozen; returned-hit lists not frozen\n"
        "Register custody: exact pre-integration context baselines; no circular live dependency\n"
        "Physical null, beta_TM score, and independent G cross-check: not claimed\n"
    )
    (HERE / "VERIFICATION.txt").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
