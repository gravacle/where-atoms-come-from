#!/usr/bin/env python3
"""Independent verifier for the HUST public calibrated source lane."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

EXPECTED_LOCAL_HASHES = {
    "THEOREM.md": "ce2177498cdae8a4fe611deb8a84526554177bd33b9fa59eff48f7b5f7391aa2",
    "CALIBRATION_FIELDS.json": "535d766d7ac69871e47471e6f49bcf02a4d96e9ba891b32c1415fe719b614493",
    "SOURCE_CUSTODY.json": "15615db166d816cad83dc9e7a7bfe6ee5253275a88c4841424f3308e96cd489d",
    "analyze_hust_public_calibrated_source.py": "ac50cc2d4330e9a6c71b74116a2d589ff15262286fcd5a37d254f7fcdf287766",
    "RESULT.json": "8b4ec6841a1032680218383778552e685f172eada6fbcd79edebfae5fc5d976f"
}

EXPECTED_SOURCE_HASHES = {
    "SOURCE/HUST_2018_main_article_public_mirror.pdf": "40756ec0fb8f00c1fde31020b294521a3b220a196bef884a2ea5f3534d77dfaa",
    "SOURCE/41586_2018_431_MOESM1_ESM.pdf": "5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb",
    "SOURCE/nature_main_table1_error_budget.html": "23436d4be7600a7a9dffa02cc4167a20b6eea032a181e77899bb57bb90aa02e9"
}

EXPECTED_MASS_SUMS = {
    "AAF-I": 22.40,
    "AAF-II": 22.89,
    "AAF-III": 22.29,
    "TOS-I-F1-first": 90.33,
    "TOS-I-F1-repeat": 90.33,
    "TOS-I-F2": 92.05,
    "TOS-I-F3-first": 91.24,
    "TOS-I-F3-repeat": 91.24,
    "TOS-II-F4-first": 88.70,
    "TOS-II-F4-repeat": 88.70,
}

EXPECTED_ANELASTIC = {
    "TOS-I-F1-first": -6.01,
    "TOS-I-F1-repeat": -6.01,
    "TOS-I-F2": -8.38,
    "TOS-I-F3-first": -5.68,
    "TOS-I-F3-repeat": -5.68,
    "TOS-II-F4-first": -6.92,
    "TOS-II-F4-repeat": -6.92,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def close(left: float, right: float, *, absolute: float = 1e-12, relative: float = 1e-12) -> bool:
    return math.isclose(left, right, abs_tol=absolute, rel_tol=relative)


def positive_definite(matrix: list[list[float]], tolerance: float = 1e-12) -> bool:
    size = len(matrix)
    lower = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            value = matrix[row][column] - sum(
                lower[row][index] * lower[column][index] for index in range(column)
            )
            if row == column:
                if value <= tolerance:
                    return False
                lower[row][column] = math.sqrt(value)
            else:
                lower[row][column] = value / lower[column][column]
    return True


class Audit:
    def __init__(self) -> None:
        self.total = 0
        self.failed: list[str] = []

    def check(self, condition: bool, label: str) -> None:
        self.total += 1
        if not condition:
            self.failed.append(label)

    def finish(self) -> None:
        if self.failed:
            print(f"FAIL {self.total - len(self.failed)}/{self.total}")
            for label in self.failed:
                print(f"FAILED: {label}")
            raise SystemExit(1)
        print(f"PASS {self.total}/{self.total}")


def main() -> None:
    audit = Audit()

    for relative, expected in EXPECTED_LOCAL_HASHES.items():
        path = HERE / relative
        audit.check(path.is_file(), f"local file exists: {relative}")
        audit.check(sha256(path) == expected, f"local hash: {relative}")
    for relative, expected in EXPECTED_SOURCE_HASHES.items():
        path = HERE / relative
        audit.check(path.is_file(), f"source exists: {relative}")
        audit.check(sha256(path) == expected, f"source hash: {relative}")

    for relative in [
        "THEOREM.md", "README.md", "CALIBRATION_FIELDS.json", "SOURCE_CUSTODY.json",
        "analyze_hust_public_calibrated_source.py", "verify_hust_public_calibrated_source.py",
        "RESULT.json",
    ]:
        data = (HERE / relative).read_bytes()
        audit.check(b"\x00" not in data, f"no NUL: {relative}")
        audit.check(b"\r" not in data, f"LF-only: {relative}")
    for relative in ["SOURCE/HUST_2018_main_article_public_mirror.pdf", "SOURCE/41586_2018_431_MOESM1_ESM.pdf"]:
        audit.check((HERE / relative).read_bytes()[:5] == b"%PDF-", f"PDF magic: {relative}")

    fields = json.loads((HERE / "CALIBRATION_FIELDS.json").read_text(encoding="utf-8"))
    custody = json.loads((HERE / "SOURCE_CUSTODY.json").read_text(encoding="utf-8"))
    result = json.loads((HERE / "RESULT.json").read_text(encoding="utf-8"))
    nominal = json.loads((ROOT / "LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001/RESULT.json").read_text())
    conditional = json.loads((ROOT / "LANE_GRA_HUST_CONDITIONAL_HOMOGENEOUS_G_CROSSCHECK_V001/RESULT.json").read_text())

    audit.check(result["accepted_or_CODATA_G_numeric_inputs"] == [], "no accepted G numeric inputs")
    audit.check(custody["input_exclusion"]["accepted_or_CODATA_G"] == "ABSENT", "custody excludes accepted G")
    audit.check("nonpublisher public" in custody["sources"][0]["custody_note"].lower(), "mirror custody disclosed")
    audit.check(result["status"].startswith("PASS__PUBLIC_CALIBRATED"), "result status")
    audit.check(result["date"] == "2026-08-27", "result date")
    audit.check(fields["date"] == "2026-08-27", "fields date")

    spec = importlib.util.spec_from_file_location(
        "hust_calibrated_analysis", HERE / "analyze_hust_public_calibrated_source.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    regenerated = module.build_result()
    audit.check(regenerated == result, "RESULT exactly regenerates")

    analyzer_text = (HERE / "analyze_hust_public_calibrated_source.py").read_text(encoding="utf-8")
    comparator_marker = analyzer_text.index("# Comparator attachment")
    audit.check(analyzer_text.index('"authors_derived_G_SI"') > comparator_marker, "authors G attached only after primary")
    audit.check(analyzer_text.index('"processed_kernel_kg_m-3"') > comparator_marker, "processed kernel attached only after primary")

    nominal_aaf = {row["id"]: row for row in nominal["AAF"]}
    nominal_tos = {row["id"]: row for row in nominal["TOS"]}
    conditional_aaf = {row["id"]: row for row in conditional["AAF"]}
    conditional_tos = {row["id"]: row for row in conditional["TOS"]}

    for row in result["AAF"]:
        row_id = row["id"]
        expected_sum = EXPECTED_MASS_SUMS[row_id]
        calibration = fields["AAF"][row_id]
        independent_sum = sum(value[0] for value in calibration["mass_model_G_corrections_ppm"].values())
        independent_u = math.sqrt(sum(value[1] ** 2 for value in calibration["mass_model_G_corrections_ppm"].values()))
        k0 = nominal_aaf[row_id]["nominal_homogeneous_coefficient_kg_m-3"]
        kpartial = k0 / (1.0 + independent_sum * 1e-6)
        source = conditional_aaf[row_id]
        g = source["response_alpha_nrad_s-2"] * 1e-9 * source["mechanical_factor_held_at_displayed_value"] / kpartial
        audit.check(close(independent_sum, expected_sum), f"AAF mass sum {row_id}")
        audit.check(close(row["public_mass_model_G_correction_sum_ppm"], independent_sum), f"AAF stored mass sum {row_id}")
        audit.check(close(row["public_mass_model_correction_standard_u_RSS_ppm"], independent_u), f"AAF mass u {row_id}")
        audit.check(close(row["public_calibrated_partial_kernel_kg_m-3"], kpartial), f"AAF partial kernel {row_id}")
        audit.check(close(row["primary_partial_G_SI"], g, absolute=1e-24), f"AAF primary G {row_id}")
        audit.check(abs(row["identity_relative_residual"]) < 5e-16, f"AAF deprocessing identity {row_id}")
        audit.check(row["identified_family"]["public_packet_owned_compact_domain"] is None, f"AAF no packet interval {row_id}")
        formula = row["acquisition_formula_check_ppm"]
        audit.check(abs(formula["data_average_calculated"] - formula["data_average_published"]) < 0.01, f"AAF averaging formula {row_id}")
        audit.check(abs(formula["numeric_derivative_calculated"] - formula["numeric_derivative_published"]) < 0.01, f"AAF derivative formula {row_id}")
        comparator = row["post_calculation_comparator"]
        audit.check("QUARANTINED" in comparator["role"], f"AAF comparator quarantine {row_id}")
        audit.check(comparator["absolute_remainder_reduction_factor"] > 1.7, f"AAF remainder reduced {row_id}")

    for row in result["TOS"]:
        row_id = row["id"]
        expected_sum = EXPECTED_MASS_SUMS[row_id]
        calibration = fields["TOS"][row_id]
        independent_sum = sum(value[0] for value in calibration["mass_model_G_corrections_ppm"].values())
        independent_u = math.sqrt(sum(value[1] ** 2 for value in calibration["mass_model_G_corrections_ppm"].values()))
        k0 = nominal_tos[row_id]["nominal_homogeneous_Delta_Cg_over_I_kg_m-3"]
        kpartial = k0 / (1.0 + independent_sum * 1e-6)
        source = conditional_tos[row_id]
        anelastic = calibration["anelastic_G_correction_ppm"][0]
        magnetic = calibration["magnetic_damper_G_correction_ppm"][0]
        g = source["response_Delta_omega2_s-2"] * (1.0 + (anelastic + magnetic) * 1e-6) / kpartial
        audit.check(close(independent_sum, expected_sum), f"ToS mass sum {row_id}")
        audit.check(close(row["public_mass_model_G_correction_sum_ppm"], independent_sum), f"ToS stored mass sum {row_id}")
        audit.check(close(row["public_mass_model_correction_standard_u_RSS_ppm"], independent_u), f"ToS mass u {row_id}")
        audit.check(close(row["public_calibrated_partial_kernel_kg_m-3"], kpartial), f"ToS partial kernel {row_id}")
        audit.check(close(row["signed_anelastic_G_correction_ppm"], EXPECTED_ANELASTIC[row_id]), f"ToS signed anelastic {row_id}")
        audit.check(row["signed_anelastic_G_correction_ppm"] < 0.0, f"ToS anelastic sign {row_id}")
        audit.check(close(row["primary_partial_G_SI"], g, absolute=1e-24), f"ToS primary G {row_id}")
        audit.check(abs(row["identity_relative_residual"]) < 5e-16, f"ToS deprocessing identity {row_id}")
        audit.check(row["identified_family"]["public_packet_owned_compact_domain"] is None, f"ToS no packet interval {row_id}")
        audit.check("removes c_f" in row["identified_family"]["reason"], f"ToS c_f closure text {row_id}")
        audit.check("independently reconstructed" in row["identified_family"]["reason"], f"ToS independent ownership text {row_id}")
        comparator = row["post_calculation_comparator"]
        audit.check("QUARANTINED" in comparator["role"], f"ToS comparator quarantine {row_id}")
        audit.check(comparator["absolute_remainder_reduction_factor"] > 8.0, f"ToS remainder reduced {row_id}")

    aaf_residuals = [abs(row["post_calculation_comparator"]["authors_processed_minus_public_partial_relative_ppm"]) for row in result["AAF"]]
    tos_residuals = [abs(row["post_calculation_comparator"]["authors_processed_minus_public_partial_relative_ppm"]) for row in result["TOS"]]
    audit.check(max(aaf_residuals) < 30.0, "AAF comparator remainder below 30 ppm")
    audit.check(max(tos_residuals) < 12.12, "ToS comparator remainder below 12.12 ppm")
    audit.check(min(tos_residuals) < 0.57, "ToS best comparator remainder below 0.57 ppm")

    density = result["material_density_inventory"]
    aaf_air = density["AAF_air_check_at_article_approximate_rho_air"]
    audit.check(7964.0 < aaf_air["rho_sphere_aggregate_kg_m-3"] < 7967.0, "AAF sphere density from public masses/diameters")
    audit.check(148.0 < aaf_air["rho_air_over_rho_sphere_ppm"] < 148.3, "approximate air-density physics check")

    covariance = result["published_category_covariance"]
    aaf_cov = covariance["AAF"]["covariance_ppm2"]
    tos_run_cov = covariance["TOS"]["run_covariance_ppm2"]
    tos_fibre_cov = covariance["TOS"]["fibre_covariance_ppm2"]
    audit.check(all(close(aaf_cov[i][j], aaf_cov[j][i]) for i in range(3) for j in range(3)), "AAF covariance symmetric")
    audit.check(all(close(tos_run_cov[i][j], tos_run_cov[j][i]) for i in range(7) for j in range(7)), "ToS run covariance symmetric")
    audit.check(all(close(tos_fibre_cov[i][j], tos_fibre_cov[j][i]) for i in range(4) for j in range(4)), "ToS fibre covariance symmetric")
    audit.check(positive_definite(aaf_cov), "AAF covariance positive definite")
    audit.check(positive_definite(tos_run_cov), "ToS run covariance positive definite")
    audit.check(positive_definite(tos_fibre_cov), "ToS fibre covariance positive definite")
    audit.check(abs(covariance["AAF"]["implied_combined_standard_u_ppm"] - 11.61) < 0.01, "AAF combined uncertainty closure")
    audit.check(abs(covariance["TOS"]["implied_combined_standard_u_ppm"] - 11.64) < 0.01, "ToS combined uncertainty closure")
    audit.check(len(covariance["TOS"]["inferred_same_fibre_shared_background"]) == 3, "three ToS shared-background inferences")
    audit.check("not a measured raw-data covariance" in covariance["ceiling"], "covariance ceiling")

    theorem = result["minimal_independent_remainder_theorem"]
    audit.check(len(theorem["row_remainders"]) == 10, "ten minimal row remainders")
    audit.check(len(set(theorem["row_remainders"])) == 10, "row remainder ids unique")
    audit.check("one scalar" in theorem["point_summary_level"], "one-scalar sufficiency stated")
    audit.check("at least one" in theorem["point_summary_level"], "one-scalar necessity stated")
    audit.check("authors' processed kernels" in theorem["ownership_qualification"].lower(), "public comparator ownership stated")
    audit.check("not asserted" in theorem["ownership_qualification"], "row-coordinate independence not overclaimed")
    audit.check("not a claim" in theorem["compact_interval_result"], "no overclaim of mathematical unboundedness")
    audit.check("authors-model conventional display band" in theorem["compact_interval_result"], "authors-model display-band boundary")
    audit.check(len(result["strict_ceilings"]) == 6, "six strict ceilings")
    audit.check(any("not raw numerator" in ceiling for ceiling in result["strict_ceilings"]), "raw numerator ceiling")
    audit.check(any("declared multiplicative" in ceiling for ceiling in result["strict_ceilings"]), "ToS deprocessing convention ceiling")
    audit.check(any("RGRL" in ceiling for ceiling in result["strict_ceilings"]), "gravity-emergence quarantine")

    audit.finish()


if __name__ == "__main__":
    main()
