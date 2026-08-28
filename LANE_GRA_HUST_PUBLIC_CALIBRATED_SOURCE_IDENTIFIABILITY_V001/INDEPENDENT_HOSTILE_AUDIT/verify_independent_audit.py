#!/usr/bin/env python3
"""Verifier for the independent HUST calibrated-source hostile audit."""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


class Checks:
    def __init__(self) -> None:
        self.total = 0
        self.failures: list[str] = []

    def check(self, condition: bool, label: str) -> None:
        self.total += 1
        if not condition:
            self.failures.append(label)

    def finish(self) -> int:
        if self.failures:
            for failure in self.failures:
                print(f"FAIL: {failure}")
            print(f"FAIL: {self.total-len(self.failures)}/{self.total} checks passed")
            return 1
        print(f"PASS: {self.total}/{self.total} independent-audit checks passed")
        return 0


def verify_manifest(checks: Checks, path: Path, base: Path) -> None:
    checks.check(path.is_file(), f"manifest exists {path.name}")
    if not path.is_file():
        return
    for line in path.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        target = base / relative
        checks.check(target.is_file(), f"manifest target exists {relative}")
        if target.is_file():
            checks.check(sha256(target) == expected, f"manifest target hash {relative}")


def main() -> int:
    checks = Checks()

    # The builder/core packet is frozen and must remain exactly as sealed.
    verify_manifest(checks, LANE / "MANIFEST.sha256", LANE)
    core_seal = (LANE / "LANE_SEAL.sha256").read_text().splitlines()
    checks.check(len(core_seal) == 2, "builder seal has two bindings")
    for line in core_seal:
        expected, relative = line.split("  ", 1)
        checks.check(sha256(LANE / relative) == expected, f"builder seal {relative}")

    audit_result = json.loads((AUDIT / "AUDIT_RESULT.json").read_text())
    transcription = json.loads((AUDIT / "INDEPENDENT_SOURCE_TRANSCRIPTION.json").read_text())
    reacquisition = json.loads((AUDIT / "SOURCE_REACQUISITION.json").read_text())
    checks.check(audit_result["arithmetic_and_custody_all_pass"] is True, "audit arithmetic/custody pass")
    checks.check(audit_result["disposition"] == "PASS__M1_N1_CORE_REPAIRS_INCORPORATED__PUBLICATION_SAFE_WITH_INTRINSIC_CEILINGS", "audit disposition")
    repair = audit_result["core_repair_audit"]
    checks.check(repair["all_pass"] is True, "all core repair checks pass")
    checks.check(repair["M1_disposition"] == "CLOSED_IN_CORE_LANGUAGE_AND_EXECUTABLE_SCHEMA", "M1 closed in core")
    checks.check(repair["N1_disposition"] == "CLOSED_IN_CORE_LANGUAGE_AND_EXECUTABLE_SCHEMA", "N1 closed in core")
    for label, value in repair["checks"].items():
        checks.check(value is True, f"core repair check {label}")
    checks.check("publication-safe unchanged" in repair["publication_disposition"], "core publication disposition")
    checks.check(len(audit_result["official_source_transcription_checks"]["manual_correction_fields"]) == 80, "80 manual correction checks")
    checks.check(len(audit_result["official_source_transcription_checks"]["official_html_error_budget_fields"]) == 42, "42 official HTML vector checks")
    for item in audit_result["official_source_transcription_checks"]["manual_correction_fields"]:
        checks.check(item["pass"], f"manual source field {item['field']}")
    for item in audit_result["official_source_transcription_checks"]["official_html_error_budget_fields"]:
        checks.check(item["pass"], f"official table field {item['field']}")

    checks.check(len(audit_result["local_source_custody"]) == 3, "three local primary sources")
    for item in audit_result["local_source_custody"]:
        checks.check(item["pass"], f"local source custody {item['id']}")
    checks.check(len(audit_result["sealed_dependency_custody"]) == 7, "seven sealed upstream dependencies")
    for item in audit_result["sealed_dependency_custody"]:
        checks.check(item["pass"], f"dependency custody {item['id']}")
    for item in reacquisition["objects"]:
        checks.check(item["sealed_sha256_match"], f"fresh reacquisition {item['id']}")

    checks.check(len(audit_result["independent_forwards"]["AAF"]) == 3, "three AAF rows")
    checks.check(len(audit_result["independent_forwards"]["TOS"]) == 7, "seven ToS rows")
    for row in audit_result["independent_forwards"]["AAF"]:
        checks.check(row["reported_match"], f"AAF forward {row['id']}")
        checks.check(row["source_response_match"], f"AAF response {row['id']}")
        checks.check(row["source_formula_matches"], f"AAF sinc formulas {row['id']}")
        checks.check(abs(row["deprocessed_identity_relative_residual"]) < 5e-16, f"AAF deprocessing identity {row['id']}")
    for row in audit_result["independent_forwards"]["TOS"]:
        checks.check(row["reported_match"], f"ToS forward {row['id']}")
        checks.check(row["source_response_match"], f"ToS response {row['id']}")
        checks.check(row["signed_anelastic_ppm"] < 0, f"ToS signed anelastic {row['id']}")
        checks.check(row["magnetic_ppm"] > 0, f"ToS positive magnetic {row['id']}")
        checks.check(abs(row["deprocessed_identity_relative_residual"]) < 5e-16, f"ToS deprocessing identity {row['id']}")
        checks.check(abs(row["multiplicative_minus_additive_ppm"]) < 0.011, f"ToS composition scale {row['id']}")

    covariance = audit_result["independent_covariance"]
    checks.check(math.isclose(covariance["AAF_combined_standard_u_ppm"], 11.61623880460628, rel_tol=1e-13), "AAF covariance arithmetic")
    checks.check(math.isclose(covariance["TOS_combined_standard_u_ppm"], 11.63746506158266, rel_tol=1e-13), "ToS covariance arithmetic")
    checks.check(covariance["AAF_positive_definite"], "AAF covariance positive definite")
    checks.check(covariance["TOS_run_positive_definite"], "ToS run covariance positive definite")
    checks.check(covariance["TOS_fibre_positive_definite"], "ToS fibre covariance positive definite")
    for label, value in audit_result["covariance_matches_frozen_result"].items():
        checks.check(value, f"frozen covariance match {label}")

    accepted = audit_result["accepted_or_CODATA_G_audit"]
    checks.check(accepted["numeric_input_used_in_independent_reconstruction"] is False, "no accepted G input")
    checks.check(accepted["frozen_result_declared_inputs"] == [], "frozen no accepted G list")
    ident = audit_result["identifiability_audit"]
    checks.check("independently owned" in ident["critical_ownership_qualification"], "remainder ownership qualification")
    checks.check("does not prove ten independent" in ident["joint_dimension_qualification"], "joint dimension qualification")
    checks.check("not sufficient" in ident["raw_reanalysis_qualification"], "raw reanalysis qualification")
    domain = audit_result["standard_uncertainty_domain_audit"]
    checks.check(domain["authors_model_compact_display_bands_possible"] is True, "authors-model compact display possible")
    checks.check(domain["independent_deterministic_admissible_domain_owned"] is False, "no independent deterministic domain")

    report = (AUDIT / "INDEPENDENT_HOSTILE_AUDIT.md").read_text()
    for phrase in (
        "one independently owned physical-harmonic",
        "Calling the scalar simply “unpublished” would be false",
        "not recovery of a raw",
        "not a source-model-free numerator",
        "deterministic admissible set",
        "not a new measurement",
        "### M1 — closed in the repaired core",
        "### N1 — closed in the repaired core",
        "publication-safe unchanged",
    ):
        checks.check(phrase in report, f"report claim boundary: {phrase}")

    replay = subprocess.run(
        [sys.executable, str(AUDIT / "audit_calibrated_source.py"), "--check"],
        cwd=AUDIT,
        text=True,
        capture_output=True,
        check=False,
    )
    checks.check(replay.returncode == 0, "audit deterministic replay")
    checks.check("PASS" in replay.stdout, "audit replay emits PASS")

    builder_verify = subprocess.run(
        [sys.executable, str(LANE / "verify_hust_public_calibrated_source.py")],
        cwd=LANE,
        text=True,
        capture_output=True,
        check=False,
    )
    checks.check(builder_verify.returncode == 0, "frozen builder verifier")
    checks.check("PASS 192/192" in builder_verify.stdout, "frozen builder verifier count")

    builder_replay = subprocess.run(
        [sys.executable, str(LANE / "analyze_hust_public_calibrated_source.py")],
        cwd=LANE,
        text=True,
        capture_output=True,
        check=False,
    )
    checks.check(builder_replay.returncode == 0, "frozen builder analyzer replay")
    checks.check(builder_replay.stdout == (LANE / "RESULT.json").read_text(), "frozen RESULT exact replay")

    for path in AUDIT.iterdir():
        if path.is_file() and path.suffix in (".py", ".md", ".json", ".sha256"):
            data = path.read_bytes()
            checks.check(not any(byte < 32 and byte not in (9, 10, 13) for byte in data), f"text control bytes {path.name}")

    verify_manifest(checks, AUDIT / "AUDIT_MANIFEST.sha256", AUDIT)
    seal_path = AUDIT / "AUDIT_SEAL.json"
    checks.check(seal_path.is_file(), "audit seal exists")
    if seal_path.is_file():
        seal = json.loads(seal_path.read_text())
        checks.check(seal["disposition"] == audit_result["disposition"], "seal disposition")
        checks.check(seal["builder_core_files_edited"] is False, "seal frozen-core declaration")
        for field, relative in {
            "audit_manifest_sha256": "AUDIT_MANIFEST.sha256",
            "audit_result_sha256": "AUDIT_RESULT.json",
            "audit_report_sha256": "INDEPENDENT_HOSTILE_AUDIT.md",
            "audit_executable_sha256": "audit_calibrated_source.py",
            "audit_verifier_sha256": "verify_independent_audit.py",
        }.items():
            checks.check(seal[field] == sha256(AUDIT / relative), f"seal hash {relative}")

    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
