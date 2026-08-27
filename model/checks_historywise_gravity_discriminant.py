#!/usr/bin/env python3
"""Sixty focused checks for the pinned formal discriminant and custody gate."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
import hashlib
from pathlib import Path
import shutil
from types import SimpleNamespace
import tempfile
from unittest import mock

import historywise_gravity_discriminant as hgd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LANE = ROOT / hgd.LANE_ID


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def refused(callable_object) -> bool:
    try:
        callable_object()
    except hgd.HistorywiseGravityDiscriminantRefusal:
        return True
    return False


def type_error(callable_object) -> bool:
    try:
        callable_object()
    except TypeError:
        return True
    return False


def set_item(mapping, key, value) -> None:
    mapping[key] = value


def main() -> int:
    checks = 0

    def check(condition: bool) -> None:
        nonlocal checks
        assert condition
        checks += 1

    check(hgd.SCHEMA == "WAC_HISTORYWISE_GRAVITY_DISCRIMINANT_CERTIFICATE_V001")  # 1
    check(hgd.CLAIM_CLASS == "FORMAL_FINITE_GROUP_DISCRIMINANT_ONLY")  # 2
    check(hgd.LANE_ID == LANE.name)  # 3
    check(hgd.MANIFEST_SHA256 == digest(LANE / "MANIFEST.sha256"))  # 4
    check(len(hgd._ARTIFACT_SHA256) == 9)  # 5
    check(set(hgd._ARTIFACT_SHA256) | {"MANIFEST.sha256"} == {p.name for p in LANE.iterdir()})  # 6

    custody = hgd._verify_custody(LANE)
    check(custody.manifest_sha256 == hgd.MANIFEST_SHA256)  # 7
    check(custody.audit_sha256 == digest(LANE / "AUDIT.md"))  # 8
    check(custody.verification_sha256 == digest(LANE / "VERIFICATION.txt"))  # 9
    check(custody.verifier_sha256 == digest(LANE / "verify_historywise_gravity_actualization.py"))  # 10
    check(len(custody.fresh_verifier_stdout_sha256) == 64)  # 11

    result = hgd.historywise_gravity_discriminant()
    check(isinstance(result, hgd.HistorywiseGravityDiscriminant))  # 12
    check(result.claim_class == hgd.CLAIM_CLASS)  # 13
    check(result.manifest_sha256 == hgd.MANIFEST_SHA256)  # 14
    certificate = result.certificate()
    check(certificate["schema"] == hgd.SCHEMA)  # 15
    check(certificate["claim_class"] == hgd.CLAIM_CLASS)  # 16
    check(set(certificate["formal_results"]) == {"HGA1", "HGA1_feedback", "HGA1a", "HGA1b", "HGA2", "HGA3"})  # 17
    check(certificate["formal_results"]["HGA1"] == "EXACT_FINITE_GROUP_ORBIT_NONSELECTION")  # 18
    check(certificate["formal_results"]["HGA1_feedback"] == "EXACT_EQUIVARIANT_FIXED_POINT_ORBIT_CLOSURE")  # 19
    check("FIXED_HISTORY" in certificate["formal_results"]["HGA1a"])  # 20
    check("NOT_A_SAMPLE" in certificate["formal_results"]["HGA1b"])  # 21
    check(certificate["formal_results"]["HGA2"] == "CONDITIONAL_FINITE_AFFINE_MEAN_FIELD_NONSELECTION")  # 22
    check("STABILIZER_CRITERION" in certificate["formal_results"]["HGA3"])  # 23
    check(certificate["discriminant"] == "ENDOGENOUS_EQUIVARIANT_SINGLETON_EXCLUDED_ON_FIXED_POINT_FREE_ORBIT")  # 24
    check(certificate["orienting_input"] == "K_SUBSET_L_MATHEMATICAL_CAPACITY_ONLY")  # 25
    check(certificate["witnesses"]["negative"].endswith("NOT_PHYSICAL_GR"))  # 26
    check(certificate["witnesses"]["positive"].endswith("NOT_A_PHYSICAL_FIELD"))  # 27

    expected_statuses = {
        "physical_GARH_D": "NOT_ADMITTED_BY_THIS_DISCRIMINANT",
        "GARH_Q": "NOT_DERIVED_NOT_FORCED_BY_THIS_DISCRIMINANT",
        "GARH_D_Q_DECISION": "NOT_MADE",
        "objective_actualization": "OPEN_IN_THIS_LANE",
        "physical_gravity": "NO_PROOF_OUTPUT",
        "record_causes_gravity": "NO_PROOF_OUTPUT",
        "Born_law": "NO_PROOF_OUTPUT",
        "general_relativity": "NO_PROOF_OUTPUT",
        "empirical_validation": "NONE",
        "caller_input_scientific_weight": "ZERO",
    }
    check(dict(certificate["statuses"]) == expected_statuses)  # 28
    check(certificate["statuses"]["physical_gravity"] == "NO_PROOF_OUTPUT")  # 29
    check(certificate["statuses"]["record_causes_gravity"] == "NO_PROOF_OUTPUT")  # 30
    check(certificate["statuses"]["Born_law"] == "NO_PROOF_OUTPUT")  # 31
    check(certificate["statuses"]["general_relativity"] == "NO_PROOF_OUTPUT")  # 32
    check(certificate["statuses"]["empirical_validation"] == "NONE")  # 33
    check(certificate["statuses"]["caller_input_scientific_weight"] == "ZERO")  # 34
    check(certificate["authorizations"] and not any(certificate["authorizations"].values()))  # 35
    check(not any(certificate["nonpromotion"].values()))  # 36
    check(certificate["executable_scope"]["statement"] == "EXECUTABLE_DOES_NOT_PROVE_GENERAL_FINITE_GROUP_THEOREM")  # 37
    check(certificate["executable_scope"]["general_finite_group_theorem_machine_proved"] is False)  # 38
    check(certificate["executable_scope"]["physical_or_empirical_proof_weight"] == "ZERO")  # 39
    check(certificate["custody"]["lane_id"] == hgd.LANE_ID)  # 40
    check(certificate["custody"]["manifest_sha256"] == hgd.MANIFEST_SHA256)  # 41
    check(certificate["custody"]["manifest_artifact_count"] == 9)  # 42
    check(certificate["custody"]["audit_disposition"] == hgd.AUDIT_DISPOSITION)  # 43
    check(certificate["custody"]["independent_audit_verdict"] == hgd.INDEPENDENT_AUDIT_VERDICT)  # 44
    check(certificate["custody"]["sealed_verifier_total"] == hgd.VERIFIER_TOTAL)  # 45
    check(certificate["custody"]["sealed_verifier_verdict"] == hgd.VERIFIER_VERDICT)  # 46

    certificate_again = result.certificate()
    check(certificate_again == certificate and certificate_again is not certificate)  # 47
    check(certificate_again["statuses"] is not certificate["statuses"])  # 48
    check(type_error(lambda: set_item(certificate, "schema", "bad")))  # 49
    check(type_error(lambda: set_item(certificate["statuses"], "physical_gravity", "PROVED")))  # 50
    try:
        result._custody = custody
    except FrozenInstanceError:
        frozen = True
    else:
        frozen = False
    check(frozen)  # 51
    check(type_error(lambda: hgd.historywise_gravity_discriminant(True)))  # 52
    check(type_error(lambda: hgd.historywise_gravity_discriminant(root=LANE)))  # 53
    check(type_error(lambda: hgd.historywise_gravity_discriminant_certificate(True)))  # 54
    check(type_error(lambda: hgd.historywise_gravity_discriminant_certificate(packet={})))  # 55

    with tempfile.TemporaryDirectory(prefix="wac-hgd-core-") as directory:
        temporary = Path(directory)
        check(refused(lambda: hgd._verify_custody(temporary / "WRONG_LANE")))  # 56
        symlink_lane = temporary / hgd.LANE_ID
        symlink_lane.symlink_to(LANE, target_is_directory=True)
        check(refused(lambda: hgd._verify_custody(symlink_lane)))  # 57

    with tempfile.TemporaryDirectory(prefix="wac-hgd-core-") as directory:
        copied = Path(directory) / hgd.LANE_ID
        shutil.copytree(LANE, copied)
        (copied / "README.md").write_text("tampered\n", encoding="utf-8")
        check(refused(lambda: hgd._verify_custody(copied)))  # 58

    bad_result = SimpleNamespace(returncode=0, stdout="different\n", stderr="")
    with mock.patch.object(hgd.subprocess, "run", return_value=bad_result):
        check(refused(lambda: hgd._run_fresh_verifier(LANE, "sealed\n")))  # 59

    check(refused(lambda: hgd._sealed_fresh_stdout("EXECUTED CHECKS\n\nTOTAL 63/64 PASS\n\nREPRODUCTION AND SCOPE CHECKS")))  # 60

    assert checks == 60
    print("HISTORYWISE_GRAVITY_DISCRIMINANT_CORE_CHECKS: 60/60 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
