#!/usr/bin/env python3
"""Default-refuted verifier for the URM world-observation input door.

This verifier audits only delegation from the public URM class to the already
separate world-observation contract.  It does not score formation, gravity, a
general theory, universality, or independent experimental reproduction.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Callable


VERIFY = Path(__file__).resolve().parent
LANE = VERIFY.parent
REPO = LANE.parent
MODEL = REPO / "model"
MANIFEST = LANE / "world_observation.json"
DATA = LANE / "world_observation.csv"
STORED_CERTIFICATE = LANE / "world_observation_certificate.json"

# These are the exact inputs commissioned for this audit.  Any change makes the
# result REFUTED until this independent verifier is deliberately recommissioned.
FROZEN_HASHES = {
    "model/project_model.py": "97059fb6d2ee107c3d03363fd2fc9f21c7273aaab9d40da1a9ec07bf537daaae",
    "model/checks_world_observation_urm.py": "b65335dfcdfe6c40c29d263abc2747cc5fdb4455cdce82e34adbe6e7c633f623",
    "model/world_observation.py": "8e511fad0862afa4d48be24b29feca688c4817418e6a83c614d44ca8e760b81f",
    "model/lakeshore_vsm.py": "68221ff4deab5442370e4be57db44916810dd534676d6ead39f378e905370bf1",
    "LANE_T53_A_WORLD_WRITE/OBSERVATION_PROTOCOL.md": "1d32ac55feced8e15d1183e8652dfe7f7d74d19aa1ac0b68fc8b288ae12c6b70",
    "LANE_T53_A_WORLD_WRITE/world_observation.csv": "82fca51d5f7923107763b01340a85842d5a06d770ecc88e8e518a696a2e3a891",
    "LANE_T53_A_WORLD_WRITE/world_observation.json": "15352eaa9c6363c5ff16c02c7c75887941999c9ab914c758c4218ff25b7e1d20",
    "LANE_T53_A_WORLD_WRITE/world_observation_certificate.json": "04beada93f26d2c7f4d3e7b88f016e4b4f27c70e38d8f8812f7e352c3eba0faf",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freeze() -> dict[str, Any]:
    actual = {name: sha256(REPO / name) for name in FROZEN_HASHES}
    return {
        "expected": FROZEN_HASHES,
        "actual": actual,
        "pass": actual == FROZEN_HASHES,
    }


def import_subjects():
    # Keep MODEL on sys.path for the contract's deliberately lazy adapter import.
    # The public check script receives the same path from Python automatically.
    sys.path.insert(0, str(MODEL))
    import project_model  # type: ignore
    import world_observation  # type: ignore
    return project_model, world_observation


def without_docstring(body: list[ast.stmt]) -> list[ast.stmt]:
    if (
        body
        and isinstance(body[0], ast.Expr)
        and isinstance(body[0].value, ast.Constant)
        and isinstance(body[0].value.value, str)
    ):
        return body[1:]
    return body


def ast_exact_delegate_gate() -> dict[str, Any]:
    """Independently constrain the two public methods to transparent delegation."""
    tree = ast.parse((MODEL / "project_model.py").read_text(encoding="utf-8"))
    classes = {
        node.name: node
        for node in tree.body
        if isinstance(node, ast.ClassDef)
    }
    base = classes["ProjectModel"]
    public = classes["URM"]
    base_methods = {
        node.name for node in base.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    if {"world_observation", "world_observation_certificate"} & base_methods:
        raise AssertionError("raw ProjectModel contains a world-observation public door")

    methods = {
        node.name: node
        for node in public.body
        if isinstance(node, ast.FunctionDef)
    }
    load = methods.get("world_observation")
    cert = methods.get("world_observation_certificate")
    if load is None or cert is None:
        raise AssertionError("URM lacks one or both commissioned public methods")
    for method in (load, cert):
        if not (
            len(method.decorator_list) == 1
            and isinstance(method.decorator_list[0], ast.Name)
            and method.decorator_list[0].id == "staticmethod"
        ):
            raise AssertionError(f"{method.name} is not exactly a staticmethod")
        if [argument.arg for argument in method.args.args] != ["manifest_path"]:
            raise AssertionError(f"{method.name} changes the public argument surface")

    load_body = without_docstring(load.body)
    if len(load_body) != 2:
        raise AssertionError("world_observation contains logic beyond import and return")
    imported, returned = load_body
    if not (
        isinstance(imported, ast.ImportFrom)
        and imported.module == "world_observation"
        and [(alias.name, alias.asname) for alias in imported.names]
        == [("load_world_observation", None)]
    ):
        raise AssertionError("world_observation does not import the frozen contract directly")
    if not (
        isinstance(returned, ast.Return)
        and isinstance(returned.value, ast.Call)
        and isinstance(returned.value.func, ast.Name)
        and returned.value.func.id == "load_world_observation"
        and len(returned.value.args) == 1
        and isinstance(returned.value.args[0], ast.Name)
        and returned.value.args[0].id == "manifest_path"
        and not returned.value.keywords
    ):
        raise AssertionError("world_observation is not an exact one-argument delegate")

    cert_body = without_docstring(cert.body)
    if len(cert_body) != 1 or not isinstance(cert_body[0], ast.Return):
        raise AssertionError("world_observation_certificate contains extra behavior")
    value = cert_body[0].value
    if not (
        isinstance(value, ast.Call)
        and not value.args
        and not value.keywords
        and isinstance(value.func, ast.Attribute)
        and value.func.attr == "certificate"
        and isinstance(value.func.value, ast.Call)
        and isinstance(value.func.value.func, ast.Attribute)
        and isinstance(value.func.value.func.value, ast.Name)
        and value.func.value.func.value.id == "URM"
        and value.func.value.func.attr == "world_observation"
        and len(value.func.value.args) == 1
        and isinstance(value.func.value.args[0], ast.Name)
        and value.func.value.args[0].id == "manifest_path"
        and not value.func.value.keywords
    ):
        raise AssertionError("world_observation_certificate is not an exact certificate delegate")
    return {
        "raw_base_methods_absent": True,
        "public_methods_static": True,
        "loader_body_is_exact_delegate": True,
        "certificate_body_is_exact_delegate": True,
        "scientific_scoring_statements": 0,
    }


def runtime_exact_delegate_gate(PM, WO) -> dict[str, Any]:
    calls: list[object] = []
    opaque_path = object()
    sentinel_certificate: dict[str, object] = {"opaque": object()}

    class SentinelObservation:
        def __init__(self) -> None:
            self.certificate_calls = 0

        def certificate(self):
            self.certificate_calls += 1
            return sentinel_certificate

    sentinel_observation = SentinelObservation()
    original = WO.load_world_observation

    def sentinel_loader(path):
        calls.append(path)
        return sentinel_observation

    WO.load_world_observation = sentinel_loader
    try:
        if PM.URM.world_observation(opaque_path) is not sentinel_observation:
            raise AssertionError("URM transformed or replaced the contract result")
        if calls != [opaque_path]:
            raise AssertionError("URM did not pass the path through exactly once")
        calls.clear()
        if PM.URM.world_observation_certificate(opaque_path) is not sentinel_certificate:
            raise AssertionError("URM transformed or replaced the certificate")
        if calls != [opaque_path] or sentinel_observation.certificate_calls != 1:
            raise AssertionError("certificate door did not use exactly one loader/certificate call")
    finally:
        WO.load_world_observation = original

    refusal = WO.ObservationRefusal("INDEPENDENT SENTINEL REFUSAL")

    def refuse(_path):
        raise refusal

    WO.load_world_observation = refuse
    try:
        for method in (PM.URM.world_observation, PM.URM.world_observation_certificate):
            try:
                method(opaque_path)
            except WO.ObservationRefusal as exc:
                if exc is not refusal:
                    raise AssertionError("URM replaced the contract's refusal object")
            else:
                raise AssertionError("URM swallowed a contract refusal")
    finally:
        WO.load_world_observation = original

    return {
        "argument_passed_by_identity": True,
        "observation_returned_by_identity": True,
        "certificate_returned_by_identity": True,
        "refusal_propagated_by_identity": True,
    }


def raw_base_gate(PM) -> dict[str, Any]:
    names = ("world_observation", "world_observation_certificate")
    if any(name in PM.ProjectModel.__dict__ for name in names):
        raise AssertionError("world input door appears in ProjectModel.__dict__")
    if any(hasattr(PM.ProjectModel, name) for name in names):
        raise AssertionError("world input door resolves on the raw ProjectModel class")
    if any(hasattr(PM.ProjectModel(), name) for name in names):
        raise AssertionError("world input door resolves on a raw ProjectModel instance")
    return {
        "project_model_class_unavailable": True,
        "project_model_instance_unavailable": True,
        "urm_class_available": all(hasattr(PM.URM, name) for name in names),
    }


def actual_bundle_gate(PM, WO) -> dict[str, Any]:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    with DATA.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, strict=True)
        if reader.fieldnames != manifest["columns"]:
            raise AssertionError("independent CSV header differs from manifest")
        raw_rows = list(reader)
    if len(raw_rows) != 760:
        raise AssertionError(f"independent CSV row count is {len(raw_rows)}, not 760")
    if len({row["row_id"] for row in raw_rows}) != 760:
        raise AssertionError("independent CSV scan found duplicate row IDs")
    if hashlib.sha256(DATA.read_bytes()).hexdigest() != manifest["data_sha256"]:
        raise AssertionError("independent data hash differs from manifest")
    for item in manifest["source_artifacts"]:
        if sha256(LANE / item["path"]) != item["sha256"]:
            raise AssertionError(f"source hash drift: {item['path']}")

    via_urm = PM.URM.world_observation(MANIFEST)
    direct = WO.load_world_observation(MANIFEST)
    if via_urm != direct:
        raise AssertionError("URM observation differs from direct contract result")
    if len(via_urm.rows) != 760:
        raise AssertionError("URM did not expose all 760 rows")
    if [row["row_id"] for row in via_urm.rows] != [row["row_id"] for row in raw_rows]:
        raise AssertionError("URM row identity/order differs from independent CSV scan")

    via_certificate = PM.URM.world_observation_certificate(MANIFEST)
    direct_certificate = direct.certificate()
    stored_certificate = json.loads(STORED_CERTIFICATE.read_text(encoding="utf-8"))
    if via_certificate != direct_certificate or via_certificate != stored_certificate:
        raise AssertionError("URM, direct, and stored certificates differ")
    if WO.certificate_json(via_urm).encode("utf-8") != STORED_CERTIFICATE.read_bytes():
        raise AssertionError("URM certificate serialization differs from frozen artifact")

    expected_nonproof = {
        "scientific_verdict": "NONE_NOT_SCORED",
        "independent_reproduction_attested": False,
        "record_formation_proof_authorized": False,
        "universal_claim_authorized": False,
        "public_urm_registration_authorized": False,
    }
    actual_nonproof = {key: via_certificate.get(key) for key in expected_nonproof}
    if actual_nonproof != expected_nonproof:
        raise AssertionError(f"certificate scientifically elevates the input: {actual_nonproof}")
    if via_certificate.get("scope_classification") != "CONFIGURATION_EVIDENCE_ONLY":
        raise AssertionError("retrospective bundle was elevated beyond configuration evidence")
    if via_certificate.get("coverage", {}).get("row_count") != 760:
        raise AssertionError("certificate does not report 760 rows")

    return {
        "independent_csv_rows": 760,
        "urm_rows": len(via_urm.rows),
        "row_identity_and_order_exact": True,
        "urm_equals_direct_contract": True,
        "stored_certificate_exact": True,
        "scope_classification": via_certificate["scope_classification"],
        **actual_nonproof,
    }


def end_to_end_refusal_gate(PM, WO) -> dict[str, Any]:
    """Use an independently damaged bundle, not the commissioned unit-check fixture."""
    original = json.loads(MANIFEST.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix="wac-urm-refusal-", dir="/private/tmp") as temp:
        root = Path(temp)
        shutil.copy2(LANE / original["protocol_path"], root / original["protocol_path"])
        shutil.copy2(LANE / original["data_file"], root / original["data_file"])
        damaged = dict(original)
        damaged["data_sha256"] = "0" * 64
        path = root / "damaged.json"
        path.write_text(json.dumps(damaged, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        outcomes: list[tuple[type[BaseException], str]] = []
        for loader in (WO.load_world_observation, PM.URM.world_observation):
            try:
                loader(path)
            except WO.ObservationRefusal as exc:
                outcomes.append((type(exc), str(exc)))
            else:
                raise AssertionError("damaged data hash was accepted")
    if len(outcomes) != 2 or outcomes[0] != outcomes[1]:
        raise AssertionError("URM refusal type/message differs from direct contract")
    if "data_sha256 mismatch" not in outcomes[0][1]:
        raise AssertionError("damaged fixture refused for an unexpected reason")

    missing = LANE / "definitely-absent" / "manifest.json"
    missing_outcomes: list[tuple[type[BaseException], str]] = []
    for loader in (WO.load_world_observation, PM.URM.world_observation):
        try:
            loader(missing)
        except WO.ObservationRefusal as exc:
            missing_outcomes.append((type(exc), str(exc)))
        else:
            raise AssertionError("missing bundle was accepted")
    if len(missing_outcomes) != 2 or missing_outcomes[0] != missing_outcomes[1]:
        raise AssertionError("missing-input refusal did not propagate exactly")
    return {
        "damaged_hash_refused": True,
        "damaged_hash_refusal_exact": True,
        "missing_manifest_refused": True,
        "missing_manifest_refusal_exact": True,
    }


def commissioned_check_gate() -> dict[str, Any]:
    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    process = subprocess.run(
        [sys.executable, "-B", str(MODEL / "checks_world_observation_urm.py")],
        cwd=REPO,
        env=env,
        capture_output=True,
        check=False,
    )
    stdout = process.stdout.decode("utf-8", errors="replace").strip()
    stderr = process.stderr.decode("utf-8", errors="replace").strip()
    expected = "URM_WORLD_OBSERVATION_CHECKS: 5/5 PASS"
    if process.returncode != 0 or stdout != expected or stderr:
        raise AssertionError(
            f"commissioned check failed: rc={process.returncode}, stdout={stdout!r}, stderr={stderr!r}"
        )
    return {"returncode": 0, "stdout": stdout, "stderr": stderr}


def run_case(name: str, operation: Callable[[], dict[str, Any]]) -> dict[str, Any]:
    try:
        detail = operation()
    except Exception as exc:  # default-refuted: every unexpected condition survives.
        return {
            "name": name,
            "pass": False,
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
    return {"name": name, "pass": True, "detail": detail}


def write_results(result: dict[str, Any]) -> None:
    (VERIFY / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    failed = [case["name"] for case in result["cases"] if not case["pass"]]
    lines = [
        f"T53A URM OBSERVATION DELEGATE: {result['verdict']}",
        f"admissible_for_input_door_integration={str(result['admissible_for_input_door_integration']).lower()}",
        f"cases={result['case_count']}",
        f"failed_cases={','.join(failed) if failed else 'NONE'}",
        "scientific_proof_weight=ZERO",
        "scientific_claim_scored=false",
        "record_formation_proof_authorized=false",
        "gravity_emergence_proof_authorized=false",
        "universal_claim_authorized=false",
        "independent_reproduction_authorized=false",
    ]
    (VERIFY / "RESULT.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    if failed:
        findings = "\n".join(
            f"- `{case['name']}`: {case.get('error_type')}: {case.get('error')}"
            for case in result["cases"]
            if not case["pass"]
        )
    else:
        findings = "- No surviving defect in the narrow delegation gate."
    report = f"""# D24 — URM world-observation delegate audit

Posture: **REFUTED by default**.

Verdict: **{result['verdict']}** for the narrow public input-door integration only.
Admissible for that narrow integration: **{str(result['admissible_for_input_door_integration']).lower()}**.

## Findings

{findings}

## Verifier correction log

- The first run removed `model/` from `sys.path` immediately after importing the top-level
  modules.  That verifier-only environment error prevented the contract's deliberately lazy
  `lakeshore_vsm` import.  The verifier now retains the same module path that the public check
  script receives automatically.  No subject, manifest, data, protocol, or certificate file
  changed; the pinned pre/post hashes show that correction is confined to this audit directory.

## Boundary

This verifier gives the result **zero scientific-proof weight**.  It establishes only that
`URM.world_observation()` and `URM.world_observation_certificate()` transparently expose the
frozen measurement contract, preserve its refusals, expose the actual 760-row bundle, remain
absent from raw `ProjectModel`, and add no scientific verdict.  It does not validate record
formation, gravity emergence, universality, physical-origin authenticity, prospective agreement,
or independent experimental reproduction.  Robustness of the underlying contract against its
full malformed-input battery is owned by the separate world-ingest verifier; this result neither
replaces nor upgrades that verdict.  Any pinned contract-byte change automatically refutes this
delegate result until a new independent audit.
"""
    (VERIFY / "D24_AUDIT.md").write_text(report, encoding="utf-8")


def main() -> int:
    initial = freeze()
    PM, WO = import_subjects()
    cases = [
        run_case("ast_exact_delegate", ast_exact_delegate_gate),
        run_case("runtime_exact_delegate_and_refusal_identity", lambda: runtime_exact_delegate_gate(PM, WO)),
        run_case("raw_project_model_unavailable", lambda: raw_base_gate(PM)),
        run_case("actual_760_row_bundle", lambda: actual_bundle_gate(PM, WO)),
        run_case("end_to_end_refusal_propagation", lambda: end_to_end_refusal_gate(PM, WO)),
        run_case("commissioned_check", commissioned_check_gate),
    ]
    final = freeze()
    all_pass = initial["pass"] and final["pass"] and all(case["pass"] for case in cases)
    result = {
        "schema": "WAC_T53A_URM_OBSERVATION_DELEGATE_VERIFY_V001",
        "default": "REFUTED",
        "verdict": "NOT_REFUTED" if all_pass else "REFUTED",
        "admissible_for_input_door_integration": all_pass,
        "scientific_proof_weight": "ZERO",
        "scientific_claim_scored": False,
        "record_formation_proof_authorized": False,
        "gravity_emergence_proof_authorized": False,
        "universal_claim_authorized": False,
        "independent_reproduction_authorized": False,
        "initial_freeze": initial,
        "final_freeze": final,
        "case_count": len(cases),
        "cases": cases,
    }
    write_results(result)
    print((VERIFY / "RESULT.txt").read_text(encoding="utf-8"), end="")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
