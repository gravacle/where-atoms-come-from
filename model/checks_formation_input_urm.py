#!/usr/bin/env python3
"""Six fixed checks for the public URM-only V002 formation delegates."""

from __future__ import annotations

import ast
from pathlib import Path
import tempfile

from checks_formation_input import write_bundle
from formation_input import (
    FormationRefusal,
    assess_validation_pair,
    attach_formation_execution,
    load_formation_input,
)
from project_model import ProjectModel, URM


def main() -> int:
    checks = 0
    with tempfile.TemporaryDirectory(prefix="wac-formation-urm-") as directory:
        manifest_path, execution_path = write_bundle(Path(directory))

        direct_input = load_formation_input(manifest_path)
        via_urm = URM.formation_input(manifest_path)
        assert via_urm == direct_input
        assert via_urm.certificate()["validation_dataset_eligible"] is True
        checks += 1

        certificate = URM.formation_input_certificate(manifest_path)
        assert certificate == direct_input.certificate()
        assert certificate["scientific_verdict"] == "NONE_NOT_SCORED"
        assert certificate["scientific_validation_authorized"] is False
        checks += 1

        direct_execution = attach_formation_execution(direct_input, execution_path)
        via_execution = URM.formation_execution(manifest_path, execution_path)
        assert via_execution == direct_execution
        assert via_execution.certificate()["generic_predicate_status"] == "ALL_FROZEN_GENERIC_PREDICATES_PASS"
        assert via_execution.certificate()["record_formation_claim_authorized"] is False
        checks += 1

        assessment = URM.formation_validation(manifest_path)
        assert assessment == assess_validation_pair(direct_input)
        assert assessment["validation_dataset_eligible"] is True
        assert assessment["scientific_validation_authorized"] is False
        assert assessment["product_reproduction_attested"] is False
        checks += 1

        try:
            URM.formation_input(Path(directory) / "absent" / "formation-input.json")
        except FormationRefusal as exc:
            assert "manifest is not readable JSON" in str(exc)
        else:
            raise AssertionError("URM accepted a nonexistent formation bundle")
        checks += 1

        method_names = {
            "formation_input",
            "formation_input_certificate",
            "formation_execution",
            "formation_validation",
        }
        assert not method_names & set(ProjectModel.__dict__)
        source_path = Path(__file__).resolve().parent / "formation_input.py"
        tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
        assert not any(
            isinstance(node, (ast.Name, ast.Attribute))
            and getattr(node, "id", getattr(node, "attr", None)) == "RecordSurface"
            for node in ast.walk(tree)
        )
        checks += 1

    assert checks == 6, checks
    print(f"FORMATION_INPUT_URM_CHECKS: {checks}/6 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
