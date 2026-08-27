#!/usr/bin/env python3
"""Six fixed checks for the public URM-only gamma-flow delegates."""

from __future__ import annotations

import ast
from pathlib import Path
import tempfile

from checks_gamma_flow import build_manifest, write_manifest
from gamma_flow import GammaFlowRefusal, load_gamma_flow
from project_model import ProjectModel, URM


def main() -> int:
    checks = 0
    with tempfile.TemporaryDirectory(prefix="wac-gamma-flow-urm-") as directory:
        root = Path(directory)
        manifest_path = write_manifest(root, build_manifest(), "gamma-flow.json")

        direct = load_gamma_flow(manifest_path)
        via_urm = URM.gamma_flow(manifest_path)
        assert via_urm == direct
        assert via_urm.certificate()["scientific_result"] == "NO_PROOF_OUTPUT"
        checks += 1

        certificate = URM.gamma_flow_certificate(manifest_path)
        assert certificate == direct.certificate()
        assert certificate["actual_platform_present"] is False
        assert certificate["synthetic_fixture_scientific_weight"] == "ZERO"
        checks += 1

        assert URM.gamma_flow_states(manifest_path) == dict(
            direct.internal_discovery_states
        )
        assert set(URM.gamma_flow_states(manifest_path).values()) <= {
            "PASS", "FAIL", "UNCLASSIFIED", "UNSCOREABLE"
        }
        checks += 1

        outputs = URM.gamma_flow_proof_outputs(manifest_path)
        assert outputs == dict(direct.proof_outputs)
        assert set(outputs.values()) == {"NO_PROOF_OUTPUT"}
        checks += 1

        try:
            URM.gamma_flow(root / "absent" / "gamma-flow.json")
        except GammaFlowRefusal as exc:
            assert "manifest is not readable JSON" in str(exc)
        else:
            raise AssertionError("URM accepted a nonexistent gamma-flow envelope")
        checks += 1

        method_names = {
            "gamma_flow",
            "gamma_flow_certificate",
            "gamma_flow_states",
            "gamma_flow_proof_outputs",
        }
        assert not method_names & set(ProjectModel.__dict__)
        assert method_names <= set(URM.__dict__)
        assert {
            "world_observation",
            "formation_input",
        } <= set(URM.__dict__)
        source_path = Path(__file__).resolve().parent / "gamma_flow.py"
        tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
        names = {
            getattr(node, "id", getattr(node, "attr", None))
            for node in ast.walk(tree)
            if isinstance(node, (ast.Name, ast.Attribute))
        }
        imports = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        assert "RecordSurface" not in names
        assert "ProjectModel" not in names
        assert "importlib" not in imports
        assert not {"__import__", "eval", "exec"} & names
        checks += 1

    assert checks == 6, checks
    print(f"GAMMA_FLOW_URM_CHECKS: {checks}/6 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
