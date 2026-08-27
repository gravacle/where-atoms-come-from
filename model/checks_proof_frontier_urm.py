#!/usr/bin/env python3
"""Six public-URM delegate and isolation checks for the proof frontier."""

from __future__ import annotations

import json
from pathlib import Path
import tempfile

from checks_proof_frontier import build_manifest, write_manifest
from proof_frontier import ProofFrontierRefusal, load_proof_frontier
from project_model import URM


def main() -> int:
    checks = 0
    with tempfile.TemporaryDirectory(prefix="wac-proof-frontier-urm-") as directory:
        root = Path(directory) / "bundle"
        root.mkdir()
        manifest_path = write_manifest(root, build_manifest(root))

        direct = load_proof_frontier(manifest_path)
        delegated = URM.proof_frontier(manifest_path)
        assert delegated.manifest_sha256 == direct.manifest_sha256
        checks += 1

        assert URM.proof_frontier_certificate(manifest_path) == direct.certificate()
        checks += 1

        assert URM.proof_frontier_proof_states(manifest_path) == {
            "URF": "BLOCKED_MISSING_DATA",
            "GE": "BLOCKED_MISSING_DATA",
        }
        checks += 1

        assert URM.proof_frontier_theory_states(manifest_path) == {
            "T_RECORD": "PARTIAL_SUPPORT_ONLY"
        }
        assert set(URM.proof_frontier_execution_frontier(manifest_path)) == {"u_cal", "g_cal"}
        checks += 1

        try:
            URM.proof_frontier(root / "absent.json")
        except ProofFrontierRefusal:
            pass
        else:
            raise AssertionError("missing frontier did not refuse")
        checks += 1

        required = {
            "proof_frontier",
            "proof_frontier_certificate",
            "proof_frontier_proof_states",
            "proof_frontier_theory_states",
            "proof_frontier_execution_frontier",
        }
        assert required <= set(dir(URM))
        source = (Path(__file__).resolve().parent / "proof_frontier.py").read_text(encoding="utf-8")
        assert "from formation_input" not in source
        assert "from gamma_flow" not in source
        assert "from world_observation" not in source
        assert json.loads(manifest_path.read_text(encoding="utf-8"))["schema"] == "WAC_PROOF_FRONTIER_V001"
        checks += 1

    assert checks == 6
    print("PROOF_FRONTIER_URM_CHECKS: 6/6 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

