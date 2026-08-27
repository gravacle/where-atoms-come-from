#!/usr/bin/env python3
"""Focused public-door checks for URM world observations.

These checks establish only that the public URM faithfully delegates to the closed
measurement contract.  They cannot authorize a scientific or universal claim.
"""

from __future__ import annotations

from pathlib import Path

from project_model import ProjectModel, URM
from world_observation import ObservationRefusal, load_world_observation


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ACTUAL_MANIFEST = ROOT / "LANE_T53_A_WORLD_WRITE" / "world_observation.json"


def main() -> int:
    checks = 0

    via_urm = URM.world_observation(ACTUAL_MANIFEST)
    direct = load_world_observation(ACTUAL_MANIFEST)
    assert via_urm == direct
    assert len(via_urm.rows) == 760
    checks += 1

    certificate = URM.world_observation_certificate(ACTUAL_MANIFEST)
    assert certificate == direct.certificate()
    assert certificate["evidence_class"] == "ACTUAL_SURFACE_MEASUREMENT"
    assert certificate["scope_classification"] == "CONFIGURATION_EVIDENCE_ONLY"
    checks += 1

    assert certificate["scientific_verdict"] == "NONE_NOT_SCORED"
    assert certificate["independent_reproduction_attested"] is False
    assert certificate["record_formation_proof_authorized"] is False
    assert certificate["universal_claim_authorized"] is False
    assert certificate["public_urm_registration_authorized"] is False
    checks += 1

    try:
        URM.world_observation(ROOT / "does-not-exist" / "manifest.json")
    except ObservationRefusal as exc:
        assert "manifest is not readable JSON" in str(exc)
    else:
        raise AssertionError("URM accepted a nonexistent observation manifest")
    checks += 1

    assert "world_observation" not in ProjectModel.__dict__
    checks += 1

    print(f"URM_WORLD_OBSERVATION_CHECKS: {checks}/5 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
