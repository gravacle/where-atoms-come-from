#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6BA hostile audit."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001"
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS + 1}] {label}")
    CHECKS += 1


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        parts = raw.split("  ", 1)
        check(len(parts) == 2, f"two-column hash row {path.name}:{lineno}")
        digest, relative = parts
        check(bool(re.fullmatch(r"[0-9a-f]{64}", digest)),
              f"digest syntax {path.name}:{lineno}")
        check(relative not in output, f"unique path {path.name}:{lineno}")
        output[relative] = digest
    return output


REQUIRED = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256", "VERIFICATION.txt",
    "independent_gl6ba_replay.py", "verify_audit_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}

AUTHOR_NAMES = {
    "COLLAR_LEDGER.json", "DEPENDENCIES.md", "DEPENDENCIES.sha256",
    "MANIFEST.sha256", "README.md", "RESULT.md", "SEAL.sha256",
    "SELF_AUDIT.md", "THEOREM.md", "VERIFICATION.txt",
    "verify_finite_mission_collar.py", "verify_packet.py",
}

DEPENDENCY_TARGETS = {
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/AUDIT.md",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/POSTFREEZE_AUDIT.md",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/MANIFEST.sha256",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/SEAL.sha256",
    "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/MANIFEST.sha256",
    "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/SEAL.sha256",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256",
}


def verify_frozen_author() -> None:
    targets = rows(HERE / "AUDITED_TARGETS.sha256")
    check(len(targets) == 12, "twelve frozen author targets")
    check({Path(relative).name for relative in targets} == AUTHOR_NAMES,
          "author target names are exact")
    for relative, expected in sorted(targets.items()):
        path = ROOT / relative
        check(path.is_file() and not path.is_symlink(), f"regular author target: {relative}")
        check(path.parent == AUTHOR, f"author target confined: {relative}")
        check(sha256(path) == expected, f"author target hash: {relative}")

    check(targets[str(AUTHOR.relative_to(ROOT) / "THEOREM.md")]
          == "d7ce0a7527a68f49e6ea2ee8edbb400a142fbb49297d8fe99cae78ffa0154ab0",
          "frozen theorem hash is exact")
    check(targets[str(AUTHOR.relative_to(ROOT) / "MANIFEST.sha256")]
          == "6e14332230f713d51e393a5889fe78964fe0e63588b4b841533fa6af7ef19103",
          "frozen author-manifest hash is exact")
    check(targets[str(AUTHOR.relative_to(ROOT) / "SEAL.sha256")]
          == "34f29b3c03d53c4dbc9736d1bf7a7785e0a49a2aad04299b69a3804290c5971e",
          "frozen author-seal-file hash is exact")

    author_manifest = rows(AUTHOR / "MANIFEST.sha256")
    check(set(author_manifest) == AUTHOR_NAMES - {"MANIFEST.sha256", "SEAL.sha256"},
          "author manifest covers the exact pre-manifest packet")
    for relative, expected in sorted(author_manifest.items()):
        check(sha256(AUTHOR / relative) == expected,
              f"author manifest resolves: {relative}")
    author_seal = rows(AUTHOR / "SEAL.sha256")
    check(author_seal == {"MANIFEST.sha256": sha256(AUTHOR / "MANIFEST.sha256")},
          "author seal resolves")


def verify_dependencies() -> None:
    dependencies = rows(AUTHOR / "DEPENDENCIES.sha256")
    check(len(dependencies) == 12, "twelve frozen dependency targets")
    check(set(dependencies) == DEPENDENCY_TARGETS,
          "frozen dependency path set is exact")
    for relative, expected in sorted(dependencies.items()):
        path = ROOT / relative
        check(path.is_file() and not path.is_symlink(), f"regular dependency: {relative}")
        check(sha256(path) == expected, f"dependency hash: {relative}")
    for marker in ("GL6AK", "GL6AZ"):
        check(any(marker in relative for relative in dependencies),
              f"dependency family present: {marker}")


def verify_audit_custody() -> None:
    check({path.name for path in HERE.iterdir() if path.is_file()} == REQUIRED,
          "audit file set is exact")
    manifest = rows(HERE / "MANIFEST.sha256")
    check(set(manifest) == REQUIRED - {"MANIFEST.sha256", "SEAL.sha256"},
          "audit manifest covers the exact pre-manifest packet")
    for relative, expected in sorted(manifest.items()):
        path = HERE / relative
        check(path.is_file() and not path.is_symlink(), f"regular audit target: {relative}")
        check(sha256(path) == expected, f"audit manifest resolves: {relative}")
    seal = rows(HERE / "SEAL.sha256")
    check(seal == {"MANIFEST.sha256": sha256(HERE / "MANIFEST.sha256")},
          "audit seal resolves")


def verify_scope() -> None:
    audit = " ".join((HERE / "AUDIT.md").read_text(encoding="utf-8").split())
    verification = " ".join((HERE / "VERIFICATION.txt").read_text(encoding="utf-8").split())
    required_audit_tokens = (
        "Hostile verdict: PASS",
        "complete finite all-formed/`MATCH` FPSS exteriors",
        "-3R-3R=-6R",
        "C_L=12(3L^2+3L+1)",
        "d_L(z,p_{\\rm inside})\\ge2L",
        "={1\\over12}",
        "T_{2L+1}(48R|s|)",
        "Taylor coefficients match through order `4L+1`",
        "order `4L+2`",
        "real two-root `M_beta`",
        "not promoted to a claim that the fully summed coefficient is nonzero",
        "At `L=0`",
        "At `R=0`",
        "exact reduction of that same state",
        "The `1/2` is exact",
        "No read outcome or flag is postselected",
        "random route mixture",
        "mathematical completion",
        "96|\\sigma_{\\rm obs}|",
        "120|\\sigma_{\\rm obs}|",
        "uses ADHH/`GL6AY`",
        "No conventional gravity premise",
    )
    for token in required_audit_tokens:
        check(token in audit, f"audit scope token: {token}")

    required_verification_tokens = (
        "PASS: 636296/636296 GL6BA constructive checks",
        "PASS: 135/135 GL6BA packet checks",
        "PASS__INDEPENDENT_GL6BA_HOSTILE_REPLAY__2008230/2008230",
        "PASS__GL6BA_HOSTILE_AUDIT_PACKET__",
    )
    for token in required_verification_tokens:
        check(token in verification, f"verification token: {token}")

    lower = audit.lower()
    forbidden = (
        "we have proved gravity",
        "this is gravity",
        "derives newton's constant",
        "selects r=2 as nature's value",
        "selects r=5/2 as nature's value",
        "controls every random route mixture",
        "bounds the full retained output distribution",
        "authenticates one infinite record",
        "one fixed collar controls infinite time",
    )
    for phrase in forbidden:
        check(phrase not in lower, f"forbidden promotion absent: {phrase}")


def main() -> None:
    verify_frozen_author()
    verify_dependencies()
    verify_audit_custody()
    verify_scope()
    print(f"PASS__GL6BA_HOSTILE_AUDIT_PACKET__{CHECKS}/{CHECKS}")


if __name__ == "__main__":
    main()
