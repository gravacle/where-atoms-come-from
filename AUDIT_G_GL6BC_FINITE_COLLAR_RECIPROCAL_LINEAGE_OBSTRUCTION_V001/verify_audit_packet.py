#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6BC hostile audit."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6BC_FINITE_COLLAR_RECIPROCAL_LINEAGE_OBSTRUCTION_V001"
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS + 1}] {label}")
    CHECKS += 1


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        pieces = raw.split("  ", 1)
        check(len(pieces) == 2, f"two-column hash row {path.name}:{line_number}")
        digest, relative = pieces
        check(bool(re.fullmatch(r"[0-9a-f]{64}", digest)),
              f"digest syntax {path.name}:{line_number}")
        check(relative not in output, f"unique path {path.name}:{line_number}")
        output[relative] = digest
    return output


REQUIRED = {
    "AUDIT.md",
    "AUDITED_TARGETS.sha256",
    "DEPENDENCIES.sha256",
    "README.md",
    "VERIFICATION.txt",
    "independent_gl6bc_replay.py",
    "verify_audit_packet.py",
    "MANIFEST.sha256",
    "SEAL.sha256",
}

AUTHOR_NAMES = {
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "MANIFEST.sha256",
    "README.md",
    "RECIPROCITY_LEDGER.json",
    "RESULT.md",
    "SEAL.sha256",
    "SELF_AUDIT.md",
    "THEOREM.md",
    "VERIFICATION.txt",
    "verify_finite_collar_reciprocal_lineage_obstruction.py",
    "verify_packet.py",
}

DEPENDENCY_TARGETS = {
    "GRAVITY_RECORD_FIRST_WORKING_THEORY_CLOSURE_V001.md",
    "GRAVITY_RGRL_LOCAL_RESPONSE_G_IDENTIFIABILITY_CEILING_V002.md",
    "GRAVITY_RGRL_ONSHELL_OFFSHELL_CLARIFICATION_ADOPTION_V001.md",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/AUDIT.md",
    "AUDIT_G_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/AUDIT.md",
    "AUDIT_G_GL6T_F3_LINEAGE_GATED_Q4_E2_RESPONSE_V001/AUDIT.md",
    "AUDIT_G_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/AUDIT.md",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6T_F3_LINEAGE_GATED_Q4_E2_RESPONSE_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001/THEOREM.md",
    "LANE_GRA_DQ_F3_RECIPROCAL_BACKPROPAGATION_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md",
    "LANE_GRA_DQ_F3_RECIPROCAL_BACKPROPAGATION_SCREEN_V001/THEOREM.md",
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

    author_prefix = str(AUTHOR.relative_to(ROOT)) + "/"
    check(targets[author_prefix + "THEOREM.md"]
          == "6355e6e1dda470e363122f5e3342c01346dc81e8992d05b1436f30b899041ea6",
          "frozen theorem hash is exact")
    check(targets[author_prefix + "MANIFEST.sha256"]
          == "614570080fd9ce3ebac3edfef4c74ede52766c96948659eba643e2b3b286ab5a",
          "frozen author-manifest hash is exact")
    check(targets[author_prefix + "SEAL.sha256"]
          == "c6ccb51c0796ef889b471fde9c691407cdda16d23239deb917cc4d6b0405cd03",
          "frozen author-seal-file hash is exact")

    manifest = rows(AUTHOR / "MANIFEST.sha256")
    check(set(manifest) == AUTHOR_NAMES - {"MANIFEST.sha256", "SEAL.sha256"},
          "author manifest covers the exact pre-manifest packet")
    for relative, expected in sorted(manifest.items()):
        check(sha256(AUTHOR / relative) == expected,
              f"author manifest resolves: {relative}")
    seal = rows(AUTHOR / "SEAL.sha256")
    check(seal == {"MANIFEST.sha256": sha256(AUTHOR / "MANIFEST.sha256")},
          "author seal resolves")


def verify_dependencies() -> None:
    audit_dependencies = rows(HERE / "DEPENDENCIES.sha256")
    author_dependencies = rows(AUTHOR / "DEPENDENCIES.sha256")
    check(len(audit_dependencies) == 13, "thirteen frozen dependency targets")
    check(set(audit_dependencies) == DEPENDENCY_TARGETS,
          "dependency path set is exact")
    check(audit_dependencies == author_dependencies,
          "audit dependency ledger exactly reproduces author dependency ledger")
    for relative, expected in sorted(audit_dependencies.items()):
        path = ROOT / relative
        check(path.is_file() and not path.is_symlink(), f"regular dependency: {relative}")
        check(sha256(path) == expected, f"dependency hash: {relative}")
    for marker in ("GL6BA", "GL6AI", "GL6T", "GL6V", "GRA_DQ"):
        normalized_marker = marker.replace("_", "-") if marker == "GRA_DQ" else marker
        if marker == "GRA_DQ":
            check(any("GRA_DQ" in relative for relative in audit_dependencies),
                  "dependency family present: GRA-DQ")
        else:
            check(any(normalized_marker in relative for relative in audit_dependencies),
                  f"dependency family present: {marker}")


def verify_audit_custody() -> None:
    check({path.name for path in HERE.iterdir() if path.is_file()} == REQUIRED,
          "audit file inventory is exact")
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
        "No author byte was edited",
        "whole-exterior all-`KEEP`/all-`BREAK` product is an explicit finite tensor-product premise",
        "p^{\\Omega,0}=p^{L,0}",
        "violates a one-tail lower bound",
        "[U[j],\\Pi_\\beta]=0",
        "Both retarded orientations therefore vanish as operator identities",
        "A route operation or writer is not an allowed `B`",
        "not the full reduced `K` density matrix",
        "Bayesian selection, not a K-changing channel",
        "bounded local cylinder observable",
        "does not create or measure one global infinite support-word projector",
        "not a selected positive number for `R=2` or `R=5/2`",
        "minimal missing object",
        "does **not** prove a unique minimal Hamiltonian term",
        "undefined",
        "no reciprocal future-record back-reaction datum",
        "no selected moderate-`R` signal",
        "Newton constant `G`",
    )
    for token in required_audit_tokens:
        check(token in audit, f"audit scope token: {token}")

    required_verification_tokens = (
        "PASS: 29552/29552 GL6BC constructive checks",
        "PASS: 137/137 GL6BC packet checks",
        "PASS__INDEPENDENT_GL6BC_HOSTILE_REPLAY__60118/60118",
        "PASS__INDEPENDENT_GL6BA_HOSTILE_REPLAY__2008230/2008230",
        "PASS__INDEPENDENT_GL6AI_REPLAY__273976/273976",
        "PASS__INDEPENDENT_GL6T_REPLAY__178/178",
        "PASS__INDEPENDENT_GL6V_REPLAY__70946/70946",
        "PASS__GL6BC_HOSTILE_AUDIT_PACKET__",
    )
    for token in required_verification_tokens:
        check(token in verification, f"verification token: {token}")

    lower = audit.lower()
    forbidden = (
        "we have proved gravity",
        "this is gravity",
        "derives newton's constant",
        "all postselected route probabilities are conserved",
        "the entire k density matrix is conserved",
        "one infinite route word is measured",
        "unique minimal hamiltonian writer",
        "future-writer response is nonzero",
        "selects r=2",
        "selects r=5/2",
    )
    for phrase in forbidden:
        check(phrase not in lower, f"forbidden promotion absent: {phrase}")


def main() -> None:
    verify_frozen_author()
    verify_dependencies()
    verify_audit_custody()
    verify_scope()
    print(f"PASS__GL6BC_HOSTILE_AUDIT_PACKET__{CHECKS}/{CHECKS}")


if __name__ == "__main__":
    main()
