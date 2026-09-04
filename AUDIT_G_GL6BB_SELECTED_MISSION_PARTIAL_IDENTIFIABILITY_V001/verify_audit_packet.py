#!/usr/bin/env python3
"""Fail-closed custody and scope verifier for the GL6BB hostile audit."""

from __future__ import annotations

import hashlib
import re
import unicodedata
from pathlib import Path, PurePosixPath


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6BB_SELECTED_MISSION_PARTIAL_IDENTIFIABILITY_V001"
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS}] {label}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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
        pure = PurePosixPath(relative)
        check(not pure.is_absolute() and all(part not in ("", ".", "..")
                                             for part in pure.parts)
              and "\\" not in relative,
              f"confined POSIX path {path.name}:{lineno}")
        check(relative not in output, f"unique path {path.name}:{lineno}")
        output[relative] = digest
    return output


REQUIRED = {
    "README.md", "AUDIT.md", "AUDITED_TARGETS.sha256", "VERIFICATION.txt",
    "independent_gl6bb_replay.py", "verify_audit_packet.py",
    "MANIFEST.sha256", "SEAL.sha256",
}

AUTHOR_NAMES = {
    "DEPENDENCIES.md", "DEPENDENCIES.sha256", "IDENTIFIABILITY_LEDGER.json",
    "MANIFEST.sha256", "README.md", "RESULT.md", "SEAL.sha256",
    "SELF_AUDIT.md", "THEOREM.md", "VERIFICATION.txt",
    "calculate_prepared_blank_collar0.py", "verify_packet.py",
    "verify_selected_mission_partial_identifiability.py",
}

DEPENDENCY_TARGETS = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md",
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MANIFEST.sha256",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/THEOREM.md",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/MANIFEST.sha256",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/INDEPENDENT_REAUDIT.md",
    "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/MANIFEST.sha256",
    "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/SEAL.sha256",
    "AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/AUDIT.md",
    "AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/MANIFEST.sha256",
    "AUDIT_G_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001/SEAL.sha256",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/MANIFEST.sha256",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/SEAL.sha256",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md",
    "LANE_CROSS_RFT_GRA_GL5ZZF_H6_DOUBLED_REAL_METRIC_SOURCE_BIND_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/AUDIT_MANIFEST.sha256",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/AUDIT.md",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/THEOREM.md",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/MANIFEST.sha256",
    "LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/SEAL.sha256",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/AUDIT.md",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/MANIFEST.sha256",
    "AUDIT_G_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/SEAL.sha256",
}


def resolve_ledger_path(parent: Path, relative: str) -> Path:
    root_candidate = ROOT / relative
    if root_candidate.exists():
        return root_candidate
    return parent / relative


def verify_frozen_author() -> None:
    check({entry.name for entry in AUTHOR.iterdir()} == AUTHOR_NAMES,
          "author directory has the exact frozen file set")
    check(all((AUTHOR / name).is_file() and not (AUTHOR / name).is_symlink()
              for name in AUTHOR_NAMES),
          "every author object is a regular nonsymlink file")

    targets = rows(HERE / "AUDITED_TARGETS.sha256")
    check(len(targets) == 13, "thirteen frozen author targets")
    check({Path(relative).name for relative in targets} == AUTHOR_NAMES,
          "author target names are exact")
    for relative, expected in sorted(targets.items()):
        path = ROOT / relative
        check(path.is_file() and not path.is_symlink(),
              f"regular author target: {relative}")
        check(path.parent == AUTHOR, f"author target confined: {relative}")
        check(sha256(path) == expected, f"author target hash: {relative}")

    prefix = str(AUTHOR.relative_to(ROOT)) + "/"
    check(targets[prefix + "THEOREM.md"]
          == "ed67f9d3bcd972fd4298c927280bbb522b9e129c5b0c15166a18c7e640ef6f88",
          "frozen theorem hash is exact")
    check(targets[prefix + "MANIFEST.sha256"]
          == "8ba5883baea82eff8e758f059e0db35011bbe5eb2c1bcae5486090d3237e53c5",
          "frozen author-manifest hash is exact")
    check(targets[prefix + "SEAL.sha256"]
          == "b6094f1a8569cc960afc10f1e2a608c194430d4c6cf0797a64a2498f0e1f0fe8",
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


def verify_dependencies() -> list[Path]:
    dependencies = rows(AUTHOR / "DEPENDENCIES.sha256")
    check(len(dependencies) == 28, "twenty-eight frozen dependency targets")
    check(set(dependencies) == DEPENDENCY_TARGETS,
          "frozen dependency path set is exact")
    paths = []
    for relative, expected in sorted(dependencies.items()):
        path = ROOT / relative
        paths.append(path)
        check(path.is_file() and not path.is_symlink(),
              f"regular dependency: {relative}")
        check(sha256(path) == expected, f"dependency hash: {relative}")

    # Check that every pinned manifest internally points to every pinned peer it
    # covers, irrespective of whether that manifest uses basename or root paths.
    dependency_paths = {str((ROOT / relative).resolve()): digest
                        for relative, digest in dependencies.items()}
    for relative in sorted(dependencies):
        if not relative.endswith(("MANIFEST.sha256", "AUDIT_MANIFEST.sha256")):
            continue
        manifest_path = ROOT / relative
        manifest_rows = rows(manifest_path)
        covered = 0
        for child_relative, child_digest in manifest_rows.items():
            child = resolve_ledger_path(manifest_path.parent, child_relative).resolve()
            key = str(child)
            if key in dependency_paths:
                covered += 1
                check(child_digest == dependency_paths[key],
                      f"dependency manifest agrees: {relative}:{child_relative}")
        check(covered >= 1, f"dependency manifest covers at least one pinned peer: {relative}")

    for relative in sorted(dependencies):
        if not relative.endswith("SEAL.sha256"):
            continue
        seal_path = ROOT / relative
        seal_rows = rows(seal_path)
        check(len(seal_rows) == 1, f"dependency seal has one row: {relative}")
        target_relative, expected = next(iter(seal_rows.items()))
        target = resolve_ledger_path(seal_path.parent, target_relative)
        check(target.is_file() and target.name == "MANIFEST.sha256",
              f"dependency seal targets a manifest: {relative}")
        check(expected == sha256(target), f"dependency seal resolves: {relative}")
    return paths


def verify_no_hidden_controls(paths: list[Path], label: str) -> None:
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="strict")
        forbidden = [character for character in text
                     if unicodedata.category(character).startswith("C")
                     and character not in "\n\t"]
        check(not forbidden, f"no hidden Unicode/control code points in {label}: {path.name}")
        check("\ufeff" not in text, f"no byte-order mark in {label}: {path.name}")


def verify_audit_custody() -> None:
    check({entry.name for entry in HERE.iterdir()} == REQUIRED,
          "audit directory has the exact sealed file set")
    check(all((HERE / name).is_file() and not (HERE / name).is_symlink()
              for name in REQUIRED),
          "every audit object is a regular nonsymlink file")
    manifest = rows(HERE / "MANIFEST.sha256")
    check(set(manifest) == REQUIRED - {"MANIFEST.sha256", "SEAL.sha256"},
          "audit manifest covers the exact pre-manifest packet")
    for relative, expected in sorted(manifest.items()):
        check(sha256(HERE / relative) == expected,
              f"audit manifest resolves: {relative}")
    seal = rows(HERE / "SEAL.sha256")
    check(seal == {"MANIFEST.sha256": sha256(HERE / "MANIFEST.sha256")},
          "audit seal resolves")


def verify_scope() -> None:
    audit = " ".join((HERE / "AUDIT.md").read_text(encoding="utf-8").split())
    verification = " ".join((HERE / "VERIFICATION.txt").read_text(
        encoding="utf-8").split())
    required_audit = (
        "Hostile verdict: PASS",
        "All three entries remain unselected",
        "exactly `[0,1]`",
        "not an authenticated-preparation claim",
        "pointwise correlation",
        "adds one `eta`, not two",
        "five-dimensional cyclic subspace",
        "(1,1/2,1/3,1/2,1)",
        "R<G>=-<sum_a X_a><=4",
        "5/6",
        "13/15",
        "11/6-e^{96|\\sigma_{\\rm obs}|}",
        "28/15-e^{120|\\sigma_{\\rm obs}|}",
        "No read outcome or retained flag is postselected",
        "sampling endpoint is immediately before",
        "only inside the explicitly prepared-blank two-member scenario",
        "single actual mission still needs",
        "fail closed",
        "No graviton, Ricci/Einstein premise, gravity identification, or `G`",
    )
    for token in required_audit:
        check(token in audit, f"audit scope token: {token}")

    required_verification = (
        "PASS: 198747 exact GL6BB checks",
        "PASS: 221/221 GL6BB packet checks",
        "PASS__INDEPENDENT_GL6BB_HOSTILE_REPLAY__42886/42886",
        "PASS__GL6BB_HOSTILE_AUDIT_PACKET__PASS",
    )
    for token in required_verification:
        check(token in verification, f"verification token: {token}")

    forbidden = (
        "we have proved gravity", "this is gravity", "derives newton's constant",
        "selects r=2 as nature's value", "selects r=5/2 as nature's value",
        "every endpoint state has an authenticated preparation",
        "one fixed collar controls unbounded time",
        "sigma is the sole datum for an actual selected mission",
    )
    lower = audit.lower()
    for phrase in forbidden:
        check(phrase not in lower, f"forbidden promotion absent: {phrase}")


def main() -> int:
    verify_frozen_author()
    dependency_paths = verify_dependencies()
    verify_audit_custody()
    verify_no_hidden_controls([AUTHOR / name for name in sorted(AUTHOR_NAMES)], "author")
    verify_no_hidden_controls(dependency_paths, "dependency")
    verify_no_hidden_controls([HERE / name for name in sorted(REQUIRED)], "audit")
    verify_scope()
    print("PASS__GL6BB_HOSTILE_AUDIT_PACKET__PASS")
    print(f"CHECKS {CHECKS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
