#!/usr/bin/env python3
"""Fail-closed custody, manifest, seal, and scope verifier for GL6BB."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
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


REQUIRED = {
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "IDENTIFIABILITY_LEDGER.json",
    "calculate_prepared_blank_collar0.py",
    "verify_selected_mission_partial_identifiability.py",
    "verify_packet.py",
    "VERIFICATION.txt",
    "MANIFEST.sha256",
    "SEAL.sha256",
}


def parse_hash_ledger(path: Path) -> dict[str, str]:
    rows: dict[str, str] = {}
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        parts = raw.split("  ", 1)
        check(len(parts) == 2, f"two-column hash row {path.name}:{lineno}")
        digest, relative = parts
        check(bool(re.fullmatch(r"[0-9a-f]{64}", digest)),
              f"valid digest {path.name}:{lineno}")
        check(relative not in rows, f"unique path {path.name}:{lineno}")
        rows[relative] = digest
    return rows


def verify_dependencies() -> None:
    rows = parse_hash_ledger(LANE / "DEPENDENCIES.sha256")
    check(len(rows) == 28, "exact dependency row count")
    for relative, expected in sorted(rows.items()):
        path = ROOT / relative
        check(path.is_file(), f"dependency exists: {relative}")
        check(sha256(path) == expected, f"dependency hash: {relative}")


def verify_manifest() -> None:
    rows = parse_hash_ledger(LANE / "MANIFEST.sha256")
    expected = REQUIRED - {"MANIFEST.sha256", "SEAL.sha256"}
    check(set(rows) == expected, "manifest covers exact pre-manifest packet")
    for relative, digest in sorted(rows.items()):
        path = LANE / relative
        check(path.is_file(), f"manifest target exists: {relative}")
        check(sha256(path) == digest, f"manifest target hash: {relative}")


def verify_seal() -> None:
    rows = parse_hash_ledger(LANE / "SEAL.sha256")
    check(set(rows) == {"MANIFEST.sha256"}, "seal covers manifest only")
    check(rows["MANIFEST.sha256"] == sha256(LANE / "MANIFEST.sha256"),
          "seal pins exact manifest")


def verify_scope() -> None:
    names = {path.name for path in LANE.iterdir() if path.is_file()}
    check(names == REQUIRED, "no missing or unsealed extra packet files")
    check(all(
        all(byte >= 32 or byte in (9, 10) for byte in (LANE / name).read_bytes())
        for name in REQUIRED
    ), "packet contains no hidden control characters")
    text = "\n".join(
        (LANE / name).read_text(encoding="utf-8", errors="strict")
        for name in sorted(REQUIRED)
        if name.endswith((".md", ".json", ".txt"))
    )
    required = [
        "does not select any one of the three inputs",
        "not a selection",
        "no finite upper bound",
        "prepared-blank state",
        "=[0,1]",
        "pointwise robust",
        "trace-distance",
        "five-state Dicke Hamiltonian",
        "1-{1\\over3R}",
        "11/6-e^{96|\\sigma_{\\rm obs}|}",
        "28/15-e^{120|\\sigma_{\\rm obs}|}",
        "one authenticated tuple",
        "required command-line input and has no physical default",
        "No graviton",
        "no gravity",
        "no G",
    ]
    for token in required:
        check(token in text, f"packet scope token: {token}")
    check("PASS:" in (LANE / "VERIFICATION.txt").read_text(encoding="utf-8"),
          "verification transcript records constructive pass")


def main() -> None:
    verify_dependencies()
    verify_manifest()
    verify_seal()
    verify_scope()
    print(f"PASS: {CHECKS}/{CHECKS} GL6BB packet checks")


if __name__ == "__main__":
    main()
