#!/usr/bin/env python3
"""Fail-closed custody, manifest, seal, and scope verifier for GL6BA."""

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
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


REQUIRED = {
    "README.md",
    "THEOREM.md",
    "RESULT.md",
    "SELF_AUDIT.md",
    "DEPENDENCIES.md",
    "DEPENDENCIES.sha256",
    "COLLAR_LEDGER.json",
    "verify_finite_mission_collar.py",
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
        digest, rel = parts
        check(bool(re.fullmatch(r"[0-9a-f]{64}", digest)),
              f"valid digest {path.name}:{lineno}")
        check(rel not in rows, f"unique path {path.name}:{lineno}")
        rows[rel] = digest
    return rows


def verify_dependencies() -> None:
    rows = parse_hash_ledger(LANE / "DEPENDENCIES.sha256")
    check(len(rows) == 12, "exact dependency row count")
    for rel, expected in sorted(rows.items()):
        path = ROOT / rel
        check(path.is_file(), f"dependency exists: {rel}")
        check(sha256(path) == expected, f"dependency hash: {rel}")


def verify_manifest() -> None:
    rows = parse_hash_ledger(LANE / "MANIFEST.sha256")
    expected = REQUIRED - {"MANIFEST.sha256", "SEAL.sha256"}
    check(set(rows) == expected, "manifest covers exact pre-manifest packet")
    for rel, digest in sorted(rows.items()):
        path = LANE / rel
        check(path.is_file(), f"manifest target exists: {rel}")
        check(sha256(path) == digest, f"manifest target hash: {rel}")


def verify_seal() -> None:
    rows = parse_hash_ledger(LANE / "SEAL.sha256")
    check(set(rows) == {"MANIFEST.sha256"}, "seal covers manifest only")
    check(rows["MANIFEST.sha256"] == sha256(LANE / "MANIFEST.sha256"),
          "seal pins exact manifest")


def verify_scope() -> None:
    names = {path.name for path in LANE.iterdir() if path.is_file()}
    check(names == REQUIRED, "no missing or unsealed extra packet files")
    text = "\n".join(
        (LANE / name).read_text(encoding="utf-8", errors="strict")
        for name in sorted(REQUIRED)
        if name.endswith((".md", ".json", ".txt"))
    )
    required = [
        "full F3 dynamics versus a spatial collar of the same",
        "C_L=36L^2+36L+12",
        "T_{2L+1}(48R|s|)",
        "ordinary nested-commutator order `4L+2`",
        "D_{\\rm TV}(p^\\Omega,p^{(L)})",
        "complete finite authenticated exterior",
        "complete finite all-formed/`MATCH` FPSS",
        "at most the `C_L` crossing terms",
        "mathematical completion",
        "not a claim that one infinite record",
        "selected-factor binary pair marginal",
        "96|\\sigma_{\\rm obs}|",
        "120|\\sigma_{\\rm obs}|",
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
    print(f"PASS: {CHECKS}/{CHECKS} GL6BA packet checks")


if __name__ == "__main__":
    main()
