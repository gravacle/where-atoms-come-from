#!/usr/bin/env python3
"""Fail-closed dependency, manifest, seal, and scope verifier for GL6BC."""

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
    "RECIPROCITY_LEDGER.json",
    "verify_finite_collar_reciprocal_lineage_obstruction.py",
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
    check(len(rows) == 13, "exact dependency row count")
    for rel, expected in sorted(rows.items()):
        path = ROOT / rel
        check(path.is_file(), f"dependency exists: {rel}")
        check(sha256(path) == expected, f"dependency hash: {rel}")


def verify_manifest() -> None:
    rows = parse_hash_ledger(LANE / "MANIFEST.sha256")
    expected = REQUIRED - {"MANIFEST.sha256", "SEAL.sha256"}
    check(set(rows) == expected, "manifest covers exact pre-manifest packet")
    for rel, expected in sorted(rows.items()):
        path = LANE / rel
        check(path.is_file(), f"manifest target exists: {rel}")
        check(sha256(path) == expected, f"manifest target hash: {rel}")


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
        "d_L>\\delta_L",
        "one-tail",
        "p^{\\Omega,0}=p^{L,0}",
        "\\chi^R_{\\Pi_\\beta,B}(s)=0",
        "\\chi^R_{B,\\Pi_\\beta}(s)=0",
        "future writer/formation channel",
        "undefined",
        "identity",
        "not a no-go",
        "No graviton",
        "no gravity",
        "no G",
    ]
    for token in required:
        check(token in text, f"packet scope token: {token}")
    verification = (LANE / "VERIFICATION.txt").read_text(encoding="utf-8")
    check("PASS: 29552/29552" in verification,
          "verification transcript records constructive pass")
    check("distinct hostile audit remains required" in verification,
          "verification preserves audit ceiling")


def main() -> None:
    verify_dependencies()
    verify_manifest()
    verify_seal()
    verify_scope()
    print(f"PASS: {CHECKS}/{CHECKS} GL6BC packet checks")


if __name__ == "__main__":
    main()

