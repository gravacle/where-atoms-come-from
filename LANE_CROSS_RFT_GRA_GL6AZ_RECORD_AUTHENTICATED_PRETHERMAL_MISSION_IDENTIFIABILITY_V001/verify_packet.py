#!/usr/bin/env python3
"""Fail-closed packet, custody, manifest, and seal verifier for GL6AZ."""

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
    "SEARCH_LEDGER.md",
    "IDENTIFIABILITY_LEDGER.json",
    "verify_record_authenticated_prethermal_mission.py",
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
    check(len(rows) == 33, "exact dependency row count")
    for rel, expected in sorted(rows.items()):
        path = ROOT / rel
        check(path.is_file(), f"dependency exists: {rel}")
        check(sha256(path) == expected, f"dependency hash: {rel}")


def verify_manifest() -> None:
    rows = parse_hash_ledger(LANE / "MANIFEST.sha256")
    expected_names = REQUIRED - {"MANIFEST.sha256", "SEAL.sha256"}
    check(set(rows) == expected_names, "manifest covers exact pre-manifest packet")
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
    names = {p.name for p in LANE.iterdir() if p.is_file()}
    check(names == REQUIRED, "no missing or unsealed extra packet files")
    text = "\n".join(
        (LANE / name).read_text(encoding="utf-8", errors="strict")
        for name in sorted(REQUIRED)
        if name.endswith((".md", ".json", ".txt"))
    )
    required_tokens = [
        "theorem-domain failure, not physical-phase failure",
        "R={U_d\\over h}",
        "D_{\\rm TV}(p^H,p^{eff})",
        "36\\pi e",
        "432\\pi e^2",
        "18665728.0078",
        "additive application-domain correction",
        "selected-factor binary pair marginal",
        "No graviton",
        "no gravity",
        "no G",
    ]
    for token in required_tokens:
        check(token in text, f"packet scope token: {token}")
    check("PASS:" in (LANE / "VERIFICATION.txt").read_text(encoding="utf-8"),
          "verification transcript records pass")


def main() -> None:
    verify_dependencies()
    verify_manifest()
    verify_seal()
    verify_scope()
    print(f"PASS: {CHECKS}/{CHECKS} GL6AZ packet checks")


if __name__ == "__main__":
    main()
