#!/usr/bin/env python3
"""Verify the sealed independent GL6AI audit packet."""

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise SystemExit("FAIL:" + label)
    checks += 1


for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(), "regular:" + relative)
    check(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
          "digest:" + relative)
seal_expected, seal_relative = (HERE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(ROOT / seal_relative == HERE / "MANIFEST.sha256", "seal-target")
check(
    hashlib.sha256((HERE / "MANIFEST.sha256").read_bytes()).hexdigest()
    == seal_expected,
    "seal-digest",
)
print(f"PASS__GL6AI_AUDIT_PACKET__{checks}/{checks}")
