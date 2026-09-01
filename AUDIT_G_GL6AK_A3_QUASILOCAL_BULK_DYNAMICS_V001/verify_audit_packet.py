#!/usr/bin/env python3
"""Verify the terminal sealed GL6AK audit packet."""

import hashlib
from pathlib import Path


here = Path(__file__).resolve().parent
root = here.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise SystemExit("FAIL:" + label)
    checks += 1


for line in (here / "MANIFEST.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = root / relative
    check(path.is_file() and not path.is_symlink(), "regular:" + relative)
    check(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
          "digest:" + relative)
expected, relative = (here / "SEAL.sha256").read_text().strip().split("  ", 1)
check(root / relative == here / "MANIFEST.sha256", "seal target")
check(hashlib.sha256((here / "MANIFEST.sha256").read_bytes()).hexdigest() == expected,
      "seal digest")
print(f"PASS__GL6AK_AUDIT_PACKET__{checks}/{checks}")
