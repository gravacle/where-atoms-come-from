#!/usr/bin/env python3
"""Verify the frozen independent GL6T audit packet."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
count = 0
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    if not path.is_file() or path.is_symlink():
        raise SystemExit(f"FAIL:regular-file:{relative}")
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f"FAIL:{relative}:{actual}!={expected}")
    count += 1
print(f"PASS__GL6T_AUDIT_PACKET__{count}/{count}")
