#!/usr/bin/env python3
"""Verify the sealed independent GL6AF hostile-audit packet."""

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
count = 0
for line in (HERE / "MANIFEST.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    if not path.is_file() or path.is_symlink():
        raise SystemExit("FAIL:regular-file:" + relative)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit("FAIL:%s:%s!=%s" % (relative, actual, expected))
    count += 1
print("PASS__GL6AF_AUDIT_PACKET__%d/%d" % (count, count))
