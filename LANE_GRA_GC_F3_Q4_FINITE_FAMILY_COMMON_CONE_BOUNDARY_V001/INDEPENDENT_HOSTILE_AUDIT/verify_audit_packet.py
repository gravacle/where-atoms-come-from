#!/usr/bin/env python3
"""Replay the sealed independent GC hostile-audit packet."""

from hashlib import sha256
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


for line in (HERE / "AUDIT_MANIFEST.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = HERE / relative
    check(path.is_file() and not path.is_symlink() and digest(path) == expected,
          f"audit manifest custody: {relative}")

seal_lines = (HERE / "AUDIT_SEAL.sha256").read_text().splitlines()
manifest_hash, manifest_name = seal_lines[0].split("  ", 1)
verification_hash, verification_name = seal_lines[1].split("  ", 1)
check(manifest_name == "AUDIT_MANIFEST.sha256" and
      manifest_hash == digest(HERE / manifest_name),
      "seal owns the audit manifest")
check(verification_name == "VERIFICATION.txt" and
      verification_hash == digest(HERE / verification_name),
      "seal owns the frozen verification transcript")
check(seal_lines[2] == "VERDICT PASS -- BOUNDED GC CLAIMS ONLY",
      "seal retains the bounded PASS disposition")

replay = subprocess.run(
    [sys.executable, str(HERE / "independent_verify_gc.py")],
    cwd=HERE.parent, capture_output=True, text=True, check=False)
check(replay.returncode == 0,
      "independent physics verifier exits successfully")
check("SUMMARY 66/66 independent hostile GC checks passed" in replay.stdout and
      "VERDICT PASS" in replay.stdout,
      "independent physics verifier reproduces the frozen PASS")
check("no matched-family pole, canonical native-source amplitude, stress Ward "
      "closure, common physical cone, Gravity Formation D1/D2/D3, gravity, or G"
      in replay.stdout,
      "independent replay retains the physical claim ceiling")

print(f"SUMMARY {checks}/{checks} sealed GC audit-packet checks passed")
