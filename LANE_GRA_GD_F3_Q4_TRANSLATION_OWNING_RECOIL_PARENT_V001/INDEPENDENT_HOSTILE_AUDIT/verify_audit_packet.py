#!/usr/bin/env python3
"""Replay the sealed independent GD hostile-audit packet."""

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
      "seal owns the frozen verification summary")
check(seal_lines[2] == "VERDICT PASS -- BOUNDED B1 ALGEBRAIC EXISTENCE ONLY",
      "seal retains the bounded B1 disposition")

independent = subprocess.run(
    [sys.executable, str(HERE / "independent_verify_gd.py")],
    cwd=HERE.parent, capture_output=True, text=True, check=False)
check(independent.returncode == 0,
      "independent physics verifier exits successfully")
check("SUMMARY 154/154 independent hostile GD checks passed" in
      independent.stdout and
      "VERDICT PASS -- BOUNDED B1 ALGEBRAIC EXISTENCE ONLY" in
      independent.stdout,
      "independent physics verifier reproduces the frozen PASS")
check("no physical diamond-space placement, spacetime source, T0j, stress "
      "Ward identity, tensor cone, gravity, or G" in independent.stdout,
      "independent replay retains the physical claim ceiling")

author = subprocess.run(
    [sys.executable, str(HERE.parent / "verify_translation_owning_recoil_parent.py")],
    cwd=HERE.parent, capture_output=True, text=True, check=False)
check(author.returncode == 0 and
      "SUMMARY 144/144 GD translation-owning recoil checks passed" in author.stdout,
      "frozen author verifier reproduces 144/144")

print(f"SUMMARY {checks}/{checks} sealed GD audit-packet checks passed")
