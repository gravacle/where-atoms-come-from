#!/usr/bin/env python3
"""Verify the distinct frozen-byte GL6AK custody and physics replay."""

import hashlib
from pathlib import Path
import subprocess
import sys


here = Path(__file__).resolve().parent
root = here.parent
target = root / "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001"
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise SystemExit("FAIL:" + label)
    checks += 1


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for line in (here / "AUDITED_TARGETS.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = root / relative
    check(path.is_file() and not path.is_symlink(), "regular:" + relative)
    check(sha(path) == expected, "digest:" + relative)

physics = subprocess.run(
    [sys.executable, "-B", str(target / "verify_a3_bulk_dynamics.py")],
    cwd=root, check=True, text=True, capture_output=True,
).stdout
check("PASS 6304/6304" in physics, "author physics")
packet = subprocess.run(
    [sys.executable, "-B", str(target / "verify_mutable_packet.py")],
    cwd=root, check=True, text=True, capture_output=True,
).stdout
check("PASS 104/104" in packet, "author frozen packet")

# Independently replay every frozen dependency and pre-screen audit pin.
for ledger_name in ("DEPENDENCIES.sha256", "PRESCREEN_AUDIT.sha256"):
    for line in (target / ledger_name).read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = root / relative.removeprefix("../")
        check(path.is_file() and not path.is_symlink(), "pin regular:" + relative)
        check(sha(path) == expected, "pin digest:" + relative)

spot = subprocess.run(
    [sys.executable, "-B", str(here / "postfreeze_gl6ak_spot.py")],
    cwd=root, check=True, text=True, capture_output=True,
).stdout
check("PASS__INDEPENDENT_GL6AK_POSTFREEZE_SPOT__33398/33398" in spot,
      "independent spot")
check("NO_GLOBAL_RECORD_STATE_SELECTION_POLE_MOMENTUM_RICCI_GRAVITY_G" in spot,
      "independent ceiling")

text = " ".join((target / "THEOREM.md").read_text().replace("**", "").split())
for token in (
    "does not turn the infinite net into one infinite record",
    "arbitrarily strong or new boundary laws are outside this statement",
    "It is not Poincaré covariance",
    "proves existence, not uniqueness, purity, a ground-state or KMS property",
    "It is not physical momentum",
    "does not say that any sector is gapless, pole-dominated, propagating, or gravitational",
):
    check(token.lower() in text.lower(), "ceiling:" + token)

# Stability check repeated at the end of the potentially long replay.
check(sha(target / "THEOREM.md") ==
      "083d5fbb8a48e27e365167075da132ffa23e395587a4c0e40cc572d8b761ad30",
      "terminal theorem stability")
check(sha(target / "MANIFEST.sha256") ==
      "d38f89c618ea6f77c7b399b005ad0f0abe04d3865e06921f8c765feb44f40620",
      "terminal manifest stability")
check(sha(target / "SEAL.sha256") ==
      "322bf51a00f8fea3f36a09656dda4ebf89ba56b9a88d60b50e9cc7ab33223987",
      "terminal seal stability")

print(f"PASS__GL6AK_POSTFREEZE_HOSTILE_AUDIT__{checks}/{checks}")
print(spot.strip())
