#!/usr/bin/env python3
"""Fast custody and claim-surface verifier for the sealed FY packet."""

from hashlib import sha256
import json
from pathlib import Path


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
EXPECTED_FILES = {
    "DEPENDENCIES.sha256", "README.md", "RESULT.json", "RESULT.md",
    "RUN_STATUS.md", "SELF_AUDIT.md", "THEOREM.md", "VERIFICATION.txt",
    "derive_native_support_m1_response.py", "verify_exact_m1_h2_lift.py",
    "verify_packet.py",
}
checks = 0


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"dependency is a regular file: {relative}")
    check(digest(path) == expected, f"dependency hash is pinned: {relative}")

result = json.loads((LANE / "RESULT.json").read_text())
check(result["status"] == "PASS" and result["full_replay"] == "72/72",
      "result records the clean 72/72 full replay")
check(result["exact_m1_diagonal_lift"] == [
          {"order": 2, "coefficient": "-1"},
          {"order": 4, "coefficient": "-37/12"},
          {"order": 6, "coefficient": "-16247/900"}],
      "result records all three exact m=1 lift coefficients")
check(all(sample["ranks"] == {
          "operator": 6, "commutator": 6, "ground_spectral": 6,
          "static_spectral_kernel": 6, "M1": 6,
          "TT_ground_image": 2}
          for sample in result["samples"]),
      "both declared samples retain the frozen rank packet")
check(all(sample["pole_ranks"] == [1, 3, 1, 1]
          for sample in result["samples"]),
      "both declared samples retain pole ranks 1,3,1,1")

joined = " ".join(" ".join((LANE / name).read_text().split())
                     for name in ("THEOREM.md", "RESULT.md",
                                  "SELF_AUDIT.md"))
for phrase in (
        "native-support", "Q(zeta_240)", "Phi_240", "m=1",
        "m=29", "720 forward", "720 reverse", "sampled finite ranks",
        "Ward", "continuum locality", "massless", "RGRL-B", "gravity",
        "Newton"):
    check(phrase in joined, f"claim surface retains: {phrase}")
for forbidden in (
        "proves a graviton", "proves continuum locality",
        "proves the Ward identity", "derives Newton's constant",
        "generic-in-x rank theorem"):
    check(forbidden not in joined, f"forbidden promotion absent: {forbidden}")

manifest = LANE / "MANIFEST.sha256"
if manifest.is_file():
    listed = set()
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        listed.add(relative)
        check(digest(LANE / relative) == expected,
              f"manifest custody: {relative}")
    check(listed == EXPECTED_FILES,
          "manifest covers exactly the eleven frozen builder files")

seal = LANE / "SEAL.sha256"
if seal.is_file():
    expected, relative = seal.read_text().strip().split("  ", 1)
    check(relative == "MANIFEST.sha256" and digest(manifest) == expected,
          "seal pins the complete builder manifest")

print(f"SUMMARY {checks}/{checks} FY packet checks passed")
print("CEILING finite Z30 m=1; sampled response; no Ward/locality/gravity claim")
