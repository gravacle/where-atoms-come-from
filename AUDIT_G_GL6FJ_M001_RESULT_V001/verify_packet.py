#!/usr/bin/env python3
"""Verify the closed GL6FJ m001 result-level audit and independently replay it."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows(path):
    result = {}
    for line in path.read_text().splitlines():
        check(line.count("  ") == 1, "checksum separator")
        value, name = line.split("  ", 1)
        pure = PurePosixPath(name)
        check(len(value) == 64 and set(value) <= set("0123456789abcdef") and
              name and not pure.is_absolute() and
              all(part not in {"", ".", ".."} for part in pure.parts) and
              name not in result, "safe unique checksum row")
        result[name] = value
    return result


def closed_audit():
    mode = os.lstat(HERE).st_mode
    check(stat.S_ISDIR(mode) and not stat.S_ISLNK(mode), "ordinary audit directory")
    files = {}
    for root, dirs, names in os.walk(HERE, followlinks=False):
        for name in dirs:
            check(not (Path(root) / name).is_symlink(), "no directory symlink")
        for name in names:
            p = Path(root) / name
            check(not p.is_symlink() and stat.S_ISREG(os.lstat(p).st_mode),
                  "ordinary audit file")
            files[p.relative_to(HERE).as_posix()] = p
    expected = {
        "README.md", "AUDIT_REPORT.md", "TARGET.sha256", "INDEPENDENT_RESULT.json",
        "independent_gl6fj_m001_result_audit.py", "verify_packet.py",
    }
    manifest = rows(HERE / "MANIFEST.sha256")
    check(set(manifest) == expected == (set(files) - {"MANIFEST.sha256", "SEAL.sha256"}),
          "exact audit payload census")
    for name, value in manifest.items():
        check(digest(files[name]) == value, "audit payload digest " + name)
    check(rows(HERE / "SEAL.sha256") == {"MANIFEST.sha256": digest(HERE / "MANIFEST.sha256")},
          "audit seal")


def main():
    closed_audit()
    stored = json.loads((HERE / "INDEPENDENT_RESULT.json").read_text())
    replay = subprocess.run([sys.executable, "-B", str(HERE / "independent_gl6fj_m001_result_audit.py")],
                            cwd=ROOT, text=True, capture_output=True, check=False, timeout=120)
    check(replay.returncode == 0 and replay.stderr == "", "independent audit replay")
    observed = json.loads(replay.stdout)
    check(observed == stored, "independent audit result exact replay")
    check(stored["physical_ward_null_claimed"] is False and stored["gravity_claimed"] is False,
          "claim ceiling retained")
    print(json.dumps({"checks": CHECKS, "result": "PASS__GL6FJ_M001_RESULT_AUDIT_V001",
                      "independent_checks": stored["independent_checks"],
                      "gravity_claimed": False}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, OSError, ValueError, json.JSONDecodeError, subprocess.SubprocessError) as error:
        print("FAIL=" + str(error), file=sys.stderr)
        raise SystemExit(1)
