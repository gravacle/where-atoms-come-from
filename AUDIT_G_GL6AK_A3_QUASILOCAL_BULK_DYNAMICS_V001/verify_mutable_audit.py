#!/usr/bin/env python3
"""Structure check for the unsealed GL6AK hostile pre-freeze audit."""

from pathlib import Path


HERE = Path(__file__).resolve().parent
required = {
    "AUDIT.md",
    "README.md",
    "VERIFICATION.txt",
    "independent_gl6ak_replay.py",
    "verify_mutable_audit.py",
}

checks = 0
for name in sorted(required):
    path = HERE / name
    assert path.is_file(), f"missing {name}"
    assert path.stat().st_size > 0, f"empty {name}"
    checks += 2

assert not (HERE / "MANIFEST.sha256").exists(), "pre-freeze audit must remain unmanifested"
assert not (HERE / "SEAL.sha256").exists(), "pre-freeze audit must remain unsealed"
checks += 2

audit = (HERE / "AUDIT.md").read_text()
for token in (
    "CLEAN_PREFREEZE",
    "AUTHOR_BYTES_NOT_FROZEN_OR_EDITED",
    "distinct post-freeze custody audit",
):
    assert token.lower() in audit.lower(), f"missing audit status token: {token}"
    checks += 1

print(f"PASS {checks}/{checks}")
