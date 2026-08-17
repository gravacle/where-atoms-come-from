"""
CUTOFF GUARD — LANE W17 / ROUTE R2 (THE DISPOSITION)

This lane may read REGISTER_V001.md only up to and including the W-06 row.
Every other script in this directory imports `CUT` from here and never touches
the raw file.  The cut is computed, not asserted: we find the line index of the
W-07 heading and truncate one line before it.
"""
import os, hashlib

REPO = "/Users/bgm/MB Work/where-atoms-come-from"
REG = os.path.join(REPO, "REGISTER_V001.md")

_raw = open(REG, encoding="utf-8").read().replace("\r\n", "\n").split("\n")

# locate the first heading strictly after W-06
_w07 = None
for i, ln in enumerate(_raw):
    if ln.startswith("## W-07"):
        _w07 = i
        break
assert _w07 is not None, "W-07 heading not found; cut cannot be established"

CUT_LAST_LINE = _w07            # 0-indexed; line _w07 itself is EXCLUDED
CUT = _raw[:CUT_LAST_LINE]      # list of lines, 0-indexed; register line N == CUT[N-1]

def line(n):
    """1-indexed register line, cutoff-enforced."""
    assert 1 <= n <= len(CUT), f"line {n} is past the W-06 cutoff (last in-cut line = {len(CUT)})"
    return CUT[n - 1]

def has(n, needle):
    return needle in line(n)

TEXT = "\n".join(CUT)

def sweep(needle):
    return TEXT.count(needle)

if __name__ == "__main__":
    print(f"REGISTER_V001.md total lines      : {len(_raw)}")
    print(f"W-07 heading at line              : {_w07 + 1}  (EXCLUDED)")
    print(f"in-cut lines (W-01 .. W-06)       : {len(CUT)}")
    print(f"in-cut bytes                      : {len(TEXT.encode())}")
    print(f"sha256 of the in-cut slice        : {hashlib.sha256(TEXT.encode()).hexdigest()}")
    # row-level pointer discipline, custody §8: every row carries
    # ruling · where the proof is · exact reopen condition
    rows, cur = [], None
    for ln in CUT:
        if ln.startswith("## "):
            cur = [ln, []]
            rows.append(cur)
        elif cur is not None:
            cur[1].append(ln)
    print("\nROW-LEVEL POINTER DISCIPLINE (custody §1 pointer rule, §8 register form)")
    print(f"{'row':<12} {'WHERE THE PROOF IS':>20} {'sha256':>8} {'REOPEN':>8}")
    for head, body in rows:
        b = "\n".join(body)
        tag = head.split("—")[0].replace("## ", "").strip()[:11]
        print(f"{tag:<12} {str('WHERE THE PROOF IS' in b):>20} "
              f"{str('sha256' in b):>8} {str('REOPEN' in b):>8}")
