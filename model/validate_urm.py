"""VALIDATE THE T-54 URM FAMILY LAYERS, THEN THE EXISTING GEOMETRY/PROJECT STACK.

The four family blocks remain independently countable and keep their sealed-lane
semantics: ARROW 27, COUNTLAW 40, CLASSES 52, WRITING 57 (176 gates total).  The
chain then runs validate_geometry.py, which in turn runs validate_project.py and its
D-25 construction/anchor gates.  Exit code is zero only when every family gate, the
declared gate counts, geometry, and project validation all pass.

Use --no-chain to run only the four T-54 family blocks.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from checks_arrow import run_arrow_checks
from checks_countlaw import run_countlaw_checks
from checks_classes import run_classes_checks
from checks_writing import run_writing_checks


n_pass = 0
n_fail = 0


def check(name, cond, detail=""):
    global n_pass, n_fail
    if cond:
        n_pass += 1
        print(f"  PASS  {name}  {detail}")
    else:
        n_fail += 1
        print(f"  FAIL  {name}  {detail}")


FAMILIES = (
    ("ARROW", 27, run_arrow_checks),
    ("COUNTLAW", 40, run_countlaw_checks),
    ("CLASSES", 52, run_classes_checks),
    ("WRITING", 57, run_writing_checks),
)

print("VALIDATE THE T-54 URM FAMILY LAYERS")
print("=" * 78)
t0 = time.time()
count_shape_ok = True
family_rows = []
for label, expected, runner in FAMILIES:
    print()
    print(f"{label}: expected {expected} gates")
    print("-" * 78)
    p0, f0 = n_pass, n_fail
    tf = time.time()
    runner(check)
    passed, failed = n_pass - p0, n_fail - f0
    observed = passed + failed
    shape_ok = observed == expected
    count_shape_ok &= shape_ok
    family_rows.append((label, passed, failed, observed, expected, time.time() - tf))
    print(f"  {label}: {passed} PASS, {failed} FAIL; gates {observed}/{expected} "
          f"({'COUNT OK' if shape_ok else 'COUNT MISMATCH'}; {time.time() - tf:.1f} s)")

print()
print("=" * 78)
print(f"  T-54 FAMILIES: {n_pass} PASS, {n_fail} FAIL; "
      f"gates {n_pass + n_fail}/176 ({'COUNT OK' if count_shape_ok else 'COUNT MISMATCH'})")

chain_ok = True
if "--no-chain" not in sys.argv:
    print()
    print("CHAIN: validate_geometry.py (which chains the 24-gate project/D-25 validator)")
    print("-" * 78)
    sys.stdout.flush()
    result = subprocess.run([sys.executable, os.path.join(HERE, "validate_geometry.py")],
                            cwd=HERE)
    chain_ok = result.returncode == 0
    print("-" * 78)
    print(f"  CHAIN validate_geometry.py: {'PASS' if chain_ok else 'FAIL'}")

print("=" * 78)
overall = n_fail == 0 and count_shape_ok and chain_ok
print(f"  URM OVERALL: {'PASS' if overall else 'FAIL'} "
      f"(families {n_pass}/{n_pass + n_fail}, "
      f"counts {'ok' if count_shape_ok else 'MISMATCH'}"
      + ("" if "--no-chain" in sys.argv else f", chain {'ok' if chain_ok else 'FAILED'}")
      + f"; {time.time() - t0:.1f} s)")
sys.exit(0 if overall else 1)
