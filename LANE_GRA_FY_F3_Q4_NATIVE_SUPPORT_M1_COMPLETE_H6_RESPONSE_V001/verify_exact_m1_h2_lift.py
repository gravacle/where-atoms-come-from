#!/usr/bin/env python3
"""Cheap exact replay of the FY H2 nonzero-momentum lift."""

from contextlib import redirect_stdout
from fractions import Fraction as F
from io import StringIO
from pathlib import Path
import runpy


LANE = Path(__file__).resolve().parent
MAIN = LANE / "derive_native_support_m1_response.py"
with redirect_stdout(StringIO()):
    fy = runpy.run_path(str(MAIN))

checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


direct = fy["direct_ledgers"]
qdiag2 = fy["qdiag2_ledgers"]
linear = fy["ledger_linear"]
freeze = fy["freeze_ledger"]
zero = fy["ZERO_ROW"]

residuals = tuple(linear((1, q2), (1, pair))
                  for q2, pair in zip(qdiag2, direct))
check(len({freeze(residual) for residual in residuals}) == 1,
      "Qdiag2+Qpair is the same native density on all six FO orbits")
residual = residuals[0]
check(all(residual.get((support, cell), zero) ==
          residual.get((support, 0), zero)
          for support in range(fy["SUPPORT_COUNT"])
          for cell in range(fy["CELL_COUNT"])),
      "the exact H2 residual is cell-uniform within every native species")
check(fy["ledger_sum_row"](residual) ==
      (F(-40), F(-40), F(-40), F(0), F(0), F(0)),
      "the uniform residual sums to the frozen homogeneous -40 identity row")
check(fy["exact_m1_relation"](qdiag2, direct, F(-1)),
      "cyclotomic reduction proves Qdiag2(m=1)=-Qpair(m=1)")

print(f"SUMMARY {checks}/{checks} exact H2 m=1 lift checks passed")
print("CLAIM Qdiag2(m=1)=-Qpair(m=1) exactly over Q(zeta_240)")
print("CEILING H4/H6 decided only by full FY replay; no generic response-rank claim")
