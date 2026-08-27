#!/usr/bin/env python3
"""Exact structural checks for the induced-EH back-reaction gate."""

from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
THEOREM = (ROOT / "THEOREM.md").read_text()
RESULT = (ROOT / "RESULT.md").read_text()
AUDIT = (ROOT / "AUDIT.md").read_text()
SOURCES = (ROOT / "PRIMARY_SOURCES.md").read_text()

checks = []


def check(condition, label):
    checks.append((bool(condition), label))


def rank(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    if not a:
        return 0
    rows, cols = len(a), len(a[0])
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        p = a[r][c]
        a[r] = [v / p for v in a[r]]
        for i in range(rows):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = [x - f * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == rows:
            break
    return r


# Proper-time power integrals on the finite fast shell
# [kappa^-2, mu^-2].
kappa = Fraction(7, 3)
mu = Fraction(5, 4)
a = 1 / kappa**2
b = 1 / mu**2
i0 = Fraction(1, 2) * (1 / a**2 - 1 / b**2)
i1 = 1 / a - 1 / b
check(i0 == (kappa**4 - mu**4) / 2, "fast-shell a0 moment")
check(i1 == kappa**2 - mu**2, "fast-shell a1 moment")

# Einstein normalization: C_R=1/(16 pi G) implies 1/(2 C_R)=8 pi G.
# Strip the common symbolic pi and check rational coefficients.
g = Fraction(5, 11)
c_over_pi = 1 / (16 * g)  # C_R / pi^{-1}; symbolic bookkeeping
rhs_over_pi = 1 / (2 * c_over_pi)
check(rhs_over_pi == 8 * g, "metric variation gives 8 pi G normalization")

# Full-rank composite map: D maps collective coordinates (columns) onto all
# three toy physical metric variations (rows).  D^T has no residual null.
d_full = [
    [1, 0, 0, 1],
    [0, 1, 0, 1],
    [0, 0, 1, 1],
]
check(rank(d_full) == 3, "toy collective-to-metric tangent is surjective")

# D^T E=0 is tested by augmenting D rows with candidate E dot columns.
for e in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, -1, 2]):
    pulled = [sum(d_full[i][j] * e[i] for i in range(3)) for j in range(4)]
    check(any(pulled), f"full-rank adjoint detects nonzero residual {e}")

# Rank-deficient countermodel: only two metric combinations vary.  E=(0,0,1)
# is a nonzero Einstein residual invisible to all collective variations.
d_def = [
    [1, 0],
    [0, 1],
    [0, 0],
]
e_hidden = [0, 0, 1]
pulled_hidden = [sum(d_def[i][j] * e_hidden[i] for i in range(3)) for j in range(2)]
check(rank(d_def) == 2, "rank-deficient tangent identified")
check(pulled_hidden == [0, 0], "nonzero residual can hide from deficient pullback")
check(any(e_hidden), "hidden residual is nonzero")

# An independently explicit collective force can balance a nonzero metric
# residual even for the full-rank map; it must vanish before E=0 follows.
e_force = [2, -1, 3]
pulled_force = [sum(d_full[i][j] * e_force[i] for i in range(3)) for j in range(4)]
f_explicit = [-x for x in pulled_force]
check(
    all(x + y == 0 for x, y in zip(pulled_force, f_explicit)),
    "explicit collective force can balance nonzero metric residual",
)
check(any(e_force), "force-balanced metric residual remains nonzero")

# Derivative hierarchy is open: choose positive C_R and bounded C_4 so there
# is a nonzero momentum interval in which C_4 k^2 / C_R < 1.
c_r = Fraction(3, 2)
c_4 = Fraction(5, 7)
k2 = Fraction(1, 10)
check(c_r > 0, "positive induced Ricci coefficient")
check(abs(c_4 * k2 / c_r) < 1, "nonempty leading-EH derivative domain")

required_theorem = [
    "full nonlinear functional",
    "does **not** generally imply",
    "range of `Dg` is dense",
    "Mere nonsurjectivity is not enough",
    "F_A^{\\rm explicit}",
    "explicit origin, sign, and normalization gate",
    "does not guarantee a nonzero positive coefficient",
    "C_R^pre=0",
    "induced contribution/renormalization",
    "Gamma can diagnose formation",
    "projected metric dynamics",
]
for phrase in required_theorem:
    check(phrase in THEOREM, f"theorem contains scope phrase: {phrase}")

required_audit = [
    "does not imply that F3 has earned a metric",
    "does not prove a common cone",
    "does not solve the cosmological-constant problem",
    "non-dense physical range",
    "Gamma is not substituted",
    "convention-only split",
]
for phrase in required_audit:
    check(phrase in AUDIT, f"audit contains ceiling: {phrase}")

check("full nonlinear" in RESULT, "result states nonlinear leading-derivative advance")
check("projected metric dynamics" in RESULT, "result preserves span ceiling")
check("gr-qc/0204062" in SOURCES, "Visser source bound")
check("gr-qc/0106002" in SOURCES, "Barcelo-Liberati-Visser source bound")
check("1007.1246" in SOURCES, "Lorentz-violation source bound")

failed = [label for ok, label in checks if not ok]
for ok, label in checks:
    print(("PASS" if ok else "FAIL"), label)
print(f"SUMMARY {len(checks)-len(failed)}/{len(checks)} PASS")
if failed:
    raise SystemExit(1)
