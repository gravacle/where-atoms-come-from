#!/usr/bin/env python3
"""Exact algebra and proof-custody replay for the GSGB bridge."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def eye(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def mm(a, b):
    bt = transpose(b)
    return [[sum((x * y for x, y in zip(row, col)), F(0))
             for col in bt] for row in a]


def mv(a, v):
    return [sum((x * y for x, y in zip(row, v)), F(0)) for row in a]


def mscale(c, a):
    return [[c * value for value in row] for row in a]


def dot(v, w):
    return sum((x * y for x, y in zip(v, w)), F(0))


def determinant(a):
    work = [row[:] for row in a]
    out = F(1)
    for col in range(len(work)):
        pivot = next((row for row in range(col, len(work))
                      if work[row][col]), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            out = -out
        value = work[col][col]
        out *= value
        for j in range(col, len(work)):
            work[col][j] /= value
        for row in range(col + 1, len(work)):
            factor = work[row][col]
            for j in range(col, len(work)):
                work[row][j] -= factor * work[col][j]
    return out


def rank(a):
    work = [row[:] for row in a]
    rows = len(work)
    cols = len(work[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][col]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [item / value for item in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            factor = work[row][col]
            work[row] = [x - factor * y
                         for x, y in zip(work[row], work[pivot_row])]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def permutation_matrix(order):
    out = [[F(0) for _ in range(4)] for _ in range(4)]
    for source, target in enumerate(order):
        out[target][source] = F(1)
    return out


def inverse_orthogonal(a):
    check(mm(transpose(a), a) == eye(len(a)), "orthogonal inverse exists")
    return transpose(a)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


# EO rank-three contrast algebra and the GSGB QFI scale lock.
I4 = eye(4)
P = [[I4[i][j] - F(1, 4) for j in range(4)] for i in range(4)]
one = [F(1)] * 4
check(mm(P, P) == P, "P is a projector")
check(transpose(P) == P, "P is symmetric")
check(mv(P, one) == [F(0)] * 4, "P removes the common direction")
check(rank(P) == 3, "P has rank three")

for order in permutations(range(4)):
    perm = permutation_matrix(order)
    check(mm(perm, P) == mm(P, perm), f"PERM4 covariance {order}")

f_qfi = F(7, 5)
a2 = F(9, 4)
ell2 = F(4) * a2 / (F(3) * f_qfi)
qfi = mscale(f_qfi, P)
s = mscale(ell2, qfi)
check(s == mscale(F(4) * a2 / F(3), P), "one-scale QFI metric lock")

vectors = [mv(P, [F(int(i == a)) for i in range(4)]) for a in range(4)]
for a in range(4):
    for b in range(4):
        value = dot(vectors[a], mv(s, vectors[b]))
        target = a2 if a == b else -a2 / F(3)
        check(value == target, f"spatial tetrahedral Gram {a},{b}")

K = [[-a2 + dot(vectors[a], mv(s, vectors[b]))
      for b in range(4)] for a in range(4)]
K_expected = mscale(F(4) * a2 / F(3),
                    [[F(int(i == j)) - F(1) for j in range(4)]
                     for i in range(4)])
check(K == K_expected, "full Lorentz Gram")
check(determinant(K) == -F(256, 27) * a2 ** 4,
      "Lorentz Gram determinant")
check(mv(K, one) == [-F(4) * a2] * 4, "timelike eigenvalue")
for contrast in ([F(1), F(-1), F(0), F(0)],
                 [F(0), F(1), F(-1), F(0)],
                 [F(0), F(0), F(1), F(-1)]):
    check(mv(K, contrast) == [F(4) * a2 * x / F(3)
                              for x in contrast],
          "three positive contrast eigenvectors")


# ER finite common-connection contrast.
R1 = [[F(0), F(0), F(1)],
      [F(1), F(0), F(0)],
      [F(0), F(1), F(0)]]
R2 = [[F(0), F(-1), F(0)],
      [F(0), F(0), F(1)],
      [F(-1), F(0), F(0)]]
R1i = inverse_orthogonal(R1)
R2i = inverse_orthogonal(R2)
P21_keep = mm(R1i, R2)
P12_keep = mm(R2i, R1)
H_keep = mm(P21_keep, inverse_orthogonal(P12_keep))
P21_break = mm(R1i, R2i)
P12_break = mm(R2, R1)
check(P21_break == P12_break, "BREAK paths agree")
H_break = mm(P21_break, inverse_orthogonal(P12_break))
check(H_keep != eye(3), "KEEP holonomy is nontrivial")
H_keep_expected = [[F(0), F(-1), F(0)],
                   [F(0), F(0), F(-1)],
                   [F(1), F(0), F(0)]]
check(H_keep == H_keep_expected, "KEEP holonomy is the ER (142) matrix")
check(H_break == eye(3), "BREAK holonomy is flat")
check(determinant(H_keep) == 1, "KEEP holonomy is proper")


# Textual type and anti-overclaim guards.
theorem = (HERE / "THEOREM.md").read_text()
for marker in (
    "gamma_Q[p_K,p_B]",
    "gamma_{\\rm state}",
    "gamma_{\\rm IF}",
    "GSGB-JOIN",
    "operatorname{ROUTE}_{\\rm ER}(R)",
    "S_4`-fixed reference point `xi_*=0 in V`",
    "Sigma_f^{\\mu\\nu}:={1\\over2}\\int_f dx^\\mu\\wedge dx^\\nu",
    "{\\cal F}_{ij}=-2",
    "ker((Dg)^*)={0}",
    "{c_*^3\\over16\\pi G_{\\rm eff}}",
    "{8\\pi G_{\\rm eff}\\over c_*^4}",
    "QFI metric = spacetime metric without G-SOLDER",
    "physical ancestry into gravity",
):
    check(marker in theorem or marker in (HERE / "RESULT.md").read_text(),
          f"required type/boundary marker: {marker}")


# Dependency custody.
dep_file = HERE / "DEPENDENCIES.sha256"
check(dep_file.exists(), "dependency manifest exists")
for line in dep_file.read_text().splitlines():
    if not line.strip():
        continue
    expected, rel = line.split("  ", 1)
    target = ROOT / rel
    check(target.is_file(), f"dependency exists: {rel}")
    check(digest(target) == expected, f"dependency hash: {rel}")

manifest_file = HERE / "MANIFEST.sha256"
check(manifest_file.exists(), "source manifest exists")
for line in manifest_file.read_text().splitlines():
    if not line.strip():
        continue
    expected, rel = line.split("  ", 1)
    target = HERE / rel
    check(target.is_file(), f"manifest source exists: {rel}")
    check(digest(target) == expected, f"manifest source hash: {rel}")

print(f"PASS {checks}/{checks}")
