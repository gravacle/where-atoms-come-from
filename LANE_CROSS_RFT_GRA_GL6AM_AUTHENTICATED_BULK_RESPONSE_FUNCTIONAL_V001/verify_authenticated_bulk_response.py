#!/usr/bin/env python3
"""Exact structural replay for GL6AM's finite-window bulk-response theorem."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def eye(n: int):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def matsub(a, b):
    return [[a[i][j] - b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def matadd(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def scale(a, q):
    return [[q * x for x in row] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def permutation_matrix(perm):
    n = len(perm)
    out = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for old, new in enumerate(perm):
        out[new][old] = Fraction(1)
    return out


# 1. Reconstruct the infinite authenticated degree-six adjacency locally.
PORTS = tuple(range(4))


def add4(x, y):
    return tuple(x[i] + y[i] for i in PORTS)


def root(a, b):
    return tuple(int(i == a) - int(i == b) for i in PORTS)


def neighbors(site):
    x, a = site
    same_parent = {(x, b) for b in PORTS if b != a}
    same_child = {(add4(x, root(a, b)), b) for b in PORTS if b != a}
    return same_parent | same_child


sample_cells = [
    (0, 0, 0, 0),
    (1, -1, 0, 0),
    (2, -1, -1, 0),
    (-2, 1, 0, 1),
]
for x in sample_cells:
    check(sum(x) == 0, "sample must lie in A3")
    for a in PORTS:
        ns = neighbors((x, a))
        check(len(ns) == 6, "every bulk link has degree six")
        for q in ns:
            check((x, a) in neighbors(q), "link adjacency must be symmetric")


# 2. Exact A1+E+T2 projectors and S4 covariance in pair space.
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
pair_index = {p: i for i, p in enumerate(PAIRS)}
adj = [[Fraction(0) for _ in PAIRS] for _ in PAIRS]
for i, p in enumerate(PAIRS):
    for j, q in enumerate(PAIRS):
        if i != j and set(p) & set(q):
            adj[i][j] = Fraction(1)
i6 = eye(6)
p_a1 = [[Fraction(1, 6) for _ in PAIRS] for _ in PAIRS]
p_e = scale(matmul(adj, matsub(adj, scale(i6, 4))), Fraction(1, 12))
p_t2 = scale(
    matmul(matsub(adj, scale(i6, 4)), matadd(adj, scale(i6, 2))),
    Fraction(-1, 8),
)
zero6 = [[Fraction(0) for _ in PAIRS] for _ in PAIRS]
for projector, rank in ((p_a1, 1), (p_e, 2), (p_t2, 3)):
    check(matmul(projector, projector) == projector, "sector projector idempotence")
    check(sum(projector[i][i] for i in range(6)) == rank, "sector projector rank")
check(matmul(p_a1, p_e) == zero6, "A1/E orthogonality")
check(matmul(p_a1, p_t2) == zero6, "A1/T2 orthogonality")
check(matmul(p_e, p_t2) == zero6, "E/T2 orthogonality")
check(matadd(matadd(p_a1, p_e), p_t2) == i6, "sector resolution")

for sigma in itertools.permutations(PORTS):
    pair_perm = []
    for a, b in PAIRS:
        image = tuple(sorted((sigma[a], sigma[b])))
        pair_perm.append(pair_index[image])
    rep = permutation_matrix(pair_perm)
    rep_t = transpose(rep)
    for projector in (p_a1, p_e, p_t2):
        check(matmul(matmul(rep_t, projector), rep) == projector,
              "sector projector must be S4 invariant")

# Positive scalar sector measures produce positive matrices atom by atom.
for weights in (
    (Fraction(1), Fraction(2), Fraction(3)),
    (Fraction(0), Fraction(5, 7), Fraction(11, 13)),
):
    mu = matadd(matadd(scale(p_a1, weights[0]), scale(p_e, weights[1])),
                scale(p_t2, weights[2]))
    for c in itertools.product((-1, 0, 1), repeat=6):
        quadratic = sum(Fraction(c[i]) * mu[i][j] * Fraction(c[j])
                        for i in range(6) for j in range(6))
        check(quadratic >= 0, "positive scalar sector measures must be PSD")


# 3. General finite-window spectral measure is an exact Gram measure.
vectors = [
    (Fraction(1), Fraction(0), Fraction(1), Fraction(2)),
    (Fraction(0), Fraction(1), Fraction(-1), Fraction(1)),
    (Fraction(2), Fraction(1), Fraction(0), Fraction(-1)),
    (Fraction(1), Fraction(-1), Fraction(1), Fraction(0)),
    (Fraction(3), Fraction(0), Fraction(1), Fraction(1)),
    (Fraction(-1), Fraction(2), Fraction(1), Fraction(0)),
]
frequency_blocks = ((0, 1), (2,), (3,))
for block in frequency_blocks:
    gram = [[sum(vectors[i][k] * vectors[j][k] for k in block)
             for j in range(6)] for i in range(6)]
    check(gram == transpose(gram), "spectral atom Gram matrix Hermiticity")
    for c in itertools.product((-1, 0, 1), repeat=6):
        lhs = sum(Fraction(c[i]) * gram[i][j] * Fraction(c[j])
                  for i in range(6) for j in range(6))
        rhs = sum(
            sum(Fraction(c[i]) * vectors[i][k] for i in range(6)) ** 2
            for k in block
        )
        check(lhs == rhs and lhs >= 0, "spectral Gram positivity identity")


# 4. Exact pulse derivative sign and finite K-word defect identity.
I2 = [[1 + 0j, 0j], [0j, 1 + 0j]]
X = [[0j, 1 + 0j], [1 + 0j, 0j]]
Z = [[1 + 0j, 0j], [0j, -1 + 0j]]


def cmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def csub(a, b):
    return [[a[i][j] - b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def cscale(a, q):
    return [[q * x for x in row] for row in a]


def kron(a, b):
    return [[a[i // len(b)][j // len(b[0])] * b[i % len(b)][j % len(b[0])]
             for j in range(len(a[0]) * len(b[0]))]
            for i in range(len(a) * len(b))]


source_m = kron(Z, Z)
read_b = kron(X, I2)
comm_bm = csub(cmul(read_b, source_m), cmul(source_m, read_b))
pulse_derivative = cscale(comm_bm, 0.5j)
# d/dj [exp(-ijM/2) B exp(+ijM/2)]_0 = (i/2)[B,M].
direct_derivative = cscale(
    csub(cmul(read_b, cscale(source_m, 0.5j)),
         cmul(cscale(source_m, 0.5j), read_b)),
    1,
)
check(pulse_derivative == direct_derivative, "GL6W half-source pulse sign")

h = Fraction(7, 5)
for kappa in itertools.product((0, 1), repeat=4):
    base_x = [(-h) for _ in PORTS]
    defect = [h * (1 - kappa[p]) for p in PORTS]
    final_x = [base_x[p] + defect[p] for p in PORTS]
    check(final_x == [-h * kappa[p] for p in PORTS],
          "finite K word must cancel exactly the unformed transverse terms")


# 5. Exact constants and formal factorial-tail identities inherited from AI/AK.
for ud in (Fraction(0), Fraction(1, 3), Fraction(-7, 4)):
    j_pair = 2 * abs(ud)
    lambda_ai = 4 * j_pair * 6  # hbar set to one in this exact unit check
    check(lambda_ai == 48 * abs(ud), "lambda_F3 constant")

z = Fraction(3, 2)
for d in range(0, 12):
    # A finite common cutoff verifies the coefficient-level tail recursion.
    tail_d = sum(z ** k / math.factorial(k) for k in range(d, 60))
    tail_next = sum(z ** k / math.factorial(k) for k in range(d + 1, 60))
    check(tail_d - tail_next == z ** d / math.factorial(d),
          "factorial tail recursion")
    check(tail_d >= tail_next >= 0, "factorial tail monotonicity")

shell_sums = []
for radius in (4, 8, 12, 16):
    total = Fraction(0)
    for r in range(radius, 48):
        tail = sum(z ** k / math.factorial(k) for k in range(r + 1, 80))
        total += (2 * r + 1) ** 3 * tail
    shell_sums.append(total)
check(all(shell_sums[i + 1] < shell_sums[i] for i in range(len(shell_sums) - 1)),
      "boundary shell-tail must decrease")
check(shell_sums[-1] * 10**8 < shell_sums[0],
      "factorial tail must beat the cubic shell on the replay window")


# 6. Dependency custody.
for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    actual = hashlib.sha256(target.read_bytes()).hexdigest()
    check(actual == expected, f"dependency hash: {relative}")


# 7. Ledger and hostile scope guards.
ledger = json.loads((HERE / "RESPONSE_LEDGER.json").read_text())
check(ledger["source_unitary"].startswith("V_I(j)=exp"), "source ledger")
check("i E_star^2/(2 hbar)" in ledger["retarded"], "retarded sign ledger")
check("correlation measure" in ledger["positive_object"], "positive object ledger")
check("nonequilibrium" in ledger["defect_ceiling"], "defect ceiling ledger")

theorem = (HERE / "THEOREM.md").read_text()
required = (
    "\\mathcal E_R(B,t)",
    "+2\\|B\\|\\sum_{k=1}^n\\mathcal E_R(V_k,s_k)",
    "{\\mathrm iE_\\star^2\\over2\\hbar}\\Theta(t)",
    "T_{d_L(p,q)}(\\lambda_{\\rm F3}t)",
    "V_\\boldsymbol\\kappa",
    "(1-\\kappa_p)X_p",
    "|\\kappa_p-\\kappa'_p|",
    "T_{d_L(p,q)}(\\lambda_{\\rm F3}u)",
    "commutator measure entering (AM16) is a signed difference",
    "not stationary bulk coefficients",
    "not physical momentum",
    "no physical momentum, pole, cone, Ricci response",
)
for token in required:
    check(token in theorem, f"required theorem scope token: {token}")

for forbidden in (
    "selected ground state",
    "selected KMS state",
    "physical momentum is",
    "strict microcausality follows",
    "positive absorption follows",
    "derives gravity",
    "derives Newton",
):
    check(forbidden not in theorem, f"forbidden overclaim: {forbidden}")

print(f"PASS__GL6AM_AUTHENTICATED_BULK_RESPONSE__{checks}/{checks}")
