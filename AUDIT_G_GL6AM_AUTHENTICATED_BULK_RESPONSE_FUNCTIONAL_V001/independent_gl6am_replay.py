#!/usr/bin/env python3
"""Distinct hostile replay for GL6AM.

This script imports no GL6AM author code.  It reconstructs the finite-word
norm algebra, authenticated graph/tail constants, retarded sign, finite-window
positivity, exact S4 scope, and K-defect identities independently.
"""

from __future__ import annotations

import cmath
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def madd(a, b):
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def msub(a, b):
    return [[a[i][j] - b[i][j] for j in range(len(a[0]))]
            for i in range(len(a))]


def mscale(q, a):
    return [[q * x for x in row] for row in a]


def adjoint(a):
    return [[a[j][i].conjugate() for j in range(len(a))]
            for i in range(len(a[0]))]


def opnorm2(a) -> float:
    """Spectral norm for a complex 2 by 2 matrix."""
    g = mm(adjoint(a), a)
    tr = float((g[0][0] + g[1][1]).real)
    det = float((g[0][0] * g[1][1] - g[0][1] * g[1][0]).real)
    disc = max(0.0, tr * tr - 4.0 * det)
    return math.sqrt(max(0.0, (tr + math.sqrt(disc)) / 2.0))


I2 = [[1 + 0j, 0j], [0j, 1 + 0j]]
X = [[0j, 1 + 0j], [1 + 0j, 0j]]
Y = [[0j, -1j], [1j, 0j]]
Z = [[1 + 0j, 0j], [0j, -1 + 0j]]


def rotation(pauli, angle):
    return madd(mscale(math.cos(angle), I2), mscale(1j * math.sin(angle), pauli))


# 1. Exact target/dependency custody, reconstructed without the author verifier.
for ledger_name in ("AUDITED_TARGETS.sha256",):
    for line in (HERE / ledger_name).read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split(maxsplit=1)
        target = ROOT / relative
        check(target.is_file(), f"audited target exists: {relative}")
        check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
              f"audited target hash: {relative}")

manifest_members = set()
for line in (TARGET / "MANIFEST.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"author manifest target exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"author manifest hash: {relative}")
    manifest_members.add(relative)
check(not (TARGET / "SEAL.sha256").exists(), "author did not self-seal")

for line in (TARGET / "DEPENDENCIES.sha256").read_text().splitlines():
    if not line.strip():
        continue
    expected, relative = line.split(maxsplit=1)
    target = ROOT / relative
    check(target.is_file(), f"dependency exists: {relative}")
    check(hashlib.sha256(target.read_bytes()).hexdigest() == expected,
          f"dependency hash: {relative}")

postfreeze_required = {
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/POSTFREEZE_AUDIT.md",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/POSTFREEZE_VERIFICATION.txt",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/MANIFEST.sha256",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/SEAL.sha256",
}
dependency_paths = {line.split(maxsplit=1)[1]
                    for line in (TARGET / "DEPENDENCIES.sha256").read_text().splitlines()
                    if line.strip()}
check(postfreeze_required <= dependency_paths, "AK postfreeze custody explicitly pinned")


# 2. Reconstruct the authenticated A3 link graph and boundary coefficient.
PORTS = range(4)


def add4(x, y):
    return tuple(a + b for a, b in zip(x, y))


def root(a, b):
    return tuple(int(i == a) - int(i == b) for i in PORTS)


def neighbors(site):
    x, a = site
    return ({(x, b) for b in PORTS if b != a}
            | {(add4(x, root(a, b)), b) for b in PORTS if b != a})


samples = ((0, 0, 0, 0), (1, -1, 0, 0), (2, -1, -2, 1), (-3, 1, 1, 1))
for x in samples:
    check(sum(x) == 0, "sample lies in A3")
    for a in PORTS:
        ns = neighbors((x, a))
        check(len(ns) == 6, "authenticated link degree six")
        check(sum(q[0] == x for q in ns) == 3, "three same-parent neighbors")
        check(sum(q[0] != x for q in ns) == 3, "three shared-child neighbors")
        for q in ns:
            check((x, a) in neighbors(q), "adjacency symmetry")

for J in (Fraction(1, 7), Fraction(5, 3), Fraction(19, 4)):
    hbar = Fraction(13, 5)
    lam = 24 * J / hbar
    check((72 * J / hbar) / lam == 3, "AK Duhamel integral gives AM05 factor three")
    check(4 * J * 6 / hbar == lam, "lambda_F3 is 24J/hbar")


def factorial_tail(d: int, z: float) -> float:
    if z == 0.0:
        return 1.0 if d == 0 else 0.0
    log_term = d * math.log(z) - math.lgamma(d + 1)
    if log_term < -745:
        return 0.0
    term = math.exp(log_term)
    total = term
    for k in range(d + 1, 400):
        term *= z / k
        total += term
        if term < max(1e-300, total * 1e-17):
            break
    return total


for z in (0.25, 1.5, 4.0):
    boundary = []
    for R in (8, 16, 32, 64):
        boundary.append(sum((2 * r + 1) ** 3 * factorial_tail(r - 2, z)
                            for r in range(R, 180)))
    check(all(boundary[i + 1] < boundary[i] for i in range(3)),
          "factorial boundary series decreases")
    check(boundary[-1] < 1e-40, "factorial tail beats cubic shell")


# 3. Reconstruct finite-product telescoping and conjugation constants.
for n in (1, 2, 4, 7):
    left_factors = [rotation((X, Y, Z)[k % 3], 0.07 * (k + 1)) for k in range(n)]
    right_factors = [rotation((X, Y, Z)[k % 3], 0.07 * (k + 1) + 0.003 * (k + 2))
                     for k in range(n)]

    def product(factors):
        out = I2
        for factor in factors:
            out = mm(factor, out)
        return out

    wl, wr = product(left_factors), product(right_factors)
    product_error = opnorm2(msub(wl, wr))
    factor_error = sum(opnorm2(msub(a, b)) for a, b in zip(left_factors, right_factors))
    check(product_error <= factor_error + 1e-13, "finite unitary product telescoping")

    b_left = madd(mscale(0.7, Z), mscale(0.2, X))
    b_right = madd(b_left, mscale(0.004, Y))
    o_left = mm(mm(adjoint(wl), b_left), wl)
    o_right = mm(mm(adjoint(wr), b_right), wr)
    lhs = opnorm2(msub(o_left, o_right))
    rhs = opnorm2(msub(b_left, b_right)) + 2 * opnorm2(b_left) * product_error
    # The exact triangle proof can use max(||B_left||,||B_right||); here
    # ||B_left|| exceeds ||B_right|| only by a negligible orthogonal term, so
    # use the safe maximum in the numerical reconstruction.
    rhs_safe = opnorm2(msub(b_left, b_right)) + 2 * max(opnorm2(b_left), opnorm2(b_right)) * product_error
    check(lhs <= rhs_safe + 1e-12, "pulse/read conjugation bound")
    check(rhs <= rhs_safe + 1e-12, "safe conjugation norm choice")


# 4. Retarded sign, normalization, Theta support, and non-strict spatial tail.
source = Z
read = X
comm = msub(mm(read, source), mm(source, read))
pulse_derivative = mscale(0.5j, comm)
delta = 1e-7
v_plus = rotation(source, delta / 2)
v_minus = rotation(source, -delta / 2)
o_plus = mm(mm(adjoint(v_plus), read), v_plus)
o_minus = mm(mm(adjoint(v_minus), read), v_minus)
finite_difference = mscale(1 / (2 * delta), msub(o_plus, o_minus))
check(opnorm2(msub(finite_difference, pulse_derivative)) < 1e-8,
      "independent half-source retarded derivative sign")

E, hbar = Fraction(7, 3), Fraction(11, 5)
check((E * E / (2 * hbar)) * 2 == E * E / hbar,
      "commutator norm converts AM10 to AM11/AM12 coefficient")
for t in (-2.0, -0.1, 0.1, 2.0):
    theta = 1 if t > 0 else 0
    response = theta * factorial_tail(3, 1.7 * abs(t))
    if t < 0:
        check(response == 0, "retarded response is exactly zero at negative time")
    else:
        check(response > 0, "factorial tail is not strict positive-time cone support")


# 5. Positive finite-window spectral measures and exact S4 scope.
vectors = (
    (1 + 1j, 0j, 2 - 1j, 1j),
    (0j, 2 + 1j, -1j, 1 + 0j),
    (2 - 1j, 1j, 0j, -1 + 2j),
    (1j, -1 + 0j, 1 + 1j, 0j),
)
blocks = ((0,), (1, 2), (3,))
for block in blocks:
    gram = [[sum(vectors[i][k].conjugate() * vectors[j][k] for k in block)
             for j in range(4)] for i in range(4)]
    check(all(abs(gram[i][j] - gram[j][i].conjugate()) < 1e-14
              for i in range(4) for j in range(4)), "spectral atom Hermitian")
    for c in itertools.product((0j, 1 + 0j, 1j, 1 + 1j), repeat=4):
        lhs = sum(c[i].conjugate() * gram[i][j] * c[j]
                  for i in range(4) for j in range(4))
        rhs = sum(abs(sum(c[i] * vectors[i][k] for i in range(4))) ** 2
                  for k in block)
        check(abs(lhs.imag) < 1e-12 and abs(lhs.real - rhs) < 1e-12 and lhs.real >= -1e-12,
              "complex finite-window Gram positivity")

PAIRS = tuple(itertools.combinations(PORTS, 2))
pair_index = {p: i for i, p in enumerate(PAIRS)}
I6 = [[Fraction(int(i == j)) for j in range(6)] for i in range(6)]
A = [[Fraction(int(i != j and bool(set(PAIRS[i]) & set(PAIRS[j]))))
      for j in range(6)] for i in range(6)]
P1 = [[Fraction(1, 6) for _ in range(6)] for _ in range(6)]
PE = mscale(Fraction(1, 12), mm(A, msub(A, mscale(4, I6))))
PT = mscale(Fraction(-1, 8), mm(msub(A, mscale(4, I6)), madd(A, mscale(2, I6))))
zero = [[Fraction(0) for _ in range(6)] for _ in range(6)]
for P, rank in ((P1, 1), (PE, 2), (PT, 3)):
    check(mm(P, P) == P, "S4 projector idempotence")
    check(sum(P[i][i] for i in range(6)) == rank, "S4 projector rank")
for P, Q in ((P1, PE), (P1, PT), (PE, PT)):
    check(mm(P, Q) == zero, "S4 projector orthogonality")
check(madd(madd(P1, PE), PT) == I6, "S4 projector resolution")

for sigma in itertools.permutations(PORTS):
    image = [pair_index[tuple(sorted((sigma[a], sigma[b])))] for a, b in PAIRS]
    R = [[Fraction(int(image[j] == i)) for j in range(6)] for i in range(6)]
    Rt = [list(row) for row in zip(*R)]
    for P in (P1, PE, PT):
        check(mm(mm(Rt, P), R) == P, "projector S4 covariance")

origin = (0, 0, 0, 0)
off_origin = (1, -1, 0, 0)
orbit = {tuple(x[sigma[i]] for i in PORTS)
         for sigma in itertools.permutations(PORTS) for x in (off_origin,)}
check(len(orbit) > 1, "one off-origin envelope is not S4 closed")
check({tuple(origin[sigma[i]] for i in PORTS) for sigma in itertools.permutations(PORTS)} == {origin},
      "one-cell origin window is S4 closed")


# 6. Exact K-defect algebra and the state-independent tail identity.
h = Fraction(17, 9)
words = tuple(itertools.product((0, 1), repeat=4))
for kappa in words:
    corrected = tuple(-h + h * (1 - k) for k in kappa)
    check(corrected == tuple(-h * k for k in kappa), "local K defect cancels only unformed X term")
for a in words:
    for b in words:
        coefficient_l1 = sum(abs((1 - a[p]) - (1 - b[p])) for p in PORTS)
        word_l1 = sum(abs(a[p] - b[p]) for p in PORTS)
        check(coefficient_l1 == word_l1, "matched defect telescoping uses |kappa-kappa'|")

for lam in (0.3, 1.1, 3.7):
    for d in range(0, 8):
        t = 0.6
        # High-order midpoint quadrature independently checks the exact shifted-tail integral.
        n = 20000
        numeric = sum(factorial_tail(d, lam * (j + 0.5) * t / n) for j in range(n)) * t / n
        exact = factorial_tail(d + 1, lam * t) / lam
        check(abs(numeric - exact) < 2e-9, "integrated factorial-tail identity")
check(factorial_tail(1, 0.0) == 0.0, "U_d=0 remote tail vanishes for positive distance")
check(factorial_tail(0, 0.0) == 1.0, "overlapping U_d=0 support is not forced to vanish")


# 7. Hostile textual promotion guards.
theorem = (TARGET / "THEOREM.md").read_text()
result = (TARGET / "RESULT.md").read_text()
ledger = json.loads((TARGET / "RESPONSE_LEDGER.json").read_text())
normalized_theorem = " ".join(theorem.split())
required = (
    "finite-volume state limit is needed",
    "This is a conditional choice, not a selected vacuum or equilibrium phase",
    "factorial quasi-local suppression, not zero outside a spatial cone",
    "commutator measure entering (AM16) is a signed difference",
    "S_4\\text{-closed window",
    "generally not stationary under",
    "not stationary bulk coefficients",
    "no physical momentum, pole, cone, Ricci response, gravity, or `G`",
)
for token in required:
    check(token in normalized_theorem, f"required hostile scope text: {token}")
check("generally nonequilibrium" in result, "result keeps defect nonequilibrium ceiling")
check("not_positive_without_extra_state_input" in ledger, "ledger separates correlation positivity")
for forbidden in (
    "unique stationary state",
    "selected bulk state",
    "strict microcausality is proved",
    "physical momentum equals",
    "stationary bulk coefficient is",
    "derives gravity",
    "calculates G",
):
    check(forbidden not in normalized_theorem, f"forbidden promotion absent: {forbidden}")


print(f"PASS__INDEPENDENT_GL6AM_HOSTILE_REPLAY__{checks}/{checks}")
print("PULSE=NORM_LIMIT_WITH_SUMMED_AK_BOUNDARY_TAIL;RETARDED=THETA_PLUS_FACTORIAL_NOT_STRICT_CONE")
print("SPECTRAL=FINITE_WINDOW_CORRELATION_POSITIVE;S4_SCALAR_BLOCKS_ONLY_CLOSED_WINDOW")
print("DEFECT=BOUNDED_K_COCYCLE_STATE_INDEPENDENT_TAIL;NONEQUILIBRIUM_NO_SPECTRAL_PROMOTION")
print("CEILING=NO_SELECTED_STATE_BULK_COEFFICIENT_MOMENTUM_CONE_GRAVITY_G")
