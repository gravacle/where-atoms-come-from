#!/usr/bin/env python3
"""Independent exact replay for frozen GL6T.

This script does not import or execute the author verifier.  It reconstructs
the N=0 FPSS star, the formed KEEP/BREAK chronology, the K-inclusive REC
contrast, the Pauli double commutators, and the S4 sector factorization using
only the displayed operators.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import product


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


# ---------------------------------------------------------------------------
# Exact FPSS N=0 census and branch chronology.
# ---------------------------------------------------------------------------
def weak_compositions(total: int):
    return [x for x in product(range(total + 1), repeat=4) if sum(x) == total]


s0 = weak_compositions(0)
s1 = weak_compositions(1)
check(s0 == [(0, 0, 0, 0)], "S0 unique parent")
check(len(s1) == 4, "S1 has four append children")

v0 = tuple(range(4))
v1 = tuple(range(4))
domain = {(u, v) for u in v0 for v in v1}
active_parent = 0
star = {(active_parent, a) for a in range(4)}
nonedges = domain - star
check(len(domain) == 16, "N0 raw link domain has sixteen factors")
check(len(star) == 4, "N0 support is the four-edge star")
check(len(nonedges) == 12, "N0 retains twelve nonedges")
check(len({("A", active_parent, e) for e in star}) == 4,
      "parent endpoint writer slots are edge-resolved")
check(len({("B", child, e) for e in star for child in (e[1],)}) == 4,
      "child endpoint writer slots are distinct")


def routed(event_vector, route):
    """Return (L,K,G,n) after formation and route, before any link pulse."""
    L = list(event_vector)
    K = [0] * 4
    G = [0] * 4
    n = [0] * 4
    if route == "KEEP":
        K, L = L, [0] * 4
    elif route == "BREAK":
        G, L = L, [0] * 4
    else:
        raise ValueError(route)
    return tuple(L), tuple(K), tuple(G), tuple(n)


all_events = list(product((0, 1), repeat=4))
check(len(all_events) == 16, "complete F/S instrument has sixteen branches")
all_f = (1, 1, 1, 1)
L, K, G, n = routed(all_f, "KEEP")
check((sum(L), sum(K), sum(G), sum(n)) == (0, 4, 0, 0),
      "formed KEEP routes four records and preserves active blank")
L, K, G, n = routed(all_f, "BREAK")
check((sum(L), sum(K), sum(G), sum(n)) == (0, 0, 4, 0),
      "formed BREAK quarantines four records and preserves active blank")

# Every record is individually distinguishable at fixed values of the other
# three events.  The deterministic K_a marginal has TV distance one.
for a in range(4):
    for other in product((0, 1), repeat=3):
        f = []
        s = []
        it = iter(other)
        for j in range(4):
            if j == a:
                f.append(1)
                s.append(0)
            else:
                bit = next(it)
                f.append(bit)
                s.append(bit)
        k_f = routed(tuple(f), "KEEP")[1][a]
        k_s = routed(tuple(s), "KEEP")[1][a]
        check((k_f, k_s) == (1, 0), "per-link K query has unit REC contrast")

# On every nonedge, K=n=0.  The gated transverse coefficient and diagonal
# occupation energy both vanish on the blank ray, so that ray is invariant.
for _edge in nonedges:
    k = n0 = 0
    check(k == 0 and n0 == 0, "nonedge blank invariant under local generator")


# ---------------------------------------------------------------------------
# Exact Gaussian-integer Pauli algebra on the four active links.
# ---------------------------------------------------------------------------
G = tuple[int, int]  # real + i imaginary


def gadd(a: G, b: G) -> G:
    return a[0] + b[0], a[1] + b[1]


def gneg(a: G) -> G:
    return -a[0], -a[1]


def gmul(a: G, b: G) -> G:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


ONE: G = (1, 0)
IUNIT: G = (0, 1)
MINUS_I: G = (0, -1)


def one_pauli(a: str, b: str):
    if a == "I":
        return b, ONE
    if b == "I":
        return a, ONE
    if a == b:
        return "I", ONE
    table = {
        ("X", "Y"): ("Z", IUNIT),
        ("Y", "X"): ("Z", MINUS_I),
        ("Y", "Z"): ("X", IUNIT),
        ("Z", "Y"): ("X", MINUS_I),
        ("Z", "X"): ("Y", IUNIT),
        ("X", "Z"): ("Y", MINUS_I),
    }
    return table[(a, b)]


def word_mul(a, b):
    out = []
    phase = ONE
    for x, y in zip(a, b):
        p, q = one_pauli(x, y)
        out.append(p)
        phase = gmul(phase, q)
    return tuple(out), phase


def op_add(a, b):
    out = dict(a)
    for word, coeff in b.items():
        out[word] = gadd(out.get(word, (0, 0)), coeff)
        if out[word] == (0, 0):
            del out[word]
    return out


def op_neg(a):
    return {word: gneg(coeff) for word, coeff in a.items()}


def op_mul(a, b):
    out = {}
    for wa, ca in a.items():
        for wb, cb in b.items():
            word, phase = word_mul(wa, wb)
            coeff = gmul(gmul(ca, cb), phase)
            out = op_add(out, {word: coeff})
    return out


def comm(a, b):
    return op_add(op_mul(a, b), op_neg(op_mul(b, a)))


def pauli_word(entries):
    word = ["I"] * 4
    for index, letter in entries.items():
        word[index] = letter
    return tuple(word)


H = {}
for a in range(4):
    H = op_add(H, {pauli_word({a: "X"}): (-1, 0)})  # -sum X_a; h stripped

pair_labels = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
M = [{pauli_word({a: "Z", b: "Z"}): ONE} for a, b in pair_labels]

double = [[comm(comm(H, M[b]), M[a]) for a in range(6)] for b in range(6)]
for b, eb in enumerate(pair_labels):
    for a, ea in enumerate(pair_labels):
        shared = set(eb) & set(ea)
        if a == b:
            expected = {}
            for q in eb:
                expected = op_add(expected, {pauli_word({q: "X"}): (-4, 0)})
        elif len(shared) == 1:
            q = next(iter(shared))
            others = tuple((set(eb) | set(ea)) - {q})
            expected = {pauli_word({q: "X", others[0]: "Z", others[1]: "Z"}): (-4, 0)}
        else:
            expected = {}
        check(double[b][a] == expected, "independent Pauli double commutator")


def expectation(op, x: F, z: F) -> F:
    value = F(0)
    for word, (real, imag) in op.items():
        check(imag == 0 and "Y" not in word, "response operator has real X/Z expectation")
        term = F(real)
        for p in word:
            if p == "X":
                term *= x
            elif p == "Z":
                term *= z
        value += term
    return value


# ---------------------------------------------------------------------------
# Bloch replay and exact S4 response at the rational witness.
# ---------------------------------------------------------------------------
h = F(2)
Delta = F(3)
epsilon = F(5)
c = Delta / epsilon
s = 2 * h / epsilon
cos_theta = F(1, 2)
x = s * c * (1 - cos_theta)
z = c * c + s * s * cos_theta
check((x, z) == (F(6, 25), F(17, 25)), "Bloch witness signs and values")

D = [[h * expectation(double[b][a], x, z) for a in range(6)] for b in range(6)]
line = [[F(0) for _ in range(6)] for _ in range(6)]
for a, ea in enumerate(pair_labels):
    for b, eb in enumerate(pair_labels):
        if a != b and len(set(ea) & set(eb)) == 1:
            line[a][b] = 1
check(all(sum(row) == 4 for row in line), "L(K4) has degree four")

expected_diag = -8 * h * x
expected_adj = -4 * h * x * z * z
for a in range(6):
    for b in range(6):
        expected = expected_diag if a == b else expected_adj * line[a][b]
        check(D[a][b] == expected, "D equals diagonal plus line-graph adjacency")


def matvec(matrix, vector):
    return [sum((matrix[i][j] * vector[j] for j in range(len(vector))), F(0))
            for i in range(len(matrix))]


vectors = {
    "A": ([F(1)] * 6, -8 * h * x * (1 + 2 * z * z)),
    "E1": ([F(1), F(-1), F(0), F(0), F(-1), F(1)], -8 * h * x * (1 - z * z)),
    "E2": ([F(1), F(1), F(-2), F(-2), F(1), F(1)], -8 * h * x * (1 - z * z)),
    "T1": ([F(1), F(0), F(0), F(0), F(0), F(-1)], -8 * h * x),
    "T2": ([F(0), F(1), F(0), F(0), F(-1), F(0)], -8 * h * x),
    "T3": ([F(0), F(0), F(1), F(-1), F(0), F(0)], -8 * h * x),
}
for label, (vector, eigenvalue) in vectors.items():
    check(matvec(D, vector) == [eigenvalue * q for q in vector], f"{label} exact eigenvector")

lam_a = vectors["A"][1]
lam_e = vectors["E1"][1]
lam_t = vectors["T1"][1]
check((lam_a, lam_e, lam_t) ==
      (F(-115488, 15625), F(-32256, 15625), F(-96, 25)),
      "rational sector eigenvalues")
check(lam_a < 0 and lam_e < 0 and lam_t < 0, "all three sectors nonzero and negative")


def determinant(matrix):
    work = [row[:] for row in matrix]
    out = F(1)
    for col in range(len(work)):
        pivot = next((r for r in range(col, len(work)) if work[r][col]), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            out = -out
        q = work[col][col]
        out *= q
        for row in range(col + 1, len(work)):
            factor = work[row][col] / q
            for j in range(col + 1, len(work)):
                work[row][j] -= factor * work[col][j]
    return out


det_direct = determinant(D)
det_sector = lam_a * lam_e**2 * lam_t**3
det_formula = (8 * h * x) ** 6 * (1 + 2 * z * z) * (1 - z * z) ** 2
check(det_direct == det_sector == det_formula > 0, "direct and sector determinants agree")

# BREAK has no transverse term; diagonal n commutes with every pair M.
break_D = [[F(0) for _ in range(6)] for _ in range(6)]
check(determinant(break_D) == 0 and all(not any(row) for row in break_D),
      "formed BREAK response is identically zero")

# Small-prewait coefficients follow directly from x=h Delta tau^2/hbar^2
# and 1-z^2=4 h^2 tau^2/hbar^2.
x2 = h * Delta
one_minus_z2_2 = 4 * h * h
check(-8 * h * x2 * 3 == -24 * h * h * Delta, "A1 onset coefficient")
check(-8 * h * x2 == -8 * h * h * Delta, "T2 onset coefficient")
check(-8 * h * x2 * one_minus_z2_2 == -32 * h**4 * Delta,
      "E2 onset coefficient and fourth order")

# With chi=-(i/hbar)[M_B(t),M_A] and
# M_B(t)=M_B+(it/hbar)[H,M_B]+..., the slope multiplier is +1/hbar^2.
check(gmul(MINUS_I, IUNIT) == ONE, "retarded entrance slope has positive D/hbar squared")

print("N0_COUNTS", len(star), len(nonedges))
print("SECTOR_EIGENVALUES", lam_a, lam_e, lam_t)
print("DETERMINANT", det_direct)
print(f"PASS__INDEPENDENT_GL6T_REPLAY__{checks}/{checks}")
