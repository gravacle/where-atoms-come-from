#!/usr/bin/env python3
"""Exact checks for the q4 pair-memory intrinsic-curvature symbol theorem."""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent
PASSES = 0


def check(condition, label):
    global PASSES
    if not condition:
        raise AssertionError(label)
    PASSES += 1
    print(f"PASS {PASSES:03d} {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def outer(a, b):
    return [[a[i] * b[j] for j in range(3)] for i in range(3)]


def add(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def sub(a, b):
    return [[a[i][j] - b[i][j] for j in range(3)] for i in range(3)]


def scale(c, a):
    return [[c * a[i][j] for j in range(3)] for i in range(3)]


def matmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), Q(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def matvec(a, v):
    return [dot(row, v) for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def trace(a):
    return sum((a[i][i] for i in range(3)), Q(0))


def contract(a, h, b):
    return dot(a, matvec(h, b))


def coords(h):
    return [h[0][0], h[1][1], h[2][2], h[0][1], h[0][2], h[1][2]]


def from_coords(c):
    xx, yy, zz, xy, xz, yz = c
    return [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]


def rank(rows):
    a = [list(map(Q, row)) for row in rows]
    if not a:
        return 0
    m, n = len(a), len(a[0])
    pivot_row = 0
    for col in range(n):
        pivot = next((row for row in range(pivot_row, m) if a[row][col]), None)
        if pivot is None:
            continue
        a[pivot_row], a[pivot] = a[pivot], a[pivot_row]
        value = a[pivot_row][col]
        a[pivot_row] = [entry / value for entry in a[pivot_row]]
        for row in range(m):
            if row == pivot_row or not a[row][col]:
                continue
            factor = a[row][col]
            a[row] = [a[row][j] - factor * a[pivot_row][j] for j in range(n)]
        pivot_row += 1
        if pivot_row == m:
            break
    return pivot_row


def determinant(square):
    a = [list(map(Q, row)) for row in square]
    n = len(a)
    value = Q(1)
    sign = 1
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        if pivot is None:
            return Q(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            sign *= -1
        pivot_value = a[col][col]
        value *= pivot_value
        for row in range(col + 1, n):
            if not a[row][col]:
                continue
            factor = a[row][col] / pivot_value
            for j in range(col, n):
                a[row][j] -= factor * a[col][j]
    return sign * value


def solve(square, vector):
    n = len(square)
    a = [list(map(Q, square[row])) + [Q(vector[row])] for row in range(n)]
    for col in range(n):
        pivot = next(row for row in range(col, n) if a[row][col])
        a[col], a[pivot] = a[pivot], a[col]
        value = a[col][col]
        a[col] = [entry / value for entry in a[col]]
        for row in range(n):
            if row == col:
                continue
            factor = a[row][col]
            a[row] = [a[row][j] - factor * a[col][j] for j in range(n + 1)]
    return [a[row][-1] for row in range(n)]


def riemann_component(i, j, k, l, h, wavevector):
    """Frozen Fourier convention from PMICS01, with d_a -> i k_a."""
    r = wavevector
    return Q(1, 2) * (
        -r[k] * r[j] * h[i][l]
        -r[l] * r[i] * h[j][k]
        +r[k] * r[i] * h[j][l]
        +r[l] * r[j] * h[i][k]
    )


def riemann_contract(a, b, c, d, h, wavevector):
    return sum((
        a[i] * b[j] * c[k] * d[l] *
        riemann_component(i, j, k, l, h, wavevector)
        for i, j, k, l in product(range(3), repeat=4)
    ), Q(0))


def ricci(h, r):
    r2 = dot(r, r)
    hr = matvec(h, r)
    trh = trace(h)
    return [[Q(1, 2) * (
        r2 * h[i][j] + r[i] * r[j] * trh
        - r[i] * hr[j] - r[j] * hr[i]
    ) for j in range(3)] for i in range(3)]


def scalar_curvature(h, r):
    return dot(r, r) * trace(h) - contract(r, h, r)


def einstein(h, r):
    ric = ricci(h, r)
    scalar = scalar_curvature(h, r)
    return [[ric[i][j] - (Q(1, 2) * scalar if i == j else Q(0))
             for j in range(3)] for i in range(3)]


DEPENDENCIES = {
    "LANE_CROSS_RFT_GRA_ES_GAMMA_SOLDERING_GRAVITY_BRIDGE_V001/THEOREM.md":
        "8721183a72f4b864d06f79f6e68405a5393e37d1a6f4d24d5f4c3c2b79a81075",
    "LANE_CROSS_RFT_GRA_ES_GAMMA_SOLDERING_GRAVITY_BRIDGE_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "b666b6ff258e3339f6cb4f30a956f0bdd814ff2bf93c899d306088c9397e62e1",
    "LANE_CROSS_RFT_GRA_EV_FIDELITY_EDGE_REGGE_SOLDERING_V001/THEOREM.md":
        "cc719b366f54585328c9816e79577b4ec82f91f30e02881c92667934613fa7b6",
    "LANE_CROSS_RFT_GRA_EV_FIDELITY_EDGE_REGGE_SOLDERING_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "2436414511894e26da74bf3b93507bd7554a3a5534fec9bb3724701befff54f0",
    "LANE_CROSS_RFT_GRA_EW_Q4_PAIR_MEMORY_METRIC_DEFORMATION_CLOSURE_V001/THEOREM.md":
        "495e4e99171f4e3e5809f24e5a9a5b68116e996f4f29115bd22f780127d4714e",
    "LANE_CROSS_RFT_GRA_EW_Q4_PAIR_MEMORY_METRIC_DEFORMATION_CLOSURE_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "71c4529413a66ebc657709807643486089a058cf0e84a7964ae153aa4da93984",
    "GRAVITY_RGRL_ADOPTION_V001.md":
        "bca6146dfa2f2a32cea42db43c85c5d5fb1ee7e6114206e321066809e7c0db1f",
    "GRAVITY_RGRL_POST_ADOPTION_STRUCTURAL_THEOREM_V001.md":
        "733b18ecaa29c7acd755db6947b790a9ae37240a3c74d199752d5e278280783d",
    "LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001/THEOREM.md":
        "c0750b7d8a6a7f1b12d3ef76e8d5a6a3754a86e714f75b7efc203a56c7cfeaf9",
    "LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001/INDEPENDENT_HOSTILE_AUDIT/INDEPENDENT_HOSTILE_AUDIT.md":
        "a6ad33451d21f94a530a8b06a89b12c4e1a085de284e353fdaa9893f5748b63e",
}

for relative, expected in DEPENDENCIES.items():
    check(digest(ROOT / relative) == expected, f"dependency hash {relative}")


# EW tetrahedral pair-memory map and its Fisher derivation.
N = ((Q(1), Q(1), Q(1)), (Q(1), Q(-1), Q(-1)),
     (Q(-1), Q(1), Q(-1)), (Q(-1), Q(-1), Q(1)))
V = tuple(tuple(entry / 2 for entry in vector) for vector in N)
PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
D_COLUMNS = tuple(add(outer(V[a], V[b]), outer(V[b], V[a])) for a, b in PAIRS)
D = [[coords(column)[row] for column in D_COLUMNS] for row in range(6)]
check(rank(D) == 6, "EW pair-memory map D has rank six")

states = tuple(product((-1, 1), repeat=4))
ys = [[state[a] * state[b] for a, b in PAIRS] for state in states]
fj = [[sum((Q(y[e] * y[f], 16) for y in ys), Q(0))
       for f in range(6)] for e in range(6)]
check(fj == [[Q(int(e == f)) for f in range(6)] for e in range(6)],
      "uniform pair-character Fisher block is identity")

for edge, (a, b) in enumerate(PAIRS):
    derivative = [[Q(0) for _ in range(3)] for _ in range(3)]
    for state, y in zip(states, ys):
        x = [sum((V[p][i] * state[p] for p in range(4)), Q(0))
             for i in range(3)]
        derivative = add(derivative, scale(Q(y[edge], 16), outer(x, x)))
    check(derivative == D_COLUMNS[edge], f"EW Fisher derivative equals D column {edge}")

i3 = [[Q(int(i == j)) for j in range(3)] for i in range(3)]
for state, y in zip(states, ys):
    x = [sum((V[p][i] * state[p] for p in range(4)), Q(0))
         for i in range(3)]
    reconstructed = [row[:] for row in i3]
    for edge in range(6):
        reconstructed = add(reconstructed, scale(y[edge], D_COLUMNS[edge]))
    check(outer(x, x) == reconstructed,
          f"statewise XX^T expectation-memory identity {state}")


# Exact FZ direction and transverse frame.
r = [Q(7), Q(15), Q(-17)]
r2 = dot(r, r)
u = [Q(15), Q(-7), Q(0)]
w = [Q(-119), Q(-255), Q(-274)]
check(r2 == 563, "FZ direction has squared norm 563")
check(dot(r, u) == dot(r, w) == dot(u, w) == 0,
      "u and w are an orthogonal transverse frame")

identity = [[Q(int(i == j)) for j in range(3)] for i in range(3)]
projector = sub(identity, scale(Q(1, r2), outer(r, r)))
check(matmul(projector, projector) == projector, "transverse projector is idempotent")
check(matvec(projector, r) == [Q(0)] * 3, "transverse projector kills k")
check(rank(projector) == 2, "transverse projector has rank two")


screen = [
    [contract(u, column, u) for column in D_COLUMNS],
    [contract(u, column, w) for column in D_COLUMNS],
    [contract(w, column, w) for column in D_COLUMNS],
]
expected_screen = [
    [88, -88, -32, -242, -88, 88],
    [-2744, 3840, 1496, -1496, -270, -826],
    [-132840, -44712, -32400, 28290, 20500, 6900],
]
check(screen == [list(map(Q, row)) for row in expected_screen],
      "exact FZ-direction curvature-block matrix reproduced")
check(rank(screen) == 3, "pair-memory curvature symbol has rank three")
minor = [[screen[row][column] for column in (0, 1, 3)] for row in range(3)]
check(determinant(minor) == -173782321152,
      "declared (12,13,23) curvature minor is nonzero")


# Independent Fourier-Riemann reconstruction of the same screen.
for edge, h in enumerate(D_COLUMNS):
    for label, a, b in (("uu", u, u), ("uw", u, w), ("ww", w, w)):
        intrinsic = riemann_contract(a, r, b, r, h, r) / r2
        expected = Q(r2, 2) * contract(a, h, b)
        check(intrinsic == expected,
              f"Fourier Riemann equals k^2/2 transverse block {edge}:{label}")


# Gauge kernel and its pullback to the pair-memory coordinates.
basis = ((Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1)))
gauge_columns = tuple(add(outer(r, e), outer(e, r)) for e in basis)
check(rank([coords(column) for column in gauge_columns]) == 3,
      "spatial pure-gradient image has rank three")
for index, h in enumerate(gauge_columns):
    check([contract(u, h, u), contract(u, h, w), contract(w, h, w)] == [0, 0, 0],
          f"curvature screen annihilates pure-gradient column {index}")

preimages = [solve(D, coords(column)) for column in gauge_columns]
check(rank(preimages) == 3, "pure-gradient preimages span three pair-memory directions")
for index, preimage in enumerate(preimages):
    check([dot(row, preimage) for row in screen] == [0, 0, 0],
          f"pulled-back gauge column lies in pair-memory kernel {index}")
check(6 - rank(screen) == rank(preimages),
      "pure-gradient image exhausts the curvature-symbol kernel by dimension")


# Intrinsic scalar-curvature and STF split.
u2, w2 = dot(u, u), dot(w, w)
trace_row = [screen[0][e] / u2 + screen[2][e] / w2 for e in range(6)]
stf_rows = [
    [screen[0][e] / u2 - screen[2][e] / w2 for e in range(6)],
    list(screen[1]),
]
check(rank([trace_row]) == 1, "transverse trace curvature has rank one")
check(rank(stf_rows) == 2, "transverse STF tidal curvature has rank two")
check(rank([trace_row] + stf_rows) == 3, "trace plus STF exhaust curvature quotient")
for edge, h in enumerate(D_COLUMNS):
    scalar = scalar_curvature(h, r)
    transverse_trace = contract(u, h, u) / u2 + contract(w, h, w) / w2
    check(scalar == r2 * transverse_trace,
          f"intrinsic scalar curvature equals k^2 transverse trace {edge}")


# Independent Ricci/Einstein contractions and injectivity on the quotient.
einstein_columns = tuple(einstein(h, r) for h in D_COLUMNS)
einstein_map = [[coords(column)[row] for column in einstein_columns] for row in range(6)]
check(rank(einstein_map) == 3, "spatial Einstein symbol has rank three")
for edge, tensor in enumerate(einstein_columns):
    check(matvec(tensor, r) == [Q(0)] * 3,
          f"spatial Einstein symbol is transverse on D column {edge}")
    h = D_COLUMNS[edge]
    transverse_trace = contract(u, h, u) / u2 + contract(w, h, w) / w2
    for label, a, b in (("uu", u, u), ("uw", u, w), ("ww", w, w)):
        metric_ab = dot(a, b)
        expected = Q(r2, 2) * (contract(a, h, b) - metric_ab * transverse_trace)
        check(contract(a, tensor, b) == expected,
              f"Einstein transverse trace reversal {edge}:{label}")


# Negative controls and anti-overclaim checks.
zero = [Q(0)] * 3
for edge, h in enumerate(D_COLUMNS):
    zero_components = [riemann_component(i, j, k, l, h, zero)
                       for i, j, k, l in product(range(3), repeat=4)]
    check(all(value == 0 for value in zero_components),
          f"uniform k=0 deformation has zero linearized curvature {edge}")
check(rank(stf_rows) == 2 and rank([trace_row] + stf_rows) == 3,
      "dropping the trace leaves only two curvature directions")
bad_projector = identity
bad_gauge = matmul(bad_projector, matmul(gauge_columns[0], bad_projector))
check(bad_gauge != [[Q(0)] * 3 for _ in range(3)],
      "nontransverse screen leaks a pure-gradient mutation")
check(rank(stf_rows) < rank(screen),
      "a rank-two TT response cannot substitute for the rank-three curvature quotient")

theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
theorem_flat = " ".join(theorem.split())
result_flat = " ".join(result.split())
required_theorem_phrases = (
    "not a four-dimensional Ricci scalar",
    "not a metric perturbation",
    "RGRL is an adopted working postulate",
    "does not prove that every record creates curvature",
    "No graviton",
    "flat-reference/principal-symbol qualifier is load bearing",
    "full finite-background linearized curvature operator",
    "observable pair-memory coordinates",
    "natural parameter `J_ab` is a preparation/control coordinate and is not itself a retained record",
)
for phrase in required_theorem_phrases:
    check(phrase in theorem_flat, f"theorem retains ceiling: {phrase}")
check("FZ cannot be substituted for the metric side" in result_flat,
      "result keeps FZ source and metric response type-separated")

print(f"SUMMARY {PASSES}/{PASSES} PASS")
print("DISPOSITION EXACT_Q4_PAIR_MEMORY_INTRINSIC_CURVATURE_SYMBOL_RANK3__PURE_GRADIENT_KERNEL_RANK3__TRACE1_PLUS_TIDAL2__RGRL_PHYSICAL_SOLDER_AXIOMATIC__F3_DYNAMICAL_JOIN_OPEN__NO_GRAVITON_OR_GRAVITY_PROMOTION")
