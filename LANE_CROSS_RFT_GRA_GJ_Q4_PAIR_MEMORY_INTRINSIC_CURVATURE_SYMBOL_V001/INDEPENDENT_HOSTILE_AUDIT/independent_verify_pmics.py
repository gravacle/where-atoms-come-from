#!/usr/bin/env python3
"""Independent exact audit of PMICS.

This audit reconstructs the algebra from the defining exponential-family
statistics and Fourier curvature convention.  It does not import the builder
verifier or any of its helper code.
"""

from fractions import Fraction as Q
from functools import reduce
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent
passed = 0


def check(statement, label):
    global passed
    if not statement:
        raise AssertionError(label)
    passed += 1
    print(f"PASS {passed:03d} {label}")


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def parity(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant(matrix):
    """Leibniz determinant, deliberately independent of row reduction."""
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("determinant requires a square matrix")
    total = Q(0)
    for sigma in permutations(range(n)):
        term = reduce(
            lambda left, right: left * right,
            (matrix[i][sigma[i]] for i in range(n)),
            Q(parity(sigma)),
        )
        total += term
    return total


def has_nonzero_minor(matrix, order):
    rows = range(len(matrix))
    columns = range(len(matrix[0]))
    for selected_rows in combinations(rows, order):
        for selected_columns in combinations(columns, order):
            minor = [
                [matrix[i][j] for j in selected_columns]
                for i in selected_rows
            ]
            if determinant(minor) != 0:
                return True
    return False


def rank_by_minors(matrix):
    for order in range(min(len(matrix), len(matrix[0])), 0, -1):
        if has_nonzero_minor(matrix, order):
            return order
    return 0


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), Q(0))


def bilinear(left, matrix, right):
    return sum(
        left[i] * matrix[i][j] * right[j]
        for i in range(3)
        for j in range(3)
    )


def symmetrized_outer(left, right):
    return [
        [left[i] * right[j] + right[i] * left[j] for j in range(3)]
        for i in range(3)
    ]


def tensor_coordinates(matrix):
    return (
        matrix[0][0], matrix[1][1], matrix[2][2],
        matrix[0][1], matrix[0][2], matrix[1][2],
    )


def matrix_add(left, right):
    return [
        [left[i][j] + right[i][j] for j in range(3)]
        for i in range(3)
    ]


def matrix_scale(number, matrix):
    return [[number * entry for entry in row] for row in matrix]


def matrix_vector(matrix, vector):
    return [dot(row, vector) for row in matrix]


def trace(matrix):
    return sum((matrix[i][i] for i in range(3)), Q(0))


def even_spin_moment(indices):
    """Uniform independent-spin moment from parity, without state enumeration."""
    return Q(int(all(indices.count(index) % 2 == 0 for index in range(4))))


def directional_fourier_riemann(a, b, c, d, h, k):
    """R(a,b,c,d) after two derivatives act as minus momentum products."""
    return Q(1, 2) * (
        -dot(c, k) * dot(b, k) * bilinear(a, h, d)
        -dot(d, k) * dot(a, k) * bilinear(b, h, c)
        +dot(c, k) * dot(a, k) * bilinear(b, h, d)
        +dot(d, k) * dot(b, k) * bilinear(a, h, c)
    )


def fourier_ricci(h, k):
    k_squared = dot(k, k)
    hk = matrix_vector(h, k)
    return [
        [
            Q(1, 2) * (
                k_squared * h[i][j]
                + k[i] * k[j] * trace(h)
                - k[i] * hk[j]
                - k[j] * hk[i]
            )
            for j in range(3)
        ]
        for i in range(3)
    ]


def fourier_scalar(h, k):
    return dot(k, k) * trace(h) - bilinear(k, h, k)


def fourier_einstein(h, k):
    ricci = fourier_ricci(h, k)
    scalar = fourier_scalar(h, k)
    return [
        [
            ricci[i][j] - Q(1, 2) * scalar * Q(int(i == j))
            for j in range(3)
        ]
        for i in range(3)
    ]


EXPECTED_TARGET = {
    "THEOREM.md": "47f155990137f7467be00e9f7d2cd80633ad35df0639583b72b3262ec44f54a5",
    "RESULT.md": "ff2f2d17aa5f3249a5b80d65719d02c0bf9041f6a696bfbf71e2bd29ddf8a8cc",
    "README.md": "6aaf636596167d6da74c0a27661e5eb457261800862e3922012f3fa86f51e7ac",
    "DEPENDENCIES.sha256": "afc7c08ecbde64342c2eda789f7e8933d068289e21985907692d7b1e1d1be34b",
    "verify_pair_memory_intrinsic_curvature_symbol.py": "417550774549d92afeccc3ac9336e1fd7805a95eab86db9e879b9145a4ac319f",
}

EXPECTED_DEPENDENCIES = {
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


# Exact target and dependency custody.
for relative, expected in EXPECTED_TARGET.items():
    check(file_hash(LANE / relative) == expected, f"target custody {relative}")
for relative, expected in EXPECTED_DEPENDENCIES.items():
    check(file_hash(ROOT / relative) == expected, f"dependency custody {relative}")

declared_lines = {
    line.strip()
    for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
}
expected_lines = {f"{digest}  {relative}" for relative, digest in EXPECTED_DEPENDENCIES.items()}
check(declared_lines == expected_lines, "declared dependency ledger is exact and complete")


# Reconstruct EW's map from tetrahedral contrasts and Walsh parity moments.
v = (
    (Q(1, 2), Q(1, 2), Q(1, 2)),
    (Q(1, 2), Q(-1, 2), Q(-1, 2)),
    (Q(-1, 2), Q(1, 2), Q(-1, 2)),
    (Q(-1, 2), Q(-1, 2), Q(1, 2)),
)
edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

direct_columns = [symmetrized_outer(v[a], v[b]) for a, b in edges]
moment_columns = []
for a, b in edges:
    derivative = []
    for i in range(3):
        row = []
        for j in range(3):
            row.append(sum(
                v[c][i] * v[d][j] * even_spin_moment([a, b, c, d])
                for c in range(4)
                for d in range(4)
            ))
        derivative.append(row)
    moment_columns.append(derivative)

check(moment_columns == direct_columns, "Walsh moments reconstruct every EW D column")
D = [[tensor_coordinates(column)[row] for column in direct_columns] for row in range(6)]
check(determinant(D) == Q(-1, 2), "EW D determinant is exactly -1/2")

identity_three = [[Q(int(i == j)) for j in range(3)] for i in range(3)]
for state in product((-1, 1), repeat=4):
    x = tuple(
        sum((v[a][i] * state[a] for a in range(4)), Q(0))
        for i in range(3)
    )
    xx = [[x[i] * x[j] for j in range(3)] for i in range(3)]
    memory_side = identity_three
    for edge, (a, b) in enumerate(edges):
        memory_side = matrix_add(
            memory_side,
            matrix_scale(Q(state[a] * state[b]), direct_columns[edge]),
        )
    check(xx == memory_side,
          f"statewise observable-memory identity holds on Walsh state {state}")

pair_fisher = [
    [even_spin_moment([*edges[e], *edges[f]]) for f in range(6)]
    for e in range(6)
]
check(pair_fisher == [[Q(int(i == j)) for j in range(6)] for i in range(6)],
      "D_J<C> and pair-score Fisher are identity at uniformity")


# Canonical Fourier-symbol reconstruction with k along the third axis.
q = Q(5)
kz = (Q(0), Q(0), q)
ex = (Q(1), Q(0), Q(0))
ey = (Q(0), Q(1), Q(0))
ez = (Q(0), Q(0), Q(1))
canonical_basis = (
    [[Q(1), Q(0), Q(0)], [Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(0)]],
    [[Q(0), Q(0), Q(0)], [Q(0), Q(1), Q(0)], [Q(0), Q(0), Q(0)]],
    [[Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(1)]],
    [[Q(0), Q(1), Q(0)], [Q(1), Q(0), Q(0)], [Q(0), Q(0), Q(0)]],
    [[Q(0), Q(0), Q(1)], [Q(0), Q(0), Q(0)], [Q(1), Q(0), Q(0)]],
    [[Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(1)], [Q(0), Q(1), Q(0)]],
)
canonical_curvature = [
    [Q(2) * directional_fourier_riemann(a, ez, b, ez, h, kz) / (q * q)
     for h in canonical_basis]
    for a, b in ((ex, ex), (ex, ey), (ey, ey))
]
check(canonical_curvature == [
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 1, 0, 0, 0, 0],
], "canonical Riemann symbol selects exactly the transverse symmetric block")
check(rank_by_minors(canonical_curvature) == 3, "canonical Riemann symbol has rank three")

canonical_einstein = [fourier_einstein(h, kz) for h in canonical_basis]
einstein_coordinates = [
    [tensor_coordinates(column)[row] for column in canonical_einstein]
    for row in range(6)
]
check(rank_by_minors(einstein_coordinates) == 3, "canonical Einstein symbol has rank three")
check(all(matrix_vector(column, kz) == [Q(0)] * 3 for column in canonical_einstein),
      "canonical Einstein symbol is transverse")
check(fourier_scalar(canonical_basis[0], kz) == q * q and
      fourier_scalar(canonical_basis[1], kz) == q * q and
      fourier_scalar(canonical_basis[2], kz) == 0,
      "scalar-curvature sign and longitudinal cancellation are correct")


# Exact FZ momentum witness, reconstructed without normalizing the frame.
r = (Q(7), Q(15), Q(-17))
u = (Q(15), Q(-7), Q(0))
w = (Q(-119), Q(-255), Q(-274))
r_squared = dot(r, r)
check(r_squared == 563, "FZ momentum norm is 563")
check(dot(r, u) == dot(r, w) == dot(u, w) == 0,
      "FZ frame is mutually orthogonal")

screen = [
    [bilinear(a, h, b) for h in direct_columns]
    for a, b in ((u, u), (u, w), (w, w))
]
declared_screen = [
    [88, -88, -32, -242, -88, 88],
    [-2744, 3840, 1496, -1496, -270, -826],
    [-132840, -44712, -32400, 28290, 20500, 6900],
]
check(screen == [[Q(entry) for entry in row] for row in declared_screen],
      "declared FZ-direction screen is reproduced exactly")
selected_minor = [[screen[i][j] for j in (0, 1, 3)] for i in range(3)]
check(determinant(selected_minor) == Q(-173782321152),
      "declared FZ-direction minor and its sign are exact")
check(rank_by_minors(screen) == 3, "FZ-direction curvature screen has rank three")

for edge, h in enumerate(direct_columns):
    for label, a, b in (("uu", u, u), ("uw", u, w), ("ww", w, w)):
        reconstructed = directional_fourier_riemann(a, r, b, r, h, r) / r_squared
        predicted = Q(r_squared, 2) * bilinear(a, h, b)
        check(reconstructed == predicted,
              f"R(a,n,b,n) sign and scale agree for edge {edge}:{label}")


# Kernel: three independent metric gradients lie in a rank-three nullspace.
unit_vectors = (ex, ey, ez)
gauge_columns = [symmetrized_outer(r, direction) for direction in unit_vectors]
gauge_coordinate_matrix = [
    [tensor_coordinates(column)[j] for column in gauge_columns]
    for j in range(6)
]
check(rank_by_minors(gauge_coordinate_matrix) == 3,
      "k symmetrized with arbitrary xi has dimension three")
check(all(
    all(bilinear(a, h, b) == 0 for a, b in ((u, u), (u, w), (w, w)))
    for h in gauge_columns
), "the transverse screen annihilates all three gradient generators")
check(6 - rank_by_minors(screen) == rank_by_minors(gauge_coordinate_matrix),
      "gradient image exhausts the screen kernel by dimension")
check(all(
    directional_fourier_riemann(a, r, b, r, h, r) == 0
    for h in gauge_columns
    for a, b in ((u, u), (u, w), (w, w))
), "Fourier Riemann annihilates every gradient generator")


# Trace plus transverse trace-free census and Einstein trace reversal.
u_squared = dot(u, u)
w_squared = dot(w, w)
trace_row = [screen[0][e] / u_squared + screen[2][e] / w_squared for e in range(6)]
stf_rows = [
    [screen[0][e] / u_squared - screen[2][e] / w_squared for e in range(6)],
    list(screen[1]),
]
check(rank_by_minors([trace_row]) == 1, "transverse trace sector has dimension one")
check(rank_by_minors(stf_rows) == 2, "transverse trace-free sector has dimension two")
check(rank_by_minors([trace_row, *stf_rows]) == 3,
      "trace plus trace-free sectors exhaust the quotient")

for edge, h in enumerate(direct_columns):
    transverse_trace = bilinear(u, h, u) / u_squared + bilinear(w, h, w) / w_squared
    check(fourier_scalar(h, r) == r_squared * transverse_trace,
          f"three-scalar curvature has the declared positive sign on edge {edge}")
    einstein = fourier_einstein(h, r)
    check(matrix_vector(einstein, r) == [Q(0)] * 3,
          f"Einstein symbol is transverse on edge {edge}")
    for label, a, b in (("uu", u, u), ("uw", u, w), ("ww", w, w)):
        expected = Q(r_squared, 2) * (
            bilinear(a, h, b) - dot(a, b) * transverse_trace
        )
        check(bilinear(a, einstein, b) == expected,
              f"Einstein trace reversal agrees on edge {edge}:{label}")


# Explicit tidal representatives and zero-momentum negative control.
h_plus = [[Q(1), Q(0), Q(0)], [Q(0), Q(-1), Q(0)], [Q(0), Q(0), Q(0)]]
h_cross = [[Q(0), Q(1), Q(0)], [Q(1), Q(0), Q(0)], [Q(0), Q(0), Q(0)]]
check(fourier_scalar(h_plus, kz) == fourier_scalar(h_cross, kz) == 0,
      "two independent tidal representatives have zero scalar curvature")
check(fourier_einstein(h_plus, kz) != [[Q(0)] * 3 for _ in range(3)] and
      fourier_einstein(h_cross, kz) != [[Q(0)] * 3 for _ in range(3)],
      "both tidal representatives retain nonzero Einstein curvature")
zero_k = (Q(0), Q(0), Q(0))
check(all(
    fourier_scalar(h, zero_k) == 0 and
    fourier_einstein(h, zero_k) == [[Q(0)] * 3 for _ in range(3)]
    for h in direct_columns
), "uniform first-order deformation has zero Fourier curvature")


# Claim ceiling: mathematical quotient and adopted-law interpretation remain separate.
theorem_text = " ".join((LANE / "THEOREM.md").read_text().split())
result_text = " ".join((LANE / "RESULT.md").read_text().split())
for phrase in (
    "RGRL is an adopted working postulate",
    "does not prove that every record creates curvature",
    "not a four-dimensional Ricci scalar",
    "No graviton",
    "A deeper F3 derivation remains open",
    "flat-reference/principal-symbol qualifier is load bearing",
    "full finite-background linearized curvature operator",
    "observable pair-memory expectation fields `C`",
    "not by the preparation coordinate `J`",
):
    check(phrase in theorem_text or phrase in result_text,
          f"claim ceiling retained: {phrase}")

print(f"SUMMARY {passed}/{passed} PASS")
print("VERDICT PASS_AT_EXACT_FLAT_BACKGROUND_FOURIER_SYMBOL_CEILING")
