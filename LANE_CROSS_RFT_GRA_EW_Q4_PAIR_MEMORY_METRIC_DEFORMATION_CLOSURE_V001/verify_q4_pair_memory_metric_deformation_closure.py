#!/usr/bin/env python3
"""Exact replay for EW q4 pair-memory metric-deformation closure."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import exp, isclose, sqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def eye(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def mm(left, right):
    right_t = transpose(right)
    return [
        [sum((x * y for x, y in zip(row, col)), F(0))
         for col in right_t]
        for row in left
    ]


def mv(matrix, vector):
    return [
        sum((x * y for x, y in zip(row, vector)), F(0))
        for row in matrix
    ]


def madd(left, right):
    return [
        [x + y for x, y in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def mscale(value, matrix):
    return [[value * x for x in row] for row in matrix]


def outer(left, right):
    return [[x * y for y in right] for x in left]


def rank(matrix):
    work = [row[:] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    pivot_row = 0
    for col in range(cols):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][col]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row:
                continue
            value = work[row][col]
            if value:
                work[row] = [
                    entry - value * pivot_entry
                    for entry, pivot_entry in zip(
                        work[row], work[pivot_row]
                    )
                ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant(matrix):
    work = [row[:] for row in matrix]
    out = F(1)
    for col in range(len(work)):
        pivot = next(
            (row for row in range(col, len(work)) if work[row][col]),
            None,
        )
        if pivot is None:
            return F(0)
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            out = -out
        value = work[col][col]
        out *= value
        work[col] = [entry / value for entry in work[col]]
        for row in range(col + 1, len(work)):
            value = work[row][col]
            if value:
                work[row] = [
                    entry - value * pivot_entry
                    for entry, pivot_entry in zip(work[row], work[col])
                ]
    return out


def vec_sym3(matrix):
    return (
        matrix[0][0], matrix[1][1], matrix[2][2],
        matrix[0][1], matrix[0][2], matrix[1][2],
    )


def pull_to_v(matrix):
    return mm(mm(transpose(BASIS), matrix), BASIS)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


I4 = eye(4)
P = [[I4[i][j] - F(1, 4) for j in range(4)] for i in range(4)]
ONE = [F(1)] * 4
EDGES = tuple(combinations(range(4), 2))
STATES = tuple(product((-1, 1), repeat=4))
BASIS = [
    [F(1), F(0), F(0)],
    [F(0), F(1), F(0)],
    [F(0), F(0), F(1)],
    [F(-1), F(-1), F(-1)],
]
VECTORS = [
    [P[row][column] for row in range(4)]
    for column in range(4)
]


# EW01--EW02: exact tetrahedral contrast frame.
check(mm(P, P) == P, "P is a projector")
check(mv(P, ONE) == [F(0)] * 4, "P removes the common direction")
check(rank(P) == 3, "P has rank three")
check(rank(transpose(VECTORS)) == 3, "tetrahedral vectors span V")
check([
    sum(VECTORS[a][i] for a in range(4)) for i in range(4)
] == [F(0)] * 4, "tetrahedral vectors sum to zero")
for a in range(4):
    for b in range(4):
        observed = sum(
            VECTORS[a][i] * VECTORS[b][i] for i in range(4)
        )
        check(observed == F(int(a == b)) - F(1, 4),
              f"tetrahedral Gram {a},{b}")


# EW05--EW10: exact SCGQA symmetric tangent no-go.
delta_basis = [
    [F(1), F(0), F(0), F(-1)],
    [F(0), F(1), F(0), F(-1)],
    [F(0), F(0), F(1), F(-1)],
]
sc_columns = []
for index, delta_p in enumerate(delta_basis):
    diagonal = [
        [delta_p[i] * F(int(i == j)) for j in range(4)]
        for i in range(4)
    ]
    ambient = mm(mm(P, diagonal), P)
    pulled = pull_to_v(ambient)
    sc_columns.append(vec_sym3(pulled))
    trace_v = sum(P[i][i] * delta_p[i] for i in range(4))
    check(trace_v == 0, f"SC tangent trace-free {index}")
check(rank([list(row) for row in zip(*sc_columns)]) == 3,
      "SC tangent map has rank three")
scale_column = vec_sym3(pull_to_v(P))
check(rank([
    list(row) for row in zip(*(sc_columns + [scale_column]))
]) == 4, "SC tangent plus scalar has rank four")

for delta_p in product(range(-2, 3), repeat=4):
    if sum(delta_p) != 0:
        continue
    root_values = [
        delta_p[a] + delta_p[b] for a, b in EDGES
    ]
    check(
        (all(value == 0 for value in root_values))
        == (all(value == 0 for value in delta_p)),
        f"SC tangent kernel control {delta_p}",
    )


# EW11--EW20: exact Walsh family, marginals, and Fisher blocks.
def x_stat(state):
    return tuple(mv(P, [F(value) for value in state]))


def y_stat(state):
    return tuple(F(state[a] * state[b]) for a, b in EDGES)


features = []
for state in STATES:
    single = (
        F(state[0] - state[3]),
        F(state[1] - state[3]),
        F(state[2] - state[3]),
    )
    features.append((F(1),) + single + y_stat(state))
check(rank([list(row) for row in features]) == 10,
      "constant plus nine sufficient statistics independent")

uniform = F(1, len(STATES))
mean_x = [
    sum(uniform * x_stat(state)[i] for state in STATES)
    for i in range(4)
]
mean_y = [
    sum(uniform * y_stat(state)[edge] for state in STATES)
    for edge in range(6)
]
check(mean_x == [F(0)] * 4, "uniform X mean zero")
check(mean_y == [F(0)] * 6, "uniform pair means zero")

f_theta = [[F(0)] * 4 for _ in range(4)]
f_pair = [[F(0)] * 6 for _ in range(6)]
f_mixed = [[F(0)] * 6 for _ in range(4)]
for state in STATES:
    x = x_stat(state)
    y = y_stat(state)
    f_theta = madd(f_theta, mscale(uniform, outer(x, x)))
    f_pair = madd(f_pair, mscale(uniform, outer(y, y)))
    f_mixed = madd(f_mixed, mscale(uniform, outer(x, y)))
check(f_theta == P, "uniform localization Fisher is P")
check(f_pair == eye(6), "uniform pair Fisher is identity")
check(f_mixed == [[F(0)] * 6 for _ in range(4)],
      "uniform mixed Fisher block vanishes")

derivative_matrices = []
for edge_index, (a, b) in enumerate(EDGES):
    derivative = [[F(0)] * 4 for _ in range(4)]
    for state in STATES:
        x = x_stat(state)
        y = y_stat(state)[edge_index]
        derivative = madd(
            derivative, mscale(uniform * y, outer(x, x))
        )
    target = madd(
        outer(VECTORS[a], VECTORS[b]),
        outer(VECTORS[b], VECTORS[a]),
    )
    check(derivative == target,
          f"pair derivative exact edge {a + 1}{b + 1}")
    derivative_matrices.append(target)

metric_columns = [
    vec_sym3(pull_to_v(matrix)) for matrix in derivative_matrices
]
check(rank([list(row) for row in zip(*metric_columns)]) == 6,
      "six pair derivatives span Sym2(V)")

# Orthogonal-complement proof checked on a general rational symmetric tensor.
test_q = [
    [F(2), F(3), F(5)],
    [F(3), F(7), F(11)],
    [F(5), F(11), F(13)],
]
test_q_ambient = mm(
    mm(BASIS, test_q),
    transpose(BASIS),
)
pairings = []
for matrix in derivative_matrices:
    pairings.append(sum(
        test_q_ambient[i][j] * matrix[i][j]
        for i in range(4) for j in range(4)
    ))
check(any(value != 0 for value in pairings),
      "nonzero symmetric tensor detected by derivative basis")


# Exact global-flip matching for a nontrivial rational pair-weight family.
coupling_weights = (F(2), F(3), F(5), F(7), F(11), F(13))
raw_weights = {}
for state in STATES:
    weight = F(1)
    for edge_index, y in enumerate(y_stat(state)):
        if y == 1:
            weight *= coupling_weights[edge_index]
    raw_weights[state] = weight
    flipped = tuple(-value for value in state)
    check(y_stat(flipped) == y_stat(state),
          f"pair statistics global-flip invariant {state}")
normalization = sum(raw_weights.values(), F(0))
probabilities = {
    state: weight / normalization for state, weight in raw_weights.items()
}
for a in range(4):
    plus = sum(
        probability for state, probability in probabilities.items()
        if state[a] == 1
    )
    minus = sum(
        probability for state, probability in probabilities.items()
        if state[a] == -1
    )
    check(plus == minus == F(1, 2),
          f"nontrivial J one-port marginal uniform a={a}")

weighted_x_mean = [
    sum(probabilities[state] * x_stat(state)[i] for state in STATES)
    for i in range(4)
]
check(weighted_x_mean == [F(0)] * 4,
      "global flip keeps localization score mean zero")
weighted_mixed = [[F(0)] * 6 for _ in range(4)]
for state in STATES:
    y_centered = [
        y_stat(state)[edge] - sum(
            probabilities[t] * y_stat(t)[edge] for t in STATES
        )
        for edge in range(6)
    ]
    weighted_mixed = madd(
        weighted_mixed,
        mscale(probabilities[state], outer(x_stat(state), y_centered)),
    )
check(weighted_mixed == [[F(0)] * 6 for _ in range(4)],
      "mixed Fisher block vanishes for nontrivial J")


# S4 covariance and exact representation character.
def permute_state(state, order):
    out = [0] * 4
    for source, target in enumerate(order):
        out[target] = state[source]
    return tuple(out)


def permute_edge(edge, order):
    return tuple(sorted((order[edge[0]], order[edge[1]])))


def cycle_type(order):
    seen = set()
    lengths = []
    for start in range(4):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = order[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


characters = {
    "A1": {(1, 1, 1, 1): 1, (2, 1, 1): 1, (2, 2): 1,
           (3, 1): 1, (4,): 1},
    "E": {(1, 1, 1, 1): 2, (2, 1, 1): 0, (2, 2): 2,
          (3, 1): -1, (4,): 0},
    "T2": {(1, 1, 1, 1): 3, (2, 1, 1): 1, (2, 2): -1,
           (3, 1): 0, (4,): -1},
    "T1": {(1, 1, 1, 1): 3, (2, 1, 1): -1, (2, 2): -1,
           (3, 1): 0, (4,): 1},
    "A2": {(1, 1, 1, 1): 1, (2, 1, 1): -1, (2, 2): 1,
           (3, 1): 1, (4,): -1},
}
edge_character = {}
for order in permutations(range(4)):
    for state in STATES:
        permuted = permute_state(state, order)
        observed_x = x_stat(permuted)
        expected_x = [F(0)] * 4
        original_x = x_stat(state)
        for source, target in enumerate(order):
            expected_x[target] = original_x[source]
        check(list(observed_x) == expected_x,
              f"X is S4 covariant {order} {state}")
        for edge_index, edge in enumerate(EDGES):
            mapped = permute_edge(edge, order)
            mapped_index = EDGES.index(mapped)
            check(y_stat(permuted)[mapped_index] == y_stat(state)[edge_index],
                  f"Y is S4 covariant {order} {state} {edge}")
    fixed_edges = sum(
        permute_edge(edge, order) == edge for edge in EDGES
    )
    fixed_points = sum(order[i] == i for i in range(4))
    squared = tuple(order[order[i]] for i in range(4))
    fixed_points_squared = sum(squared[i] == i for i in range(4))
    sym2_character = (
        (fixed_points - 1) ** 2 + (fixed_points_squared - 1)
    ) // 2
    check(fixed_edges == sym2_character,
          f"edge and Sym2(V) characters agree {order}")
    edge_character[order] = fixed_edges

for name, character in characters.items():
    multiplicity = sum(
        edge_character[order] * character[cycle_type(order)]
        for order in permutations(range(4))
    ) // 24
    expected = 1 if name in {"A1", "E", "T2"} else 0
    check(multiplicity == expected,
          f"edge representation multiplicity {name}")

# Opposite-pair coordinates form an exact six-coordinate transform.
edge_order = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
opposite_transform = [
    [F(1), F(0), F(0), F(0), F(0), F(1)],
    [F(0), F(1), F(0), F(0), F(1), F(0)],
    [F(0), F(0), F(1), F(1), F(0), F(0)],
    [F(1), F(0), F(0), F(0), F(0), F(-1)],
    [F(0), F(1), F(0), F(0), F(-1), F(0)],
    [F(0), F(0), F(1), F(-1), F(0), F(0)],
]
check(edge_order == EDGES, "edge order matches theorem")
check(rank(opposite_transform) == 6,
      "opposite-pair sum/difference coordinates invertible")


# EW27--EW31: exact edge-quadratic map and inverse.
vertices = (
    (F(1), F(1), F(1)),
    (F(1), F(-1), F(-1)),
    (F(-1), F(1), F(-1)),
    (F(-1), F(-1), F(1)),
)
edge_matrix = []
for a, b in EDGES:
    d = [vertices[a][i] - vertices[b][i] for i in range(3)]
    edge_matrix.append([
        d[0] * d[0], d[1] * d[1], d[2] * d[2],
        2 * d[0] * d[1], 2 * d[0] * d[2], 2 * d[1] * d[2],
    ])
check(determinant(edge_matrix) == -F(2) ** 19,
      "edge-quadratic determinant is -2^19")

for q_components in (
    (F(1), F(2), F(3), F(4), F(5), F(6)),
    (F(-2), F(7), F(11), F(-3), F(5, 2), F(13)),
):
    y = mv(edge_matrix, list(q_components))
    y12, y13, y14, y23, y24, y34 = y
    a_sum = y14 + y23
    b_sum = y13 + y24
    c_sum = y12 + y34
    recovered = (
        (a_sum + b_sum - c_sum) / 16,
        (a_sum + c_sum - b_sum) / 16,
        (b_sum + c_sum - a_sum) / 16,
        (y14 - y23) / 16,
        (y13 - y24) / 16,
        (y12 - y34) / 16,
    )
    check(recovered == q_components,
          f"edge-quadratic inverse {q_components}")


# Numerical diagonal-state fidelity witness for the full family.
theta = (0.17, -0.08, 0.03, -0.12)
j_values = (0.11, -0.07, 0.05, 0.02, -0.09, 0.04)
theta_2 = (-0.04, 0.09, -0.02, -0.03)
j_values_2 = (-0.02, 0.03, 0.08, -0.05, 0.01, -0.06)


def model_probability(theta_value, j_value):
    raw = []
    for state in STATES:
        exponent_value = sum(
            theta_value[i] * float(x_stat(state)[i]) for i in range(4)
        ) + sum(
            j_value[edge] * float(y_stat(state)[edge])
            for edge in range(6)
        )
        raw.append(exp(exponent_value))
    total = sum(raw)
    return [value / total for value in raw]


p = model_probability(theta, j_values)
q = model_probability(theta_2, j_values_2)
check(isclose(sum(p), 1.0, abs_tol=2e-14), "first model normalized")
check(isclose(sum(q), 1.0, abs_tol=2e-14), "second model normalized")
check(all(value > 0 for value in p + q), "finite family full support")
bc = sum(sqrt(left * right) for left, right in zip(p, q))
gamma_query = bc * bc
gamma_diagonal_state = bc * bc
check(isclose(gamma_query, gamma_diagonal_state, abs_tol=1e-15),
      "complete query equals diagonal state fidelity")
check(0.0 < gamma_query < 1.0, "distinct parameters have subunit gamma")


# Product sufficient statistics and additive Fisher scaling.
fiber_weights = {}
for left in STATES:
    for right in STATES:
        statistic = tuple(
            x_stat(left)[i] + x_stat(right)[i] for i in range(4)
        ) + tuple(
            y_stat(left)[edge] + y_stat(right)[edge] for edge in range(6)
        )
        weight = raw_weights[left] * raw_weights[right]
        if statistic in fiber_weights:
            check(fiber_weights[statistic] == weight,
                  f"product likelihood constant on statistic fiber {statistic}")
        else:
            fiber_weights[statistic] = weight
for copies in range(1, 7):
    theta_sum = [[F(0)] * 4 for _ in range(4)]
    pair_sum = [[F(0)] * 6 for _ in range(6)]
    for _ in range(copies):
        theta_sum = madd(theta_sum, f_theta)
        pair_sum = madd(pair_sum, f_pair)
    check(mscale(F(copies), f_theta) == theta_sum,
          f"localization Fisher adds N={copies}")
    check(mscale(F(copies), f_pair) == pair_sum,
          f"pair Fisher adds N={copies}")


# Text, dependency, manifest, and byte-custody guards.
required = {
    "THEOREM.md": (
        "MUTABLE_PRESCREEN_READY__BUILDER_REPLAY_PASS__NOT_SEALED",
        "\\delta{\\cal F}_{\\rm SC}\\big|_V",
        "\\operatorname{im}L_{\\rm SC}=T_2",
        "p_{\\theta,J}(s)",
        "\\Pr_J(s_a=+1)=\\Pr_J(s_a=-1)=\\frac12",
        "{\\partial{\\cal F}_\\theta\\over\\partial J_{ab}}",
        "D_J{\\cal F}_\\theta(0)",
        "\\mathbb R^{{\\cal E}_4}=A_1\\oplus E\\oplus T_2",
        "\\det{\\cal Q}_{\\rm edge}=-2^{19}",
        "This is a strong one-body match,\nnot a complete stress match.",
        "Six spatial symmetric variations do not by\n   themselves supply lapse",
    ),
    "README.md": (
        "two-dimensional \\(E\\) shear sector absent",
        "Independent hostile prescreen is\nrequired before any freeze.",
    ),
    "RESULT.md": (
        "SCGQA_SYMMETRIC_TANGENT_NO_GO_EXACT",
        "ACTUAL_WORLD_SOLDERING_AND_GRAVITY_OPEN",
    ),
    "AUDIT.md": (
        "not\nindependent hostile review",
        "cannot seal the lane.",
    ),
    "PRESCREEN_REQUEST.md": (
        "REQUESTED__MUTABLE_SOURCES__NO_FREEZE_AUTHORIZED",
        "Do not edit builder sources.",
    ),
    "DECISIVE_NEXT_TEST.md": (
        "same-parent pair-memory deformation intervention",
        "complete collateral match",
    ),
}
for name, tokens in required.items():
    text = (HERE / name).read_text()
    for token in tokens:
        check(token in text, f"{name} required token: {token}")

dependency_lines = [
    line for line in (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
    if line.strip()
]
check(len(dependency_lines) == 7, "dependency ledger has seven entries")
for line in dependency_lines:
    expected, relative = line.split(None, 1)
    target = (HERE / relative.strip()).resolve()
    check(target.is_file(), f"dependency exists {relative.strip()}")
    check(digest(target) == expected, f"dependency hash {relative.strip()}")

manifest_path = HERE / "MANIFEST.sha256"
if manifest_path.exists():
    expected_names = {
        "README.md", "THEOREM.md", "RESULT.md", "AUDIT.md",
        "PRESCREEN_REQUEST.md", "DECISIVE_NEXT_TEST.md",
        "DEPENDENCIES.sha256", "VERIFICATION.txt",
        "verify_q4_pair_memory_metric_deformation_closure.py",
    }
    lines = [
        line for line in manifest_path.read_text().splitlines() if line.strip()
    ]
    check(len(lines) == len(expected_names), "manifest has nine entries")
    seen = set()
    for line in lines:
        expected, name = line.split(None, 1)
        name = name.strip()
        seen.add(name)
        check((HERE / name).is_file(), f"manifest target exists {name}")
        check(digest(HERE / name) == expected, f"manifest hash {name}")
    check(seen == expected_names, "manifest exact source census")

for path in sorted(HERE.iterdir()):
    if not path.is_file():
        continue
    data = path.read_bytes()
    check(b"\x00" not in data, f"no NUL bytes {path.name}")
    check(b"\r" not in data, f"no CR bytes {path.name}")
    check(data.endswith(b"\n"), f"terminal newline {path.name}")


print("SCGQA_SYMMETRIC_A1_PLUS_T2_NO_GO_EXACT")
print("PAIR_WALSH_FISHER_AND_UNIFORM_MARGINALS_EXACT")
print("SAME_FAMILY_SPATIAL_METRIC_TANGENT_RANK_SIX_EXACT")
print("EDGE_QUADRATIC_RECONSTRUCTION_INVERTIBLE_EXACT")
print("PHYSICAL_SOLDERING_TIME_STRESS_RIEHB_AND_GRAVITY_OPEN")
print(f"PASS {checks}/{checks}")
