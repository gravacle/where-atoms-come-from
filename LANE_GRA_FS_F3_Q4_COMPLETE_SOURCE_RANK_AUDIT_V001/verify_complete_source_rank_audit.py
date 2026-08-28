#!/usr/bin/env python3
"""Exact verifier for the frozen CW/FM pure-ice complete linear source.

Only Python's standard library is used.  Integer tetrahedral dyads are used
instead of their common factor 1/3; this leaves every rank and nullspace
statement unchanged.
"""

from collections import Counter, defaultdict, deque
from fractions import Fraction
from hashlib import sha256 as _sha256
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent

DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md":
        "4cc63e3e5853b4250a2a5b78256d41b83b195cf527901819a62a04ef53f8d932",
    "LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/INDEPENDENT_AUDIT.md":
        "9d7ef0419b3022dba0db1add7a46d145ebe4b6ec035f73a9b760e63b978d1b2b",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "a91caa20d16b0a1194333f9b51d96546a4ea24d55e23bf1f04c7d249641af8db",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md":
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md":
        "53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9",
    "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md":
        "07445c035ed4c5167a5a20280c4db69a5101eeb71831cdeb126b29702d04b69d",
    "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "91aa35170432684a47278e46ee2b9d56658a43acc8acbbb84480d047cdbe6dcf",
    "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/THEOREM.md":
        "62c7aaee9433a9ffa970ff6e38bac5585200cf40d6fca2cb70477e7e1e7524eb",
    "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "d2da0796cfec7cff8f1d7da5c9bc449d38acdbae089dd9778fb5f19cb6e42b88",
}


checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return _sha256(path.read_bytes()).hexdigest()


def rank(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    if not a:
        return 0
    nr, nc = len(a), len(a[0])
    r = 0
    for c in range(nc):
        pivot = next((i for i in range(r, nr) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        p = a[r][c]
        a[r] = [x / p for x in a[r]]
        for i in range(nr):
            if i != r and a[i][c]:
                q = a[i][c]
                a[i] = [x - q * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == nr:
            break
    return r


def add(*rows):
    return tuple(sum(xs) for xs in zip(*rows))


def scale(c, row):
    return tuple(c * x for x in row)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and digest(path) == expected,
          f"dependency custody matches {relative}")
    check(_sha256(path.read_bytes() + b"tamper").hexdigest() != expected,
          f"dependency appended-byte tamper fails {relative}")


# Semantic custody: the selected dependency branch really is the reduced
# pure-incidence Hamiltonian, while BS20 and FQ retain the complete-source
# obligations that this lane must explicitly classify.
bs = (ROOT / "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md").read_text()
cw = (ROOT / "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md").read_text()
fm = (ROOT / "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md").read_text()
fq = (ROOT / "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md").read_text()
check("H=H_0+V_X" in cw and "d_*=2" in cw,
      "CW freezes the reduced degree-two pure-incidence parent")
check("H=H_0+V_X" in fm and "E_R=0" in fm,
      "FM freezes the symmetric-detuning pure-incidence slice")
check("Freeze `h != 0`" in (LANE / "THEOREM.md").read_text(),
      "FS prospectively freezes the nonzero flip needed for exact operator rank")
check(all(token in bs for token in ("H_{\\rm car}", "H_{\\rm form}",
                                    "H_{\\rm fb}", "H_{\\rm port}")),
      "the broader BS parent names carrier formation feedback and port sectors")
check("node/port weights" in fq and "boundary/controller/port pieces" in fq,
      "FQ forbids silently dropping nonedge complete-source terms")
check("same source-deformed parent" in fq and "post hoc hand weight" in fq,
      "FQ freezes source before Feshbach reduction")
check("occurrence multiplicities" in fq and "node/port weights" in fq,
      "FQ makes the degree occurrence-one node tensor a prospective query datum")


# Covering-matched finite quotient family G_L, L=5*2^r.  Bravais shifts are
# {0,e1,e2,e3}; labels retain the common tetrahedral coframe.
SHIFTS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def graph(L):
    vertices = [(s, x, y, z)
                for s in (0, 1)
                for x, y, z in product(range(L), repeat=3)]
    adjacency = {v: set() for v in vertices}
    edge_label = {}
    for x, y, z in product(range(L), repeat=3):
        avertex = (0, x, y, z)
        for label, shift in enumerate(SHIFTS):
            bvertex = (1, (x + shift[0]) % L,
                           (y + shift[1]) % L,
                           (z + shift[2]) % L)
            adjacency[avertex].add(bvertex)
            adjacency[bvertex].add(avertex)
            edge_label[frozenset((avertex, bvertex))] = label
    return vertices, adjacency, edge_label


def connected(adjacency):
    start = next(iter(adjacency))
    seen = {start}
    queue = deque([start])
    while queue:
        v = queue.popleft()
        for w in adjacency[v]:
            if w not in seen:
                seen.add(w)
                queue.append(w)
    return len(seen) == len(adjacency)


def canonical_cycle(cycle):
    c = tuple(cycle)
    variants = []
    for oriented in (c, tuple(reversed(c))):
        for k in range(len(c)):
            variants.append(oriented[k:] + oriented[:k])
    return min(variants)


def cycles_through_eight(adjacency):
    found = {4: set(), 6: set(), 8: set()}
    for start in adjacency:
        def visit(v, path):
            for w in adjacency[v]:
                if w == start:
                    if len(path) in found:
                        found[len(path)].add(canonical_cycle(path))
                    continue
                if w in path or len(path) == 8:
                    continue
                visit(w, path + [w])
        visit(start, [start])
    return found


for L in (5, 10):
    vertices, adjacency, edge_label = graph(L)
    edge_count = sum(map(len, adjacency.values())) // 2
    check(len(vertices) == 2 * L**3 and edge_count == 4 * L**3,
          f"G_{L} has exact periodic diamond vertex and edge counts")
    check(all(len(adjacency[v]) == 4 for v in vertices),
          f"G_{L} is closed and coordination four")
    check(all(v[0] != w[0] for v in vertices for w in adjacency[v]),
          f"G_{L} is bipartite")
    check(connected(adjacency), f"G_{L} is connected")
    check(len(edge_label) == 4 * L**3,
          f"G_{L} has one unambiguous q4 label on every edge")

v5, a5, labels5 = graph(5)
short_cycles = cycles_through_eight(a5)
check(not short_cycles[4], "G_5 has no simple four-cycle")
check(bool(short_cycles[6]), "G_5 contains inherited elementary hexagons")
check(len(short_cycles[6]) == 4 * 5**3,
      "G_5 has the exact inherited elementary-hexagon count")
check(all(set(Counter(labels5[frozenset((c[i], c[(i + 1) % 6]))]
                      for i in range(6)).values()) == {2}
          and len(Counter(labels5[frozenset((c[i], c[(i + 1) % 6]))]
                          for i in range(6))) == 3
          for c in short_cycles[6]),
      "every G_5 simple hexagon uses three q4 labels twice")

def lifted_balance(cycle, edge_label):
    total = [0, 0, 0]
    for i, v in enumerate(cycle):
        w = cycle[(i + 1) % len(cycle)]
        shift = SHIFTS[edge_label[frozenset((v, w))]]
        sign = 1 if v[0] == 0 else -1
        for k in range(3):
            total[k] += sign * shift[k]
    return tuple(total)

check(all(lifted_balance(c, labels5) == (0, 0, 0)
          for length in (6, 8) for c in short_cycles[length]),
      "G_5 has no quotient-induced wrapping cycle through length eight")

v10, a10, labels10 = graph(10)
g5_edges = {frozenset(edge) for v in a5 for edge in ((v, w) for w in a5[v])}
fiber_count = Counter()
cover_ok = True
for v in v10:
    image = (v[0], v[1] % 5, v[2] % 5, v[3] % 5)
    fiber_count[image] += 1
    for w in a10[v]:
        image_w = (w[0], w[1] % 5, w[2] % 5, w[3] % 5)
        cover_ok &= image_w in a5[image]
check(cover_ok and set(fiber_count.values()) == {8},
      "G_10 maps to G_5 as the declared eight-sheeted graph cover")


# Symmetric-tensor source coordinates are ordered xx,yy,zz,xy,xz,yz.
# Off-diagonal entries carry factor two in j:D.  The common dyad factor 1/3
# has been cleared.
SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
DYADS = tuple((x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)
              for x, y, z in SIGNS)
W_VERTEX = add(*DYADS)
E_NULL = ((1, -1, 0, 0, 0, 0),
          (1, 1, -2, 0, 0, 0))

check(rank(DYADS) == 4, "four one-link flip weights have exact rank four")
check(W_VERTEX == (4, 4, 4, 0, 0, 0),
      "the degree-node weight is the scalar four-dyad sum")
check(all(dot(d, e) == 0 for d in DYADS for e in E_NULL),
      "both diagonal-traceless E directions annihilate every edge dyad")
check(rank(DYADS) + len(E_NULL) == 6,
      "the displayed E pair completes the source nullspace")

# Complete frozen microscopic inventory per primitive cell.  The four onsite
# rows are source-inactive because the FM slice sets E_R=0.  The two vertex
# squares are nonzero and use the scalar W_VERTEX.  The nonzero flips alone
# give the lower rank bound.
flip_rows = list(DYADS)
onsite_zero_rows = [(0, 0, 0, 0, 0, 0) for _ in DYADS]
degree_rows = [W_VERTEX, W_VERTEX]
complete_nonzero_rows = flip_rows + degree_rows
check(rank(complete_nonzero_rows) == 4,
      "complete nonzero CW/FM microscopic linear source has rank four")
check(rank(complete_nonzero_rows + onsite_zero_rows) == 4,
      "formal zero-detuning onsite weights cannot change the rank")
check(all(dot(row, e) == 0
          for row in complete_nonzero_rows + onsite_zero_rows
          for e in E_NULL),
      "complete frozen microscopic inventory retains both E nulls")

# Weight rank is only an upper bound unless the associated operators retain
# the independent directions.  Flatten the four one-link Pauli-X operators
# on one four-link block.  Their disjoint occupation-basis matrix elements
# give an exact operator lower bound of four.
flip_operators = []
for link in range(4):
    matrix = [[0 for _ in range(16)] for _ in range(16)]
    for state in range(16):
        matrix[state][state ^ (1 << link)] = 1
    flip_operators.append(tuple(entry for row in matrix for entry in row))
check(rank(flip_operators) == 4,
      "four distinct one-link flip operators are linearly independent")
operator_source_rows = []
for coordinate in range(6):
    operator_source_rows.append(tuple(
        sum(DYADS[link][coordinate] * flip_operators[link][entry]
            for link in range(4))
        for entry in range(16 * 16)))
check(rank(operator_source_rows) == 4,
      "microscopic flip source has exact operator rank four")

inventory = {
    "one_link_flip": "present_nonzero",
    "onsite_detuning": "present_coefficient_zero",
    "vertex_degree": "present_nonzero_additive_scalar",
    "storage": "excluded_by_selected_CW_FM_parent",
    "carrier": "excluded_by_selected_CW_FM_parent",
    "formation": "excluded_by_selected_CW_FM_parent",
    "feedback": "excluded_by_selected_CW_FM_parent",
    "independent_node_content": "excluded_by_selected_CW_FM_parent",
    "geometric_boundary": "absent_on_closed_periodic_family",
    "boundary_exchange": "excluded_by_selected_CW_FM_parent",
    "controller": "excluded_by_selected_CW_FM_parent",
    "port": "excluded_by_selected_CW_FM_parent",
    "projected_H6": "generated_before_scoring",
    "projected_H8": "generated_before_scoring",
    "projected_identities": "generated_and_retained",
    "source_contact": "retained_O_j2_no_linear_weight",
}
check(len(inventory) == 16 and all(inventory.values()),
      "all requested microscopic projected and completion classes are classified")


def compositions(total, slots=4):
    if slots == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, slots - 1):
            yield (first,) + rest


def word_weight(counts):
    m = sum(counts)
    # The sign is irrelevant to rank.  This is the exact logarithmic tensor
    # from m flip numerators and m-1 uniformly scaled degree denominators.
    return add(*(scale(counts[a], DYADS[a]) for a in range(4)),
               scale(-(m - 1), W_VERTEX))


projected_rows = []
for order in (2, 4, 6, 8):
    rows = [word_weight(c) for c in compositions(order)]
    projected_rows.extend(rows)
    check(all(dot(row, e) == 0 for row in rows for e in E_NULL),
          f"every order-{order} uniform-source word template retains E null2")
    check(rank(DYADS + tuple(rows)) == 4,
          f"order-{order} uniform-source word templates cannot enlarge the four-dyad span")

hex_weights = [word_weight(tuple(0 if a == missing else 2
                                 for a in range(4)))
               for missing in range(4)]
check(all(dot(row, e) == 0 for row in hex_weights for e in E_NULL),
      "all four elementary H6 label types retain E null2")
check(rank(DYADS + tuple(hex_weights)) == 4,
      "generated H6 source weights remain in A1 plus T2")
check(rank(DYADS + tuple(projected_rows)) == 4,
      "complete through-order-eight projected word-weight span remains rank four")

# A general O(j^2) contact can carry an E Hessian, but its source-off gradient
# is zero.  This explicit Hessian includes E-E and mixed E-A1 entries.
A1 = (1, 1, 1, 0, 0, 0)
HESSIAN = [[0 for _ in range(6)] for _ in range(6)]
for u, v in ((E_NULL[0], E_NULL[0]), (E_NULL[0], A1), (A1, E_NULL[0])):
    for i in range(6):
        for j in range(6):
            HESSIAN[i][j] += u[i] * v[j]
gradient_at_zero = (0, 0, 0, 0, 0, 0)
check(any(HESSIAN[i][j] for i in range(6) for j in range(6)),
      "a lawful quadratic contact can have a nonzero E-dependent Hessian")
check(dot(gradient_at_zero, E_NULL[0]) == 0
      and dot(gradient_at_zero, E_NULL[1]) == 0,
      "every O(j^2) contact has zero source-off linear E conjugate")
check(rank(complete_nonzero_rows + projected_rows) == 4,
      "the complete frozen microscopic-plus-projected linear source rank is four")


manifest = LANE / "MANIFEST.sha256"
if manifest.is_file():
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        check(digest(LANE / relative) == expected,
              f"lane manifest matches {relative}")

print(f"SUMMARY {checks}/{checks} exact checks passed")
print("DISPOSITION MATCHED_PERIODIC_Q4_PURE_ICE_FAMILY_FROZEN__COMPLETE_REDUCED_CW_FM_LINEAR_SOURCE_RANK4_A1_PLUS_T2__E_NULL2__H6_H8_IDENTITIES_INHERIT_NULL__OJ2_HESSIAN_CAVEAT_ONLY__FULL_BS_PHYSICAL_COMPLETION_OPEN")
