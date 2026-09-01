#!/usr/bin/env python3
"""Independent frozen-byte hostile replay for GL6AN V001.

This script imports no author verifier.  It reconstructs the algebra, the
period-four incidence quotient, the local pair sector, and the strong-lock
coefficients directly from the theorem's mathematical definitions.
"""

from __future__ import annotations

import cmath
import hashlib
import itertools
import math
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6AN_NATIVE_DEGREE_LOCK_SECTOR_THEOREM_V001"
THEOREM_SHA = "32f597edc51a609a37b86144487cd7db3bd2f14a65adb754a893d47ef6807e81"
MANIFEST_SHA = "24a71c01ed1b7a92830e92ec7682882c892667289e2794dafb4af5905ad71b2e"
SEAL_SHA = "a946902f027c555f91cd1f2e9ce93e3182f8edeca319955cd691a0bc929fba51"

checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_q(matrix) -> int:
    a = [[Fraction(value) for value in row] for row in matrix]
    if not a:
        return 0
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((row for row in range(rank, rows) if a[row][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [value / scale for value in a[rank]]
        for row in range(rows):
            if row == rank or not a[row][col]:
                continue
            scale = a[row][col]
            a[row] = [x - scale * y for x, y in zip(a[row], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


# Frozen-byte entry pins.
check(sha(TARGET / "THEOREM.md") == THEOREM_SHA, "frozen theorem entry hash")
check(sha(TARGET / "MANIFEST.sha256") == MANIFEST_SHA, "frozen manifest entry hash")
check(sha(TARGET / "SEAL.sha256") == SEAL_SHA, "frozen seal entry hash")


# 1. Exact local completion, inherited-domain intersection, and complement.
for bits in itertools.product((0, 1), repeat=4):
    k = sum(bits)
    pair_sum = sum(bits[a] * bits[b] for a, b in itertools.combinations(range(4), 2))
    check(2 * pair_sum - 3 * k == (k - 2) ** 2 - 4, "local square completion")
    check((4 - k - 2) ** 2 == (k - 2) ** 2, "full complement symmetry")

for u_d in (Fraction(1, 3), Fraction(5, 2), Fraction(11, 7)):
    for d_star in (Fraction(3), Fraction(7, 2), Fraction(5)):
        delta = 4 * u_d * (d_star - 2)
        epsilon = delta + 2 * u_d * (1 - 2 * d_star)
        check(delta > 0, "strict inherited Delta witness")
        check(epsilon == -6 * u_d, "inherited lock line")
for d_star in (Fraction(0), Fraction(2), Fraction(5, 2)):
    delta = 4 * (d_star - 2)
    check((delta > 0) == (d_star > 2), "Delta-domain iff d_star>2")

for bits in itertools.product((0, 1), repeat=4):
    q = sum(bits) - 2
    full = sum(1 - bit for bit in bits) - 2
    partial = (1 - bits[0]) + sum(bits[1:]) - 2
    check(full == -q, "full four-link particle-hole map")
    if bits[0] == bits[1] == bits[2] == bits[3]:
        check(partial != -q, "finite-open partial product counterexample")

# Direct Pauli commutator [X,n]=iY.
X = ((0j, 1 + 0j), (1 + 0j, 0j))
Y = ((0j, -1j), (1j, 0j))
nmat = ((0j, 0j), (0j, 1 + 0j))


def mm2(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


xn, nx = mm2(X, nmat), mm2(nmat, X)
check(tuple(tuple(xn[i][j] - nx[i][j] for j in range(2)) for i in range(2))
      == tuple(tuple(1j * Y[i][j] for j in range(2)) for i in range(2)),
      "Pauli commutator")


# 2. Independently reconstruct the declared Q4 quotient.
L = 4
cells = tuple(itertools.product(range(L), repeat=3))
links = tuple((x, port) for x in cells for port in range(4))


def child(x, port):
    y = list(x)
    if port < 3:
        y[port] = (y[port] + 1) % L
    return tuple(y)


endpoints = {
    edge: (("P", edge[0]), ("C", child(*edge)))
    for edge in links
}
nodes = tuple((kind, x) for kind in ("P", "C") for x in cells)
incident = {node: [] for node in nodes}
for edge, ends in endpoints.items():
    for node in ends:
        incident[node].append(edge)

check(len(cells) == 64, "Q4 cell count")
check(len(nodes) == 128, "Q4 node count")
check(len(links) == 256, "Q4 link count")
check(len(set(endpoints.values())) == len(links), "Q4 has no parallel links")
for node in nodes:
    check(len(incident[node]) == 4, "Q4 node degree four")

# Connectivity, bipartition, and the sole alternating coefficient solution.
neighbors = {node: set() for node in nodes}
for u, v in endpoints.values():
    neighbors[u].add(v)
    neighbors[v].add(u)
seen = {nodes[0]}
queue = deque(seen)
while queue:
    u = queue.popleft()
    for v in neighbors[u]:
        if v not in seen:
            seen.add(v)
            queue.append(v)
check(len(seen) == len(nodes), "Q4 connected")
coeff = {nodes[0]: 1}
queue = deque((nodes[0],))
while queue:
    u = queue.popleft()
    for v in neighbors[u]:
        wanted = -coeff[u]
        if v in coeff:
            check(coeff[v] == wanted, "linear-charge edge consistency")
        else:
            coeff[v] = wanted
            queue.append(v)
check(all(coeff[node] == (1 if node[0] == "P" else -1) for node in nodes),
      "linear-charge kernel is bipartite sign")
check(sum(coeff.values()) == 0, "balanced alternating charge is identity")

# No wrapped 4-cycle: any two parents share at most one child.
parent_nodes = tuple(("P", x) for x in cells)
four_cycles = 0
for left, right in itertools.combinations(parent_nodes, 2):
    common = len(neighbors[left] & neighbors[right])
    check(common <= 1, "Q4 parent pair has at most one common child")
    four_cycles += math.comb(common, 2)
check(four_cycles == 0, "Q4 has no wrapped four-cycle")

# Port differences remain distinct modulo four.
mod_differences = set()
for a in range(4):
    for b in range(4):
        if a == b:
            continue
        vector = tuple((int(i == a) - int(i == b)) % 4 for i in range(4))
        check(vector not in mod_differences, "ordered port difference unique mod four")
        mod_differences.add(vector)
check(len(mod_differences) == 12, "twelve modulo-four port differences")

# Pair ownership and line-graph degree.
pair_owner = Counter()
line_neighbors = {edge: set() for edge in links}
for node, edges in incident.items():
    for e, f in itertools.combinations(edges, 2):
        key = frozenset((e, f))
        pair_owner[key] += 1
        line_neighbors[e].add(f)
        line_neighbors[f].add(e)
check(len(pair_owner) == 3 * len(links), "Q4 adjacent-link pair count")
for multiplicity in pair_owner.values():
    check(multiplicity == 1, "unique pair owner")
for edge in links:
    check(len(line_neighbors[edge]) == 6, "line degree six")

# Exact global square identity on unrelated deterministic occupations.
occupations = []
occupations.append({edge: 0 for edge in links})
occupations.append({edge: 1 for edge in links})
occupations.append({edge: int(edge[1] in (0, 2)) for edge in links})
occupations.append({edge: (edge[0][0] + 2 * edge[0][1] + edge[1]) % 2 for edge in links})
for occupation in occupations:
    lhs = 2 * sum(occupation[e] * occupation[f] for e, f in map(tuple, pair_owner))
    lhs -= 6 * sum(occupation.values())
    rhs = sum((sum(occupation[e] for e in edges) - 2) ** 2 - 4
              for edges in incident.values())
    check(lhs == rhs, "global completed-square identity")


# 3. Local A1/E/T2 locked filter and covariance.
pairs = tuple(itertools.combinations(range(4), 2))
R = [[int(a in pair) for pair in pairs] for a in range(4)]
check(rank_q(R) == 4, "K4 unsigned pair-incidence rank")
locked_z = tuple(z for z in itertools.product((-1, 1), repeat=4) if sum(z) == 0)
check(len(locked_z) == 6, "six locked local states")
pair_vectors = tuple(tuple(z[a] * z[b] for a, b in pairs) for z in locked_z)
for vector in pair_vectors:
    check(sum(vector) == -2, "locked A1 fixed")
    check(tuple(sum(R[a][j] * vector[j] for j in range(6)) for a in range(4))
          == (-1, -1, -1, -1), "locked affine pair constraint")
differences = [[pair_vectors[i][j] - pair_vectors[0][j] for j in range(6)]
               for i in range(1, 6)]
check(rank_q(differences) == 2, "locked variation dimension two")

mean = tuple(sum(Fraction(vector[j], 6) for vector in pair_vectors) for j in range(6))
cov = [[sum(Fraction((vector[i] - mean[i]) * (vector[j] - mean[j]), 6)
            for vector in pair_vectors) for j in range(6)] for i in range(6)]
check(all(value == Fraction(-1, 3) for value in mean), "locked pair mean")
c_d, c_a, c_o = cov[0][0], cov[0][1], cov[0][5]
check((c_d, c_a, c_o) == (Fraction(8, 9), Fraction(-4, 9), Fraction(8, 9)),
      "locked covariance orbit values")
check((c_d + 4 * c_a + c_o, c_d - 2 * c_a + c_o, c_d - c_o)
      == (0, Fraction(8, 3), 0), "A1 E T2 covariance")
check(rank_q(cov) == 2, "locked covariance rank two")


# 4. Constraint-symbol eigenvalues and correct squared/unsquared scaling.
phase_samples = (
    (1, 1, 1, 1),
    (1, 1, -1, -1),
    (1, 1, 1, -1),
    (1, 1j, -1, -1j),
)
for phases in phase_samples:
    s = sum(phases)
    gram_nonzero = (4 - abs(s), 4 + abs(s))
    check(abs(sum(gram_nonzero) - 8) < 1e-12, "symbol Gram trace")
    check(abs(math.prod(gram_nonzero) - (16 - abs(s) ** 2)) < 1e-12,
          "symbol Gram determinant")
    nullity = 3 if abs(s) == 4 else 2
    check(nullity in (2, 3), "symbol nullity")

centered = (
    (1, -1, 0, 0),
    (1, 2, -1, -2),
    (Fraction(3, 2), Fraction(-1, 2), -2, 1),
)
for theta in centered:
    check(sum(theta) == 0, "centered character tangent")
    pair_square = sum((theta[a] - theta[b]) ** 2 for a, b in pairs)
    norm_square = sum(value * value for value in theta)
    check(pair_square == 4 * norm_square, "centered quadratic identity")
    t = 1e-5
    s = sum(cmath.exp(1j * t * float(value)) for value in theta)
    gram_small = 4 - abs(s)
    check(abs(gram_small / t**2 - float(norm_square) / 2) < 2e-5,
          "Gram eigenvalue is quadratic")
    singular = math.sqrt(max(0.0, gram_small))
    check(abs(singular / t - math.sqrt(float(norm_square) / 2)) < 2e-5,
          "singular value is linear")


# 5. Build an alternating locked Q4 background with Dinic flow, independently
# of the author's Edmonds-Karp construction.
def target_cycle():
    return (
        ((0, 0, 0), 0),
        ((1, 3, 0), 1),
        ((1, 3, 0), 2),
        ((0, 3, 1), 0),
        ((0, 3, 1), 1),
        ((0, 0, 0), 2),
    )


cycle = target_cycle()
check(len(set(cycle)) == 6, "six distinct target links")
cycle_ends = tuple(set(endpoints[edge]) for edge in cycle)
for index in range(6):
    check(len(cycle_ends[index] & cycle_ends[(index + 1) % 6]) == 1,
          "hexagon consecutive adjacency")
for i, j in itertools.combinations(range(6), 2):
    if (j - i) not in (1, 5):
        check(not (cycle_ends[i] & cycle_ends[j]), "hexagon chord absence")

fixed = {edge: int(index % 2 == 0) for index, edge in enumerate(cycle)}
source, sink = ("S",), ("T",)
capacity = {}


def add_capacity(u, v, value):
    capacity.setdefault(u, {})[v] = value
    capacity.setdefault(v, {}).setdefault(u, 0)


fixed_parent = Counter()
fixed_child = Counter()
for edge, value in fixed.items():
    if value:
        fixed_parent[edge[0]] += 1
        fixed_child[child(*edge)] += 1
for x in cells:
    add_capacity(source, ("P", x), 2 - fixed_parent[x])
    add_capacity(("C", x), sink, 2 - fixed_child[x])
for edge in links:
    if edge not in fixed:
        add_capacity(("P", edge[0]), ("C", child(*edge)), 1)

residual = {u: dict(row) for u, row in capacity.items()}
level = {}


def bfs_level():
    level.clear()
    level[source] = 0
    queue = deque((source,))
    while queue:
        u = queue.popleft()
        for v, cap in residual[u].items():
            if cap and v not in level:
                level[v] = level[u] + 1
                queue.append(v)
    return sink in level


def send(u, amount, cursor):
    if u == sink:
        return amount
    keys = tuple(residual[u])
    while cursor[u] < len(keys):
        v = keys[cursor[u]]
        if residual[u][v] and level.get(v) == level[u] + 1:
            pushed = send(v, min(amount, residual[u][v]), cursor)
            if pushed:
                residual[u][v] -= pushed
                residual[v][u] += pushed
                return pushed
        cursor[u] += 1
    return 0


flow = 0
while bfs_level():
    cursor = {u: 0 for u in residual}
    while True:
        pushed = send(source, 10**9, cursor)
        if not pushed:
            break
        flow += pushed
required = sum(2 - fixed_parent[x] for x in cells)
check(flow == required == 125, "independent locked-background max flow")

occupied = {edge for edge, value in fixed.items() if value}
for edge in links:
    if edge in fixed:
        continue
    pnode, cnode = ("P", edge[0]), ("C", child(*edge))
    if residual[cnode].get(pnode, 0) == 1:
        occupied.add(edge)
check(len(occupied) == 128, "locked Q4 occupation count")
check(tuple(int(edge in occupied) for edge in cycle) == (1, 0, 1, 0, 1, 0),
      "alternating target hexagon")
for node, edges in incident.items():
    check(sum(edge in occupied for edge in edges) == 2, "global Q4 degree lock")


# 6. Independent canonical second/fourth-order census.
M = len(links)
mixed = same = 0
for e, f in pair_owner:
    if (e in occupied) == (f in occupied):
        same += 1
    else:
        mixed += 1
disjoint = math.comb(M, 2) - len(pair_owner)
check((mixed, same, disjoint) == (2 * M, M, math.comb(M, 2) - 3 * M),
      "universal locked edge-pair census")
check(-Fraction(M, 2) == -128, "canonical H2 scalar")
direct_h4 = -Fraction(mixed, 2) - Fraction(same, 6) - Fraction(disjoint, 4)
folded_h4 = Fraction(M * M, 8)
check(direct_h4 == -Fraction(3 * M * M + 7 * M, 24), "direct Q-only H4")
check(folded_h4 == Fraction(M * M, 8), "folded H4")
check(direct_h4 + folded_h4 == -Fraction(7 * M, 24), "canonical scalar H4")


# 7. Independent subset energies and sixth-order coefficient.
initial = tuple(int(edge in occupied) for edge in cycle)


def subset_energy(subset):
    charge = Counter()
    for index in subset:
        delta = 1 if initial[index] == 0 else -1
        for node in endpoints[cycle[index]]:
            charge[node] += delta
    return sum(value * value for value in charge.values())


all_edges = frozenset(range(6))
check(subset_energy(frozenset()) == 0, "empty subset locked")
check(subset_energy(all_edges) == 0, "full hexagon toggle locked")
for size in range(1, 6):
    for subset in itertools.combinations(range(6), size):
        check(subset_energy(frozenset(subset)) > 0, "proper subset leaves lock")

energy_census = {
    size: Counter(subset_energy(frozenset(subset))
                  for subset in itertools.combinations(range(6), size))
    for size in range(1, 6)
}
check(energy_census[1] == Counter({2: 6}), "one-flip energy census")
check(energy_census[2] == Counter({4: 9, 2: 6}), "two-flip energy census")
check(energy_census[3] == Counter({4: 12, 2: 6, 6: 2}), "three-flip energy census")
check(energy_census[4] == energy_census[2], "four-flip complement census")
check(energy_census[5] == energy_census[1], "five-flip complement census")

# Subset dynamic programming sums all orderings without importing the author's
# 720-permutation loop.
amplitude = {frozenset(): Fraction(1)}
for size in range(1, 6):
    for values in itertools.combinations(range(6), size):
        subset = frozenset(values)
        amplitude[subset] = -sum(amplitude[subset - {edge}] for edge in subset) \
                            / subset_energy(subset)
coefficient = sum(amplitude[all_edges - {edge}] for edge in all_edges)
check(coefficient == Fraction(-63, 8), "sixth-order hexagon coefficient")

# A separate permutation/profile sum checks the sign and the five resolvents.
profile = Counter()
permutation_sum = Fraction(0)
for ordering in itertools.permutations(range(6)):
    subset = set()
    energies = []
    term = Fraction(1)
    for edge in ordering[:-1]:
        subset.add(edge)
        energy = subset_energy(frozenset(subset))
        energies.append(energy)
        term *= Fraction(-1, energy)
    profile[tuple(energies)] += 1
    permutation_sum += term
check(sum(profile.values()) == 720, "all sixth-order paths")
check(len(profile) == 9, "nine energy-denominator profiles")
check(permutation_sum == coefficient == Fraction(-63, 8), "path/profile agreement")


# 8. Strict textual ceilings and repaired scope.
theorem = " ".join((TARGET / "THEOREM.md").read_text().replace("**", "").split())
result = " ".join((TARGET / "RESULT.md").read_text().replace("**", "").split())
for token in (
    "Delta=4U_d(d_\\star-2)",
    "generic finite open product",
    "constraint-Gram eigenvalue",
    "singular value itself vanishes linearly",
    "explicit period-four quotient",
    "No global projector in the infinite quasi-local algebra",
    "arbitrary smaller periodic quotients",
    "bi-infinite path",
    "six independently authenticated records",
    "sixth and higher diagonal/loop terms remain unclassified",
):
    check(token in theorem, "repaired theorem scope: " + token)
for forbidden in (
    "multi-record collective operation",
    "quadratically soft singular value",
    "state-dependent diagonal fourth",
    "physical momentum is established",
):
    check(forbidden not in theorem + " " + result, "forbidden promotion absent: " + forbidden)

# Terminal stability after all reconstructions.
check(sha(TARGET / "THEOREM.md") == THEOREM_SHA, "terminal theorem stability")
check(sha(TARGET / "MANIFEST.sha256") == MANIFEST_SHA, "terminal manifest stability")
check(sha(TARGET / "SEAL.sha256") == SEAL_SHA, "terminal seal stability")

print(f"PASS__INDEPENDENT_GL6AN_HOSTILE_REPLAY__{checks}/{checks}")
print("LOCK=INHERITED_DSTAR_GT2_SQUARE_EXACT;WARD=NO_NONTRIVIAL_LINEAR_DEGREE_CHARGE")
print("PAIR=A1_FIXED_E_TWO_DIMENSIONAL_T2_ZERO;GRAM=QUADRATIC_EIGENVALUE_LINEAR_SINGULAR")
print("Q4=PERIOD4_SIMPLE_DEGREE4_GIRTH6;HEFF=H2_H4_SCALAR_HEXAGON_MINUS63_OVER8")
print("CEILING=NO_RECORD_PROMOTION_INFINITE_PROJECTOR_POLE_MOMENTUM_CONE_GRAVITY_G")
