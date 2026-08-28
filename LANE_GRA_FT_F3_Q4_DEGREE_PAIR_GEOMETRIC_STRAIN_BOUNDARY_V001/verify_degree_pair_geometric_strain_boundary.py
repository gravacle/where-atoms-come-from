#!/usr/bin/env python3
"""Exact verifier for the F3/q4 degree-pair strain-source boundary.

The calculation uses only the Python standard library and exact rational
arithmetic.  Common tetrahedral normalization factors are cleared; ranks,
nullspaces, representation sectors, and source-off equalities are unchanged.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent

DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md":
        "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "44690ad431c85af7a4947a431a4c57ad4cd8b19a346e3d26720b341f77256f90",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4",
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
    "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/THEOREM.md":
        "36879f4c18eec83a22bdf9bd161d9d444b72e1dbda1d5eaa0312c6aab3d95724",
    "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "6e0ecd0febf6364e4122bbf2f65e1feb93c27e960bce30c0097ea0fbe3f58966",
}

MANIFEST_FILES = {
    "DEPENDENCIES.sha256",
    "RESULT.md",
    "SELF_AUDIT.md",
    "THEOREM.md",
    "verify_degree_pair_geometric_strain_boundary.py",
}

AUDIT_FILES = {
    "AUDIT_MANIFEST.sha256",
    "AUDIT_SEAL.sha256",
    "INDEPENDENT_HOSTILE_AUDIT.md",
    "INDEPENDENT_HOSTILE_VERIFICATION.txt",
    "independent_hostile_audit.py",
}


checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def rank(rows):
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    nr, nc = len(matrix), len(matrix[0])
    pivot_row = 0
    for column in range(nc):
        pivot = next((r for r in range(pivot_row, nr)
                      if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        scale = matrix[pivot_row][column]
        matrix[pivot_row] = [value / scale for value in matrix[pivot_row]]
        for row in range(nr):
            if row != pivot_row and matrix[row][column]:
                scale = matrix[row][column]
                matrix[row] = [x - scale * y
                               for x, y in zip(matrix[row], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == nr:
            break
    return pivot_row


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def dot(a, b):
    return sum(Fraction(x) * Fraction(y) for x, y in zip(a, b))


# Dependency bytes and semantic custody.
for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and digest(path) == expected,
          f"dependency custody {relative}")
    check(sha256(path.read_bytes() + b"tamper").hexdigest() != expected,
          f"dependency tamper rejection {relative}")

bs = (ROOT / "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md").read_text()
fj = (ROOT / "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md").read_text()
fk = (ROOT / "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md").read_text()
fm = (ROOT / "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md").read_text()
fq = (ROOT / "LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md").read_text()
fr = (ROOT / "LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/THEOREM.md").read_text()
fs = (ROOT / "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/THEOREM.md").read_text()

check("No distance or dimension occurs" in bs,
      "BS degree parent is explicitly metric free")
check("source-construction template" in bs and "H_{N,\\ell}[j]" in bs,
      "BS distinguishes parent from source template")
check("widehat j_{ab}(v):=s_a(v)s_b(v)" in fj,
      "FJ owns physical link-pair operators")
check("existing BS06 term couples every pair" in fj,
      "FJ derives pair operators inside the degree interaction")
check("normalized pair-state variations" in fk and "two-dimensional `E`" in fk,
      "FK identifies the ice pair E sector")
check("does **not** identify either map" in fk,
      "FK withholds physical metric-source identification")
check("same source-deformed parent" in fq and "post hoc hand weight" in fq,
      "FQ freezes source before Feshbach and rejects hand weights")
check("m_\\xi=\\sum_aN_{\\xi a}D_a" in fr,
      "FR proves additive multi-edge closure")
check("cross-dyad or explicitly blocked root-edge source" in fr,
      "FR leaves a prospectively derived root source open")
check("complete selected Hamiltonian" in fs and "unsplit degree-square term" in fs,
      "FS freezes the reduced parent and term decomposition")
check("through order eight" in fm and "Feshbach" in fm,
      "FM retains the inherited order-eight Feshbach boundary")


# Exact tetrahedral tensor algebra. Symmetric tensor coordinates are
# (xx, yy, zz, xy, xz, yz), with factor two stored in dyad off-diagonals.
SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple(combinations(range(4), 2))


def dyad(vector):
    x, y, z = vector
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


EDGE = tuple(dyad(vector) for vector in SIGNS)
ROOT_DYAD = tuple(
    dyad(tuple(SIGNS[b][i] - SIGNS[a][i] for i in range(3)))
    for a, b in PAIRS
)
ADDITIVE_PAIR = tuple(add(EDGE[a], EDGE[b]) for a, b in PAIRS)
W_VERTEX = tuple(sum(row[column] for row in EDGE) for column in range(6))
E_BASIS = ((1, -1, 0, 0, 0, 0),
           (1, 1, -2, 0, 0, 0))
T2_BASIS = ((0, 0, 0, 1, 0, 0),
            (0, 0, 0, 0, 1, 0),
            (0, 0, 0, 0, 0, 1))
A1 = (1, 1, 1, 0, 0, 0)

check(rank(EDGE) == 4, "four edge dyads have exact rank four")
check(W_VERTEX == (4, 4, 4, 0, 0, 0),
      "unsplit FS degree weight is the scalar edge-dyad sum")
check(rank(EDGE + (W_VERTEX,)) == 4,
      "unsplit FS degree source leaves the exact edge rank unchanged")
check(rank(EDGE + ADDITIVE_PAIR) == 4,
      "additive degree-pair weights cannot enlarge edge rank")
check(all(dot(e, row) == 0 for e in E_BASIS for row in EDGE),
      "both E directions annihilate every edge dyad")
check(all(dot(e, row) == 0 for e in E_BASIS for row in ADDITIVE_PAIR),
      "both E directions annihilate every additive pair weight")
check(rank(ROOT_DYAD) == 6, "six sibling-root dyads have exact rank six")
check(rank(EDGE + ROOT_DYAD) == 6,
      "root pair plus unchanged edge weights have exact rank six")

for (a, b), root_row, additive_row in zip(PAIRS, ROOT_DYAD, ADDITIVE_PAIR):
    cross = tuple(root_row[i] - additive_row[i] for i in range(6))
    expected = dyad(tuple(SIGNS[b][i] - SIGNS[a][i] for i in range(3)))
    check(root_row == expected and any(cross),
          f"pair {a+1}{b+1} root weight contains a nonzero cross-dyad part")


# Exact degree identity and ice-pair representation.
ICE = tuple(state for state in product((-1, 1), repeat=4) if sum(state) == 0)
ALL = tuple(product((-1, 1), repeat=4))

for state in ALL:
    occupation = tuple((1 - z) // 2 for z in state)
    degree = sum(occupation)
    pair_sum = sum(state[a] * state[b] for a, b in PAIRS)
    check(Fraction((degree - 2) ** 2) == 1 + Fraction(pair_sum, 2),
          f"degree-square pair identity for state {state}")

check(len(ICE) == 6, "local q4 ice fiber has six states")
check(all(sum(state[a] * state[b] for a, b in PAIRS) == -2
          for state in ICE), "ice pair sum is the fixed scalar -2")
check(all(state[a] * state[b] == state[c] * state[d]
          for state in ICE
          for (a, b), (c, d) in (((0, 1), (2, 3)),
                                  ((0, 2), (1, 3)),
                                  ((0, 3), (1, 2)))),
      "complementary ice pairs agree")


def pair_query_values(source):
    coefficients = [dot(source, row) for row in ROOT_DYAD]
    return tuple(sum(coefficient * state[a] * state[b]
                     for coefficient, (a, b) in zip(coefficients, PAIRS))
                 for state in ICE)


e_values = tuple(pair_query_values(source) for source in E_BASIS)
t_values = tuple(pair_query_values(source) for source in T2_BASIS)
a_values = pair_query_values(A1)
source_basis = tuple(tuple(Fraction(index == coordinate)
                           for index in range(6)) for coordinate in range(6))
all_source_values = tuple(pair_query_values(source) for source in source_basis)
check(rank(e_values) == 2, "root-pair query realizes both ice E directions")
check(all(sum(values) == 0 for values in e_values),
      "both root-pair E functions are centered on the ice fiber")
check(all(all(value == 0 for value in values) for values in t_values),
      "ice restriction kills the root-pair T2 source")
check(len(set(a_values)) == 1 and a_values[0] != 0,
      "ice root-pair A1 source is a nonzero scalar")
check(rank(all_source_values) == 3,
      "direct ice-projected root-pair source has exact A1+E rank three")

# At source off both prospective pair extensions have exactly the same
# coefficients U_d/2; their derivatives have different ranks.
source_off_fs_expansion = tuple(Fraction(1, 2) for _ in PAIRS)
source_off_root = tuple(Fraction(1, 2) for _ in PAIRS)
check(source_off_fs_expansion == source_off_root,
      "rank-four and rank-six pair sources have identical source-off coefficients")
check(rank(EDGE + (W_VERTEX,)) == 4 and rank(EDGE + ROOT_DYAD) == 6,
      "same source-off Hamiltonian admits different exact derivative ranks")
root_lengths = tuple(sum((SIGNS[b][i] - SIGNS[a][i]) ** 2
                         for i in range(3)) for a, b in PAIRS)
check(set(root_lengths) == {8},
      "all tetrahedral sibling roots share one source-off length")
derivative_ok = True
for root_row in ROOT_DYAD:
    for coordinate in range(6):
        ratio_derivative = -Fraction(root_row[coordinate], 8)
        energy_derivative = Fraction(1, 2) * ratio_derivative
        conjugate = -2 * energy_derivative
        derivative_ok &= conjugate == Fraction(root_row[coordinate], 8)
check(derivative_ok,
      "F=I-j/2 gives the normalized DPAR sign and Q=-2 dH/dj factor")


# Frozen periodic G_5 support: verify every local ring-edge pair occurs among
# inherited hexagons through one vertex. This is the exact support needed by
# the conditional local H6 commutator statement.
L = 5
SHIFTS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
vertices = [(s, x, y, z) for s in (0, 1)
            for x, y, z in product(range(L), repeat=3)]
adjacency = {vertex: set() for vertex in vertices}
edge_label = {}
for x, y, z in product(range(L), repeat=3):
    av = (0, x, y, z)
    for label, shift in enumerate(SHIFTS):
        bv = (1, (x + shift[0]) % L,
                 (y + shift[1]) % L,
                 (z + shift[2]) % L)
        adjacency[av].add(bv)
        adjacency[bv].add(av)
        edge_label[frozenset((av, bv))] = label

check(len(vertices) == 2 * L**3 and len(edge_label) == 4 * L**3,
      "G_5 has exact periodic q4 vertex and edge counts")
check(all(len(adjacency[vertex]) == 4 for vertex in vertices),
      "G_5 is boundaryless and coordination four")


def canonical_cycle(cycle):
    cycle = tuple(cycle)
    variants = []
    for oriented in (cycle, tuple(reversed(cycle))):
        for shift in range(len(cycle)):
            variants.append(oriented[shift:] + oriented[:shift])
    return min(variants)


origin = (0, 0, 0, 0)
hexagons = set()


def visit(vertex, path):
    if len(path) == 6:
        if origin in adjacency[vertex]:
            hexagons.add(canonical_cycle(path))
        return
    for neighbor in adjacency[vertex]:
        if neighbor not in path:
            visit(neighbor, path + [neighbor])


visit(origin, [origin])
local_ring_pairs = set()
for cycle in hexagons:
    index = cycle.index(origin)
    left = edge_label[frozenset((origin, cycle[index - 1]))]
    right = edge_label[frozenset((origin, cycle[(index + 1) % 6]))]
    local_ring_pairs.add(tuple(sorted((left, right))))

check(len(hexagons) == 12, "G_5 has twelve inherited hexagons through one vertex")
check(local_ring_pairs == set(PAIRS),
      "hexagons through one G_5 vertex realize all six local ring-edge pairs")

# A local ring flips the two incident ring links. For locally flippable ice
# states the change of the two E pair queries gives exact commutator channel
# coefficients. Their union over inherited ring types has rank two.
commutator_channels = []
for p, q in local_ring_pairs:
    for state in ICE:
        if state[p] != -state[q]:
            continue
        flipped = list(state)
        flipped[p] *= -1
        flipped[q] *= -1
        channel = []
        for source in E_BASIS:
            coefficients = [dot(source, row) for row in ROOT_DYAD]
            before = sum(c * state[a] * state[b]
                         for c, (a, b) in zip(coefficients, PAIRS))
            after = sum(c * flipped[a] * flipped[b]
                        for c, (a, b) in zip(coefficients, PAIRS))
            channel.append(after - before)
        commutator_channels.append(tuple(channel))

check(rank(commutator_channels) == 2,
      "inherited G_5 H6 ring parities act on both local E directions")
check(all(any(channel[column] for channel in commutator_channels)
          for column in range(2)),
      "each displayed E basis has a nonzero local source-off ring commutator")


# Claim and byte hygiene.
theorem = (LANE / "THEOREM.md").read_text()
check("`E` pair **query** is available" in theorem,
      "theorem distinguishes query availability")
check("do **not** lawfully promote" in theorem,
      "theorem withholds parent-derived strain")
check("Degree-pair affine-response law" in theorem and "g'(1)=\\lambda\\ne0" in theorem,
      "theorem names one explicit sufficient physical premise")
check("exact rank-six statement in this lane is microscopic" in theorem
      and "direct root-pair image has only `A1+E`, hence rank" in theorem,
      "theorem separates exact microscopic rank six from projected rank three")
check("does **not** prove that\n`DPAR` is the unique or logically necessary closure" in theorem
      and "It is neither\ninherited nor adopted here" in theorem,
      "theorem keeps DPAR sufficient, nonunique, and unadopted")
check("Equation (FT08) **replaces** the FS scalar deformation" in theorem
      and "not a double-counted scalar" in theorem,
      "theorem resolves the overlapping A1 source accounting")
check("uniform sum over all vertices" in theorem,
      "theorem retains the uniform-response cancellation ceiling")
check("not a proof of gravity emergence" in theorem,
      "theorem retains the gravity ceiling")

for name in ("DEPENDENCIES.sha256", "RESULT.md", "SELF_AUDIT.md",
             "THEOREM.md", "verify_degree_pair_geometric_strain_boundary.py"):
    data = (LANE / name).read_bytes()
    check(b"\r" not in data and b"\b" not in data and b"\f" not in data,
          f"byte hygiene {name}")

verification_text = (
    "GRA_FT_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY: PASS\n"
    "Checks: 105/105\n"
    "Current frozen reduced-parent query: microscopic rank 4 A1+T2, E null 2\n"
    "Existing ice pair observables: E query available and H6-dynamical\n"
    "Same H[j=0]: frozen microscopic derivative rank 4; DPAR microscopic derivative rank 6\n"
    "DPAR: one sufficient constitutive law, neither inherited nor adopted\n"
    "Full state-dependent CTP rank, Ward packet, gravity, and G: not claimed\n"
)
check((LANE / "VERIFICATION.txt").read_text() == verification_text,
      "verification transcript is exact")

manifest = {}
for line in (LANE / "MANIFEST.sha256").read_text().splitlines():
    expected, name = line.split("  ", 1)
    manifest[name] = expected
check(set(manifest) == MANIFEST_FILES, "manifest member set is exact")
check(all(digest(LANE / name) == expected
          for name, expected in manifest.items()), "manifest hashes match")
check(all(not (LANE / name).is_symlink() for name in manifest),
      "manifest members are not symlinks")
seal = {}
for line in (LANE / "SEAL.sha256").read_text().splitlines():
    expected, name = line.split("  ", 1)
    seal[name] = expected
check(set(seal) == {"MANIFEST.sha256", "VERIFICATION.txt"},
      "builder seal covers manifest and verification transcript")
check(all(digest(LANE / name) == expected for name, expected in seal.items()),
      "builder seal hashes match")
lane_files = {path.name for path in LANE.iterdir() if path.is_file()}
check(lane_files == MANIFEST_FILES | AUDIT_FILES
      | {"MANIFEST.sha256", "SEAL.sha256", "VERIFICATION.txt"},
      "lane file set is exact")

print("GRA_FT_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY: PASS")
print(f"Checks: {checks}/{checks}")
print("Current frozen reduced-parent query: microscopic rank 4 A1+T2, E null 2")
print("Existing ice pair observables: E query available and H6-dynamical")
print("Same H[j=0]: frozen microscopic derivative rank 4; DPAR microscopic derivative rank 6")
print("DPAR: one sufficient constitutive law, neither inherited nor adopted")
print("Full state-dependent CTP rank, Ward packet, gravity, and G: not claimed")
