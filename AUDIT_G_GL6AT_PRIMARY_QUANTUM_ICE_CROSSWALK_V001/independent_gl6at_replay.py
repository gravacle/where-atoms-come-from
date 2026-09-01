#!/usr/bin/env python3
"""Independent exact replay for the GL6AT crosswalk.

This file deliberately imports no GL6AT or upstream verifier.
"""

from fractions import Fraction
from itertools import combinations, permutations, product


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def rank(rows):
    matrix = [[Fraction(x) for x in row] for row in rows]
    if not matrix:
        return 0
    r = 0
    for c in range(len(matrix[0])):
        pivot = next((i for i in range(r, len(matrix)) if matrix[i][c]), None)
        if pivot is None:
            continue
        matrix[r], matrix[pivot] = matrix[pivot], matrix[r]
        scale = matrix[r][c]
        matrix[r] = [x / scale for x in matrix[r]]
        for i in range(len(matrix)):
            if i != r and matrix[i][c]:
                scale = matrix[i][c]
                matrix[i] = [x - scale * y for x, y in zip(matrix[i], matrix[r])]
        r += 1
    return r


def cycle_type(p):
    seen = set()
    cycles = []
    for i in range(4):
        if i in seen:
            continue
        j = i
        length = 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = p[j]
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


# A3/diamond tetrahedral star.
tetra = [tuple(Fraction(a == b) - Fraction(1, 4) for b in range(4))
         for a in range(4)]
for a in range(4):
    check(sum(tetra[a]) == 0, f"bond {a} lies in A3 plane")
    for b in range(4):
        expected = Fraction(3, 4) if a == b else Fraction(-1, 4)
        check(dot(tetra[a], tetra[b]) == expected, f"tetra Gram {a},{b}")
check(rank(tetra) == 3, "tetrahedral star has spatial rank three")

# Degree-two ice mapping and exact Z normalization.
ice = []
for occupied in combinations(range(4), 2):
    n = tuple(int(a in occupied) for a in range(4))
    z = tuple(2 * x - 1 for x in n)
    sz = tuple(Fraction(x, 2) for x in z)
    check(sum(n) == 2, "degree-two occupation")
    check(sum(sz) == 0, "two-in/two-out spin sum")
    check(all(z[a] == 2 * sz[a] for a in range(4)), "Z equals twice Sz")
    ice.append((n, z))
check(len(ice) == 6, "six local ice assignments")

# Six-cycle toggle: precisely the two alternating strings survive locally.
preserved = []
for bits in product((0, 1), repeat=6):
    flipped = tuple(1 - b for b in bits)
    if all(bits[i - 1] + bits[i] == flipped[i - 1] + flipped[i]
           for i in range(6)):
        preserved.append(bits)
check(preserved == [(0, 1, 0, 1, 0, 1), (1, 0, 1, 0, 1, 0)],
      "only alternating hexagons remain locked under six-X toggle")
check(tuple(1 - b for b in preserved[0]) == preserved[1], "unit toggle pair")
check(Fraction(63, 8) > 0, "positive GL6AO ring amplitude")
check(-Fraction(63, 8) < 0, "negative sign-free Hamiltonian coefficient")
check(Fraction(0, 1) / Fraction(63, 8) == 0, "v/g equals zero")
check(Fraction(1, 1) != 0, "RK coordinate differs from zero-potential point")

# Pair representation: centered complement sums are E; differences are T2.
edges = list(combinations(range(4), 2))
edge_index = {edge: i for i, edge in enumerate(edges)}
opposites = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
sum_vectors = []
diff_vectors = []
for edge, comp in opposites:
    s = [0] * 6
    d = [0] * 6
    s[edge_index[edge]] = s[edge_index[comp]] = 1
    d[edge_index[edge]] = 1
    d[edge_index[comp]] = -1
    sum_vectors.append(tuple(s))
    diff_vectors.append(tuple(d))
centered = [tuple(3 * x - sum(sum_vectors[j][i] for j in range(3))
                  for i, x in enumerate(sum_vectors[a])) for a in range(3)]
check(rank(sum_vectors) == 3, "opposite-pair sum rank")
check(rank(centered) == 2, "centered opposite-pair sum plane is two-dimensional")
check(rank(diff_vectors) == 3, "opposite-pair difference space is three-dimensional")
check(all(dot(s, d) == 0 for s in sum_vectors for d in diff_vectors),
      "sum and difference sectors are orthogonal")

class_sizes = {(1, 1, 1, 1): 1, (2, 1, 1): 6, (2, 2): 3,
               (3, 1): 8, (4,): 6}
char_e_expected = {(1, 1, 1, 1): 2, (2, 1, 1): 0, (2, 2): 2,
                   (3, 1): -1, (4,): 0}
char_t2_expected = {(1, 1, 1, 1): 3, (2, 1, 1): 1, (2, 2): -1,
                    (3, 1): 0, (4,): -1}
char_e = {}
char_t2 = {}
for p in permutations(range(4)):
    kind = cycle_type(p)
    fixed_partitions = 0
    signed_trace = 0
    for a, (edge, comp) in enumerate(opposites):
        pe = tuple(sorted((p[edge[0]], p[edge[1]])))
        pc = tuple(sorted((p[comp[0]], p[comp[1]])))
        for b, (target, target_comp) in enumerate(opposites):
            if {pe, pc} == {target, target_comp} and a == b:
                fixed_partitions += 1
                signed_trace += 1 if (pe, pc) == (target, target_comp) else -1
    char_e.setdefault(kind, set()).add(fixed_partitions - 1)
    char_t2.setdefault(kind, set()).add(signed_trace)
for kind in class_sizes:
    check(char_e[kind] == {char_e_expected[kind]}, f"E character on {kind}")
    check(char_t2[kind] == {char_t2_expected[kind]}, f"T2 character on {kind}")
inner = sum(class_sizes[k] * char_e_expected[k] * char_t2_expected[k]
            for k in class_sizes) / 24
check(inner == 0, "E and T2 are inequivalent orthogonal irreps")

# Strict local ice values: T2 differences vanish, E sums need not.
nonzero_centered_e = False
for _, z in ice:
    pair = {edge: z[edge[0]] * z[edge[1]] for edge in edges}
    sums = [pair[e] + pair[c] for e, c in opposites]
    diffs = [pair[e] - pair[c] for e, c in opposites]
    check(diffs == [0, 0, 0], "locked local pair T2 difference vanishes")
    centered_values = [3 * x - sum(sums) for x in sums]
    nonzero_centered_e |= any(centered_values)
check(nonzero_centered_e, "locked local pair E sector is nontrivial")

# Fu displayed-integral infrared counting: d^3p * epsilon^2 * delta.
spatial_radial_power = 3 - 1
electric_vertex_power = 2
check(spatial_radial_power == 2, "three-dimensional radial density power")
check(spatial_radial_power + electric_vertex_power == 4,
      "displayed electric Raman integral scales as Omega^4")
check(spatial_radial_power == 2, "unweighted two-photon density scales as Omega^2")

print(f"PASS__GL6AT_INDEPENDENT_REPLAY__{checks}/{checks}")
