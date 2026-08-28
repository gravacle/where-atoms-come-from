#!/usr/bin/env python3
"""Algorithmically independent hostile checks for the repaired FZ theorem.

This verifier does not import the FZ author verifier or FY source machinery.
It imports only the frozen FO finite parent, reconstructs the pair source and
one complete 720+720 ring entry from the microscopic support rules, and then
checks the exact TT and finite Liouvillian/Ward algebra independently.
"""

from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
from itertools import permutations
from pathlib import Path
import runpy


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent
FO_SCRIPT = (ROOT / "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001" /
             "verify_finite_tt_four_point.py")

checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def verify_hash_list(list_path, base):
    count = 0
    for line in list_path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        path = base / relative
        check(path.is_file() and not path.is_symlink() and
              digest(path) == expected,
              f"custody {list_path.name}: {relative}")
        count += 1
    return count


# Freeze the exact repaired bytes audited here.
check(verify_hash_list(AUDIT / "TARGET_CUSTODY.sha256", LANE) == 6,
      "target custody contains the six repaired FZ core files")

# Independently replay FZ's declared direct parent custody.
check(verify_hash_list(LANE / "DEPENDENCIES.sha256", ROOT) == 6,
      "FZ direct dependency list has six hash-pinned files")


def verify_parent_seal(name):
    parent = ROOT / name
    manifest = parent / "MANIFEST.sha256"
    seal = parent / "SEAL.sha256"
    count = verify_hash_list(manifest, parent)
    seal_lines = seal.read_text().splitlines()
    expected_manifest, manifest_name = seal_lines[0].split("  ", 1)
    check(manifest_name == "MANIFEST.sha256" and
          expected_manifest == digest(manifest),
          f"{name}: seal owns the replayed manifest")
    return count


for parent_name in (
    "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001",
    "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001",
):
    check(verify_parent_seal(parent_name) >= 5,
          f"{parent_name}: sealed parent packet replays")


# Read the repaired claim surfaces before doing any physics inference.
theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
verifier = (LANE / "verify_continuity_contact_ward_boundary.py").read_text()
theorem_flat = " ".join(theorem.split())
result_flat = " ".join(result.split())
self_audit_flat = " ".join(self_audit.split())

check("admits a zero-current representation, but it does not\nselect one" in
      theorem and "divergence-free or circulating projected currents" in
      theorem,
      "theorem does not infer a unique or absent current from zero density")
check("does not exclude divergence-free/circulating currents" in result_flat and
      "no bond current is derived" in result_flat,
      "result preserves the zero-current-representation boundary")
check("physical discrete divergence has not been derived" in theorem_flat and
      "future `Delta_m` need not equal `k_i`" in theorem_flat,
      "the supplied embedding contraction is not typed as physical divergence")
check("under the supplied embedding contraction" in self_audit_flat and
      "unknown physical divergence is retained as undecided" in self_audit_flat,
      "self-audit retains the repaired Ward-scope boundary")
check("full Ward identity remains undecided rather than failed" in result_flat,
      "result does not promote nontransversality to Ward failure")
check("no derived \"\n      \"physical discrete divergence" in verifier and
      "Ward \"\n      \"closure, continuum, gravity, or G" in verifier,
      "author verifier prints the narrowed physical ceiling")


# Load only the frozen finite graph/component.  All source reconstruction
# below is new code and does not call FY or FZ source helpers.
with redirect_stdout(StringIO()):
    fo = runpy.run_path(str(FO_SCRIPT))

states = tuple(fo["states"])
state_index = dict(fo["state_index"])
edges = tuple(tuple(edge) for edge in fo["edges"])
edge_labels = tuple(fo["edge_labels"])
incidence = tuple(tuple(fo["incidence"][vertex]) for vertex in range(60))
ring_patterns = tuple(fo["ring_patterns"])
hexagons = tuple(tuple(row) for row in fo["hexagons"])
translation_orbits = tuple(tuple(row) for row in fo["translation_orbits"])

CELL_COUNT = 30
VERTEX_COUNT = 60
EDGE_COUNT = 120
SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple((a, b) for a in range(4) for b in range(a + 1, 4))
SUPPORT_EXPONENT = (0, 10, 5, 9, 25, 201)
ZERO6 = (F(0),) * 6

check((len(states), len(edges), len(incidence)) ==
      (180, EDGE_COUNT, VERTEX_COUNT),
      "independent replay uses the exact FO 180/120/60 finite parent")


def degree_row(state):
    return tuple(sum((state >> edge) & 1 for _, edge in incidence[vertex])
                 for vertex in range(VERTEX_COUNT))


degree_rows = tuple(degree_row(state) for state in states)
check(all(row == (2,) * VERTEX_COUNT for row in degree_rows),
      "all projected incidence-charge eigenvalues d_v-2 vanish exactly")

# Rebuild the ring-transition graph from masks rather than trusting FO's H.
transition_cycle = {}
transition_hits = 0
for row, state in enumerate(states):
    for (mask, first, second), cycle in zip(ring_patterns, hexagons):
        if (state & mask) not in (first, second):
            continue
        target = state ^ mask
        if target not in state_index:
            raise AssertionError("ring endpoint left C_180")
        transition_hits += 1
        key = tuple(sorted((row, state_index[target])))
        previous = transition_cycle.setdefault(key, cycle)
        if set(previous) != set(cycle):
            raise AssertionError("ring transition has inconsistent owner")

check(transition_hits == 840,
      "all 840 directed ring hits remain in C_180 with consistent owners")
check(len(transition_cycle) == 420,
      "independent ring inventory has 420 undirected H6 transitions")
check(all(degree_rows[row] == degree_rows[column]
          for row, column in transition_cycle),
      "every reconstructed H6 transition preserves every projected charge")


# -------------------------------------------------------------------------
# Independent exact Q(zeta_240) arithmetic and pair-source reconstruction.


def dyad(vector):
    x, y, z = vector
    return (x * x, y * y, z * z,
            2 * x * y, 2 * x * z, 2 * y * z)


def add_rows(*rows):
    return tuple(sum(values, F(0)) for values in zip(*rows))


def scale_row(scale, row):
    return tuple(F(scale) * value for value in row)


def ledger_add(ledger, key, row, scale=F(1)):
    updated = add_rows(ledger.get(key, ZERO6), scale_row(scale, row))
    if any(updated):
        ledger[key] = updated
    elif key in ledger:
        del ledger[key]


def vertex_support(vertex):
    return (0, vertex) if vertex < CELL_COUNT else (1, vertex - CELL_COUNT)


def edge_support(edge):
    return (2 + edge_labels[edge], edge // 4)


def vertex_z(state):
    result = []
    for vertex in range(VERTEX_COUNT):
        row = [None] * 4
        for _, edge in incidence[vertex]:
            row[edge_labels[edge]] = 1 - 2 * ((state >> edge) & 1)
        if not all(value in (-1, 1) for value in row):
            raise AssertionError("q4 vertex lacks one signed value per label")
        result.append(tuple(row))
    return tuple(result)


def direct_pair_ledger(state):
    ledger = {}
    for vertex, z in enumerate(vertex_z(state)):
        row = ZERO6
        for a, b in PAIRS:
            root = tuple(SIGNS[b][axis] - SIGNS[a][axis]
                         for axis in range(3))
            # lambda=-1/2 and Rhat=dyad(root)/8.
            row = add_rows(row,
                           scale_row(F(-z[a] * z[b], 16), dyad(root)))
        ledger[vertex_support(vertex)] = row
    return ledger


PHI240 = [F(0)] * 65
for power, coefficient in ((0, 1), (8, 1), (24, -1), (32, -1),
                           (40, -1), (56, 1), (64, 1)):
    PHI240[power] = F(coefficient)
PHI240 = tuple(PHI240)


def trim(polynomial):
    result = list(polynomial)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def reduce_phi(polynomial):
    result = trim(tuple(F(value) for value in polynomial))
    while len(result) >= len(PHI240):
        shift = len(result) - len(PHI240)
        factor = result[-1]
        for index, value in enumerate(PHI240):
            result[index + shift] -= factor * value
        result = trim(result)
    return tuple(result + [F(0)] * (64 - len(result)))


x240_minus_one = [F(0)] * 241
x240_minus_one[0] = F(-1)
x240_minus_one[240] = F(1)
check(not any(reduce_phi(x240_minus_one)),
      "hard-coded degree-64 Phi_240 exactly divides z^240-1")


def component_polynomial(ledger, component):
    polynomial = [F(0)] * 240
    for (support, cell), row in ledger.items():
        exponent = (8 * cell + SUPPORT_EXPONENT[support]) % 240
        polynomial[exponent] += row[component]
    return tuple(polynomial)


def linear_polynomial(ledger, covector):
    polynomial = [F(0)] * 240
    for component, coefficient in enumerate(covector):
        if not coefficient:
            continue
        source = component_polynomial(ledger, component)
        for power, value in enumerate(source):
            polynomial[power] += F(coefficient) * value
    return tuple(polynomial)


LONGITUDINAL_2 = (
    (14, 0, 0, 15, -17, 0),
    (0, 30, 0, 7, 0, -17),
    (0, 0, -34, 0, 7, 15),
)

pair_ledger = direct_pair_ledger(translation_orbits[0][0])
pair_remainders = tuple(reduce_phi(linear_polynomial(pair_ledger, covector))
                        for covector in LONGITUDINAL_2)
pair_terms = tuple((power, value) for power, value in
                   enumerate(pair_remainders[0]) if value)
EXPECTED_PAIR_TERMS = (
    (0, F(-14)), (2, F(28)), (8, F(14)), (10, F(14)),
    (24, F(14)), (26, F(-14)), (34, F(-28)), (42, F(-14)),
    (48, F(-14)), (56, F(-14)), (58, F(28)),
)
check(pair_terms == EXPECTED_PAIR_TERMS,
      "independent pair construction reproduces exact FZ05 remainder")


# -------------------------------------------------------------------------
# Independent naive 720+720 ring-history derivative reconstruction.


def local_eprime8(z):
    total = [0] * 6
    for a, b in PAIRS:
        difference = tuple(SIGNS[b][axis] - SIGNS[a][axis]
                           for axis in range(3))
        row = dyad(difference)
        for component, value in enumerate(row):
            total[component] += z[a] * z[b] * value // 4
    return tuple(total)


def prefix_gap_data(state, z_initial, toggled):
    degree_delta = {}
    affected = set()
    for edge in toggled:
        occupation_delta = -1 if (state >> edge) & 1 else 1
        for vertex in edges[edge]:
            affected.add(vertex)
            degree_delta[vertex] = (degree_delta.get(vertex, 0) +
                                    occupation_delta)
    gap = sum(value * value for value in degree_delta.values())
    if gap <= 0:
        raise AssertionError("proper ring prefix returned to ice")
    gp8 = {}
    for vertex in affected:
        before = z_initial[vertex]
        after = list(before)
        for edge in toggled:
            if vertex in edges[edge]:
                after[edge_labels[edge]] *= -1
        row = tuple(new - old for new, old in
                    zip(local_eprime8(tuple(after)),
                        local_eprime8(before)))
        if any(row):
            gp8[vertex_support(vertex)] = row
    return F(gap), gp8


def naive_oriented_ring_derivative(state, cycle):
    z_initial = vertex_z(state)
    total = {}
    path_count = 0
    for order in permutations(tuple(cycle)):
        path_count += 1
        path = []
        for prefix_length in range(1, 6):
            path.append(prefix_gap_data(
                state, z_initial, order[:prefix_length]))
        p0 = F(-1)
        for gap, _ in path:
            p0 /= gap

        # d product(hop)/dj = -sum D_label/2 = -sum dyad/6.
        for edge in order:
            ledger_add(total, edge_support(edge),
                       dyad(SIGNS[edge_labels[edge]]), F(-1, 6) * p0)

        # d[-prod(1/g)] = p0[-sum (dg/g)]; gp8 stores 8*dg.
        for gap, gp8 in path:
            for support, row in gp8.items():
                ledger_add(total, support, row, -p0 / (8 * gap))
    check(path_count == 720,
          "one oriented six-link ring sums all 6!=720 histories")
    return total


ring_key = min(transition_cycle)
ring_cycle = transition_cycle[ring_key]
forward = naive_oriented_ring_derivative(states[ring_key[0]], ring_cycle)
reverse = naive_oriented_ring_derivative(states[ring_key[1]], ring_cycle)
ring_ledger = {}
for support, row in forward.items():
    ledger_add(ring_ledger, support, row, F(-8, 63))
for support, row in reverse.items():
    ledger_add(ring_ledger, support, row, F(-8, 63))

missing_label = next(iter(set(range(4)) -
                          {edge_labels[edge] for edge in ring_cycle}))
expected_ring_sum = add_rows(
    (F(-31, 6), F(-31, 6), F(-31, 6), F(0), F(0), F(0)),
    scale_row(F(3, 2), dyad(SIGNS[missing_label])),
)
check(add_rows(*ring_ledger.values()) == expected_ring_sum,
      "naive 720+720 derivative recovers the frozen homogeneous ring row")

ring_remainders = tuple(reduce_phi(linear_polynomial(ring_ledger, covector))
                        for covector in LONGITUDINAL_2)
ring_first = next(remainder for remainder in ring_remainders if any(remainder))
ring_terms = tuple((power, value) for power, value in enumerate(ring_first)
                   if value)
EXPECTED_RING_TERMS = (
    (2, F(-20, 3)), (10, F(-20, 3)), (13, F(10, 3)),
    (17, F(-32, 3)), (18, F(-25, 3)), (25, F(-50, 3)),
    (34, F(20, 3)), (41, F(50, 3)), (42, F(20, 3)),
    (48, F(20, 3)), (49, F(50, 3)), (53, F(-10, 3)),
    (56, F(25, 3)), (57, F(14, 3)), (58, F(5, 3)),
    (61, F(-10, 3)),
)
check(ring_terms == EXPECTED_RING_TERMS,
      "independent ring construction reproduces exact ring obstruction")


# The sample coefficients are exact and the diagonal/off-diagonal split makes
# cancellation impossible entry by entry.
for x, rho, expected_f, expected_product in (
    (F(2, 5), F(15625, 504), F(2415673, 3515625),
     F(2415673, 113400)),
    (F(1, 2), F(512, 63), F(15853, 57600), F(31706, 14175)),
):
    f_e = (1 - x**2 - F(37, 12) * x**4 -
           F(16247, 900) * x**6)
    check(f_e == expected_f and rho * f_e == expected_product,
          f"x={x}: exact diagonal coefficient matches FZ07 and is nonzero")


# -------------------------------------------------------------------------
# Exact TT projector, with a separately constructed rational image basis.


def matmul(first, second):
    return tuple(tuple(sum(first[row][inner] * second[inner][column]
                           for inner in range(len(second)))
                       for column in range(len(second[0])))
                 for row in range(len(first)))


def matadd(first, second):
    return tuple(tuple(first[row][column] + second[row][column]
                       for column in range(len(first[0])))
                 for row in range(len(first)))


def matscale(scale, matrix):
    return tuple(tuple(F(scale) * value for value in row) for row in matrix)


def outer(first, second):
    return tuple(tuple(F(a) * F(b) for b in second) for a in first)


def rank(matrix):
    work = [list(map(F, row)) for row in matrix]
    row = 0
    for column in range(len(work[0])):
        pivot = next((candidate for candidate in range(row, len(work))
                      if work[candidate][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        divisor = work[row][column]
        work[row] = [value / divisor for value in work[row]]
        for other in range(len(work)):
            if other == row or not work[other][column]:
                continue
            factor = work[other][column]
            work[other] = [left - factor * right
                           for left, right in zip(work[other], work[row])]
        row += 1
    return row


r = tuple(map(F, (7, 15, -17)))
identity3 = tuple(tuple(F(int(i == j)) for j in range(3)) for i in range(3))
P = matadd(identity3, matscale(F(-1, 563), outer(r, r)))
check(matmul(P, P) == P and
      all(sum(P[i][j] * r[j] for j in range(3)) == 0 for i in range(3)) and
      sum(P[i][i] for i in range(3)) == 2,
      "independent exact transverse projector obeys P^2=P, Pr=0, trP=2")


def unpack(column):
    values = [F(0)] * 6
    values[column] = F(1)
    return ((values[0], values[3] / 2, values[4] / 2),
            (values[3] / 2, values[1], values[5] / 2),
            (values[4] / 2, values[5] / 2, values[2]))


def pack(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            2 * matrix[0][1], 2 * matrix[0][2], 2 * matrix[1][2])


def tt(matrix):
    pap = matmul(matmul(P, matrix), P)
    trace_pa = sum(P[i][j] * matrix[i][j]
                   for i in range(3) for j in range(3))
    return matadd(pap, matscale(-trace_pa / 2, P))


columns = tuple(pack(tt(unpack(column))) for column in range(6))
TT = tuple(tuple(columns[column][row] for column in range(6))
           for row in range(6))
check(matmul(TT, TT) == TT and rank(TT) == 2,
      "independent six-coordinate TT projector is idempotent and rank two")

# Construct two exact TT tensors without using the projector columns.
u = tuple(map(F, (15, -7, 0)))
v = tuple(map(F, (-119, -255, -274)))  # r cross u
A = matadd(outer(u, u), matscale(F(-1, 563), outer(v, v)))
B = matadd(outer(u, v), outer(v, u))
check(tt(A) == A and tt(B) == B and rank((pack(A), pack(B))) == 2,
      "two independently constructed exact tensors span the TT image")


# -------------------------------------------------------------------------
# Degeneracy-aware inverse Liouvillian and retarded-sign checks.


energies = (0, 0, 2, 5)
L = [[0j for _ in energies] for _ in energies]
for a, b, value in ((0, 2, 3), (1, 3, 2), (2, 3, -4)):
    L[a][b] = complex(value)
    L[b][a] = complex(value)
P_solution = [[0j for _ in energies] for _ in energies]
for a, energy_a in enumerate(energies):
    for b, energy_b in enumerate(energies):
        if energy_a != energy_b:
            P_solution[a][b] = 1j * L[a][b] / (energy_a - energy_b)
check(all(abs(1j * (energies[a] - energies[b]) * P_solution[a][b] +
                  L[a][b]) == 0
          for a in range(4) for b in range(4)),
      "FZ11 has the correct sign and solves i[H,P]+L=0 off energy blocks")

L_blocked = [row[:] for row in L]
L_blocked[0][1] = L_blocked[1][0] = 1
check(any(L_blocked[a][b] != 0 and energies[a] == energies[b]
          for a in range(4) for b in range(4)),
      "a nonzero complete degenerate-energy block is an exact ad_H obstruction")


def spectral_chi(Aop, Bop, omega):
    gap = 2.0
    return (Aop[0][1] * Bop[1][0] / (omega - gap) -
            Bop[0][1] * Aop[1][0] / (omega + gap))


Ptoy = ((0j, 1 + 0j), (1 + 0j, 0j))
Btoy = ((0j, -1j), (1j, 0j))
Ltoy = ((0j, 2j), (-2j, 0j))  # -i[diag(0,2),P]
omega = 0.7 + 0.3j
chi_p = spectral_chi(Ptoy, Btoy, omega)
chi_l = spectral_chi(Ltoy, Btoy, omega)
commutator_ground = (Ptoy[0][1] * Btoy[1][0] -
                     Btoy[0][1] * Ptoy[1][0])
check(abs((chi_l - 1j * omega * chi_p) -
          (-1j * commutator_ground)) < 1e-12,
      "FZ12 retarded sign matches e^{+i omega t} and the equal-time contact")


print("PAIR_WITNESS", [(power, str(value)) for power, value in pair_terms])
print("RING_WITNESS", [(power, str(value)) for power, value in ring_terms])
print(f"SUMMARY {checks}/{checks} independent hostile FZ checks passed")
print("VERDICT PASS")
print("CEILING exact projected-charge identity, supplied-embedding "
      "longitudinal diagnostics, exact TT algebra, and algebraic Ward "
      "dependency only; no physical discrete divergence, current, Ward "
      "closure, continuum, gravity, or G")
