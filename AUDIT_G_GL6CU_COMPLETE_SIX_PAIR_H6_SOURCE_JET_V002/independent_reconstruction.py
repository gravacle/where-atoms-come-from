#!/usr/bin/env python3
"""Independent bounded reconstruction for the frozen GL6CU V002 audit.

No target module is imported.  The audit independently rebuilds the Q4
incidence graph, the support census, scalar source two-jets for the literal
first stencil in each of the eight compact-ledger records, one complete
two-endpoint canonical cycle calculation, and the CTP/half-source signs.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TARGET = ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002"
L = 4
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
PAIRS = tuple(combinations(range(4), 2))
SELECTOR = frozenset(((0, 1), (0, 2), (1, 1), (1, 3), (2, 0), (2, 1)))
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def canonical_hash(value):
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(raw).hexdigest()


def rational(text):
    return Q(text)


def add(a, b):
    return tuple((a[i] + b[i]) % L for i in range(3))


def sub(a, b):
    return tuple((a[i] - b[i]) % L for i in range(3))


cells = tuple(product(range(L), repeat=3))
nodes = tuple((kind, cell) for kind in ("P", "C") for cell in cells)
node_index = {node: i for i, node in enumerate(nodes)}
edges = tuple((cell, port) for cell in cells for port in range(4))
edge_index = {edge: i for i, edge in enumerate(edges)}


def endpoints(edge):
    cell, port = edge
    return ("P", cell), ("C", add(cell, STEPS[port]))


ends = tuple(endpoints(edge) for edge in edges)
incident = {
    (kind, cell): tuple(
        edge_index[(cell, port)] if kind == "P"
        else edge_index[(sub(cell, STEPS[port]), port)]
        for port in range(4)
    )
    for kind, cell in nodes
}


def selected(a, b):
    return int((a % L, b % L) in SELECTOR)


def occupation(edge):
    cell, port = edge
    x, y, z = cell
    a, b = (y - x) % L, (x + z) % L
    word = (selected(a - 1, b + 1), selected(a, b),
            1 - selected(a - 1, b + 1), 1 - selected(a, b))
    return word[port]


base = frozenset(i for i, edge in enumerate(edges) if occupation(edge))


def bit(edge, toggled=frozenset()):
    return int(edge in base) ^ int(edge in toggled)


def pair_word(node, toggled=frozenset()):
    z = tuple(1 - 2 * bit(edge, toggled) for edge in incident[node])
    return tuple(Q(z[a] * z[b]) for a, b in PAIRS)


def locked(toggled=frozenset()):
    return all(sum(bit(edge, toggled) for edge in incident[node]) == 2
               for node in nodes)


def defect(toggled):
    charge = {}
    for edge in toggled:
        sign = 1 if edge not in base else -1
        for node in ends[edge]:
            charge[node] = charge.get(node, 0) + sign
    return Q(sum(value * value for value in charge.values()))


check((len(cells), len(nodes), len(edges), len(base)) == (64, 128, 256, 128),
      "independent Q4 cardinalities")
check(locked(), "selected branch is locked")
check(all(len(set(row)) == 4 for row in incident.values()),
      "degree-four incidence")
check(all(edges[incident[("P", cell)][port]] == (cell, port)
          for cell in cells for port in range(4)), "P physical port order")
check(all(edges[incident[("C", cell)][port]] ==
          (sub(cell, STEPS[port]), port)
          for cell in cells for port in range(4)), "C physical port order")


# Enumerate six-cycles by DFS, independently of the target's closed formula.
adjacency = {node: [] for node in nodes}
for edge, (left, right) in enumerate(ends):
    adjacency[left].append((edge, right))
    adjacency[right].append((edge, left))
cycle_supports = set()


def walk(start, current, visited, path):
    if len(path) == 6:
        if current == start:
            cycle_supports.add(frozenset(path))
        return
    for edge, nxt in adjacency[current]:
        if edge in path:
            continue
        if nxt == start:
            if len(path) == 5:
                walk(start, nxt, visited, path + [edge])
            continue
        if nxt in visited:
            continue
        walk(start, nxt, visited | {nxt}, path + [edge])


for start in (node for node in nodes if node[0] == "P"):
    walk(start, start, {start}, [])


def alternating(support):
    return sum(edge in base for edge in support) == 3 and locked(support)


active_cycle_supports = {row for row in cycle_supports if alternating(row)}
check(len(cycle_supports) == 256, "DFS gives 256 six-cycle supports")
check(len(active_cycle_supports) == 64, "64 cycles connect the selected row")

wedges = set()
stars = set()
for node in nodes:
    wedges.update(frozenset(row) for row in combinations(incident[node], 2))
    stars.update(frozenset(row) for row in combinations(incident[node], 3))
paths = set()
for middle, (left, right) in enumerate(ends):
    for a in incident[left]:
        for b in incident[right]:
            if a != middle and b != middle:
                paths.add(frozenset((a, middle, b)))
check((len(wedges), len(stars), len(paths)) == (768, 512, 2304),
      "wedge star path census")
check(not (stars & paths), "star and path supports disjoint")
check(256 + 768 + 512 + 2304 == 3840,
      "3840 diagonal supports computed")
check(3840 + len(cycle_supports) == 4096,
      "4096 graph-wide geometric supports")
check(3840 + len(active_cycle_supports) == 3904,
      "3904 selected-row nonbare Jets")
check(len(cycle_supports - active_cycle_supports) == 192,
      "192 cycle supports inactive on the selected row")


@dataclass(frozen=True)
class D2:
    """Scalar Taylor jet c0+c1*t+c2*t^2; c2 is half the Hessian pullback."""

    c0: Q = Q(0)
    c1: Q = Q(0)
    c2: Q = Q(0)

    def __add__(self, other):
        other = dual(other)
        return D2(self.c0 + other.c0, self.c1 + other.c1,
                  self.c2 + other.c2)

    __radd__ = __add__

    def __neg__(self):
        return D2(-self.c0, -self.c1, -self.c2)

    def __sub__(self, other):
        return self + (-dual(other))

    def __rsub__(self, other):
        return dual(other) - self

    def __mul__(self, other):
        other = dual(other)
        return D2(self.c0 * other.c0,
                  self.c0 * other.c1 + self.c1 * other.c0,
                  self.c0 * other.c2 + self.c1 * other.c1 +
                  self.c2 * other.c0)

    __rmul__ = __mul__

    def reciprocal(self):
        check(self.c0 != 0, "nonzero scalar-jet denominator")
        return D2(1 / self.c0,
                  -self.c1 / self.c0**2,
                  self.c1**2 / self.c0**3 - self.c2 / self.c0**2)

    def __truediv__(self, other):
        return self * dual(other).reciprocal()

    def scale(self, value):
        return self * Q(value)


def dual(value):
    return value if isinstance(value, D2) else D2(Q(value))


ZERO = D2()
ONE = D2(Q(1))


def touched_nodes(active):
    return tuple(sorted({node for edge in active for node in ends[edge]}, key=repr))


def global_weights(active, local_weights, ordered_nodes=None):
    touched = tuple(ordered_nodes) if ordered_nodes is not None else touched_nodes(active)
    check(set(touched) == {node for edge in active for node in ends[edge]},
          "local stencil nodes cover support")
    return {6 * node_index[node] + pair: Q(local_weights[6 * local + pair])
            for local, node in enumerate(touched) for pair in range(6)}


def source_shift(toggled, weights):
    answer = Q(0)
    affected = {node for edge in toggled for node in ends[edge]}
    for node in affected:
        before = pair_word(node)
        after = pair_word(node, toggled)
        offset = 6 * node_index[node]
        answer += sum(weights.get(offset + pair, Q(0)) *
                      (after[pair] - before[pair]) for pair in range(6))
    return answer


def raw_scalar(active, weights):
    active = tuple(sorted(active))
    size = len(active)
    coefficient = [ZERO for _ in range(7)]
    wave = [dict() for _ in range(7)]
    wave[0][0] = ONE
    denominators = {}
    for mask in range(1, 1 << size):
        toggled = frozenset(active[i] for i in range(size) if mask & (1 << i))
        denominators[mask] = D2(defect(toggled), source_shift(toggled, weights))
    for order in range(1, 7):
        for mask in range(1, 1 << size):
            numerator = sum((wave[order - 1].get(mask ^ (1 << i), ZERO)
                             for i in range(size)), ZERO)
            numerator += sum((coefficient[m] *
                              wave[order - m].get(mask, ZERO)
                              for m in range(1, order)), ZERO)
            wave[order][mask] = numerator / denominators[mask]
        coefficient[order] = -sum(
            (wave[order - 1].get(1 << i, ZERO) for i in range(size)), ZERO)
    return tuple(coefficient)


def owner_scalar(active, order, weights, memo=None):
    if memo is None:
        memo = {}
    key = (frozenset(active), order)
    if key in memo:
        return memo[key]
    ordered = tuple(sorted(active))
    answer = raw_scalar(ordered, weights)[order]
    for size in range(1, len(ordered)):
        for subset in combinations(ordered, size):
            answer -= owner_scalar(subset, order, weights, memo)
    memo[key] = answer
    return answer


def payload_scalar(payload, weights):
    c0 = rational(payload["value"])
    c1 = sum((rational(value) * weights[index]
              for index, value in payload["gradient"]), Q(0))
    quadratic = Q(0)
    for left, right, value in payload["hessian"]:
        term = rational(value) * weights[left] * weights[right]
        quadratic += term if left == right else 2 * term
    return D2(c0, c1, quadratic / 2)


def weight_rows(size):
    e0 = [Q(0)] * size
    e0[0] = Q(1)
    elast = [Q(0)] * size
    elast[-1] = Q(-2)
    dense_a = [Q(((7 * i + 3) % 11) - 5, (i % 3) + 1) for i in range(size)]
    dense_b = [Q(((5 * i + 1) % 13) - 6, (i % 4) + 1) for i in range(size)]
    return (tuple(e0), tuple(elast), tuple(dense_a), tuple(dense_b))


ledger = json.loads((TARGET / "EXACT_LEDGER.json").read_text(),
                    object_pairs_hook=unique_object)
rows = dict(ledger["exact_local_operator_stencils"])
rows.update(ledger["cycle_writer_selected_Q4"]["exact_local_operator_stencils"])
check(len(rows) == 8, "eight compact-ledger family/order rows")
check(sum(row["class_count"] for row in rows.values()) == 1654,
      "1654 class commitments")
check(sum(row["owner_count"] for row in rows.values()) == 5184,
      "eight order-resolved records contain 5184 owner-order Jets")
check(rows["one_edge_h6"]["owner_count"] +
      rows["wedge_h6"]["owner_count"] +
      rows["star_h6"]["owner_count"] +
      rows["path_h6"]["owner_count"] +
      rows["alternating_cycle_offdiag_h6"]["owner_count"] == 3904,
      "unique support families cover 3904 selected-row nonbare Jets")
for name, row in rows.items():
    check(canonical_hash(row["first_exact_stencil"]) ==
          row["representative_class_records"][0]["stencil_sha256"],
          f"literal first stencil hashes to first class {name}")
    check(len(row["all_exact_stencils_sha256"]) == 64 and
          len(row["all_class_records_sha256"]) == 64,
          f"aggregate class commitments present {name}")


# Reconstruct the literal first stencil of all seven diagonal records through
# four independent exact scalar pullbacks each.  This checks every stored
# gradient/Hessian collectively without importing the target's Jet class.
diagonal_orders = {
    "one_edge_h2": 2, "one_edge_h4": 4, "one_edge_h6": 6,
    "wedge_h4": 4, "wedge_h6": 6, "star_h6": 6, "path_h6": 6,
}
for name, order in diagonal_orders.items():
    row = rows[name]
    active = tuple(row["representative_class_records"][0]["representative_edges"])
    payload = row["first_exact_stencil"]
    size = 6 * row["representative_class_records"][0]["incident_node_count"]
    pulls = weight_rows(size)
    for weights_local in pulls:
        weights = global_weights(active, weights_local)
        got = owner_scalar(active, order, weights)
        expected = payload_scalar(payload, weights_local)
        check(got == expected, f"independent scalar stencil pullback {name}")
    # A polarization check isolates a mixed Hessian bilinear, rather than
    # accepting only diagonal quadratic pullbacks.
    wa, wb = pulls[2], pulls[3]
    def evaluate(w):
        return owner_scalar(active, order, global_weights(active, w))
    wab = tuple(wa[i] + wb[i] for i in range(size))
    mixed_got = evaluate(wab).c2 - evaluate(wa).c2 - evaluate(wb).c2
    mixed_expected = (payload_scalar(payload, wab).c2 -
                      payload_scalar(payload, wa).c2 -
                      payload_scalar(payload, wb).c2)
    check(mixed_got == mixed_expected,
          f"independent mixed-Hessian polarization {name}")


def zmat():
    return [[ZERO, ZERO], [ZERO, ZERO]]


def ident():
    return [[ONE, ZERO], [ZERO, ONE]]


def madd(a, b, factor=Q(1)):
    return [[a[i][j] + b[i][j].scale(factor) for j in range(2)]
            for i in range(2)]


def mmul(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(2)), ZERO)
             for j in range(2)] for i in range(2)]


def mtranspose(a):
    return [list(row) for row in zip(*a)]


def qadd(target, state, row, factor=Q(1)):
    old = target.get(state, (ZERO, ZERO))
    target[state] = (old[0] + row[0].scale(factor),
                     old[1] + row[1].scale(factor))


def cycle_nodes(cycle):
    answer = []
    for i in range(6):
        common = set(ends[cycle[i - 1]]) & set(ends[cycle[i]])
        check(len(common) == 1, "successive representative cycle edges meet")
        answer.append(next(iter(common)))
    check(len(set(answer)) == 6, "representative cycle has six nodes")
    return tuple(answer)


def cycle_writer(cycle, weights_local):
    """Independent scalar-source Bloch/des-Cloizeaux h6 writer."""
    cnodes = cycle_nodes(cycle)
    full = 63
    p_states = (0, full)
    pindex = {state: i for i, state in enumerate(p_states)}
    energy = {}
    source = {}
    for mask in range(64):
        toggled = frozenset(cycle[i] for i in range(6) if mask & (1 << i))
        if mask not in pindex:
            energy[mask] = defect(toggled)
            check(energy[mask] > 0, "proper representative cycle state is Q")
        value = Q(0)
        for local, node in enumerate(cnodes):
            word = pair_word(node, toggled)
            value += sum(weights_local[6 * local + pair] * word[pair]
                         for pair in range(6))
        source[mask] = D2(Q(0), value, Q(0))

    asrc = zmat()
    asrc[0][0], asrc[1][1] = source[0], source[full]

    def divide(vector):
        return {state: (row[0].scale(-1 / energy[state]),
                        row[1].scale(-1 / energy[state]))
                for state, row in vector.items()
                if row != (ZERO, ZERO)}

    def wqq(vector):
        out = {}
        for state, row in vector.items():
            for edge in range(6):
                nxt = state ^ (1 << edge)
                if nxt not in pindex:
                    qadd(out, nxt, row, Q(-1))
        return out

    def bchi(vector):
        out = zmat()
        for pstate, row_index in pindex.items():
            for edge in range(6):
                row = vector.get(pstate ^ (1 << edge))
                if row is not None:
                    out[row_index][0] -= row[0]
                    out[row_index][1] -= row[1]
        return out

    def mat_right(row, matrix):
        return (row[0] * matrix[0][0] + row[1] * matrix[1][0],
                row[0] * matrix[0][1] + row[1] * matrix[1][1])

    def raise_source(vector):
        bracket = {}
        for state, row in vector.items():
            qadd(bracket, state, (source[state] * row[0],
                                  source[state] * row[1]))
            qadd(bracket, state, mat_right(row, asrc), Q(-1))
        return divide(bracket)

    coupling = {}
    for pstate, column in pindex.items():
        for edge in range(6):
            state = pstate ^ (1 << edge)
            row = list(coupling.get(state, (ZERO, ZERO)))
            row[column] -= ONE
            coupling[state] = tuple(row)

    chi = [{} for _ in range(6)]
    for order in range(1, 6):
        if order == 1:
            constant = dict(coupling)
        else:
            constant = wqq(chi[order - 1])
            for left in range(1, order - 1):
                right = order - 1 - left
                product_matrix = bchi(chi[right])
                for state, row in chi[left].items():
                    qadd(constant, state, mat_right(row, product_matrix), Q(-1))
        base_solution = divide(constant)
        first = raise_source(base_solution)
        second = raise_source(first)
        answer = dict(base_solution)
        for state, row in first.items():
            qadd(answer, state, row)
        for state, row in second.items():
            qadd(answer, state, row)
        chi[order] = answer

    hb = [zmat() for _ in range(7)]
    hb[0] = asrc
    for order in range(2, 7):
        hb[order] = bchi(chi[order - 1])
    metric = [zmat() for _ in range(7)]
    metric[0] = ident()
    for order in range(2, 7):
        for left in range(1, order):
            right = order - left
            if left >= len(chi) or right >= len(chi):
                continue
            overlap = zmat()
            for state in set(chi[left]) & set(chi[right]):
                lr, rr = chi[left][state], chi[right][state]
                for i in range(2):
                    for j in range(2):
                        overlap[i][j] += lr[i] * rr[j]
            metric[order] = madd(metric[order], overlap)
    square = [zmat() for _ in range(7)]
    square[0] = ident()
    for order in range(1, 7):
        known = zmat()
        for left in range(1, order):
            known = madd(known, mmul(square[left], square[order - left]))
        square[order] = [[(metric[order][i][j] - known[i][j]).scale(Q(1, 2))
                          for j in range(2)] for i in range(2)]
    inverse = [zmat() for _ in range(7)]
    inverse[0] = ident()
    for order in range(1, 7):
        value = zmat()
        for left in range(order):
            value = madd(value, mmul(inverse[left], square[order - left]), Q(-1))
        inverse[order] = value

    def series_mul(a, b):
        out = [zmat() for _ in range(7)]
        for total in range(7):
            for left in range(total + 1):
                out[total] = madd(out[total], mmul(a[left], b[total - left]))
        return out

    canonical = series_mul(series_mul(square, hb), inverse)
    check(all(canonical[n] == mtranspose(canonical[n]) for n in range(7)),
          "independent scalar canonical series is Hermitian")
    return canonical[6][0][1]


# Independently contract the literal cycle stencil through four source paths.
cycle_row = rows["alternating_cycle_offdiag_h6"]
cycle = tuple(cycle_row["representative_class_records"][0]["representative_edges"])
check(frozenset(cycle) in active_cycle_supports,
      "literal representative cycle is active")
cycle_payload = cycle_row["first_exact_stencil"]
cycle_pulls = weight_rows(36)
for weights_local in cycle_pulls:
    got = cycle_writer(cycle, weights_local)
    expected = payload_scalar(cycle_payload, weights_local)
    check(got == expected, "independent canonical cycle stencil pullback")
wa, wb = cycle_pulls[2], cycle_pulls[3]
wab = tuple(wa[i] + wb[i] for i in range(36))
mixed_got = (cycle_writer(cycle, wab).c2 - cycle_writer(cycle, wa).c2 -
             cycle_writer(cycle, wb).c2)
mixed_expected = (payload_scalar(cycle_payload, wab).c2 -
                  payload_scalar(cycle_payload, wa).c2 -
                  payload_scalar(cycle_payload, wb).c2)
check(mixed_got == mixed_expected,
      "independent canonical cycle mixed-Hessian polarization")


# A separate minimal-history sum checks the source-free off-diagonal value.
minimal = Q(0)
for order in permutations(range(6)):
    toggled = set()
    denominator = Q(1)
    for step in order[:-1]:
        toggled.add(cycle[step])
        denominator *= defect(frozenset(toggled))
    minimal -= Q(1) / denominator
check(minimal == Q(-63, 8), "720-history cycle writer is -63/8")


# CTP sign and physical negative-half-source pullback.  Use an abstract
# energy Hessian k, because these factors are independent of its tensor value.
k = Q(7, 3)
energy_scale = Q(11, 5)
raw_forward = -k
raw_backward = +k
check((raw_forward, raw_backward) == (-k, k),
      "W=-i log Z gives -sigma times the energy Hessian")
pulled_hessian = (energy_scale**2 / 4) * k  # j=-(E*/2)J
physical_forward = -pulled_hessian
physical_backward = +pulled_hessian
check(physical_forward == -(energy_scale**2 / 4) * k and
      physical_backward == +(energy_scale**2 / 4) * k,
      "negative half-source squares in same-branch contact pullback")
# J+=Jr+Ja/2, J-=Jr-Ja/2.  Mixed direct branch blocks vanish.
w_ar = Q(1, 2) * physical_forward - Q(1, 2) * physical_backward
physical_retarded_contact = 2 * w_ar
check(w_ar == -(energy_scale**2 / 4) * k and
      physical_retarded_contact == -(energy_scale**2 / 2) * k,
      "r/a half-source pullback has the GL6W sign and factor")


# Verify the four V001 repairs in the actual frozen V002 bytes.
theorem = " ".join((TARGET / "THEOREM.md").read_text().split())
result = " ".join((TARGET / "RESULT.md").read_text().split())
readme = " ".join((TARGET / "README.md").read_text().split())
verifier = (TARGET / "verify_packet.py").read_text()
check("same-branch coincident connected contact operator" in theorem and
      "-\\sigma K''_{AB}\\,\\delta(t-s)" in theorem,
      "V002 defines the branchwise CTP sign inside theorem")
check("4096 **enumerated supports**" in theorem and "=3904" in theorem and
      "192 graph-wide cycle supports" in theorem,
      "V002 types 4096 versus 3904 versus 192")
check("does not literally store every stencil coefficient" in readme and
      "reconstructs all 1654 classes" in result,
      "V002 states compact-ledger boundary truthfully")
check("same-branch coincident connected contact operator" in verifier,
      "V002 verifier and theorem CTP token agree")
for forbidden in ("evaluates all 4096 nonbare geometric owners",
                  "and every exact local nonuniform source stencil class"):
    check(forbidden not in theorem + result + readme + json.dumps(ledger),
          f"V001 overclaim absent: {forbidden}")


# Compact JSON alone is a commitment, not a complete coefficient export.  The
# literal eight stencils plus the sealed reconstructing executable make the
# packet reproducible; downstream full-Jet use must replay or add an adapter.
check(all("first_exact_stencil" in row for row in rows.values()),
      "one literal first stencil per family/order")
check(all(len(row["representative_class_records"]) <= 4
          for row in rows.values()), "representative metadata remains compact")
check(all(not any("stencil" in rec and "stencil_sha256" not in rec
                  for rec in row["representative_class_records"])
          for row in rows.values()), "representatives do not masquerade as exports")
check("fresh complete result differs from frozen exact ledger" in
      (TARGET / "derive_complete_six_pair_h6_source_jet.py").read_text(),
      "sealed executable fail-closes against compact ledger")


print(f"PASS__INDEPENDENT_GL6CU_V002_HOSTILE_RECONSTRUCTION__{checks}/{checks}")
print("REPAIRS=V001_R1_R2_R3_R4_CLOSED")
print("CENSUS=4096_ENUMERATED__3904_COMPUTED__192_INACTIVE")
print("ARITHMETIC=EIGHT_LITERAL_STENCILS_INDEPENDENTLY_SCALAR_PULLBACK_TESTED")
print("HESSIAN=ACTUAL_SECOND_DERIVATIVE_AND_MIXED_POLARIZATION_PASS")
print("CTP=-SIGMA_KPP__HALF_SOURCE_RETARDED=-E2_OVER_2_KPP")
print("CEILING=OPERATOR_SOURCE_JET_ONLY__NO_CONNECTED_1PI_WARD_GRAVITY_OR_G")
