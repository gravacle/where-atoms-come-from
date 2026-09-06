#!/usr/bin/env python3
"""GL6CU V002: exact complete raw six-pair source 2-jet through h^6.

This executable never imports a frozen predecessor.  It constructs the
selected-branch and incident-cycle source-before-projection canonical
effective-Hamiltonian inventory from exact rational arithmetic.  Diagonal terms are assigned to their distinct-edge support by
Möbius inversion; the only configuration-changing owner through sixth order
is an alternating six-cycle, treated with the full two-state Bloch/des-
Cloizeaux canonicalization.

The source coordinates are the six raw pair observables at every constraint
node.  A Jet stores the value, first derivative, and actual mixed second
derivative (not a Taylor coefficient divided by two) at source zero.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PERIOD = 4
STEPS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
PAIR_ORDER = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIR_ORDER)}
SELECTOR = frozenset(((0, 1), (0, 2), (1, 1),
                      (1, 3), (2, 0), (2, 1)))
RAW_DIRECTIONS = (
    (1, 1, 1, 1, 1, 1),
    (1, 1, -2, -2, 1, 1),
    (1, -1, 0, 0, -1, 1),
    (1, 0, 0, 0, 0, -1),
    (0, 1, 0, 0, -1, 0),
    (0, 0, 1, -1, 0, 0),
)
SOURCE_NAMES = ("A", "Ea", "Eb", "T1", "T2", "T3")
T_DIRECTIONS = RAW_DIRECTIONS[3:]
CHECKS: list[str] = []


def check(condition, label):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def qtext(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest_payload(value):
    return sha256(canonical_json(value).encode()).hexdigest()


def file_digest(path):
    return sha256(path.read_bytes()).hexdigest()


def clean(mapping):
    return {key: F(value) for key, value in mapping.items() if value}


def hp(left, right):
    return (left, right) if left <= right else (right, left)


@dataclass(frozen=True)
class Jet:
    """Sparse exact source 2-jet with actual mixed Hessian entries."""

    value: F = F(0)
    grad: dict[int, F] | None = None
    hess: dict[tuple[int, int], F] | None = None

    def __post_init__(self):
        object.__setattr__(self, "value", F(self.value))
        object.__setattr__(self, "grad", clean(self.grad or {}))
        object.__setattr__(self, "hess", clean(self.hess or {}))

    @staticmethod
    def linear(gradient):
        return Jet(F(0), gradient, {})

    def __add__(self, other):
        other = as_jet(other)
        grad = dict(self.grad)
        for key, value in other.grad.items():
            grad[key] = grad.get(key, F(0)) + value
        hess = dict(self.hess)
        for key, value in other.hess.items():
            hess[key] = hess.get(key, F(0)) + value
        return Jet(self.value + other.value, grad, hess)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value,
                   {key: -value for key, value in self.grad.items()},
                   {key: -value for key, value in self.hess.items()})

    def __sub__(self, other):
        return self + (-as_jet(other))

    def __rsub__(self, other):
        return as_jet(other) - self

    def __mul__(self, other):
        other = as_jet(other)
        grad = {}
        for key in set(self.grad) | set(other.grad):
            grad[key] = self.grad.get(key, F(0)) * other.value + self.value * other.grad.get(key, F(0))
        hess = {}
        for key in set(self.hess) | set(other.hess):
            hess[key] = (self.hess.get(key, F(0)) * other.value
                         + self.value * other.hess.get(key, F(0)))
        for left, lvalue in self.grad.items():
            for right, rvalue in other.grad.items():
                key = hp(left, right)
                hess[key] = hess.get(key, F(0)) + (2 if left == right else 1) * lvalue * rvalue
        return Jet(self.value * other.value, grad, hess)

    __rmul__ = __mul__

    def reciprocal(self):
        if not self.value:
            raise ZeroDivisionError("jet reciprocal at zero value")
        grad = {key: -value / self.value ** 2
                for key, value in self.grad.items()}
        keys = set(self.hess)
        for left in self.grad:
            for right in self.grad:
                keys.add(hp(left, right))
        hess = {}
        for left, right in keys:
            hess[(left, right)] = (
                2 * self.grad.get(left, F(0)) * self.grad.get(right, F(0)) / self.value ** 3
                - self.hess.get((left, right), F(0)) / self.value ** 2
            )
        return Jet(1 / self.value, grad, hess)

    def __truediv__(self, other):
        return self * as_jet(other).reciprocal()


def as_jet(value):
    return value if isinstance(value, Jet) else Jet(F(value))


class JetAccumulator:
    """Mutable exact sum used only at the global owner-aggregation boundary."""

    def __init__(self):
        self.value = F(0)
        self.grad = {}
        self.hess = {}

    def add(self, jet):
        self.value += jet.value
        for key, value in jet.grad.items():
            updated = self.grad.get(key, F(0)) + value
            if updated:
                self.grad[key] = updated
            else:
                self.grad.pop(key, None)
        for key, value in jet.hess.items():
            updated = self.hess.get(key, F(0)) + value
            if updated:
                self.hess[key] = updated
            else:
                self.hess.pop(key, None)

    def freeze(self):
        return Jet(self.value, self.grad, self.hess)


ZERO = Jet()
ONE = Jet(F(1))


# ---------------------------------------------------------------- Q4 parent
def add(left, right):
    return tuple((left[i] + right[i]) % PERIOD for i in range(3))


def sub(left, right):
    return tuple((left[i] - right[i]) % PERIOD for i in range(3))


CELLS = tuple(product(range(PERIOD), repeat=3))
EDGES = tuple((cell, port) for cell in CELLS for port in range(4))
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
NODES = tuple((kind, cell) for kind in ("P", "C") for cell in CELLS)
NODE_INDEX = {node: index for index, node in enumerate(NODES)}


def endpoints(edge):
    cell, port = edge
    return ("P", cell), ("C", add(cell, STEPS[port]))


ENDS = tuple(endpoints(edge) for edge in EDGES)
# The tuple position is the physical port label.  In particular, quotient
# edge-enumeration order is not port order at a C node.
INCIDENT = {
    (kind, cell): tuple(
        EDGE_INDEX[(cell, port)] if kind == "P" else
        EDGE_INDEX[(sub(cell, STEPS[port]), port)]
        for port in range(4)
    )
    for kind, cell in NODES
}


def selector(a, b):
    return int((a % PERIOD, b % PERIOD) in SELECTOR)


def background_occupation(edge):
    cell, port = edge
    x, y, z = cell
    a, b = (y - x) % PERIOD, (x + z) % PERIOD
    word = (selector(a - 1, b + 1), selector(a, b),
            1 - selector(a - 1, b + 1), 1 - selector(a, b))
    return word[port]


BASE = frozenset(index for index, edge in enumerate(EDGES)
                 if background_occupation(edge))


def bit(occupied, edge_index, toggled=frozenset()):
    return int(edge_index in occupied) ^ int(edge_index in toggled)


def pair_word(occupied, node, toggled=frozenset()):
    z = tuple(1 - 2 * bit(occupied, edge_index, toggled)
              for edge_index in INCIDENT[node])
    return tuple(F(z[a] * z[b]) for a, b in PAIR_ORDER)


def locked(occupied):
    return all(sum(int(edge in occupied) for edge in INCIDENT[node]) == 2
               for node in NODES)


def defect_energy(occupied, toggled):
    charge = {}
    for edge_index in toggled:
        sign = 1 if edge_index not in occupied else -1
        for node in ENDS[edge_index]:
            charge[node] = charge.get(node, 0) + sign
    return F(sum(value * value for value in charge.values()))


def source_delta(occupied, toggled):
    affected = {node for edge in toggled for node in ENDS[edge]}
    gradient = {}
    for node in affected:
        before = pair_word(occupied, node)
        after = pair_word(occupied, node, frozenset(toggled))
        base = 6 * NODE_INDEX[node]
        for pair in range(6):
            value = after[pair] - before[pair]
            if value:
                gradient[base + pair] = value
    return gradient


check(len(EDGES) == 256 and len(NODES) == 128,
      "Q4 has 256 links and 128 nodes")
check(len(BASE) == 128 and locked(BASE),
      "selected GL6CC Q4 state is degree-two locked")
check(all(len(row) == 4 for row in INCIDENT.values()),
      "Q4 parent is degree four")
check(all(EDGES[INCIDENT[("P", cell)][port]] == (cell, port)
          for cell in CELLS for port in range(4)),
      "every P-node incidence tuple position is its physical port label")
check(all(EDGES[INCIDENT[("C", cell)][port]] ==
          (sub(cell, STEPS[port]), port)
          for cell in CELLS for port in range(4)),
      "every C-node incidence tuple position is its physical port label")


# ------------------------------------------- scalar linked-owner recurrence
def remap_jet(jet, coordinate_map):
    return Jet(
        jet.value,
        {coordinate_map[key]: value for key, value in jet.grad.items()},
        {hp(coordinate_map[left], coordinate_map[right]): value
         for (left, right), value in jet.hess.items()},
    )


def scalar_cluster_problem(active, occupied):
    """Return a translation-free local problem and its global coordinate map."""
    active = tuple(sorted(active))
    touched = tuple(sorted({node for edge in active for node in ENDS[edge]}, key=repr))
    local_node = {node: index for index, node in enumerate(touched)}
    coordinate_map = {
        6 * local_node[node] + pair: 6 * NODE_INDEX[node] + pair
        for node in touched for pair in range(6)
    }
    rows = []
    for mask in range(1, 1 << len(active)):
        toggled = frozenset(active[i] for i in range(len(active)) if mask & (1 << i))
        energy = defect_energy(occupied, toggled)
        if energy <= 0:
            raise AssertionError(f"early locked return in scalar owner {active}/{mask}")
        global_delta = source_delta(occupied, toggled)
        local_delta = tuple(sorted(
            (6 * local_node[NODES[key // 6]] + key % 6, value)
            for key, value in global_delta.items()
        ))
        rows.append((energy, local_delta))
    return (len(active), tuple(rows)), coordinate_map


LOCAL_SCALAR_CACHE = {}


def scalar_cluster_coefficients(active, occupied=BASE):
    """Intermediate-normalized branch coefficients e_0..e_6 as Jets.

    Active is a tuple of distinct physical edge indices.  For |active|<=3,
    girth>=6 ensures that mask zero is the only locked configuration.  Exact
    translated copies share one local source-jet calculation and are remapped
    afterward; no source coordinate is averaged or identified.
    """
    active = tuple(sorted(active))
    signature, coordinate_map = scalar_cluster_problem(active, occupied)
    if signature in LOCAL_SCALAR_CACHE:
        return tuple(remap_jet(jet, coordinate_map)
                     for jet in LOCAL_SCALAR_CACHE[signature])
    size, rows = signature
    states = range(1 << size)
    deltas = {}
    for mask, (energy, gradient) in enumerate(rows, start=1):
        deltas[mask] = Jet(energy, dict(gradient), {})

    coefficient = [ZERO for _ in range(7)]
    wave = [dict() for _ in range(7)]
    wave[0] = {0: ONE}
    for order in range(1, 7):
        current = {}
        for mask in states:
            if not mask:
                continue
            # W=-sum X.  The numerator is -W c_(n-1)+sum e_m c_(n-m).
            numerator = ZERO
            for local_edge in range(size):
                numerator += wave[order - 1].get(mask ^ (1 << local_edge), ZERO)
            for lower in range(1, order):
                numerator += coefficient[lower] * wave[order - lower].get(mask, ZERO)
            current[mask] = numerator / deltas[mask]
        wave[order] = current
        coefficient[order] = -sum(
            (wave[order - 1].get(1 << local_edge, ZERO)
             for local_edge in range(size)), ZERO)
    LOCAL_SCALAR_CACHE[signature] = tuple(coefficient)
    return tuple(remap_jet(jet, coordinate_map) for jet in coefficient)


class LinkedOwners:
    def __init__(self, occupied):
        self.occupied = occupied
        self.raw = {}
        self.owner = {}

    def raw_coefficients(self, active):
        key = frozenset(active)
        if key not in self.raw:
            self.raw[key] = scalar_cluster_coefficients(tuple(key), self.occupied)
        return self.raw[key]

    def coefficient(self, active, order):
        key = frozenset(active)
        cache_key = (key, order)
        if cache_key in self.owner:
            return self.owner[cache_key]
        answer = self.raw_coefficients(key)[order]
        ordered = tuple(sorted(key))
        for size in range(1, len(ordered)):
            for subset in combinations(ordered, size):
                answer -= self.coefficient(subset, order)
        self.owner[cache_key] = answer
        return answer

    def release(self, active, order):
        """Release a completed high-cardinality owner after its ledger use."""
        key = frozenset(active)
        self.owner.pop((key, order), None)
        self.raw.pop(key, None)


# -------------------------------------------------------- exact owner census
WEDGES = set()
STARS = set()
for node in NODES:
    for pair in combinations(INCIDENT[node], 2):
        WEDGES.add(frozenset(pair))
    for triple in combinations(INCIDENT[node], 3):
        STARS.add(frozenset(triple))

PATHS = set()
for middle, (left_node, right_node) in enumerate(ENDS):
    for left in INCIDENT[left_node]:
        if left == middle:
            continue
        for right in INCIDENT[right_node]:
            if right == middle:
                continue
            PATHS.add(frozenset((left, middle, right)))


def canonical_cycle(cell, ports):
    a, b, c = ports
    ab = add(sub(cell, STEPS[b]), STEPS[a])
    cb = add(sub(cell, STEPS[b]), STEPS[c])
    return ((cell, a), (ab, b), (ab, c),
            (cb, a), (cb, b), (cell, c))


CYCLES = []
cycle_seen = set()
for cell in CELLS:
    for ports in combinations(range(4), 3):
        row = tuple(EDGE_INDEX[edge] for edge in canonical_cycle(cell, ports))
        key = frozenset(row)
        check(len(key) == 6 and key not in cycle_seen,
              "canonical Q4 six-cycle owner is unique")
        cycle_seen.add(key)
        CYCLES.append(row)

check(len(WEDGES) == 768, "Q4 has 768 unique node-owned wedges")
check(len(STARS) == 512, "Q4 has 512 unique three-edge stars")
check(len(PATHS) == 2304, "Q4 has 2304 unique length-three paths")
check(len(CYCLES) == 256, "Q4 has 256 unique elementary six-cycles")
check(not (STARS & PATHS), "star and path owner classes are disjoint")


def is_alternating(occupied, cycle):
    word = tuple(int(edge in occupied) for edge in cycle)
    return word in ((1, 0, 1, 0, 1, 0), (0, 1, 0, 1, 0, 1))


ACTIVE_CYCLES = tuple(cycle for cycle in CYCLES if is_alternating(BASE, cycle))
check(len(ACTIVE_CYCLES) == 64,
      "selected Q4 state has 64 active alternating cycle transitions")


# ---------------------------------- two-state canonical six-cycle Jet engine
def jzero_matrix(size=2):
    return [[ZERO for _ in range(size)] for _ in range(size)]


def jidentity(size=2):
    answer = jzero_matrix(size)
    for i in range(size):
        answer[i][i] = ONE
    return answer


def jmadd(left, right, factor=F(1)):
    return [[left[i][j] + factor * right[i][j]
             for j in range(len(left[0]))] for i in range(len(left))]


def jmmul(left, right):
    return [[sum((left[i][k] * right[k][j]
                  for k in range(len(right))), ZERO)
             for j in range(len(right[0]))] for i in range(len(left))]


def jmtranspose(matrix):
    return [list(row) for row in zip(*matrix)]


def series_product(left, right):
    output = [jzero_matrix() for _ in range(7)]
    for total in range(7):
        for first in range(total + 1):
            output[total] = jmadd(
                output[total], jmmul(left[first], right[total - first]))
    return output


def qvadd(target, source, matrix=None, factor=F(1)):
    for state, row in source.items():
        if matrix is None:
            addition = tuple(factor * value for value in row)
        else:
            addition = tuple(
                factor * sum((row[k] * matrix[k][column]
                              for k in range(2)), ZERO)
                for column in range(2))
        prior = target.get(state, (ZERO, ZERO))
        value = (prior[0] + addition[0], prior[1] + addition[1])
        if value == (ZERO, ZERO):
            target.pop(state, None)
        else:
            target[state] = value


LOCAL_CYCLE_CACHE = {}


def remap_matrix_series(series, coordinate_map):
    return [
        [[remap_jet(entry, coordinate_map) for entry in row]
         for row in matrix]
        for matrix in series
    ]


def cycle_canonical_jet(cycle, occupied=BASE):
    """Canonical 2x2 effective Hamiltonian Jets for one alternating cycle."""
    cycle = tuple(cycle)
    if not is_alternating(occupied, cycle):
        raise AssertionError("cycle owner requires an alternating locked endpoint")
    full = (1 << 6) - 1
    p_states = (0, full)
    p_index = {state: index for index, state in enumerate(p_states)}
    q_states = tuple(state for state in range(1 << 6) if state not in p_index)
    q_energy = {}
    source = {}
    cycle_nodes = []
    for index in range(6):
        common = set(ENDS[cycle[index - 1]]) & set(ENDS[cycle[index]])
        check(len(common) == 1,
              "successive canonical cycle edges meet at one constraint node")
        cycle_nodes.append(next(iter(common)))
    check(len(set(cycle_nodes)) == 6,
          "alternating owner has six distinct constraint nodes")
    local_node = {node: index for index, node in enumerate(cycle_nodes)}
    local_to_global = {
        6 * local_node[node] + pair: 6 * NODE_INDEX[node] + pair
        for node in cycle_nodes for pair in range(6)
    }
    for mask in range(1 << 6):
        toggled = frozenset(cycle[i] for i in range(6) if mask & (1 << i))
        if mask in p_index:
            check(defect_energy(occupied, toggled) == 0,
                  "cycle endpoints are locked")
        else:
            q_energy[mask] = defect_energy(occupied, toggled)
            check(q_energy[mask] > 0, "proper cycle subset remains in Q")
        gradient = {}
        for node in cycle_nodes:
            word = pair_word(occupied, node, toggled)
            base = 6 * local_node[node]
            for pair, value in enumerate(word):
                if value:
                    gradient[base + pair] = value
        source[mask] = Jet.linear(gradient)

    # This is a proof-complete exact local problem signature: it retains the
    # defect energy and all six raw source eigenvalues for every one of the 64
    # cycle subsets.  Equal signatures therefore have identical canonical
    # Jets, while their independent global source coordinates are remapped.
    signature = tuple(
        (q_energy.get(mask, F(0)), tuple(sorted(source[mask].grad.items())))
        for mask in range(1 << 6)
    )
    if signature in LOCAL_CYCLE_CACHE:
        return remap_matrix_series(LOCAL_CYCLE_CACHE[signature], local_to_global), tuple(cycle_nodes)

    a_source = jzero_matrix()
    for index, state in enumerate(p_states):
        a_source[index][index] = source[state]

    def divide(vector):
        return {state: tuple(-value / q_energy[state] for value in row)
                for state, row in vector.items() if row != (ZERO, ZERO)}

    def apply_wqq(vector):
        output = {}
        for state, row in vector.items():
            for local_edge in range(6):
                target = state ^ (1 << local_edge)
                if target in p_index:
                    continue
                qvadd(output, {target: row}, factor=F(-1))
        return output

    def b_times(vector):
        output = jzero_matrix()
        for pstate, row_index in p_index.items():
            for local_edge in range(6):
                state = pstate ^ (1 << local_edge)
                values = vector.get(state)
                if values is None:
                    continue
                for column in range(2):
                    output[row_index][column] -= values[column]
        return output

    def source_raise(vector):
        bracket = {}
        for state, row in vector.items():
            qvadd(bracket, {state: (source[state] * row[0],
                                    source[state] * row[1])})
        qvadd(bracket, vector, a_source, factor=F(-1))
        return divide(bracket)

    coupling = {}
    for pstate, column in p_index.items():
        for local_edge in range(6):
            state = pstate ^ (1 << local_edge)
            row = list(coupling.get(state, (ZERO, ZERO)))
            row[column] -= ONE
            coupling[state] = tuple(row)

    chi = [{} for _ in range(6)]
    for order in range(1, 6):
        if order == 1:
            constant = coupling
        else:
            constant = apply_wqq(chi[order - 1])
            for left in range(1, order - 1):
                right = order - 1 - left
                qvadd(constant, chi[left], b_times(chi[right]), factor=F(-1))
        base_solution = divide(constant)
        first_raise = source_raise(base_solution)
        second_raise = source_raise(first_raise)
        answer = dict(base_solution)
        qvadd(answer, first_raise)
        qvadd(answer, second_raise)
        chi[order] = answer

    hb = [jzero_matrix() for _ in range(7)]
    hb[0] = a_source
    for order in range(2, 7):
        hb[order] = b_times(chi[order - 1])

    metric = [jzero_matrix() for _ in range(7)]
    metric[0] = jidentity()
    for order in range(2, 7):
        for left in range(1, order):
            right = order - left
            if left >= len(chi) or right >= len(chi):
                continue
            overlap = jzero_matrix()
            common = set(chi[left]) & set(chi[right])
            for state in common:
                lrow, rrow = chi[left][state], chi[right][state]
                for i in range(2):
                    for j in range(2):
                        overlap[i][j] += lrow[i] * rrow[j]
            metric[order] = jmadd(metric[order], overlap)

    square = [jzero_matrix() for _ in range(7)]
    square[0] = jidentity()
    for order in range(1, 7):
        known = jzero_matrix()
        for left in range(1, order):
            known = jmadd(known, jmmul(square[left], square[order - left]))
        square[order] = [[F(1, 2) * (metric[order][i][j] - known[i][j])
                          for j in range(2)] for i in range(2)]

    inverse = [jzero_matrix() for _ in range(7)]
    inverse[0] = jidentity()
    for order in range(1, 7):
        value = jzero_matrix()
        for left in range(order):
            value = jmadd(value, jmmul(inverse[left], square[order - left]),
                          factor=F(-1))
        inverse[order] = value

    canonical = series_product(series_product(square, hb), inverse)
    for order in range(7):
        check(canonical[order] == jmtranspose(canonical[order]),
              f"cycle canonical h{order} Jet matrix is Hermitian")
    LOCAL_CYCLE_CACHE[signature] = canonical
    return remap_matrix_series(canonical, local_to_global), tuple(cycle_nodes)


# ------------------------------------------------------- exact contractions
def jet_payload(jet):
    return {
        "value": qtext(jet.value),
        "gradient": [[key, qtext(value)] for key, value in sorted(jet.grad.items())],
        "hessian": [[left, right, qtext(value)]
                    for (left, right), value in sorted(jet.hess.items())],
    }


def jet_hash(jet):
    return digest_payload(jet_payload(jet))


def localize_owner_jet(jet, active, ordered_nodes=None):
    """Express an owner Jet in canonical local-node/raw-pair coordinates."""
    incident = {node for edge in active for node in ENDS[edge]}
    touched = (tuple(ordered_nodes) if ordered_nodes is not None else
               tuple(sorted(incident, key=repr)))
    check(len(touched) == len(incident) and set(touched) == incident,
          "local stencil node order covers every incident node exactly once")
    coordinate_map = {
        6 * NODE_INDEX[node] + pair: 6 * local + pair
        for local, node in enumerate(touched) for pair in range(6)
    }
    check(set(jet.grad).issubset(coordinate_map),
          "owner gradient is supported only on incident constraint nodes")
    check(all(left in coordinate_map and right in coordinate_map
              for left, right in jet.hess),
          "owner Hessian is supported only on incident constraint nodes")
    return remap_jet(jet, coordinate_map), touched


def record_stencil(groups, family, order, active, jet, ordered_nodes=None):
    """Classify a raw nonuniform owner stencil without identifying sources."""
    local, touched = localize_owner_jet(jet, active, ordered_nodes)
    payload = jet_payload(local)
    key = canonical_json(payload)
    bucket = groups.setdefault((family, order), {})
    if key not in bucket:
        bucket[key] = {
            "count": 0,
            "representative_edges": list(active),
            "incident_node_count": len(touched),
            "stencil": payload,
        }
    bucket[key]["count"] += 1


def stencil_ledger(groups):
    answer = {}
    for (family, order), buckets in sorted(groups.items()):
        rows = sorted(buckets.values(), key=lambda row: digest_payload(row["stencil"]))
        class_rows = []
        for row in rows:
            stencil = row["stencil"]
            class_rows.append({
                "count": row["count"],
                "representative_edges": row["representative_edges"],
                "incident_node_count": row["incident_node_count"],
                "value": stencil["value"],
                "gradient_entries": len(stencil["gradient"]),
                "hessian_entries": len(stencil["hessian"]),
                "stencil_sha256": digest_payload(stencil),
            })
        multiplicities = {}
        for row in class_rows:
            key = str(row["count"])
            multiplicities[key] = multiplicities.get(key, 0) + 1
        answer[f"{family}_h{order}"] = {
            "class_count": len(rows),
            "owner_count": sum(row["count"] for row in rows),
            "class_multiplicity_histogram": multiplicities,
            "gradient_entry_count_classes": sorted({
                row["gradient_entries"] for row in class_rows}),
            "hessian_entry_count_classes": sorted({
                row["hessian_entries"] for row in class_rows}),
            "all_class_records_sha256": digest_payload(class_rows),
            "all_exact_stencils_sha256": digest_payload(
                [row["stencil"] for row in rows]),
            "representative_class_records": (
                class_rows if len(class_rows) <= 4 else
                class_rows[:2] + class_rows[-2:]
            ),
            "first_exact_stencil": rows[0]["stencil"],
        }
    return answer


def matrix_from_hessian(jet, left_node, right_node):
    left_base, right_base = 6 * NODE_INDEX[left_node], 6 * NODE_INDEX[right_node]
    matrix = [[F(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(6):
            matrix[i][j] = jet.hess.get(hp(left_base + i, right_base + j), F(0))
    return matrix


def matvec(matrix, vector):
    return [sum((matrix[i][j] * F(vector[j]) for j in range(6)), F(0))
            for i in range(6)]


def field_bilinear(jet, left_field, right_field):
    answer = F(0)
    for (left, right), value in jet.hess.items():
        if left == right:
            answer += F(left_field.get(left, 0)) * value * F(right_field.get(right, 0))
        else:
            answer += value * (
                F(left_field.get(left, 0)) * F(right_field.get(right, 0))
                + F(left_field.get(right, 0)) * F(right_field.get(left, 0)))
    return answer


def uniform_field(direction):
    return {6 * node + pair: F(direction[pair])
            for node in range(len(NODES)) for pair in range(6)
            if direction[pair]}


def project_t_gradient(jet, node):
    base = 6 * NODE_INDEX[node]
    return tuple(sum(F(direction[pair]) * jet.grad.get(base + pair, F(0))
                     for pair in range(6)) for direction in T_DIRECTIONS)


def frozen_regressions(diagonal):
    dependencies = {
        "CF": (ROOT / "LANE_CROSS_RFT_GRA_GL6CF_DENSE_PARENT_GLOBAL_H2_CONTACT_SYMBOL_V001/EXACT_LEDGER.json",
               "bbe21a2739c8b6848eb0069b6c3e7b14577af87c0db524c3da6db87275dd92a1"),
        "CG": (ROOT / "LANE_CROSS_RFT_GRA_GL6CG_DENSE_PARENT_GLOBAL_H4_SOURCE_OPERATOR_V001/EXACT_LEDGER.json",
               "ab38bdac8fac569ae890e20a79d24a4e7e75aa639fb81621462421de2ac33c21"),
        "CG_STENCIL": (ROOT / "LANE_CROSS_RFT_GRA_GL6CG_DENSE_PARENT_GLOBAL_H4_SOURCE_OPERATOR_V001/OPERATOR_STENCIL.json",
                       "2fb7835f6880631cbb0ddd29a23b8f72e2b192a6a939aa42bd054e1a45f8a561"),
        "CH": (ROOT / "LANE_CROSS_RFT_GRA_GL6CH_GLOBAL_H6_TENSOR_WRITER_V001/EXACT_LEDGER.json",
               "b934da4d79fee1ccdab92886c2fdbc2749a2470b156ec66ee63cb9bcdfe9a7a0"),
        "CN": (ROOT / "LANE_CROSS_RFT_GRA_GL6CN_COMPLETE_DIAGONAL_H6_T2_FIRST_SOURCE_V001/EXACT_LEDGER.json",
               "a950362c4c23abdad52786b858395c8140af9f6b11511c20538c34ae14a5a3da"),
        "BV": (ROOT / "LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001/VERIFICATION.txt",
               "f612288cf1c4fec3c366765f3112cdfe19f7ae03e3be847f21e74ac9e42d409a"),
        "BV_THEOREM": (ROOT / "LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001/THEOREM.md",
                       "a8daf338600135af58cc65f76612c91d3efdf26117af5fbd9ae94d6e018d63db"),
        "BV_MANIFEST": (ROOT / "LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001/MANIFEST.sha256",
                        "5e6ba6ba609a5627bf0d9fdac91464be1e7727b9b5ef718dd8e891266f0b7abf"),
        "BV_SEAL": (ROOT / "LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001/SEAL.sha256",
                    "a1c6f8f537474219be0f3dda7696da0ed290b6d001f367bde6aab04100c8472e"),
        "AO": (ROOT / "LANE_CROSS_RFT_GRA_GL6AO_COMPLETE_SIXTH_ORDER_LOCKED_HAMILTONIAN_V001/VERIFICATION.txt",
               "fdab6654cf9bb1ead27f9117dfd5ed1692a022145f24a1061d0e694a7df66fd1"),
        "CU_V001_THEOREM": (ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/THEOREM.md",
                            "22a23eb8c1d4ad2c1f49eb7bb37046ba73df39bdfb34c2b77a90a3258d8647a1"),
        "CU_V001_MANIFEST": (ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/MANIFEST.sha256",
                             "a0b79af8608ae8f151dbfcc62c63ec581643cc4c411c789c71859fee7b86a5c3"),
        "CU_V001_SEAL": (ROOT / "LANE_CROSS_RFT_GRA_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/SEAL.sha256",
                         "359b55882ccc4d015be1592d2f31dc9e1556bc6627323df31480935b189674d5"),
        "CU_V001_AUDIT_REPORT": (ROOT / "AUDIT_G_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/AUDIT_REPORT.md",
                                 "479d64500cc334c47f19ef0336f1e11de8e7e9004810583247a4681006de4e6c"),
        "CU_V001_AUDIT_MANIFEST": (ROOT / "AUDIT_G_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/MANIFEST.sha256",
                                   "c3a0ae90e13fdd1900c493b1409fb8eaee82d50b4348e80c34016c6f93692e2c"),
        "CU_V001_AUDIT_SEAL": (ROOT / "AUDIT_G_GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V001/SEAL.sha256",
                               "57b2127d8f96935b8beaa964f8e095751a2f5106a245a2f7637c527cdb292797"),
    }
    for label, (path, expected) in dependencies.items():
        check(path.is_file() and file_digest(path) == expected,
              f"frozen {label} regression input has pinned bytes")

    cf = json.loads(dependencies["CF"][0].read_text())
    ch = json.loads(dependencies["CH"][0].read_text())
    check(ch["direct_hexagon"]["amplitude"] == "-63/8" and
          ch["direct_hexagon"]["canonical_pair_gradient_coefficient"] == "105/8",
          "frozen CH regression constants are typed as expected")

    # CF common uniform source block, per 64 parent cells.
    fields = [uniform_field(direction) for direction in RAW_DIRECTIONS]
    common = [[field_bilinear(diagonal[2], fields[i], fields[j]) / 64
               for j in range(6)] for i in range(6)]
    expected = [[F(value) for value in row] for row in cf["symbol"]["k0"]]
    check(common == expected,
          "complete arbitrary-source h2 Hessian contracts to frozen CF k0 symbol")

    # BV onsite tight-frame block in every raw T direction.
    for node in NODES:
        block = matrix_from_hessian(diagonal[2], node, node)
        for direction in T_DIRECTIONS:
            check(matvec(block, direction) == [-2 * F(value) for value in direction],
                  "h2 onsite tensor contact equals -2 P_T")

    # CG and CN first-source regressions, pointwise at every node.
    for node in NODES:
        memory = pair_word(BASE, node)
        base = 6 * NODE_INDEX[node]
        got2 = tuple(diagonal[2].grad.get(base + pair, F(0)) for pair in range(6))
        got4 = tuple(diagonal[4].grad.get(base + pair, F(0)) for pair in range(6))
        check(got2 == tuple(-value for value in memory),
              "complete h2 first source equals -M pointwise")
        check(got4 == tuple(-F(4, 9) - F(37, 12) * value for value in memory),
              "complete h4 first source matches frozen CG pointwise operator")
        check(project_t_gradient(diagonal[6], node) == (F(0), F(0), F(0)),
              "complete diagonal h6 first source is T2-dark as frozen CN requires")
    return {label: expected for label, (_path, expected) in dependencies.items()}


def run():
    started = time.monotonic()
    engine = LinkedOwners(BASE)
    diagonal_sums = {2: JetAccumulator(), 4: JetAccumulator(), 6: JetAccumulator()}
    stencil_groups = {}
    bare_gradient = {}
    bare_word_census = {}
    for node in NODES:
        word = pair_word(BASE, node)
        bare_word_census[word] = bare_word_census.get(word, 0) + 1
        base = 6 * NODE_INDEX[node]
        for pair, value in enumerate(word):
            bare_gradient[base + pair] = value
    bare = Jet.linear(bare_gradient)
    check(len(bare.grad) == 6 * len(NODES) and not bare.hess,
          "bare selected-branch P source retains all 768 raw coordinates and has zero Hessian")

    for edge in range(len(EDGES)):
        owner = (edge,)
        for order in (2, 4, 6):
            contribution = engine.coefficient(owner, order)
            diagonal_sums[order].add(contribution)
            record_stencil(stencil_groups, "one_edge", order, owner, contribution)
    for owner_set in sorted(WEDGES, key=lambda row: tuple(sorted(row))):
        owner = tuple(owner_set)
        for order in (4, 6):
            contribution = engine.coefficient(owner, order)
            diagonal_sums[order].add(contribution)
            record_stencil(stencil_groups, "wedge", order, owner, contribution)
    for owner_set in sorted(STARS | PATHS, key=lambda row: tuple(sorted(row))):
        owner = tuple(owner_set)
        family = "star" if owner_set in STARS else "path"
        contribution = engine.coefficient(owner, 6)
        diagonal_sums[6].add(contribution)
        record_stencil(stencil_groups, family, 6, owner, contribution)
        engine.release(owner, 6)

    diagonal = {order: accumulator.freeze()
                for order, accumulator in diagonal_sums.items()}

    check(diagonal[2].value == -F(256, 2),
          "source-off h2 diagonal sum matches AO")
    check(diagonal[4].value == -F(7 * 256, 24),
          "source-off h4 diagonal sum matches AO")
    check(diagonal[6].value == -F(893 * 256, 1080),
          "source-off h6 diagonal sum matches AO")

    # Disconnected cumulants are zero as full Jets, not only source-free.
    disjoint_pair = next(
        frozenset((left, right))
        for left, right in combinations(range(len(EDGES)), 2)
        if frozenset((left, right)) not in WEDGES)
    check(engine.coefficient(tuple(disjoint_pair), 6) == ZERO,
          "representative disjoint two-edge owner vanishes as a complete source Jet")
    adjacent = next(iter(WEDGES))
    isolated = next(edge for edge in range(len(EDGES))
                    if all(frozenset((edge, member)) not in WEDGES
                           for member in adjacent))
    one_adjacent = tuple(adjacent | {isolated})
    check(engine.coefficient(one_adjacent, 6) == ZERO,
          "representative wedge-plus-isolated h6 owner vanishes as a complete source Jet")
    matching = None
    for trial in combinations(range(len(EDGES)), 3):
        if all(frozenset(pair) not in WEDGES for pair in combinations(trial, 2)):
            matching = trial
            break
    check(matching is not None and engine.coefficient(matching, 6) == ZERO,
          "representative three-edge matching h6 owner vanishes as a complete source Jet")

    dependency_hashes = frozen_regressions(diagonal)

    # The canonical cycle computation is exact for arbitrary nonuniform raw
    # source coordinates.  Run every active Q4 transition and freeze hashes.
    cycle_rows = []
    cycle_groups = {}
    for cycle_index, cycle in enumerate(ACTIVE_CYCLES):
        canonical, cycle_nodes = cycle_canonical_jet(cycle, BASE)
        off = canonical[6][0][1]
        check(off.value == -F(63, 8),
              "every active cycle has frozen source-off writer coefficient -63/8")
        expected_gradient = {}
        for node in cycle_nodes:
            local_ports = sorted(EDGES[edge][1] for edge in cycle if node in ENDS[edge])
            check(len(local_ports) == 2,
                  "cycle meets every cycle node in two distinct ports")
            pair = PAIR_INDEX[tuple(local_ports)]
            expected_gradient[6 * NODE_INDEX[node] + pair] = F(105, 8)
        check(off.grad == expected_gradient,
              "every active cycle has the complete frozen CH raw-pair first writer")
        check(bool(off.hess),
              "canonical h6 cycle writer has a nonzero mixed second-source Jet")
        record_stencil(cycle_groups, "alternating_cycle_offdiag", 6, cycle, off,
                       ordered_nodes=cycle_nodes)
        cycle_rows.append({
            "index": cycle_index,
            "cycle": list(cycle),
            "nodes": [NODE_INDEX[node] for node in cycle_nodes],
            "offdiag_value": qtext(off.value),
            "gradient_entries": len(off.grad),
            "hessian_entries": len(off.hess),
            "offdiag_jet_sha256": jet_hash(off),
        })

    cycle_hessian_counts = sorted({row["hessian_entries"] for row in cycle_rows})
    diagonal_stencils = stencil_ledger(stencil_groups)
    cycle_stencils = stencil_ledger(cycle_groups)
    diagonal_owners_computed = 256 + 768 + 512 + 2304
    selected_row_jets_computed = diagonal_owners_computed + len(cycle_rows)
    graph_wide_geometric_supports = diagonal_owners_computed + len(CYCLES)
    inactive_cycle_supports = len(CYCLES) - len(cycle_rows)
    all_stencil_records = list(diagonal_stencils.values()) + [
        cycle_stencils["alternating_cycle_offdiag_h6"]
    ]
    committed_class_count = sum(row["class_count"] for row in all_stencil_records)
    check(graph_wide_geometric_supports == 4096,
          "graph-wide geometric support census is 4096")
    check(selected_row_jets_computed == 3904,
          "computed selected-row nonbare source-Jet census is 3904")
    check(inactive_cycle_supports == 192,
          "192 graph-wide cycles are enumerated but inactive on selected row")
    check(committed_class_count == 1654,
          "compact ledger commits 1654 reconstructed stencil classes")
    check(len(all_stencil_records) == 8,
          "eight family/order records each retain one literal first stencil")
    result = {
        "schema": "GL6CU_COMPLETE_SIX_PAIR_H6_SOURCE_JET_V002",
        "status": "PASS",
        "claim": (
            "complete selected-branch and incident-cycle source-before-projection canonical "
            "raw-six-pair first and mixed-second source operator through h6 on the finite "
            "simple degree-four bipartite girth>=6 domain; literal full selected-Q4 "
            "connected-owner and active-cycle regression"
        ),
        "source_convention": {
            "hamiltonian": "H=U_d D+hW+sum_(v,ab) j_(v,ab) M_(v,ab)",
            "dimensionless_source": "eta=j/U_d",
            "fixed_objects": "P, z, M, and the occupation basis are source-independent",
            "retained_space": (
                "selected Q4 locked branch for diagonal owners; local two-endpoint P span "
                "for each incident active cycle; not the full global locked-manifold matrix"
            ),
            "pair_order": [list(pair) for pair in PAIR_ORDER],
            "jet": "value, first derivative, actual mixed second derivative at eta=0",
            "physical_first_scale": "h^(2m)/U_d^(2m) for K_(2m)'",
            "physical_second_scale": "h^(2m)/U_d^(2m+1) for K_(2m)''",
            "downstream_writer_sign": "for H=H0-JY, Y_A=-V_A where V_A=partial H_can/partial j_A",
            "canonical_chain_rule": (
                "full D2F[G_A,G_B]+DF[G_AB] of the affine microscopic source family; "
                "bare G_AB=0, induced denominator/matching/fold/canonical terms retained"
            ),
            "ctp_definition": (
                "Z[j+,j-]=Tr(T exp[-i integral H(j+)dt] rho antiT exp[+i integral H(j-)dt]); "
                "W=-i log Z; sigma=+1 forward and -1 backward"
            ),
            "ctp_branch_contact": (
                "direct same-branch connected contact is -sigma K''_AB delta(t-s); "
                "forward ++ is -K'' delta, backward -- is +K'' delta, mixed direct contacts vanish"
            ),
        },
        "owner_census": {
            "bare_node_sources": 128,
            "one_edge": 256,
            "two_edge_wedge": 768,
            "three_edge_star": 512,
            "three_edge_path": 2304,
            "elementary_six_cycle": 256,
            "total_nonbare_geometric_supports": 4096,
            "graph_wide_geometric_supports_enumerated": graph_wide_geometric_supports,
            "diagonal_nonbare_operator_jets_computed": diagonal_owners_computed,
            "selected_row_nonbare_operator_jets_computed": selected_row_jets_computed,
            "inactive_six_cycle_supports_enumerated_not_differentiated": inactive_cycle_supports,
            "active_alternating_cycle_transitions_in_selected_Q4_state": 64,
            "ownership": "Möbius inversion by distinct flipped-edge support; one undirected cycle owns both Hermitian directions",
        },
        "bare_endpoint_source_selected_Q4": {
            "value": qtext(bare.value),
            "gradient_entries": len(bare.grad),
            "hessian_entries": len(bare.hess),
            "jet_sha256": jet_hash(bare),
            "local_pair_word_census": [
                {"word": [qtext(value) for value in word], "count": count}
                for word, count in sorted(bare_word_census.items())
            ],
        },
        "diagonal_selected_Q4": {
            str(order): {
                "value": qtext(diagonal[order].value),
                "gradient_entries": len(diagonal[order].grad),
                "hessian_entries": len(diagonal[order].hess),
                "jet_sha256": jet_hash(diagonal[order]),
            } for order in (2, 4, 6)
        },
        "exact_local_operator_stencils": diagonal_stencils,
        "cycle_writer_selected_Q4": {
            "active_count": len(cycle_rows),
            "source_free_coefficient": "-63/8",
            "raw_first_coefficient_per_cycle_node": "(105/8)e_ab for the two cycle ports a,b",
            "hessian_entry_count_classes": cycle_hessian_counts,
            "rows_sha256": digest_payload(cycle_rows),
            "first_row": cycle_rows[0],
            "exact_local_operator_stencils": cycle_stencils,
        },
        "stencil_storage_contract": {
            "family_order_records": len(all_stencil_records),
            "all_reconstructed_local_stencil_classes_committed": committed_class_count,
            "literal_first_exact_stencils_stored": len(all_stencil_records),
            "representative_records_per_family_order_maximum": 4,
            "custody": (
                "compact ledger stores counts aggregate hashes representative metadata and one literal first exact stencil "
                "per family/order; sealed executable reconstructs all 1654 classes and verifies aggregate hashes"
            ),
        },
        "disconnected_owner_tests": {
            "two_edge_matching": "zero complete Jet",
            "wedge_plus_isolated": "zero complete Jet",
            "three_edge_matching": "zero complete Jet",
        },
        "frozen_dependency_sha256": dependency_hashes,
        "ceilings": [
            "operator/source/contact/writer census, not a stationary-state expectation",
            "complete around the selected Q4 branch, not a materialization of every row of the global locked manifold",
            "no connected-to-1PI Schur quotient, Ward identity, Ricci tensor, gravity, or G",
            "a later nonlinear physical-source embedding with its own nonzero G_AB is not included",
            "Q4 is the literal regulator; the general-domain proof uses only support parity, degree four, simplicity, and girth at least six",
        ],
        "checks": len(CHECKS),
    }
    result["result_sha256_without_runtime"] = digest_payload(result)
    elapsed = time.monotonic() - started
    return result, elapsed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result, elapsed = run()
    frozen = HERE / "EXACT_LEDGER.json"
    if frozen.is_file():
        if json.loads(frozen.read_text()) != result:
            raise AssertionError("fresh complete result differs from frozen exact ledger")
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        print(f"PASS__GL6CU_V002_COMPLETE_SIX_PAIR_H6_SOURCE_JET__{len(CHECKS)}/{len(CHECKS)}")
        print("CENSUS=4096_GEOMETRIC_SUPPORTS_ENUMERATED;3904_SELECTED_ROW_NONBARE_JETS_COMPUTED;192_INACTIVE_CYCLES")
        print("OWNERS=128_BARE+256_LINK+768_WEDGE+512_STAR+2304_PATH+64_ACTIVE_CYCLE_JETS")
        print("STENCILS=1654_CLASSES_HASH_COMMITTED_AND_EXECUTABLY_RECONSTRUCTED;8_LITERAL_FIRST_STENCILS_STORED")
        print("DIAGONAL_VALUES=" + ",".join(
            f"h{order}:{qtext(result['diagonal_selected_Q4'][str(order)]['value'])}"
            for order in (2, 4, 6)))
        print("ACTIVE_CYCLE_WRITERS=64;H6=-63/8;FIRST_RAW=105/8;SECOND=NONZERO")
        print("CTP=DEFINED_Z_PLUS_MINUS_AND_W_MINUS_I_LOG_Z;DIRECT_CONTACT=-SIGMA_ENERGY_HESSIAN_DELTA")
        print(f"RUNTIME_SECONDS={elapsed:.3f}")


if __name__ == "__main__":
    main()
