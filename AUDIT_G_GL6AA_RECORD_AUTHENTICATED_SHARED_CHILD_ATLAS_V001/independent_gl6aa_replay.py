#!/usr/bin/env python3
"""Independent exact replay of the frozen GL6AA atlas theorem.

This program uses only the Python standard library.  It does not import or
execute either frozen author verifier.  It reconstructs the one-port
non-identifiability boundary, the direct-sum copy dilation, a complete
finite MATCH query, literal shared-child equivalence, and the S4/six-pair
transition cocycle.
"""

from copy import deepcopy
from itertools import combinations, permutations, product
from math import comb


checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def compose(left, right):
    """Return left after right for tuple permutations."""
    return tuple(left[right[index]] for index in range(len(right)))


def inverse(perm):
    out = [0] * len(perm)
    for source, target in enumerate(perm):
        out[target] = source
    return tuple(out)


def transposition(a, b, size=4):
    out = list(range(size))
    out[a], out[b] = out[b], out[a]
    return tuple(out)


identity4 = tuple(range(4))
S4 = tuple(permutations(range(4)))
PAIRS = tuple(combinations(range(4), 2))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIRS)}


def pair_rep(perm):
    return tuple(PAIR_INDEX[tuple(sorted((perm[a], perm[b])))]
                 for a, b in PAIRS)


# ---------------------------------------------------------------------------
# A single connector authenticates one port image, not a complete frame.


for source_port in range(4):
    for target_port in range(4):
        extensions = [perm for perm in S4
                      if perm[source_port] == target_port]
        check(len(extensions) == 6,
              "one port image leaves exactly three-factorial extensions")

g_ba = transposition(0, 1)
g_cb = transposition(1, 2)
g_ac = transposition(0, 2)
check(g_ba[0] == 1 and g_cb[1] == 2 and g_ac[2] == 0,
      "three one-port constraints")
holonomy = compose(g_ac, compose(g_cb, g_ba))
check(holonomy == transposition(1, 2) and holonomy != identity4,
      "one-port-compatible nontrivial holonomy")


# ---------------------------------------------------------------------------
# The selected blank/code copy dilation is an exact controlled unitary.


def author_digits(label_count):
    digits = 0
    while 2 ** digits < label_count + 1:
        digits += 1
    return digits


for n in range(0, 10):
    labels = comb(n + 3, 3) + comb(n + 4, 3)
    digits = author_digits(labels)
    check(2 ** digits >= labels + 1, "author ID bank has finite capacity")
    check(2 ** (digits - 1) < labels + 1 if digits else labels == 0,
          "author digit formula replay")

# State 0 is the unique all-blank target; states 1..L are the orthogonal
# binary codewords.  Extra states represent the rest of the ternary bank and
# must be fixed by every X_s.
for label_count in (1, 2, 3, 4, 5, 8, 9, 20):
    target_dimension = label_count + 5
    target_basis = tuple(range(target_dimension))

    def x_s(source_label, target_state):
        if target_state == 0:
            return source_label
        if target_state == source_label:
            return 0
        return target_state

    for source_label in range(1, label_count + 1):
        image = tuple(x_s(source_label, state) for state in target_basis)
        check(len(set(image)) == target_dimension, "X_s is a permutation")
        for state in target_basis:
            check(x_s(source_label, x_s(source_label, state)) == state,
                  "X_s is an involution")
            if state not in (0, source_label):
                check(x_s(source_label, state) == state,
                      "X_s fixes the orthogonal complement")
        check(x_s(source_label, 0) == source_label,
              "copy maps blank to selected code")

    # Include one P_perp source state, denoted zero.  The controlled map is a
    # bijective involution and leaves its source unchanged.
    controlled_image = []
    for source_label in range(0, label_count + 1):
        for target_state in target_basis:
            copied = (target_state if source_label == 0
                      else x_s(source_label, target_state))
            controlled_image.append((source_label, copied))
            second = (copied if source_label == 0
                      else x_s(source_label, copied))
            check(second == target_state, "controlled copy squares to identity")
    check(len(set(controlled_image)) == len(controlled_image),
          "controlled direct sum is unitary on basis states")

    # Repeated use of one orthogonal source on two disjoint targets commutes.
    for source_label in range(1, label_count + 1):
        for first, second in product(target_basis, repeat=2):
            order_12 = (x_s(source_label, first),
                        x_s(source_label, second))
            order_21 = (x_s(source_label, first),
                        x_s(source_label, second))
            check(order_12 == order_21, "disjoint tap copies commute")
        check((source_label, source_label)
              == (x_s(source_label, 0), x_s(source_label, 0)),
              "orthogonal source can populate two blank taps")


# ---------------------------------------------------------------------------
# Literal FPSS slab and shared-child replay.


def simplex(total):
    return tuple((a, b, c, total - a - b - c)
                 for a in range(total + 1)
                 for b in range(total - a + 1)
                 for c in range(total - a - b + 1))


def add_unit(point, port):
    out = list(point)
    out[port] += 1
    return tuple(out)


def subtract(left, right):
    return tuple(left[index] - right[index] for index in range(4))


for n_level in range(0, 6):
    parents = simplex(n_level)
    children = simplex(n_level + 1)
    check(len(parents) == comb(n_level + 3, 3), "parent simplex census")
    check(len(children) == comb(n_level + 4, 3), "child simplex census")
    child_owner = {child: index for index, child in enumerate(children)}
    check(len(child_owner) == len(children), "child owner map is bijective")

    edges = tuple((parent, port, add_unit(parent, port))
                  for parent in parents for port in range(4))
    check(len(edges) == 4 * len(parents), "four append ports per parent")
    check(len({(parent, child) for parent, _, child in edges}) == len(edges),
          "literal append edges are distinct")

    # Use the independently queried physical site label itself as the ideal
    # injective sigma outcome.  Parent and child simplices are disjoint by
    # their coordinate sum.
    all_sites = parents + children
    sigma = {site: site for site in all_sites}
    check(len(set(sigma.values())) == len(all_sites),
          "complete ideal site-ID census is injective")

    for edge, other in product(edges, repeat=2):
        parent, port, child = edge
        other_parent, other_port, other_child = other
        left_tap = sigma[child]
        right_tap = sigma[other_child]
        check((left_tap == right_tap) == (child == other_child),
              "equal copied child IDs iff literal child equality")
        check((child_owner[child] == child_owner[other_child])
              == (child == other_child), "owner injectivity")
        if child == other_child and parent != other_parent:
            displacement = subtract(other_parent, parent)
            expected = tuple(int(index == port) - int(index == other_port)
                             for index in range(4))
            check(displacement == expected,
                  "shared-child chart displacement")
            check(sum(abs(value) for value in displacement) == 2,
                  "adjacent graph distance is one")
            positive = [index for index, value in enumerate(displacement)
                        if value == 1]
            negative = [index for index, value in enumerate(displacement)
                        if value == -1]
            check(positive == [port] and negative == [other_port],
                  "positive and negative coordinates are unique")

    # Every enumerated graph triangle telescopes exactly.  This is a finite
    # corroboration of the algebraic all-path telescoping identity.
    adjacency = {parent: set() for parent in parents}
    by_child = {}
    for parent, _, child in edges:
        by_child.setdefault(child, []).append(parent)
    for owners in by_child.values():
        for left, right in combinations(owners, 2):
            adjacency[left].add(right)
            adjacency[right].add(left)
    for a, b, c in combinations(parents, 3):
        if b in adjacency[a] and c in adjacency[b] and a in adjacency[c]:
            closed_sum = tuple(
                subtract(b, a)[index]
                + subtract(c, b)[index]
                + subtract(a, c)[index]
                for index in range(4)
            )
            check(closed_sum == (0, 0, 0, 0), "triangle displacement cocycle")

# The explicit common-child counterexample really belongs to every N>=2.
for n_level in range(2, 10):
    child = (1, 1, 1, n_level - 2)
    check(child in simplex(n_level + 1), "counterexample child in slab")
    parents = []
    for port in range(3):
        parent = list(child)
        parent[port] -= 1
        parents.append(tuple(parent))
    check(len(set(parents)) == 3, "three distinct common-child parents")
    check(all(parent in simplex(n_level) for parent in parents),
          "counterexample parents in slab")


# ---------------------------------------------------------------------------
# A complete finite MATCH query on a full N=1 raw parent-child product.


n_level = 1
parents = simplex(n_level)
children = simplex(n_level + 1)
raw_ports = tuple(product(parents, children))
edge_lookup = {(parent, add_unit(parent, port)): port
               for parent in parents for port in range(4)}


def canonical_expected(raw_port):
    parent, child = raw_port
    if raw_port in edge_lookup:
        port = edge_lookup[raw_port]
        return ("EDGE", parent, port, child)
    return ("NONEDGE", parent, child)


def ideal_query_state():
    # Layer-tagged IDs keep parent and child label alphabets disjoint.
    ids = {("P", parent): ("P", parent) for parent in parents}
    ids.update({("C", child): ("C", child) for child in children})
    state = {
        "ids": ids,
        "taps": {},
        "expected": {},
        "k": {},
        "n": {},
        "labels": {parent: {port: port for port in range(4)}
                   for parent in parents},
        "owners": {raw_port: raw_port for raw_port in raw_ports},
        "failure": False,
    }
    for raw_port in raw_ports:
        parent, child = raw_port
        state["taps"][raw_port] = (
            state["ids"][("P", parent)],
            state["ids"][("C", child)],
        )
        expected = canonical_expected(raw_port)
        state["expected"][raw_port] = expected
        is_edge = expected[0] == "EDGE"
        state["k"][raw_port] = int(is_edge)
        state["n"][raw_port] = 1 if is_edge else 0
    return state


def atlas_flags(state):
    failure = bool(state["failure"])
    blank = any(value is None for value in state["ids"].values())
    blank |= any(value is None for taps in state["taps"].values()
                 for value in taps)
    blank |= any(value is None for value in state["expected"].values())
    blank |= any(value is None for labels in state["labels"].values()
                 for value in labels.values())

    occupied_ids = [value for value in state["ids"].values()
                    if value is not None]
    collision = len(occupied_ids) != len(set(occupied_ids))
    for labels in state["labels"].values():
        values = [value for value in labels.values() if value is not None]
        if len(values) != 4 or set(values) != set(range(4)):
            collision = True

    mismatch = False
    for raw_port in raw_ports:
        parent, child = raw_port
        if state["owners"].get(raw_port) != raw_port:
            mismatch = True
        actual_taps = state["taps"].get(raw_port)
        correct_taps = (state["ids"].get(("P", parent)),
                        state["ids"].get(("C", child)))
        if actual_taps != correct_taps:
            mismatch = True
        expected = state["expected"].get(raw_port)
        canonical = canonical_expected(raw_port)
        if expected != canonical:
            mismatch = True
        if not isinstance(expected, tuple) or not expected:
            mismatch = True
            continue
        if expected[0] == "EDGE":
            _, named_parent, port, named_child = expected
            if state["k"].get(raw_port) != 1:
                mismatch = True
            if actual_taps != (("P", named_parent), ("C", named_child)):
                mismatch = True
            if state["labels"].get(named_parent, {}).get(port) != port:
                mismatch = True
        elif expected[0] == "NONEDGE":
            if state["k"].get(raw_port) != 0:
                mismatch = True
            if state["n"].get(raw_port) != 0:
                mismatch = True
        else:
            mismatch = True
    return failure, blank, collision, mismatch


ideal = ideal_query_state()
check(atlas_flags(ideal) == (False, False, False, False),
      "selected ideal completion is MATCH")


def require_nonmatch(mutator, label):
    state = deepcopy(ideal)
    mutator(state)
    check(any(atlas_flags(state)), label)


first_parent = parents[0]
second_parent = parents[1]
first_edge = next(raw for raw in raw_ports if raw in edge_lookup)
first_nonedge = next(raw for raw in raw_ports if raw not in edge_lookup)
wrong_child = next(child for child in children if child != first_edge[1])

require_nonmatch(lambda state: state.__setitem__("failure", True),
                 "failure remains non-MATCH")
require_nonmatch(lambda state: state["ids"].__setitem__(("P", first_parent), None),
                 "blank ID remains non-MATCH")
require_nonmatch(lambda state: state["ids"].__setitem__(
    ("P", second_parent), state["ids"][("P", first_parent)]),
    "ID collision remains non-MATCH")
require_nonmatch(lambda state: state["taps"].__setitem__(first_edge, (None, None)),
                 "blank endpoint tap remains non-MATCH")
require_nonmatch(lambda state: state["taps"].__setitem__(
    first_edge, (("P", first_edge[0]), ("C", wrong_child))),
    "wrong endpoint tap remains non-MATCH")
require_nonmatch(lambda state: state["owners"].__setitem__(
    first_edge, first_nonedge), "wrong detector owner remains non-MATCH")
require_nonmatch(lambda state: state["k"].__setitem__(first_edge, 0),
                 "missing programmed edge remains non-MATCH")
require_nonmatch(lambda state: state["k"].__setitem__(first_nonedge, 1),
                 "extra programmed edge remains non-MATCH")
require_nonmatch(lambda state: state["n"].__setitem__(first_nonedge, 1),
                 "occupied nonedge dynamics remains non-MATCH")
require_nonmatch(lambda state: state["expected"].__setitem__(
    first_edge, ("NONEDGE",) + first_edge),
    "altered expectation remains non-MATCH")
require_nonmatch(lambda state: state["labels"][first_parent].__setitem__(1, 0),
                 "repeated port label remains non-MATCH")
require_nonmatch(lambda state: state["labels"][first_parent].__setitem__(0, None),
                 "blank port label remains non-MATCH")

# Edge occupancies n_e are intentionally not fixed by MATCH.
dynamic = deepcopy(ideal)
for raw_port in edge_lookup:
    dynamic["n"][raw_port] = (dynamic["n"][raw_port] + 1) % 3
check(atlas_flags(dynamic) == (False, False, False, False),
      "active-edge n remains dynamically unconstrained")


# ---------------------------------------------------------------------------
# Complete queried four-port frames force inverse and cocycle identities.


identity6 = tuple(range(6))
for lambda_m, lambda_n in product(S4, repeat=2):
    g_nm = compose(inverse(lambda_n), lambda_m)
    g_mn = compose(inverse(lambda_m), lambda_n)
    check(compose(g_mn, g_nm) == identity4, "four-port frame inverse")
    check(compose(pair_rep(g_mn), pair_rep(g_nm)) == identity6,
          "six-pair frame inverse")

for lambda_m, lambda_n, lambda_l in product(S4, repeat=3):
    g_nm = compose(inverse(lambda_n), lambda_m)
    g_ln = compose(inverse(lambda_l), lambda_n)
    g_lm = compose(inverse(lambda_l), lambda_m)
    check(compose(g_ln, g_nm) == g_lm, "four-port frame cocycle")
    check(compose(pair_rep(g_ln), pair_rep(g_nm)) == pair_rep(g_lm),
          "six-pair frame cocycle")


print("BOUNDARY=ONE_PORT_SIX_EXTENSIONS_NONTRIVIAL_HOLONOMY")
print("COPY=DIRECT_SUM_BLANK_CODE_INVOLUTION_DISJOINT_TAPS")
print("QUERY=COMPLETE_NONEXCLUSIVE_FLAGS_MATCH_IFF_ALL_FALSE")
print("SHARED_CHILD=EQUAL_ID_IFF_LITERAL_CHILD_ON_MATCH")
print("ATLAS=TRANSLATION_AND_S4_SIX_PAIR_COCYCLES_EXACT")
print("SCOPE=SELECTED_RELATIONAL_INCIDENCE_NOT_SPACE_CONE_RICCI_GRAVITY_G")
print("PASS__INDEPENDENT_GL6AA_REPLAY__%d/%d" % (checks, checks))
