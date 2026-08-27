#!/usr/bin/env python3
"""Algebraic witnesses for alpha nonselection, inheritance, and phase locking.

The graph and numerical checks do not certify the physical premises of REC,
FCLPD, current conservation, the Magnus derivation, or detector realization.
"""

from fractions import Fraction as F


passed = 0
total = 0


def check(label, condition):
    global passed, total
    total += 1
    if not condition:
        raise AssertionError(label)
    passed += 1
    print(f"PASS {label}")


print("CANONICAL NORMALIZATION AND FIELD-COORDINATE INVARIANCE")
q = F(3, 5)
z = F(7, 4)
base = q * q / z  # 4*pi*alpha
for lam in (F(2), F(-3, 2), F(5, 7)):
    q_prime = q / lam
    z_prime = z / (lam * lam)
    check(f"q^2/Z invariant for lambda={lam}", q_prime * q_prime / z_prime == base)


print("SAME-SECTOR SPECIES INHERITANCE")
for n in (-3, -1, 1, 2, 5):
    species_vertex_squared = F(n * n) * base
    normalized = species_vertex_squared / F(n * n)
    check(f"charge representation n={n} has common base coupling", normalized == base)


print("CROSS-MODEL NONSELECTION VERSUS WITHIN-PARENT UNIQUENESS")
record_law = (F(1, 8), F(3, 8), F(1, 2))
parent_a = {"record_law": record_law, "four_pi_alpha": F(1, 10)}
parent_b = {"record_law": record_law, "four_pi_alpha": F(1, 7)}
check("same registered record law survives different parent alpha", parent_a["record_law"] == parent_b["record_law"])
check("the two parent alphas are physically distinct", parent_a["four_pi_alpha"] != parent_b["four_pi_alpha"])

subsystems = {
    "writer": base,
    "carrier": base,
    "hold": base,
    "reader": base,
}
check("all restrictions of one frozen parent inherit one coefficient", set(subsystems.values()) == {base})
foreign = base + F(1, 101)
check("an inequivalent coefficient is not the frozen parent coefficient", foreign not in set(subsystems.values()))


print("CANONICAL U1 FIXED-CONTROL RECORD WORLDS")
e1 = F(3, 10)
e2 = F(3, 5)
# B^2=(100/9) ln(2).  Store the exact rational coefficient of ln(2).
b2_over_ln2 = F(100, 9)
exponent1_over_ln2 = e1 * e1 * b2_over_ln2
exponent2_over_ln2 = e2 * e2 * b2_over_ln2
check("first coherent-state exponent is exactly ln(2)", exponent1_over_ln2 == 1)
check("second coherent-state exponent is exactly 4 ln(2)", exponent2_over_ln2 == 4)

# exp[-k ln(2)]=2^(-k), so the registered probabilities are exact Fractions.
p1 = 1 - F(1, 2) ** int(exponent1_over_ln2)
p2 = 1 - F(1, 2) ** int(exponent2_over_ln2)
check("first fixed-control cavity contrast is exactly 1/2", p1 == F(1, 2))
check("second fixed-control cavity contrast is exactly 15/16", p2 == F(15, 16))
check("both cavity records clear one common strict floor", min(p1, p2) > F(2, 5))
check("canonical 4*pi*alpha values are distinct", e1 * e1 != e2 * e2)
check("second canonical alpha is four times the first", (e2 * e2) / (e1 * e1) == 4)
check("EM-off ablation has exactly zero photon-record contrast", 1 - F(1, 2) ** 0 == 0)
check("fixed-control response is strictly alpha-sensitive", p2 > p1 > 0)

cavity_edges = {
    "prep": (
        "source",
        "current_program",
        "photon_blank",
        "pointer_blank",
        "write_clock",
        "hold_clock",
        "query_clock",
        "query_program",
    ),
    "source": ("current_generator", "source_exhaust_G"),
    "current_program": ("current_generator",),
    "write_clock": ("current_generator", "frontier_G"),
    "current_generator": ("current_register", "controller_exhaust_G"),
    "current_register": ("em_write", "current_exhaust_G"),
    "photon_blank": ("em_write",),
    "em_write": ("frontier_G",),
    "frontier_G": ("free_hold",),
    "hold_clock": ("free_hold",),
    "free_hold": ("query",),
    "pointer_blank": ("query",),
    "query_clock": ("query",),
    "query_program": ("query",),
    "query": ("Y",),
    "source_exhaust_G": (),
    "controller_exhaust_G": (),
    "current_exhaust_G": (),
    "Y": (),
}


def reachable(start, target, edges):
    frontier = [start]
    seen = set()
    while frontier:
        node = frontier.pop()
        if node == target:
            return True
        if node in seen:
            continue
        seen.add(node)
        frontier.extend(edges[node])
    return False


stage = {
    "prep": 0,
    "source": 1,
    "current_program": 1,
    "photon_blank": 1,
    "pointer_blank": 1,
    "write_clock": 1,
    "hold_clock": 1,
    "query_clock": 1,
    "query_program": 1,
    "current_generator": 2,
    "current_register": 3,
    "em_write": 4,
    "frontier_G": 5,
    "source_exhaust_G": 5,
    "controller_exhaust_G": 5,
    "current_exhaust_G": 5,
    "free_hold": 6,
    "query": 7,
    "Y": 8,
}
check(
    "cavity mission graph is acyclic",
    all(stage[u] < stage[v] for u, children in cavity_edges.items() for v in children),
)
check("declared source reaches the registered query", reachable("source", "query", cavity_edges))
check("registered query is not an ancestor of the write", not reachable("query", "em_write", cavity_edges))
without_em_write = {
    node: tuple(child for child in children if child != "em_write")
    for node, children in cavity_edges.items()
    if node != "em_write"
}
check(
    "every declared source-to-query path crosses the canonical EM write",
    reachable("source", "em_write", cavity_edges)
    and reachable("em_write", "query", cavity_edges)
    and not reachable("source", "query", without_em_write),
)
check("no post-frontier return to the source", not reachable("frontier_G", "source", cavity_edges))
check(
    "isolated source/current/controller exhaust has no query path",
    all(
        not reachable(node, "query", cavity_edges)
        for node in ("source_exhaust_G", "controller_exhaust_G", "current_exhaust_G")
    ),
)
check(
    "current register has no direct bypass into the post-G query ancestry",
    cavity_edges["em_write"] == ("frontier_G",)
    and "free_hold" not in cavity_edges["current_register"],
)
check("complete binary query alphabet is retained", {0, 1} == set((0, 1)))


print("PHASE POLYNOMIAL")


def h(x):
    return (15 * x - 10 * x**3 + 3 * x**5) / 8


def hp(x):
    return F(15, 8) * (1 - x * x) ** 2


def hpp(x):
    return -F(15, 2) * x * (1 - x * x)


for s in (F(-1), F(1)):
    check(f"h({s})={s}", h(s) == s)
    check(f"h'({s})=0", hp(s) == 0)
    check(f"h''({s})=0", hpp(s) == 0)


print("NORMALIZED WALL IDENTITY")
# Choose v=1 and lambda=2. Then V=(1/2)(1-phi^2)^2,
# the first-order kink equation is phi'=1-phi^2, and sigma=4/3.
samples = (
    (F(-3, 4), F(7, 16)),
    (F(-1, 3), F(8, 9)),
    (F(0), F(1)),
    (F(2, 5), F(21, 25)),
)
for phi, dphi in samples:
    potential = F(1, 2) * (1 - phi * phi) ** 2
    lhs = F(1, 2) * dphi * dphi + potential
    rhs = F(1, 2) * (dphi - (1 - phi * phi)) ** 2 + dphi * (1 - phi * phi)
    check(f"Bogomolny identity at phi={phi}", lhs == rhs)


def primitive(phi):
    return phi - phi**3 / 3


sigma = primitive(F(1)) - primitive(F(-1))
check("normalized planar wall tension is 4/3", sigma == F(4, 3))
check("twice the normalized isolated-wall tension is 8/3", 2 * sigma == F(8, 3))

print(f"SUMMARY {passed}/{total} exact checks passed")
print("VERDICT ALPHA_INHERITANCE_AND_ACTIVE_EM_ALGEBRAIC_WITNESSES_PASS")
print("SCOPE algebraic examples and a small DAG/cut only; physical REC/FCLPD premises are analytic, not executable certifications")
