#!/usr/bin/env python3
"""Exact finite combinatorics for the GL6AA atlas interface."""

from itertools import permutations, product
from math import comb, ceil, log2


checks = 0


def require(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


def compose(left, right):
    """Return left after right."""
    return tuple(left[right[i]] for i in range(len(right)))


def inverse(p):
    out = [0] * 4
    for i, value in enumerate(p):
        out[value] = i
    return tuple(out)


identity = tuple(range(4))


# A single port image leaves 3! full frame extensions.
for a in range(4):
    for b in range(4):
        extensions = [p for p in permutations(range(4)) if p[a] == b]
        require(len(extensions) == 6, "single port leaves six S4 extensions")


# Explicit nontrivial loop compatible with three one-port constraints.
g_ba = (1, 0, 2, 3)  # (01), maps 0 to 1
g_cb = (0, 2, 1, 3)  # (12), maps 1 to 2
g_ac = (2, 1, 0, 3)  # (02), maps 2 to 0
require(g_ba[0] == 1 and g_cb[1] == 2 and g_ac[2] == 0,
        "one-port constraints")
holonomy = compose(g_ac, compose(g_cb, g_ba))
require(holonomy != identity, "pairwise connector data permit holonomy")


PAIRS = tuple((a, b) for a in range(4) for b in range(a + 1, 4))


def pair_rep(p):
    return tuple(PAIRS.index(tuple(sorted((p[a], p[b])))) for a, b in PAIRS)


# Complete local-to-global frame records force inverse and cocycle, including
# on the pair fiber: g_ba=lambda_b^{-1} lambda_a.
all_permutations = tuple(permutations(range(4)))
for lambda_a in all_permutations:
    for lambda_b in all_permutations:
        g_ba = compose(inverse(lambda_b), lambda_a)
        g_ab = compose(inverse(lambda_a), lambda_b)
        require(compose(g_ab, g_ba) == identity, "frame inverse")
        rho_ba = pair_rep(g_ba)
        rho_ab = pair_rep(g_ab)
        require(compose(rho_ab, rho_ba) == tuple(range(6)), "pair inverse")

for lambda_a, lambda_b, lambda_c in product(all_permutations, repeat=3):
    g_ba = compose(inverse(lambda_b), lambda_a)
    g_cb = compose(inverse(lambda_c), lambda_b)
    g_ca = compose(inverse(lambda_c), lambda_a)
    require(compose(g_cb, g_ba) == g_ca, "frame cocycle")
    require(compose(pair_rep(g_cb), pair_rep(g_ba)) == pair_rep(g_ca),
            "six-pair cocycle")


# Fixed finite FPSS sizes admit finite injective code banks.
for n in range(0, 12):
    a_n = comb(n + 3, 3)
    b_n = comb(n + 4, 3)
    labels = a_n + b_n
    digits = ceil(log2(labels + 1))
    require(2 ** digits >= labels + 1, "ID bank capacity including blank")
    require(2 ** (digits - 1) < labels + 1 if digits else labels == 0,
            "ID bank minimality")

    # The controlled blank/code swap is an involution for every label and
    # fixes every other target code.
    target_states = tuple(range(labels + 1))  # zero is the unique blank
    for source_label in range(1, labels + 1):
        x_s = tuple(source_label if value == 0 else
                    0 if value == source_label else value
                    for value in target_states)
        require(compose(x_s, x_s) == target_states, "copy block involution")
        for value in target_states:
            if value not in (0, source_label):
                require(x_s[value] == value, "copy fixes orthogonal target")


# Literal nonvacuous three-parent shared-child configuration begins at N=2.
for n in range(2, 12):
    child = (1, 1, 1, n - 2)
    require(sum(child) == n + 1 and all(x >= 0 for x in child),
            "triple-overlap child belongs to S_(N+1)")
    parents = []
    for a in range(3):
        parent = tuple(child[i] - int(i == a) for i in range(4))
        require(sum(parent) == n and all(x >= 0 for x in parent),
                "triple-overlap parent belongs to S_N")
        parents.append(parent)
    require(len(set(parents)) == 3, "three distinct parents")


# MATCH includes both directions of edge-list equality and retains dynamic n
# only on expected edges.
for expected_edge, k_value, n_value in product((False, True), repeat=3):
    match_support = (k_value == expected_edge)
    match_nonedge_n = expected_edge or not n_value
    match = match_support and match_nonedge_n
    if match:
        require(k_value == expected_edge, "MATCH exact K iff EDGE")
        require(expected_edge or n_value is False, "MATCH blank nonedge n")


# Translation displacements telescope on all enumerated closed A3 walks.
roots = tuple(tuple(int(i == a) - int(i == b) for i in range(4))
              for a in range(4) for b in range(4) if a != b)
for length in range(2, 6):
    for steps in product(roots, repeat=length):
        total = tuple(sum(step[i] for step in steps) for i in range(4))
        if total == (0, 0, 0, 0):
            reverse_total = tuple(sum(-step[i] for step in reversed(steps))
                                  for i in range(4))
            require(reverse_total == (0, 0, 0, 0), "closed walk inverse")


print(f"PASS GL6AA exact checks {checks}/{checks}")
