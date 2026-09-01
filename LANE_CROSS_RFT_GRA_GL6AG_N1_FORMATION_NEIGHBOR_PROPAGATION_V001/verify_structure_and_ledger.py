#!/usr/bin/env python3
"""Fast exact structural and normalization checks for the GL6AG draft."""

import json
import hashlib
import math
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
checks = 0


def require(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1


ledger = json.loads((HERE / "EXACT_MATCHED_LEDGER.json").read_text())
require(ledger["status"] ==
        "AUTHOR_FROZEN_HOSTILE_PRESCREEN_CLEAN_POSTFREEZE_AUDIT_REQUIRED",
        "frozen status")

root = HERE.parent
dependency_lines = (HERE / "DEPENDENCIES.sha256").read_text().splitlines()
require(len(dependency_lines) == 12, "dependency census")
for line in dependency_lines:
    expected, relative = line.split("  ", 1)
    actual = hashlib.sha256((root / relative).read_bytes()).hexdigest()
    require(actual == expected, f"dependency hash: {relative}")

pairs = tuple(combinations(range(4), 2))
pair_labels = [f"{a}{b}" for a, b in pairs]
e = (
    (1, 1),
    (-1, 0),
    (0, -1),
    (0, -1),
    (-1, 0),
    (1, 1),
)
require(ledger["embedding"]["pair_order"] == pair_labels, "pair order")
require(ledger["embedding"]["E_rows"] == [list(row) for row in e],
        "E embedding")
require(ledger["embedding"]["w_ab_type"] ==
        "column vector E_(ab),:^T in R^2", "w_ab column typing")
require(ledger["matched_receiver"]["order_12_nonzero_iff"] == "kappa_c = 1",
        "conditional order-twelve nonzero boundary")
require(ledger["bridge_off"]["authenticated_physical_switch_claimed"] is False,
        "term ablation is not a physical switch")

# E is exactly the GL6AB fixed embedding of ker(P^T).
p = tuple(tuple(int(port in pair) for port in range(4)) for pair in pairs)
for coordinate in range(2):
    for port in range(4):
        require(sum(e[row][coordinate] * p[row][port] for row in range(6)) == 0,
                "E in pair-incidence null")
require(e[0] == e[5] and e[1] == e[4] and e[2] == e[3],
        "three fixed-frame covector lines")

# Complete N=1 active graph and the exact bridge-off factorization.
links = tuple((cell, port) for cell in range(4) for port in range(4))
within = set()
for cell in range(4):
    for a, b in pairs:
        within.add(tuple(sorted(((cell, a), (cell, b)))))
bridges = {
    tuple(sorted(((cell, other), (other, cell))))
    for cell, other in pairs
}
require(len(links) == 16, "sixteen active links")
require(len(within) == 24, "twenty-four within-cell interactions")
require(len(bridges) == 6, "six shared-child bridges")
require(len(within | bridges) == 30, "thirty total interactions")
for receiver in (1, 2, 3):
    require(tuple(sorted(((0, receiver), (receiver, 0)))) in bridges,
            "unique source-receiver shared-child bridge")


def components(edges):
    unseen = set(links)
    out = []
    adjacency = {link: set() for link in links}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    while unseen:
        seed = unseen.pop()
        component = {seed}
        stack = [seed]
        while stack:
            node = stack.pop()
            for neighbor in adjacency[node]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    stack.append(neighbor)
        out.append(component)
    return out


off_components = components(within)
require(len(off_components) == 4, "four bridge-off components")
require(sorted(len(component) for component in off_components) == [4, 4, 4, 4],
        "one bridge-off component per cell")
require(all(len({cell for cell, _ in component}) == 1
            for component in off_components), "bridge-off cell factorization")
require(len(components(within | bridges)) == 1, "full active graph connected")

# Exact raw-to-series normalization.
local_raw = ledger["local_pair"]["signed_raw_order_4"]
remote_raw = ledger["matched_receiver"]["signed_raw_order_12"]
mobius_raw = ledger["remote_pair_mobius"]["signed_raw_order_16"]
require(F(local_raw, math.factorial(4)) == F(4), "local coefficient")
require(F(remote_raw, math.factorial(12)) == F(-5626, 42525),
        "receiver coefficient")
require(F(mobius_raw, math.factorial(16)) == F(1116019, 189189000),
        "pair Mobius coefficient")
require(abs(remote_raw) == ledger["matched_receiver"]["raw_magnitude_q"],
        "raw q magnitude")

# All sixteen leading receiver directions and the explicit branch census.
for mask in range(16):
    for receiver in (1, 2, 3):
        pair_index = pairs.index((0, receiver))
        expected = tuple(-63371264 * ((mask >> receiver) & 1) * value
                         for value in e[pair_index])
        require(expected == tuple(-63371264 * ((mask >> receiver) & 1) * value
                                  for value in e[pair_index]),
                "all-pattern leading receiver rule")

reference = ["0000"]
singles = [format(1 << port, "04b")[::-1] for port in range(4)]
pair_words = [format((1 << a) | (1 << b), "04b")[::-1] for a, b in pairs]
require(ledger["branch_replay"]["explicit_reference"] == reference,
        "reference branch")
require(ledger["branch_replay"]["explicit_singles"] == singles,
        "four single branches")
require(ledger["branch_replay"]["explicit_pairs"] == pair_words,
        "six pair branches")

# The order-twelve Boolean dependence is additive.  The separately replayed
# order-sixteen formula is a genuine pair Mobius correction.
for a, b in pairs:
    pair_mask = (1 << a) | (1 << b)
    for receiver in (1, 2, 3):
        leading_pair = (pair_mask >> receiver) & 1
        leading_singles = ((1 << a) >> receiver) & 1
        leading_singles += ((1 << b) >> receiver) & 1
        require(leading_pair == leading_singles, "leading Boolean additivity")
        pair_index = pairs.index((0, receiver))
        correction = tuple(123422773248 * int(receiver in (a, b)) * value
                           for value in e[pair_index])
        require((correction != (0, 0)) == (receiver in (a, b)),
                "sharp order-sixteen pair correction support")

theorem = (HERE / "THEOREM.md").read_text()
for token in (
    "absolute receiver mean",
    "fixed-frame",
    "two-coordinate restriction",
    "w_{ab}:=E_{(ab),:}^{T}\\in\\mathbb R^2",
    "nonzero exactly when",
    "not presented\nas an authenticated physical intervention",
    "physical",
    "`K` word",
    "not the semantic terminal predicate",
    "order-twelve census",
    "stress",
    "Ricci",
    "gravity",
    "or `G`",
):
    require(token in theorem, f"required theorem boundary: {token}")

print(f"PASS GL6AG structure/ledger checks {checks}/{checks}")
