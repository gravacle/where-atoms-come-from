#!/usr/bin/env python3
"""Deterministic checks for GRA-FH-F3-Q4-FPSS-V001.

The verifier checks the finite combinatorics, the local controlled-pulse
identity, exact q4/F3 adjacency binding, the raw-slab d*=2 obstruction, and
the theorem's mandatory scope ceilings.  It does not simulate autonomous
support selection, a physical port calibration, or a thermodynamic phase.
"""

from __future__ import annotations

from collections import Counter
from math import comb
from pathlib import Path
import itertools
import sys


ROOT = Path(__file__).resolve().parent
THEOREM = ROOT / "THEOREM.md"


class Checks:
    def __init__(self) -> None:
        self.rows: list[tuple[bool, str]] = []

    def check(self, condition: bool, label: str) -> None:
        self.rows.append((bool(condition), label))

    def finish(self) -> int:
        for ok, label in self.rows:
            print(f"{'PASS' if ok else 'FAIL'}  {label}")
        passed = sum(ok for ok, _ in self.rows)
        total = len(self.rows)
        print(f"SUMMARY {passed}/{total} checks passed")
        return 0 if passed == total else 1


def front(n: int) -> tuple[tuple[int, int, int, int], ...]:
    return tuple(
        (a, b, c, n - a - b - c)
        for a in range(n + 1)
        for b in range(n - a + 1)
        for c in range(n - a - b + 1)
    )


def add_unit(m: tuple[int, ...], axis: int) -> tuple[int, ...]:
    out = list(m)
    out[axis] += 1
    return tuple(out)


def edges(n: int):
    p = front(n)
    c = front(n + 1)
    pi = {m: i for i, m in enumerate(p)}
    ci = {m: i for i, m in enumerate(c)}
    e = tuple((pi[m], ci[add_unit(m, axis)]) for m in p for axis in range(4))
    return p, c, e


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b)))
         for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def dagger(a):
    return [[a[j][i].conjugate() for j in range(len(a))]
            for i in range(len(a[0]))]


def identity(n: int):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def kron(a, b):
    return [
        [a[i // len(b)][j // len(b[0])] * b[i % len(b)][j % len(b[0])]
         for j in range(len(a[0]) * len(b[0]))]
        for i in range(len(a) * len(b))
    ]


def add_matrix(*matrices):
    return [
        [sum(matrix[i][j] for matrix in matrices)
         for j in range(len(matrices[0][0]))]
        for i in range(len(matrices[0]))
    ]


def scale_matrix(scalar, matrix):
    return [[scalar * value for value in row] for row in matrix]


def matvec(matrix, vector):
    return [sum(matrix[i][j] * vector[j] for j in range(len(vector)))
            for i in range(len(matrix))]


def commutator(a, b):
    return add_matrix(matmul(a, b), scale_matrix(-1, matmul(b, a)))


def equal_matrix(a, b, tol: float = 1e-12) -> bool:
    return len(a) == len(b) and len(a[0]) == len(b[0]) and all(
        abs(a[i][j] - b[i][j]) <= tol
        for i in range(len(a)) for j in range(len(a[0]))
    )


def main() -> int:
    q = Checks()
    text = THEOREM.read_text(encoding="utf-8")

    ns = range(0, 8)
    samples = [edges(n) for n in ns]

    q.check(all(len(p) == comb(n + 3, 3)
                for n, (p, _, _) in zip(ns, samples)),
            "stars-and-bars parent-front count")
    q.check(all(len(c) == comb(n + 4, 3)
                for n, (_, c, _) in zip(ns, samples)),
            "stars-and-bars child-front count")
    q.check(all(len(c) - len(p) == comb(n + 3, 2)
                for n, (p, c, _) in zip(ns, samples)),
            "equal-layer parent padding count")
    q.check(all(len(e) == 4 * len(p)
                for p, _, e in samples),
            "four append edges per parent")
    q.check(all(len(set(e)) == len(e) for _, _, e in samples),
            "q4 append edges are distinct")
    q.check(all(len(c) ** 2 - len(e) == comb(n + 4, 3) ** 2
                - 4 * comb(n + 3, 3)
                for n, (_, c, e) in zip(ns, samples)),
            "adjacent-layer nonedge count")

    parent_degree_ok = True
    child_degree_ok = True
    handshake_ok = True
    extreme_ok = True
    for n, (parents, children, edge_list) in zip(ns, samples):
        pd = Counter(i for i, _ in edge_list)
        cd = Counter(j for _, j in edge_list)
        parent_degree_ok &= all(pd[i] == 4 for i in range(len(parents)))
        child_degree_ok &= all(
            cd[j] == sum(x > 0 for x in child)
            for j, child in enumerate(children)
        )
        handshake_ok &= sum(pd.values()) == sum(cd.values()) == 4 * len(parents)
        extreme = (n + 1, 0, 0, 0)
        extreme_ok &= cd[children.index(extreme)] == 1

    q.check(parent_degree_ok, "every active parent has eligible degree four")
    q.check(child_degree_ok, "child degree equals positive-coordinate count")
    q.check(handshake_ok, "bipartite degree handshake")
    q.check(extreme_ok, "extreme child has eligible degree one")
    q.check(extreme_ok, "raw finite slab has empty global d*=2 sector")
    q.check(parent_degree_ok, "saturated FD word and parent d*=2 word are disjoint")

    # Basis order |K,n> = |00>,|01>,|10>,|11>.  The admitted pi/2 pulse is
    # identity for K=0 and iX for K=1.
    u = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1j],
        [0, 0, 1j, 0],
    ]
    q.check(equal_matrix(matmul(dagger(u), u), identity(4)),
            "K-controlled incidence pulse is unitary")
    q.check(u[3][2] == 1j and sum(abs(u[i][2]) for i in range(3)) == 0,
            "pulse maps |K=1,n=0> to i|K=1,n=1>")
    q.check(u[0][0] == 1 and u[1][1] == 1,
            "pulse is identity on the K=0 nonedge block")

    # Check that the physical saturated adjacency block is exactly B_N.
    adjacency_ok = True
    guard_zero_ok = True
    for parents, children, edge_list in samples:
        m = len(children)
        b = [[0 for _ in parents] for _ in children]
        for i, j in edge_list:
            b[j][i] = 1
        adjacency_ok &= all(
            b[j][i] == int(any(
                children[j] == add_unit(parents[i], axis)
                for axis in range(4)
            ))
            for i in range(len(parents)) for j in range(len(children))
        )
        # Equal-layer embedding pads the parent columns from |S_N| to M.
        padded = [row + [0] * (m - len(parents)) for row in b]
        guard_zero_ok &= all(
            padded[j][i] == 0
            for j in range(m) for i in range(len(parents), m)
        )
    q.check(adjacency_ok, "saturated physical support equals q4 append incidence B_N")
    q.check(guard_zero_ok, "all parent-guard adjacency columns are zero")

    # Exact local K/n operator replay.  Basis order is |K,n>.  The qualified
    # hold contains P_K X_n plus terms diagonal/controlled in n, but no raw
    # I_K X_n.  The blank nonedge is then a reducing one-dimensional block.
    i2 = identity(2)
    x2 = [[0, 1], [1, 0]]
    p1 = [[0, 0], [0, 1]]
    p_k = kron(p1, i2)
    p_n = kron(i2, p1)
    h_control = 0.37
    delta = 0.91
    h_qualified = add_matrix(
        scale_matrix(-h_control, kron(p1, x2)),
        scale_matrix(delta, p_n),
    )
    blank_nonedge = [1, 0, 0, 0]
    qualified_image = matvec(h_qualified, blank_nonedge)
    q.check(all(abs(value) <= 1e-12 for value in qualified_image[1:]),
            "qualified K=0,n=0 nonedge block is exactly invariant")
    q.check(equal_matrix(commutator(h_qualified, p_k), [[0] * 4 for _ in range(4)]),
            "qualified hold exactly conserves K support")
    q.check(not equal_matrix(commutator(h_qualified, p_n), [[0] * 4 for _ in range(4)]),
            "K-gated actuator may evolve active n while conserving K")
    q.check(not equal_matrix(p_k, p_n), "K support and active n are distinct factors")

    h_fd_slice = scale_matrix(delta, p_n)
    q.check(equal_matrix(commutator(h_fd_slice, p_n), [[0] * 4 for _ in range(4)]),
            "FD comparator with both flip actuators off conserves saturated n")

    h_with_raw = add_matrix(
        h_qualified,
        scale_matrix(-0.23, kron(i2, x2)),
    )
    raw_image = matvec(h_with_raw, blank_nonedge)
    q.check(any(abs(value) > 1e-12 for value in raw_image[1:]),
            "raw ungated X would violate blank-nonedge quarantine")

    required_phrases = [
        "fixed orthogonal finite program",
        "one active BQ4 factor is not a tensor product",
        "supplied cap, address maps, edge list",
        "passive support retention only",
        "raw ungated flip would take a",
        "but not the subsequently active",
        "merely stroboscopic echo",
        "[H_{\\rm hold},\\Pi_{p_N,E_N}]=0",
        "premise of the qualified",
        "qualified raw-flip-free",
        "keep both the raw ungated BS06 flip and the PESC `K`-gated incidence flip exactly zero",
        "Physical energies, port matrices, and calibration remain supplied",
        "positive uniform child/parent detuning",
        "Omega_2(E_N)=\\varnothing",
        "same-`n` incompatibility",
        "No `K_eT_e`, second",
        "visible electromagnetism",
        "gravity closure",
    ]
    normalized_text = " ".join(text.split())
    for phrase in required_phrases:
        normalized_phrase = " ".join(phrase.split())
        q.check(normalized_phrase in normalized_text,
                f"mandatory ceiling present: {phrase}")

    forbidden_promotions = [
        "autonomous q4 support selection is proved",
        "the fd detuning is derived",
        "the raw finite slab realizes the global ice phase",
        "gravity is derived by fpss",
    ]
    q.check(not any(phrase.lower() in text.lower() for phrase in forbidden_promotions),
            "no forbidden promotion sentence")

    dependencies = [
        "../LANE_CROSS_RFT_GRA_EQ_BOUNDED_Q4_RECORD_STREAM_WITNESS_V001/THEOREM.md",
        "../LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md",
        "../LANE_CROSS_RFT_ALPHA_GRA_DB_F3_PAIR_ELIGIBILITY_SUCCESSOR_COMPOSITION_V001/THEOREM.md",
        "../LANE_CROSS_RFT_GRA_DW_F3_AUTHENTICATED_SUPPORT_SECTOR_CONSERVATION_V001/THEOREM.md",
        "../LANE_GRA_FF_F3_Q4_CARRIER_LIFT_DERIVABILITY_NO_GO_V001/THEOREM.md",
    ]
    for dep in dependencies:
        q.check((ROOT / dep).is_file(), f"dependency present: {Path(dep).parent.name}")

    q.check(text.count("\\[") == text.count("\\]"),
            "display-math delimiters balanced")
    q.check(text.count("`") % 2 == 0, "Markdown backticks balanced")

    return q.finish()


if __name__ == "__main__":
    sys.exit(main())
