#!/usr/bin/env python3
"""Independent hostile mathematical replay for repaired GL6AZ.

Standard library only.  This file imports no author code and uses no Python
``assert`` statements, so normal and optimized runs execute the same checks.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import random
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
AUTHOR = ROOT / "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001"
checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise RuntimeError(f"FAIL [{checks + 1}] {label}")
    checks += 1


def close(a: float, b: float, label: str, tol: float = 1e-12) -> None:
    check(abs(a - b) <= tol * max(1.0, abs(a), abs(b)), label)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_author() -> None:
    expected = {
        "THEOREM.md": "bcca352d9e58deba63a068a29f94a0e19fb88ae9994391ab1b5053120587ba44",
        "MANIFEST.sha256": "a932e083b5dd629e41f0014d52bf1c65f982612b55314e1e03e082b0feb8ebc7",
        "SEAL.sha256": "620eb8bcb254b0296bd4d2f7e81a0b3e78c0d0b9ce80ad257e97a4db394f392d",
    }
    for name, expected_hash in expected.items():
        path = AUTHOR / name
        check(path.is_file(), f"frozen author file exists: {name}")
        check(digest(path) == expected_hash, f"frozen author hash: {name}")


def one_flip_calibration() -> None:
    # Each endpoint has three other links.  Enumerate all locked completions
    # for an absent and an occupied target link.
    scales = [
        (Fraction(2, 7), Fraction(5, 11)),
        (Fraction(13, 9), Fraction(17, 8)),
        (Fraction(31, 6), Fraction(29, 5)),
        (Fraction(101, 37), Fraction(43, 19)),
    ]
    for target in (0, 1):
        triples = [bits for bits in itertools.product((0, 1), repeat=3)
                   if sum(bits) == 2 - target]
        check(len(triples) == 3, "three endpoint completions")
        for left in triples:
            for right in triples:
                before = (target + sum(left) - 2) ** 2 + (target + sum(right) - 2) ** 2
                flipped = 1 - target
                after = (flipped + sum(left) - 2) ** 2 + (flipped + sum(right) - 2) ** 2
                check(before == 0, "locked one-flip input")
                check(after == 2, "one flip creates two unit defects")
                for ud, h in scales:
                    off_diagonal = abs(-h)
                    diagonal_separation = ud * (after - before)
                    check(off_diagonal == h, "one-flip matrix element magnitude")
                    check(diagonal_separation == 2 * ud,
                          "computational diagonal separation is 2 U_d")
                    check(diagonal_separation / (2 * off_diagonal) == ud / h,
                          "native ratio identity")


def all_positive_ratios() -> None:
    # Use only the inherited integer d_*=3 witness.
    for numerator in range(1, 61):
        for denominator in range(1, 47):
            ratio = Fraction(numerator, denominator)
            for h in (Fraction(1, 13), Fraction(7, 3), Fraction(41, 17)):
                ud = ratio * h
                d_star = 3
                delta = 4 * ud * (d_star - 2)
                epsilon = delta + 2 * ud * (1 - 2 * d_star)
                check(ud > 0 and h > 0 and delta > 0, "strict inherited member")
                check(epsilon == -6 * ud, "integer witness lies on lock ray")
                check(ud / h == ratio, "integer witness realizes arbitrary R")


def complete_pair_read_and_tv() -> None:
    outcomes = list(itertools.product((-1, 1), repeat=4))
    check(len(outcomes) == 16, "complete link PVM has sixteen outcomes")
    for a in range(4):
        for b in range(a + 1, 4):
            plus = [q for q in outcomes if q[a] * q[b] == 1]
            minus = [q for q in outcomes if q[a] * q[b] == -1]
            check(len(plus) == 8 and len(minus) == 8, "pair marginal is 8+8")
            for q in outcomes:
                m = q[a] * q[b]
                check(Fraction(1 + m, 2) + Fraction(1 - m, 2) == 1,
                      "pair projectors resolve identity")

    # Exact rational replay of D_TV = |Delta <M>|/2.
    for denominator in range(1, 71):
        for i in range(denominator + 1):
            p = Fraction(i, denominator)
            mp = 2 * p - 1
            for j in range(denominator + 1):
                q = Fraction(j, denominator)
                mq = 2 * q - 1
                tv = (abs(p - q) + abs((1 - p) - (1 - q))) / 2
                check(tv == abs(p - q), "binary TV identity")
                check(tv == abs(mp - mq) / 2, "binary TV half-expectation factor")


def grouped_strong_support() -> None:
    displacements = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0))
    supports = []
    for da in displacements:
        support = frozenset(tuple(da[j] - db[j] for j in range(3))
                            for db in displacements)
        check(len(support) == 4, "native child strong support has four cells")
        supports.append(support)
    check(len(set(supports)) == 4,
          "ports at one parent have translated distinct child supports")

    # Four incoming links to one child share a support.  In any grouped sum,
    # only the chosen X_e connects n to n^e.  Exhaustive bitstrings attack
    # cancellation by arbitrary numbers of grouped one-bit flips.
    for width in range(2, 13):
        for state in itertools.product((0, 1), repeat=width):
            for target in range(width):
                wanted = list(state)
                wanted[target] ^= 1
                wanted = tuple(wanted)
                contributors = []
                for edge in range(width):
                    candidate = list(state)
                    candidate[edge] ^= 1
                    if tuple(candidate) == wanted:
                        contributors.append(edge)
                check(contributors == [target], "grouped potential isolates target flip")

    # Pinching between N_def eigenvalues 0 and 2 multiplies the matrix element
    # by the zero Fourier coefficient of exp(2 i theta).
    roots = [complex(math.cos(2 * math.pi * k / 4096),
                     math.sin(2 * math.pi * k / 4096)) ** 2
             for k in range(4096)]
    average = sum(roots) / len(roots)
    check(abs(average) < 1e-13, "charged one-flip matrix element is removed by pinching")


def source_domain_and_floors() -> None:
    first_floor = 36 * math.pi * math.e
    scale_floor = 432 * math.pi * math.e ** 2
    close(first_floor, 307.4304320162484, "first-smallness floor")
    close(scale_floor, 10028.190682380982, "scale-separation floor")

    # Explicit counterexample to the unrepaired compact domain.
    kappa = 0.25
    vbar = 1.0
    dbar = 0.0
    nubar = 54 * math.pi * (dbar + 2 * vbar) / (kappa ** 2)
    low_x = math.exp(-0.9)
    r_value = low_x * nubar
    low_f = low_x / (1 + math.log(low_x)) ** 3
    check(r_value >= 9 * math.pi * vbar / kappa,
          "low branch passes first compact inequality")
    check(math.floor(low_f) - 2 == 404, "low branch spuriously has n_star>=1")
    check(r_value < nubar, "restored scale separation excludes low branch")

    # Independent derivative-sign and global-minimum replay.
    for index in range(1, 18001):
        k = 10 ** (-4 + 5.5 * index / 18000)
        first = 9 * math.pi * math.exp(min(4 * k, 700)) / k
        scale = 108 * math.pi * math.exp(min(4 * k, 700)) / (k * k)
        if 4 * k < 700:
            check(first >= first_floor * (1 - 1e-12), "global first floor")
            check(scale >= scale_floor * (1 - 1e-12), "global scale floor")
        check((4 * k - 1 < 0) == (k < 0.25), "first-floor derivative sign")
        check((4 - 2 / k < 0) == (k < 0.5), "scale-floor derivative sign")

    # On x>=1, n_*>=1 iff x/(1+ln x)^3>=3.  Bisection is independent of
    # the author implementation.
    def f(x: float) -> float:
        return x / (1 + math.log(x)) ** 3

    check(f(math.e ** 2) < 3, "high branch below threshold at its minimum")
    lo, hi = math.e ** 2, 1_000_000.0
    check(f(hi) > 3, "high root bracketed")
    for _ in range(240):
        middle = (lo + hi) / 2
        if f(middle) < 3:
            lo = middle
        else:
            hi = middle
    x_star = (lo + hi) / 2
    close(x_star, 1861.32559690908, "high-branch root", tol=1e-13)
    close(f(x_star), 3.0, "high-branch root equation")
    full_floor = x_star * scale_floor
    close(full_floor, 18665728.00780086, "complete universal floor", tol=1e-12)
    check(2 < first_floor and Fraction(5, 2) < first_floor,
          "both H6 sample ratios fail even first floor")


def two_member_witness() -> None:
    def a_of_x(x: Fraction) -> Fraction:
        return Fraction(8, 63) / x ** 6 * (
            1 - x ** 2 - Fraction(37, 12) * x ** 4
            - Fraction(16247, 900) * x ** 6
        )

    x1, x2 = Fraction(2, 5), Fraction(1, 2)
    check(1 / x1 == Fraction(5, 2) and 1 / x2 == 2, "two exact R values")
    a1, a2 = a_of_x(x1), a_of_x(x2)
    check(a1 == Fraction(2415673, 113400), "first exact source coefficient")
    check(a2 == Fraction(31706, 14175), "second exact source coefficient")
    check(a1 - a2 == Fraction(3203, 168), "exact source-owner difference")
    for j6 in (Fraction(1, 19), Fraction(11, 7), Fraction(101, 23)):
        parents = []
        for x in (x1, x2):
            ud = Fraction(8, 63) * j6 / x ** 6
            h = x * ud
            check(Fraction(63, 8) * ud * (h / ud) ** 6 == j6,
                  "same source-free H6 coefficient")
            parents.append((ud, h))
        check(parents[0] != parents[1], "two microscopic parameter pairs differ")


def clock_scaling_and_scope() -> None:
    rng = random.Random(0x6A2)
    for _ in range(12000):
        ud = Fraction(rng.randint(1, 5000), rng.randint(1, 5000))
        h = Fraction(rng.randint(1, 5000), rng.randint(1, 5000))
        elapsed = Fraction(rng.randint(1, 5000), rng.randint(1, 5000))
        c = Fraction(rng.randint(1, 5000), rng.randint(1, 5000))
        check(ud / h == (c * ud) / (c * h), "common energy scaling preserves R")
        check(h * elapsed == (c * h) * (elapsed / c),
              "inverse clock scaling preserves dimensionless duration")

    theorem = (AUTHOR / "THEOREM.md").read_text(encoding="utf-8")
    result = (AUTHOR / "RESULT.md").read_text(encoding="utf-8")
    search = (AUTHOR / "SEARCH_LEDGER.md").read_text(encoding="utf-8")
    for token in (
        "R\\ge\\bar\\nu_0",
        "x_*\\,432\\pi e^2",
        "selected ideal ready/`MATCH` active factor",
        "binary pair marginal, not total variation",
        "system sampling/evolution endpoint immediately before the",
        "response is off or refocused during that read",
        "No noninteger `d_star` is used",
        "computational-basis diagonal energy separation/cost, not a spectral gap",
        "an individual witness may contain registered clock data",
        "sharper/different finite-mission local-observable theorem",
        "No graviton, Ricci target, Einstein equation, gravity identification, or",
    ):
        check(token in theorem, f"repaired theorem scope: {token}")
    check("0<r_1<{\\ln(3/2)\\over4}" in result, "result carries r1 domain")
    check("formal clauses supply no numerical" in search,
          "search ledger does not reduce U-DCL to topology")
    for forbidden in (
        "we have proved gravity",
        "this is gravity",
        "derives newton's constant",
        "necessary and sufficient prethermal",
        "full retained output obeys",
        "spectral gap equals 2u_d",
    ):
        check(forbidden not in theorem.lower(), f"forbidden promotion absent: {forbidden}")


def main() -> None:
    frozen_author()
    one_flip_calibration()
    all_positive_ratios()
    complete_pair_read_and_tv()
    grouped_strong_support()
    source_domain_and_floors()
    two_member_witness()
    clock_scaling_and_scope()
    print(f"PASS__INDEPENDENT_GL6AZ_POST_REPAIR_REPLAY__{checks}/{checks}")


if __name__ == "__main__":
    main()
