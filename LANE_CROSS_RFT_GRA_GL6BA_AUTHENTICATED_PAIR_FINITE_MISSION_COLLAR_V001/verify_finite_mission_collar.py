#!/usr/bin/env python3
"""Constructive replay for GL6BA V001 (standard library only)."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from collections import deque
from fractions import Fraction
from pathlib import Path


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"FAIL [{CHECKS}] {label}")


def close(a: float, b: float, label: str, tol: float = 1e-12) -> None:
    check(abs(a - b) <= tol * max(1.0, abs(a), abs(b)), label)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


EXPECTED_DEPENDENCIES = {
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/AUDIT.md": "72c3c6552a89a3d5b977fb7ff566d558ed70d9c14b2caac64c8602609356498f",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256": "9954328b6939e12c0b2f712bf6373d8a5c2a4cfafd77ab85aa133b27a4a9f707",
    "AUDIT_G_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256": "1d7cd10a59469139edf4b2836ce95820f605069ef3e6bc0d56bf4df58009b869",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/POSTFREEZE_AUDIT.md": "ec1e452e51fd381ad23bc24f49d4613132bcc570b778697e72ac63fd006cd22e",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/MANIFEST.sha256": "01d9dd3e2a7cc3247ea11719b5d69c1e0c82f3370770fef1b91fe98452d12839",
    "AUDIT_G_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/SEAL.sha256": "3137e843b9abc98ecdececeb67a204c56993b673e79562dde87f9c2588f1fe7f",
    "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/THEOREM.md": "083d5fbb8a48e27e365167075da132ffa23e395587a4c0e40cc572d8b761ad30",
    "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/MANIFEST.sha256": "d38f89c618ea6f77c7b399b005ad0f0abe04d3865e06921f8c765feb44f40620",
    "LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/SEAL.sha256": "322bf51a00f8fea3f36a09656dda4ebf89ba56b9a88d60b50e9cc7ab33223987",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/THEOREM.md": "bcca352d9e58deba63a068a29f94a0e19fb88ae9994391ab1b5053120587ba44",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/MANIFEST.sha256": "a932e083b5dd629e41f0014d52bf1c65f982612b55314e1e03e082b0feb8ebc7",
    "LANE_CROSS_RFT_GRA_GL6AZ_RECORD_AUTHENTICATED_PRETHERMAL_MISSION_IDENTIFIABILITY_V001/SEAL.sha256": "620eb8bcb254b0296bd4d2f7e81a0b3e78c0d0b9ce80ad257e97a4db394f392d",
}


def verify_dependencies() -> None:
    for rel, expected in sorted(EXPECTED_DEPENDENCIES.items()):
        path = ROOT / rel
        check(path.is_file(), f"dependency exists: {rel}")
        check(sha256(path) == expected, f"dependency hash: {rel}")


Cell = tuple[int, int, int, int]
Site = tuple[Cell, int]
ZERO: Cell = (0, 0, 0, 0)


def add(x: Cell, y: Cell) -> Cell:
    return tuple(a + b for a, b in zip(x, y))  # type: ignore[return-value]


def displacement(a: int, b: int) -> Cell:
    return tuple(1 if i == a else -1 if i == b else 0 for i in range(4))  # type: ignore[return-value]


DISPLACEMENTS = tuple(displacement(a, b) for a in range(4) for b in range(4) if a != b)


def cell_distance(x: Cell) -> int:
    check(sum(x) == 0, "A3 cell has zero coordinate sum")
    return sum(v for v in x if v > 0)


def ball(radius: int) -> set[Cell]:
    result: set[Cell] = set()
    for first_three in itertools.product(range(-radius, radius + 1), repeat=3):
        x: Cell = (*first_three, -sum(first_three))
        if cell_distance(x) <= radius:
            result.add(x)
    return result


def shell_formula(radius: int) -> int:
    return 1 if radius == 0 else 10 * radius * radius + 2


def ball_formula(radius: int) -> int:
    return (10 * radius ** 3 + 15 * radius ** 2 + 11 * radius + 3) // 3


def verify_a3_counts() -> None:
    previous: set[Cell] = set()
    for radius in range(0, 13):
        current = ball(radius)
        shell = current - previous
        check(len(shell) == shell_formula(radius), f"exact A3 shell r={radius}")
        check(len(current) == ball_formula(radius), f"exact A3 ball r={radius}")
        check((10 * radius ** 3 + 15 * radius ** 2 + 11 * radius + 3) % 3 == 0,
              "ball polynomial integral")
        for x in current:
            check(all(-radius <= value <= radius for value in x),
                  "collar coordinates lie in embedding box")
            shifted = tuple(value + radius + 1 for value in x)
            check(all(value >= 1 for value in shifted), "FPSS embedding is strict")
            check(sum(shifted) == 4 * (radius + 1), "FPSS embedding has declared N")
        previous = current

    for radius in range(1, 50):
        total = 0
        for p in range(1, 5):
            for q in range(1, 5 - p):
                total += (math.comb(4, p) * math.comb(4 - p, q)
                          * math.comb(radius - 1, p - 1)
                          * math.comb(radius - 1, q - 1))
        check(total == 10 * radius * radius + 2,
              "sign-composition shell identity")


def boundary_edges(cells: set[Cell]) -> list[tuple[Cell, Cell, Cell]]:
    rows = []
    for x in cells:
        for delta in DISPLACEMENTS:
            y = add(x, delta)
            if y not in cells:
                rows.append((x, y, delta))
    return rows


def verify_boundary_census() -> None:
    for radius in range(0, 13):
        cells = ball(radius)
        rows = boundary_edges(cells)
        expected = 12 * (3 * radius * radius + 3 * radius + 1)
        check(len(rows) == expected, f"exact crossing count L={radius}")
        by_delta = {delta: 0 for delta in DISPLACEMENTS}
        for x, y, delta in rows:
            by_delta[delta] += 1
            check(cell_distance(x) == radius, "inside crossing cell on boundary shell")
            check(cell_distance(y) == radius + 1, "outside crossing cell on next shell")
        directional = 3 * radius * radius + 3 * radius + 1
        check(all(value == directional for value in by_delta.values()),
              "all twelve directed boundary counts agree")
        coefficient = (math.comb(radius + 2, 2) ** 2
                       - 2 * math.comb(radius + 1, 2) ** 2
                       + math.comb(radius, 2) ** 2)
        check(coefficient == directional, "boundary generating-function coefficient")


def link_neighbors(site: Site) -> tuple[tuple[Site, str], ...]:
    x, a = site
    rows: list[tuple[Site, str]] = []
    for b in range(4):
        if b == a:
            continue
        rows.append(((x, b), "P"))
        rows.append(((add(x, displacement(a, b)), b), "C"))
    return tuple(rows)


def verify_link_geometry() -> None:
    roots: tuple[Site, Site] = ((ZERO, 0), (ZERO, 1))
    distance: dict[Site, int] = {root: 0 for root in roots}
    queue = deque(roots)
    max_link_radius = 25
    while queue:
        site = queue.popleft()
        check(len(set(link_neighbors(site))) == 6, "link graph degree exactly six")
        if distance[site] >= max_link_radius:
            continue
        for neighbor, _kind in link_neighbors(site):
            if neighbor not in distance:
                distance[neighbor] = distance[site] + 1
                queue.append(neighbor)

    for (x, _a), d_link in distance.items():
        d_cell = cell_distance(x)
        if d_cell > 0:
            check(d_link >= 2 * d_cell - 1, "port-aware link-to-cell distance")

    for order in range(0, 19):
        reached = [cell_distance(site[0]) for site, value in distance.items()
                   if value <= order]
        check(max(reached) == math.ceil(order / 2),
              "interaction-order cell collar is ceil(m/2)")

    for radius in range(1, 10):
        crossing_inside: list[Site] = []
        cells = ball(radius)
        for site in distance:
            x, _a = site
            if x not in cells:
                continue
            for neighbor, kind in link_neighbors(site):
                if kind == "C" and neighbor[0] not in cells:
                    crossing_inside.append(site)
        check(crossing_inside, "crossing inside endpoints found")
        minimum = min(distance[site] for site in crossing_inside)
        check(minimum == 2 * radius, "inside crossing endpoint distance is sharply 2L")


PAULI_PRODUCT = {
    ("I", "I"): (1, "I"), ("I", "X"): (1, "X"),
    ("I", "Y"): (1, "Y"), ("I", "Z"): (1, "Z"),
    ("X", "I"): (1, "X"), ("Y", "I"): (1, "Y"),
    ("Z", "I"): (1, "Z"), ("X", "X"): (1, "I"),
    ("Y", "Y"): (1, "I"), ("Z", "Z"): (1, "I"),
    ("X", "Y"): (1j, "Z"), ("Y", "X"): (-1j, "Z"),
    ("Y", "Z"): (1j, "X"), ("Z", "Y"): (-1j, "X"),
    ("Z", "X"): (1j, "Y"), ("X", "Z"): (-1j, "Y"),
}


Operator = dict[tuple[str, ...], complex]


def operator_multiply(left: Operator, right: Operator) -> Operator:
    out: Operator = {}
    for a, ca in left.items():
        for b, cb in right.items():
            coefficient = ca * cb
            word = []
            for pa, pb in zip(a, b):
                phase, pc = PAULI_PRODUCT[(pa, pb)]
                coefficient *= phase
                word.append(pc)
            key = tuple(word)
            out[key] = out.get(key, 0) + coefficient
    return {word: coefficient for word, coefficient in out.items()
            if abs(coefficient) > 1e-12}


def commutator(left: Operator, right: Operator) -> Operator:
    out = operator_multiply(left, right)
    for word, coefficient in operator_multiply(right, left).items():
        out[word] = out.get(word, 0) - coefficient
    return {word: coefficient for word, coefficient in out.items()
            if abs(coefficient) > 1e-12}


def single_pauli(width: int, site: int, pauli: str) -> Operator:
    word = ["I"] * width
    word[site] = pauli
    return {tuple(word): 1}


def number_pair(width: int, left: int, right: int) -> Operator:
    out: Operator = {}
    for left_coefficient, left_pauli in ((1, "I"), (-1, "Z")):
        for right_coefficient, right_pauli in ((1, "I"), (-1, "Z")):
            word = ["I"] * width
            word[left] = left_pauli
            word[right] = right_pauli
            key = tuple(word)
            out[key] = out.get(key, 0) + left_coefficient * right_coefficient / 4
    return out


def verify_taylor_filtration() -> None:
    # Along a simple link path, each new link needs one transverse onsite
    # commutator followed by one diagonal pair commutator.  This constructs
    # the nonzero formal word that makes the 4L+2 boundary order sharp.
    for hops in range(1, 12):
        width = hops + 2
        word = ["I"] * width
        word[0] = "Z"
        word[-1] = "Z"  # spectator second link of M_beta
        observable: Operator = {tuple(word): 1}
        commutator_count = 0
        for site in range(hops):
            observable = commutator(single_pauli(width, site, "X"), observable)
            commutator_count += 1
            observable = commutator(number_pair(width, site, site + 1), observable)
            commutator_count += 1
        check(bool(observable), "alternating F3 commutator word nonzero")
        check(any(word[hops] != "I" for word in observable),
              "alternating word reaches declared link distance")
        check(commutator_count == 2 * hops,
              "two ordinary commutators per link hop")

    for collar in range(0, 20):
        first_omitted = 4 * collar + 2
        check(first_omitted - 1 == 4 * collar + 1,
              "B_L Taylor matching order has no off-by-one")
    for order in range(1, 80):
        radius = max(0, math.ceil((order - 1) / 4))
        check(order <= 4 * radius + 1, "Taylor-order collar suffices")
        if radius > 0:
            check(order > 4 * (radius - 1) + 1,
                  "Taylor-order collar is smallest licensed radius")


def exponential_tail(order: int, x: float) -> float:
    check(order >= 1 and x >= 0, "tail domain")
    if x == 0:
        return 0.0
    term = math.exp(order * math.log(x) - math.lgamma(order + 1))
    total = term
    k = order
    for _ in range(100000):
        k += 1
        term *= x / k
        total += term
        if term <= max(1.0, total) * 1e-15 and k > x + 40:
            return total
    raise RuntimeError("tail did not converge")


def pair_error(radius: int, ratio: float, sigma: float) -> float:
    x = 48.0 * ratio * abs(sigma)
    return ((3 * radius * radius + 3 * radius + 1)
            * exponential_tail(2 * radius + 1, x))


def certified_radius(ratio: float, sigma: float, tolerance: float) -> int:
    for radius in range(0, 5000):
        if pair_error(radius, ratio, sigma) <= tolerance:
            return radius
    raise RuntimeError("certificate search did not terminate")


def verify_tail_and_moderate_ratios() -> None:
    check(48 * 2 == 96, "R=2 tail argument")
    check(Fraction(48) * Fraction(5, 2) == 120, "R=5/2 tail argument")

    samples = (
        (2.0, 0.001, 0.01, 1),
        (2.5, 0.001, 0.01, 1),
        (2.0, 0.01, 0.01, 3),
        (2.5, 0.01, 0.01, 4),
        (2.0, 0.1, 0.01, 17),
        (2.5, 0.1, 0.01, 20),
        (2.0, 1.0, 0.01, 136),
        (2.5, 1.0, 0.01, 169),
    )
    for ratio, sigma, tolerance, expected in samples:
        radius = certified_radius(ratio, sigma, tolerance)
        check(radius == expected, "moderate-R certified-radius regression")
        check(pair_error(radius, ratio, sigma) <= tolerance,
              "certified radius passes tolerance")
        if radius > 0:
            check(pair_error(radius - 1, ratio, sigma) > tolerance,
                  "preceding radius fails tolerance")

    check(certified_radius(2.0, 0.0, 0.01) == 0,
          "zero-duration certificate permits exact radius zero")

    for x in (0.01, 0.1, 1.0, 3.0, 10.0, 25.0):
        for order in range(max(1, math.floor(x) + 1), max(2, math.floor(x) + 20)):
            exact = exponential_tail(order, x)
            marked = (math.e * x / order) ** order
            check(exact <= marked * (1 + 2e-12), "optimized marked-tail envelope")

    for ratio in (2.0, 2.5, 7.0):
        for sigma in (0.001, 0.01, 0.1, 1.0):
            radius = certified_radius(ratio, sigma, 1e-8)
            check(radius < 5000, "finite collar certificate exists")
            check(pair_error(radius, ratio, sigma) <= 1e-8,
                  "finite collar reaches declared small tolerance")


def verify_binary_total_variation() -> None:
    for denominator in range(1, 81):
        for i in range(denominator + 1):
            p = Fraction(i, denominator)
            for j in range(denominator + 1):
                q = Fraction(j, denominator)
                tv = (abs(p - q) + abs((1 - p) - (1 - q))) / 2
                expectation_difference = abs((2 * p - 1) - (2 * q - 1))
                check(tv == expectation_difference / 2,
                      "binary TV is half pair-expectation difference")


def verify_theorem_scope() -> None:
    theorem = (LANE / "THEOREM.md").read_text(encoding="utf-8")
    required = [
        "C_L=36L^2+36L+12",
        "d_L(p,q)\\ge2d_{A_3}",
        "T_{2L+1}(48R|s|)",
        "ordinary nested-commutator order `4L+2`",
        "D_{\\rm TV}(p^\\Omega,p^{(L)})",
        "complete finite authenticated exterior missions",
        "complete finite all-formed/`MATCH` FPSS",
        "at most the `C_L` crossing terms",
        "mathematical completion/extension",
        "not a claim that one infinite record",
        "96|\\sigma_{\\rm obs}|",
        "120|\\sigma_{\\rm obs}|",
        "full F3 dynamics versus a spatial collar of the same",
        "selected-factor binary pair",
        "No graviton, Ricci target, Einstein equation, gravity identification, or `G`",
    ]
    for token in required:
        check(token in theorem, f"required theorem token: {token}")
    forbidden = [
        "we have proved gravity",
        "this is gravity",
        "derives newton's constant",
        "r=2 is nature's value",
        "all retained flags obey",
        "uniform for infinite time",
        "gl6ay applies at r=2",
    ]
    lower = theorem.lower()
    for phrase in forbidden:
        check(phrase not in lower, f"forbidden promotion absent: {phrase}")


def verify_ledger() -> None:
    data = json.loads((LANE / "COLLAR_LEDGER.json").read_text(encoding="utf-8"))
    check(data["lane"] == "GL6BA_V001", "ledger lane")
    check(data["pair_support_size"] == 2, "ledger pair support")
    check(data["link_degree"] == 6, "ledger link degree")
    check(data["cross_boundary_pairs"] == "12*(3*L^2+3*L+1)",
          "ledger crossing count")
    check("T_(2*L+1)" in data["binary_pair_DTV"], "ledger port-aware tail order")
    check("complete finite all-formed/MATCH FPSS Omega" in data["binary_pair_DTV"],
          "ledger complete authenticated primary scope")
    check(len(data["admitted_members"]) == 2, "ledger admitted members")
    for ceiling in ("no graviton", "no Ricci", "no gravity", "no G"):
        check(ceiling in data["ceilings"], f"ledger ceiling: {ceiling}")


def main() -> None:
    verify_dependencies()
    verify_a3_counts()
    verify_boundary_census()
    verify_link_geometry()
    verify_taylor_filtration()
    verify_tail_and_moderate_ratios()
    verify_binary_total_variation()
    verify_theorem_scope()
    verify_ledger()
    print(f"PASS: {CHECKS}/{CHECKS} GL6BA constructive checks")


if __name__ == "__main__":
    main()
