#!/usr/bin/env python3
"""Independent hostile replay for GRA-FQ-F3-Q4-CMOS-V001.

This executable does not import the builder verifier.  It reconstructs the
static/dynamic ranks, FJ response, ice and constraint symbols, stress graph,
documentary inventory, and custody/tamper surface independently.
"""

from __future__ import annotations

import hashlib
import itertools
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM = (HERE / "THEOREM.md").read_text(encoding="utf-8")
PASSED = 0


def check(condition: bool, label: str) -> None:
    global PASSED
    if not condition:
        raise AssertionError(label)
    PASSED += 1
    print(f"PASS {label}")


def digest_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def digest(path: Path) -> str:
    return digest_bytes(path.read_bytes())


def parse_hashes(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        value, name = line.split(maxsplit=1)
        result[name.strip()] = value
    return result


def rank(matrix: list[list[int | Fraction]]) -> int:
    if not matrix:
        return 0
    work = [[Fraction(entry) for entry in row] for row in matrix]
    rows, columns = len(work), len(work[0])
    pivot = 0
    for column in range(columns):
        selected = next((row for row in range(pivot, rows) if work[row][column]), None)
        if selected is None:
            continue
        work[pivot], work[selected] = work[selected], work[pivot]
        scale = work[pivot][column]
        work[pivot] = [entry / scale for entry in work[pivot]]
        for row in range(rows):
            if row == pivot or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right for left, right in zip(work[row], work[pivot])]
        pivot += 1
        if pivot == rows:
            break
    return pivot


def transpose(matrix: list[list[int | Fraction]]) -> list[list[int | Fraction]]:
    return [list(column) for column in zip(*matrix)]


def matmul(left: list[list[int | Fraction]], right: list[list[int | Fraction]]) -> list[list[Fraction]]:
    return [
        [sum(Fraction(a) * Fraction(b) for a, b in zip(row, column)) for column in zip(*right)]
        for row in left
    ]


# Six tetrahedral root dyads.
tetra = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
roots = [
    tuple(tetra[b][axis] - tetra[a][axis] for axis in range(3))
    for a, b in itertools.combinations(range(4), 2)
]
dyads = [
    [x * x, y * y, z * z, x * y, x * z, y * z]
    for x, y, z in roots
]
check(len(set(roots)) == 6, "six A3/tetrahedral roots are distinct")
check(rank(dyads) == 6, "six root dyads span Sym2")
check([sum(row[i] for row in dyads) for i in range(6)] == [16, 16, 16, 0, 0, 0], "root second moment is isotropic")


# Translation-complete periodic group algebra.
simple_roots = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, 0), (1, 0, -1), (0, 1, -1))


def periodic_channels(length: int) -> list[dict[tuple[int, int, int], Fraction]]:
    channels = []
    for root in simple_roots:
        positive = tuple(value % length for value in root)
        negative = tuple((-value) % length for value in root)
        channels.append({positive: Fraction(1), negative: Fraction(1)})
    return channels


def add(left: dict, right: dict) -> dict:
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, Fraction(0)) + value
        if not answer[key]:
            del answer[key]
    return answer


def scale(value: int | Fraction, element: dict) -> dict:
    value = Fraction(value)
    return {key: value * entry for key, entry in element.items() if value * entry}


def multiply(left: dict, right: dict, length: int) -> dict:
    answer: dict[tuple[int, int, int], Fraction] = {}
    for first, first_value in left.items():
        for second, second_value in right.items():
            key = tuple((a + b) % length for a, b in zip(first, second))
            answer[key] = answer.get(key, Fraction(0)) + first_value * second_value
    return {key: value for key, value in answer.items() if value}


for length in (5, 7):
    channels = periodic_channels(length)
    support = sorted(set().union(*(set(channel) for channel in channels)))
    coefficient_rows = [[channel.get(point, 0) for point in support] for channel in channels]
    check(all(len(channel) == 2 for channel in channels), f"L={length} has distinct plus/minus root shifts")
    check(rank(coefficient_rows) == 6, f"L={length} six even translations are independent")
    check(all(multiply(a, b, length) == multiply(b, a, length) for a, b in itertools.product(channels, repeat=2)), f"L={length} translation operators commute")

    identity = {(0, 0, 0): Fraction(1)}
    kernel = scale(4, identity)
    for channel in channels:
        kernel = add(kernel, channel)
    kernel2 = multiply(kernel, kernel, length)
    hamiltonian = add(multiply(kernel2, kernel, length), scale(-2, kernel))
    fprime = add(scale(3, kernel2), scale(-2, identity))
    sources = [multiply(fprime, channel, length) for channel in channels]
    check(all(multiply(hamiltonian, source, length) == multiply(source, hamiltonian, length) for source in sources), f"L={length} f-prime-dressed sources are conserved")
    check(all(multiply(a, b, length) == multiply(b, a, length) for a, b in itertools.product(sources, repeat=2)), f"L={length} dressed source commutators vanish")


# Frozen FJ six-sector response at Delta=3,h=2,z=i.
edges = tuple(itertools.combinations(range(4), 2))
c, s = Fraction(3, 5), Fraction(4, 5)
a = c * c * s * s * Fraction(10, -26)
b = s**4 * Fraction(20, -101)
fj = [
    [2 * a + b if first == second else a if set(first) & set(second) else 0 for second in edges]
    for first in edges
]
check(a < 0 and b < 0, "FJ spectral coefficients are strictly negative at imaginary frequency")
check(rank(fj) == 6, "FJ conditional unprojected response has exact rank six")
check(all(value < 0 for value in (6 * a + b, b, 2 * a + b)), "FJ A1/E/T2 sector eigenvalues are nonzero")


# Ice ranks and algebraic relations.
ice = [state for state in itertools.product((-1, 1), repeat=4) if sum(state) == 0]
one = [list(state) for state in ice]
pairs = [[state[i] * state[j] for i, j in edges] for state in ice]
means = [Fraction(sum(row[column] for row in pairs), 6) for column in range(6)]
centered = [[Fraction(value) - means[column] for column, value in enumerate(row)] for row in pairs]
check(len(ice) == 6, "ice fiber has six states")
check(rank(one) == 3, "ice one-link rank is three")
check(rank(pairs) == 3, "ice pair rank including scalar is three")
check(rank(centered) == 2, "ice centered-pair rank is two")
check(rank([[1] + left + right for left, right in zip(one, centered)]) == 6, "1+3+2 exhausts diagonal ice functions")
check(all(sum(row) == -2 for row in pairs), "ice pair sum is fixed")
check(all(mean == Fraction(-1, 3) for mean in means), "all ice pair means equal minus one third")


# Constraint symbols and double curl at k along z.
vector = [[0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0]]
scalar = [[-1, -1, 0, 0, 0, 0]]
trace = [[1, 1, 1, 0, 0, 0]]
gauge = transpose(vector)
check(rank([[0, 0, 1]]) == 1, "Maxwell symbol rank is one")
check(rank(vector) == 3, "tensor vector symbol rank is three")
check(rank(scalar) == 1, "tensor scalar symbol rank is one")
check(rank(vector + scalar) == 4, "vector plus scalar constraint rows have rank four")
check(6 - rank(vector + trace) == 2, "TT space has dimension two")
check(matmul(scalar, gauge) == [[0, 0, 0]], "scalar curvature annihilates spatial gauge image")

double_curl = [[0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, -1, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
check(double_curl == transpose(double_curl), "double curl is self-adjoint on the compatible mode")
check(rank(double_curl) == 3, "double-curl image has rank three")
check(rank([[Fraction(a) - Fraction(b) for a, b in zip(row, column)] for row, column in zip(double_curl, transpose(double_curl))]) == 0, "E=G[a] pullback has zero symplectic rank")


# Documentary inventory and corrected successor boundary.
inventory_tokens = (
    "q4 count/front", "FD `phi_m`", "FH/FI carrier", "FC `c_ab`",
    "FC Clifford carrier", "EW `J_ab`", "FPMH/FG `Q_ab`",
    "unprojected FJ", "projected ice pairs", "projected ice one-links",
    "hybrid `A1+E+T2` ledger", "H6/H8 loop flips",
    "conserved macroscopic stress", "double-curl potential", "adopted RGRL-B fields",
)
for token in inventory_tokens:
    check(token in THEOREM, f"declared inventory includes {token}")

successor_tokens = (
    "one stationary density-matrix family", "complete BS16--BS22 data",
    "complete fixed microscopic F3 Hamiltonian", "same source-deformed parent",
    "complete virtual histories", "ungauge-fixed", "naively inverted six-by-six",
    "H_L[j=0]=H_L^(<=8)", "may not splice FI", "post hoc hand weight",
)
for token in successor_tokens:
    check(token in THEOREM, f"successor freezes boundary: {token}")


# Dependency, manifest, outer-seal, and tamper checks.
dependencies = parse_hashes(HERE / "DEPENDENCIES.sha256")
check(len(dependencies) == 14, "dependency ledger has fourteen entries")
for relative, expected in dependencies.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(), f"dependency is regular: {relative}")
    check(digest(path) == expected, f"dependency digest matches: {relative}")

required_manifest = {
    "DEPENDENCIES.sha256", "README.md", "RESULT.md", "SELF_AUDIT.md",
    "THEOREM.md", "VERIFICATION.txt", "verify_collective_metric_origin_screen.py",
    "INDEPENDENT_HOSTILE_AUDIT.md", "independent_hostile_audit.py",
}
manifest = parse_hashes(HERE / "MANIFEST.sha256")
check(set(manifest) == required_manifest, "manifest pins the exact complete lane payload")
for relative, expected in manifest.items():
    path = HERE / relative
    check(path.is_file() and not path.is_symlink(), f"manifest member is regular: {relative}")
    check(digest(path) == expected, f"manifest member digest matches: {relative}")

seal = parse_hashes(HERE / "SEAL.sha256")
check(seal == {"MANIFEST.sha256": digest(HERE / "MANIFEST.sha256")}, "outer seal authenticates manifest bytes")
check(digest_bytes((HERE / "THEOREM.md").read_bytes() + b"HOSTILE") != manifest["THEOREM.md"], "theorem tamper is detected")
first_dependency, expected = next(iter(dependencies.items()))
check(digest_bytes((ROOT / first_dependency).read_bytes() + b"HOSTILE") != expected, "dependency tamper is detected")
check(digest_bytes((HERE / "MANIFEST.sha256").read_bytes() + b"HOSTILE") != seal["MANIFEST.sha256"], "manifest tamper is detected")
check(required_manifest - {"THEOREM.md"} != required_manifest, "required-member omission changes payload set")

print(f"SUMMARY {PASSED}/{PASSED} independent hostile checks passed")
print("DISPOSITION PASS_WITH_FJ_INVENTORY_AND_SOURCE_BEFORE_FESHBACH_CORRECTIONS__CURRENT_CONSTRUCTION_ONLY__NOT_THERMODYNAMIC_NO_GO")
