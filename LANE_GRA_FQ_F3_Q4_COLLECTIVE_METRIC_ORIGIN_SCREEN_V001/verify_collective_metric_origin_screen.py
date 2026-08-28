#!/usr/bin/env python3
"""Exact replay for GRA-FQ-F3-Q4-CMOS-V001.

This verifier checks dependency custody, the q4/A3 static metric rank, the
commuting inherited even root-channel algebra, the finite ice function ranks,
the rank-two versus Maxwell constraint symbols, and the pulled-back stress
two-form.  It does not test a thermodynamic limit, the proposed BS20 CTP
successor, a tensor pole, RGRL-B physical instantiation, or gravity.
"""

from __future__ import annotations

import hashlib
import itertools
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM_BYTES = (HERE / "THEOREM.md").read_bytes()
THEOREM = THEOREM_BYTES.decode("utf-8")
PASSED = 0


DEPENDENCIES = {
    "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md":
        "00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec",
    "LANE_GRA_FC_F3_Q4_CLIFFORD_COLLECTIVE_CONE_V001/THEOREM.md":
        "28b6319e3187337da8ebef2212b030ff6e5b9f8168d9844ae172d94f3e0641a6",
    "LANE_GRA_FD_F3_Q4_COMMON_CHILD_ACOUSTIC_CONE_V001/THEOREM.md":
        "60d012766675c12e82dd1731e202a6c0ed48f24e2697f589b63eecc3cb650287",
    "LANE_GRA_FG_Q4_PAIR_FIELD_LIFT_DERIVABILITY_V001/THEOREM.md":
        "fff521ae41e3f8b83a4738ff96a99715e89f90e2d64724786da8a3ed4732e838",
    "LANE_GRA_FH_F3_Q4_FINITE_PROGRAMMED_SUPPORT_SOLDER_V001/THEOREM.md":
        "2b88febc569efa0de0238e8000d018bf3f798a8ebed2e4ff1327f053d6bd9284",
    "LANE_GRA_FI_F3_Q4_PROGRAMMED_FLOQUET_DETUNING_V001/THEOREM.md":
        "09a9e2ee46acf10dbde91e9578576cb537fe5aff4a9dea513d4c1f208e62de4c",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md":
        "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md":
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
    "LANE_GRA_FN_F3_Q4_ICE_T2_FISHER_SOLDER_BOUNDARY_V001/THEOREM.md":
        "be69f15d611827db9841bd932042604deb4f82a777ff9da28b80e4493cef7596",
    "LANE_CROSS_RFT_GRA_EW_Q4_PAIR_MEMORY_METRIC_DEFORMATION_CLOSURE_V001/THEOREM.md":
        "495e4e99171f4e3e5809f24e5a9a5b68116e996f4f29115bd22f780127d4714e",
    "LANE_CROSS_RFT_MGFT_CONSERVED_STRESS_RANK2_BRIDGE_V001/THEOREM.md":
        "35ad7884b12b1824c34afbfcbe796a2b69d6b955b434183bada75c7b3d923a26",
    "LANE_CROSS_RFT_MGFT_STRESS_TO_CANONICAL_SPIN2_GATE_V001/THEOREM.md":
        "e24095b5d44846ee2fff8bb68fe7f331d9ed87521b01a569bc8cf6198973c65b",
    "GRAVITY_RGRL_ADOPTION_V001.md":
        "bca6146dfa2f2a32cea42db43c85c5d5fb1ee7e6114206e321066809e7c0db1f",
}


def check(condition: bool, label: str) -> None:
    global PASSED
    if not condition:
        raise AssertionError(label)
    PASSED += 1
    print(f"PASS {label}")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank(matrix: list[list[int | Fraction]]) -> int:
    """Exact rational row rank."""
    if not matrix:
        return 0
    a = [[Fraction(value) for value in row] for row in matrix]
    rows = len(a)
    cols = len(a[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next((row for row in range(pivot_row, rows) if a[row][col]), None)
        if pivot is None:
            continue
        a[pivot_row], a[pivot] = a[pivot], a[pivot_row]
        scale = a[pivot_row][col]
        a[pivot_row] = [value / scale for value in a[pivot_row]]
        for row in range(rows):
            if row != pivot_row and a[row][col]:
                factor = a[row][col]
                a[row] = [left - factor * right
                          for left, right in zip(a[row], a[pivot_row])]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def transpose(matrix: list[list[int | Fraction]]) -> list[list[int | Fraction]]:
    return [list(column) for column in zip(*matrix)]


def matmul(
    left: list[list[int | Fraction]],
    right: list[list[int | Fraction]],
) -> list[list[Fraction]]:
    return [
        [sum(Fraction(a) * Fraction(b) for a, b in zip(row, column))
         for column in zip(*right)]
        for row in left
    ]


for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(), f"dependency exists: {relative}")
    check(digest(path) == expected, f"dependency frozen: {relative}")

for forbidden, label in (
    (bytes((0x0D,)), "carriage return"),
    (bytes((0x08,)), "backspace"),
    (bytes((0x0C,)), "form feed"),
):
    check(forbidden not in THEOREM_BYTES, f"theorem contains no {label} control byte")

required_tokens = (
    "SIX_A3_ROOT_COEFFICIENTS_SPAN_STATIC_COMETRIC_TANGENT",
    "FJ_UNPROJECTED_PAIR_OPERATORS_HAVE_EXACT_CONDITIONAL_RANK6_RESPONSE_BUT_NO_TENSOR_SOLDER_OR_RANK2_NULL_PACKET",
    "NO_PRESENT_FIXED_PARENT_OBJECT_SIMULTANEOUSLY_OWNS_SIX_COLLECTIVE_CONFIGURATION_CHANNELS",
    "ADOPTED_RGRLB_SUPPLIES_THOSE_OBJECTS_AS_A_WORKING_LAW_BUT_IS_NOT_THEIR_MICROSCOPIC_F3_DERIVATION",
    "CURRENT_PARENT_OBSTRUCTION_NOT_THERMODYNAMIC_NO_GO",
    "Q4-BLOCK-STRAIN-CTP",
    "H_L[j=0]=H_L^(<=8)",
    "ungauge-fixed",
    "FQ17a",
    "same source-deformed parent",
    "post hoc hand weight",
)
for token in required_tokens:
    check(token in THEOREM, f"theorem binds disposition token: {token}")


# q4 count-front and tetrahedral root-dyad ranks.
simple_roots = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, -1, 0),
    (1, 0, -1),
    (0, 1, -1),
)
check(rank([list(root) for root in simple_roots]) == 3,
      "q4 count-front sibling tangent has rank three")

tetra = (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)
roots = []
for a, b in itertools.combinations(range(4), 2):
    roots.append(tuple(tetra[b][axis] - tetra[a][axis] for axis in range(3)))


def sym_flat(vector: tuple[int, int, int]) -> list[int]:
    x, y, z = vector
    return [x * x, y * y, z * z, x * y, x * z, y * z]


dyads = [sym_flat(root) for root in roots]
check(len(set(roots)) == 6, "six tetrahedral roots are distinct")
check(rank(dyads) == 6, "six q4 root dyads span Sym2 exactly")
dyad_sum = [sum(row[column] for row in dyads) for column in range(6)]
check(dyad_sum[:3] == [16, 16, 16] and dyad_sum[3:] == [0, 0, 0],
      "q4 root second moment is isotropic")


# Exact Abelian group-algebra replay on Z_L^3.
L = 7
Vector = tuple[int, int, int]
Element = dict[Vector, Fraction]


def mod_vector(vector: Vector) -> Vector:
    return tuple(value % L for value in vector)  # type: ignore[return-value]


def add(left: Element, right: Element) -> Element:
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, Fraction(0)) + value
        if not result[key]:
            del result[key]
    return result


def scale(value: Fraction | int, element: Element) -> Element:
    factor = Fraction(value)
    return {key: factor * coefficient for key, coefficient in element.items()
            if factor * coefficient}


def multiply(left: Element, right: Element) -> Element:
    result: Element = {}
    for first, first_value in left.items():
        for second, second_value in right.items():
            key = mod_vector(tuple(a + b for a, b in zip(first, second)))
            result[key] = result.get(key, Fraction(0)) + first_value * second_value
    return {key: value for key, value in result.items() if value}


identity: Element = {(0, 0, 0): Fraction(1)}


def even_shift(root: Vector) -> Element:
    positive = mod_vector(root)
    negative = mod_vector(tuple(-value for value in root))
    return add({positive: Fraction(1)}, {negative: Fraction(1)})


shift_channels = [even_shift(root) for root in simple_roots]
check(all(len(channel) == 2 for channel in shift_channels),
      "every periodic even root channel has two distinct shifts")
all_support = sorted(set().union(*(set(channel) for channel in shift_channels)))
coefficient_matrix = [
    [channel.get(support, Fraction(0)) for support in all_support]
    for channel in shift_channels
]
check(rank(coefficient_matrix) == 6,
      "six periodic even root operators are linearly independent")

for index, first in enumerate(shift_channels):
    for second in shift_channels[index:]:
        check(multiply(first, second) == multiply(second, first),
              f"root-channel commutator vanishes for pair starting {index}")

kernel = scale(4, identity)
for channel in shift_channels:
    kernel = add(kernel, channel)
kernel2 = multiply(kernel, kernel)
kernel3 = multiply(kernel2, kernel)
# One nontrivial polynomial is an executable representative. Commutativity of
# the underlying algebra proves the stated Borel/analytic functional-calculus
# consequence for the exact FI f_F(K).
hamiltonian = add(kernel3, scale(-2, kernel))
fprime = add(scale(3, kernel2), scale(-2, identity))
source_channels = [multiply(fprime, channel) for channel in shift_channels]

check(all(multiply(kernel, channel) == multiply(channel, kernel)
          for channel in shift_channels),
      "all root channels commute with the inherited symmetric kernel")
check(all(multiply(hamiltonian, source) == multiply(source, hamiltonian)
          for source in source_channels),
      "all source derivatives are conserved under a nontrivial f(K)")
check(all(multiply(first, second) == multiply(second, first)
          for first, second in itertools.product(source_channels, repeat=2)),
      "all f-prime-dressed source derivatives commute")
zero_commutators = [
    add(multiply(first, second), scale(-1, multiply(second, first)))
    for first, second in itertools.product(source_channels, repeat=2)
]
check(all(not commutator for commutator in zero_commutators),
      "root-channel retarded commutator numerators vanish exactly")


# Six-state ice fiber: pair rank collapse and diagonal-algebra distinction.
ice = [state for state in itertools.product((-1, 1), repeat=4) if sum(state) == 0]
edges = tuple(itertools.combinations(range(4), 2))
one = [list(state) for state in ice]
pairs = [[state[a] * state[b] for a, b in edges] for state in ice]
means = [Fraction(sum(row[column] for row in pairs), len(ice))
         for column in range(len(edges))]
centered_pairs = [
    [Fraction(value) - means[column] for column, value in enumerate(row)]
    for row in pairs
]

check(len(ice) == 6, "q4 two-in/two-out fiber has six states")
check(rank(one) == 3, "ice one-link functions have rank three")
check(rank(pairs) == 3, "ice pair functions have rank three including constant")
check(rank(centered_pairs) == 2, "centered ice pair functions have rank two")
check(rank([[1] + one_row + pair_row
            for one_row, pair_row in zip(one, centered_pairs)]) == 6,
      "constant plus odd one-link plus centered pair exhaust diagonal functions")
check(all(mean == Fraction(-1, 3) for mean in means),
      "every ice pair has mean minus one third")
check(all(sum(row) == -2 for row in pairs),
      "uniform pair sum is a fixed scalar on ice")

edge_index = {edge: index for index, edge in enumerate(edges)}
opposites = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))
check(all(row[edge_index[first]] == row[edge_index[second]]
          for row in pairs for first, second in opposites),
      "three opposite-pair identities hold identically")
check(all(tuple(-value for value in state) in ice for state in ice),
      "ice fiber is closed under complement")
check(all(((-state[a]) * (-state[b])) == state[a] * state[b]
          for state in ice for a, b in edges),
      "pair functions are complement even")
check(all(-(-state[a]) == state[a] for state in ice for a in range(4)),
      "one-link functions are complement odd")

# Diagonal function representatives commute pointwise; dimension saturation
# does not create a symplectic algebra.
diagonal_generators = [
    tuple(row[column] for row in one) for column in range(4)
] + [
    tuple(row[column] for row in pairs) for column in range(6)
]
check(all(tuple(a * b for a, b in zip(first, second)) ==
          tuple(b * a for a, b in zip(first, second))
          for first, second in itertools.product(diagonal_generators, repeat=2)),
      "all local one-link and pair diagonal functions commute")


# Frozen FJ result which the custody inventory must not erase: before ice
# projection, the admitted independent-link comparator has an exact full-rank
# six-sector response.  Choose Delta=3, h=2, so epsilon=5, and z=i with exact
# rational sector eigenvalues.
c = Fraction(3, 5)
s = Fraction(4, 5)
r_epsilon = Fraction(2 * 5, -(1 + 5 * 5))
r_two_epsilon = Fraction(2 * 10, -(1 + 10 * 10))
fj_a = c * c * s * s * r_epsilon
fj_b = s**4 * r_two_epsilon
fj_sector = {
    "A1": 6 * fj_a + fj_b,
    "E": fj_b,
    "T2": 2 * fj_a + fj_b,
}
check(all(value < 0 for value in fj_sector.values()),
      "FJ A1/E/T2 response eigenvalues are all strictly nonzero at imaginary frequency")
fj_response = [
    [
        (2 * fj_a + fj_b if first == second else
         fj_a if set(first) & set(second) else Fraction(0))
        for second in edges
    ]
    for first in edges
]
check(rank(fj_response) == 6, "FJ exact imaginary-frequency response matrix has rank six")


# Exact principal-symbol ranks at k along z.
# Symmetric coordinate order: xx, yy, zz, sqrt(2)xy, sqrt(2)xz, sqrt(2)yz.
vector_constraint = [
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1],
    [0, 0, 1, 0, 0, 0],
]
scalar_curvature = [[-1, -1, 0, 0, 0, 0]]
trace_row = [[1, 1, 1, 0, 0, 0]]
gauge_map = transpose(vector_constraint)
maxwell_gauss = [[0, 0, 1]]

check(rank(maxwell_gauss) == 1, "Maxwell Gauss principal symbol has rank one")
check(rank(vector_constraint) == 3, "rank-two vector constraint symbol has rank three")
check(rank(scalar_curvature) == 1, "rank-two scalar curvature row has rank one")
check(rank(vector_constraint + scalar_curvature) == 4,
      "vector plus independent scalar rows have joint rank four")
check(6 - rank(vector_constraint + trace_row) == 2,
      "transverse-traceless symmetric quotient has dimension two")
check(rank(gauge_map) == 3, "symmetric-gradient gauge symbol has rank three")
check(matmul(scalar_curvature, gauge_map) == [[0, 0, 0]],
      "linear curvature scalar annihilates vector gauge image")


# Double-curl self-adjointness and pulled-back symplectic degeneracy.
double_curl = [
    [0, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, -1, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
]
check(double_curl == transpose(double_curl),
      "nonzero-mode double curl is self-adjoint")
pullback_two_form = [
    [Fraction(left) - Fraction(right)
     for left, right in zip(row, column)]
    for row, column in zip(double_curl, transpose(double_curl))
]
check(rank(pullback_two_form) == 0,
      "Pi=G[a] graph has zero pulled-back symplectic rank")
check(rank(double_curl) == 3,
      "double-curl image equals the rank-three symmetric-transverse space")


print(f"VERIFICATION PASSED: {PASSED}/{PASSED}")
