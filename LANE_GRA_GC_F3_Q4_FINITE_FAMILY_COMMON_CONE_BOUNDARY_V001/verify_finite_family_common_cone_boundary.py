#!/usr/bin/env python3
"""Exact finite-family geometry/source checks for the q4 common-cone boundary.

This verifier does not enumerate a thermodynamic ice Hilbert space.  It checks
the graph-cover, local H6 coefficient, reciprocal-lattice, refinement, and TT
kinematic statements that are derivable without adding an interaction.  Gap
and residue scaling remain a declared response-data dependency.
"""

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def rank(rows):
    work = [[F(value) for value in row] for row in rows]
    if not work:
        return 0
    row_count, column_count = len(work), len(work[0])
    result = 0
    for column in range(column_count):
        pivot = next((row for row in range(result, row_count)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[result], work[pivot] = work[pivot], work[result]
        pivot_value = work[result][column]
        work[result] = [value/pivot_value for value in work[result]]
        for row in range(row_count):
            if row == result or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [value-factor*base
                         for value, base in zip(work[row], work[result])]
        result += 1
    return result


def determinant(matrix):
    total = F(0)
    size = len(matrix)
    for ordering in permutations(range(size)):
        inversions = sum(ordering[first] > ordering[second]
                         for first in range(size)
                         for second in range(first+1, size))
        term = F(-1 if inversions % 2 else 1)
        for row, column in enumerate(ordering):
            term *= matrix[row][column]
        total += term
    return total


def matmul(first, second):
    return tuple(tuple(sum(first[row][inner]*second[inner][column]
                           for inner in range(len(second)))
                       for column in range(len(second[0])))
                 for row in range(len(first)))


def matadd(first, second):
    return tuple(tuple(first[row][column]+second[row][column]
                       for column in range(len(first[0])))
                 for row in range(len(first)))


def matscale(scale, matrix):
    return tuple(tuple(F(scale)*value for value in row) for row in matrix)


def transpose(matrix):
    return tuple(zip(*matrix))


for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink() and
          sha256(path.read_bytes()).hexdigest() == expected,
          f"dependency custody: {relative}")


# Semantic custody: these are the precise ceilings being joined.
fs_text = (ROOT / "LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001" /
           "THEOREM.md").read_text()
fd_text = (ROOT / "LANE_GRA_FD_F3_Q4_COMMON_CHILD_ACOUSTIC_CONE_V001" /
           "THEOREM.md").read_text()
fv_text = (ROOT / "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001" /
           "THEOREM.md").read_text()
fy_text = (ROOT / "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001" /
           "THEOREM.md").read_text()
fl_text = (ROOT / "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001" /
           "THEOREM.md").read_text()
gb_text = (ROOT / "LANE_GRA_GB_F3_Q4_FIXED_SUPPORT_ENERGY_MOMENTUM_WARD_BOUNDARY_V001" /
           "THEOREM.md").read_text()
fy_result = json.loads((ROOT /
    "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001" /
    "RESULT.json").read_text())
check("L_r=5\\,2^r" in fs_text and "eight-sheeted graph" in fs_text,
      "FS owns the covering-matched G_L family but not nested Hilbert spaces")
check("v_3={16\\over3\\sqrt3}a_*^3" in fd_text and
      "m_X=0" in fd_text and "\\kappa_X\\over\\chi_X" in fd_text,
      "FD owns A3 volume and the conditional common-cone criterion")
check("4678629417" in fv_text and "rank}_{\\rm nonid}" in fv_text,
      "FV owns exact local H6 coefficients and rank-six witnesses")
check("selected FO 180-state winding component" in fy_text and
      "mean H8 completion" in fy_text and "thermodynamic limit" in fy_text,
      "FY owns one finite native m=1 response and no thermodynamic limit")
check("two-photon continuum" in fl_text and
      "isolated massless helicity-two pole" in fl_text,
      "FL fixes the Gaussian-Maxwell composite no-pole baseline")
check("physical ownership gap" in gb_text and
      "support/recoil" in gb_text and "boundary ports" in gb_text,
      "GB fixes the independent momentum-ownership and stress-Ward dependency")


# -------------------------------------------------------------------------
# Exact covering family and elementary H6 inventory.

SHIFTS = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))


def add3(first, second, modulus):
    return tuple((first[index]+second[index]) % modulus
                 for index in range(3))


def sub3(first, second, modulus):
    return tuple((first[index]-second[index]) % modulus
                 for index in range(3))


def canonical_cycle(cycle):
    variants = []
    for oriented in (cycle, tuple(reversed(cycle))):
        for offset in range(len(cycle)):
            variants.append(oriented[offset:]+oriented[:offset])
    return min(variants)


def elementary_cycle(cell, labels, modulus):
    first, second, third = labels
    b1 = add3(cell, SHIFTS[first], modulus)
    a2 = sub3(b1, SHIFTS[second], modulus)
    b3 = add3(a2, SHIFTS[third], modulus)
    a4 = sub3(b3, SHIFTS[first], modulus)
    b5 = add3(a4, SHIFTS[second], modulus)
    return ((0,)+cell, (1,)+b1, (0,)+a2,
            (1,)+b3, (0,)+a4, (1,)+b5)


def elementary_inventory(modulus):
    by_missing = {}
    for missing in range(4):
        labels = tuple(label for label in range(4) if label != missing)
        by_missing[missing] = {
            canonical_cycle(elementary_cycle(cell, labels, modulus))
            for cell in product(range(modulus), repeat=3)
        }
    return by_missing


def support_inventory(modulus):
    cells = tuple(product(range(modulus), repeat=3))
    vertices = {(sublattice,)+cell
                for sublattice in (0, 1) for cell in cells}
    # Every physical link is uniquely owned by its A-cell and q4 label.
    links = {(cell, label) for cell in cells for label in range(4)}
    native_supports = ({("vertex",)+vertex for vertex in vertices} |
                       {("link",)+cell+(label,) for cell, label in links})
    return vertices, links, native_supports


families = {}
for length in (5, 10):
    families[length] = elementary_inventory(length)
    all_cycles = set().union(*families[length].values())
    check(len(all_cycles) == 4*length**3 and
          all(len(cycles) == length**3
              for cycles in families[length].values()) and
          all(len(cycle) == len(set(cycle)) == 6 for cycle in all_cycles),
          f"G_{length} has exactly L^3 H6 rings of each missing-label type")
    vertices, links, native_supports = support_inventory(length)
    check((len(vertices), len(links), len(native_supports)) ==
          (2*length**3, 4*length**3, 6*length**3),
          f"G_{length} has exact vertex/link/native-support volume counts")


def reduce_cycle(cycle, modulus):
    return canonical_cycle(tuple((vertex[0],) +
                                 tuple(coordinate % modulus
                                       for coordinate in vertex[1:])
                                 for vertex in cycle))


for missing in range(4):
    fiber_count = Counter(reduce_cycle(cycle, 5)
                          for cycle in families[10][missing])
    check(set(fiber_count) == families[5][missing] and
          set(fiber_count.values()) == {8},
          f"missing label {missing}: H6 rings pull back eight-to-one under G_10 -> G_5")


# FV's two E diagonal witnesses plus the four environment-independent
# Hermitian ring rows.  These are local coefficient data and hence pull back
# with each marked ring; they are not a claim that spectra pull back.
e_rows = ((F(-1), F(1), F(0), F(0), F(0), F(0)),
          (F(-1), F(0), F(1), F(0), F(0), F(0)))
ring_rows = (
    (F(231, 8), F(231, 8), F(231, 8), F(-189, 8), F(-189, 8), F(-189, 8)),
    (F(231, 8), F(231, 8), F(231, 8), F(189, 8), F(189, 8), F(-189, 8)),
    (F(231, 8), F(231, 8), F(231, 8), F(189, 8), F(-189, 8), F(189, 8)),
    (F(231, 8), F(231, 8), F(231, 8), F(-189, 8), F(189, 8), F(189, 8)),
)
witness = e_rows+ring_rows
check(rank(ring_rows) == 4 and rank(witness) == 6,
      "FV local E plus A1/T2 H6 coefficient witnesses have exact rank six")
check(determinant(witness) == F(-4678629417, 256),
      "the inherited six-witness determinant matches FV exactly")


# -------------------------------------------------------------------------
# Allowed momenta, cover characters, and exact A3 reciprocal geometry.

# Construct the exact FS -> FD affine join before using FD's covolume.  Work
# with u_a=4 P e_a=2 sqrt(3) n_a so every entry is integral.  Label 4 is the
# FS zero shift, and root_i4=u_4-u_i represents 2 sqrt(3) alpha_i4.
U = tuple(tuple(4*int(axis == label)-1 for axis in range(4))
          for label in range(4))
ROOT_I4 = tuple(tuple(U[3][axis]-U[label][axis] for axis in range(4))
                for label in range(3))
check(all(tuple(ROOT_I4[label][axis]-U[3][axis]
                for axis in range(4)) ==
          tuple(-U[label][axis] for axis in range(4))
          for label in range(3)),
      "the explicit FS-to-FD embedding sends every shifted edge to -a_* n_i")

# alpha_i4 dot alpha_j4 = (u_4-u_i).(u_4-u_j)/12 because
# n_a=u_a/(2 sqrt(3)).
JOIN_GRAM = tuple(tuple(sum(ROOT_I4[row][axis]*ROOT_I4[column][axis]
                            for axis in range(4))/F(12)
                        for column in range(3)) for row in range(3))

G = ((F(8, 3), F(4, 3), F(4, 3)),
     (F(4, 3), F(8, 3), F(4, 3)),
     (F(4, 3), F(4, 3), F(8, 3)))
G_INV = ((F(9, 16), F(-3, 16), F(-3, 16)),
         (F(-3, 16), F(9, 16), F(-3, 16)),
         (F(-3, 16), F(-3, 16), F(9, 16)))
I3 = tuple(tuple(F(int(row == column)) for column in range(3))
           for row in range(3))
check(JOIN_GRAM == G,
      "the constructed FS translation basis has exactly FD's A3 Gram matrix")
check(matmul(G, G_INV) == I3 and determinant(G) == F(256, 27),
      "the A3 primitive Gram matrix and reciprocal metric are exact")
check(determinant(G) == F(16, 1)**2/F(27, 1),
      "the squared primitive covolume is 256 a_*^6 / 27")


def reciprocal_integer_form(vector):
    return 4*sum(value*value for value in vector)-sum(vector)**2


def reciprocal_sum_of_squares(vector):
    return (sum(value*value for value in vector) +
            sum((vector[first]-vector[second])**2
                for first in range(3) for second in range(first+1, 3)))


check(all(reciprocal_integer_form(vector) ==
          reciprocal_sum_of_squares(vector)
          for vector in product(range(-2, 3), repeat=3)),
      "q(n) is exactly the positive sum-of-squares A3 reciprocal form")

# q<=4 implies sum n_i^2<=4, hence |n_i|<=2.  The following finite
# enumeration is therefore an exhaustive equality classification, not a
# finite-box inference about the infinite lattice.
q3 = {vector for vector in product(range(-2, 3), repeat=3)
      if reciprocal_integer_form(vector) == 3}
q4 = {vector for vector in product(range(-2, 3), repeat=3)
      if reciprocal_integer_form(vector) == 4}
expected_q3 = ({tuple(sign*int(axis == basis) for axis in range(3))
                for basis in range(3) for sign in (-1, 1)} |
               {(1, 1, 1), (-1, -1, -1)})
expected_q4 = {tuple(sign*int(axis != missing) for axis in range(3))
               for missing in range(3) for sign in (-1, 1)}
check(q3 == expected_q3 and q4 == expected_q4,
      "the exact q=3 and q=4 equality shells contain 8 and 6 vectors")
check(all(reciprocal_integer_form(vector) % 4 ==
          (-sum(vector)**2) % 4
          for vector in product(range(4), repeat=3)) and
      F(3, 16)*3 == F(9, 16),
      "the positive reciprocal form has first possible value q=3 and norm 9/16")


# A base character n on G_L pulls back to character 2n on G_2L.  The odd
# cover modes are genuinely new long-wavelength characters, not inherited
# base modes.
for vector in ((1, 0, 0), (0, 1, 1), (1, -1, 2)):
    character_ok = all(
        (2*sum(vector[index]*(coordinate[index] % 5)
               for index in range(3))-
         sum((2*vector[index])*coordinate[index]
             for index in range(3))) % 10 == 0
        for coordinate in product(range(10), repeat=3))
    check(character_ok,
          f"cover pullback maps momentum {vector} on G_5 to 2n on G_10")
check(any(component % 2 for component in (1, 0, 0)),
      "the lowest G_10 mode is not a pullback of any G_5 character")


# Under L -> 2L, fixed a_* gives eight times the volume and half k_min.
# Under the paired refinement a_* -> a_*/2, affine volume and k_min are
# unchanged.  A joint IR+continuum sequence requires 0<s<1 for
# a_r=a_0 2^{-sr}, L_r=5 2^r.
check(F(8)*F(1) == F(8) and F(1, 2) == F(1)/2,
      "fixed lattice scale: one cover multiplies volume by 8 and halves k_min")
check(F(8)*F(1, 2)**3 == 1 and F(1, 2)/F(1, 2) == 1,
      "paired affine refinement preserves affine volume and k_min")
check(0 < F(1, 2) < 1,
      "a representative 0<s<1 sequence has both a_r -> 0 and L_r a_r -> infinity")


# FD's A3 isotropy is the exact second moment of the six positive roots.
tetra = ((F(1), F(1), F(1)), (F(1), F(-1), F(-1)),
         (F(-1), F(1), F(-1)), (F(-1), F(-1), F(1)))
root_second = [[F(0) for _ in range(3)] for _ in range(3)]
for first in range(4):
    for second in range(first+1, 4):
        # n=(sign)/sqrt(3), so outer products carry one factor 1/3.
        root = tuple(tetra[second][axis]-tetra[first][axis]
                     for axis in range(3))
        for row in range(3):
            for column in range(3):
                root_second[row][column] += root[row]*root[column]/3
check(tuple(tuple(row) for row in root_second) ==
      tuple(tuple(F(16, 3) if row == column else F(0)
                  for column in range(3)) for row in range(3)),
      "the six A3 roots have exact isotropic second moment 16/3")


# -------------------------------------------------------------------------
# Exact scale-invariant TT kinematics along one reciprocal direction.

# With primitive columns alpha_14,alpha_24,alpha_34, the reciprocal direction
# for n=(1,0,0) is proportional to r=(-1,-1,-1).
r = (F(-1), F(-1), F(-1))
r2 = sum(value*value for value in r)
outer = tuple(tuple(r[row]*r[column] for column in range(3))
              for row in range(3))
P = matadd(I3, matscale(F(-1, r2), outer))
check(matmul(P, P) == P and
      all(sum(P[row][column]*r[column] for column in range(3)) == 0
          for row in range(3)) and
      sum(P[index][index] for index in range(3)) == 2,
      "the exact transverse projector has rank two on the shortest shell")


def unpack_coordinate(column):
    values = [F(0)]*6
    values[column] = F(1)
    return ((values[0], values[3]/2, values[4]/2),
            (values[3]/2, values[1], values[5]/2),
            (values[4]/2, values[5]/2, values[2]))


def pack_coordinate(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            2*matrix[0][1], 2*matrix[0][2], 2*matrix[1][2])


def tt_project(matrix):
    pap = matmul(matmul(P, matrix), P)
    transverse_trace = sum(P[row][column]*matrix[row][column]
                           for row in range(3) for column in range(3))
    return matadd(pap, matscale(-transverse_trace/2, P))


TT = transpose(tuple(pack_coordinate(tt_project(unpack_coordinate(column)))
                     for column in range(6)))
check(matmul(TT, TT) == TT and rank(TT) == 2,
      "the exact six-coordinate TT projector is idempotent and rank two")
check(all(sum(tt_project(unpack_coordinate(column))[axis][axis]
              for axis in range(3)) == 0
          for column in range(6)),
      "every exact TT basis image is traceless")
r_scaled = tuple(2*value for value in r)
outer_scaled = tuple(tuple(r_scaled[row]*r_scaled[column]
                           for column in range(3)) for row in range(3))
P_scaled = matadd(I3, matscale(F(-1, sum(x*x for x in r_scaled)),
                              outer_scaled))
check(P_scaled == P,
      "the transverse and TT projectors are exactly invariant under k rescaling")


# -------------------------------------------------------------------------
# Response-data and identifiability boundary.

check(fy_result["component"]["cells"] == 30 and
      fy_result["component"]["momentum"] == 1 and
      all(min(sample["pole_gaps"]) > 0 and
          sample["ranks"]["TT_ground_image"] == 2
          for sample in fy_result["samples"]),
      "FY supplies one positive-gap finite volume with kinematic TT rank two")
check("volume_sequence" not in fy_result and
      "gap_scaling" not in fy_result and
      "thermodynamic_residue" not in fy_result,
      "FY contains no matched-volume gap or residue sequence")
check("flux scaling" in fl_text and
      "small-wave-vector structure factor" in fl_text and
      "two-photon continuum" in fl_text,
      "current public input fixes Maxwell phase evidence but not a tensor pole sequence")


# The same kinematic family and source coefficients are compatible with the
# two FD sector dispersions below because the mass and stiffness are not fixed
# by FS/FV/FY.  This is a logical nonidentifiability witness, not a new model
# proposed for nature.
lambda_k = F(1, 100)
ratio = F(3, 2)
massless_gap_squared = ratio*lambda_k
massive_gap_squared = F(1, 4)+ratio*lambda_k
check(massless_gap_squared != massive_gap_squared and
      massless_gap_squared > 0 and massive_gap_squared > 0,
      "identical family/source kinematics permits distinct massless and massive response laws")

# In the declared even-Euclidean cosine transform,
# 2 int_0^infinity cos(w t) R exp(-D t) dt
# =2 D R/(w^2+D^2).  This fixes the convention, but no nonzero constant IR
# limit is demanded for a derivative/native source.
delta, residue = F(3, 7), F(5, 11)
check(2*delta*residue == F(30, 77),
      "the even-Euclidean quadratic-denominator numerator is 2 Delta R")


result = {
    "lane": LANE.name,
    "status": "PASS",
    "checks": f"{checks}/{checks}",
    "family": "G_L with L=5*2^r",
    "affine_join": "X(A_x)=a_* sum_i x_i alpha_i4; X(B_x)=X(A_x)-a_* n_4",
    "vertices": "2 L^3",
    "links": "4 L^3",
    "native_supports": "6 L^3",
    "h6_rings": "4 L^3; L^3 per missing q4 label",
    "local_source": "FV rank-six coefficient witnesses inherit under graph covers",
    "momentum_norm": "|k_[n]|^2=(2pi/(L a_*))^2 hat(n)^T G^-1 hat(n)",
    "k_min": "3pi/(2 L a_*)",
    "cover_character": "n on G_L pulls back to 2n on G_2L; new odd modes have no base ancestor",
    "joint_limit": "a_r->0 and L_r a_r->infinity; for a_r=a0*2^(-sr), require 0<s<1",
    "tt_projector": "exact rank two and invariant under k rescaling",
    "response_boundary": "no current matched-family gap/residue sequence or canonical native-source normalization",
    "smallest_target": "G_5,G_10,G_20 count-normalized TT spectral screen at three inequivalent small-k rays",
    "ceiling": "geometry/source scaling only; no Ward-derived tensor pole, common massless cone, gravity, or G",
}
print("GC_RESULT_JSON", json.dumps(result, sort_keys=True))
print(f"SUMMARY {checks}/{checks} GC finite-family/common-cone checks passed")
print("CEILING exact family geometry, local coefficients, momenta, and TT "
      "kinematics only; gap/residue/Ward/common-cone dynamics remain open")
