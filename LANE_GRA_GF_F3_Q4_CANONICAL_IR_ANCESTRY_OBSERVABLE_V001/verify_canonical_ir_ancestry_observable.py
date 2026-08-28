#!/usr/bin/env python3
"""Verify the prospective GF/CIAO observable contract.

All finite matrices below are algebraic contract self-tests.  They are not
spectral data from G_L and must not be reported as evidence for a pole.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import product
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


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink() and
          digest(path) == expected,
          f"dependency custody: {relative}")

check(len((LANE / "DEPENDENCIES.sha256").read_text().splitlines()) == 25,
      "twenty-five shortcut/FD/GC/FY/FX/GD dependencies are frozen")


# -------------------------------------------------------------------------
# Machine-readable claim surface.

contract = json.loads((LANE / "OBSERVABLE_CONTRACT.json").read_text())
result = json.loads((LANE / "RESULT.json").read_text())
theorem = " ".join((LANE / "THEOREM.md").read_text().split())
self_audit = " ".join((LANE / "SELF_AUDIT.md").read_text().split())

check(contract["version"] == "V005" and
      contract["status"] ==
      "REPAIRED_V005_AWAITING_INDEPENDENT_AUDIT_PHYSICS_UNEXECUTED" and
      result["status"] ==
      "REPAIRED_AFTER_V004_REJECTION_AWAITING_INDEPENDENT_SEAL_DESIGN_ONLY",
      "machine status marks V005 repaired after the V004 rejection")
check(contract["ceiling"] ==
      "no matched-family pole data, no Poincare/helicity result, no positive G2, no G3 soft theorem, no gravity, and no G",
      "contract retains the no-pole/no-helicity/no-G2 ceiling")
check(result["ceiling"] ==
      "no numerical pole evidence, no helicity result, no positive G2, no G3, no gravity, and no G",
      "result retains the no-numerical-evidence ceiling")
check(contract["source"]["post_result_rescaling"] == "forbidden" and
      "without such D is INDETERMINATE" in result["residue"],
      "post-result pole normalization is forbidden")
check(contract["momenta"]["raw_all_classes_required"] is True and
      contract["momenta"]["conjugates_required"] is True,
      "complete raw momentum classes and conjugates are mandatory")
check(contract["momenta"]["tt_domain"] == "all nonzero characters only" and
      "Pi_TT(0) is undefined and forbidden" in
      contract["momenta"]["zero_character"],
      "TT domain excludes k=0 and types homogeneous data separately")
check(contract["momenta"]["rays"] == [[1, 0, 0], [1, 1, 0], [1, -1, 0]] and
      contract["momenta"]["ray_q"] == [3, 4, 8],
      "the three inequivalent low-momentum rays are frozen")
check("complete exact energy-and-quantum-number projectors" in
      contract["cluster"]["projector"] and
      "no individual eigenvector" in contract["cluster"]["basis_rule"],
      "degenerate spectral data are projector typed")
check("liminf eta_min>0" in contract["ancestry"]["pass"] and
      "generalized eigenvalues" in contract["ancestry"]["normalized_matrix"],
      "normalized ancestry has a basis-independent positive gate")
check("ancestry" not in contract["canonical_source_binding"]["fail"] and
      "reported independently" in
      contract["canonical_source_binding"]["separation_rule"] and
      "separate outputs" in result["decision_separation"],
      "canonical amplitude and normalized ancestry are disjoint classifiers")
check("finite sigma_max" in contract["canonical_source_binding"]["pass"] and
      "upper-bound divergence" in
      contract["canonical_source_binding"]["fail"],
      "canonical amplitude PASS has an executable finite upper-bound complement")
check("rho0=P0/rank(P0)" in contract["ground"]["allowed"] and
      "never load bearing" in contract["ground"]["basis_rule"],
      "ground query is rank one or basis-invariant on a degenerate ground space")
check("independently audited and sealed FD certificate" in
      contract["family"]["physical_binding"] and
      "equal slopes are not" in contract["common_cone"]["ceiling"],
      "physical scale binding and common-cone ceiling are explicit")
check("helicities +2,-2" in contract["representation"]["little_group"] and
      "scalar doublet" in contract["representation"]["fail"],
      "Poincare little-group gate rejects the scalar-doublet counterexample")
check("at least one route" in contract["factorization"]["pass"] and
      "not longitudinal decoupling" in contract["factorization"]["ceiling"],
      "factorization readiness is separated from G3")
check(len(contract["spectral_pass"]) == 8,
      "spectral pass surface includes representation and canonical outcome")
check(len(contract["payload"]["required"]) == 10,
      "the executable G2 payload has ten frozen fields")
check("Changing any of those after examining a spectrum" in theorem and
      "three-size/ray packet is a screen" in self_audit,
      "the prose packet forbids post-result rescue and finite-size promotion")
check("P_{C,L}(k)Q^{TT\\dagger}_{L,\\rm vol}(k)P_{0,L}" in theorem and
      "ground-to-pole map" in theorem,
      "the ancestry identity is typed as a ground-to-pole map")
check("Frobenius inner product" in theorem and
      "arbitrary nonorthogonal" in theorem,
      "canonical square-root invariance is limited to the physical TT metric")


# -------------------------------------------------------------------------
# Exact Q(sqrt(3)) volume/source conversion at a_*=1.

def q3_add(first, second):
    return (first[0] + second[0], first[1] + second[1])


def q3_mul(first, second):
    # a+b sqrt(3)
    return (first[0] * second[0] + 3 * first[1] * second[1],
            first[0] * second[1] + first[1] * second[0])


v3 = (F(0), F(16, 9))          # 16 sqrt(3)/9 = 16/(3 sqrt(3))
volume_factor_squared = (F(0), F(3, 8))  # 2/v3 = 3 sqrt(3)/8
check(q3_mul(v3, volume_factor_squared) == (F(2), F(0)),
      "Q_aff/Q_cnt squared is exactly 2/v3")
source_factor_squared = (F(0), F(8, 9))  # v3/2 = 8 sqrt(3)/9
check(q3_mul(volume_factor_squared, source_factor_squared) == (F(1), F(0)),
      "the conjugate source conversion is exactly inverse")
for length in (5, 10, 20):
    cells = length ** 3
    # Squared coefficient identity:
    # (2/v3)*(1/(2 L^3)) = 1/(L^3 v3).
    left = q3_mul(volume_factor_squared, (F(1, 2 * cells), F(0)))
    inverse_volume = q3_mul((F(0), F(3, 16)),
                            (F(1, cells), F(0)))
    check(left == inverse_volume,
          f"L={length}: affine-volume mode equals V_L^(-1/2) cell sum")


# -------------------------------------------------------------------------
# Exact A3 momentum registry.

def qform(vector):
    return 4 * sum(value * value for value in vector) - sum(vector) ** 2


def canonical_representative(residue, length):
    # Searching shifts [-2,2]^3 is global.  A coordinate outside this shift
    # box has norm >2L, hence q>=|n|^2>4L^2, while a centered candidate has
    # q<=4|n|^2<=3L^2.
    candidates = []
    for shift in product(range(-2, 3), repeat=3):
        vector = tuple(residue[axis] + length * shift[axis]
                       for axis in range(3))
        candidates.append((qform(vector), vector))
    return min(candidates)[1]


registries = {}
for length in (5, 10, 20):
    registry = {}
    class_membership_ok = True
    for residue in product(range(length), repeat=3):
        representative = canonical_representative(residue, length)
        class_membership_ok = (class_membership_ok and
                               all((representative[axis] - residue[axis]) %
                                   length == 0 for axis in range(3)))
        registry[residue] = representative
    registries[length] = registry
    check(class_membership_ok,
          f"L={length} every representative remains in its momentum class")
    check(len(registry) == length ** 3,
          f"L={length} registry contains every character exactly once")
    check(registry[(0, 0, 0)] == (0, 0, 0) and
          all(qform(vector) > 0 for residue, vector in registry.items()
              if residue != (0, 0, 0)),
          f"L={length} only the zero character has zero reciprocal norm")
    conjugates_ok = True
    for residue, representative in registry.items():
        negative_class = tuple((-entry) % length for entry in residue)
        negative_representative = registry[negative_class]
        conjugates_ok = (conjugates_ok and
                         all((residue[axis] + negative_class[axis]) % length
                             == 0 for axis in range(3)) and
                         qform(negative_representative) == qform(representative) and
                         all((-representative[axis] - negative_class[axis]) % length
                             == 0 for axis in range(3)))
    check(conjugates_ok,
          f"L={length} every character has an exact equal-norm conjugate")

rays = ((1, 0, 0), (1, 1, 0), (1, -1, 0))
check(qform((0, 0, 0)) == 0 and all(qform(ray) > 0 for ray in rays),
      "k=0 is excluded from TT while every frozen TT ray is admissible")
check(tuple(qform(ray) for ray in rays) == (3, 4, 8),
      "the frozen rays have exact reciprocal forms 3,4,8")
check(min(qform(vector) for vector in registries[5].values()
          if vector != (0, 0, 0)) == 3,
      "the exact shortest nonzero shell has q=3")
check(F(3, 4) * 3 == F(9, 4),
      "q=3 gives |k_min|^2=9 pi^2/(4 L^2 a_*^2)")

pullback_ok = True
for residue, representative in registries[5].items():
    pulled = tuple((2 * entry) % 10 for entry in residue)
    # The physical vector is represented exactly by 2*hat(n) on G_2L,
    # regardless of the minimum representative selected for reporting.
    pullback_ok = (pullback_ok and
                   all((2 * representative[axis] - pulled[axis]) % 10 == 0
                       for axis in range(3)))
check(pullback_ok, "cover pullback [n]->[2n] is exact for every G_5 class")
check(any(all(entry % 2 for entry in residue)
          for residue in registries[10]),
      "G_10 contains new odd momentum classes without G_5 ancestors")


# -------------------------------------------------------------------------
# Exact Cartesian reciprocal directions and TT projectors.

def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def matmul(first, second):
    return tuple(tuple(sum(first[row][inner] * second[inner][column]
                           for inner in range(len(second)))
                       for column in range(len(second[0])))
                 for row in range(len(first)))


def matadd(first, second):
    return tuple(tuple(first[row][column] + second[row][column]
                       for column in range(len(first[0])))
                 for row in range(len(first)))


def matscale(scale, matrix):
    return tuple(tuple(scale * value for value in row) for row in matrix)


def identity(size):
    return tuple(tuple(F(int(row == column)) for column in range(size))
                 for row in range(size))


def inverse(matrix):
    size = len(matrix)
    work = [[F(value) for value in row] + list(identity(size)[index])
            for index, row in enumerate(matrix)]
    for column in range(size):
        pivot = next(row for row in range(column, size)
                     if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        factor = work[column][column]
        work[column] = [value / factor for value in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [value - factor * base
                         for value, base in zip(work[row], work[column])]
    return tuple(tuple(row[size:]) for row in work)


def matrix_rank(matrix):
    work = [[F(value) for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work))
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        factor = work[rank][column]
        work[rank] = [value / factor for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [value - factor * base
                         for value, base in zip(work[row], work[rank])]
        rank += 1
    return rank


B = ((F(-2), F(-2), F(0)),
     (F(-2), F(0), F(-2)),
     (F(0), F(2), F(2)))
BinvT = transpose(inverse(B))


def matvec(matrix, vector):
    return tuple(sum(row[column] * vector[column]
                     for column in range(len(vector))) for row in matrix)


def outer(first, second):
    return tuple(tuple(a * b for b in second) for a in first)


def projector_for_ray(ray):
    vector = matvec(BinvT, tuple(F(value) for value in ray))
    norm = sum(value * value for value in vector)
    return matadd(identity(3), matscale(-F(1, norm), outer(vector, vector)))


SYM_BASIS = (
    ((1, 0, 0), (0, 0, 0), (0, 0, 0)),
    ((0, 0, 0), (0, 1, 0), (0, 0, 0)),
    ((0, 0, 0), (0, 0, 0), (0, 0, 1)),
    ((0, 1, 0), (1, 0, 0), (0, 0, 0)),
    ((0, 0, 1), (0, 0, 0), (1, 0, 0)),
    ((0, 0, 0), (0, 0, 1), (0, 1, 0)),
)


def symmetric_coordinates(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            matrix[0][1], matrix[0][2], matrix[1][2])


def tt_matrix(ray):
    p = projector_for_ray(ray)
    columns = []
    for basis in SYM_BASIS:
        ptp = matmul(matmul(p, basis), p)
        contraction = sum(p[row][column] * basis[row][column]
                          for row in range(3) for column in range(3))
        projected = matadd(ptp, matscale(-contraction / 2, p))
        columns.append(symmetric_coordinates(projected))
    return transpose(tuple(columns))


for index, ray in enumerate(rays):
    vector = matvec(BinvT, tuple(F(value) for value in ray))
    check(16 * sum(value * value for value in vector) == qform(ray),
          f"ray {index}: Cartesian reciprocal norm matches q/16")
    p = projector_for_ray(ray)
    check(matmul(p, p) == p and sum(p[i][i] for i in range(3)) == 2,
          f"ray {index}: transverse vector projector is exact rank two")
    tt = tt_matrix(ray)
    check(matmul(tt, tt) == tt and matrix_rank(tt) == 2,
          f"ray {index}: symmetric TT projector is idempotent rank two")


# -------------------------------------------------------------------------
# Degenerate/coalescing projectors and normalized ancestry self-tests.

P1 = ((F(1), F(0), F(0), F(0)),
      (F(0), F(0), F(0), F(0)),
      (F(0), F(0), F(0), F(0)),
      (F(0), F(0), F(0), F(0)))
P2 = ((F(0), F(0), F(0), F(0)),
      (F(0), F(1), F(0), F(0)),
      (F(0), F(0), F(0), F(0)),
      (F(0), F(0), F(0), F(0)))
PC = matadd(P1, P2)
rotation = ((F(3, 5), F(-4, 5), F(0), F(0)),
            (F(4, 5), F(3, 5), F(0), F(0)),
            (F(0), F(0), F(1), F(0)),
            (F(0), F(0), F(0), F(1)))
check(matmul(matmul(rotation, PC), transpose(rotation)) == PC,
      "complete degenerate projector is invariant under internal basis rotation")
check(matadd(P1, P2) == matadd(P2, P1) and matrix_rank(PC) == 2,
      "coalescing-projector sum is order independent and rank two")

# A degenerate two-ground-state sector is queried with rho0=P0/2.  Its trace
# Gram is invariant under internal ground-basis rotation, unlike a selected
# pure ground vector.
rho0 = ((F(1, 2), F(0)), (F(0), F(1, 2)))
Qg = ((F(1), F(0)), (F(0), F(2)))
Og = ((F(3, 5), F(-4, 5)), (F(4, 5), F(3, 5)))
Qg_rot = matmul(matmul(Og, Qg), transpose(Og))


def matrix_trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


mixed_gram = matrix_trace(matmul(rho0, matmul(Qg, transpose(Qg))))
mixed_gram_rot = matrix_trace(
    matmul(rho0, matmul(Qg_rot, transpose(Qg_rot))))
pure_gram = matmul(Qg, transpose(Qg))[0][0]
pure_gram_rot = matmul(Qg_rot, transpose(Qg_rot))[0][0]
check(mixed_gram == F(5, 2) and mixed_gram_rot == mixed_gram,
      "normalized complete-ground density is basis invariant")
check(pure_gram != pure_gram_rot,
      "regression: a selected vector in a degenerate ground space is rejected")

# F maps two TT source directions into a four-state Hilbert space.
Fmap = ((F(1), F(0)),
        (F(0), F(1)),
        (F(1), F(0)),
        (F(0), F(2)))
R = matmul(matmul(transpose(Fmap), PC), Fmap)
S = matmul(transpose(Fmap), Fmap)
check(R == ((F(1), F(0)), (F(0), F(1))) and
      S == ((F(2), F(0)), (F(0), F(5))),
      "synthetic source and cluster Grams are positive with 0<R<S")


def det2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def trace2(matrix):
    return matrix[0][0] + matrix[1][1]


SinvR = matmul(inverse(S), R)
check(trace2(SinvR) == F(7, 10) and det2(SinvR) == F(1, 10),
      "normalized ancestry generalized eigenvalues are 1/2 and 1/5")
C = ((F(2), F(1)), (F(1), F(1)))
R2 = matmul(matmul(transpose(C), R), C)
S2 = matmul(matmul(transpose(C), S), C)
SinvR2 = matmul(inverse(S2), R2)
check(trace2(SinvR2) == trace2(SinvR) and
      det2(SinvR2) == det2(SinvR),
      "ancestry generalized spectrum is invariant under TT source-basis change")

R1 = matmul(matmul(transpose(Fmap), P1), Fmap)
R2pole = matmul(matmul(transpose(Fmap), P2), Fmap)
numerator = matadd(matscale(F(2), R1), matscale(F(4), R2pole))
check(numerator == ((F(2), F(0)), (F(0), F(4))) and
      matrix_rank(numerator) == 2,
      "coalescing 2 Delta R numerator is positive rank two without eigenvector choice")


def canonical_outcome(raw_rank, controlled_power, map_frozen, lower_bound,
                      upper_bounded):
    if raw_rank < 2:
        return "FAIL"
    if not map_frozen:
        return "INDETERMINATE"
    if lower_bound > 0 and upper_bounded:
        return "PASS"
    return "FAIL"


check(canonical_outcome(2, False, True, F(1, 3), True) == "PASS",
      "derived full-rank canonical vertex yields PASS")
check(canonical_outcome(1, False, True, F(1, 3), True) == "FAIL",
      "persistent source rank loss yields FAIL")
check(canonical_outcome(2, True, False, F(0), False) == "INDETERMINATE",
      "controlled raw power without a frozen map is INDETERMINATE")
check(canonical_outcome(2, False, False, F(1), True) == "INDETERMINATE",
      "nonvanishing full-rank raw vertex without a canonical map is INDETERMINATE")
check(canonical_outcome(2, False, True, F(1), False) == "FAIL",
      "unbounded canonical upper singular value yields FAIL")

# Exact V002 overlap regression: canonical amplitude stays positive while the
# independently normalized Phase-A ancestry vanishes.
for L in (5, 10, 20, 40):
    delta_overlap = F(1, L)
    residue_overlap = F(L, 2)
    source_gram_overlap = F(L * L)
    numerator_overlap = 2 * delta_overlap * residue_overlap
    ancestry_overlap = residue_overlap / source_gram_overlap
    check(numerator_overlap == 1 and ancestry_overlap == F(1, 2 * L),
          f"L={L}: overlap regression has Z_can=1 and eta=1/(2L)")
check(canonical_outcome(2, False, True, F(1), True) == "PASS" and
      F(1, 2 * 40) < F(1, 50),
      "V002 counterexample is amplitude PASS with vanishing ancestry")
check(result["decision_separation"].endswith("overall G2 FAIL"),
      "amplitude PASS plus ancestry FAIL is overall G2 FAIL")


# -------------------------------------------------------------------------
# Shared-cone and factorization algebra self-tests (not parent evidence).

shared_slope2 = F(25, 9)
sector_slopes = {sector: shared_slope2 for sector in ("TT", "EM", "matter")}
check(len(set(sector_slopes.values())) == 1 and shared_slope2 > 0,
      "common-cone decision rule accepts one shared positive parent-unit slope")
sector_slopes["matter"] += F(1, 100)
check(len(set(sector_slopes.values())) != 1,
      "common-cone decision rule rejects a sector-dependent slope")


def helicity_two_representation(covariant, helicities):
    return covariant and tuple(sorted(helicities)) == (-2, 2)


check(helicity_two_representation(True, (-2, 2)),
      "little-group gate accepts a covariant helicity +/-2 pair")
check(not helicity_two_representation(True, (0, 0)),
      "hostile regression: degenerate factorizing scalar doublet is rejected")
check(not helicity_two_representation(False, (-2, 2)),
      "helicity labels without pole-bundle covariance are rejected")

Z = ((F(2), F(0)), (F(0), F(3)))
Zinv = inverse(Z)
amplitude = ((F(5), F(-1)), (F(2), F(7)))
residue = matmul(Z, amplitude)
check(matmul(Zinv, residue) == amplitude,
      "pole amputation recovers a finite source-independent vertex")
analytic_contact = ((F(11), F(13)), (F(17), F(19)))
check(matmul(Zinv, residue) == amplitude and analytic_contact != residue,
      "an analytic contact changes the regular term but not the pole residue")


print(f"SUMMARY {checks}/{checks} GF canonical-IR contract checks passed")
print("STATUS REPAIRED_V005_AWAITING_INDEPENDENT_AUDIT_PHYSICS_UNEXECUTED")
print("CEILING no matched-family pole data, Poincare/helicity result, positive G2, G3 soft theorem, gravity, or G")
