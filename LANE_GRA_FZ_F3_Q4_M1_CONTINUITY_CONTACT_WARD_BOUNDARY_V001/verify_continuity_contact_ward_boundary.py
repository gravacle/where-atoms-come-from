#!/usr/bin/env python3
"""Exact finite checks for the FZ continuity/contact/Ward boundary.

The default replay is deliberately inexpensive.  It reuses FY only in its
structural mode, proves an exact nonzero supplied-embedding longitudinal witness in
Q(zeta_240), checks the projected-ice charge identity, and constructs the
exact rational TT projector for the frozen m=1 direction.

``--full-liouvillian`` additionally rebuilds FY's complete H6 native source
and tests whether its embedding-longitudinal operator lies in the algebraic range of
ad_H at the two frozen samples.  Such a range solution is only a finite,
generally nonlocal operator; it is not a physical current or a Ward closure.
"""

from argparse import ArgumentParser
from contextlib import redirect_stdout
from fractions import Fraction as F
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path
import runpy
import sys

import numpy as np


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
FY_DIR = ROOT / "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001"
FY_SCRIPT = FY_DIR / "derive_native_support_m1_response.py"

parser = ArgumentParser()
parser.add_argument("--full-liouvillian", action="store_true")
args = parser.parse_args()

checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


for dependency_line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = dependency_line.split("  ", 1)
    dependency = ROOT / relative
    check(dependency.is_file() and not dependency.is_symlink() and
          sha256(dependency.read_bytes()).hexdigest() == expected,
          f"dependency custody: {relative}")

fy_theorem = (FY_DIR / "THEOREM.md").read_text()
fu_theorem = (ROOT /
    "LANE_GRA_FU_F3_Q4_PAIR_RESOLVED_MAXWELL_DPAR_DERIVATION_V001" /
    "THEOREM.md").read_text()
fv_theorem = (ROOT /
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001" /
    "THEOREM.md").read_text()
check("FY has not constructed the temporal density,\ncurrent operators, "
      "contact terms" in fy_theorem,
      "FY explicitly leaves temporal/current/contact operators unconstructed")
check("[Q_R,T_{a,\\pm}]=\\pm2q_*T_{a,\\pm}" in fu_theorem and
      "[Q_{\\rm tot},\\widetilde X_a]=0" in fu_theorem,
      "FU09b explicitly owns the required conserved dressed-transfer algebra")
check("S10 / FV-PURE" in fv_theorem and
      "No additional nonidentity source operator is present" in fv_theorem,
      "FV-PURE freezes the complete projected nonidentity spatial source")


# Load the exact FY ledgers without triggering its long H4/H6 enumeration.
saved_argv = sys.argv[:]
try:
    sys.argv = [str(FY_SCRIPT)]
    with redirect_stdout(StringIO()):
        fy = runpy.run_path(str(FY_SCRIPT))
finally:
    sys.argv = saved_argv

fo = fy["fo"]
states = fo["states"]
incidence = fo["incidence"]
H = np.asarray(fo["hamiltonian"], dtype=complex)
dimension = len(states)

check((fo["CELL_COUNT"], fo["VERTEX_COUNT"], len(fo["edges"]), dimension) ==
      (30, 60, 120, 180),
      "FZ uses exactly the frozen FO Z30/180-state parent")


# -------------------------------------------------------------------------
# Exact projected incidence-charge identity.

degree_rows = tuple(tuple(sum((state >> edge) & 1 for _, edge in
                              incidence[vertex])
                          for vertex in range(fo["VERTEX_COUNT"]))
                    for state in states)
check(all(all(degree == 2 for degree in row) for row in degree_rows),
      "every projected state has exact local incidence charge d_v-2=0")
check(all(sum((state >> edge) & 1 for edge in range(len(fo["edges"]))) == 60
          for state in states),
      "every projected state has exactly sixty occupied links")
check(np.count_nonzero(H) == 840 and np.array_equal(H, H.T),
      "the nontrivial H6 ring generator is the frozen 420-edge Hermitian graph")

# P G_v P is the zero matrix for every vertex, hence both [H,G_v] and its
# Heisenberg derivative vanish exactly.  This is a genuine but trivial charge
# continuity law; it supplies no temporal partner for the stress source.
projected_charges = tuple(
    np.diag([row[vertex] - 2 for row in degree_rows])
    for vertex in range(fo["VERTEX_COUNT"])
)
check(all(not np.any(charge) for charge in projected_charges),
      "all sixty projected Gauss-charge operators vanish exactly")
check(all(degree_rows[row] == degree_rows[column]
          for row, column in zip(*np.nonzero(H))),
      "every H6 transition preserves all sixty local incidence charges")

# Minimal concrete representation of the already-required FU09b dressed
# transfer (q_*=h=hbar=1).  It is not installed in H6; it verifies the exact
# current seed identified as the missing physical dependency.
z_link = np.diag((-1, 1)).astype(complex)
sigma_plus = np.array(((0, 0), (1, 0)), dtype=complex)
sigma_minus = sigma_plus.T
q_reservoir = np.diag((1, -1)).astype(complex)
t_minus = np.array(((0, 0), (1, 0)), dtype=complex)
t_plus = t_minus.T
identity2 = np.eye(2, dtype=complex)
a_transfer = np.kron(sigma_plus, t_minus)
b_transfer = np.kron(sigma_minus, t_plus)
x_dressed = a_transfer + b_transfer
q_link_full = np.kron(z_link, identity2)
q_reservoir_full = np.kron(identity2, q_reservoir)
q_total = q_link_full + q_reservoir_full
h_flip = -x_dressed
q_link_dot = 1j*(h_flip@q_link_full-q_link_full@h_flip)
q_reservoir_dot = 1j*(h_flip@q_reservoir_full-
                       q_reservoir_full@h_flip)
check(np.array_equal(q_total@x_dressed, x_dressed@q_total),
      "FU09b dressed transfer commutes exactly with total charge")
check(np.array_equal(q_link_dot, 2j*(a_transfer-b_transfer)) and
      np.array_equal(q_reservoir_dot, -q_link_dot) and
      np.any(q_link_dot),
      "FU09b supplies an exact nonzero compensating transfer-current seed")


# -------------------------------------------------------------------------
# Exact nonzero longitudinal witness over Q(zeta_240).

# The frozen shortest reciprocal alias has k proportional to (7,15,-17).
# Coordinates use (xx,yy,zz,2xy,2xz,2yz); multiplying by two removes halves.
K_DIRECTION = (7, 15, -17)
LONGITUDINAL_COVECTORS_2 = (
    (14, 0, 0, 15, -17, 0),
    (0, 30, 0, 7, 0, -17),
    (0, 0, -34, 0, 7, 15),
)
check(sum(value * value for value in K_DIRECTION) == 563,
      "the exact m=1 direction has squared integer norm 563")
frozen_q, frozen_wavevector, _, _, _, _ = fo["polarization_data"](1)
check(np.allclose(frozen_q, np.array((1, 5, -11))/30) and
      np.allclose(frozen_wavevector,
                  np.pi*np.sqrt(3)/60*np.array(K_DIRECTION)),
      "the FO embedding gives exactly k=(pi sqrt(3)/60)(7,15,-17)")


def add_polynomials(terms):
    length = max(len(polynomial) for _, polynomial in terms)
    result = [F(0)] * length
    for scale, polynomial in terms:
        for power, value in enumerate(polynomial):
            result[power] += F(scale) * value
    return tuple(result)


direct = fy["direct_ledgers"][0]
coordinate_polynomials = tuple(
    fy["ledger_component_polynomial"](direct, component)
    for component in range(6)
)
longitudinal_remainders = []
for covector in LONGITUDINAL_COVECTORS_2:
    polynomial = add_polynomials(tuple(
        (coefficient, coordinate_polynomials[component])
        for component, coefficient in enumerate(covector)
        if coefficient
    ))
    longitudinal_remainders.append(
        fy["polynomial_reduce_phi240"](polynomial))
check(any(any(remainder) for remainder in longitudinal_remainders),
      "Qpair(m=1) has an exact nonzero supplied-embedding contraction")

# Freeze a compact exact checksum for the first nonzero component.
first_nonzero = next(remainder for remainder in longitudinal_remainders
                     if any(remainder))
first_terms = tuple((power, value) for power, value in enumerate(first_nonzero)
                    if value)
check(first_terms and all(isinstance(value, F) for _, value in first_terms),
      "the longitudinal witness is exact rational cyclotomic data")

# FY's structural replay already contains one complete 720+720 Hermitian ring
# ledger.  Its off-diagonal longitudinal contraction can be decided exactly
# without rerunning the exhaustive fourteen-orbit inventory.
ring_coordinate_polynomials = tuple(
    fy["ledger_component_polynomial"](fy["_witness_q"], component)
    for component in range(6)
)
ring_longitudinal_remainders = []
for covector in LONGITUDINAL_COVECTORS_2:
    polynomial = add_polynomials(tuple(
        (coefficient, ring_coordinate_polynomials[component])
        for component, coefficient in enumerate(covector)
        if coefficient
    ))
    ring_longitudinal_remainders.append(
        fy["polynomial_reduce_phi240"](polynomial))
check(any(any(remainder) for remainder in ring_longitudinal_remainders),
      "one complete ring entry has exact nonzero embedding contraction")
ring_first_nonzero = next(remainder for remainder in
                          ring_longitudinal_remainders if any(remainder))
ring_first_terms = tuple((power, value) for power, value in
                         enumerate(ring_first_nonzero) if value)
check(ring_first_terms and all(isinstance(value, F)
                               for _, value in ring_first_terms),
      "the ring longitudinal witness is exact rational cyclotomic data")

coefficients = (F(-1), F(-37, 12), F(-16247, 900))
samples = ((F(2, 5), F(15625, 504), F(2415673, 3515625)),
           (F(1, 2), F(512, 63), F(15853, 57600)))
for x, rho, expected_f in samples:
    f_e = 1 + coefficients[0]*x**2 + coefficients[1]*x**4 + coefficients[2]*x**6
    check(f_e == expected_f and rho*f_e != 0,
          f"x={x}: complete diagonal H6 source has exact nonzero rho*f_E")
check(fy["exact_m1_relation"](fy["qdiag2_ledgers"],
                              fy["direct_ledgers"], coefficients[0]),
      "the inexpensive exact H2 ledger lifts the first f_E coefficient")

# FY proves that the ring source is strictly off diagonal.  Therefore it
# cannot cancel the nonzero diagonal longitudinal witness at either sample.
fy_result = json.loads((FY_DIR / "RESULT.json").read_text())
check(fy_result["exact_m1_diagonal_lift"] == [
          {"order": 2, "coefficient": "-1"},
          {"order": 4, "coefficient": "-37/12"},
          {"order": 6, "coefficient": "-16247/900"}],
      "FY custody fixes all three exact m=1 diagonal lift coefficients")
check(fy_result["exact_ring_result"] ==
      "nonzero independent off-diagonal m=1 source",
      "FY custody fixes the independent off-diagonal ring source")


# -------------------------------------------------------------------------
# Exact rational transverse-traceless projector for k || (7,15,-17).


def matrix_multiply(first, second):
    return tuple(tuple(sum(first[row][inner] * second[inner][column]
                           for inner in range(len(second)))
                       for column in range(len(second[0])))
                 for row in range(len(first)))


def matrix_add(first, second):
    return tuple(tuple(first[row][column] + second[row][column]
                       for column in range(len(first[0])))
                 for row in range(len(first)))


def matrix_scale(scale, matrix):
    return tuple(tuple(F(scale) * value for value in row) for row in matrix)


def transpose(matrix):
    return tuple(zip(*matrix))


identity3 = tuple(tuple(F(int(row == column)) for column in range(3))
                  for row in range(3))
r = tuple(F(value) for value in K_DIRECTION)
outer_rr = tuple(tuple(r[row] * r[column] for column in range(3))
                 for row in range(3))
P = matrix_add(identity3, matrix_scale(F(-1, 563), outer_rr))
check(matrix_multiply(P, P) == P,
      "the exact spatial transverse projector is idempotent")
check(all(sum(P[row][column] * r[column] for column in range(3)) == 0
          for row in range(3)),
      "the exact transverse projector annihilates the m=1 direction")
check(sum(P[index][index] for index in range(3)) == 2,
      "the transverse plane has exact rank two")


def unpack_coordinate(column):
    values = [F(0)] * 6
    values[column] = F(1)
    return ((values[0], values[3]/2, values[4]/2),
            (values[3]/2, values[1], values[5]/2),
            (values[4]/2, values[5]/2, values[2]))


def pack_coordinate(matrix):
    return (matrix[0][0], matrix[1][1], matrix[2][2],
            2*matrix[0][1], 2*matrix[0][2], 2*matrix[1][2])


def tt_project(matrix):
    pap = matrix_multiply(matrix_multiply(P, matrix), P)
    transverse_trace = sum(P[row][column] * matrix[row][column]
                           for row in range(3) for column in range(3))
    return matrix_add(pap, matrix_scale(-transverse_trace/2, P))


tt_columns = tuple(pack_coordinate(tt_project(unpack_coordinate(column)))
                   for column in range(6))
TT = transpose(tt_columns)
check(matrix_multiply(TT, TT) == TT,
      "the exact six-coordinate TT projector is idempotent")
for column in range(6):
    projected = tt_project(unpack_coordinate(column))
    check(sum(projected[index][index] for index in range(3)) == 0,
          f"TT basis image {column} is exactly traceless")
    check(all(sum(r[row] * projected[row][column2] for row in range(3)) == 0
              for column2 in range(3)),
          f"TT basis image {column} is exactly transverse")


def rational_rank(matrix):
    work = [list(row) for row in matrix]
    rows, columns = len(work), len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [value/scale for value in work[rank]]
        for row in range(rows):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [value - factor*pivot_value
                         for value, pivot_value in zip(work[row], work[rank])]
        rank += 1
    return rank


check(rational_rank(TT) == 2,
      "the exact m=1 TT projector has coordinate rank two")

# FY's sampled response supplies both directions in that exact kinematic
# quotient, but all poles are finite.  This is response evidence, not Ward
# ancestry or a massless graviton.
for sample in fy_result["samples"]:
    check(sample["ranks"]["TT_ground_image"] == 2 and
          min(sample["TT_gram_eigenvalues"]) > 0,
          f"x={sample['x']}: sampled complete source excites two TT images")
    check(min(sample["pole_gaps"]) > 0,
          f"x={sample['x']}: sampled finite response has no zero pole")


# -------------------------------------------------------------------------
# Optional complete Liouvillian range diagnostic.

liouvillian_results = []
if args.full_liouvillian:
    with redirect_stdout(StringIO()):
        _, qdiag_ledgers = fy["enumerate_diagonal_ledgers"]()
        _, _, ring_ledgers = fy["enumerate_ring_ledgers"]()
    direct_m1 = fy["diagonal_source_from_ledgers"](
        fy["direct_ledgers"], 1)
    qdiag_m1 = tuple(fy["diagonal_source_from_ledgers"](ledgers, 1)
                     for ledgers in qdiag_ledgers)
    ring_m1 = fy["ring_source_from_ledgers"](ring_ledgers, 1)
    _, wavevector, _, _, _, _ = fo["polarization_data"](1)
    eigenvalues = np.asarray(fo["eigenvalues"], dtype=float)
    eigenvectors = np.asarray(fo["eigenvectors"], dtype=complex)
    delta = eigenvalues[:, None] - eigenvalues[None, :]
    degenerate = np.abs(delta) < 2e-10
    for x, rho, _ in samples:
        coordinate_source = float(rho) * (
            direct_m1 + float(x**2)*qdiag_m1[0] +
            float(x**4)*qdiag_m1[1] + float(x**6)*qdiag_m1[2]) + ring_m1
        longitudinal = fy["longitudinal_sources"](coordinate_source,
                                                   wavevector)
        transformed = np.asarray([
            eigenvectors.conj().T @ operator @ eigenvectors
            for operator in longitudinal
        ])
        blocked_norm = float(np.linalg.norm(transformed[:, degenerate]))
        total_norm = float(np.linalg.norm(transformed))
        relative = blocked_norm/max(total_norm, 1e-300)
        admits = relative < 2e-10
        liouvillian_results.append({
            "x": str(x),
            "degenerate_block_relative_norm": relative,
            "embedding_longitudinal_algebraic_adH_range": admits,
        })
        check(np.isfinite(relative),
              f"x={x}: Liouvillian compatibility diagnostic is finite")

print("EXACT_LONGITUDINAL_WITNESS", [(power, str(value))
      for power, value in first_terms])
print("EXACT_RING_LONGITUDINAL_WITNESS", [(power, str(value))
      for power, value in ring_first_terms])
if liouvillian_results:
    print("LIouvillian_RESULT", json.dumps(liouvillian_results,
                                           sort_keys=True))
print(f"SUMMARY {checks}/{checks} FZ continuity/contact/Ward checks passed")
print("MODE", "FULL_LIOUVILLIAN" if args.full_liouvillian else "EXACT_FAST")
print("CEILING supplied-embedding contraction and TT algebra only; no derived "
      "physical discrete divergence, temporal/current/contact source, Ward "
      "closure, continuum, gravity, or G")
