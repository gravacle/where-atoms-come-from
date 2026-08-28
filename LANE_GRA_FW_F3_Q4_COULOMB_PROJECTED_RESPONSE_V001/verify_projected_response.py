#!/usr/bin/env python3
"""Finite projected-response screen for the FV-WITNESS source.

The script replays the frozen FO 180-state translation-closed H6 ice sector,
constructs the two rank-closing source pieces displayed by FV-PURE -- the
order-zero direct pair term and irreducible differentiated H6 ring term --
and separates four notions that must not be conflated:

* operator rank modulo the identity;
* rank of the source-to-commutator map ``Q -> [H,Q]``;
* ground-state spectral/retarded rank; and
* the rank of the first nonzero ground-state commutator moment.

Units are J6=1 and hbar=1.  The direct Coulomb-DPAR pair source is measured
in units of Ud, so its physical coefficient relative to the ring source is
rho=Ud/J6.  The executable certificate uses rho=1 and verifies the separate
direct/ring factors from which the displayed all-rho formulas follow.
"""

from contextlib import redirect_stdout
from hashlib import sha256
from io import StringIO
from itertools import combinations
from pathlib import Path
import runpy

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
LANE = Path(__file__).resolve().parent
FO_SCRIPT = (ROOT / "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001" /
             "verify_finite_tt_four_point.py")

DEPENDENCIES = {
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/THEOREM.md":
        "44fc28edc9820d2b4ea67cef9f83beef60e53bfe6407b582ffbeccfe42f756c5",
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/INDEPENDENT_AUDIT.md":
        "84d8c02c3e560198f6a9fae04f5ee81bc72c98354e02f1f58d506a8c3171c453",
    "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/verify_finite_tt_four_point.py":
        "fb44d45290c0530098c0e8f9593dff1c0f8149d42598f842374b61004a8ff6c2",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/THEOREM.md":
        "6fc221a31151340b91a946d33e442971c1373500e067c354b6c610e3964edb1c",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/RESULT.md":
        "b5d4c3de99aa4e100519c19a9b74de487b47c1a2d3671204e77740bd9094771a",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/verify_projected_source_rank.py":
        "0e93d84f9eb7cf7fdd62b5a14d5c6705c74841899dd1676bdd7e7a41eb971a00",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/MANIFEST.sha256":
        "651a66b9afd7545b04aa80e5f90952fda9327d011ecd19b973aa80ff51a739f3",
    "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001/SEAL.sha256":
        "8301da6bbc026d0e14d985592c5dabe3d91072957c4aaa4b1bebf1f45aadd894",
}

MANIFEST_FILES = {
    "DEPENDENCIES.sha256", "README.md", "RESULT.md", "SELF_AUDIT.md",
    "THEOREM.md", "VERIFICATION.txt", "verify_projected_response.py",
}


checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


for relative, expected in DEPENDENCIES.items():
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink(),
          f"dependency is a regular file: {relative}")
    check(digest(path) == expected, f"dependency hash is pinned: {relative}")


# Reuse the already audited FO quotient/sector constructor byte-for-byte.  Its
# own verification output is suppressed here; every object consumed below is
# checked again at the interface.
captured = StringIO()
with redirect_stdout(captured):
    fo = runpy.run_path(str(FO_SCRIPT))

states = fo["states"]
state_index = fo["state_index"]
incidence = fo["incidence"]
edge_labels = fo["edge_labels"]
hexagons = fo["hexagons"]
ring_patterns = fo["ring_patterns"]
H = fo["hamiltonian"]
eigenvalues = fo["eigenvalues"]
eigenvectors = fo["eigenvectors"]
ground = fo["ground"]
ground_energy = fo["ground_energy"]
translation_permutation = fo["translation_permutation"]
translation_orbits = fo["translation_orbits"]
dimension = len(states)

check(dimension == 180, "FO interface supplies exactly the 180-state sector")
check(H.shape == (180, 180) and np.array_equal(H, H.T),
      "FO interface supplies the exact real-symmetric H6 matrix")
check(len(hexagons) == 120 and len(ring_patterns) == 120,
      "FO interface supplies all 120 elementary hexagons")
check(abs(ground_energy + 2 + 2*np.sqrt(2)) < 2e-12,
      "FO ground energy retains its exact algebraic checksum")


# Tensor coordinates are (xx, yy, zz, 2xy, 2xz, 2yz), exactly as in FV.
SIGNS = np.array(((1, 1, 1), (1, -1, -1),
                  (-1, 1, -1), (-1, -1, 1)), dtype=int)
PAIRS = tuple(combinations(range(4), 2))


def dyad_coordinates(vector, denominator):
    x, y, z = vector
    return np.array((x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z),
                    dtype=float) / denominator


EDGE_DYADS = tuple(dyad_coordinates(vector, 3) for vector in SIGNS)
ROOT_DYADS = {
    (a, b): dyad_coordinates(SIGNS[b] - SIGNS[a], 8)
    for a, b in PAIRS
}

# Direct source in Ud units: lambda=-1/2 times every normalized root dyad.
Q_direct = np.zeros((6, dimension, dimension), dtype=float)
for row, state in enumerate(states):
    for vertex in range(fo["VERTEX_COUNT"]):
        z = {edge_labels[edge]: 1 - 2*((state >> edge) & 1)
             for _, edge in incidence[vertex]}
        for a, b in PAIRS:
            Q_direct[:, row, row] += (-0.5 * z[a] * z[b] *
                                      ROOT_DYADS[(a, b)])

# Differentiated H6 ring source in J6 units.  Multiplying FV11 by -8/63
# gives -31 I/6 + 9 D_d/2 for a ring whose missing label is d.
Q_ring = np.zeros((6, dimension, dimension), dtype=float)
isotropic = np.array((1, 1, 1, 0, 0, 0), dtype=float)
for row, state in enumerate(states):
    for (mask, first, second), cycle in zip(ring_patterns, hexagons):
        if (state & mask) not in (first, second):
            continue
        column = state_index[state ^ mask]
        present = {edge_labels[edge] for edge in cycle}
        missing = next(iter(set(range(4)) - present))
        Q_ring[:, row, column] = (-31/6 * isotropic +
                                  9/2 * EDGE_DYADS[missing])

check(np.array_equal(Q_direct, Q_direct.transpose(0, 2, 1)),
      "direct Coulomb-DPAR source is exactly diagonal Hermitian")
check(np.array_equal(Q_ring, Q_ring.transpose(0, 2, 1)),
      "Hermitian FV ring source is exactly real symmetric")
check(np.count_nonzero(Q_ring[0]) == 840,
      "ring-source support is the two orientations of all 420 transitions")

# Independent analytic zero-momentum reduction.  The free translation action
# has six length-30 orbits, so every homogeneous source reduces to a 6 x 6
# matrix.  In this basis H0 is integer, the ground vector is closed-form, and
# the projectors onto the two responding excited energies are rational
# polynomials in H0.  This avoids inferring the algebraic residues solely from
# floating-point eigenvector choices inside degenerate spaces.
U0 = np.zeros((dimension, 6), dtype=float)
for orbit_index, orbit in enumerate(translation_orbits):
    for state in orbit:
        U0[state_index[state], orbit_index] = 1/np.sqrt(30)
check(np.linalg.norm(U0.T@U0 - np.eye(6)) < 2e-15,
      "translation-orbit zero-momentum isometry is orthonormal")
H0 = U0.T@H@U0
EXPECTED_H0 = np.array((
    (0, -1, -1, -1, -1, -2),
    (-1, 0, -1, -1, 0, -1),
    (-1, -1, 0, 0, -1, -1),
    (-1, -1, 0, 0, -1, -1),
    (-1, 0, -1, -1, 0, -1),
    (-2, -1, -1, -1, -1, 0),
), dtype=float)
check(np.linalg.norm(H0 - EXPECTED_H0) < 2e-15,
      "zero-momentum Hamiltonian is the exact integer 6x6 block")
ground0 = np.array((1/2, 1/(2*np.sqrt(2)), 1/(2*np.sqrt(2)),
                    1/(2*np.sqrt(2)), 1/(2*np.sqrt(2)), 1/2))
check(abs(ground0@ground0 - 1) < 2e-15 and
      np.linalg.norm(H0@ground0 - (-2-2*np.sqrt(2))*ground0) < 2e-15,
      "zero-momentum ground vector is certified in closed algebraic form")
P_ZERO = ((H0@H0 + 4*H0 - 4*np.eye(6)) @ (H0 - 2*np.eye(6))) / 8
P_TWO = ((H0@H0 + 4*H0 - 4*np.eye(6)) @ H0) / 16
check(np.linalg.norm(P_ZERO@P_ZERO-P_ZERO) < 2e-14 and
      np.linalg.norm(P_TWO@P_TWO-P_TWO) < 2e-14,
      "rational polynomial projectors onto energies zero and two are idempotent")
check(round(np.trace(P_ZERO)) == 2 and round(np.trace(P_TWO)) == 2 and
      np.linalg.norm(P_ZERO@P_TWO) < 2e-14,
      "both responding zero-momentum energy projectors have exact rank two")

# Uniform source must remain in the zero cyclic-momentum block.
for source in tuple(Q_direct) + tuple(Q_ring):
    check(np.array_equal(source,
                         source[np.ix_(translation_permutation,
                                       translation_permutation)]),
          "each homogeneous tensor source commutes with cyclic translation")


# Frobenius-orthonormal Sym^2 basis: A1, E1, E2, Txy, Txz, Tyz.
SQRT2 = np.sqrt(2.0)
SQRT3 = np.sqrt(3.0)
SQRT6 = np.sqrt(6.0)
BASIS_COORDINATES = np.array((
    (1/SQRT3, 1/SQRT3, 1/SQRT3, 0, 0, 0),
    (1/SQRT2, -1/SQRT2, 0, 0, 0, 0),
    (1/SQRT6, 1/SQRT6, -2/SQRT6, 0, 0, 0),
    (0, 0, 0, 1/SQRT2, 0, 0),
    (0, 0, 0, 0, 1/SQRT2, 0),
    (0, 0, 0, 0, 0, 1/SQRT2),
), dtype=float)
NAMES = ("A1", "E1", "E2", "Txy", "Txz", "Tyz")


def contract_sources(coordinate_sources):
    return np.einsum("ac,cij->aij", BASIS_COORDINATES,
                     coordinate_sources)


D = contract_sources(Q_direct)
R = contract_sources(Q_ring)

# The uniform orbit sums cancel the 1/30 normalization.  Direct entries lie
# on the quarter lattice and ring entries on the third lattice, hence the
# exact common denominator is twelve.  Recover those exact 6 x 6 coordinate
# blocks before introducing the orthonormal radical tensor basis.
D0_COORD_RAW = np.array([U0.T@source@U0 for source in Q_direct])
R0_COORD_RAW = np.array([U0.T@source@U0 for source in Q_ring])
D0_COORD = np.rint(12*D0_COORD_RAW)/12
R0_COORD = np.rint(12*R0_COORD_RAW)/12
check(np.linalg.norm(D0_COORD_RAW-D0_COORD) < 2e-12 and
      np.linalg.norm(R0_COORD_RAW-R0_COORD) < 2e-12,
      "zero-momentum source blocks lie on the exact denominator-twelve lattice")

check(np.linalg.norm(D[3:]) < 1e-13,
      "direct projected pair source contains no T2 operator")
check(np.linalg.norm(R[1:3]) < 1e-13,
      "differentiated ring source contains no E operator")


def matrix_rank_psd(matrix, tolerance=1e-9):
    return int(np.count_nonzero(np.linalg.eigvalsh(
        (matrix + matrix.T.conj()) / 2) > tolerance))


def response_packet(rho):
    """Return all state-specific response objects for positive rho=Ud/J6."""
    Q = rho*D + R
    identity = np.eye(dimension)

    # Exact finite-sector identities and conservation laws.
    check(np.linalg.norm(SQRT3*Q[0] -
                         (60*rho*identity + 11*H)) < 2e-12,
          f"A1 source equals (60 rho I + 11 H)/sqrt(3) at rho={rho:g}")

    E_conserved = SQRT3/2*Q[1] + 0.5*Q[2]
    E_active = -0.5*Q[1] + SQRT3/2*Q[2]
    T_plus = (Q[3] + Q[4]) / SQRT2
    T_minus = (Q[3] - Q[4]) / SQRT2
    check(np.linalg.norm(E_conserved - 16*SQRT6*rho*identity) < 3e-12,
          f"one E direction is an identity on this component at rho={rho:g}")
    check(np.linalg.norm(H@Q[0] - Q[0]@H) < 2e-12,
          f"A1 is conserved at homogeneous momentum for rho={rho:g}")
    check(np.linalg.norm(H@Q[5] - Q[5]@H) < 2e-12,
          f"one T2 direction is exactly conserved for rho={rho:g}")
    check(np.linalg.norm((Q[5]@ground) -
                         (3*SQRT2 - 6)*ground) < 3e-12,
          f"conserved Tyz has the exact ground eigenvalue at rho={rho:g}")

    # Off-shell operator rank after identities are quotiented.
    centered_trace = np.array([
        operator - np.trace(operator)/dimension*identity for operator in Q
    ])
    operator_gram = np.einsum("aij,bij->ab", centered_trace, centered_trace)
    check(matrix_rank_psd(operator_gram) == 5,
          f"uniform 180-state operator rank modulo identity is five at rho={rho:g}")

    # Source-to-commutator map, a state-independent dynamical precursor.
    commutators = np.array([H@operator - operator@H for operator in Q])
    commutator_gram = np.einsum("aij,bij->ab", commutators, commutators)
    expected_commutator_gram = np.zeros((6, 6))
    expected_commutator_gram[1:3, 1:3] = rho*rho*np.array((
        (960, -960*SQRT3), (-960*SQRT3, 2880)))
    expected_commutator_gram[3, 3] = 25920
    expected_commutator_gram[4, 4] = 25920
    check(np.linalg.norm(commutator_gram - expected_commutator_gram) < 2e-9,
          f"Hilbert-Schmidt commutator Gram has its exact block form at rho={rho:g}")
    check(matrix_rank_psd(commutator_gram) == 3,
          f"source-to-commutator map has rank three at rho={rho:g}")

    # Ground-state spectral amplitudes.
    excited_values = eigenvalues[1:]
    excited_vectors = eigenvectors[:, 1:]
    gaps = excited_values - ground_energy
    amplitudes = np.stack([excited_vectors.T@operator@ground
                           for operator in Q], axis=1)
    delta_1 = 2 + 2*SQRT2
    delta_2 = 4 + 2*SQRT2
    mask_1 = np.abs(gaps - delta_1) < 2e-10
    mask_2 = np.abs(gaps - delta_2) < 2e-10
    residue_1 = amplitudes[mask_1].T@amplitudes[mask_1]
    residue_2 = amplitudes[mask_2].T@amplitudes[mask_2]
    residue_all = amplitudes.T@amplitudes

    vector_1 = np.array((0, rho/SQRT2, -rho*np.sqrt(3/2),
                         -3/SQRT2, -3/SQRT2, 0))
    vector_2 = np.array((0, 0, 0, 3/SQRT2, -3/SQRT2, 0))
    exact_residue_1 = np.outer(vector_1, vector_1)
    exact_residue_2 = np.outer(vector_2, vector_2)
    check(np.linalg.norm(residue_1 - exact_residue_1) < 2e-11,
          f"first pole has the exact rank-one residue at rho={rho:g}")
    check(np.linalg.norm(residue_2 - exact_residue_2) < 2e-11,
          f"second pole has the exact rank-one residue at rho={rho:g}")
    check(np.linalg.norm(residue_all - residue_1 - residue_2) < 2e-11,
          f"no other ground-state spectral support survives at rho={rho:g}")
    check(matrix_rank_psd(residue_1) == 1 and
          matrix_rank_psd(residue_2) == 1 and
          matrix_rank_psd(residue_all) == 2,
          f"ground-state spectral rank is exactly two at rho={rho:g}")
    check(abs(np.trace(residue_1@residue_2)) < 2e-11,
          f"the two rank-one residue directions are orthogonal at rho={rho:g}")

    # Independent algebraic replay using only the displayed closed-form
    # ground vector and rational polynomial projectors of the exact H0 block.
    Q0 = np.einsum("ac,cij->aij", BASIS_COORDINATES,
                   rho*D0_COORD + R0_COORD)
    source_vectors0 = np.array([operator@ground0 for operator in Q0])
    algebraic_residue_1 = source_vectors0@P_ZERO@source_vectors0.T
    algebraic_residue_2 = source_vectors0@P_TWO@source_vectors0.T
    check(np.linalg.norm(algebraic_residue_1-exact_residue_1) < 5e-11,
          f"polynomial-projector replay certifies the first algebraic residue at rho={rho:g}")
    check(np.linalg.norm(algebraic_residue_2-exact_residue_2) < 5e-11,
          f"polynomial-projector replay certifies the second algebraic residue at rho={rho:g}")
    expectations0 = np.array([ground0@operator@ground0 for operator in Q0])
    centered_vectors0 = source_vectors0 - expectations0[:, None]*ground0
    algebraic_remainder = centered_vectors0 - (
        centered_vectors0@P_ZERO + centered_vectors0@P_TWO)
    check(np.linalg.norm(algebraic_remainder) < 5e-12,
          f"analytic zero-momentum source vectors have no third-pole remainder at rho={rho:g}")

    # A nonconserved source can still be dark on one selected state.  This is
    # the fourth ground-response null and must not be called a Ward identity.
    dark = 3*E_active - SQRT2*rho*T_plus
    dark_centered = dark - (ground@dark@ground)*identity
    check(np.linalg.norm(dark_centered@ground) < 3e-12,
          f"one nonconserved source is exactly ground-dark at rho={rho:g}")
    check(np.linalg.norm(H@dark - dark@H) > 1.0,
          f"ground-dark source is not an operator conservation law at rho={rho:g}")

    # First commutator moments as defined in FQ18a.
    moment_0 = np.empty((6, 6), dtype=float)
    moment_1 = np.empty((6, 6), dtype=float)
    for a in range(6):
        ad = H@Q[a] - Q[a]@H
        for b in range(6):
            moment_0[a, b] = ground@(Q[a]@Q[b] - Q[b]@Q[a])@ground
            moment_1[a, b] = ground@(ad@Q[b] - Q[b]@ad)@ground
    exact_moment_1 = -2*(delta_1*exact_residue_1 +
                         delta_2*exact_residue_2)
    check(np.linalg.norm(moment_0) < 2e-12,
          f"equal-time ground commutator moment M0 vanishes at rho={rho:g}")
    check(np.linalg.norm(moment_1 - exact_moment_1) < 2e-10,
          f"M1 equals minus twice the gap-weighted residue at rho={rho:g}")
    check(matrix_rank_psd(-moment_1) == 2,
          f"first nonzero ground commutator moment has rank two at rho={rho:g}")

    # Positive static Kubo matrix and one off-axis retarded-frequency replay.
    kubo = 2*exact_residue_1/delta_1 + 2*exact_residue_2/delta_2
    expected_kubo_eigenvalues = sorted((
        2*(2*rho*rho + 9)/delta_1,
        18/delta_2,
    ))
    observed_kubo_eigenvalues = [value for value in np.linalg.eigvalsh(kubo)
                                 if value > 1e-9]
    check(np.allclose(observed_kubo_eigenvalues,
                      expected_kubo_eigenvalues, atol=2e-12, rtol=0),
          f"static Kubo matrix has exactly two positive eigenvalues at rho={rho:g}")

    z = 1.25 + 0.2j
    spectral_retarded = np.zeros((6, 6), dtype=complex)
    for gap, amplitude_row in zip(gaps, amplitudes):
        spectral_retarded += np.outer(amplitude_row, amplitude_row) * (
            1/(z-gap) - 1/(z+gap))
    two_pole_retarded = exact_residue_1*(1/(z-delta_1)-1/(z+delta_1))
    two_pole_retarded += exact_residue_2*(1/(z-delta_2)-1/(z+delta_2))
    check(np.linalg.norm(spectral_retarded - two_pole_retarded) < 3e-12,
          f"full Lehmann sum equals the two-pole retarded formula at rho={rho:g}")

    return {
        "Q": Q,
        "operator_rank": matrix_rank_psd(operator_gram),
        "commutator_rank": matrix_rank_psd(commutator_gram),
        "spectral_rank": matrix_rank_psd(residue_all),
        "moment_rank": matrix_rank_psd(-moment_1),
        "delta_1": delta_1,
        "delta_2": delta_2,
        "residue_1_eigenvalue": 2*rho*rho + 9,
        "residue_2_eigenvalue": 9,
        "kubo_eigenvalues": expected_kubo_eigenvalues,
    }


# rho=1 is the displayed normalized certificate.  rho=2 is an independent
# scale replay proving that the direct/ring separation, rather than a tuned
# equality of Ud and J6, drives the formulas.
packet_1 = response_packet(1.0)
packet_2 = response_packet(2.0)


# Documentary scope and the new FV-PURE conditionality are executable gates.
theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
joined = " ".join((theorem + result + self_audit).split())
for phrase in (
    "FV-PURE",
    "FV-WITNESS",
    "S1`--`S10",
    "Q_diag^(2,4,6)",
    "off-shell projected operator rank",
    "source-to-commutator rank",
    "ground-state retarded rank",
    "not a Ward identity",
    "finite-sector, homogeneous-source result",
    "does not contradict FV",
    "thermodynamic massless tensor",
    "gravity",
):
    check(phrase in joined, f"scope text retains: {phrase}")

for forbidden in (
    "FV rank six implies retarded rank six",
    "the rank-two response proves a graviton",
    "the ground-dark source is a Ward identity",
    "this finite sector excludes a thermodynamic tensor phase",
    "Newton's constant is calculated",
):
    check(forbidden not in joined, f"forbidden promotion absent: {forbidden}")

# Activated only after the builder payload is frozen.  The verifier never
# rewrites its own custody files.
manifest = LANE / "MANIFEST.sha256"
if manifest.is_file():
    listed = set()
    for line in manifest.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        listed.add(relative)
        check(digest(LANE / relative) == expected,
              f"manifest custody {relative}")
    check(listed == MANIFEST_FILES,
          "manifest lists exactly the seven frozen builder files")


print(f"SUMMARY {checks}/{checks} projected-response checks passed")
print("SECTOR states=180 transitions=420 source_momentum=0")
print("RANKS offshell_mod_identity=5 adH=3 ground_retarded=2 M1=2")
print("POLES delta1=2+2sqrt2 delta2=4+2sqrt2; each residue rank1")
print("CEILING FV-WITNESS only; Qdiag(2,4,6)/Ward/RGRLB/gravity/G open")
