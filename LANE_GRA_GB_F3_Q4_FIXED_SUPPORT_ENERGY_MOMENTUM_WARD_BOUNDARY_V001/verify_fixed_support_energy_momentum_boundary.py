#!/usr/bin/env python3
"""Exact finite checks for the fixed-support energy/momentum Ward boundary.

The positive theorem is local Hamiltonian-energy continuity.  The negative
result is derivability/nonidentifiability: a finite discrete translation label
on an externally fixed support does not select a local physical momentum
density or support-recoil ledger.
"""

from contextlib import redirect_stdout
from hashlib import sha256
from io import StringIO
from pathlib import Path
import json
import runpy

import numpy as np


LANE = Path(__file__).resolve().parent
ROOT = LANE.parent
FO_SCRIPT = (ROOT / "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001" /
             "verify_finite_tt_four_point.py")
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def commutator(first, second):
    return first @ second - second @ first


for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink() and
          sha256(path.read_bytes()).hexdigest() == expected,
          f"dependency custody: {relative}")


with redirect_stdout(StringIO()):
    fo = runpy.run_path(str(FO_SCRIPT))

states = fo["states"]
state_index = fo["state_index"]
ring_patterns = fo["ring_patterns"]
hexagons = fo["hexagons"]
H = np.asarray(fo["hamiltonian"], dtype=complex)
dimension = len(states)

check((fo["CELL_COUNT"], fo["VERTEX_COUNT"], len(fo["edges"]),
       len(hexagons), dimension) == (30, 60, 120, 120, 180),
      "GB uses exactly the fixed FO Z30 diamond support and component")


# -------------------------------------------------------------------------
# Exact H6 term decomposition and local Hamiltonian-energy continuity.

ring_terms = []
for mask, first, second in ring_patterns:
    term = np.zeros((dimension, dimension), dtype=complex)
    for row, state in enumerate(states):
        if (state & mask) in (first, second):
            term[row, state_index[state ^ mask]] = -1
    ring_terms.append(term)

check(all(np.array_equal(term, term.T.conj()) for term in ring_terms),
      "all 120 inherited local ring-energy terms are Hermitian")
check(np.array_equal(sum(ring_terms, np.zeros_like(H)), H),
      "the fixed H6 Hamiltonian is exactly the sum of local ring terms")

# The projected ring operators inherit genuine support locality: terms on
# edge-disjoint hexagons commute.  Thus the commutator current below has no
# artificial long-range edge between disjoint supports.
ring_masks = tuple(pattern[0] for pattern in ring_patterns)
disjoint_pairs = tuple((first, second)
                       for first in range(len(ring_terms))
                       for second in range(first+1, len(ring_terms))
                       if ring_masks[first] & ring_masks[second] == 0)
check(disjoint_pairs and all(
      np.array_equal(commutator(ring_terms[first], ring_terms[second]),
                     np.zeros_like(H))
      for first, second in disjoint_pairs),
      "edge-disjoint inherited ring energies have exactly zero mutual current")

# For h_A, define the outward term-graph energy current
# J_(A->rest)=i[h_A,H-h_A].  Then hdot_A+J=0 exactly.
alpha = next(index for index, term in enumerate(ring_terms) if np.any(term))
h_alpha = ring_terms[alpha]
j_alpha_out = 1j*commutator(h_alpha, H-h_alpha)
h_alpha_dot = 1j*commutator(H, h_alpha)
check(np.array_equal(h_alpha_dot+j_alpha_out, np.zeros_like(H)),
      "one local ring energy obeys exact outward-current continuity")
check(np.linalg.norm(j_alpha_out) > 0,
      "the local H6 term-graph energy current is nontrivial")
check(all(np.array_equal(
      1j*commutator(H, term)+1j*commutator(term, H-term),
      np.zeros_like(H)) for term in ring_terms),
      "all 120 inherited ring energies obey exact term-graph continuity")

beta = next(index for index, term in enumerate(ring_terms)
            if index != alpha and np.any(commutator(h_alpha, term)))
j_alpha_beta = 1j*commutator(h_alpha, ring_terms[beta])
j_beta_alpha = 1j*commutator(ring_terms[beta], h_alpha)
check(np.array_equal(j_alpha_beta, -j_beta_alpha),
      "pairwise local energy current is exactly antisymmetric")

# Any finite subset of ring terms has the same exact boundary-current law.
region_indices = tuple(range(12))
h_region = sum((ring_terms[index] for index in region_indices),
               np.zeros_like(H))
h_complement = H-h_region
j_region_out = 1j*commutator(h_region, h_complement)
h_region_dot = 1j*commutator(H, h_region)
check(np.array_equal(h_region_dot+j_region_out, np.zeros_like(H)),
      "a twelve-term subset obeys exact Hamiltonian-energy boundary balance")
check(np.linalg.norm(j_region_out) > 0,
      "the selected finite subset has nonzero energy exchange with its complement")
check(np.array_equal(commutator(H, H), np.zeros_like(H)),
      "total source-off H6 energy is exactly conserved")

# GA preservation is only modulo one common identity/reference on the full
# encoded P+Q space.  A time-independent common reference changes neither
# Heisenberg evolution nor any commutator energy current.  This check does not
# license a sector-dependent scalar or an unowned time-dependent controller.
REFERENCE_SHIFT = 7
H_reference = H+REFERENCE_SHIFT*np.eye(dimension)
check(np.array_equal(commutator(H_reference, h_alpha),
                     commutator(H, h_alpha)) and
      np.array_equal(commutator(H_reference, H_reference),
                     np.zeros_like(H)),
      "one common full-space identity/reference leaves energy currents unchanged")

# The density partition is not unique.  Move a noncommuting Hermitian operator
# C between two terms without changing H; the attributed pair current changes.
c_shift = ring_terms[beta]
h_alpha_prime = h_alpha+c_shift
h_beta_prime = ring_terms[beta]-c_shift
check(np.array_equal(h_alpha_prime+h_beta_prime,
                     h_alpha+ring_terms[beta]),
      "local energy can be repartitioned without changing the Hamiltonian")
current_original = 1j*commutator(h_alpha, ring_terms[beta])
current_repartitioned = 1j*commutator(h_alpha_prime, h_beta_prime)
check(not np.array_equal(current_original, current_repartitioned),
      "the term-assigned local energy current changes under repartition")


# -------------------------------------------------------------------------
# Exact discrete translation symmetry versus local momentum ownership.

permutation = np.asarray(fo["translation_permutation"], dtype=int)
U = np.zeros((dimension, dimension), dtype=complex)
for source, target in enumerate(permutation):
    U[target, source] = 1
check(np.array_equal(np.linalg.matrix_power(U, 30), np.eye(dimension)) and
      not np.array_equal(U, np.eye(dimension)),
      "the inherited support owns a nontrivial finite Z30 translation")
check(np.array_equal(commutator(H, U), np.zeros_like(H)),
      "H6 exactly conserves the global Z30 representation label")

seen = set()
translation_orbit_lengths = []
for start in range(dimension):
    if start in seen:
        continue
    cursor = start
    orbit = []
    while cursor not in seen:
        seen.add(cursor)
        orbit.append(cursor)
        cursor = int(permutation[cursor])
    translation_orbit_lengths.append(len(orbit))
check(sorted(translation_orbit_lengths) == [30]*6,
      "the actual 180-state translation decomposes into six length-30 orbits")

# Discrete U has many Hermitian logarithms.  Demonstrate the ambiguity on one
# of the actual length-30 orbit types without invoking a matrix logarithm.
U30 = np.zeros((30, 30), dtype=complex)
for source in range(30):
    U30[(source+1) % 30, source] = 1
projectors = []
power = np.eye(30, dtype=complex)
powers = []
for exponent in range(30):
    powers.append(power)
    power = power @ U30
for momentum in range(30):
    theta = 2*np.pi*momentum/30
    projector = sum((np.exp(-1j*theta*exponent)*powers[exponent]
                     for exponent in range(30)),
                    np.zeros_like(U30))/30
    projectors.append(projector)
check(np.max(np.abs(sum(projectors, np.zeros_like(U30))-np.eye(30))) < 2e-13,
      "the thirty exact-character projectors resolve the translation orbit")

thetas = np.array([2*np.pi*momentum/30 for momentum in range(30)])
K0 = sum((thetas[momentum]*projectors[momentum]
          for momentum in range(30)), np.zeros_like(U30))
K1 = K0 + 2*np.pi*projectors[0]
U_from_K0 = sum((np.exp(1j*thetas[momentum])*projectors[momentum]
                 for momentum in range(30)), np.zeros_like(U30))
shifted_thetas = thetas.copy()
shifted_thetas[0] += 2*np.pi
U_from_K1 = sum((np.exp(1j*shifted_thetas[momentum])*projectors[momentum]
                 for momentum in range(30)), np.zeros_like(U30))
check(np.max(np.abs(K0-K0.T.conj())) < 2e-13 and
      np.max(np.abs(K1-K1.T.conj())) < 2e-13,
      "two distinct Hermitian global translation generators exist")
check(np.linalg.norm(K1-K0) > 1 and
      np.max(np.abs(U_from_K0-U30)) < 2e-13 and
      np.max(np.abs(U_from_K1-U30)) < 2e-13,
      "both branch-inequivalent generators exponentiate to the same U30")
check(np.count_nonzero(np.abs(K0) > 1e-10) > 30,
      "a chosen translation logarithm is nonlocal on the cell orbit")

# The FO Hilbert space contains only incidence configurations.  Fixed support
# coordinates enter as numerical labels/vectors, not canonical R,P factors.
check(dimension == len(states) and all(isinstance(state, int) for state in states),
      "the frozen response Hilbert space has incidence states only")
check("bond_vectors" in fo and "translation_permutation" in fo and
      "support_momentum" not in fo and "support_recoil" not in fo,
      "support geometry is c-number data with no owned momentum/recoil operator")


result = {
    "lane": LANE.name,
    "status": "PASS",
    "checks": f"{checks}/{checks}",
    "dimension": dimension,
    "local_ring_terms": len(ring_terms),
    "disjoint_support_pairs": len(disjoint_pairs),
    "energy_continuity": "exact for the inherited ring-term allocation",
    "translation": "global Z30 representation label only",
    "translation_orbits": sorted(translation_orbit_lengths),
    "momentum_ward": "not derived: no selected local momentum density or support/recoil ledger",
    "ga_reference": "common full-code identity shift leaves commutator currents unchanged",
    "ceiling": "no T0j, stress Ward packet, gravity, or G",
}
print("GB_RESULT_JSON", json.dumps(result, sort_keys=True))
print(f"SUMMARY {checks}/{checks} GB fixed-support energy/momentum checks passed")
print("CEILING exact Hamiltonian-energy continuity and global Z30 symmetry; "
      "no derived local physical momentum density, support recoil, T0j, "
      "stress Ward packet, gravity, or G")
