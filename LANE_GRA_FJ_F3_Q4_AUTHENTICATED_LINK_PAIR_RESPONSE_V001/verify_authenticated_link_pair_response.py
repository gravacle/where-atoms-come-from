#!/usr/bin/env python3
"""Finite replay for GRA-FJ-F3-Q4-ALPR-V001.

The replay checks the physical four-link Walsh realization, exact conditional
local and cross-cell spectral response, gated nonedge quarantine, and the
first nonzero two-step operator-spreading coefficient.  It does not test a
PMMDC physical-port solder, record qualification of the pair observable,
state preparation, a continuum phase, tensor gravity, or G.
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
THEOREM = (HERE / "THEOREM.md").read_text(encoding="utf-8")

passed = 0


def check(condition: bool, label: str) -> None:
    global passed
    if not condition:
        raise AssertionError(label)
    passed += 1
    print(f"PASS {label}")


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
n_op = (I - Z) / 2


def kron_all(factors: list[np.ndarray]) -> np.ndarray:
    result = np.array([[1.0 + 0.0j]])
    for factor in factors:
        result = np.kron(result, factor)
    return result


def embed(local: np.ndarray, site: int, count: int) -> np.ndarray:
    factors = [I] * count
    factors[site] = local
    return kron_all(factors)


def commutator(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a @ b - b @ a


def spectral_response(
    hamiltonian: np.ndarray,
    a_op: np.ndarray,
    b_op: np.ndarray,
    z_energy: complex,
) -> complex:
    energies, vectors = np.linalg.eigh(hamiltonian)
    ground = vectors[:, 0]
    e0 = energies[0]
    total = 0.0 + 0.0j
    for index in range(1, len(energies)):
        excited = vectors[:, index]
        gap = energies[index] - e0
        gan = np.vdot(ground, a_op @ excited)
        nbg = np.vdot(excited, b_op @ ground)
        gbn = np.vdot(ground, b_op @ excited)
        nag = np.vdot(excited, a_op @ ground)
        total += gan * nbg / (z_energy - gap)
        total -= gbn * nag / (z_energy + gap)
    return total


# Four physical link factors and their six Walsh-pair operators.
edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
z_links = [embed(Z, site, 4) for site in range(4)]
pair_ops = [z_links[a] @ z_links[b] for a, b in edges]

gram = np.array(
    [[np.trace(a.conj().T @ b) / 16 for b in pair_ops] for a in pair_ops]
)
check(np.allclose(gram, np.eye(6)), "six physical Walsh pair operators are linearly independent")
check(all(np.allclose(op, op.conj().T) for op in pair_ops), "all pair operators are self-adjoint")
check(
    all(np.allclose(commutator(a, b), 0) for a in pair_ops for b in pair_ops),
    "pair operators commute at equal time",
)

# S4 permutes the exact six-edge set.
edge_set = set(edges)
s4_closed = True
for permutation in itertools.permutations(range(4)):
    images = {
        tuple(sorted((permutation[a], permutation[b]))) for a, b in edges
    }
    s4_closed &= images == edge_set
check(s4_closed, "all 24 S4 permutations close the six physical pair labels")

# Smallest raw append slab with an interior degree-four child: N=3 and
# c=(1,1,1,1) in S_4.  Its four distinct parents are in S_3, and every one of
# those parents has four append children.
interior_child = (1, 1, 1, 1)
interior_parents = []
for label in range(4):
    parent = list(interior_child)
    parent[label] -= 1
    interior_parents.append(tuple(parent))
check(
    len(set(interior_parents)) == 4
    and all(sum(parent) == 3 for parent in interior_parents)
    and sum(coordinate > 0 for coordinate in interior_child) == 4,
    "N=3 raw slab contains an adjacent degree-four parent/child pair",
)

# Line graph L(K4) and its exact A1/E/T2 spectrum.
line_adjacency = np.zeros((6, 6))
for row, first in enumerate(edges):
    for column, second in enumerate(edges):
        if row != column and set(first).intersection(second):
            line_adjacency[row, column] = 1.0

check(np.allclose(line_adjacency.sum(axis=1), 4), "L(K4) is four-regular")
check(
    np.allclose(np.linalg.eigvalsh(line_adjacency), [-2, -2, 0, 0, 0, 4]),
    "L(K4) spectrum is 4,-2,-2,0,0,0",
)

opposite = np.zeros((6, 6))
for first, second in ((0, 5), (1, 4), (2, 3)):
    opposite[first, second] = 1
    opposite[second, first] = 1
ones = np.ones((6, 1))
p_a = ones @ ones.T / 6
p_t = (np.eye(6) - opposite) / 2
p_e = (np.eye(6) + opposite) / 2 - p_a
check(tuple(np.linalg.matrix_rank(p) for p in (p_a, p_e, p_t)) == (1, 2, 3), "A1/E/T2 ranks are 1,2,3")
check(np.allclose(p_a + p_e + p_t, np.eye(6)), "A1/E/T2 projectors resolve identity")

# Exact admitted independent-link comparator.
h = 0.37
delta = 1.10
kappa = 0.29
z_energy = 1j * kappa
epsilon = np.sqrt(delta**2 + 4 * h**2)
c = delta / epsilon
s = 2 * h / epsilon

check(np.allclose(c**2 + s**2, 1), "single-link rotation is normalized")
check(h > 0 and delta > 0, "response point lies in positive h,Delta domain")

h4 = sum(-h * embed(X, site, 4) + delta * embed(n_op, site, 4) for site in range(4))
energies4, vectors4 = np.linalg.eigh(h4)
ground4 = vectors4[:, 0]
check(np.count_nonzero(np.isclose(energies4, energies4[0], atol=1e-11)) == 1, "four-link comparator ground state is unique")
check(all(np.allclose(np.vdot(ground4, op @ ground4), c) for op in z_links), "ground link polarization equals Delta/epsilon")

comm_expected = 2j * h * (
    embed(Y, 0, 4) @ z_links[1] + z_links[0] @ embed(Y, 1, 4)
)
check(np.allclose(commutator(h4, pair_ops[0]), comm_expected), "pair nonconservation commutator is exact")
check(np.linalg.norm(comm_expected) > 1e-10, "physical pair direction is nonconserved")


def response_factor(energy: float) -> complex:
    return 1 / (z_energy - energy) - 1 / (z_energy + energy)


a_coeff = c**2 * s**2 * response_factor(epsilon)
b_coeff = s**4 * response_factor(2 * epsilon)
analytic = a_coeff * (2 * np.eye(6) + line_adjacency) + b_coeff * np.eye(6)
numeric = np.empty((6, 6), dtype=complex)
for row, a_op in enumerate(pair_ops):
    for column, b_op in enumerate(pair_ops):
        numeric[row, column] = spectral_response(h4, a_op, b_op, z_energy)

check(np.allclose(numeric, analytic, atol=2e-11), "full six-by-six spectral response matches closed form")
check(np.allclose(numeric, numeric.T), "spectral response is reciprocal in comparator")

sector_values = (6 * a_coeff + b_coeff, b_coeff, 2 * a_coeff + b_coeff)
for projector, value, name in zip((p_a, p_e, p_t), sector_values, ("A1", "E", "T2")):
    check(np.allclose(numeric @ projector, value * projector, atol=2e-11), f"{name} response eigenvalue is exact")
    check(abs(value) > 1e-10, f"{name} response is nonzero at imaginary energy")

# Adjacent cells share one literal physical link; disjoint pairs do not.
h3 = sum(-h * embed(X, site, 3) + delta * embed(n_op, site, 3) for site in range(3))
z3 = [embed(Z, site, 3) for site in range(3)]
pair_left = z3[0] @ z3[1]
pair_right = z3[0] @ z3[2]
cross = spectral_response(h3, pair_left, pair_right, z_energy)
check(np.allclose(cross, a_coeff, atol=2e-11), "shared-link adjacent-cell response equals a(z)")
check(abs(cross) > 1e-10, "shared-link adjacent-cell response is nonzero")
disjoint = spectral_response(h4, pair_ops[0], pair_ops[5], z_energy)
check(np.allclose(disjoint, 0, atol=2e-11), "disjoint-cell pair response vanishes in independent-link comparator")

# The K-gated actuator protects a blank nonedge; the raw actuator would not.
q_k = np.diag([0.0, 1.0])
h_gated = -h * np.kron(q_k, X) + delta * np.kron(I, n_op)
blank_nonedge = np.array([1, 0, 0, 0], dtype=complex)
check(np.allclose(h_gated @ blank_nonedge, 0), "K=0,n=0 nonedge is invariant under gated response block")
h_raw = -h * np.kron(I, X) + delta * np.kron(I, n_op)
check(np.linalg.norm(h_raw @ blank_nonedge) > 1e-10, "raw ungated X would violate nonedge quarantine")

# In the K=0 block the link Hamiltonian is diagonal, so Z-pair response is
# exactly zero.  This is an operator-block statement, not a matched physical
# KEEP/BREAK experiment.
h_diag2 = sum(delta * embed(n_op, site, 2) for site in range(2))
z_diag2 = embed(Z, 0, 2) @ embed(Z, 1, 2)
check(
    np.allclose(spectral_response(h_diag2, z_diag2, z_diag2, z_energy), 0, atol=2e-11),
    "K=0 diagonal block has zero Z-pair spectral response",
)

# Exact fifth-order operator-spreading coefficient on a two-step line-graph path.
u_degree = 0.46
j_coupling = u_degree / 2
h_chain = sum(-h * embed(X, site, 3) + delta * embed(n_op, site, 3) for site in range(3))
h_chain += j_coupling * (z3[0] @ z3[1] + z3[1] @ z3[2])

ad = z3[0]
for order in range(1, 5):
    ad = commutator(h_chain, ad)
    check(np.allclose(commutator(ad, z3[2]), 0, atol=1e-10), f"no two-step response below nested order {order + 1}")
ad = commutator(h_chain, ad)
fifth_cross = commutator(ad, z3[2])
xxx = embed(X, 0, 3) @ embed(X, 1, 3) @ embed(X, 2, 3)
xxx_coefficient = np.trace(xxx.conj().T @ fifth_cross) / 8
check(np.allclose(xxx_coefficient, -64 * h**3 * j_coupling**2, atol=2e-10), "fifth-order spreading coefficient is -64 h^3 J^2")
check(np.linalg.norm(fifth_cross) > 1e-10, "degree interaction produces nonzero two-step operator spreading")

# The minimal fifth-order Pauli coefficient is independent of one-link
# diagonal fields because its five commutators are already exhausted by three
# transverse flips and two path couplings.
delta_alt = 0.23
h_chain_alt = sum(
    -h * embed(X, site, 3) + delta_alt * embed(n_op, site, 3)
    for site in range(3)
)
h_chain_alt += j_coupling * (z3[0] @ z3[1] + z3[1] @ z3[2])
ad_alt = z3[0]
for _ in range(5):
    ad_alt = commutator(h_chain_alt, ad_alt)
fifth_cross_alt = commutator(ad_alt, z3[2])
xxx_coefficient_alt = np.trace(xxx.conj().T @ fifth_cross_alt) / 8
check(
    np.allclose(xxx_coefficient_alt, -64 * h**3 * j_coupling**2, atol=2e-10),
    "fifth-order coefficient is independent of one-link diagonal detuning",
)

# Claim and type ceilings must remain explicit.
required_phrases = (
    "These types must not be collapsed",
    "This is not a new interaction",
    "unlike a mere equality of six dimensions",
    "does **not** identify the four-link query with PMMDC's physical four-port",
    "not the full PMMDC/Q4 physical solder demanded by `PFLD`",
    "For the adjacent-cell and spreading statements choose `N>=3`",
    "does not prove that the supplied FPSS address program",
    "supplied stationary ground state",
    "not yet a matched KEEP/BREAK response intervention",
    "not yet a propagating continuum cone",
    "does not prove a massless pole",
    "complete the Gravity Formation proof",
)
normalized_theorem = " ".join(THEOREM.split())
for phrase in required_phrases:
    normalized_phrase = " ".join(phrase.split())
    check(normalized_phrase in normalized_theorem, f"claim/type ceiling: {phrase}")

check(
    sum(line.strip() == r"\[" for line in THEOREM.splitlines())
    == sum(line.strip() == r"\]" for line in THEOREM.splitlines()),
    "display-math delimiters are balanced",
)

print(f"SUMMARY {passed}/{passed} PASS")
print("VERDICT FPMH_SUPPORTED_Q4_WALSH_OPERATOR_AND_CONDITIONAL_FINITE_RESPONSE_PASS__PMMDC_SOLDER_RECORD_QUALIFICATION_STATE_PREPARATION_MATCHED_BREAK_PORT_PHASE_CONTINUUM_AND_GRAVITY_OPEN")
