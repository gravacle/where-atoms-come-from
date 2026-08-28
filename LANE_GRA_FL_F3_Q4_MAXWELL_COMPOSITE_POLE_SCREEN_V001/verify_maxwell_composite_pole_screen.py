#!/usr/bin/env python3
"""Replay GRA-FL-F3-Q4-MCPS-V001.

The finite replay verifies the exact ice/S4/parity identities, the vector and
TT projector algebra, a Fock-space one- versus two-particle selection-rule
witness, and the massless two-particle threshold.  It also freezes the FJ and
CROSS-CW dependencies.  The thermodynamic Maxwell phase itself is an
explicit imported premise and is not numerically re-proved by this script.
"""

from __future__ import annotations

import hashlib
import itertools
import re
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
THEOREM = (HERE / "THEOREM.md").read_text(encoding="utf-8")
THEOREM_COMPACT = " ".join(THEOREM.split())

passed = 0


def check(condition: bool, label: str) -> None:
    global passed
    if not condition:
        raise AssertionError(label)
    passed += 1
    print(f"PASS {label}")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# Frozen theorem and evidence custody.
dependencies = {
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md":
        "05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769",
    "LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "44690ad431c85af7a4947a431a4c57ad4cd8b19a346e3d26720b341f77256f90",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/PRIMARY_SOURCES.md":
        "4ee84b4f9b78003cdc5ce80a86cba6cbab618feb1fcd78d25903cb5e97c42a62",
}
for relative, expected in dependencies.items():
    path = ROOT / relative
    check(path.is_file(), f"dependency exists: {relative}")
    check(sha256(path) == expected, f"dependency hash frozen: {relative}")
    check(expected in THEOREM, f"dependency hash printed in theorem: {relative}")


# Enumerate the exact two-in/two-out local fiber.
states = np.array(
    [s for s in itertools.product((-1, 1), repeat=4) if sum(s) == 0],
    dtype=int,
)
edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
pairs = np.array(
    [[s[a] * s[b] for a, b in edges] for s in states],
    dtype=int,
)

check(states.shape == (6, 4), "ice fiber has exactly six states")
check(np.all(states.sum(axis=1) == 0), "every enumerated state obeys the ice rule")
check(np.linalg.matrix_rank(states) == 3, "nonconstant one-link span has rank three")
check(np.linalg.matrix_rank(pairs) == 3, "pair span including its constant has rank three")

center = np.eye(6) - np.ones((6, 6)) / 6
centered_pairs = center @ pairs
check(np.linalg.matrix_rank(centered_pairs) == 2, "connected pair span has rank two")
check(np.all(pairs.sum(axis=1) == -2), "pair A1 sum is the fixed value minus two")

opposite_indices = ((0, 5), (1, 4), (2, 3))
check(
    all(np.array_equal(pairs[:, a], pairs[:, b]) for a, b in opposite_indices),
    "all three opposite-edge pair identities hold on ice",
)

state_index = {tuple(s): index for index, s in enumerate(states)}
complement = np.array([state_index[tuple(-s)] for s in states])
check(np.array_equal(states[complement], -states), "complement makes every one-link value odd")
check(np.array_equal(pairs[complement], pairs), "complement makes every pair value even")


# Exact A1/E/T2 projectors for the six tetrahedral edges.
opposite = np.zeros((6, 6))
for a, b in opposite_indices:
    opposite[a, b] = opposite[b, a] = 1
ones6 = np.ones((6, 1))
p_a1 = ones6 @ ones6.T / 6
p_t2 = (np.eye(6) - opposite) / 2
p_e = (np.eye(6) + opposite) / 2 - p_a1

check(
    tuple(np.linalg.matrix_rank(p) for p in (p_a1, p_e, p_t2)) == (1, 2, 3),
    "edge A1/E/T2 projector ranks are one, two, three",
)
check(np.allclose(p_a1 + p_e + p_t2, np.eye(6)), "edge projectors resolve identity")
check(
    all(np.allclose(p @ p, p) for p in (p_a1, p_e, p_t2)),
    "edge projectors are idempotent",
)
check(
    all(np.allclose(p @ q, 0) for p, q in ((p_a1, p_e), (p_a1, p_t2), (p_e, p_t2))),
    "edge projectors are mutually orthogonal",
)
check(np.allclose(pairs @ p_t2, 0), "ice projection kills pair T2 exactly")
check(np.allclose(centered_pairs @ p_a1, 0), "connected response kills constant pair A1")
check(np.allclose(centered_pairs @ p_e, centered_pairs), "all connected pair variation lies in E")


# All S4 permutations preserve the projectors.  The vertex sum-zero module and
# pair T2 module have the same character, fixing the convention called T2.
p_vertex_t2 = np.eye(4) - np.ones((4, 4)) / 4
character_match = True
commutes = True
for permutation in itertools.permutations(range(4)):
    r4 = np.zeros((4, 4))
    for old in range(4):
        r4[permutation[old], old] = 1
    r6 = np.zeros((6, 6))
    for old, (a, b) in enumerate(edges):
        image = tuple(sorted((permutation[a], permutation[b])))
        r6[edges.index(image), old] = 1
    commutes &= np.allclose(r4 @ p_vertex_t2, p_vertex_t2 @ r4)
    commutes &= all(np.allclose(r6 @ p, p @ r6) for p in (p_a1, p_e, p_t2))
    character_match &= np.isclose(np.trace(p_vertex_t2 @ r4), np.trace(p_t2 @ r6))

check(commutes, "all 24 S4 actions preserve the stated sector projectors")
check(character_match, "one-link sum-zero and pair T2 modules have identical S4 characters")


# Exact tetrahedral frame and lattice electric-flux reconstruction.
frame = np.array(
    [
        (1, 1, 1),
        (1, -1, -1),
        (-1, 1, -1),
        (-1, -1, 1),
    ],
    dtype=float,
) / np.sqrt(3)

expected_gram = np.full((4, 4), -1 / 3)
np.fill_diagonal(expected_gram, 1)
check(np.allclose(frame.sum(axis=0), 0), "tetrahedral bond vectors sum to zero")
check(np.allclose(frame @ frame.T, expected_gram), "tetrahedral Gram matrix is exact")
check(np.allclose(frame.T @ frame, 4 * np.eye(3) / 3), "tetrahedral frame is complete")

fluxes = (3 / 4) * states @ frame
reconstructed = fluxes @ frame.T
check(np.allclose(reconstructed, states), "oriented epsilon equals e dot local electric flux on every ice state")
check(np.linalg.matrix_rank(fluxes) == 3, "local electric-flux image spans three dimensions")

# Across a shared edge, both the outward frame and incidence-corrected scalar
# reverse, so their vector product is endpoint independent.
epsilon_plus = states
epsilon_minus = -states
frame_plus = frame
frame_minus = -frame
check(
    np.allclose(
        epsilon_plus[:, :, None] * frame_plus,
        epsilon_minus[:, :, None] * frame_minus,
    ),
    "incidence-corrected flux glues consistently across the two sublattices",
)

# The full tetrahedral point group T_d is isomorphic to S4 and contains
# improper operations.  Its polar-vector representation is the vertex T2,
# while polar symmetric-traceless rank two restricts as E+T2.  This is a
# finite-group character statement, not an identification of continuum spin.
polar_actions_match = True
rank_two_characters_match = True
determinants_match_permutation_parity = True
for permutation in itertools.permutations(range(4)):
    r4 = np.zeros((4, 4))
    for old in range(4):
        r4[permutation[old], old] = 1
    r6 = np.zeros((6, 6))
    for old, (a, b) in enumerate(edges):
        image = tuple(sorted((permutation[a], permutation[b])))
        r6[edges.index(image), old] = 1

    # R e_a = e_{permutation(a)} for the polar O(3) realization of T_d.
    r3 = (3 / 4) * frame.T @ r4 @ frame
    polar_actions_match &= np.allclose((r3 @ frame.T).T, frame[list(permutation)])

    # Character of Sym^2(R) with the scalar trace removed.
    chi_l2 = 0.5 * (np.trace(r3) ** 2 + np.trace(r3 @ r3)) - 1
    chi_e_plus_t2 = np.trace((p_e + p_t2) @ r6)
    rank_two_characters_match &= np.isclose(chi_l2, chi_e_plus_t2)

    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    permutation_sign = -1 if inversions % 2 else 1
    determinants_match_permutation_parity &= np.isclose(np.linalg.det(r3), permutation_sign)

check(polar_actions_match, "all 24 S4 permutations realize the full tetrahedral polar O(3) action")
check(rank_two_characters_match, "polar symmetric-traceless rank two restricts as E plus T2")
check(
    determinants_match_permutation_parity,
    "odd S4 permutations are the improper classes of full T_d",
)


# Maxwell vector projectors and the complete four-link photon residue.
momenta = (
    np.array((1.0, 2.0, 3.0)),
    np.array((-2.0, 0.5, 1.25)),
    np.array((0.0, 0.0, 2.0)),
)
for number, momentum in enumerate(momenta, start=1):
    unit = momentum / np.linalg.norm(momentum)
    transverse = np.eye(3) - np.outer(unit, unit)
    check(np.allclose(transverse, transverse.T), f"sample {number} transverse vector projector is symmetric")
    check(np.allclose(transverse @ transverse, transverse), f"sample {number} transverse vector projector is idempotent")
    check(np.linalg.matrix_rank(transverse) == 2, f"sample {number} vector photon space has rank two")
    check(np.allclose(transverse @ unit, 0), f"sample {number} longitudinal vector is removed")
    link_residue = frame @ transverse @ frame.T
    check(np.linalg.matrix_rank(link_residue) == 2, f"sample {number} complete link residue has rank two")
    check(np.all(np.linalg.eigvalsh(link_residue) > -1e-11), f"sample {number} link residue is positive semidefinite")

# Canonical Maxwell normalization: for one transverse oscillator,
# |<0|E|1>|^2 = omega_k/(2 chi), exactly the positive-frequency residue of
# c^2 k^2/[chi(z^2-c^2 k^2)].  The omega^2 numerator differs by a contact 1.
chi = 1.7
c_speed = 0.83
k_abs = 0.41
omega_k = c_speed * k_abs
electric_matrix_element_sq = omega_k / (2 * chi)
electric_pole_residue = (c_speed**2 * k_abs**2 / chi) / (2 * omega_k)
check(
    np.isclose(electric_pole_residue, electric_matrix_element_sq),
    "electric-field positive-frequency pole residue is linear in momentum",
)
complex_frequencies = (0.2 + 0.3j, 1.1 + 0.4j)
check(
    all(
        np.isclose(
            z**2 / (z**2 - omega_k**2)
            - omega_k**2 / (z**2 - omega_k**2),
            1,
        )
        for z in complex_frequencies
    ),
    "omega-squared and k-squared Maxwell numerators differ by an analytic contact",
)


# A genuine continuum helicity-two pole has a TT residue.  Its rotations about
# k carry angle 2 theta, unlike the vector pole's angle theta.
unit = np.array((0.0, 0.0, 1.0))
transverse = np.eye(3) - np.outer(unit, unit)
tt = np.empty((3, 3, 3, 3))
for i, j, k, l in itertools.product(range(3), repeat=4):
    tt[i, j, k, l] = 0.5 * (
        transverse[i, k] * transverse[j, l]
        + transverse[i, l] * transverse[j, k]
        - transverse[i, j] * transverse[k, l]
    )
tt_matrix = tt.reshape(9, 9)
check(np.allclose(tt_matrix, tt_matrix.T), "TT projector is symmetric")
check(np.allclose(tt_matrix @ tt_matrix, tt_matrix), "TT projector is idempotent")
check(np.linalg.matrix_rank(tt_matrix) == 2, "TT helicity space has rank two")
check(np.allclose(np.einsum("iikl->kl", tt), 0), "TT projector removes the tensor trace")
check(np.allclose(np.einsum("i,ijkl->jkl", unit, tt), 0), "TT projector is transverse")

theta = np.pi / 5
rotation = np.array(
    [
        (np.cos(theta), -np.sin(theta), 0),
        (np.sin(theta), np.cos(theta), 0),
        (0, 0, 1),
    ]
)
tensor_rotation = np.empty((3, 3, 3, 3))
for i, j, k, l in itertools.product(range(3), repeat=4):
    tensor_rotation[i, j, k, l] = rotation[i, k] * rotation[j, l]
vector_character = np.trace(transverse @ rotation)
tensor_character = np.trace(tt_matrix @ tensor_rotation.reshape(9, 9))
check(np.isclose(vector_character, 2 * np.cos(theta)), "vector pole transforms with helicity angle theta")
check(np.isclose(tensor_character, 2 * np.cos(2 * theta)), "TT pole transforms with helicity angle two theta")
check(not np.isclose(vector_character, tensor_character), "spin-one and helicity-two residues are not the same representation")


# One oscillator mode supplies an exact selection-rule witness: an odd field
# reaches one particle, while its normal-ordered square reaches two but not one.
dimension = 6
annihilation = np.zeros((dimension, dimension))
for n in range(1, dimension):
    annihilation[n - 1, n] = np.sqrt(n)
creation = annihilation.T
field = annihilation + creation
normal_square = annihilation @ annihilation + creation @ creation + 2 * creation @ annihilation
number_parity = np.diag([(-1) ** n for n in range(dimension)])
vacuum = np.eye(dimension)[:, 0]
one_particle = np.eye(dimension)[:, 1]
two_particle = np.eye(dimension)[:, 2]

check(np.allclose(number_parity @ field @ number_parity, -field), "linear Maxwell field is complement odd")
check(np.allclose(number_parity @ normal_square @ number_parity, normal_square), "normal-ordered bilinear is complement even")
check(not np.isclose(vacuum @ field @ one_particle, 0), "odd field has a nonzero one-particle matrix element")
check(np.isclose(vacuum @ normal_square @ one_particle, 0), "even bilinear has zero one-particle matrix element")
check(not np.isclose(vacuum @ normal_square @ two_particle, 0), "even bilinear has a nonzero two-particle matrix element")


# Exact kinematic threshold for two linearly dispersing massless particles.
momentum = np.array((1.2, -0.5, 0.7))
momentum_norm = np.linalg.norm(momentum)
rng = np.random.default_rng(20260827)
samples = rng.normal(size=(2048, 3))
energies = np.linalg.norm(samples, axis=1) + np.linalg.norm(momentum - samples, axis=1)
check(np.all(energies >= momentum_norm - 2e-14), "sampled two-photon energies obey the triangle threshold")

unit = momentum / momentum_norm
t_values = np.linspace(0, 6, 121)
continuum_energies = np.array(
    [np.linalg.norm(-t * unit) + np.linalg.norm(momentum + t * unit) for t in t_values]
)
check(np.allclose(continuum_energies, momentum_norm + 2 * t_values), "two-photon kinematics continuously fills energies above threshold")
check(np.isclose(continuum_energies[0], momentum_norm), "two-photon continuum begins on the Maxwell light cone")

# No universal threshold-suppression exponent follows merely from using field
# strengths.  At collinear threshold, two x-polarized photons moving along z
# have a nonzero vacuum-to-two-particle matrix element for :E_x E_x:.
split = 0.37
first_energy = split * momentum_norm
second_energy = (1 - split) * momentum_norm
collinear_electric_bilinear = np.sqrt(first_energy * second_energy) / (2 * chi)
check(
    collinear_electric_bilinear > 0,
    "field-strength bilinear onset is operator dependent, not universally power suppressed",
)


# Claim boundaries and document integrity.
required_phrases = (
    "exact conditional deduction",
    "zero-charge lattice Gauss law with correct incidence signs",
    "primary-source custody is arXiv v3",
    "not an internal finite-F3 theorem",
    "Source-free Gaussian Maxwell duality exchanges those names",
    "one linearly dispersing transverse spin-one photon",
    "positive-frequency complex pole residue",
    "no one-photon pole, by (FL16)",
    "two-photon continuum beginning at `|omega|=c|k|`",
    "no universal threshold exponent is asserted",
    "not a helicity-two pole",
    "T_d \\simeq S_4",
    "proper-rotation subgroup alone is `A_4`",
    "Finite-group labels therefore do not determine continuum spin",
    "This is a no-go for the **direct Gaussian composite route**",
    "complete F3 all-orders phase and volume-uniform stability",
    "No `j-j` attraction or fitted tensor interaction may be inserted",
    "rank-two Gauss/diffeomorphism-type constraint and Ward identity",
    "they may not be renamed as the graviton",
    "separately owned local ice-sector theorem is a consistency cross-check",
    "dependency hash will be pinned only after its hostile audit is final",
    "do not form a circular proof chain",
)
for phrase in required_phrases:
    check(phrase in THEOREM_COMPACT, f"required claim boundary present: {phrase}")

for forbidden in (
    "the F3 photon is visible electromagnetism",
    "the pair observable is the graviton",
    "unconditional helicity-two pole",
    "MAXWELL-IR is an internal theorem",
    "universal threshold suppression",
):
    check(forbidden not in THEOREM, f"forbidden promotion absent: {forbidden}")

tags = re.findall(r"\\tag\{(FL\d{2})\}", THEOREM)
check(tags == [f"FL{number:02d}" for number in range(1, 24)], "equation tags FL01-FL23 are unique and ordered")
check(THEOREM.count(r"\[") == THEOREM.count(r"\]"), "display-math delimiters are balanced")
check(THEOREM.count("`") % 2 == 0, "Markdown backticks are balanced")
check(not any(line.endswith((" ", "\t")) for line in THEOREM.splitlines()), "theorem has no trailing whitespace")
check(not any(ord(char) < 32 and char not in "\n" for char in THEOREM), "theorem has no forbidden control characters")

print()
print(f"Maxwell composite pole screen verification: {passed} passed, 0 failed")
print("Exact sector ranks: one-link=3, pair-total=3, pair-connected=2, pair-T2=0")
print("Conditional IR spectrum: spin-1 pole in odd link channel; two-photon continuum in even pair channel")
print("Helicity-2 isolated pole: absent at the Gaussian Maxwell fixed point")
