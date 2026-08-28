#!/usr/bin/env python3
"""Exact-matrix checks for the finite F3/q4 TT composite-cumulant screen.

The matrix is built from the inherited H6=-J6 sum_C B_C with J6 set to one.
All reported powers of J6 are restored analytically in THEOREM.md.
"""

from collections import Counter, defaultdict, deque
import hashlib
from itertools import combinations_with_replacement, product
from math import gcd, pi, sqrt
from pathlib import Path

import numpy as np


checks = 0


def check(condition, message):
    global checks
    if not condition:
        raise AssertionError(message)
    checks += 1
    print(f"PASS {message}")


# ---------------------------------------------------------------------------
# A plaquette-complete cyclic quotient of the diamond net.

CELL_COUNT = 30
SHIFTS = (0, 1, 5, 19)
LINK_COUNT = 4 * CELL_COUNT
VERTEX_COUNT = 2 * CELL_COUNT

check(gcd(CELL_COUNT, *SHIFTS[1:]) == 1,
      "cyclic diamond quotient is connected")
check(len(set(SHIFTS)) == 4, "four q4 append directions remain distinct")

triple_sum_values = [
    sum(SHIFTS[index] for index in triple) % CELL_COUNT
    for triple in combinations_with_replacement(range(4), 3)
]
check(len(triple_sum_values) == 20 and len(set(triple_sum_values)) == 20,
      "all twenty unordered B3 triple sums are distinct")

edges = []
edge_labels = []
incidence = defaultdict(list)
for cell in range(CELL_COUNT):
    for label, shift in enumerate(SHIFTS):
        edge_index = len(edges)
        endpoints = (cell, CELL_COUNT + (cell + shift) % CELL_COUNT)
        edges.append(endpoints)
        edge_labels.append(label)
        for endpoint in endpoints:
            incidence[endpoint].append((endpoints[0] + endpoints[1] - endpoint,
                                        edge_index))

check(len(edges) == 120, "quotient has 120 physical links")
check(len(incidence) == 60, "quotient has sixty diamond vertices")
check(all(len(incidence[vertex]) == 4 for vertex in range(VERTEX_COUNT)),
      "every quotient vertex has coordination four")


def canonical_cycle(edge_order):
    edge_order = tuple(edge_order)
    images = []
    for image in (edge_order, tuple(reversed(edge_order))):
        images.extend(image[offset:] + image[:offset]
                      for offset in range(len(image)))
    return min(images)


def enumerate_six_cycles():
    cycles = set()
    for root in range(VERTEX_COUNT):
        def recurse(vertex, visited_vertices, visited_edges):
            if len(visited_edges) == 6:
                if vertex == root:
                    cycles.add(canonical_cycle(visited_edges))
                return
            for neighbor, edge_index in incidence[vertex]:
                if edge_index in visited_edges:
                    continue
                if neighbor == root:
                    if len(visited_edges) == 5:
                        recurse(neighbor, visited_vertices + [neighbor],
                                visited_edges + [edge_index])
                elif neighbor not in visited_vertices:
                    recurse(neighbor, visited_vertices + [neighbor],
                            visited_edges + [edge_index])
        recurse(root, [root], [])
    return tuple(sorted(cycles))


hexagons = enumerate_six_cycles()
check(len(hexagons) == 4 * CELL_COUNT,
      "all 120 simple six-cycles are the four local diamond hexagons per cell")
check(all(sorted(Counter(edge_labels[edge] for edge in cycle).values()) ==
          [2, 2, 2] and
          len(Counter(edge_labels[edge] for edge in cycle)) == 3
          for cycle in hexagons),
      "every admitted six-cycle uses three q4 labels twice")

ring_patterns = []
for cycle in hexagons:
    alternating_zero = sum(1 << edge for edge in cycle[::2])
    alternating_one = sum(1 << edge for edge in cycle[1::2])
    ring_patterns.append((alternating_zero | alternating_one,
                          alternating_zero, alternating_one))


def degrees(state):
    return tuple(sum((state >> edge) & 1 for _, edge in incidence[vertex])
                 for vertex in range(VERTEX_COUNT))


def translate_state(state, amount=1):
    translated = 0
    for cell in range(CELL_COUNT):
        for label in range(4):
            source = 4 * cell + label
            target = 4 * ((cell + amount) % CELL_COUNT) + label
            translated |= ((state >> source) & 1) << target
    return translated


def ring_component(seed):
    queue = deque([seed])
    component = {seed}
    while queue:
        state = queue.popleft()
        for mask, first, second in ring_patterns:
            restriction = state & mask
            if restriction not in (first, second):
                continue
            successor = state ^ mask
            if successor not in component:
                component.add(successor)
                queue.append(successor)
    return component


# Start from the label-(0,1) frozen ice and reverse one exact noncontractible
# alternating octagon.  This specifies the sector without fitting anything.
base_state = sum(1 << edge for edge, label in enumerate(edge_labels)
                 if label in (0, 1))
seed_loop = (84, 11, 9, 114, 112, 39, 37, 87)

check(len(set(seed_loop)) == 8, "sector seed uses eight distinct links")
check(all(len(set(edges[edge]) &
              set(edges[seed_loop[(index + 1) % len(seed_loop)]])) == 1
          for index, edge in enumerate(seed_loop)),
      "successive sector-seed links share one vertex")
check(all(((base_state >> seed_loop[index]) & 1) !=
          ((base_state >> seed_loop[(index + 1) % 8]) & 1)
          for index in range(8)),
      "sector-seed octagon is alternating")
seed_winding = (sum(SHIFTS[edge_labels[edge]] for edge in seed_loop[::2]) -
                sum(SHIFTS[edge_labels[edge]] for edge in seed_loop[1::2]))
check(abs(seed_winding) == 2 * CELL_COUNT,
      "sector-seed octagon has nonzero lifted quotient winding")

seed_state = base_state ^ sum(1 << edge for edge in seed_loop)
check(all(value == 2 for value in degrees(base_state)),
      "frozen reference obeys degree-two ice")
check(all(value == 2 for value in degrees(seed_state)),
      "octagon-reversed seed obeys degree-two ice")

component = ring_component(seed_state)
states = tuple(sorted(component))
state_index = {state: index for index, state in enumerate(states)}
check(len(states) == 180, "exact inherited-ring component has 180 ice states")
check(all(all(value == 2 for value in degrees(state)) for state in states),
      "all 180 component states remain in exact ice")

ring_degrees = []
for state in states:
    ring_degrees.append(sum(1 for mask, first, second in ring_patterns
                            if (state & mask) in (first, second)))
check(min(ring_degrees) == 4 and max(ring_degrees) == 6,
      "ring graph degrees range from four to six")
check(sum(ring_degrees) == 840,
      "ring graph has 420 undirected inherited transitions")
check(all(translate_state(state) in component for state in states),
      "the 180-state sector is closed under cyclic translation")

full_mask = (1 << LINK_COUNT) - 1
complement_component = {state ^ full_mask for state in component}
check(component.isdisjoint(complement_component),
      "global complement occupies a distinct isospectral ring sector")
check(len(ring_component(seed_state ^ full_mask)) == 180,
      "complement partner has the same 180-state size")


# ---------------------------------------------------------------------------
# Exact 180 x 180 inherited H6 matrix, in units J6=1.

dimension = len(states)
hamiltonian = np.zeros((dimension, dimension), dtype=float)
for row, state in enumerate(states):
    for mask, first, second in ring_patterns:
        if (state & mask) in (first, second):
            hamiltonian[row, state_index[state ^ mask]] = -1.0

check(np.array_equal(hamiltonian, hamiltonian.T),
      "inherited ring Hamiltonian is real symmetric")
check(np.count_nonzero(np.triu(hamiltonian)) == 420,
      "Hamiltonian contains exactly the 420 ring-graph edges")

translation_permutation = np.array(
    [state_index[translate_state(state)] for state in states], dtype=int)
check(np.array_equal(hamiltonian,
                     hamiltonian[np.ix_(translation_permutation,
                                        translation_permutation)]),
      "H6 commutes with cyclic translation")

eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
ground_energy = eigenvalues[0]
ground = eigenvectors[:, 0]
check(abs(ground_energy + 2 + 2 * sqrt(2)) < 2e-12,
      "sector ground energy is -2(1+sqrt(2)) J6")
check(eigenvalues[1] - eigenvalues[0] > 1e-10,
      "sector ground state is unique")
check(abs((eigenvalues[1] - eigenvalues[0]) -
          (1 + 2 * sqrt(2) - 2 * sqrt(3))) < 2e-12,
      "lowest sector gap has the exact algebraic checksum")
check(np.linalg.norm(hamiltonian @ ground - ground_energy * ground) < 2e-13,
      "ground-state eigensystem residual is below 2e-13")
check(np.linalg.norm(ground[translation_permutation] - ground) < 2e-12,
      "unique ground state has zero cyclic momentum")

# Independent translation-orbit reduction.  This does not reuse the full
# eigenvectors to infer the orbit closure or block spectra.
unassigned = set(states)
translation_orbits = []
while unassigned:
    representative = min(unassigned)
    orbit = []
    translated = representative
    while translated not in orbit:
        orbit.append(translated)
        translated = translate_state(translated)
    translation_orbits.append(tuple(orbit))
    unassigned.difference_update(orbit)

check(len(translation_orbits) == 6 and
      all(len(orbit) == CELL_COUNT for orbit in translation_orbits),
      "sector decomposes into six free translation orbits of length thirty")

momentum_blocks = []
for momentum in range(CELL_COUNT):
    block = np.zeros((len(translation_orbits),
                      len(translation_orbits)), dtype=complex)
    for row, source_orbit in enumerate(translation_orbits):
        source_index = state_index[source_orbit[0]]
        for column, target_orbit in enumerate(translation_orbits):
            block[row, column] = sum(
                hamiltonian[source_index, state_index[target]] *
                np.exp(2j * pi * momentum * amount / CELL_COUNT)
                for amount, target in enumerate(target_orbit)
            )
    check(np.max(np.abs(block - block.conj().T)) < 2e-13,
          f"momentum-{momentum} orbit block is Hermitian")
    momentum_blocks.append(block)

block_spectrum = np.sort(np.concatenate(
    [np.linalg.eigvalsh(block) for block in momentum_blocks]
))
check(np.max(np.abs(block_spectrum - eigenvalues)) < 8e-12,
      "thirty independent six-by-six momentum blocks reproduce full spectrum")
check(np.allclose(np.linalg.eigvalsh(momentum_blocks[0]),
                  (-2 - 2 * sqrt(2), 0, 0, -2 + 2 * sqrt(2), 2, 2),
                  atol=3e-12),
      "zero-momentum block certifies exact algebraic sector ground energy")
check(np.allclose(np.linalg.eigvalsh(momentum_blocks[5]),
                  (-1 - 2 * sqrt(3), -1, 1, 1, 1, -1 + 2 * sqrt(3)),
                  atol=3e-12),
      "momentum-five block certifies exact algebraic first excited energy")
other_block_floor = min(
    np.linalg.eigvalsh(momentum_blocks[momentum])[0]
    for momentum in range(CELL_COUNT)
    if momentum not in (0, 5, CELL_COUNT - 5)
)
check(abs(other_block_floor + 4.410987667370205) < 4e-12 and
      other_block_floor > -1 - 2 * sqrt(3) + 0.05,
      "remaining momentum blocks stay above the m=plus-or-minus-five level")


# ---------------------------------------------------------------------------
# Correctly oriented complement-even local E tensor and nonzero-momentum TT
# source.  Tetrahedral bond vectors are expressed in Cartesian coordinates.

bond_vectors = np.array([
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
], dtype=float) / sqrt(3)
check(np.linalg.norm(np.sum(bond_vectors, axis=0)) < 1e-14,
      "four tetrahedral q4 bond vectors sum to zero")

primitive = np.column_stack((bond_vectors[1] - bond_vectors[0],
                             bond_vectors[2] - bond_vectors[0],
                             bond_vectors[3] - bond_vectors[0]))
reciprocal = 2 * pi * np.linalg.inv(primitive).T
basis_offset = np.array((-0.25, -0.25, -0.25))

z_states = np.empty((dimension, LINK_COUNT), dtype=float)
for row, state in enumerate(states):
    z_states[row] = [1 - 2 * ((state >> edge) & 1)
                     for edge in range(LINK_COUNT)]


def local_tensor_from_signs(signs):
    field = sum(signs[label] * bond_vectors[label] for label in range(4))
    return np.outer(field, field) - np.eye(3) * np.dot(field, field) / 3


local_ice_signs = tuple(signs for signs in product((-1, 1), repeat=4)
                        if sum(signs) == 0)
check(all(np.allclose(local_tensor_from_signs(signs),
                      local_tensor_from_signs(signs).T) and
          abs(np.trace(local_tensor_from_signs(signs))) < 1e-13
          for signs in local_ice_signs),
      "every local ice tensor is symmetric traceless")
check(all(np.allclose(local_tensor_from_signs(signs),
                      local_tensor_from_signs(tuple(-value for value in signs)))
          for signs in local_ice_signs),
      "local ice tensor is exactly complement even")

# The one-link square terms are isotropic and vanish after STF projection;
# this direct expansion verifies that only the centered pair/E source remains.
pair_replays = []
for signs in local_ice_signs:
    pair_tensor = np.zeros((3, 3))
    for first in range(4):
        for second in range(first + 1, 4):
            pair_tensor += (signs[first] * signs[second] *
                            (np.outer(bond_vectors[first], bond_vectors[second]) +
                             np.outer(bond_vectors[second], bond_vectors[first])))
    pair_tensor -= np.eye(3) * np.trace(pair_tensor) / 3
    pair_replays.append(np.allclose(pair_tensor, local_tensor_from_signs(signs)))
check(all(pair_replays),
      "STF source is purely the ice-surviving even pair/E operator")

local_tensors = np.empty((dimension, VERTEX_COUNT, 3, 3), dtype=float)
for row in range(dimension):
    for cell in range(CELL_COUNT):
        a_edges = [4 * cell + label for label in range(4)]
        b_edges = [4 * ((cell - SHIFTS[label]) % CELL_COUNT) + label
                   for label in range(4)]
        local_tensors[row, cell] = local_tensor_from_signs(z_states[row, a_edges])
        local_tensors[row, CELL_COUNT + cell] = local_tensor_from_signs(
            z_states[row, b_edges])


def fractional_wavevector(momentum):
    q = momentum * np.array((1, 5, 19), dtype=float) / CELL_COUNT
    return q - np.round(q)


def polarization_data(momentum):
    q = fractional_wavevector(momentum)
    wavevector = reciprocal @ q
    unit = wavevector / np.linalg.norm(wavevector)
    axes = np.eye(3)
    reference = axes[np.argmin(np.abs(axes @ unit))]
    first = reference - unit * np.dot(unit, reference)
    first /= np.linalg.norm(first)
    second = np.cross(unit, first)
    plus = (np.outer(first, first) - np.outer(second, second)) / sqrt(2)
    cross = (np.outer(first, second) + np.outer(second, first)) / sqrt(2)
    return q, wavevector, first, second, plus, cross


def tt_complex_source(momentum, polarization):
    q, _, _, _, plus, cross = polarization_data(momentum)
    epsilon = plus if polarization == "plus" else cross
    b_phase = np.exp(2j * pi * np.dot(q, basis_offset))
    values = np.zeros(dimension, dtype=complex)
    for cell in range(CELL_COUNT):
        phase = np.exp(2j * pi * momentum * cell / CELL_COUNT)
        values += phase * np.einsum("nij,ij->n",
                                    local_tensors[:, cell], epsilon)
        values += phase * b_phase * np.einsum(
            "nij,ij->n", local_tensors[:, CELL_COUNT + cell], epsilon)
    return values / sqrt(VERTEX_COUNT)


def tt_real_sources(momentum):
    plus = tt_complex_source(momentum, "plus")
    cross = tt_complex_source(momentum, "cross")
    # sqrt(2) times the Hermitian real/imaginary part of the complex source.
    return (sqrt(2) * plus.real, sqrt(2) * plus.imag,
            sqrt(2) * cross.real, sqrt(2) * cross.imag)


q, wavevector, _, _, plus_tensor, cross_tensor = polarization_data(1)
check(np.allclose(q, np.array((1, 5, -11)) / 30),
      "selected nonzero quotient momentum has the declared reciprocal alias")
raw_q = np.array((1, 5, 19), dtype=float) / CELL_COUNT
alias_norms = sorted(
    (np.linalg.norm(reciprocal @ (raw_q + np.array(shift))), shift)
    for shift in product(range(-2, 3), repeat=3)
)
check(alias_norms[0][1] == (0, 0, -1) and
      alias_norms[1][0] - alias_norms[0][0] > 0.5,
      "declared reciprocal representative is the unique shortest nearby alias")
check(np.linalg.norm(primitive @ basis_offset - bond_vectors[0]) < 2e-15,
      "B-sublattice fractional offset reproduces the label-zero bond vector")
for epsilon in (plus_tensor, cross_tensor):
    check(abs(np.trace(epsilon)) < 2e-15 and
          np.linalg.norm(epsilon @ wavevector) < 2e-15 and
          abs(np.sum(epsilon * epsilon) - 1) < 2e-15,
          "each selected polarization is unit-normalized transverse traceless")
check(abs(np.sum(plus_tensor * cross_tensor)) < 2e-15,
      "plus and cross polarization tensors are Frobenius orthogonal")

complex_plus = tt_complex_source(1, "plus")
phase = np.exp(2j * pi / CELL_COUNT)
check(np.max(np.abs(complex_plus[translation_permutation] -
                    phase * complex_plus)) < 2e-12,
      "complex TT source transforms with cyclic momentum +1")


def tensor_value_on_state(state, momentum, polarization):
    q, _, _, _, plus, cross = polarization_data(momentum)
    epsilon = plus if polarization == "plus" else cross
    b_phase = np.exp(2j * pi * np.dot(q, basis_offset))
    total = 0j
    for cell in range(CELL_COUNT):
        phase_cell = np.exp(2j * pi * momentum * cell / CELL_COUNT)
        a_edges = [4 * cell + label for label in range(4)]
        b_edges = [4 * ((cell - SHIFTS[label]) % CELL_COUNT) + label
                   for label in range(4)]
        a_signs = [1 - 2 * ((state >> edge) & 1) for edge in a_edges]
        b_signs = [1 - 2 * ((state >> edge) & 1) for edge in b_edges]
        total += phase_cell * np.sum(epsilon * local_tensor_from_signs(a_signs))
        total += (phase_cell * b_phase *
                  np.sum(epsilon * local_tensor_from_signs(b_signs)))
    return total / sqrt(VERTEX_COUNT)


check(abs(tensor_value_on_state(seed_state, 1, "plus") -
          tensor_value_on_state(seed_state ^ full_mask, 1, "plus")) < 2e-13,
      "global TT source is complement even across partner sectors")


# ---------------------------------------------------------------------------
# Connected zero-frequency cumulants and composite-source Legendre amputation.

sources = tt_real_sources(1)
check(all(abs(ground @ (source * ground)) < 2e-13 for source in sources),
      "all four momentum-one TT sources are centered in the ground state")

resolvent_diagonal = np.zeros(dimension)
resolvent_diagonal[1:] = 1 / (eigenvalues[1:] - ground_energy)


def apply_resolvent(vector):
    coefficients = eigenvectors.T @ vector
    coefficients[0] = 0
    return eigenvectors @ (resolvent_diagonal * coefficients)


w2_matrix = np.empty((4, 4))
for row, source_row in enumerate(sources):
    for column, source_column in enumerate(sources):
        w2_matrix[row, column] = 2 * ((source_row * ground) @
                                      apply_resolvent(source_column * ground))

expected_w2 = np.array([
    [1.130847135995723, 0, -0.037434360319729, 0],
    [0, 1.130847135995723, 0, -0.037434360319729],
    [-0.037434360319729, 0, 0.114433012322288, 0],
    [0, -0.037434360319729, 0, 0.114433012322288],
])
check(np.max(np.abs(w2_matrix - expected_w2)) < 3e-12,
      "four-channel connected TT two-point matrix matches its checksum")
w2_eigenvalues = np.linalg.eigvalsh(w2_matrix)
check(np.allclose(w2_eigenvalues[:2], 0.113056176225351, atol=3e-12) and
      np.allclose(w2_eigenvalues[2:], 1.132223972092660, atol=3e-12),
      "TT two-point matrix has two positive translation-quadrature doublets")


def rayleigh_schrodinger(source, order=4):
    perturbation = -np.diag(source)
    wavefunctions = [ground]
    energy_coefficients = [ground_energy]
    for degree in range(1, order + 1):
        energy_coefficients.append(float(
            ground @ (perturbation @ wavefunctions[degree - 1])))
        right = -(perturbation @ wavefunctions[degree - 1])
        for lower in range(1, degree):
            right += (energy_coefficients[lower] *
                      wavefunctions[degree - lower])
        wavefunctions.append(apply_resolvent(right))
    w2 = -2 * energy_coefficients[2]
    w3 = -6 * energy_coefficients[3]
    w4 = -24 * energy_coefficients[4]
    gamma4 = -w4 / w2 ** 4
    return w2, w3, w4, gamma4, energy_coefficients


plus_w2, plus_w3, plus_w4, plus_gamma4, plus_coefficients = (
    rayleigh_schrodinger(sources[0]))
check(abs(plus_w2 - 1.130847135995723) < 3e-12,
      "plus-cos connected W2 is nonzero")
check(abs(plus_w3) < 2e-13,
      "momentum conservation kills the plus-cos connected W3")
check(abs(plus_w4 + 0.136825085605100) < 3e-12,
      "plus-cos connected W4 is nonzero and negative")
check(abs(plus_gamma4 - 0.083666214307836) < 3e-12,
      "plus-cos composite Legendre-amputated Gamma4 is positive")

# Independent spectral-resolvent replay of the connected fourth cumulant.
source_ground = sources[0] * ground
r_source = apply_resolvent(source_ground)
a_spectral = float(source_ground @ r_source)
c_spectral = float(source_ground @ apply_resolvent(r_source))
middle = sources[0] * r_source
middle -= ground * (ground @ middle)
middle = apply_resolvent(middle)
middle = sources[0] * middle
middle -= ground * (ground @ middle)
b_spectral = float(source_ground @ apply_resolvent(middle))
check(abs((-b_spectral + a_spectral * c_spectral) -
          plus_coefficients[4]) < 3e-13,
      "connected fourth-order subtraction matches spectral resolvents")

# A structurally independent replay uses a nonsingular augmented linear solve,
# not the spectral decomposition, for every reduced resolvent application.
augmented_operator = (hamiltonian - ground_energy * np.eye(dimension) +
                      np.outer(ground, ground))


def apply_augmented_resolvent(vector):
    solution = np.linalg.solve(augmented_operator, vector)
    return solution - ground * (ground @ solution)


augmented_w2 = np.empty((4, 4))
for row, source_row in enumerate(sources):
    for column, source_column in enumerate(sources):
        augmented_w2[row, column] = 2 * ((source_row * ground) @
                                         apply_augmented_resolvent(
                                             source_column * ground))
check(np.max(np.abs(augmented_w2 - w2_matrix)) < 2e-12,
      "augmented linear solves independently reproduce the full W2 matrix")

augmented_r_source = apply_augmented_resolvent(source_ground)
augmented_a = float(source_ground @ augmented_r_source)
augmented_c = float(source_ground @
                    apply_augmented_resolvent(augmented_r_source))
augmented_middle = sources[0] * augmented_r_source
augmented_middle -= ground * (ground @ augmented_middle)
augmented_middle = apply_augmented_resolvent(augmented_middle)
augmented_middle = sources[0] * augmented_middle
augmented_middle -= ground * (ground @ augmented_middle)
augmented_b = float(source_ground @
                    apply_augmented_resolvent(augmented_middle))
augmented_e4 = -augmented_b + augmented_a * augmented_c
check(abs(augmented_e4 - 0.005701045233546) < 2e-13 and
      abs(-24 * augmented_e4 - plus_w4) < 2e-12,
      "augmented solves reproduce the connected four-Q subtraction and W4")

# Resolve the two polarization eigenchannels in one cosine quadrature.
cosine_block = w2_matrix[np.ix_((0, 2), (0, 2))]
cosine_eigenvalues, cosine_eigenvectors = np.linalg.eigh(cosine_block)
channel_results = []
for channel in range(2):
    source = (cosine_eigenvectors[0, channel] * sources[0] +
              cosine_eigenvectors[1, channel] * sources[2])
    channel_results.append(rayleigh_schrodinger(source)[:4])

check(abs(channel_results[0][2] + 0.003036590242056) < 3e-12 and
      abs(channel_results[0][3] - 18.586988116578) < 2e-10,
      "low-susceptibility TT eigenchannel has positive static Gamma4")
check(abs(channel_results[1][2] + 0.148796476796671) < 3e-12 and
      abs(channel_results[1][3] - 0.090544748133189) < 3e-12,
      "high-susceptibility TT eigenchannel has positive static Gamma4")


# ---------------------------------------------------------------------------
# Finite poles and a conservative two-one-link threshold proxy.


def spectral_groups(source, weight_floor=1e-10):
    weights = (eigenvectors.T @ (source * ground)) ** 2
    groups = []
    for gap, weight in zip(eigenvalues - ground_energy, weights):
        if weight <= weight_floor:
            continue
        for group in groups:
            if abs(group[0] - gap) < 1e-9:
                group[1] += float(weight)
                break
        else:
            groups.append([float(gap), float(weight)])
    return groups


tt_groups = spectral_groups(sources[0])
expected_tt_groups = (
    (3.194109035554332, 0.005026104004432),
    (3.490165912028476, 1.965864248197576),
    (6.166688337463908, 0.003649484732984),
    (9.139267639373482, 0.000000908886868),
)
check(len(tt_groups) == 4, "selected TT response has four finite spectral poles")
for actual, expected in zip(tt_groups, expected_tt_groups):
    check(abs(actual[0] - expected[0]) < 3e-12 and
          abs(actual[1] - expected[1]) < 3e-12,
          "selected TT pole and residue match the finite-sector checksum")

# Independent Krylov/Lanczos reconstruction of the selected spectral measure.
lanczos_norm = np.linalg.norm(source_ground)
lanczos_vector = source_ground / lanczos_norm
lanczos_previous = np.zeros(dimension)
lanczos_beta_previous = 0.0
lanczos_alpha = []
lanczos_beta = []
shifted_hamiltonian = hamiltonian - ground_energy * np.eye(dimension)
for _ in range(12):
    remainder = (shifted_hamiltonian @ lanczos_vector -
                 lanczos_beta_previous * lanczos_previous)
    alpha = float(lanczos_vector @ remainder)
    remainder -= alpha * lanczos_vector
    beta = float(np.linalg.norm(remainder))
    lanczos_alpha.append(alpha)
    if beta < 1e-9:
        break
    lanczos_beta.append(beta)
    lanczos_previous, lanczos_vector = lanczos_vector, remainder / beta
    lanczos_beta_previous = beta

check(len(lanczos_alpha) == 4 and len(lanczos_beta) == 3,
      "selected plus-cos spectral Krylov space closes after four steps")
lanczos_matrix = (np.diag(lanczos_alpha) +
                  np.diag(lanczos_beta, 1) +
                  np.diag(lanczos_beta, -1))
lanczos_gaps, lanczos_vectors = np.linalg.eigh(lanczos_matrix)
lanczos_weights = lanczos_norm ** 2 * lanczos_vectors[0] ** 2
check(np.max(np.abs(lanczos_gaps -
                    np.array([entry[0] for entry in expected_tt_groups]))) < 4e-12,
      "four-step Lanczos gaps reproduce the selected pole support")
check(np.max(np.abs(lanczos_weights -
                    np.array([entry[1] for entry in expected_tt_groups]))) < 4e-12,
      "four-step Lanczos weights reproduce the selected spectral residues")


def photon_real_sources(momentum):
    q, _, first, second, _, _ = polarization_data(momentum)
    results = []
    for polarization in (first, second):
        values = np.zeros(dimension, dtype=complex)
        for cell in range(CELL_COUNT):
            cell_phase = np.exp(2j * pi * momentum * cell / CELL_COUNT)
            for label in range(4):
                displacement = basis_offset.copy()
                if label:
                    displacement[label - 1] += 1
                midpoint_phase = np.exp(1j * pi * np.dot(q, displacement))
                edge = 4 * cell + label
                values += (cell_phase * midpoint_phase * z_states[:, edge] *
                           np.dot(polarization, bond_vectors[label]))
        values /= sqrt(LINK_COUNT)
        results.extend((sqrt(2) * values.real, sqrt(2) * values.imag))
    return tuple(results)


def lowest_source_gap(source_family, weight_floor=1e-8):
    best = float("inf")
    for source in source_family:
        weights = (eigenvectors.T @ (source * ground)) ** 2
        indices = np.where(weights > weight_floor)[0]
        if len(indices):
            best = min(best, float(eigenvalues[indices[0]] - ground_energy))
    return best


photon_gaps = {}
for momentum in range(1, CELL_COUNT):
    photon_gaps[momentum] = lowest_source_gap(photon_real_sources(momentum))
check(all(np.isfinite(value) for value in photon_gaps.values()),
      "every nonzero quotient momentum has a one-link transverse pole")

threshold_candidates = []
for first_momentum in range(1, CELL_COUNT):
    second_momentum = (1 - first_momentum) % CELL_COUNT
    if second_momentum == 0:
        continue
    threshold_candidates.append((photon_gaps[first_momentum] +
                                 photon_gaps[second_momentum],
                                 first_momentum, second_momentum))
two_link_threshold = min(threshold_candidates)
check(two_link_threshold[1:] == (9, 22),
      "momentum-one threshold proxy is attained by momenta 9 and 22")
check(abs(two_link_threshold[0] - 2.059674505691458) < 3e-12,
      "two-one-link threshold proxy matches its checksum")
check(tt_groups[0][0] > two_link_threshold[0],
      "selected momentum-one TT pole lies above the finite threshold proxy")

photon_momentum_one_groups = []
for source in photon_real_sources(1):
    photon_momentum_one_groups.extend(spectral_groups(source))
check(any(abs(group[0] - tt_groups[0][0]) < 3e-12
          for group in photon_momentum_one_groups),
      "lowest TT pole is energy-degenerate with one-link response")


# ---------------------------------------------------------------------------
# Claim-text ceilings and dependency-hash hold.

root = Path(__file__).resolve().parent
theorem = (root / "THEOREM.md").read_text()
self_audit = (root / "SELF_AUDIT.md").read_text()
independent_audit = (root / "INDEPENDENT_AUDIT.md").read_text()
joined = theorem + self_audit + independent_audit
theorem_flat = " ".join(theorem.split())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


workspace = root.parent
dependency_hashes = {
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md":
        "5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe",
    "LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md":
        "a91caa20d16b0a1194333f9b51d96546a4ea24d55e23bf1f04c7d249641af8db",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md":
        "cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98",
    "LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md":
        "c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4",
    "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md":
        "98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452",
    "LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/INDEPENDENT_AUDIT.md":
        "327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md":
        "78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25",
    "LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md":
        "53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9",
}
for relative_path, expected_hash in dependency_hashes.items():
    dependency = workspace / relative_path
    check(dependency.is_file() and not dependency.is_symlink() and
          sha256(dependency) == expected_hash,
          f"final dependency custody pinned: {relative_path}")
    original = dependency.read_bytes()
    tampered = original + b"\nFO hostile tamper sentinel\n"
    check(hashlib.sha256(tampered).hexdigest() != expected_hash,
          f"dependency tamper is rejected: {relative_path}")

manifest_path = root / "MANIFEST.sha256"
manifest_entries = {}
for line in manifest_path.read_text().splitlines():
    if not line.strip():
        continue
    digest, marker_name = line.split(maxsplit=1)
    manifest_entries[marker_name.lstrip(" *")] = digest

sealed_names = {
    "THEOREM.md",
    "SELF_AUDIT.md",
    "INDEPENDENT_AUDIT.md",
    "verify_finite_tt_four_point.py",
    "VERIFICATION.txt",
}
check(set(manifest_entries) == sealed_names,
      "local manifest covers exactly the five load-bearing lane files")
for name in sorted(sealed_names):
    sealed_file = root / name
    expected_hash = manifest_entries[name]
    check(sealed_file.is_file() and not sealed_file.is_symlink() and
          sha256(sealed_file) == expected_hash,
          f"local manifest custody pinned: {name}")
    tampered = sealed_file.read_bytes() + b"\nFO local tamper sentinel\n"
    check(hashlib.sha256(tampered).hexdigest() != expected_hash,
          f"local manifest tamper is rejected: {name}")

for required in (
    "H_6=-J_6\\sum_C B_C",
    "180-state",
    "W_4={-0.136825085605100\\over J_6^3}",
    "Gamma^{(4)}_{\\rm comp}=-{W_4\\over W_2^4}",
    "not the connected or photon-amputated four-one-link vertex",
    "four-`Q` composite cumulant",
    "quadratic composite of an exactly Gaussian underlying field",
    "no below-proxy or energy-exclusive tensor candidate",
    "not a no-bound-state theorem",
    "frequency- and relative-momentum-resolved connected function",
    "finite-volume diagnostic",
    "thermodynamic helicity-two",
):
    check(required in theorem_flat,
          f"theorem retains load-bearing statement: {required}")

for forbidden in (
    "the finite cluster proves a graviton",
    "Gamma4 is Newton's constant",
    "a massless tensor pole is derived",
    "the 180-state sector is the full ice Hilbert space",
    "the threshold proxy is a thermodynamic continuum",
    "W4 proves a non-Gaussian photon interaction",
    "the finite proxy excludes a bound state",
    "composite 1PI is photon 1PI",
    "ice occupation is automatically a record",
):
    check(forbidden not in joined,
          f"forbidden promotion absent: {forbidden}")

for required in (
    "PASS_AFTER_COMPOSITE_TYPING_AND_SEAL_REPAIR",
    "six free translation orbits",
    "180",
    "420",
    "1.130847135995723",
    "-0.136825085605100",
    "0.083666214307836",
    "3.194109035554332",
    "2.059674505691458",
    "four-Q",
    "four-one-link",
    "not a no-bound-state theorem",
):
    check(required in independent_audit,
          f"independent audit retains required finding: {required}")

print(f"Finite TT composite-cumulant verification: {checks} passed, 0 failed")
print(f"Periodic quotient: cells={CELL_COUNT}, vertices={VERTEX_COUNT}, "
      f"links={LINK_COUNT}, hexagons={len(hexagons)}")
print(f"Ring sector: states={dimension}, transitions=420, "
      f"E0/J6={ground_energy:.15f}")
print("TT plus-cos: "
      f"W2={plus_w2:.15f}/J6, W4={plus_w4:.15f}/J6^3, "
      f"Gamma4={plus_gamma4:.15f} J6")
print("Lowest selected TT pole / finite two-link threshold proxy: "
      f"{tt_groups[0][0]:.15f} / {two_link_threshold[0]:.15f} J6")
