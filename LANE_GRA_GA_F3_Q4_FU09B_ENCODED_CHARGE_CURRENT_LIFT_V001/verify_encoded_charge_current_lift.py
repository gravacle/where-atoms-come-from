#!/usr/bin/env python3
"""Exact finite verifier for the minimal FU09b encoded-current lift.

This is a U(1) charge/current ancestry calculation.  It does not construct
T^{0j}, a metric Ward identity, continuum locality, or gravity.
"""

import json
from hashlib import sha256
from pathlib import Path

import numpy as np


LANE = Path(__file__).resolve().parent
checks = 0


def check(condition, label):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS {label}")


def commutator(first, second):
    return first @ second - second @ first


def kron_all(factors):
    result = np.array(((1,),), dtype=complex)
    for factor in factors:
        result = np.kron(result, factor)
    return result


def embed_local(operator, site, count, local_identity):
    return kron_all(tuple(operator if index == site else local_identity
                          for index in range(count)))


for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = LANE.parent / relative
    check(path.is_file() and not path.is_symlink() and
          sha256(path.read_bytes()).hexdigest() == expected,
          f"dependency custody: {relative}")


# -------------------------------------------------------------------------
# One link plus one paired reservoir.

I2 = np.eye(2, dtype=complex)
Z = np.diag((-1, 1)).astype(complex)
X = np.array(((0, 1), (1, 0)), dtype=complex)
sigma_plus = np.array(((0, 0), (1, 0)), dtype=complex)
sigma_minus = sigma_plus.T

# Reservoir basis is |+q_R>,|-q_R>, so T_- maps +q_R to -q_R.
QR = np.diag((1, -1)).astype(complex)
T_MINUS = np.array(((0, 0), (1, 0)), dtype=complex)
T_PLUS = T_MINUS.T

Z_LIFT = np.kron(Z, I2)
QR_LIFT = np.kron(I2, QR)
A_TRANSFER = np.kron(sigma_plus, T_MINUS)
B_TRANSFER = np.kron(sigma_minus, T_PLUS)
X_DRESSED = A_TRANSFER + B_TRANSFER
Q_TOTAL = Z_LIFT + QR_LIFT

# V|-q> = |-q,+q_R>, V|+q> = |+q,-q_R>.
V = np.zeros((4, 2), dtype=complex)
V[0, 0] = 1
V[3, 1] = 1
P_ENC = V @ V.T.conj()

check(np.array_equal(V.T.conj() @ V, I2) and
      np.array_equal(P_ENC @ P_ENC, P_ENC),
      "the fixed-total-charge encoder is an exact two-dimensional isometry")
check(np.array_equal(Q_TOTAL @ V, np.zeros_like(V)),
      "the encoded link-reservoir pair has exact total charge zero")
check(np.array_equal(Z_LIFT @ V, V @ Z) and
      np.array_equal(X_DRESSED @ V, V @ X),
      "dressed Z and X exactly intertwine with the inherited link operators")
check(np.array_equal(commutator(Q_TOTAL, X_DRESSED),
                     np.zeros((4, 4), dtype=complex)),
      "the FU09b dressed flip conserves total charge on the full pair space")
check(np.array_equal(commutator(P_ENC, X_DRESSED),
                     np.zeros((4, 4), dtype=complex)),
      "the fixed-total encoded subspace is invariant under the dressed flip")

# With q_*=h=hbar=1 and H_flip=-X_dressed.
H_FLIP = -X_DRESSED
Q_LINK_DOT = 1j * commutator(H_FLIP, Z_LIFT)
Q_RES_DOT = 1j * commutator(H_FLIP, QR_LIFT)
I_LINK_TO_RES = -Q_LINK_DOT
check(np.array_equal(Q_LINK_DOT, 2j*(A_TRANSFER-B_TRANSFER)) and
      np.array_equal(Q_RES_DOT, -Q_LINK_DOT) and
      np.any(I_LINK_TO_RES),
      "the dressed flip has an exact nonzero compensating U(1) current")
check(np.array_equal(Q_LINK_DOT + I_LINK_TO_RES,
                     np.zeros((4, 4), dtype=complex)) and
      np.array_equal(Q_RES_DOT - I_LINK_TO_RES,
                     np.zeros((4, 4), dtype=complex)),
      "link and paired-reservoir charge continuity closes exactly")


# -------------------------------------------------------------------------
# Four-link q4 node: exact tensor-product equivalence and spatial source.

LINKS = 4
I4 = np.eye(4, dtype=complex)
V4 = kron_all((V,)*LINKS)
check(V4.shape == (256, 16) and
      np.array_equal(V4.T.conj() @ V4, np.eye(16)),
      "four paired reservoirs encode the complete sixteen-state q4 node")

original_z = tuple(embed_local(Z, site, LINKS, I2) for site in range(LINKS))
original_x = tuple(embed_local(X, site, LINKS, I2) for site in range(LINKS))
lifted_z = tuple(embed_local(Z_LIFT, site, LINKS, I4)
                 for site in range(LINKS))
lifted_x = tuple(embed_local(X_DRESSED, site, LINKS, I4)
                 for site in range(LINKS))
lifted_qr = tuple(embed_local(QR_LIFT, site, LINKS, I4)
                  for site in range(LINKS))

check(all(np.array_equal(lifted_z[site] @ V4,
                         V4 @ original_z[site]) and
          np.array_equal(lifted_x[site] @ V4,
                         V4 @ original_x[site])
          for site in range(LINKS)),
      "all four dressed link algebras intertwine simultaneously")

sum_z = sum(original_z, np.zeros((16, 16), dtype=complex))
sum_z_lift = sum(lifted_z, np.zeros((256, 256), dtype=complex))
H_DEGREE4 = sum_z @ sum_z                    # 4(d-2)^2
H_DEGREE4_LIFT = sum_z_lift @ sum_z_lift
H_NODE = H_DEGREE4 - sum(original_x)
H_NODE_LIFT = H_DEGREE4_LIFT - sum(lifted_x)
check(np.array_equal(H_NODE_LIFT @ V4, V4 @ H_NODE),
      "the dressed source-off degree-plus-flip node is exactly inherited")

# A permitted reference shift must be one common scalar on the complete
# encoded Hilbert space, including both the ice P sector and virtual/off-ice
# Q sector.  It need not be an identity on the unused dark complement of the
# physical 256-state lift.  The exact encoded statement is therefore modulo
# P_CODE = V4 V4^dagger, not literal equality of the Hamiltonians.
P_CODE4 = V4 @ V4.T.conj()
COMMON_SCALAR = 7
H_NODE_REFERENCE_LIFT = H_NODE_LIFT + COMMON_SCALAR*P_CODE4
check(np.array_equal(H_NODE_REFERENCE_LIFT @ V4,
                     V4 @ (H_NODE+COMMON_SCALAR*np.eye(16))),
      "one full-code scalar gives the exact modulo-reference intertwiner")

SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))
PAIRS = tuple((first, second) for first in range(4)
              for second in range(first+1, 4))


def dyad_coordinates(vector):
    x, y, z = vector
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


# 16 Q_pair/U_d at lambda=-1/2 and 3 Q_flip/h are integer matrices.
q_pair = [np.zeros((16, 16), dtype=complex) for _ in range(6)]
q_pair_lift = [np.zeros((256, 256), dtype=complex) for _ in range(6)]
for first, second in PAIRS:
    root = tuple(SIGNS[second][axis]-SIGNS[first][axis]
                 for axis in range(3))
    coordinates = dyad_coordinates(root)
    pair = original_z[first] @ original_z[second]
    pair_lift = lifted_z[first] @ lifted_z[second]
    for component, coefficient in enumerate(coordinates):
        q_pair[component] -= coefficient * pair
        q_pair_lift[component] -= coefficient * pair_lift

q_flip = [np.zeros((16, 16), dtype=complex) for _ in range(6)]
q_flip_lift = [np.zeros((256, 256), dtype=complex) for _ in range(6)]
for link, signs in enumerate(SIGNS):
    coordinates = dyad_coordinates(signs)
    for component, coefficient in enumerate(coordinates):
        q_flip[component] -= coefficient * original_x[link]
        q_flip_lift[component] -= coefficient * lifted_x[link]

check(all(np.array_equal(q_pair_lift[component] @ V4,
                         V4 @ q_pair[component])
          for component in range(6)),
      "the complete six-coordinate Coulomb pair source intertwines exactly")
check(all(np.array_equal(q_flip_lift[component] @ V4,
                         V4 @ q_flip[component])
          for component in range(6)),
      "the complete six-coordinate hopping-numerator source intertwines exactly")

source_inventory = tuple(q_pair) + tuple(q_flip)
source_inventory_lift = tuple(q_pair_lift) + tuple(q_flip_lift)
source_gram = np.array([[np.trace(first.T.conj() @ second)
                         for second in source_inventory]
                        for first in source_inventory])
source_gram_lift = np.array([[np.trace(
    (V4.T.conj() @ first @ V4).T.conj() @
    (V4.T.conj() @ second @ V4))
    for second in source_inventory_lift] for first in source_inventory_lift])
check(np.array_equal(source_gram_lift, source_gram),
      "encoded restriction preserves the full pair/flip source Gram matrix")

# A source derivative of the common scalar is a declared reference/identity
# source.  It may change the uncentered source, but it cancels exactly from
# the centered nonidentity operator used for ranks, commutators, and connected
# response.
IDENTITY_SOURCE = -5
q0_reference_lift = q_pair_lift[0] + IDENTITY_SOURCE*P_CODE4
q0_reference = q_pair[0] + IDENTITY_SOURCE*np.eye(16)
q0_restricted = V4.T.conj() @ q0_reference_lift @ V4
q0_centered = q0_reference - np.trace(q0_reference)/16*np.eye(16)
q0_original_centered = q_pair[0] - np.trace(q_pair[0])/16*np.eye(16)
check(np.array_equal(q0_reference_lift @ V4, V4 @ q0_reference) and
      np.array_equal(q0_restricted, q0_reference) and
      np.array_equal(q0_centered, q0_original_centered),
      "a full-code identity source changes only the declared reference part")

# Feshbach check on the inherited sixteen-state q4 node.  P is the six-state
# sum-Z=0 sector and Q is its ten-state virtual/off-ice complement.  Under a
# common scalar c on P+Q and the consistent reference shift z -> z+c, the
# Q-resolvent denominator is literally unchanged and the Schur/Feshbach
# operator changes only by c I_P.
p_indices = tuple(index for index in range(16)
                  if sum_z[index, index] == 0)
q_indices = tuple(index for index in range(16)
                  if index not in p_indices)
h_pp = H_NODE[np.ix_(p_indices, p_indices)]
h_pq = H_NODE[np.ix_(p_indices, q_indices)]
h_qp = H_NODE[np.ix_(q_indices, p_indices)]
h_qq = H_NODE[np.ix_(q_indices, q_indices)]
z_reference = 100
denominator = z_reference*np.eye(len(q_indices))-h_qq
h_shifted = H_NODE+COMMON_SCALAR*np.eye(16)
h_shifted_pp = h_shifted[np.ix_(p_indices, p_indices)]
h_shifted_pq = h_shifted[np.ix_(p_indices, q_indices)]
h_shifted_qp = h_shifted[np.ix_(q_indices, p_indices)]
h_shifted_qq = h_shifted[np.ix_(q_indices, q_indices)]
shifted_denominator = ((z_reference+COMMON_SCALAR)*np.eye(len(q_indices))-
                       h_shifted_qq)
h_eff = h_pp+h_pq @ np.linalg.inv(denominator) @ h_qp
h_shifted_eff = (h_shifted_pp+h_shifted_pq @
                 np.linalg.inv(shifted_denominator) @ h_shifted_qp)
check(len(p_indices) == 6 and len(q_indices) == 10 and
      np.array_equal(shifted_denominator, denominator) and
      np.max(np.abs(h_shifted_eff-
                    (h_eff+COMMON_SCALAR*np.eye(len(p_indices))))) < 1e-13,
      "a common P+Q scalar cancels from virtual gaps and shifts Feshbach by identity")


# -------------------------------------------------------------------------
# Native m=1 exchange continuity on the four link-midpoint supports.

SUPPORT_EXPONENTS_240 = (5, 9, 25, 201)
phases = tuple(np.exp(2j*np.pi*exponent/240)
               for exponent in SUPPORT_EXPONENTS_240)
H_FLIPS4 = -sum(lifted_x)
rho_m = sum((phase * lifted_z[link]
             for link, phase in enumerate(phases)),
            np.zeros((256, 256), dtype=complex))
rho_dot_m = 1j*commutator(H_FLIPS4, rho_m)
currents4 = tuple(-1j*commutator(H_FLIPS4, lifted_z[link])
                  for link in range(LINKS))
i_port_m = sum((phase * currents4[link]
                for link, phase in enumerate(phases)),
               np.zeros((256, 256), dtype=complex))
check(np.max(np.abs(rho_dot_m+i_port_m)) < 2e-13,
      "co-located m=1 link-density plus reservoir-port current closes")
check(np.linalg.norm(i_port_m) > 0,
      "the co-located nonzero-momentum port-current operator is nontrivial")

# Placement is not fixed by FU09b.  Co-locating each paired reservoir with its
# link makes the full Fourier density conserved.  Moving all reservoirs to one
# phase produces a nonzero finite-m exchange which needs a connector/boundary
# current; global m=0 charge remains conserved either way.
rho_full_colocated = sum((phase*(lifted_z[link]+lifted_qr[link])
                          for link, phase in enumerate(phases)),
                         np.zeros((256, 256), dtype=complex))
rho_full_displaced = (
    sum((phase*lifted_z[link] for link, phase in enumerate(phases)),
        np.zeros((256, 256), dtype=complex)) +
    sum(lifted_qr, np.zeros((256, 256), dtype=complex)))
rho_full_m0 = sum((lifted_z[link]+lifted_qr[link]
                   for link in range(LINKS)),
                  np.zeros((256, 256), dtype=complex))
check(np.max(np.abs(commutator(H_FLIPS4, rho_full_colocated))) < 2e-13,
      "co-located link/reservoir allocation conserves full m=1 charge")
check(np.linalg.norm(commutator(H_FLIPS4, rho_full_displaced)) > 1e-8,
      "a displaced reservoir allocation requires a nonzero connector current")
check(np.array_equal(commutator(H_FLIPS4, rho_full_m0),
                     np.zeros((256, 256), dtype=complex)),
      "global m=0 total charge is allocation independent")


# -------------------------------------------------------------------------
# Explicit outer-port boundary term and its encoded-hold obstruction.

# External port uses the same +/- charge basis.  The reservoir/external
# exchange conserves QR+Qext but changes QR without changing the link.
Q_EXT = QR.copy()
R_PLUS = T_PLUS
R_MINUS = T_MINUS
E_PLUS = T_PLUS
E_MINUS = T_MINUS
H_BOUNDARY_RE = -(np.kron(R_MINUS, E_PLUS) +
                  np.kron(R_PLUS, E_MINUS))
H_FLIP_LRE = np.kron(H_FLIP, I2)
H_BOUNDARY_LRE = np.kron(I2, H_BOUNDARY_RE)
Q_LINK_LRE = np.kron(Z_LIFT, I2)
Q_RES_LRE = np.kron(QR_LIFT, I2)
Q_EXT_LRE = np.kron(np.eye(4), Q_EXT)
H_LRE = H_FLIP_LRE + H_BOUNDARY_LRE

link_dot = 1j*commutator(H_LRE, Q_LINK_LRE)
res_dot = 1j*commutator(H_LRE, Q_RES_LRE)
ext_dot = 1j*commutator(H_LRE, Q_EXT_LRE)
i_link_res = -link_dot
i_res_ext = ext_dot
check(np.array_equal(link_dot+i_link_res, np.zeros((8, 8))) and
      np.array_equal(res_dot-i_link_res+i_res_ext, np.zeros((8, 8))) and
      np.array_equal(ext_dot-i_res_ext, np.zeros((8, 8))),
      "link/reservoir/exterior continuity retains both port boundary terms")
check(np.array_equal(link_dot+res_dot+ext_dot, np.zeros((8, 8))),
      "the complete three-factor U(1) charge is exactly conserved")
P_ENC_LRE = np.kron(P_ENC, I2)
check(np.any(commutator(P_ENC_LRE, H_BOUNDARY_LRE)),
      "an active reservoir-only outer port leaks from the encoded hold subspace")


# -------------------------------------------------------------------------
# Exact source-preservation and failure classification.

# A local quadratic reservoir charging term is identity; a linear bias is not.
check(np.array_equal(QR_LIFT @ QR_LIFT, np.eye(4)),
      "local paired-reservoir charging Q_R^2 is an identity source term")
check(np.array_equal(V.T.conj() @ QR_LIFT @ V, -Z),
      "a source-independent reservoir bias descends to nonidentity link Z")

# A shared two-reservoir coupling descends to a cross-link pair operator.
V2 = np.kron(V, V)
shared_res_pair = np.kron(QR_LIFT, QR_LIFT)
check(np.array_equal(V2.T.conj() @ shared_res_pair @ V2, np.kron(Z, Z)),
      "shared reservoir charging descends to a nonidentity cross-link pair")

# The source-blind closed hold is an exact positive witness.  Active port,
# bias, or shared charging terms are exact counterexamples to unconditional
# FV-PURE inheritance.
check(np.array_equal(P_ENC @ H_FLIP, H_FLIP @ P_ENC) and
      np.any(commutator(P_ENC_LRE, H_BOUNDARY_LRE)),
      "closed encoded hold preserves equivalence while an active outer port does not")

result = {
    "lane": LANE.name,
    "status": "PASS",
    "checks": None,
    "encoded_link_dimension": 2,
    "physical_pair_dimension": 4,
    "q4_encoded_dimension": 16,
    "q4_physical_dimension": 256,
    "unitary_equivalence_on_fixed_total_charge": True,
    "m1_allocation": "co-located witness only; displaced reservoir requires connector current",
    "outer_port_boundary": "rhoR_dot - I_link_to_res + I_res_to_ext = 0",
    "fv_pure_fy": "nonidentity results preserved only modulo one common full-code identity/reference shift across encoded P+Q",
    "counterexamples": [
        "active reservoir-only outer port leaks from code",
        "reservoir linear bias descends to nonidentity Z",
        "shared reservoir charging descends to cross-link ZZ",
    ],
    "ceiling": "U(1) charge/current ancestry only; no spatial bond current, T0j, metric Ward, gravity, or G",
}
result["checks"] = f"{checks}/{checks}"
print("GA_RESULT_JSON", json.dumps(result, sort_keys=True))
print(f"SUMMARY {checks}/{checks} GA encoded charge-current checks passed")
print("CEILING U(1) link/reservoir/port current only; no spatial bond current, "
      "T0j, metric Ward, gravity, or G")
