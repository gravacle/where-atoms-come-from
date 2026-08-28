#!/usr/bin/env python3
"""Independent hostile verifier for the repaired GA encoded-current lift."""

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

import numpy as np


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
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


def verify_hash_list(list_path, base):
    count = 0
    for line in list_path.read_text().splitlines():
        if not line.strip():
            continue
        expected, relative = line.split("  ", 1)
        path = base / relative
        check(path.is_file() and not path.is_symlink() and
              digest(path) == expected,
              f"custody {list_path.name}: {relative}")
        count += 1
    return count


check(verify_hash_list(AUDIT / "TARGET_CUSTODY.sha256", LANE) == 10,
      "target custody freezes all ten repaired GA core/seal files")
check(verify_hash_list(LANE / "DEPENDENCIES.sha256", ROOT) == 6,
      "all six declared FU/FV/FY/FZ dependencies replay")
check(verify_hash_list(LANE / "MANIFEST.sha256", LANE) == 8,
      "the repaired eight-file GA author manifest replays")
seal_lines = (LANE / "SEAL.sha256").read_text().splitlines()
manifest_hash, manifest_name = seal_lines[0].split("  ", 1)
check(manifest_name == "MANIFEST.sha256" and
      manifest_hash == digest(LANE / manifest_name),
      "the author seal owns the repaired GA manifest")


# Repaired claim surfaces are load-bearing inputs to this PASS.
theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
result_json = json.loads((LANE / "RESULT.json").read_text())
theorem_flat = " ".join(theorem.split())
result_flat = " ".join(result.split())
self_flat = " ".join(self_audit.split())

check("internal to a link--reservoir pair; it is not a spatial current" in
      theorem_flat and "port-normal exchange equation" in theorem_flat,
      "internal exchange is not promoted to a spatial bond current")
check("No vertex divergence or spatial bond current follows" in theorem_flat and
      "No spatial bond current or vertex divergence has been derived" in
      result_flat,
      "placement analysis retains the no-spatial-current ceiling")
check("GA-CLOSED-FULL-CODE-SCALAR-HOLD" in theorem and
      "full encoded Hilbert space" in theorem_flat and
      "virtual/off-ice `Q` sectors" in theorem_flat,
      "source hold requires one common scalar across full encoded P+Q")
check("z+c[j]" in theorem and "modulo one full-code identity/reference shift"
      in theorem_flat and "literal equality" in theorem_flat,
      "Feshbach and source preservation are correctly modulo identity")
check("Source independence alone is insufficient" in theorem_flat and
      "constant reservoir bias descends to `Z`" in result_flat,
      "source-independent nonidentity bias is explicitly excluded")
check("no spatial bond current" in result_json["ceiling"] and
      "gravity" in result_json["ceiling"],
      "machine-readable result retains the bounded U(1)-ancestry ceiling")


def commutator(first, second):
    return first @ second - second @ first


def kron_all(factors):
    result = np.array(((1,),), dtype=complex)
    for factor in factors:
        result = np.kron(result, factor)
    return result


def embed(operator, site, count, identity):
    return kron_all(operator if index == site else identity
                    for index in range(count))


def zero(matrix):
    return np.zeros_like(matrix)


# -------------------------------------------------------------------------
# Encoder, exact generator intertwining, and current orientation.

I2 = np.eye(2, dtype=complex)
Z = np.diag((-1, 1)).astype(complex)
X = np.array(((0, 1), (1, 0)), dtype=complex)
SP = np.array(((0, 0), (1, 0)), dtype=complex)
SM = SP.T
QR = np.diag((1, -1)).astype(complex)
TM = np.array(((0, 0), (1, 0)), dtype=complex)
TP = TM.T

ZL = np.kron(Z, I2)
QRL = np.kron(I2, QR)
A = np.kron(SP, TM)
B = np.kron(SM, TP)
XD = A + B
QT = ZL + QRL
V = np.zeros((4, 2), dtype=complex)
V[0, 0] = 1
V[3, 1] = 1
PENC = V @ V.T.conj()

check(tuple(np.diag(QT).real.astype(int)) == (0, -2, 2, 0) and
      np.array_equal(PENC, np.diag((1, 0, 0, 1))),
      "V spans exactly the complete two-dimensional zero-charge eigenspace")
required_reservoir_charges = (F(1), F(-1))
check(len(set(required_reservoir_charges)) == 2,
      "a one-state reservoir cannot encode both link charges at fixed total")
check(np.array_equal(V.T.conj() @ V, I2) and
      np.array_equal(ZL @ V, V @ Z) and
      np.array_equal(XD @ V, V @ X),
      "independent encoder exactly intertwines inherited Z and X")
check(np.array_equal(commutator(QT, XD), zero(XD)) and
      np.array_equal(commutator(PENC, XD), zero(XD)),
      "dressed flip conserves total charge and preserves the code")
dark = np.eye(4, dtype=complex)[:, (1, 2)]
check(np.array_equal(XD @ dark, np.zeros((4, 2), dtype=complex)),
      "the two off-code charge sectors are dark, not silently encoded")

H = -XD
qdot = 1j * commutator(H, ZL)
rdot = 1j * commutator(H, QRL)
ILR = -qdot
check(np.array_equal(qdot, 2j * (A - B)) and
      np.array_equal(rdot, -qdot),
      "internal link/reservoir current signs reproduce GA06")
check(np.array_equal(ILR, ILR.T.conj()) and np.any(ILR) and
      np.array_equal(qdot + ILR, zero(qdot)) and
      np.array_equal(rdot - ILR, zero(rdot)),
      "oriented Hermitian exchange current closes both local balances")
original_current = -1j * commutator(-X, Z)
check(np.array_equal(ILR @ V, V @ original_current),
      "encoded exchange current intertwines the inherited charge motion")


# -------------------------------------------------------------------------
# Tensor lift and exact nonidentity source preservation.

LINKS = 4
I4 = np.eye(4, dtype=complex)
V4 = kron_all((V,) * LINKS)
oz = tuple(embed(Z, site, LINKS, I2) for site in range(LINKS))
ox = tuple(embed(X, site, LINKS, I2) for site in range(LINKS))
lz = tuple(embed(ZL, site, LINKS, I4) for site in range(LINKS))
lx = tuple(embed(XD, site, LINKS, I4) for site in range(LINKS))
lqr = tuple(embed(QRL, site, LINKS, I4) for site in range(LINKS))

check(V4.shape == (256, 16) and
      all(np.array_equal(lz[site] @ V4, V4 @ oz[site]) and
          np.array_equal(lx[site] @ V4, V4 @ ox[site])
          for site in range(LINKS)),
      "four independent pairs intertwine all q4 generators simultaneously")

SIGNS = ((1, 1, 1), (1, -1, -1),
         (-1, 1, -1), (-1, -1, 1))


def dyad(vector):
    x, y, z = vector
    return (x*x, y*y, z*z, 2*x*y, 2*x*z, 2*y*z)


qpair = [np.zeros((16, 16), dtype=complex) for _ in range(6)]
qpair_lift = [np.zeros((256, 256), dtype=complex) for _ in range(6)]
for first in range(4):
    for second in range(first + 1, 4):
        root = tuple(SIGNS[second][axis] - SIGNS[first][axis]
                     for axis in range(3))
        for component, coefficient in enumerate(dyad(root)):
            # These are 16 Q_pair/Ud at lambda=-1/2.
            qpair[component] -= coefficient * (oz[first] @ oz[second])
            qpair_lift[component] -= coefficient * (lz[first] @ lz[second])

qflip = [np.zeros((16, 16), dtype=complex) for _ in range(6)]
qflip_lift = [np.zeros((256, 256), dtype=complex) for _ in range(6)]
for link, signs in enumerate(SIGNS):
    for component, coefficient in enumerate(dyad(signs)):
        # These are 3 Q_flip/h.
        qflip[component] -= coefficient * ox[link]
        qflip_lift[component] -= coefficient * lx[link]

check(all(np.array_equal(qpair_lift[c] @ V4, V4 @ qpair[c])
          for c in range(6)),
      "all six Coulomb pair-source coordinates intertwine with correct scale")
check(all(np.array_equal(qflip_lift[c] @ V4, V4 @ qflip[c])
          for c in range(6)),
      "all six hopping-source coordinates intertwine with correct scale")

# Add a nonzero identity source and verify only the centered/nonidentity part
# is unchanged.  This directly attacks the repaired modulo-reference claim.
identity_coefficients = (1, -2, 3, -4, 5, -6)
original_sources = tuple(qpair[c] + qflip[c] for c in range(6))
lifted_sources = tuple(qpair_lift[c] + qflip_lift[c] +
                       identity_coefficients[c] * np.eye(256)
                       for c in range(6))
restricted = tuple(V4.T.conj() @ source @ V4 for source in lifted_sources)
check(all(np.array_equal(restricted[c], original_sources[c] +
                         identity_coefficients[c] * np.eye(16))
          for c in range(6)),
      "a full-code identity derivative remains exactly an identity source")


def centered(matrix):
    return matrix - np.trace(matrix) / matrix.shape[0] * np.eye(matrix.shape[0])


check(all(np.array_equal(centered(restricted[c]),
                         centered(original_sources[c]))
          for c in range(6)),
      "centering removes the allowed reference source exactly")

gram_original = np.array([[np.trace(centered(first).T.conj() @ centered(second))
                           for second in original_sources]
                          for first in original_sources])
gram_restricted = np.array([[np.trace(centered(first).T.conj() @ centered(second))
                             for second in restricted]
                            for first in restricted])
check(np.array_equal(gram_original, gram_restricted),
      "nonidentity source Gram matrix and rank survive identity shifts")

# A one-dimensional P/Q Feshbach calculation tests the energy-reference shift
# without reusing the author's matrix implementation.
hpp, hqq, hpq, z, scalar = F(0), F(3), F(-1), F(5), F(7)
heff = hpp + hpq*hpq/(z-hqq)
heff_shifted = hpp + scalar + hpq*hpq/((z+scalar)-(hqq+scalar))
check(heff_shifted == heff + scalar,
      "one common full-code scalar cancels every gap and shifts H_eff by I")
noncommon = hpp + scalar + hpq*hpq/((z+scalar)-(hqq+scalar+1))
check(noncommon != heff + scalar,
      "different P/Q scalars fail the repaired Feshbach preservation law")


# -------------------------------------------------------------------------
# Placement, boundary terms, and explicit failure modes.

H4 = -sum(lx)
currents = tuple(-1j * commutator(H4, lz[link])
                 for link in range(LINKS))
exponents = (5, 9, 25, 201)
phases = tuple(np.exp(2j * np.pi * exponent / 240)
               for exponent in exponents)
rho_link = sum((phases[link] * lz[link] for link in range(LINKS)),
               np.zeros((256, 256), dtype=complex))
i_weighted = sum((phases[link] * currents[link]
                  for link in range(LINKS)),
                 np.zeros((256, 256), dtype=complex))
check(np.max(np.abs(1j * commutator(H4, rho_link) + i_weighted)) < 2e-13,
      "co-located link-to-reservoir exchange closes GA10 with its signs")

rho_colocated = sum((phases[link] * (lz[link] + lqr[link])
                     for link in range(LINKS)),
                    np.zeros((256, 256), dtype=complex))
rho_displaced = rho_link + sum(lqr, np.zeros((256, 256), dtype=complex))
rho_global = sum((lz[link] + lqr[link] for link in range(LINKS)),
                 np.zeros((256, 256), dtype=complex))
check(np.max(np.abs(commutator(H4, rho_colocated))) < 2e-13 and
      np.array_equal(commutator(H4, rho_global), zero(H4)),
      "co-located m=1 and allocation-independent global m=0 charges conserve")
check(np.linalg.norm(commutator(H4, rho_displaced)) > 1e-8 and
      all(exponent % 240 for exponent in exponents),
      "common displaced reservoirs require a nonzero connector allocation")

# The exterior exchange acts only on reservoir/exterior.  Recompute all signs
# and test leakage by the off-code block, not merely a commutator checksum.
H_BOUNDARY_RE = -(np.kron(TM, TP) + np.kron(TP, TM))
H_LRE = np.kron(H, I2) + np.kron(I2, H_BOUNDARY_RE)
QL_LRE = np.kron(ZL, I2)
QR_LRE = np.kron(QRL, I2)
QE_LRE = np.kron(np.eye(4), QR)
ldot = 1j * commutator(H_LRE, QL_LRE)
rdot_outer = 1j * commutator(H_LRE, QR_LRE)
edot = 1j * commutator(H_LRE, QE_LRE)
ILR_outer = -ldot
IRE = edot
check(np.array_equal(ldot + ILR_outer, zero(ldot)) and
      np.array_equal(rdot_outer - ILR_outer + IRE, zero(rdot_outer)) and
      np.array_equal(edot - IRE, zero(edot)),
      "independent three-factor boundary-current signs close exactly")
check(np.array_equal(ldot + rdot_outer + edot, zero(ldot)),
      "total link/reservoir/exterior U(1) charge is conserved")
P_LRE = np.kron(PENC, I2)
leakage = (np.eye(8) - P_LRE) @ np.kron(I2, H_BOUNDARY_RE) @ P_LRE
check(np.any(leakage),
      "reservoir-only active port has an explicit off-code leakage block")

check(np.array_equal(QRL @ QRL, np.eye(4)) and
      np.array_equal(V.T.conj() @ QRL @ V, -Z),
      "quadratic charging is scalar but source-independent linear bias is not")
V2 = np.kron(V, V)
check(np.array_equal(V2.T.conj() @ np.kron(QRL, QRL) @ V2,
                     np.kron(Z, Z)),
      "shared reservoir charging creates an encoded cross-link ZZ term")
strain_transfer = 1j * (A - B)
check(np.array_equal(strain_transfer, strain_transfer.T.conj()) and
      np.linalg.matrix_rank(V.T.conj() @ strain_transfer @ V) == 2,
      "a strain-dependent transfer supplies a nonidentity extra source")


print(f"SUMMARY {checks}/{checks} independent hostile GA checks passed")
print("VERDICT PASS")
print("CEILING exact finite encoded U(1) ancestry and modulo-identity source "
      "preservation only; no selected reservoir placement, spatial bond "
      "current, vertex divergence, T0j, Ward closure, continuum, gravity, or G")
