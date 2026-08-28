#!/usr/bin/env python3
"""Independent hostile verifier for the frozen GB Ward-boundary lane."""

from contextlib import redirect_stdout
from hashlib import sha256
from io import StringIO
from pathlib import Path
import json
import runpy

import numpy as np


AUDIT = Path(__file__).resolve().parent
LANE = AUDIT.parent
ROOT = LANE.parent
FO_SCRIPT = (ROOT / "LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001" /
             "verify_finite_tt_four_point.py")
BS_ACTION = (ROOT / "LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001" /
             "MICRO_ACTION.md")
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
      "target custody freezes all ten GB core/seal files")
check(verify_hash_list(LANE / "DEPENDENCIES.sha256", ROOT) == 7,
      "all seven declared FO/FY/FZ/GA/BS dependencies replay")
check(verify_hash_list(LANE / "MANIFEST.sha256", LANE) == 8,
      "the frozen eight-file GB author manifest replays")
seal_lines = (LANE / "SEAL.sha256").read_text().splitlines()
manifest_hash, manifest_name = seal_lines[0].split("  ", 1)
check(manifest_name == "MANIFEST.sha256" and
      manifest_hash == digest(LANE / manifest_name),
      "the author seal owns the frozen GB manifest")


# Claim-surface gates enforce the negative theorem's exact logical type.
theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
result_json = json.loads((LANE / "RESULT.json").read_text())
theorem_flat = " ".join(theorem.split())
result_flat = " ".join(result.split())
self_flat = " ".join(self_audit.split())

check("does not prove that no enlarged or continuum completion could possess"
      in theorem_flat and "derivability boundary" in theorem_flat,
      "translation theorem is underdetermination, not generator nonexistence")
check("not overclaimed as a theorem that every possible enlarged generator is "
      "nonlocal" in self_flat,
      "self-audit retains the discrete-logarithm nonselection ceiling")
check("not a derived spatial vector component" in theorem_flat and
      "not a continuum energy-momentum Ward identity" in theorem_flat,
      "term energy current is not promoted to spatial stress flux")
check("A charge current is not thereby `T^{0j}`" in theorem_flat and
      "U(1) charge current from GA is not identified with `T^{0j}`" in
      self_flat,
      "energy and U(1) continuity do not manufacture T0j")
check("neither proves nor falsifies local spatial-momentum conservation" in
      theorem_flat and "undecided rather than failed" in result_flat,
      "missing local stress identity remains undecided")
check("required ownership ledger, not a result" in theorem_flat and
      "no T0j" in result_json["ceiling"],
      "momentum/recoil completion and machine-readable ceiling are bounded")


# Load the finite parent only; all term-current work below uses a sparse exact
# operator implementation independent of the author's dense-matrix method.
with redirect_stdout(StringIO()):
    fo = runpy.run_path(str(FO_SCRIPT))

states = tuple(fo["states"])
state_index = dict(fo["state_index"])
patterns = tuple(fo["ring_patterns"])
hexagons = tuple(tuple(row) for row in fo["hexagons"])
dimension = len(states)
check((fo["CELL_COUNT"], fo["VERTEX_COUNT"], len(fo["edges"]),
       len(patterns), dimension) == (30, 60, 120, 120, 180),
      "independent replay uses the exact FO Z30/180-state component")


def op_clean(operator):
    return {key: value for key, value in operator.items() if value != 0}


def op_add(*operators):
    result_op = {}
    for operator in operators:
        for key, value in operator.items():
            result_op[key] = result_op.get(key, 0) + value
    return op_clean(result_op)


def op_scale(scale, operator):
    return op_clean({key: scale * value for key, value in operator.items()})


def op_product(first, second):
    second_rows = {}
    for (row, column), value in second.items():
        second_rows.setdefault(row, []).append((column, value))
    result_op = {}
    for (row, inner), left in first.items():
        for column, right in second_rows.get(inner, ()):
            key = (row, column)
            result_op[key] = result_op.get(key, 0) + left * right
    return op_clean(result_op)


def op_commutator(first, second):
    return op_add(op_product(first, second),
                  op_scale(-1, op_product(second, first)))


def hermitian_real(operator):
    return all(operator.get((column, row), 0) == value
               for (row, column), value in operator.items())


terms = []
for mask, first, second in patterns:
    operator = {}
    for row, state in enumerate(states):
        if (state & mask) in (first, second):
            operator[(row, state_index[state ^ mask])] = -1
    terms.append(operator)

check(all(hermitian_real(term) for term in terms),
      "all 120 independently rebuilt ring terms are Hermitian")
H = op_add(*terms)
fo_H = np.asarray(fo["hamiltonian"])
expected_H = {(int(row), int(column)): int(fo_H[row, column])
              for row, column in zip(*np.nonzero(fo_H))}
check(H == expected_H and len(H) == 840,
      "sparse term sum exactly recovers all 420 Hermitian H6 edges")

# Every term obeys [H,h_a]+[h_a,H-h_a]=0.  The commutators are kept without
# the common i/hbar so the verification is exact integer arithmetic.
nonzero_currents = []
for index, term in enumerate(terms):
    hdot_over_i = op_commutator(H, term)
    complement = op_add(H, op_scale(-1, term))
    jout_over_i = op_commutator(term, complement)
    if jout_over_i:
        nonzero_currents.append(index)
    if op_add(hdot_over_i, jout_over_i):
        raise AssertionError(f"term {index} violates current sign")
    if any(jout_over_i.get((column, row), 0) != -value
           for (row, column), value in jout_over_i.items()):
        raise AssertionError(f"term {index} current is not anti-Hermitian/i")
check(len(nonzero_currents) > 0,
      "all 120 exact sign identities hold and some term currents are nonzero")

ring_masks = tuple(pattern[0] for pattern in patterns)
disjoint = []
for first in range(120):
    for second in range(first + 1, 120):
        if ring_masks[first] & ring_masks[second]:
            continue
        disjoint.append((first, second))
        if op_commutator(terms[first], terms[second]):
            raise AssertionError("edge-disjoint ring terms failed to commute")
check(len(disjoint) == 5700,
      "all 5,700 edge-disjoint ring pairs have exactly zero current")

alpha = nonzero_currents[0]
beta = next(index for index in range(120)
            if op_commutator(terms[alpha], terms[index]))
pair_current = op_commutator(terms[alpha], terms[beta])
check(pair_current == op_scale(-1,
      op_commutator(terms[beta], terms[alpha])),
      "pairwise energy-current orientation is exactly antisymmetric")

# Repartition over two necessarily overlapping rings.  The alternate support
# is their bounded union (at most eleven edges), not an arbitrary long-range
# operator, while the attributed pair current changes to zero.
alpha_prime = op_add(terms[alpha], terms[beta])
beta_prime = {}
check(op_add(alpha_prime, beta_prime) ==
      op_add(terms[alpha], terms[beta]) and pair_current and
      not op_commutator(alpha_prime, beta_prime),
      "a bounded two-ring repartition fixes H but changes attributed current")
check((ring_masks[alpha] & ring_masks[beta]) != 0 and
      bin(ring_masks[alpha] | ring_masks[beta]).count("1") <= 11,
      "the hostile repartition witness remains on one overlapping-ring union")

# Choose a geometric term region by one shared physical edge rather than the
# author's first-twelve index set.
edge_zero_region = tuple(index for index, mask in enumerate(ring_masks)
                         if mask & 1)
h_region = op_add(*(terms[index] for index in edge_zero_region))
h_outside = op_add(H, op_scale(-1, h_region))
region_dot = op_commutator(H, h_region)
region_out = op_commutator(h_region, h_outside)
check(edge_zero_region and not op_add(region_dot, region_out),
      "an independently selected edge-star region obeys boundary continuity")

identity = {(index, index): 7 for index in range(dimension)}
check(op_commutator(identity, H) == {} and
      op_commutator(op_add(H, identity), terms[alpha]) ==
      op_commutator(H, terms[alpha]),
      "a common time-independent identity leaves every commutator current")


# -------------------------------------------------------------------------
# Finite translation, branch ambiguity, and fixed-support ownership.

permutation = tuple(int(value) for value in fo["translation_permutation"])
visited = set()
orbits = []
for seed in range(dimension):
    if seed in visited:
        continue
    orbit = []
    current = seed
    while current not in orbit:
        orbit.append(current)
        visited.add(current)
        current = permutation[current]
    orbits.append(tuple(orbit))
check(sorted(map(len, orbits)) == [30] * 6,
      "the actual translation is six free length-30 orbits")
check(all(H.get((permutation[row], permutation[column]), 0) == value
          for (row, column), value in H.items()),
      "H6 is exactly invariant under the global Z30 permutation")

# Independent spectral logarithms on one orbit.  Shift character 7 rather
# than the author's character 0 branch.
U30 = np.zeros((30, 30), dtype=complex)
for source in range(30):
    U30[(source + 1) % 30, source] = 1
powers = [np.linalg.matrix_power(U30, exponent) for exponent in range(30)]
projectors = []
for momentum in range(30):
    theta = 2 * np.pi * momentum / 30
    projectors.append(sum((np.exp(-1j * theta * exponent) * powers[exponent]
                           for exponent in range(30)),
                          np.zeros_like(U30)) / 30)
thetas = tuple(2 * np.pi * momentum / 30 for momentum in range(30))
K0 = sum((thetas[momentum] * projectors[momentum]
          for momentum in range(30)), np.zeros_like(U30))
K7 = K0 + 2 * np.pi * projectors[7]
U0 = sum((np.exp(1j * thetas[momentum]) * projectors[momentum]
          for momentum in range(30)), np.zeros_like(U30))
shifted = list(thetas)
shifted[7] += 2 * np.pi
U7 = sum((np.exp(1j * shifted[momentum]) * projectors[momentum]
          for momentum in range(30)), np.zeros_like(U30))
check(np.linalg.norm(K7 - K0) > 1 and
      np.max(np.abs(U0 - U30)) < 3e-13 and
      np.max(np.abs(U7 - U30)) < 3e-13,
      "two branch-distinct Hermitian generators exponentiate to the same U")
check(np.max(np.abs(K0 - K0.T.conj())) < 3e-13 and
      np.count_nonzero(np.abs(K0) > 1e-10) > 30,
      "one permitted global logarithm is Hermitian and cell-basis nonlocal")

# This is an ownership statement about the frozen response parent, not a
# theorem that no enlarged parent could ever carry momentum.
fo_source = FO_SCRIPT.read_text()
bs_source = BS_ACTION.read_text()
check(all(isinstance(state, int) for state in states) and
      "bond_vectors" in fo and "support_momentum" not in fo and
      "support_recoil" not in fo,
      "FO response states own incidence plus c-number geometry, not recoil")
check("support/recoil" in bs_source and
      "symbolic slot alone does not prove that completion" in bs_source,
      "BS owns only an uninstantiated support/recoil completion slot")
check("support_momentum" not in fo_source and "support_recoil" not in fo_source,
      "the frozen FO executable defines no hidden support momentum operator")


print(f"SUMMARY {checks}/{checks} independent hostile GB checks passed")
print("VERDICT PASS")
print("CEILING exact inherited term-energy continuity and global Z30 symmetry "
      "only; no unique local energy density, derived physical momentum/recoil, "
      "T0j, stress Ward packet, continuum locality, gravity, or G")
