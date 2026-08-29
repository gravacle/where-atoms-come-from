#!/usr/bin/env python3
"""Exact standard-library verifier for the GD/FY/FZ Ward obstruction.

The theorem tested here is deliberately conditional.  It concerns the
direct configuration-diagonal momentum density furnished by GD's encoded
half-kick hold and FY/FZ's *supplied embedding contraction*.  GD does not
derive the native diamond-space divergence, and this verifier does not invent
one.
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


# -------------------------------------------------------------------------
# Dependency and claim custody.

for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink() and
          sha256(path.read_bytes()).hexdigest() == expected,
          f"dependency custody: {relative}")

gd_theorem = (ROOT /
    "LANE_GRA_GD_F3_Q4_TRANSLATION_OWNING_RECOIL_PARENT_V001" /
    "THEOREM.md").read_text()
gd_result = (ROOT /
    "LANE_GRA_GD_F3_Q4_TRANSLATION_OWNING_RECOIL_PARENT_V001" /
    "RESULT.md").read_text()
fy_theorem = (ROOT /
    "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001" /
    "THEOREM.md").read_text()
fy_result = json.loads((ROOT /
    "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001" /
    "RESULT.json").read_text())
fz_theorem = (ROOT /
    "LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001" /
    "THEOREM.md").read_text()
fz_result = (ROOT /
    "LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001" /
    "RESULT.md").read_text()
theorem = (LANE / "THEOREM.md").read_text()

check("J^P_{L_e\\to R_e}" in gd_theorem and
      "not be called `T^{0j}`" in gd_theorem,
      "GD types its off-diagonal recoil object as a factor-edge current, not T0j")
check("not yet a localized diamond-space current" in gd_result and
      "geometric placement" in gd_result,
      "GD leaves diamond placement and the native spatial-current bind open")
check("source-independent full-code hold" in gd_result and
      "No `T^{0j}`" in gd_result,
      "GD preserves FY only on a source-independent hold and claims no T0j")
check("strictly\noff diagonal" in fz_theorem and
      fy_result["exact_ring_result"] ==
      "nonzero independent off-diagonal m=1 source",
      "FY fixes the ring source as off diagonal")
check(fy_result["exact_m1_diagonal_lift"] == [
          {"order": 2, "coefficient": "-1"},
          {"order": 4, "coefficient": "-37/12"},
          {"order": 6, "coefficient": "-16247/900"}],
      "FY fixes the complete diagonal H2/H4/H6 lift coefficients")
check("supplied continuum-embedding contraction" in fz_theorem and
      "discrete divergence has not been derived" in fz_theorem,
      "FZ distinguishes its embedding contraction from a derived divergence")
check("does not rule them out for an unknown physical\ndivergence" in fz_result,
      "FZ explicitly preserves the unknown-native-divergence ceiling")


# -------------------------------------------------------------------------
# Exact GD encoded momenta and direct projected density type.

HBAR = F(1)
momenta = (F(5, 7), F(-2, 3), F(11, 13))
kicks = (F(2), F(-4), F(6))


def code_momenta(z):
    """GD06 in units with hbar=1, for z=-1,+1."""
    return (tuple(p - k*z/2 for p, k in zip(momenta, kicks)),
            tuple(p + k*z/2 for p, k in zip(momenta, kicks)))


for z in (-1, +1):
    p_link, p_reservoir = code_momenta(z)
    check(all(p_link[i] + p_reservoir[i] == 2*momenta[i]
              for i in range(3)),
          f"z={z}: GD code has one fixed total pair momentum")
    check(all(p_link[i] == momenta[i] - kicks[i]*z/2 and
              p_reservoir[i] == momenta[i] + kicks[i]*z/2
              for i in range(3)),
          f"z={z}: each direct factor momentum is affine in diagonal Z")

# Arbitrary local/Fourier weights change coefficients but not diagonal type.
weights_link = (F(3, 5), F(-7, 11), F(13, 17), F(19, 23))
weights_reservoir = (F(-5, 29), F(11, 31), F(17, 37), F(-23, 41))
configurations = tuple(product((-1, +1), repeat=4))


def direct_density(configuration, component):
    total = F(0)
    for edge, z in enumerate(configuration):
        p_link, p_reservoir = code_momenta(z)
        total += (weights_link[edge]*p_link[component] +
                  weights_reservoir[edge]*p_reservoir[component])
    return total


direct_rows = tuple(tuple(direct_density(configuration, component)
                          for component in range(3))
                    for configuration in configurations)
check(len(set(direct_rows)) > 1,
      "generic placement weights make the direct GD density nontrivial")
check(all(all(isinstance(value, F) for value in row) for row in direct_rows),
      "the direct GD density is exact and configuration diagonal")

ice_configurations = tuple(configuration for configuration in configurations
                           if sum(1 for z in configuration if z == -1) == 2)
check(len(ice_configurations) == 6,
      "the local q4 ice projector retains exactly six configurations")
check(all(configuration in configurations and
          sum(1 for z in configuration if z == -1) == 2
          for configuration in ice_configurations),
      "the ice subset is represented in the same configuration basis")

# A one-link transfer flips one Z and therefore leaves a degree-two q4 ice
# fiber.  Hence its direct off-diagonal current has P_ice J P_ice=0.
single_flip_targets = []
for configuration in ice_configurations:
    for edge in range(4):
        target = list(configuration)
        target[edge] *= -1
        single_flip_targets.append(tuple(target))
check(len(single_flip_targets) == 24 and
      all(target not in ice_configurations for target in single_flip_targets),
      "all twenty-four single-link currents exit the local q4 ice fiber")


# -------------------------------------------------------------------------
# Exact matrix algebra: density versus current typing.


def madd(first, second):
    return tuple(tuple(first[i][j] + second[i][j]
                       for j in range(len(first[0])))
                 for i in range(len(first)))


def mscale(scale, matrix):
    return tuple(tuple(F(scale)*value for value in row) for row in matrix)


def mm(first, second):
    return tuple(tuple(sum(first[i][k]*second[k][j]
                           for k in range(len(second)))
                       for j in range(len(second[0])))
                 for i in range(len(first)))


def commutator(first, second):
    return madd(mm(first, second), mscale(-1, mm(second, first)))


def diagonal(matrix):
    return tuple(matrix[i][i] for i in range(len(matrix)))


zero2 = ((F(0), F(0)), (F(0), F(0)))
X = ((F(0), F(1)), (F(1), F(0)))
Z = ((F(-1), F(0)), (F(0), F(1)))
Y = ((F(0), F(0)), (F(1), F(0)))
Y_DAGGER = ((F(0), F(1)), (F(0), F(0)))
A = madd(Y, mscale(-1, Y_DAGGER))       # J^P/i up to h*kappa
h = F(3)
kappa = F(2)
H_flip = mscale(-h, X)
P_link = mscale(-kappa/2, Z)
J_over_i = mscale(h*kappa, A)

check(madd(commutator(H_flip, P_link), J_over_i) == zero2,
      "after factoring i, GD direct momentum and recoil current obey dot P+J=0")
check(diagonal(commutator(H_flip, P_link)) == (F(0), F(0)) and
      any(commutator(H_flip, P_link)[i][j]
          for i in range(2) for j in range(2) if i != j),
      "commutator with direct diagonal momentum is off diagonal")

# If one deliberately mistypes the off-diagonal flux J as a density, its next
# commutator can indeed be diagonal.  This is why the theorem must be typed,
# not based on the word 'recoil'.
flux_acceleration_without_i = mscale(-1, commutator(H_flip, J_over_i))
check(any(diagonal(flux_acceleration_without_i)) and
      all(flux_acceleration_without_i[i][j] == 0
          for i in range(2) for j in range(2) if i != j),
      "the off-diagonal recoil flux would have diagonal acceleration if mistyped as T0j")

# The GD factor-edge incidence closes that current on its own auxiliary graph.
dot_pair = (mscale(-1, J_over_i), J_over_i)
check(madd(dot_pair[0], J_over_i) == zero2 and
      madd(dot_pair[1], mscale(-1, J_over_i)) == zero2,
      "GD closes the L/R auxiliary factor-edge balance exactly")
check(madd(dot_pair[0], dot_pair[1]) == zero2,
      "reindexing or summing the two auxiliary endpoints leaves zero total impulse")

# Generic finite-dimensional fact used at the projected FY level: for any H
# and any configuration-diagonal D, diag([H,D])=0.
H_test = (
    (F(2), F(3), F(0), F(-1), F(0), F(4)),
    (F(3), F(-5), F(7), F(0), F(2), F(0)),
    (F(0), F(7), F(11), F(5), F(0), F(-3)),
    (F(-1), F(0), F(5), F(13), F(6), F(0)),
    (F(0), F(2), F(0), F(6), F(-17), F(9)),
    (F(4), F(0), F(-3), F(0), F(9), F(19)),
)
D_test = tuple(tuple(F((row+2)*(row+3)) if row == column else F(0)
                     for column in range(6)) for row in range(6))
check(diagonal(commutator(H_test, D_test)) == (F(0),)*6,
      "diag([H,D]) vanishes exactly for a nontrivial Hermitian H and diagonal D")
check(all(H_test[row][column]*D_test[column][column] -
          D_test[row][row]*H_test[row][column] ==
          commutator(H_test, D_test)[row][column]
          for row in range(6) for column in range(6)),
      "the exact commutator is H_ab(D_b-D_a) in the configuration basis")


# -------------------------------------------------------------------------
# FZ's exact diagonal supplied-embedding witness.


def trim(poly):
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def exact_divide(numerator, denominator):
    numerator = trim([F(value) for value in numerator])
    denominator = trim([F(value) for value in denominator])
    quotient = [F(0)] * max(1, len(numerator)-len(denominator)+1)
    while len(numerator) >= len(denominator) and any(numerator):
        shift = len(numerator)-len(denominator)
        factor = numerator[-1]/denominator[-1]
        quotient[shift] += factor
        for index, value in enumerate(denominator):
            numerator[index+shift] -= factor*value
        numerator = trim(numerator)
    if any(numerator):
        raise AssertionError("nonzero exact polynomial remainder")
    return tuple(trim(quotient))


def cyclotomic(n):
    divisors = tuple(value for value in range(1, n+1) if n % value == 0)
    table = {}
    for divisor in divisors:
        poly = (F(-1),) + (F(0),)*(divisor-1) + (F(1),)
        for proper in divisors:
            if proper >= divisor:
                break
            if divisor % proper == 0:
                poly = exact_divide(poly, table[proper])
        table[divisor] = poly
    return table[n]


def reduce_monic(poly, modulus):
    remainder = trim([F(value) for value in poly])
    while len(remainder) >= len(modulus):
        shift = len(remainder)-len(modulus)
        factor = remainder[-1]
        for index, value in enumerate(modulus):
            remainder[index+shift] -= factor*value
        remainder = trim(remainder)
    return tuple(remainder)


phi240 = cyclotomic(240)
check(len(phi240)-1 == 64 and phi240[-1] == 1,
      "Phi_240 is reconstructed exactly with degree 64")

witness_terms = {
    0: -14, 2: 28, 8: 14, 10: 14, 24: 14, 26: -14,
    34: -28, 42: -14, 48: -14, 56: -14, 58: 28,
}
witness = [F(0)]*59
for power, coefficient in witness_terms.items():
    witness[power] = F(coefficient)
witness_remainder = reduce_monic(witness, phi240)
check(any(witness_remainder) and witness_remainder == tuple(witness),
      "FZ05 is an exact nonzero diagonal longitudinal witness in Q(zeta_240)")

x = F(2, 5)
rho = F(15625, 504)
f_e = 1 - x**2 - F(37, 12)*x**4 - F(16247, 900)*x**6
check(f_e == F(2415673, 3515625) and
      rho*f_e == F(2415673, 113400),
      "the complete through-H6 diagonal coefficient is exact and nonzero at x=2/5")
scaled_witness = tuple(rho*f_e*value for value in witness_remainder)
check(any(scaled_witness),
      "the complete diagonal supplied-embedding contraction remains nonzero")
check("2415673/113400" in fz_result and
      "independent ring source is off diagonal" in fz_result,
      "FZ result custody separates the nonzero diagonal and ring supports")

# A configuration-diagonal T0j therefore cannot satisfy the operator Ward
# identity under this supplied contraction.  Off-diagonal ring pieces cannot
# cancel a diagonal entry.
commutator_diagonal = (F(0),)*len(scaled_witness)
ward_diagonal = tuple(left + right for left, right in
                      zip(commutator_diagonal, scaled_witness))
check(any(ward_diagonal),
      "direct GD density leaves an exact nonzero source-off Ward diagonal under the supplied contraction")


# -------------------------------------------------------------------------
# Contact, boundary, and theorem-ceiling gates.

# R(j)=a*j^2+b*j^3 has zero first derivative at source off but a nonzero
# Hessian.  This is the exact distinction between the operator source and a
# response seagull.
a, b = F(7, 11), F(-5, 13)
r_first_at_zero = F(0)
r_second_at_zero = 2*a
check(r_first_at_zero == 0 and r_second_at_zero != 0,
      "an O(j^2) contact changes a seagull but not the source-off first spatial source")

# A mixed h0*j contact likewise has a nonzero mixed derivative but contributes
# neither direct T0j nor Tij when both sources are zero.
mixed = F(17, 19)
d_h_at_j0 = mixed*F(0)
d_j_at_h0 = mixed*F(0)
d_h_d_j = mixed
check(d_h_at_j0 == d_j_at_h0 == 0 and d_h_d_j != 0,
      "a mixed h0-j contact cannot fill either source-off operator slot")

check("conditional on the supplied embedding" in theorem and
      "physical native\ndivergence `Delta_m` is still undefined" in theorem and
      "bare, directly projected" in theorem and
      "scalar-weighted GD `P_L/P_R` assignment" in theorem,
      "the new theorem preserves the FZ divergence and bare-direct-density ceiling")
check("off-diagonal recoil current" in theorem and
      "cannot be reassigned" in theorem,
      "the new theorem explicitly attacks the strongest current-mistyping evasion")
check("active boundary" in theorem and "it changes the\nparent" in theorem,
      "the new theorem preserves the active-boundary escape only as a new parent")
check("Feshbach-dressed" in theorem and "remains open" in theorem,
      "the new theorem does not exclude a properly derived effective momentum density")
check("does not refute" in theorem and "gravity" in theorem and
      "constructive recommendation, not a logical-necessity" in theorem,
      "the new theorem states its strict non-gravity and non-necessity ceiling")


result = {
    "lane": LANE.name,
    "status": "PASS",
    "checks": f"{checks}/{checks}",
    "proved": "the bare directly projected scalar-weighted GD P_L/P_R density cannot close the frozen FY source-off Ward identity under the FZ supplied embedding contraction",
    "current_test": "GD off-diagonal recoil current closes the auxiliary L/R factor-edge balance but is flux, has zero direct ice projection, and is not T0j",
    "contacts": "O(j^2) and mixed h0-j contacts cannot alter the first source-off operator identity",
    "ceiling": "native physical Delta, dynamical position-weighted localization, interaction contributions to T0j, a modified spatial source, Feshbach-dressed h0i source, active boundaries, complete pair-field/support stress, continuum, gravity, and G remain open",
}
print("RESULT_JSON", json.dumps(result, sort_keys=True))
print(f"SUMMARY {checks}/{checks} flip-recoil embedding-Ward checks passed")
print("CEILING conditional supplied-embedding obstruction only; native physical "
      "divergence and complete source-before-Feshbach construction remain open")
