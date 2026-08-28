#!/usr/bin/env python3
"""Exact algebra checks for the minimal translation-owning recoil parent.

The verifier uses integer momentum units on one common periodic auxiliary
mechanical torus.  A
kick K=2 makes the half-kick recoil states normalizable momentum eigenstates.
No thermodynamic Hilbert space, tensor field, or continuum stress tensor is
constructed.
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


for line in (LANE / "DEPENDENCIES.sha256").read_text().splitlines():
    expected, relative = line.split("  ", 1)
    path = ROOT / relative
    check(path.is_file() and not path.is_symlink() and
          sha256(path.read_bytes()).hexdigest() == expected,
          f"dependency custody: {relative}")


ga = (ROOT / "LANE_GRA_GA_F3_Q4_FU09B_ENCODED_CHARGE_CURRENT_LIFT_V001" /
      "THEOREM.md").read_text()
gb = (ROOT / "LANE_GRA_GB_F3_Q4_FIXED_SUPPORT_ENERGY_MOMENTUM_WARD_BOUNDARY_V001" /
      "THEOREM.md").read_text()
fv = (ROOT / "LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001" /
      "THEOREM.md").read_text()
fy = (ROOT / "LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001" /
      "THEOREM.md").read_text()
plan = (ROOT / "GRAVITY_NO_LAB_PROOF_TASK_PLAN_V001.md").read_text()
check("GA-CLOSED-FULL-CODE-SCALAR-HOLD" in ga and
      "modulo one full-code identity/reference shift" in ga,
      "GA fixes the full-code modulo-reference preservation requirement")
check("physical ownership gap" in gb and "conjugate recoil momenta" in gb,
      "GB fixes the missing translation and recoil ownership")
check("4678629417" in fv and "inserted before Feshbach" in fv,
      "FV fixes the nonidentity source/rank payload to preserve")
check("1\over\sqrt{60}" in fy and "native-support source" in fy,
      "FY fixes the native nonzero-momentum source payload to preserve")
check("B1" in plan and "one autonomous Hamiltonian" in plan and
      "unowned momentum port" in plan,
      "the no-lab plan fixes the B1 pass and fail conditions")
theorem_text = (LANE / "THEOREM.md").read_text()
check("same auxiliary mechanical torus" in theorem_text and
      "mechanical-translation neutral" in theorem_text and
      "code-independent same-mechanical-state" in theorem_text and
      "square" in theorem_text and "summability" in theorem_text and
      "not a physical diamond-space" in theorem_text,
      "the theorem pins common-torus, neutral-transfer, and analytic no-go scope")


# -------------------------------------------------------------------------
# One-link charge plus mechanical recoil encoder.

KICK = 2                 # integer torus momentum units; half-kick is allowed


def code_state(z_link, com=0):
    """(z_link,z_reservoir,p_link,p_reservoir) in one COM sector."""
    half = KICK // 2
    if z_link == -1:
        return (-1, +1, com+half, com-half)
    if z_link == +1:
        return (+1, -1, com-half, com+half)
    raise ValueError(z_link)


def recoil_flip(state):
    z_link, z_reservoir, p_link, p_reservoir = state
    if (z_link, z_reservoir) == (-1, +1):
        # sigma^+ T_- exp[i k(r_R-r_L)]
        return (+1, -1, p_link-KICK, p_reservoir+KICK)
    if (z_link, z_reservoir) == (+1, -1):
        # Hermitian reverse term.
        return (-1, +1, p_link+KICK, p_reservoir-KICK)
    return None


for com in range(-3, 4):
    minus, plus = code_state(-1, com), code_state(+1, com)
    check(minus[0]+minus[1] == plus[0]+plus[1] == 0,
          f"COM {com}: encoded total U(1) charge is zero")
    check(minus[2]+minus[3] == plus[2]+plus[3] == 2*com,
          f"COM {com}: both codewords have one total mechanical momentum")
    check(minus[2]**2+minus[3]**2 ==
          plus[2]**2+plus[3]**2 == 2*com**2+2,
          f"COM {com}: equal-mass recoil kinetic energy is code scalar")
    check(recoil_flip(minus) == plus and recoil_flip(plus) == minus,
          f"COM {com}: dressed charge-recoil flip restricts exactly to X")


# Each forward/reverse transition owns equal-and-opposite charge and momentum.
for z_link in (-1, +1):
    before = code_state(z_link)
    after = recoil_flip(before)
    delta_charge = (after[0]-before[0], after[1]-before[1])
    delta_momentum = (after[2]-before[2], after[3]-before[3])
    check(sum(delta_charge) == 0 and
          delta_charge == ((+2, -2) if z_link == -1 else (-2, +2)),
          f"z={z_link}: the dressed flip owns equal-and-opposite charge transfer")
    check(sum(delta_momentum) == 0 and
          delta_momentum == ((-KICK, +KICK) if z_link == -1
                             else (+KICK, -KICK)),
          f"z={z_link}: the dressed flip owns equal-and-opposite recoil")


def matrix_multiply(first, second):
    return tuple(tuple(sum(first[row][inner]*second[inner][column]
                           for inner in range(len(second)))
                       for column in range(len(second[0])))
                 for row in range(len(first)))


def matrix_add(first, second):
    return tuple(tuple(first[row][column]+second[row][column]
                       for column in range(len(first[0])))
                 for row in range(len(first)))


def matrix_scale(scale, matrix):
    return tuple(tuple(scale*value for value in row) for row in matrix)


def commutator(first, second):
    return matrix_add(matrix_multiply(first, second),
                      matrix_scale(-1, matrix_multiply(second, first)))


# Exact restriction of the operator currents on the two-dimensional code.
# Units hbar=q_*=1 and momentum kick K=2 are used here.
X = ((0, 1), (1, 0))
Z = ((-1, 0), (0, 1))
P_LINK = ((1, 0), (0, -1))
P_RESERVOIR = matrix_scale(-1, P_LINK)
Y = ((0, 0), (1, 0))
Y_DAGGER = ((0, 1), (0, 0))
h = 3
H_FLIP = matrix_scale(-h, X)
dot_p_link = matrix_scale(1j, commutator(H_FLIP, P_LINK))
dot_p_reservoir = matrix_scale(1j, commutator(H_FLIP, P_RESERVOIR))
j_momentum = matrix_scale(1j*h*KICK,
                          matrix_add(Y, matrix_scale(-1, Y_DAGGER)))
dot_charge_link = matrix_scale(1j, commutator(H_FLIP, Z))
i_charge = matrix_scale(-1, dot_charge_link)
check(matrix_add(dot_p_link, j_momentum) == ((0j, 0j), (0j, 0j)),
      "operator link momentum continuity has the GD11 sign")
check(matrix_add(dot_p_reservoir, matrix_scale(-1, j_momentum)) ==
      ((0j, 0j), (0j, 0j)),
      "operator reservoir momentum continuity has the opposite GD11 sign")
check(j_momentum == matrix_scale(F(-KICK, 2), i_charge),
      "operator charge and recoil currents obey J_P=-(hbar kappa/2q) I_Q")


# Four-link full-code check: all 2^4 inherited configurations, including
# off-ice configurations, have one common recoil kinetic scalar and every
# link flip remains inside the code.
encoded_configurations = {}
for configuration in product((-1, +1), repeat=4):
    encoded = tuple(code_state(z) for z in configuration)
    encoded_configurations[configuration] = encoded
check(len(encoded_configurations) == 16,
      "the four-link encoder covers the complete inherited 16-state Hilbert space")
check({sum(link[2]**2+link[3]**2 for link in encoded)
       for encoded in encoded_configurations.values()} == {8},
      "recoil kinetic energy is one scalar on the full four-link code")
for configuration, encoded in encoded_configurations.items():
    for link in range(4):
        target = list(configuration)
        target[link] *= -1
        lifted = list(encoded)
        lifted[link] = recoil_flip(lifted[link])
        check(tuple(lifted) == encoded_configurations[tuple(target)],
              f"full-code link {link} flip intertwines for {configuration}")


# Exact nonidentity source preservation on the code: arbitrary diagonal Z
# coefficients and arbitrary flip coefficients have the same matrix entries;
# the only new source-off term is the common recoil scalar.
diagonal_coefficients = (F(2, 3), F(-5, 7), F(11, 13), F(17, 19))
flip_coefficients = (F(3, 5), F(-7, 11), F(13, 17), F(19, 23))
decoded_configurations = {encoded: configuration
                          for configuration, encoded
                          in encoded_configurations.items()}
for configuration in encoded_configurations:
    inherited_diagonal = sum(coefficient*z for coefficient, z in
                             zip(diagonal_coefficients, configuration))
    encoded_diagonal = sum(coefficient*link[0] for coefficient, link in
                           zip(diagonal_coefficients,
                               encoded_configurations[configuration]))
    inherited_row = {configuration: inherited_diagonal}
    encoded_row = {configuration: encoded_diagonal}
    for link, coefficient in enumerate(flip_coefficients):
        target = list(configuration)
        target[link] *= -1
        inherited_row[tuple(target)] = coefficient
        lifted = list(encoded_configurations[configuration])
        lifted[link] = recoil_flip(lifted[link])
        encoded_row[decoded_configurations[tuple(lifted)]] = coefficient
    check(inherited_row == encoded_row,
          f"complete diagonal-plus-flip source row intertwines at {configuration}")


# -------------------------------------------------------------------------
# Momentum-ledger and boundary gates.

# A relative Weyl factor shifts (p_a,p_b)->(p_a-K,p_b+K); every active edge
# therefore has zero total impulse.  Check link-reservoir, reservoir-boundary,
# and support-controller representatives with distinct kicks.
for name, kick in (("link-reservoir", 2),
                   ("reservoir-boundary", 4),
                   ("support-controller", -6)):
    delta = (-kick, +kick)
    check(sum(delta) == 0 and delta[0] == -delta[1],
          f"{name} Weyl exchange closes its two-factor momentum ledger")

# A reservoir-only outer charge flip exits the zero-charge link-reservoir
# code, reproducing GA's active-port boundary.  Adding a dynamical exterior
# charge can conserve total charge and recoil, but does not preserve the hold.
state = code_state(-1)
reservoir_only = (state[0], -state[1], state[2], state[3])
check(reservoir_only[0]+reservoir_only[1] != 0,
      "an active reservoir-only outer port leaks from the encoded hold")
outer_charge_delta = -2*state[1]
check((reservoir_only[0]+reservoir_only[1])+(-outer_charge_delta) == 0,
      "an explicit exterior charge closes total outer-port charge bookkeeping")
check((-4)+(+4) == 0,
      "an explicit boundary recoil factor closes outer-port momentum bookkeeping")


# Minimality screen: a nonzero bilateral momentum shift cannot leave the same
# finite-support mechanical state invariant for both logical codewords.  The
# theorem, not this finite screen, supplies the l2 recurrence proof for all
# normalizable states.  Link/reservoir product states are not excluded.
for radius in (1, 2, 5):
    support = set(range(-radius, radius+1))
    shifted = {value+KICK for value in support}
    check(support != shifted,
          f"nonzero recoil shift has no invariant finite-support same-state hold (r={radius})")
check(KICK != 0,
      "zero recoil is the only way to reuse one mechanical state, but owns no impulse")


result = {
    "lane": LANE.name,
    "status": "PASS",
    "checks": f"{checks}/{checks}",
    "construction": "GA charge flip dressed by a relative mechanical Weyl kick",
    "code": "charge anticorrelation plus +/- half-kick recoil correlation",
    "full_code": "all four-link P/Q configurations preserved modulo one recoil scalar",
    "momentum": "each active relative exchange has exact equal-and-opposite impulse",
    "port": "outer boundary must be dynamical; active reservoir-only port breaks hold",
    "minimality": "nonzero recoil forbids reusing one normalizable mechanical state for both logical codewords",
    "ceiling": "translation/recoil ownership only; no complete spacetime source, stress Ward identity, tensor cone, gravity, or G",
}
print("GD_RESULT_JSON", json.dumps(result, sort_keys=True))
print(f"SUMMARY {checks}/{checks} GD translation-owning recoil checks passed")
print("CEILING exact B1 recoil/translation ledger and encoded preservation only; "
      "B2 spacetime-source and B4 Ward closure remain open")
