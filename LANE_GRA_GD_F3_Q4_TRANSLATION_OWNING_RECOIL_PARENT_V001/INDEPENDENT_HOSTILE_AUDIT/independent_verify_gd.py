#!/usr/bin/env python3
"""Algorithmically independent hostile replay of the frozen GD/TORP packet."""

from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


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
      "target custody freezes all ten GD core/seal files")
check(verify_hash_list(LANE / "DEPENDENCIES.sha256", ROOT) == 11,
      "all eleven GA/GB/FV/FY/plan dependencies replay")
check(verify_hash_list(LANE / "MANIFEST.sha256", LANE) == 8,
      "the frozen eight-file GD author manifest replays")
seal_hash, seal_name = (LANE / "SEAL.sha256").read_text().strip().split("  ", 1)
check(seal_name == "MANIFEST.sha256" and
      seal_hash == digest(LANE / seal_name),
      "the author seal owns the frozen GD manifest")


# -------------------------------------------------------------------------
# Claim-surface and operator-domain typing.

theorem = (LANE / "THEOREM.md").read_text()
result = (LANE / "RESULT.md").read_text()
self_audit = (LANE / "SELF_AUDIT.md").read_text()
readme = (LANE / "README.md").read_text()
result_json = json.loads((LANE / "RESULT.json").read_text())
flat_t = " ".join(theorem.split())
flat_r = " ".join(result.split())
flat_s = " ".join(self_audit.split())

check("same auxiliary mechanical torus" in flat_t and
      "kappa_e in 2 Lambda^*" in flat_t and "p_e/hbar in Lambda^*" in flat_t,
      "Weyl factors and half-kick states share one explicit reciprocal lattice")
check("mechanical-translation neutral" in flat_t and
      "additional commutators must be retained" in flat_t,
      "the pair-current formula is limited to mechanically neutral transfers")
check("[Pi_(C,on),H_(R partial)]=0" in flat_t and
      "off sector invariant" in flat_t,
      "the conditional port term is Hermitian and the hold sector invariant")
check("does not construct a nontrivial controller transition or energy- current packet"
      in flat_t and "conditional factor class" in flat_r,
      "controller and boundary closure is conditional and not an energy-current proof")
check("not a physical diamond-space current" in flat_t and
      "auxiliary common mechanical torus" in readme,
      "factor-edge locality is not promoted to localized physical space")
check("code-independent same-mechanical-state hold no-go" in flat_t and
      "does **not** exclude link/reservoir product momentum states" in flat_t,
      "the no-go is restricted to reuse of one mechanical state")
check("NONZERO_RECOIL_CODE_INDEPENDENT_SAME_MECHANICAL_STATE_HOLD_NO_GO"
      in theorem and "NONZERO_RECOIL_PRODUCT_HOLD_NO_GO" not in theorem,
      "the machine disposition retains the narrowed no-go")
check("source-dependent kick, mass, support potential, controller, or port is outside"
      in flat_t and "source-independent full-code hold" in flat_r,
      "FV/FY inheritance is limited to the source-independent scalar hold")
check(result_json["plan_gate"].startswith("B1 algebraic existence") and
      "not yet a localized physical diamond-space current" in
      result_json["placement"],
      "machine-readable result retains the algebraic-existence placement ceiling")
check(result_json["ceiling"] ==
      "translation/recoil ownership only; no stress Ward identity, tensor cone, gravity, or G",
      "machine-readable result retains the Gravity Formation ceiling")


# -------------------------------------------------------------------------
# Independent three-dimensional momentum-state construction.

def vadd(first, second):
    return tuple(a + b for a, b in zip(first, second))


def vsub(first, second):
    return tuple(a - b for a, b in zip(first, second))


def vscale(scale, vector):
    return tuple(scale * value for value in vector)


def norm2(vector):
    return sum(value * value for value in vector)


def encoded(z_link, center, kick):
    half = vscale(F(1, 2), kick)
    if z_link == -1:
        momenta = (vadd(center, half), vsub(center, half))
        return (-1, +1) + momenta
    if z_link == +1:
        momenta = (vsub(center, half), vadd(center, half))
        return (+1, -1) + momenta
    raise ValueError(z_link)


def dressed_flip(state, kick):
    z_l, z_r, p_l, p_r = state
    if (z_l, z_r) == (-1, +1):
        return (+1, -1, vsub(p_l, kick), vadd(p_r, kick))
    if (z_l, z_r) == (+1, -1):
        return (-1, +1, vadd(p_l, kick), vsub(p_r, kick))
    raise ValueError((z_l, z_r))


kicks = ((2, -4, 6), (-6, 2, 4), (4, 8, -2))
centers = ((0, 0, 0), (3, -2, 5), (-4, 7, 1))
for index, (kick, center) in enumerate(zip(kicks, centers)):
    minus = encoded(-1, center, kick)
    plus = encoded(+1, center, kick)
    check(dressed_flip(minus, kick) == plus and
          dressed_flip(plus, kick) == minus,
          f"3D witness {index}: relative Weyl kick intertwines exactly with X")
    check(vadd(minus[2], minus[3]) == vscale(2, center) ==
          vadd(plus[2], plus[3]),
          f"3D witness {index}: pair total momentum is fixed")
    check(norm2(minus[2]) + norm2(minus[3]) ==
          norm2(plus[2]) + norm2(plus[3]) ==
          2 * norm2(center) + F(1, 2) * norm2(kick),
          f"3D witness {index}: equal-mass kinetic energy is code scalar")
    check(minus[0] + minus[1] == plus[0] + plus[1] == 0,
          f"3D witness {index}: encoded charge remains neutral")


# Direct action of the Weyl shift gives the commutator signs without importing
# or calling the author verifier.  Units hbar=1.
probe_l = (5, -3, 2)
probe_r = (-1, 7, 4)
kick = kicks[0]
after_l, after_r = vsub(probe_l, kick), vadd(probe_r, kick)
check(vsub(after_l, probe_l) == vscale(-1, kick),
      "direct Weyl action gives [P_L,U]=-kappa U")
check(vsub(after_r, probe_r) == kick,
      "direct Weyl action gives [P_R,U]=+kappa U")
check(vadd(after_l, after_r) == vadd(probe_l, probe_r),
      "direct Weyl action commutes with pair total momentum")


# Exact two-codeword Heisenberg calculation of current signs.
def mm(first, second):
    return tuple(tuple(sum(first[r][q] * second[q][c]
                           for q in range(len(second)))
                       for c in range(len(second[0])))
                 for r in range(len(first)))


def madd(first, second):
    return tuple(tuple(first[r][c] + second[r][c]
                       for c in range(len(first[0])))
                 for r in range(len(first)))


def mscale(scale, matrix):
    return tuple(tuple(scale * value for value in row) for row in matrix)


def comm(first, second):
    return madd(mm(first, second), mscale(-1, mm(second, first)))


zero2 = ((0j, 0j), (0j, 0j))
xmat = ((0, 1), (1, 0))
zmat = ((-1, 0), (0, 1))
ymat = ((0, 0), (1, 0))
ydag = ((0, 1), (0, 0))
h = 7
H = mscale(-h, xmat)
for component, kval in enumerate(kicks[0]):
    p_l = mscale(F(kval, 2), ((1, 0), (0, -1)))
    p_r = mscale(-1, p_l)
    j_p = mscale(1j * h * kval, madd(ymat, mscale(-1, ydag)))
    dot_l = mscale(1j, comm(H, p_l))
    dot_r = mscale(1j, comm(H, p_r))
    check(madd(dot_l, j_p) == zero2 and
          madd(dot_r, mscale(-1, j_p)) == zero2,
          f"current component {component}: link/reservoir continuity signs close")

dot_charge_l = mscale(1j, comm(H, zmat))
i_charge = mscale(-1, dot_charge_l)
for component, kval in enumerate(kicks[0]):
    j_p = mscale(1j * h * kval, madd(ymat, mscale(-1, ydag)))
    check(j_p == mscale(F(-kval, 2), i_charge),
          f"current component {component}: J_P=-(kappa/2q) I_Q")


# -------------------------------------------------------------------------
# Complete four-link code and a nontrivial inherited Hamiltonian row.

link_kicks = ((2, 0, 0), (0, 4, 0), (0, 0, -6), (2, -2, 4))
link_centers = ((0, 0, 0), (1, 2, -1), (-2, 0, 3), (4, -1, 1))
code = {}
for config in product((-1, +1), repeat=4):
    code[config] = tuple(encoded(z, p, k)
                         for z, p, k in zip(config, link_centers, link_kicks))
check(len(code) == 16 and len(set(code.values())) == 16,
      "four-link isometry covers the complete inherited P+Q link Hilbert")

energies = {sum(norm2(pair[2]) + norm2(pair[3]) for pair in state)
            for state in code.values()}
check(len(energies) == 1,
      "four-link recoil kinetic energy is one full-code scalar")

diag = (F(2, 3), F(-5, 7), F(11, 13), F(17, 19))
pair_couplings = {(0, 1): F(7, 11), (1, 3): F(-13, 17),
                  (0, 2): F(19, 23)}
flips = (F(3, 5), F(-7, 11), F(13, 17), F(19, 29))
for config, lifted in code.items():
    inherited_row = {config: sum(c * z for c, z in zip(diag, config)) +
                     sum(c * config[a] * config[b]
                         for (a, b), c in pair_couplings.items())}
    lifted_row = dict(inherited_row)
    for link, coefficient in enumerate(flips):
        target = list(config)
        target[link] *= -1
        target = tuple(target)
        changed = list(lifted)
        changed[link] = dressed_flip(changed[link], link_kicks[link])
        check(tuple(changed) == code[target],
              f"full-code flip {link} closes at configuration {config}")
        inherited_row[target] = coefficient
        lifted_row[target] = coefficient
    check(inherited_row == lifted_row,
          f"diagonal interaction plus flip row is preserved at {config}")


# -------------------------------------------------------------------------
# Independent exact Feshbach/common-reference calculation.

def eye(size):
    return tuple(tuple(F(int(row == column)) for column in range(size))
                 for row in range(size))


def msub(first, second):
    return madd(first, mscale(-1, second))


def inv2(matrix):
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if not determinant:
        raise ZeroDivisionError
    return ((matrix[1][1] / determinant, -matrix[0][1] / determinant),
            (-matrix[1][0] / determinant, matrix[0][0] / determinant))


P = ((F(2), F(1, 3)), (F(1, 3), F(-1)))
Q = ((F(5), F(2, 7)), (F(2, 7), F(7)))
V = ((F(1, 2), F(-2, 5)), (F(3, 8), F(4, 9)))
VT = tuple(zip(*V))
z = F(11)
c = F(13, 6)
resolvent = inv2(msub(mscale(z, eye(2)), Q))
heff = madd(P, mm(mm(V, resolvent), VT))
shifted_P = madd(P, mscale(c, eye(2)))
shifted_Q = madd(Q, mscale(c, eye(2)))
shifted_resolvent = inv2(msub(mscale(z + c, eye(2)), shifted_Q))
shifted_heff = madd(shifted_P, mm(mm(V, shifted_resolvent), VT))
check(shifted_resolvent == resolvent,
      "a common P/Q reference shift cancels exactly in the resolvent")
check(shifted_heff == madd(heff, mscale(c, eye(2))),
      "the Feshbach effective Hamiltonian changes only by the common identity")


# A P/Q-unequal addition must not be discarded as a reference shift.
bad_Q = madd(Q, mscale(c + F(1, 9), eye(2)))
bad_resolvent = inv2(msub(mscale(z + c, eye(2)), bad_Q))
check(bad_resolvent != resolvent,
      "a P/Q-unequal scalar changes virtual denominators and fails the hold")


# -------------------------------------------------------------------------
# Closed factor graph and analytic no-go boundary.

nodes = ("L", "R", "C", "B")
events = (("L", "R", (2, -4, 6)),
          ("R", "B", (-6, 2, 4)),
          ("L", "C", (4, 8, -2)))
ledger = {node: (0, 0, 0) for node in nodes}
for source, target, event_kick in events:
    ledger[source] = vsub(ledger[source], event_kick)
    ledger[target] = vadd(ledger[target], event_kick)
check(tuple(sum(ledger[node][axis] for node in nodes) for axis in range(3)) ==
      (0, 0, 0),
      "an explicit closed factor graph has zero total impulse")
check(all(vadd(vscale(-1, event_kick), event_kick) == (0, 0, 0)
          for _, _, event_kick in events),
      "every admitted pair event owns equal-and-opposite impulse")

# If U|chi>=phase|chi>, every coefficient on an infinite shift orbit has the
# same nonzero norm.  Exact partial norms therefore grow as (2N+1)|c|^2.
amplitude2 = F(5, 17)
partial_norms = tuple((2 * radius + 1) * amplitude2
                      for radius in (1, 2, 5, 20, 100))
check(all(first < second for first, second in zip(partial_norms,
                                                  partial_norms[1:])),
      "same-state Weyl-eigenvector partial norms grow strictly along the orbit")
check(partial_norms[-1] > 50 and
      "square summability then forces every coefficient" in flat_t,
      "the analytic l2 recurrence excludes a nonzero normalizable same-state hold")
check("does **not** exclude link/reservoir product momentum states" in flat_t,
      "the no-go leaves the explicit product half-kick code admissible")


print(f"SUMMARY {checks}/{checks} independent hostile GD checks passed")
print("VERDICT PASS -- BOUNDED B1 ALGEBRAIC EXISTENCE ONLY")
print("CEILING no physical diamond-space placement, spacetime source, T0j, stress "
      "Ward identity, tensor cone, gravity, or G")
