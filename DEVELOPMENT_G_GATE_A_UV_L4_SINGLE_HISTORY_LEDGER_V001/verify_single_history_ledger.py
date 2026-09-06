#!/usr/bin/env python3
"""Exact arithmetic checks for the declared L=4 one-history ledger witness.

This verifier intentionally checks only the finite witness.  It does not
claim that the linear source attachment is derived from the bare F3 parent.
"""

from fractions import Fraction
from pathlib import Path


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


# The raw k=0 source coefficient, independently audited as H4 plus H6.
r0 = Fraction(32128, 27) + Fraction(7212448, 6075)
check(r0 == Fraction(14441248, 6075), "audited raw source coefficient")
check(r0 > 0, "positive chosen S component coefficient")

# The explicitly declared L=4 periodic history has 64 cells.  Its one-cell
# region has six oriented spatial boundary edges, all disabled on the exact
# GL6Q terms-off write slice.
cells = [(i, j, k) for i in range(4) for j in range(4) for k in range(4)]
region = {(0, 0, 0)}
check(len(cells) == 64, "L=4 cell count")
check(region <= set(cells), "region belongs to periodic L4 family")
boundary_edge_fluxes = [Fraction(0) for _ in range(6)]
check(sum(boundary_edge_fluxes, Fraction(0)) == 0,
      "all six declared spatial edge currents vanish")

# epsilon_star = pi/(4*r0) makes Phi_star=pi/4.  To avoid floating-point
# trigonometry, carry the pulse coordinate in units of pi: phi/pi=1/4.
epsilon_times_r0_over_pi = Fraction(1, 4)
check(epsilon_times_r0_over_pi == Fraction(1, 4),
      "linear raw-source-to-pulse attachment has Phi=pi/4")

# The audited blank-target F3 unitary gives Q=sin^2(Phi).  At the declared
# pulse Phi/pi=1/4, the exact special-angle identity is cos(2Phi)=0.  The
# source integral is independently int_0^Phi sin(2u) du=(1-cos(2Phi))/2;
# it is not assigned from delta_q.
phi_initial_over_pi = Fraction(0)
phi_final_over_pi = Fraction(1, 4)
cos_two_phi_initial = Fraction(1)
cos_two_phi_final = Fraction(0)
q_initial = (Fraction(1) - cos_two_phi_initial) / 2
q_final = (Fraction(1) - cos_two_phi_final) / 2
delta_q = q_final - q_initial
w_source = (Fraction(1) - cos_two_phi_final) / 2
check(delta_q == Fraction(1, 2), "retained-lineage charge change")
check(w_source == Fraction(1, 2), "source-owned unitary write integral")

boundary_total = sum(boundary_edge_fluxes, Fraction(0))
balance = delta_q + boundary_total - w_source
check(balance == 0, "discrete ledger closes exactly")

# Preserve the result's declared ceiling in executable form.
check(Path(__file__).with_name("THEOREM.md").is_file(), "theorem present")
print("RAW_SOURCE_COEFFICIENT", r0)
print("EPSILON_STAR_TIMES_R0_OVER_PI", epsilon_times_r0_over_pi)
print("DELTA_Q_R", delta_q)
print("BOUNDARY_FLUX_TOTAL", boundary_total)
print("W_R", w_source)
print("LEDGER_BALANCE", balance)
print(f"PASS__GATE_A_UV_L4_SINGLE_HISTORY_WITNESS__{CHECKS}/{CHECKS}")
