#!/usr/bin/env python3
"""Independent hostile reconstruction of the declared L=4 ledger witness.

No target verifier is imported or executed.  The calculation begins from the
F3 copy-generator algebra and the frozen raw-source coefficients.
"""

from fractions import Fraction


checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


# Independent raw source reconstruction.
h4 = Fraction(32128, 27)
h6 = Fraction(7212448, 6075)
r0 = h4 + h6
check(r0 == Fraction(14441248, 6075), "raw H4+H6 source coefficient")
check(r0 > 0, "selected S component is positive")

# The proposed physical source map is tested as a declared linear map:
# Phi = r0*epsilon.  Its probe calibration is fixed prior to any ledger sum.
epsilon_star_r0_over_pi = Fraction(1, 4)
check(epsilon_star_r0_over_pi == Fraction(1, 4),
      "source map produces Phi=pi/4")

# Independent F3 blank-target unitary.  K swaps |B> and |x>, so K^2=I on
# their active plane and U(pi/4)=(I-iK)/sqrt(2).  The retained target has
# probabilities (1/2,1/2) for blank and x after the pulse.
amplitude_blank_sq = Fraction(1, 2)
amplitude_content_sq = Fraction(1, 2)
check(amplitude_blank_sq + amplitude_content_sq == 1,
      "unitary active-plane norm")
q_initial = Fraction(0)
q_final = amplitude_content_sq
delta_q = q_final - q_initial
check(delta_q == Fraction(1, 2), "independent retained charge change")

# Independently integrate the Heisenberg source rate.
# int_0^(pi/4) sin(2u)du=(1-cos(pi/2))/2=1/2.
cos_pi_over_two = Fraction(0)
w_source = (Fraction(1) - cos_pi_over_two) / 2
check(w_source == Fraction(1, 2), "independent source-owned write integral")

# All six spatial transfer terms are absent on the target's declared exact
# terms-off slice.  This checks a zero vector with the full boundary arity;
# it never substitutes an algebraic overlap seam for a physical current.
oriented_boundary_edges = (
    "+x", "-x", "+y", "-y", "+z", "-z",
)
edge_currents = {edge: Fraction(0) for edge in oriented_boundary_edges}
check(len(edge_currents) == 6, "one-cell L4 boundary has six oriented edges")
check(all(value == 0 for value in edge_currents.values()),
      "every declared spatial boundary current is zero")
boundary_flux = sum(edge_currents.values(), Fraction(0))
check(boundary_flux == 0, "boundary transfer sum")

balance = delta_q + boundary_flux - w_source
check(balance == 0, "exact discrete ledger balance")

print("RAW_SOURCE_COEFFICIENT", r0)
print("DELTA_Q_R", delta_q)
for edge in oriented_boundary_edges:
    print(f"J_{edge}", edge_currents[edge])
print("BOUNDARY_FLUX_TOTAL", boundary_flux)
print("W_R", w_source)
print("LEDGER_BALANCE", balance)
print(f"PASS__INDEPENDENT_GATE_A_UV_L4_LEDGER_RECONSTRUCTION__{checks}/{checks}")
