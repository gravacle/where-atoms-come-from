#!/usr/bin/env python3
"""
ARM 3 (mechanical / soft-matter Kramers) -- arithmetic on the pinned numbers.

Nothing here is a fit. Every input is a number quoted in text or a caption of an
openly-loadable source, pinned in ARM3_KRAMERS.md with its access evidence.
The point of this file is to show, margin-free, WHICH of the two record laws a
given pinned datum can actually test, and which it merely re-derives.

Laws under test (C-69, C-70):
    LIFETIME : tau      = exp(E_b/kT) / (2 f0 cosh(dE/2kT))     [= 1/(g_up+g_dn)]
    STEADY   : <R>_ss   = tanh(dE/2kT)                          [= (t_dn - t_up)/(t_dn + t_up)]
    RATIO    : tau_up/tau_dn = exp(-dE/kT)                      [detailed balance]

Run:  python3 arm3_numbers.py > arm3_numbers.txt
"""

import math

out = []


def say(s=""):
    out.append(s)


def rule(t):
    say()
    say("=" * 78)
    say(t)
    say("=" * 78)


# ----------------------------------------------------------------------------
rule("1. THE TAUTOLOGY TEST -- what a bare pair of dwell times can and cannot do")
# ----------------------------------------------------------------------------
say("""
Given only tau_up and tau_dn from one trace, the program's STEADY law is an
IDENTITY, not a measurement. Define dE from the ratio (the only thing a bare
trace offers):
        dE/kT := ln(tau_dn / tau_up)
Then tanh(dE/2kT) == (tau_dn - tau_up)/(tau_dn + tau_up) algebraically, for ANY
two positive numbers. Demonstrated below on the Bercy-Bockelmann dwell times.
A pair of dwell times therefore tests C-70 ONLY IF dE is known INDEPENDENTLY --
from a Boltzmann-inverted potential, or from a control parameter (force x
distance) whose mapping to dE has been earned.
""")


def steady_check(label, t_up, t_dn):
    dE = math.log(t_dn / t_up)                      # in units of kT, by construction
    R_from_times = (t_dn - t_up) / (t_dn + t_up)    # <R>_ss from occupancies
    R_from_law = math.tanh(dE / 2.0)                # the C-70 law
    say(f"  {label}")
    say(f"    tau_up (unfolded) = {t_up:.6g} s   tau_dn (folded) = {t_dn:.6g} s")
    say(f"    dE/kT  inferred from the ratio      = {dE:.6f}")
    say(f"    <R>_ss from the dwell times         = {R_from_times:.9f}")
    say(f"    tanh(dE/2kT) from the C-70 law      = {R_from_law:.9f}")
    say(f"    residual                            = {abs(R_from_times - R_from_law):.3e}"
        "   <-- machine zero: an identity, not a test")
    say()
    return dE


say("Bercy & Bockelmann 2015 (NAR 43:9928, CC-BY), Fig. 7 caption, T = 29 C:")
dE_dna = steady_check("DNA10, extension where both states are ~equally occupied",
                      0.11, 0.25)
dE_rna = steady_check("RNA10, extension where both states are ~equally occupied",
                      1.77, 3.1)
say(f"  Inferred asymmetries: DNA10 dE = {dE_dna:.3f} kT, RNA10 dE = {dE_rna:.3f} kT")
say("  At T = 302.15 K (29 C), kT = 26.04 meV, so:")
kT_meV = 1.380649e-23 * 302.15 / 1.602176634e-19 * 1e3
say(f"    kT              = {kT_meV:.4f} meV")
say(f"    DNA10 dE        = {dE_dna * kT_meV:.4f} meV")
say(f"    RNA10 dE        = {dE_rna * kT_meV:.4f} meV")
say("""
  SEMANTICS: these dE values are DERIVED FROM the dwell times. They are not an
  independent measurement of the asymmetry, so they close a circle. They are
  still useful as a CONSISTENCY FLOOR: any pipeline that reads these traces and
  does NOT reproduce the identity to machine precision is broken.
""")

# ----------------------------------------------------------------------------
rule("2. RONDIN 2017 -- the one pinned system whose measured object IS C-69's tau")
# ----------------------------------------------------------------------------
say("""
Rondin et al. measure R by fitting an exponential to the autocorrelation of the
binary population operator h_A. Their own text: "R = R_AC + R_CA", the sum of the
forward and backward rate constants. That is exactly the program's

        tau = 1/(g_up + g_dn)

i.e. the CORRECTED C-69 object, not the naive exp(E_b/kT)/f0. The correspondence
is structural, not numerical luck: both are the relaxation eigenvalue of a
two-state Liouvillian.

Pinned numbers (arXiv:1703.07699v2, main text + Fig. 3):
    T          = 300 K            (stated: "T = 300 K is the temperature of the gas")
    U_A        ~ 4 k_B T          (barrier from well A, "we measure the energy barriers")
    U_C        ~ 5 k_B T          (barrier from well C)
    R_peak     ~ 2.4e3 s^-1       (turnover maximum, READ FROM Fig. 3 -- digitized)
    turnover   at P_gas ~ 1200-1600 Pa (predicted), observed max near that range
""")
U_A, U_C, T_rondin = 4.0, 5.0, 300.0
dE_rondin = U_C - U_A           # well C deeper by (U_C - U_A) if the saddle is shared
say(f"  Well-depth asymmetry implied by the two barriers: dE = U_C - U_A = {dE_rondin:.2f} kT")
say(f"  (shared saddle B is stated: both barriers are measured to the SAME saddle point)")
say(f"  2 cosh(dE/2kT)      = {2*math.cosh(dE_rondin/2):.6f}")
say(f"  tanh(dE/2kT)        = {math.tanh(dE_rondin/2):.6f}   <-- C-70 prediction for <R>_ss")
R_peak = 2.4e3
say(f"  tau at turnover     = 1/R_peak = {1/R_peak:.4e} s   (R_peak digitized from Fig. 3)")
say("""
  WHAT THIS CAN TEST: the cosh factor. With U_A and U_C both stated in kT and a
  shared saddle, dE is INDEPENDENT of the dwell times -- the circle of section 1
  is broken. This is the only pinned system in Arm 3 where the measured quantity
  is literally 1/(g_up+g_dn) AND the asymmetry is independently stated.

  WHAT IT CANNOT TEST: f0. The swept parameter is gas pressure, i.e. the damping
  Gamma -- which is precisely the prefactor. The paper's whole point is that the
  prefactor is NOT constant: it rises with Gamma, peaks, and falls. C-69 assumes
  a single f0. Rondin is therefore a DIRECT WARNING to C-69's convention. The
  barriers are held fixed across the sweep, so rate ratio == prefactor ratio:
  their measured points run from ~1.2e3 to ~2.4e3 s^-1 (Fig. 3, DIGITIZED), a
  factor of ~2 in f0 across the measured decades, and the two asymptotic branches
  they plot fall away without bound outside it. Any C-69 anchor that treats f0 as
  a material constant inherits an unstated damping-regime assumption.
""")

# ----------------------------------------------------------------------------
rule("3. ZIJLSTRA 2020 -- the text-anchored consistency pair")
# ----------------------------------------------------------------------------
say("""
Zijlstra/Schuler PRL 125, 146001 (2020), Fig. 1c caption, one example particle:
    autocorrelation relaxation time      tau_relax = 4.4 +/- 0.1 s
    inverse sum of the rate coefficients 1/(k1+k2) = 3.9 +/- 0.3 s
    T = 295 K (stated in Methods, Eq. 1)
""")
t_relax, t_relax_e = 4.4, 0.1
t_invsum, t_invsum_e = 3.9, 0.3
diff = abs(t_relax - t_invsum)
comb = math.hypot(t_relax_e, t_invsum_e)
say(f"  difference          = {diff:.3f} s")
say(f"  combined 1-sigma    = {comb:.3f} s")
say(f"  discrepancy         = {diff/comb:.2f} sigma")
say("""
  SEMANTICS: this is the same identity C-69 rests on -- the population relaxation
  eigenvalue equals the inverse sum of the two rate constants -- measured two
  independent ways on the same trace and agreeing at ~1.6 sigma. It is a
  MEASURED check of the instrument's central assumption, on a real thermal
  system, at a stated temperature. It is one point, not a sweep.
""")

# ----------------------------------------------------------------------------
rule("4. THE FORCE -> dE MAPPING -- earned in one modality, not in the other")
# ----------------------------------------------------------------------------
say("""
For a hairpin held at CONSTANT FORCE F (passive force clamp), the free-energy
difference between folded and unfolded is

        dG(F) = dG_0 - F * dx        (dx = extension change on unfolding)

and F is a genuine, externally-set tilt of a two-state landscape. The mapping
dE := dG(F) is EARNED, subject to two named corrections:
    (a) handle/linker compliance -- dx is the molecular extension change, and the
        DNA handles stretch too; dx must be the handle-corrected value.
    (b) the folded/unfolded states must be the only two occupied states.

For a hairpin held at CONSTANT TRAP SEPARATION (extension clamp / "high
stiffness" modality), the trap itself stores elastic energy when the hairpin
opens, so the tilt is NOT F*dx and the effective landscape is the molecule PLUS
the trap spring. There, dE := dG(F) is NOT EARNED and must not be asserted.

This distinction is load-bearing for the Lyons/Woodside deposit, which ships BOTH
modalities in separately-named folders ("CF" and "HS"), and for Bercy &
Bockelmann, whose measurement is explicitly constant-extension.

Numerical scale of the tilt, to show the regime: a 10 nm extension change at
1 pN of force, in units of kT at 295 K.
""")
kT_J = 1.380649e-23 * 295.0
for F_pN in (1.0, 5.0, 10.0, 16.0):
    for dx_nm in (10.0, 20.0):
        E = F_pN * 1e-12 * dx_nm * 1e-9 / kT_J
        say(f"    F = {F_pN:5.1f} pN, dx = {dx_nm:4.1f} nm  ->  F*dx = {E:8.3f} kT")
say("""
  CONSEQUENCE: at the forces these experiments run (5-16 pN) and the extension
  changes hairpins show (10-20 nm), the force term is TENS of kT. The asymmetry
  dE swings through zero across a range of well under 1 pN. This is why hairpin
  work is always reported near F_1/2: a fraction of a piconewton is the whole
  dynamic range of dE. Any pinned force must therefore carry its calibration
  uncertainty, or the dE it implies is worthless.
""")

# ----------------------------------------------------------------------------
rule("5. FLOOR: what a candidate must supply to be EXECUTABLE against C-69/C-70")
# ----------------------------------------------------------------------------
say("""
    C-70 (steady/ratio):  tau_up, tau_dn  AND  dE from a source that is not the
                          dwell times themselves. Without the second, section 1
                          shows the check is vacuous.
    C-69 (lifetime):      additionally E_b and f0, each independently determined.
                          A fitted f0 makes the law an identity of the fit (this
                          is the exact failure the solidity review found).
    Both:                 a stated bath temperature with an uncertainty, and a
                          statement of what "the two states" are physically.
""")

print("\n".join(out))
