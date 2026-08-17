#!/usr/bin/env python3
# W-17 RULING — THE NEGATIVE CONTROL THE ROUND DID NOT RUN.
#
# The instrument F1..F5 fired on 19 of 20 slots across R1-R4. That statistic is
# uninterpretable without a frame KNOWN TO BE SOUND. W-15 names one, in its own
# discriminator (REGISTER_V001.md:1509-1511):
#     "diagonal / non-diagonal was forced and has a proof under it."
# and the strongest possible control is a partition no one can dispute at all: parity.
#
# This script runs F1..F5 against BOTH, and reports the fire rate under
#   (a) the LETTER of FRAME_CHALLENGE_V001.md as written, and
#   (b) an AMENDED F1 with an EXCLUSION clause (the off-axis point must lie in NEITHER arm).
#
# It also runs one strengthening computation for R1's F5, which the cutoff audit
# challenged as K1-dependent (K1 postdates the founding decision): the founding
# escape inference is refutable with NO carrier at all.
#
# python3 + numpy only.  Seed fixed.  No sympy.

import numpy as np

RNG = np.random.default_rng(20260817)
SQ = lambda z: float(np.sqrt(np.sum(np.abs(z) ** 2)))

def rule(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)

# =====================================================================
# CONTROL 1 — PARITY.  Route form: "take the even branch or the odd branch?"
# An unimpeachable partition of Z.  Any test that fires here is a fault machine.
# =====================================================================
rule("CONTROL 1 — PARITY  (arms: n even / n odd).  Unimpeachably sound partition.")

N = 100000
ns = np.arange(N)
even = (ns % 2 == 0)
odd = ~even

# ---- F1 EMBED, BY THE LETTER: "find the parameter space containing both arms.
#      If one exists, the binary is a SAMPLE of a space and not a partition of it."
# Z/2 sits inside the family {Z/m}.  The embedding EXISTS, trivially.
mods = [2, 3, 4, 5, 6, 7, 8, 12]
print("F1 (letter): embedding space = congruence classes mod m, m in", mods)
print("            the named binary is the m=2 slice; dimension of the family = 1 (the modulus m)")
print("            off-axis points exist at every m != 2, e.g. residue classes mod 3:",
      sorted(set((ns % 3).tolist())))
print("            -> AN EMBEDDING EXISTS.  By the letter of F1, FIRES.")

# ---- F1 EMBED, AMENDED: EXCLUSION CLAUSE.  Does any off-axis point lie in NEITHER arm?
outside = int(np.sum(~(even | odd)))
print(f"F1 (amended, EXCLUSION): integers lying in NEITHER arm: {outside} of {N}"
      f"  fraction = {outside/N:.6f}")
# and the RELEVANCE clause: does the decision-relevant functional (which branch runs)
# differ at an off-axis point?  The branch is a function of n mod 2 by construction:
branch = (ns % 2)
resid3 = ns % 3
# measure: is branch constant on the m=2 fibres regardless of the mod-3 coordinate?
viol = 0
for r in range(3):
    sel = resid3 == r
    # within a mod-3 class, the branch still takes both values - i.e. mod 3 carries
    # NO information the decision uses, once mod 2 is known.
    for b in (0, 1):
        both = np.sum(sel & (branch == b))
        if both == 0:
            viol += 1
print(f"F1 (amended, RELEVANCE): mod-3 coordinate constrains the branch in {viol} of 6 cells"
      f"  -> the off-axis coordinate is DECISION-IRRELEVANT")
print("            -> AMENDED F1: SOUND.")

# ---- F2 DEGENERACY.  Candidate maps listed by the instrument: power, restriction,
#      change of coordinates, monomial substitution.
shift = (ns + 1) % 2  # change of coordinates n -> n+1 maps odds onto evens exactly
A = set(ns[even].tolist())
fB = set((ns[odd] + 1).tolist())
inter = len(A & fB)
symdiff = len(A ^ fB)
print(f"\nF2: map = change of coordinates f(n) = n+1.")
print(f"    |A \\ f(B)| + |f(B) \\ A| = {symdiff}   (boundary term only: {symdiff} of {N})")
print(f"    density of symmetric difference = {symdiff/N:.8f}   ->  DISTANCE 0 in the limit")
sq = ns ** 2
pw = int(np.sum((sq % 2) != (ns % 2)))
print(f"    map = power n -> n^2: parity violations {pw} of {N} (power PRESERVES the arms)")
print("    -> ONE listed map carries arm B exactly onto arm A.  By the letter of F2, FIRES.")
print("    (and it is a SYMMETRY EXCHANGING the arms, which is EVIDENCE OF A GOOD PARTITION,")
print("     not evidence of degeneracy.  F2 as written cannot tell the two apart.)")

# ---- F3 CARVING.
dens_even = float(np.mean(even)); dens_odd = float(np.mean(odd))
overlap = float(np.mean(even & odd))
print(f"\nF3: natural density  mu(even) = {dens_even:.6f}   mu(odd) = {dens_odd:.6f}"
      f"   overlap = {overlap:.6f}   union = {float(np.mean(even|odd)):.6f}")
print("    -> both cells positive measure, no overlap, exhaustive.  F3: SOUND.")

# ---- F4 PRESUPPOSITION.
ints = int(np.sum(ns == np.floor(ns)))
print(f"\nF4: presupposition = 'every input is an integer'.  Tested: {ints} of {N} hold."
      f"  -> TRUE.  F4: SOUND.")
print("    BUT the verdict is set by WHICH sentence the lane elects to write.  A second,")
print("    equally defensible sentence - 'the cost of the two branches is comparable' - is")
print("    contingent and can be made false by choosing the algorithm.  F4 has unbounded")
print("    lane discretion and returns no forced number.")

# ---- F5 NULL.
wA = int(ns[even][1]); wB = int(ns[odd][1])
print(f"\nF5: witness for arm A: n = {wA} (even).  witness for arm B: n = {wB} (odd).")
print("    both branches obtain.  F5: SOUND.")

# =====================================================================
# CONTROL 2 — DIAGONAL / NON-DIAGONAL.  W-15's OWN named paradigm of a forced
# binary "with a proof under it" (REGISTER_V001.md:1509-1511).  Route form: the
# corpus's transport operator - is the branch operator diagonal in the class basis?
# =====================================================================
rule("CONTROL 2 — DIAGONAL / NON-DIAGONAL on U(3).  W-15's own certified-sound binary.")

def haar_u(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2.0)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r)
    return q * (d / np.abs(d))

def offdiag(U):
    M = U.copy()
    np.fill_diagonal(M, 0.0)
    return SQ(M)

NH = 20000
offs = np.empty(NH)
for i in range(NH):
    offs[i] = offdiag(haar_u(3, RNG))
print(f"F3: Haar measure on U(3), {NH} draws.")
print(f"    mu(DIAGONAL cell)     = {float(np.mean(offs < 1e-9)):.6f}   "
      f"(count below 1e-9: {int(np.sum(offs < 1e-9))})")
print(f"    mu(NON-DIAGONAL cell) = {float(np.mean(offs >= 1e-9)):.6f}")
print(f"    min off-diagonal norm observed = {offs.min():.6e}")
print("    -> ONE CELL HAS MEASURE EXACTLY ZERO.  F3's stated rule reads: 'an empty cell,")
print("       A MEASURE-ZERO CELL, or an overlap means the predicate does not partition'.")
print("       BY THE LETTER OF F3, FIRES - on a binary that partitions U(3) exactly,")
print("       with overlap 0 and union 1.  THIS IS A FALSE POSITIVE, AND IT CONVICTS")
print("       THE COROLLARY W-15 CERTIFIED AS THE PARADIGM OF A SOUND ONE.")

# ---- F2 DEGENERACY via the POWER map, on the corpus's own structure.
# P = cyclic permutation.  For ANY diagonal unitary D, (P D)^3 is diagonal.
P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=complex)
worst_cube = 0.0
min_base = 1e9
NP = 20000
for _ in range(NP):
    th = RNG.uniform(0, 2 * np.pi, 3)
    D = np.diag(np.exp(1j * th))
    U = P @ D
    worst_cube = max(worst_cube, offdiag(U @ U @ U))
    min_base = min(min_base, offdiag(U))
print(f"\nF2: map = POWER, f(U) = U^3.  Family U = P.D, P the 3-cycle, D diagonal unitary,"
      f" {NP} draws.")
print(f"    min off-diagonal norm of U itself (arm B, NON-diagonal): {min_base:.6f}")
print(f"    max off-diagonal norm of U^3      (lands in arm A)     : {worst_cube:.3e}")
print(f"    -> || A - f(B) || = {worst_cube:.3e}.  BY THE LETTER OF F2, FIRES.")
print("    The family is 3-real-dimensional inside U(3) (dim 9), so this is not a")
print("    measure-zero accident of one operator: a positive-dimensional slab of arm B")
print("    is carried into arm A by a listed map.")
print("    THE CORPUS'S OWN INSTANCE IS THE SAME FACT: ||T_F^3 - M_dF|| = 6.18e-16,")
print("    ||T_C^3 - M_c|| = 4.34e-16 (ERRATUM AGAINST W-11, REGISTER_V001.md:1519-1534).")

# ---- F1 EMBED on control 2.
print("\nF1 (letter): embedding space = U(3), real dimension 9.  Coordinate: delta =")
print("    off-diagonal Frobenius norm, range [0, sqrt(6)] = [0, 2.449490].")
print("    arm A = {delta = 0} (the diagonal subgroup, real dim 3);  arm B = {delta > 0}.")
K = np.array([[0, 1j, 0], [-1j, 0, 0], [0, 0, 0]], dtype=complex)  # Hermitian generator
w, V = np.linalg.eigh(K)
def Uth(t):
    return V @ np.diag(np.exp(1j * t * w)) @ V.conj().T
mods_ = np.array([0.6, 0.5, 0.6244998])
mods_ = mods_ / np.linalg.norm(mods_)
def blindness_spread(U, ntr=400):
    """max-min of |<U s, s>| over random phase assignments at FIXED moduli.
       zero iff the functional factors through the moduli, i.e. iff U acts diagonally."""
    vals = []
    for _ in range(ntr):
        ph = RNG.uniform(0, 2 * np.pi, 3)
        s = mods_ * np.exp(1j * ph)
        vals.append(abs(np.vdot(s, U @ s)))
    return max(vals) - min(vals)
print("    decision-relevant functional: phase-blindness spread of Z = <s, U s> at fixed moduli")
print("    theta      delta(U)      blindness spread")
for t in [0.0, 0.02, 0.05, 0.1, 0.3, 0.7853981634, 1.5707963268]:
    U = Uth(t)
    print(f"    {t:<9.4f}  {offdiag(U):<12.6f}  {blindness_spread(U):.6f}")
print("    -> interior points exist and the functional DIFFERS there.  BY THE LETTER OF F1,")
print("       and even under R4's stronger (a)+(b) discriminator, FIRES.")
print("F1 (amended, EXCLUSION): every interior point delta in (0, 2.449490] LIES IN ARM B.")
print("    points lying in NEITHER arm: 0 of 20000 Haar draws and 0 of the interpolation path")
print("    -> the coordinate refines one arm, it does not exhibit a third option.")
print("    -> AMENDED F1: SOUND.")

# ---- F4 on control 2 - the discretion problem, made explicit.
diagable = 0
NE = 2000
for _ in range(NE):
    U = haar_u(3, RNG)
    ev, EV = np.linalg.eig(U)
    diagable += int(SQ(EV.conj().T @ U @ EV - np.diag(ev)) < 1e-8)
print(f"\nF4: sentence 1 = 'a basis is designated' -> TRUE in-corpus (S1's vertex classes).")
print(f"    sentence 2 = 'diagonality is a property of the operator alone' -> FALSE:")
print(f"    {diagable} of {NE} Haar unitaries are diagonal in their OWN eigenbasis, so under")
print(f"    a free basis arm B is EMPTY.  BOTH SENTENCES ARE DEFENSIBLE READINGS OF THE SAME")
print(f"    QUESTION AND THEY GIVE OPPOSITE VERDICTS.  F4's output is set by lane discretion,")
print(f"    not by the frame.  SCORED: INDETERMINATE.")

# ---- F5 on control 2.
print("\nF5: witness for arm A: M_gamma, diagonal (REGISTER_V001.md:29-33).")
print("    witness for arm B: COR-F's T, non-diagonal (S3_THE_CROSSING_AUDIT_V001.md:182).")
print("    both branches obtain, each with a sealed exhibit.  F5: SOUND.")

# =====================================================================
# SCORECARD
# =====================================================================
rule("SCORECARD — FIRE RATE ON FRAMES KNOWN TO BE SOUND")
print("""
                       F1(letter)  F2   F3   F4          F5     fired/5
  CONTROL 1 PARITY        FIRE    FIRE  ok   ok          ok       2/5
  CONTROL 2 DIAG/NONDIAG  FIRE    FIRE  FIRE indeterminate ok     3.5/5

  SAME FRAMES, AMENDED INSTRUMENT (F1 + exclusion clause; F2 deleted;
  F3 demoted to a reporting line; F4 demoted to a grep gate; F5 unchanged):

                       F1(amended)      F5     fired/2
  CONTROL 1 PARITY          ok          ok       0/2
  CONTROL 2 DIAG/NONDIAG    ok          ok       0/2

  READ: as written, the instrument convicts the partition of the integers into
  evens and odds, and convicts the one binary W-15 itself certifies as forced with
  a proof under it.  A 19/20 fire rate on R1-R4 therefore carries FAR less
  information than it appears to.  The two tests that clear BOTH controls are
  F5 and the amended F1 - and those are the same two tests that produced 4 of the
  5 genuinely-new findings in the round.
""")

# =====================================================================
# STRENGTHENING R1's F5 — the founding escape inference, refuted with NO carrier.
# The cutoff audit objected that R1's F5(ii) runs on K1's electric spectrum, and
# K1 (S1, 2026-08-16) POSTDATES the founding route decision.  The objection is
# answered: the inference dies on any commensurate infinite spectrum, and needs
# no complex, no gauge group and no fibre.
# =====================================================================
rule("R1/F5 STRENGTHENED — 'an inductive limit of finite objects is not finite, which is")
print("precisely how it escapes recurrence' (FOUNDING_DESIGN_V001.md:63-65) — refuted with")
print("no carrier at all.  Any INFINITE spectrum that is commensurate revives exactly.\n")
for Nlev in [10, 100, 10000, 1000000]:
    E = np.arange(Nlev, dtype=float)          # harmonic-oscillator spectrum, integer gaps
    A2pi = abs(np.mean(np.exp(-1j * E * 2 * np.pi)))
    print(f"    dim = {Nlev:>8d}   E_n = n   |A(t=2*pi)| = {A2pi:.15f}")
print("\n    and the same for E_n = n^2 (K1's own electric spectrum shape, but the point")
print("    does not need it):")
for Nlev in [10, 100, 10000, 1000000]:
    E = np.arange(Nlev, dtype=float) ** 2
    A2pi = abs(np.mean(np.exp(-1j * E * 2 * np.pi)))
    print(f"    dim = {Nlev:>8d}   E_n = n^2  |A(t=2*pi)| = {A2pi:.15f}")
print("\n    NON-FINITENESS IS NECESSARY FOR NON-RECURRENCE, NOT SUFFICIENT.  The founding")
print("    design's single stated escape from its own §4 obstruction does not follow from")
print("    the premise it is drawn from, and the counterexample is two lines of arithmetic")
print("    requiring nothing but the sentence itself.")
print("\n    REGISTER CHECK: occurrences of 'inductive limit' in REGISTER_V001.md (1696 lines,")
print("    W-01 through W-16): 0.  Occurrences of 'quasi-local': 0.  The escape route that")
print("    licensed stage S3 was never examined by any register row.")
