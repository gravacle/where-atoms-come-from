# LANE W-11 R/C — LEG 1 — THE TWO-CLOCK LATTICE.
# CLAIM UNDER TEST (the registrar's): "invisibility holds exactly where both branch operators are
# DIAGONAL -- by construction under the circuit convention, and only on the sublattice
# n = 0 mod lcm under the edge one.  So carrier-independence is a restatement of the convention."
#
# LEG 1 ASKS A PRIOR QUESTION: is there one transport in play, or two?
# M_gamma = T^L is verified in the registrar's own leg A.  If so, the CIRCUIT convention is not a
# rival OPERATOR; it is a SAMPLING of COR-F's own operator.  Then B1-vs-B2 does not move "the
# transport convention" at all -- it moves the CLOCK.  Settle it on the full lattice
#        Z(mF, mC) = <T_F^{mF} s, T_C^{mC} s>
# of which BOTH conventions are rays:  EDGE (n,n);  CIRCUIT (L_F k, L_C k).
#
# ISOLATION LEDGER (leg 1)
#   HELD FIXED: carrier, connection a, the observable <branch_F,branch_C>, the state family,
#               the code path (one function Zlat for every cell of the lattice), the seed.
#               pi is held fixed BY CONSTRUCTION across the whole state family.
#   MOVED, ONE THING: the lattice point (mF, mC) -- i.e. the clock, and nothing else.
#   There is NO arm in this leg in which the operator changes: M never appears in the sweep.
import numpy as np
from w11c_lib import (K1, B0b, ops, pi_of, states_same_pi, Zlat, arms_differ, classes)

rng = np.random.default_rng(20260817)          # registrar's seed, so figures are comparable
np.set_printoptions(linewidth=200)

def report(K, a, pi, NST=64, MMAX=25):
    TF, TC, MF, MC, WF, WC = ops(K, a)
    print(f"\n================ {K.name}   |gamma_F| = {K.LF}   |gamma_C| = {K.LC} ================")
    print(f"  class multiset {K.multiset()}    pi = {np.round(pi,6)}")
    print(f"  || T_F^{K.LF} - M_dF || = {np.linalg.norm(np.linalg.matrix_power(TF,K.LF)-MF):.2e}"
          f"    || T_C^{K.LC} - M_c || = {np.linalg.norm(np.linalg.matrix_power(TC,K.LC)-MC):.2e}"
          f"   [THEOREM, not a control]")

    # --- 1a  DIAGONALITY of T^m is exact combinatorics: T is a cyclic shift with phases.
    print("\n  1a  DIAGONALITY OF THE BRANCH OPERATOR AS A FUNCTION OF THE TICK  (exact, structural)")
    dF = [m for m in range(0, 2*K.LF*K.LC+1) if np.allclose(np.linalg.matrix_power(TF,m),
                                np.diag(np.diag(np.linalg.matrix_power(TF,m))), atol=1e-12)]
    dC = [m for m in range(0, 2*K.LF*K.LC+1) if np.allclose(np.linalg.matrix_power(TC,m),
                                np.diag(np.diag(np.linalg.matrix_power(TC,m))), atol=1e-12)]
    print(f"      T_F^m diagonal at m = {dF}   <=>  {K.LF} | m   : {dF == [m for m in range(0,2*K.LF*K.LC+1) if m % K.LF == 0]}")
    print(f"      T_C^m diagonal at m = {dC}   <=>  {K.LC} | m   : {dC == [m for m in range(0,2*K.LF*K.LC+1) if m % K.LC == 0]}")

    # --- 1b  the state family.  ALL SAME pi.  Guard against a zero-variable control.
    S = states_same_pi(K, pi, NST, rng)
    P = np.array([pi_of(K, s) for s in S])
    print(f"\n  1b  {NST} ready states, IDENTICAL pi (max deviation {np.abs(P-pi).max():.2e}),")
    print(f"      differing only in within-class weight and phase.")
    arms_differ("state family, first three", S[0], S[1], S[2])

    # --- 1c  the lattice sweep.  ONE VARIABLE: (mF, mC).
    print(f"\n  1c  SPREAD OF |Z(mF,mC)| ACROSS THE {NST} SAME-pi STATES, OVER THE LATTICE")
    sp = np.zeros((MMAX+1, MMAX+1))
    for mF in range(MMAX+1):
        AF = [np.linalg.matrix_power(TF, mF) @ s for s in S]
        for mC in range(MMAX+1):
            AC = [np.linalg.matrix_power(TC, mC) @ s for s in S]
            v = np.array([abs(np.vdot(x, y)) for x, y in zip(AF, AC)])
            sp[mF, mC] = v.max() - v.min()
    zero = (sp < 1e-12)
    pred = np.zeros_like(zero)
    for mF in range(MMAX+1):
        for mC in range(MMAX+1):
            pred[mF, mC] = (mF % K.LF == 0) and (mC % K.LC == 0)
    print(f"      cells with spread < 1e-12 : {zero.sum()} of {zero.size}")
    print(f"      predicted set  (L_F | mF) and (L_C | mC) : {pred.sum()}")
    print(f"      SETS AGREE EXACTLY: {np.array_equal(zero, pred)}      "
          f"max spread OFF the sublattice = {sp[~pred].max():.3e}")
    print(f"      max spread ON the sublattice = {sp[pred].max():.3e}")
    # the three mixed quadrants, separately -- one branch closed, the other not
    q_bothopen = (~pred) & np.array([[ (mF % K.LF != 0) and (mC % K.LC != 0)
                                       for mC in range(MMAX+1)] for mF in range(MMAX+1)])
    q_Fopen    = np.array([[ (mF % K.LF != 0) and (mC % K.LC == 0)
                             for mC in range(MMAX+1)] for mF in range(MMAX+1)])
    q_Copen    = np.array([[ (mF % K.LF == 0) and (mC % K.LC != 0)
                             for mC in range(MMAX+1)] for mF in range(MMAX+1)])
    for nm, q in (("BOTH branches mid-loop", q_bothopen), ("only F mid-loop", q_Fopen),
                  ("only C mid-loop", q_Copen)):
        if q.sum():
            print(f"        {nm:<24} cells {q.sum():>4}   min spread {sp[q].min():.3e}   max {sp[q].max():.3e}")

    # --- 1d  the two conventions AS RAYS in that lattice
    print("\n  1d  THE TWO CONVENTIONS ARE RAYS IN THE SAME LATTICE")
    circ = [(K.LF*k, K.LC*k) for k in range(0, MMAX//max(K.LF,K.LC)+1)]
    edge = [(n, n) for n in range(0, MMAX+1)]
    print(f"      CIRCUIT ray (L_F k, L_C k) = {circ}")
    print(f"        every point in the invisibility sublattice? "
          f"{all(pred[p] for p in circ)}   -> BY CONSTRUCTION")
    inside = [p for p in edge if pred[p]]
    print(f"      EDGE ray (n,n): points inside the sublattice = {inside}   "
          f"(= multiples of lcm({K.LF},{K.LC}) = {np.lcm(K.LF,K.LC)})")
    both = [p for p in circ if p in edge and p != (0, 0)]
    print(f"      points on BOTH rays (n>0): {both}")
    print(f"      -> the operative variable between B1 and B2 is the RAY, not the operator:")
    print(f"         the sweep above never used M at all, and it reproduces both conventions.")
    return sp, pred, S, TF, TC, MF, MC

# ---------------------------------------------------------------------------- K1
K = K1()
aK = np.array([1.0, 0.37, 0.91, 2**0.5, 0.23, 1.77])     # registrar's leg-B connection
piK = np.array([0.0, 0.30, 0.30, 0.40])
spK, predK, SK, TFK, TCK, MFK, MCK = report(K, aK, piK)

# reproduction of the registrar's B3 identity, marked as a THEOREM not a control
print("\n  1e  REGISTRAR'S B3, RE-DERIVED:  Z_edge(3k) = Z_circuit(k) on K1 -- IDENTITY, NOT EVIDENCE")
w = max(abs(Zlat(s, TFK, TCK, 3*k, 3*k) - np.vdot(np.linalg.matrix_power(MFK,k)@s,
                                                  np.linalg.matrix_power(MCK,k)@s))
        for s in SK[:8] for k in range(1, 8))
print(f"      max |Z_edge(3k) - Z_circuit(k)| = {w:.2e}   [T]  follows from T^3 = M; could not")
print( "      have failed, and that is no charge against it because it is a theorem.")

# ---------------------------------------------------------------------------- B0b
B = B0b()
aB = np.random.default_rng(20260817).uniform(0, 2*np.pi, 18)
wB = np.array([.10, .12, .09, .14, .11, .11, .11, .11, .11]); wB /= wB.sum()
piB = pi_of(B, np.sqrt(wB)+0j)
spB, predB, SB, TFB, TCB, MFB, MCB = report(B, aB, piB)

print("\n================ THE ANSWER LEG 1 SETTLES ================")
print("  There is ONE transport semigroup in play, generated by COR-F's T.  M_gamma is T^L.")
print("  The CIRCUIT convention is a SCHEDULE on that generator, not a rival operator.")
print("  The invisibility set is a SUBLATTICE of the clock lattice: {L_F | mF} x {L_C | mC}.")
print("  The circuit ray lies inside it; the edge ray meets it only at multiples of lcm(L_F,L_C).")
print("  OPERATIVE VARIABLE, NAMED: the CLOCK RAY -- equivalently, whether the record is read")
print("  while a branch is MID-LOOP.  NOT 'the transport convention'.")
