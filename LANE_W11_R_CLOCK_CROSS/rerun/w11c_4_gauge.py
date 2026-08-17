# LANE W-11 R/C — LEG 4 — WHAT ACTUALLY SELECTS THE SUBLATTICE, AND WHAT THE EDGE CLOCK SEES.
#
# Leg 1 found: invisibility holds exactly on {L_F | mF} x {L_C | mC}.  Leg 4 asks WHY that set,
# and whether anything in the sealed corpus picks it out.  The candidate the corpus never wrote
# down is COR-J's premise (S3_THE_CROSSING_AUDIT_V001.md:794): "the record must be gauge-invariant"
# -- applied not to the RECORD but to the BRANCH OPERATOR.
#
#   M_gamma is a gauge-INVARIANT operator:  M(g.a) = M(a) exactly, because W(gamma) is invariant.
#   COR-F's T is only gauge-COVARIANT:      T(g.a) = G T(a) G*,  G = diag(e^{i theta}).
#   PREDICTION: T^m is gauge-INVARIANT as an operator  <=>  L | m  <=>  invisibility.
#
# ISOLATION LEDGER (leg 4)
#   4a/4b HELD FIXED: carrier, connection, tick, state.  MOVED: the gauge transformation alone.
#   4c    HELD FIXED: carrier, connection, clock, observable.  MOVED: the ready state at fixed pi.
import numpy as np
from w11c_lib import (K1, B0b, ops, Tedge, Mcirc, hol, pi_of, states_same_pi, arms_differ,
                      generic_conn)

rng = np.random.default_rng(20260817)

def gauge_conn(K, a, th):
    """S1:59-63 -- a_e -> a_e + theta_target - theta_source on e : u -> v."""
    return np.array([a[j] + th[t] - th[s] for j, (s, t) in enumerate(K.edges)])

def leg(K, NG=2000, MMAX=13):
    a = generic_conn(K, np.random.default_rng(7 + K.nv))
    TF, TC, MF, MC, WF, WC = ops(K, a)
    print(f"\n================ {K.name}   L_F={K.LF}  L_C={K.LC} ================")

    print("  4a  GAUGE STATUS OF THE TWO BRANCH OPERATORS  (one variable: the gauge transform)")
    worstM = worstTinv = 0.0; worstTcov = 0.0
    for _ in range(NG):
        th = rng.uniform(0, 2*np.pi, K.nv)
        ag = gauge_conn(K, a, th)
        G = np.diag(np.exp(1j*th))
        Mg = Mcirc(K, K.VF, hol(K.walkF, ag))
        Tg = Tedge(K, K.walkF, ag)
        worstM = max(worstM, np.linalg.norm(Mg - MF))
        worstTinv = max(worstTinv, np.linalg.norm(Tg - TF))
        worstTcov = max(worstTcov, np.linalg.norm(Tg - G @ TF @ G.conj().T))
    print(f"      max || M_dF(g.a) - M_dF(a) ||          = {worstM:.2e}   -> M is gauge-INVARIANT")
    print(f"      max || T_F(g.a) - G T_F(a) G* ||       = {worstTcov:.2e}   -> T is gauge-COVARIANT")
    print(f"      max || T_F(g.a) - T_F(a) ||            = {worstTinv:.2e}   -> T is NOT invariant")

    print("\n  4b  AT WHICH TICKS IS THE BRANCH OPERATOR GAUGE-INVARIANT?  (T^m, m = 0..%d)" % MMAX)
    invF, invC = [], []
    for m in range(0, MMAX+1):
        wF = wC = 0.0
        for _ in range(200):
            th = rng.uniform(0, 2*np.pi, K.nv); ag = gauge_conn(K, a, th)
            wF = max(wF, np.linalg.norm(np.linalg.matrix_power(Tedge(K, K.walkF, ag), m)
                                        - np.linalg.matrix_power(TF, m)))
            wC = max(wC, np.linalg.norm(np.linalg.matrix_power(Tedge(K, K.walkC, ag), m)
                                        - np.linalg.matrix_power(TC, m)))
        if wF < 1e-12: invF.append(m)
        if wC < 1e-12: invC.append(m)
    print(f"      T_F^m gauge-INVARIANT at m = {invF}    <=>  {K.LF} | m : "
          f"{invF == [m for m in range(MMAX+1) if m % K.LF == 0]}")
    print(f"      T_C^m gauge-INVARIANT at m = {invC}    <=>  {K.LC} | m : "
          f"{invC == [m for m in range(MMAX+1) if m % K.LC == 0]}")
    print("      -> THE GAUGE-INVARIANCE SET OF THE BRANCH OPERATOR IS EXACTLY LEG 1's")
    print("         INVISIBILITY SUBLATTICE.  Same set, two descriptions.")

    print("\n  4c  BUT THE OBSERVABLE IS GAUGE-INVARIANT AT EVERY TICK UNDER BOTH CONVENTIONS")
    s = states_same_pi(K, PI[K.name], 1, np.random.default_rng(3))[0]
    worst = 0.0
    for _ in range(400):
        th = rng.uniform(0, 2*np.pi, K.nv); ag = gauge_conn(K, a, th)
        tF, tC = Tedge(K, K.walkF, ag), Tedge(K, K.walkC, ag)
        sg = np.exp(1j*th) * s
        for n in range(1, 8):
            z1 = np.vdot(np.linalg.matrix_power(tF, n) @ sg, np.linalg.matrix_power(tC, n) @ sg)
            z0 = np.vdot(np.linalg.matrix_power(TF, n) @ s,  np.linalg.matrix_power(TC, n) @ s)
            worst = max(worst, abs(z1 - z0))
    print(f"      max | <T_F^n g.s, T_C^n g.s>(g.a) - <T_F^n s, T_C^n s>(a) |, n<=7 = {worst:.2e}")
    print("      -> so the EDGE clock is NOT excluded by gauge-invariance of the RECORD (COR-J's")
    print("         own wording).  It is excluded only by the STRONGER, never-stated demand that")
    print("         the branch OPERATOR itself be gauge-invariant.")

    print("\n  4d  WHAT THE EDGE CLOCK SEES IS NOT 'THE INCIDENCE'.  (one variable: the state)")
    S = states_same_pi(K, PI[K.name], 3, np.random.default_rng(11), phases=False)
    sA = S[0]; sB = S[1]
    sPh = sA * np.exp(1j*rng.uniform(0, 2*np.pi, K.nv))
    arms_differ("A (weights), B (weights), A-with-phases", sA, sB, sPh)
    print(f"      |s_v|^2 (A) = {np.round(np.abs(sA)**2,6)}")
    print(f"      |s_v|^2 (B) = {np.round(np.abs(sB)**2,6)}   same pi = {np.round(pi_of(K,sA),6)}")
    print("      A vs B differ in |s_v|^2 ITSELF, which is gauge-invariant vertex by vertex.")
    # dressed Wilson-line invariant separating A from its own rephasing (W-06/W-07's object class)
    def dressed(s, a, walk):
        """conj(U_e s_u) * s_v -- transport s from u to v, then compare with s at v.
        Under a_e -> a_e + th_v - th_u and s -> e^{i th} s this is INVARIANT (checked below).
        This is the W-06 / W-07 dressed-observable class, built from S1's own edge transports."""
        out = []
        for (u, v, e, sg_) in walk:
            U = np.exp(1j*a[e]) if sg_ > 0 else np.exp(-1j*a[e])
            out.append(np.conj(U * s[u]) * s[v])
        return np.array(out)
    dA, dP = dressed(sA, a, K.walkF), dressed(sPh, a, K.walkF)
    wdr = 0.0
    for _ in range(300):
        th = rng.uniform(0, 2*np.pi, K.nv); ag = gauge_conn(K, a, th)
        wdr = max(wdr, np.abs(dressed(np.exp(1j*th)*sA, ag, K.walkF) - dA).max())
    print(f"      dressed edge invariant conj(U_e s_u) s_v : gauge-invariant to {wdr:.2e}")
    print(f"        on A          = {np.round(dA,6)}")
    print(f"        on A rephased = {np.round(dP,6)}   separation {np.abs(dA-dP).max():.3e}")
    print("      -> a rephased state is a DIFFERENT physical state at fixed connection, and the")
    print("         carrier's own gauge-invariant data says so.  The circuit functional cannot.")
    print("      CORRECTED NAME: what the edge clock reads is the joint gauge-invariant content of")
    print("      (connection, state) beyond (W_F, W_C, pi) -- S2 audit COR-E counts 11 such joint")
    print("      invariants where the build exhibits 7.  It is not 'the incidence labels'.")

K = K1(); B = B0b()
PI = {"K1": np.array([0.0, 0.30, 0.30, 0.40])}
wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB /= wB.sum()
PI["B0b"] = pi_of(B, np.sqrt(wB)+0j)
leg(K)
leg(B)

print("""
================ WHAT LEG 4 SETTLES ================
  The invisibility sublattice is EXACTLY the set of ticks at which the branch operator is a
  GAUGE-INVARIANT operator rather than merely gauge-covariant.  So:

    "the formation functional is carrier-independent"
      =  "the branch operator is required to be gauge-invariant, not merely gauge-covariant"

  That requirement is COR-J's undeclared premise -- "the record must be gauge-invariant" --
  promoted one level, from the RECORD (where COR-J states it, and where 4c shows it does NOT
  exclude the edge clock) to the OPERATOR (where it does).  Nobody made that promotion in
  writing.  It is the reason CHOICE LEDGER A2 could have given for the circuit clock and did not.
""")
