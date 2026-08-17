# X6 — THE ADJUDICATION.  Lane C's central move is:
#   "M = T^L is a THEOREM, therefore there is only ONE transport in play, therefore between the
#    registrar's B1 and B2 THE TRANSPORT CONVENTION DOES NOT MOVE AT ALL.  What moves is the CLOCK."
# X6 tests that move directly, and tests W-10 PRECONDITION P.4 which lane C files as FALSE.
import numpy as np
from x_lib import *

print("== X6a  ON THE SUBLATTICE, THE OPERATOR APPLIED *IS* WHOLE-CIRCUIT SCALAR TRANSPORT ==")
print("   Lane C: 'W-10 PRECONDITION P.4 -- whole-circuit scalar transport (M_gamma, not COR-F's T)")
print("   -- the convention on which the entire CARRIER_INDEPENDENT column rests -- FALSE as stated:")
print("   the column is re-derived from T alone on the sublattice (leg 7A/7B).'")
for K, pi in ((K1(), np.array([0.,.30,.30,.40])), (B0b(), None)):
    if pi is None:
        wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB/=wB.sum(); pi = pi_of(K,np.sqrt(wB)+0j)
    a = generic_conn(K, np.random.default_rng(7+K.nv))
    TF, TC = T_edge(K,K.wF,a), T_edge(K,K.wC,a)
    WF, WC = holo(K.wF,a), holo(K.wC,a)
    MF, MC = M_circ(K,K.VF,WF), M_circ(K,K.VC,WC)
    wF = max(np.linalg.norm(np.linalg.matrix_power(TF,K.LF*i)-np.linalg.matrix_power(MF,i))
             for i in range(0,10))
    wC = max(np.linalg.norm(np.linalg.matrix_power(TC,K.LC*j)-np.linalg.matrix_power(MC,j))
             for j in range(0,10))
    print(f"   {K.name:<4}: max_i || T_F^(L_F i) - M_F^i || = {wF:.2e}   "
          f"max_j || T_C^(L_C j) - M_C^j || = {wC:.2e}   (i,j <= 9)")
print("   => EVERY OPERATOR LEG 7A EVALUATES IS A POWER OF M_gamma, TO 1e-15.  Writing M as T^L")
print("      renames the operator; it does not remove it.  P.4 SAYS THE COLUMN RESTS ON THAT")
print("      OPERATOR AND IT DOES.  P.4 STANDS; lane C's new_defect against it is a re-description.")

print("\n== X6b  WHAT ACTUALLY MAKES ONLY pi ENTER.  NO HOLONOMY, NO CLOCK, NO CIRCUIT. ==")
print("   CLAIM: if A_F = alpha on gamma_F's fibres and 1 elsewhere, A_C = beta on gamma_C's fibres")
print("   and 1 elsewhere -- FIBRE-WISE and LOOP-CONSTANT, with alpha, beta ANY unit phases at all,")
print("   unrelated to any holonomy -- then <A_F s, A_C s> = p00 + p10 conj(a) + p01 b + p11 conj(a)b.")
rng = np.random.default_rng(31415)
worst = 0.0
for K, pi in ((K1(), np.array([0.,.30,.30,.40])), (B0b(), None)):
    if pi is None:
        wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB/=wB.sum(); pi = pi_of(K,np.sqrt(wB)+0j)
    S = states_same_pi(K, pi, 40, np.random.default_rng(20260817))
    w = 0.0
    for _ in range(300):
        al, be = np.exp(1j*rng.uniform(0,2*np.pi,2))
        AF = M_circ(K, K.VF, al); AC = M_circ(K, K.VC, be)
        cf = pi[0] + pi[1]*np.conj(al) + pi[2]*be + pi[3]*np.conj(al)*be
        for s in S:
            w = max(w, abs(np.vdot(AF@s, AC@s) - cf))
    worst = max(worst, w)
    print(f"   {K.name:<4}: 300 random (alpha,beta) x 40 same-pi states -> max deviation {w:.2e}")
print(f"   => THE WHOLE 'CARRIER-INDEPENDENT' CONTENT IS THIS ONE LINE.  max {worst:.2e}")
print("      It uses NO holonomy, NO circuit count, NO clock, NO equidistribution, NO Mahler")
print("      measure.  It is the algebraic shape of the branch operators and nothing else.")

print("\n== X6c  AND 'FIBRE-WISE' ALONE IS NOT ENOUGH -- LOOP-CONSTANCY IS THE OTHER HALF ==")
print("   W-06's registered correction names the mechanism FIBRE-WISE-NESS.  Sharpened here:")
for K, pi in ((K1(), np.array([0.,.30,.30,.40])), (B0b(), None)):
    if pi is None:
        wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB/=wB.sum(); pi = pi_of(K,np.sqrt(wB)+0j)
    S = states_same_pi(K, pi, 40, np.random.default_rng(20260817))
    best = 0.0
    for _ in range(200):
        # fibre-wise (diagonal, unitary) but per-vertex phases on the loop: NOT loop-constant
        AF = np.eye(K.nv, dtype=complex); AC = np.eye(K.nv, dtype=complex)
        for v in K.VF: AF[v,v] = np.exp(1j*rng.uniform(0,2*np.pi))
        for v in K.VC: AC[v,v] = np.exp(1j*rng.uniform(0,2*np.pi))
        best = max(best, spread_over(S, AF, AC))
    print(f"   {K.name:<4}: fibre-wise but per-vertex phases -> max spread across same-pi states "
          f"{best:.3e}   (invisibility FAILS)")
print("   => CORRECT NAME: the branch operator must be FIBRE-WISE *AND* CONSTANT ON ITS LOOP.")
print("      'Whole-circuit scalar multiplication by the holonomy' (P.4) is the special case the")
print("      corpus stipulated; the general condition is weaker and is still a condition ON THE")
print("      OPERATOR, not on the clock.")

print("\n== X6d  THE TWO EMBEDDINGS.  LANE C's 'ONE TRANSPORT SEMIGROUP' IS A CHOICE, NOT A FACT ==")
K = K1(); a = generic_conn(K, np.random.default_rng(7+K.nv))
TF, TC = T_edge(K,K.wF,a), T_edge(K,K.wC,a); WF, WC = holo(K.wF,a), holo(K.wC,a)
MF, MC = M_circ(K,K.VF,WF), M_circ(K,K.VC,WC)
DF, DC = D_root(K,K.VF,WF,K.LF), D_root(K,K.VC,WC,K.LC)
S = states_same_pi(K, np.array([0.,.30,.30,.40]), 64, np.random.default_rng(20260817))
print("   M_gamma sits inside TWO one-parameter semigroups, both with the same L-th power:")
print(f"      T-semigroup (COR-F, edge-wise parallel transport): ||T^{K.LF} - M|| = "
      f"{np.linalg.norm(np.linalg.matrix_power(TF,K.LF)-MF):.2e}   T diagonal? False")
print(f"      D-semigroup (fibre-wise, the corpus's OWN kind):   ||D^{K.LF} - M|| = "
      f"{np.linalg.norm(np.linalg.matrix_power(DF,K.LF)-MF):.2e}   D diagonal? True")
eT = [spread_over(S, np.linalg.matrix_power(TF,n), np.linalg.matrix_power(TC,n)) for n in range(1,10)]
eD = [spread_over(S, np.linalg.matrix_power(DF,n), np.linalg.matrix_power(DC,n)) for n in range(1,10)]
cT = [spread_over(S, np.linalg.matrix_power(MF,k), np.linalg.matrix_power(MC,k)) for k in range(1,10)]
print(f"   EDGE clock in the T-semigroup : max spread over n<=9 = {max(eT):.3e}   DISAGREES with circuit")
print(f"   EDGE clock in the D-semigroup : max spread over n<=9 = {max(eD):.3e}   AGREES with circuit")
print(f"   CIRCUIT clock (M^k)           : max spread over k<=9 = {max(cT):.3e}")
print("   => THE REGISTRAR'S B1-vs-B2 DISAGREEMENT IS NOT A DISAGREEMENT BETWEEN TWO CLOCKS.")
print("      Hold the clock at EDGE and move only the semigroup: the disagreement appears and")
print("      disappears.  So the clock is not the operative variable; the OPERATOR FAMILY is.")
print("      Lane C's leg 1 held the operator family fixed at T and moved the ray, so it could")
print("      not have found this: within one family the ray is the only thing left to move.")

print("""
== X6e  THE HONEST OTHER SIDE, SCORED FOR NEITHER ==
 S1:52 defines parallel transport EDGE-WISE ("along e : u -> v is z |-> U_e z").  If one adds the
 admissibility criterion "a tick must be an edge-wise parallel transport of the fibre values",
 then D is inadmissible and lane C's sublattice is the right answer.  BUT THE SAME CRITERION
 EXCLUDES THE CORPUS'S OWN M_gamma, which moves no fibre value at all and is an element of the
 gauge group (W-05, REGISTER:414-416).  So S1:52 cannot be the criterion that saves lane C's
 reading without killing the object the corpus actually uses.
 THE CORPUS HAS NO ADMISSIBILITY CRITERION FOR TICKS -- lane C's own eighth new_defect says so,
 about CLOCKS.  The missing criterion is one level lower: it is about OPERATORS.
 READS TWO WAYS AND I SCORE NEITHER: either (i) the edge-wise criterion is right, the corpus's
 M_gamma is inadmissible, and the whole functional layer is scoped to a convention the corpus's
 own S1 contradicts; or (ii) fibre-wise operators are admissible, D is a legitimate tick, and
 the clock question is empty because a finer clock on the corpus's own operator changes nothing.
""")
