# X2 — THE TEST BOTH THE REGISTRAR AND LANE C DECLINED: RIVAL TICKS WITH THE SAME L-th POWER.
#
# LANE C's SELF-FLAG SAYS, IN ONE SENTENCE, BOTH OF THESE:
#   (i)  "My lattice framing makes uniqueness LESS load-bearing (ANY RIVAL TICK WITH THE SAME
#         L-TH POWER GIVES THE SAME SUBLATTICE)"
#   (ii) "a rival tick that is not a cyclic shift could have a different diagonality set."
# (i) and (ii) contradict each other, and LANE C's HEADLINE NAME -- "the operative variable is
# whether the record is read while a branch is MID-LOOP" -- rests on (i).  X2 settles it.
#
# STRUCTURE THEOREM (one line, proved not asserted).  On the loop, M_gamma = W.I_L.  Every
# unitary U with U^L = M_gamma is U = W^{1/L} V with V^L = I.  So
#       {m : U^m is DIAGONAL}  =  {m : V^m is diagonal},
# which is {L | m} for V a full cyclic shift (COR-F's T) and ANYTHING ELSE for other V.
#
# TWO EXPLICIT RIVALS, BOTH UNITARY, BOTH WITH U^L = M_gamma EXACTLY:
#   D  = diag(W^{1/L}) on the loop, identity off it.        DIAGONALITY SET = ALL m.
#        FIBRE-WISE, so it lies in the local gauge group (W-06's corrected N4 mechanism) and is
#        GAUGE-INVARIANT AT EVERY POWER.  It is the corpus's OWN kind of operator -- scalar
#        multiplication of the loop's fibres -- ticking once per EDGE instead of once per circuit.
#   U2 = W^{-1/L} T^2  on a loop of EVEN length.             DIAGONALITY SET = EVEN m.
#        Not fibre-wise; built from COR-F's own T.  Exists on B0b's gamma_F (L=4), not on K1.
#
# ISOLATION LEDGER (X2).  HELD FIXED: carrier, connection a, ready-state family (same 64 states,
#   same seed), pi by construction, the observable |<branch_F s, branch_C s>|, the evaluator, and
#   THE CLOCK (the same tick index m).  MOVED, EXACTLY ONE THING: WHICH UNITARY IS THE TICK.
#   Arms are hashed as OUTPUT vectors, not inputs (W-10 N-6).
import hashlib
import numpy as np
from x_lib import *

def oh(x): return hashlib.sha256(np.ascontiguousarray(np.asarray(x)).tobytes()).hexdigest()[:12]

def run(K, a, pi, NST=64, MM=13):
    TF, TC = T_edge(K,K.wF,a), T_edge(K,K.wC,a)
    WF, WC = holo(K.wF,a), holo(K.wC,a)
    MF, MC = M_circ(K,K.VF,WF), M_circ(K,K.VC,WC)
    DF, DC = D_root(K,K.VF,WF,K.LF), D_root(K,K.VC,WC,K.LC)
    S = states_same_pi(K, pi, NST, np.random.default_rng(20260817))
    print(f"\n================ {K.name}   L=({K.LF},{K.LC})   pi={np.round(pi,6)} ================")
    print(f"  D unitary?  ||D*D-I|| = {np.linalg.norm(DF.conj().T@DF-np.eye(K.nv)):.2e}"
          f"    ||D_F^{K.LF} - M_F|| = {np.linalg.norm(np.linalg.matrix_power(DF,K.LF)-MF):.2e}"
          f"    ||D_C^{K.LC} - M_C|| = {np.linalg.norm(np.linalg.matrix_power(DC,K.LC)-MC):.2e}")
    print(f"  T unitary?  ||T*T-I|| = {np.linalg.norm(TF.conj().T@TF-np.eye(K.nv)):.2e}"
          f"    ||T_F^{K.LF} - M_F|| = {np.linalg.norm(np.linalg.matrix_power(TF,K.LF)-MF):.2e}")
    print(f"  D is FIBRE-WISE (diagonal): {np.allclose(DF, np.diag(np.diag(DF)))}"
          f"   -> lies in the local gauge group U(1)^V, exactly as M_gamma does (W-05 REGISTER:414-416)")

    print(f"\n  2a  THE SAME CLOCK, TWO TICKS WITH THE SAME L-th POWER.  ONE VARIABLE: THE TICK.")
    print(f"      {'m':>3} | {'spread under COR-F T':>21} | {'spread under D':>15} |  branch state at m")
    zT, zD = [], []
    for m in range(1, MM+1):
        sT = spread_over(S, np.linalg.matrix_power(TF,m), np.linalg.matrix_power(TC,m))
        sD = spread_over(S, np.linalg.matrix_power(DF,m), np.linalg.matrix_power(DC,m))
        zT.append(sT); zD.append(sD)
        tag = ("BOTH CLOSED" if m%K.LF==0 and m%K.LC==0 else
               "F mid-loop" if m%K.LF and m%K.LC==0 else
               "C mid-loop" if m%K.LC and m%K.LF==0 else "BOTH mid-loop")
        print(f"      {m:>3} | {sT:>21.3e} | {sD:>15.3e} |  {tag}")
    print(f"      ARMS-DIFF  T-arm out#{oh(np.array(zT))}   D-arm out#{oh(np.array(zD))}   "
          f"{'DISTINCT' if oh(np.array(zT))!=oh(np.array(zD)) else '*** BYTE-IDENTICAL ***'}")
    print(f"      max spread under D over ALL m<= {MM} : {max(zD):.2e}    "
          f"max under T at m not divisible by L : "
          f"{max(z for m,z in zip(range(1,MM+1),zT) if m%K.LF or m%K.LC):.3e}")
    print("      => AT EVERY TICK WHERE A BRANCH IS MID-LOOP, INVISIBILITY HOLDS UNDER D AND FAILS")
    print("         UNDER T.  'READ WHILE MID-LOOP' IS THEREFORE NOT THE OPERATIVE VARIABLE.")

    # 2b — the non-fibre-wise rival, only on an even loop
    if K.LF % 2 == 0:
        print(f"\n  2b  A NON-FIBRE-WISE RIVAL ON gamma_F (L_F = {K.LF} is even):  U = W_F^(-1/{K.LF}) T_F^2")
        U = np.exp(-1j*np.angle(WF)/K.LF) * (TF @ TF)
        # restore identity off the loop (the scalar must not touch off-loop fibres)
        off = [v for v in range(K.nv) if v not in K.VF]
        U[np.ix_(off,off)] = np.eye(len(off), dtype=complex)
        for v in off: U[v,:] = 0; U[:,v] = 0; U[v,v] = 1
        print(f"      ||U*U - I|| = {np.linalg.norm(U.conj().T@U-np.eye(K.nv)):.2e}"
              f"    ||U^{K.LF} - M_F|| = {np.linalg.norm(np.linalg.matrix_power(U,K.LF)-MF):.2e}"
              f"    U diagonal? {np.allclose(U,np.diag(np.diag(U)))}")
        dia = [m for m in range(0,2*K.LF+1)
               if np.allclose(np.linalg.matrix_power(U,m), np.diag(np.diag(np.linalg.matrix_power(U,m))), atol=1e-12)]
        print(f"      diagonality set of U : {dia}   vs COR-F's T : "
              f"{[m for m in range(0,2*K.LF+1) if m%K.LF==0]}   SAME SET: {dia==[m for m in range(0,2*K.LF+1) if m%K.LF==0]}")
        print("      => a NON-DIAGONAL unitary with the SAME L-th power and a STRICTLY LARGER")
        print("         diagonality set.  Lane C's parenthetical (i) is FALSE.")

    # 2c — N1 itself, under D, on the EDGE clock
    print(f"\n  2c  N1 UNDER THE FIBRE-WISE EDGE TICK D.  Closed form and rate, EVERY tick.")
    worst = 0.0
    for m in range(0, 12):
        for n in range(0, 12):
            AF, AC = np.linalg.matrix_power(DF,m), np.linalg.matrix_power(DC,n)
            uu = np.conj(np.exp(1j*np.angle(WF)*m/K.LF)); vv = np.exp(1j*np.angle(WC)*n/K.LC)
            cf = pi[0] + pi[1]*uu + pi[2]*vv + pi[3]*uu*vv
            for s in S[:12]:
                worst = max(worst, abs(np.vdot(AF@s, AC@s) - cf))
    print(f"      max | <D_F^m s, D_C^n s> - [p00 + p10 conj(W_F)^(m/L_F) + p01 W_C^(n/L_C) + p11 ...] |"
          f" = {worst:.2e}")
    print("      -> N1's POLYNOMIAL, at EVERY (m,n), with FRACTIONAL winding.  No sublattice needed.")
    mP = m_jensen(pi)
    for N in (10**4, 10**5, 10**6):
        n = np.arange(1, N+1)
        u = np.exp(-1j*np.angle(WF)*n/K.LF); v = np.exp(1j*np.angle(WC)*n/K.LC)
        z = np.abs(pi[0] + pi[1]*u + pi[2]*v + pi[3]*u*v)
        r_closed = np.log(np.maximum(z,1e-300)).mean()
        print(f"      EDGE-clock rate under D, N={N:>7}: {r_closed:.9f}    m(P) = {mP:.9f}"
              f"    diff {r_closed-mP:+.2e}")
    # and the same by BRUTE FORCE from the operator, no closed form, to kill the substitution charge
    xF = S[0].copy(); xC = S[0].copy(); tot = 0.0
    for t in range(1, 200001):
        xF = DF@xF; xC = DC@xC
        tot += np.log(max(abs(np.vdot(xF,xC)), 1e-300))
    print(f"      SAME NUMBER FROM THE OPERATOR ITSELF (no closed form), N=200000: {tot/200000:.9f}")
    print("      => N1'S RATE IS m(P) ON THE *EDGE* CLOCK UNDER A FIBRE-WISE TICK.")

K = K1(); B = B0b()
wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB /= wB.sum()
run(K,  generic_conn(K, np.random.default_rng(7+K.nv)),  np.array([0.,.30,.30,.40]))
run(B,  generic_conn(B, np.random.default_rng(7+B.nv)),  pi_of(B, np.sqrt(wB)+0j))

print("""
================ WHAT X2 SETTLES ================
 1. LANE C's PARENTHETICAL "any rival tick with the same L-th power gives the same sublattice"
    IS FALSE, exhibited twice: D (all m) and U = W^{-1/L}T^2 on an even loop (even m).
 2. LANE C's NAME FOR THE OPERATIVE VARIABLE -- "whether the record is read while a branch is
    MID-LOOP" -- IS THEREFORE WRONG.  Mid-loop is invisible under D at every tick.
 3. THE CORRECT NAME IS THE ONE LANE C's OWN LEG 4 COMPUTES AND THEN DOES NOT ADOPT: whether the
    BRANCH OPERATOR APPLIED AT THE TICK READ IS FIBRE-WISE-AND-LOOP-CONSTANT -- i.e. scalar
    multiplication of the loop's fibres by one phase.  That is W-01/S3's transport convention with
    "the whole-circuit holonomy" weakened to "any loop-constant phase", and it is W-06's already
    registered correction to N4: FIBRE-WISE-NESS, REGISTER_V001.md (W-06, 'N4's MECHANISM RESTORED
    CORRECTED').  The name was in the register before either lane ran.
 4. CONSEQUENCE FOR THE QUESTION: the clock is not the stipulation.  A finer clock on the corpus's
    OWN kind of operator preserves N1 entirely.  The stipulation is on the OPERATOR.
""")
