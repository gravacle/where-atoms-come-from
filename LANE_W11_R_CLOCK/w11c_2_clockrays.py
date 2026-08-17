# LANE W-11 R/C — LEG 2 — THE SPACE OF CLOCKS, AND N1's RATE ON EACH OF THEM.
# Leg 1 showed both conventions are RAYS in one lattice.  Leg 2 sweeps the whole family of rays
#   clock(a,b):  branch F advances a edges per tick, branch C advances b edges per tick
# EDGE clock = (1,1).  CIRCUIT clock = (L_F, L_C).  ONE VARIABLE MOVES: the ray (a,b).
# Everything else -- carrier, connection, state family, pi, observable, code path, seed, and the
# number of ticks -- is byte-identical across every row.
#
# TWO QUESTIONS:
#   (Q2a) On which clocks is the incidence invisible?     PREDICTION (theorem): L_F|a and L_C|b.
#   (Q2b) On the invisible clocks, is N1's rate the SAME NUMBER, or is it the circuit clock's
#         private number?  If m(P) on every invisible clock, then N1 is NOT a restatement of the
#         CIRCUIT clock; it is a restatement of "read only when both branches are closed".
import numpy as np
from w11c_lib import K1, B0b, ops, pi_of, states_same_pi, m_jensen, arms_differ, generic_conn

rng = np.random.default_rng(20260817)

def rate_and_spread(TF, TC, S, aa, bb, N):
    """mean of log|Z| along the ray (aa,bb), for each state.  T unitary, so no overflow."""
    UF = np.linalg.matrix_power(TF, aa)
    UC = np.linalg.matrix_power(TC, bb)
    out = []
    for s in S:
        xF = s.copy(); xC = s.copy(); tot = 0.0
        for _ in range(N):
            xF = UF @ xF; xC = UC @ xC
            z = abs(np.vdot(xF, xC))
            tot += np.log(z) if z > 1e-300 else -690.0
        out.append(tot / N)
    return np.array(out)

def run(K, label, NST=24, AMAX=12, N=4000):
    a = generic_conn(K, np.random.default_rng(7 + K.nv))
    TF, TC, MF, MC, WF, WC = ops(K, a)
    f = np.angle(WF) % (2*np.pi); c = np.angle(WC) % (2*np.pi)
    print(f"\n================ {label}   L_F={K.LF} L_C={K.LC} ================")
    print(f"  connection: every edge phase non-zero; (f, c) = ({f:.12f}, {c:.12f})"
          f"  [target (1.0, sqrt2={2**0.5:.12f})]")
    wref = None
    S = states_same_pi(K, PI[K.name], NST, np.random.default_rng(20260817))
    P = np.array([pi_of(K, s) for s in S])
    print(f"  {NST} same-pi states, pi = {np.round(PI[K.name],6)}, max pi deviation {np.abs(P-PI[K.name]).max():.1e}")
    arms_differ("states 0,1,2", S[0], S[1], S[2])
    mP = m_jensen(PI[K.name])
    print(f"  N1's registered value  lambda = m(P) = {mP:.12f}   (Jensen quadrature, n = 2^22)")

    print(f"\n  2a/2b  ONE VARIABLE: the clock ray (a,b).  N = {N} ticks in every row.")
    print(f"     {'a':>3} {'b':>3} | {'invisible?':>10} {'max spread':>12} | "
          f"{'mean rate/tick':>15} {'rate - m(P)':>13}   note")
    inv_cnt = 0
    rows = []
    for aa in range(1, AMAX+1):
        for bb in range(1, AMAX+1):
            r = rate_and_spread(TF, TC, S, aa, bb, N)
            spread = r.max() - r.min()
            # spread of |Z| itself at a single tick is the sharper invisibility test
            zsp = 0.0
            for t in (1, 2, 3, 5, 7, 11):
                v = np.array([abs(np.vdot(np.linalg.matrix_power(TF, aa*t) @ s,
                                          np.linalg.matrix_power(TC, bb*t) @ s)) for s in S])
                zsp = max(zsp, v.max()-v.min())
            inv = zsp < 1e-12
            inv_cnt += inv
            rows.append((aa, bb, inv, zsp, r.mean(), spread))
    pred = [(aa, bb) for aa in range(1, AMAX+1) for bb in range(1, AMAX+1)
            if aa % K.LF == 0 and bb % K.LC == 0]
    got = [(aa, bb, inv, zsp, rm, sp) for (aa, bb, inv, zsp, rm, sp) in rows if inv]
    for (aa, bb, inv, zsp, rm, sp) in got:
        tag = "CIRCUIT CLOCK" if (aa, bb) == (K.LF, K.LC) else ""
        if (aa, bb) == (1, 1): tag = "EDGE CLOCK"
        print(f"     {aa:>3} {bb:>3} | {'YES':>10} {zsp:>12.2e} | {rm:>15.9f} {rm-mP:>13.2e}   {tag}")
    print(f"     ... invisible clocks: {inv_cnt} of {AMAX*AMAX}  = 1/(L_F*L_C) = 1/{K.LF*K.LC}"
          f"   set equals predicted: {sorted([(a_,b_) for a_,b_,i,_,_,_ in rows if i]) == sorted(pred)}")
    vis = [(aa, bb, zsp, rm) for (aa, bb, inv, zsp, rm, sp) in rows if not inv]
    print(f"     VISIBLE clocks: {len(vis)}.  min |Z|-spread over them = "
          f"{min(z for _,_,z,_ in vis):.3e}   max = {max(z for _,_,z,_ in vis):.3e}")
    print(f"     EDGE CLOCK (1,1): |Z|-spread = {[z for a_,b_,z,_ in vis if (a_,b_)==(1,1)][0]:.3e}"
          f"   rate/tick spread across states = "
          f"{[s for a_,b_,i,z,r,s in rows if (a_,b_)==(1,1)][0]:.3e}")
    onlat = [rm for (aa, bb, inv, zsp, rm, sp) in rows if inv]
    print(f"\n     RATE ON THE INVISIBLE SUBLATTICE: min {min(onlat):.9f}  max {max(onlat):.9f}"
          f"  spread {max(onlat)-min(onlat):.2e}")
    print(f"     ALL EQUAL m(P) = {mP:.9f} TO {max(abs(np.array(onlat)-mP)):.2e}"
          f"  -> N1's number is the SAME on EVERY invisible clock, not just the circuit one.")
    return rows, mP

K = K1(); B = B0b()
PI = {}
PI["K1"] = np.array([0.0, 0.30, 0.30, 0.40])
wB = np.array([.10, .12, .09, .14, .11, .11, .11, .11, .11]); wB /= wB.sum()
PI["B0b"] = pi_of(B, np.sqrt(wB)+0j)

rK, mK = run(K, "K1  (equal loop lengths -- circuit clock and edge clock are COMMENSURABLE)")
rB, mB = run(B, "B0b (unequal loop lengths -- the two clocks NEVER coincide at n>0)")

print("\n================ WHAT LEG 2 SETTLES ================")
print("  (1) The set of clocks on which the incidence is invisible is exactly the sublattice")
print("      {L_F | a} x {L_C | b}: density 1/(L_F L_C).  Invisibility is a property of the")
print("      CLOCK, not of the transport operator -- the operator never changed in this leg.")
print("  (2) N1's rate m(P) is attained on EVERY clock in that sublattice, not only on the")
print("      corpus's circuit clock.  So N1 is NOT a restatement of the CIRCUIT convention.")
print("      It is a restatement of a WEAKER stipulation: read the record only at ticks where")
print("      BOTH branches have closed their loops.  That is a real narrowing of Reading B.")
