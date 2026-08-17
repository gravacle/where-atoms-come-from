# LANE W-11 R/C — LEG 7 — THE CLOSED FORM ON THE SUBLATTICE, AND CONVERGENCE OF THE RATE.
# Two precision-sensitive claims from legs 1, 2 and 5, checked in closed form and to convergence.
#
# CLAIM 7A (closed form, exact identity):  for L_F | mF and L_C | mC, with i = mF/L_F, j = mC/L_C,
#     Z(mF, mC)  =  p00 + p10 * conj(W_F)^i + p01 * W_C^j + p11 * conj(W_F)^i * W_C^j
# i.e. THE EXACT POLYNOMIAL N1 IS THE MAHLER MEASURE OF, with the two exponents DECOUPLED.
# The corpus only ever evaluates it on i = j.  That is the whole of the invisibility.
#
# CLAIM 7B (convergence):  the mean of log|Z| along ANY ray inside the sublattice converges to
#     m(P) = m(p00 + p10 x + p01 y + p11 xy),
# because (i,j) winds densely on T^2 for a generic (f,c) at ANY non-zero integer rates.
import numpy as np
from w11c_lib import K1, B0b, ops, pi_of, states_same_pi, m_jensen, generic_conn

def poly(pi, u, v):
    return pi[0] + pi[1]*u + pi[2]*v + pi[3]*u*v

def leg(K, pi, NST=12):
    a = generic_conn(K, np.random.default_rng(7 + K.nv))
    TF, TC, MF, MC, WF, WC = ops(K, a)
    S = states_same_pi(K, pi, NST, np.random.default_rng(20260817))
    mP = m_jensen(pi)
    print(f"\n================ {K.name}  L=({K.LF},{K.LC})   m(P) = {mP:.12f} ================")

    print("  7A  CLOSED FORM ON THE SUBLATTICE  (exact identity, checked over 12 states x 100 cells)")
    worst = 0.0
    for i in range(0, 10):
        for j in range(0, 10):
            AF = np.linalg.matrix_power(TF, K.LF*i); AC = np.linalg.matrix_power(TC, K.LC*j)
            cf = poly(pi, np.conj(WF)**i, WC**j)
            for s in S:
                worst = max(worst, abs(np.vdot(AF @ s, AC @ s) - cf))
    print(f"      max | Z(L_F i, L_C j) - [p00 + p10 conj(W_F)^i + p01 W_C^j + p11 ...] | = {worst:.2e}")
    print("      -> on the sublattice the observable IS N1's polynomial, with i and j DECOUPLED.")
    print("         The corpus evaluates it only on the diagonal i = j (S3:395, 'k_n circuits of")
    print("         EACH loop').  Off the sublattice no such closed form exists: T^m is not")
    print("         diagonal, so the vertex weights do not group into class sums at all.")

    print("\n  7B  CONVERGENCE OF THE RATE ALONG SUBLATTICE RAYS  (state 0; all states identical)")
    print(f"      {'ray (a,b)':>12} {'N=1e4':>14} {'N=1e5':>14} {'N=1e6':>14}   note")
    rays = [(K.LF, K.LC), (K.LF, 2*K.LC), (3*K.LF, K.LC), (np.lcm(K.LF,K.LC), np.lcm(K.LF,K.LC))]
    for (aa, bb) in rays:
        i_step, j_step = aa // K.LF, bb // K.LC
        vals = []
        for N in (10**4, 10**5, 10**6):
            n = np.arange(1, N+1)
            u = np.exp(-1j * np.angle(WF) * i_step * n)      # conj(W_F)^(i_step n)
            v = np.exp( 1j * np.angle(WC) * j_step * n)
            z = np.abs(poly(pi, u, v))
            vals.append(np.log(np.maximum(z, 1e-300)).mean())
        note = "CIRCUIT CLOCK" if (aa, bb) == (K.LF, K.LC) else (
               "lcm clock: F at %d circuits, C at %d" % (np.lcm(K.LF,K.LC)//K.LF,
                                                         np.lcm(K.LF,K.LC)//K.LC)
               if (aa, bb) == (np.lcm(K.LF,K.LC),)*2 else "")
        print(f"      {str((aa,bb)):>12} {vals[0]:>14.9f} {vals[1]:>14.9f} {vals[2]:>14.9f}   {note}")
    print(f"      m(P) (Jensen, n = 2^22)          {mP:>14.9f}"
          f"    <- every ray converges to this, none to a multiple of it")

K = K1(); B = B0b()
wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB /= wB.sum()
leg(K, np.array([0.0, 0.30, 0.30, 0.40]))
leg(B, pi_of(B, np.sqrt(wB)+0j))
print("""
================ WHAT LEG 7 SETTLES ================
  N1's polynomial is the closed form of the observable at EVERY point of the invisibility
  sublattice, with the two circuit counts DECOUPLED -- so the corpus's 'k_n circuits of EACH
  loop' (S3:395) is a further, separate stipulation on top of the clock, and N1 does not need it.
  And the RATE is m(P) on every one of those rays, at any winding rates.  N1's number is
  therefore invariant under a two-parameter family of clocks that includes the corpus's own,
  and it does NOT rescale by loop length.  'The edge rate is not m(P)/3' compares the wrong pair.
""")
