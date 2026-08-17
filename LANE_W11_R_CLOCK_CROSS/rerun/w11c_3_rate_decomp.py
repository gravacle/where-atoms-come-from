# LANE W-11 R/C — LEG 3 — THE RATE, AND THE REGISTRAR'S D1 "RESCALED x3" ROW.
# D1 compares  (circuit rate per circuit)  against  (edge rate per tick) x 3  and reports
# "the EDGE rate is NOT m(P)/3 and is NOT the same number for the three states."
# The second half is a finding.  The first half is a UNIT ERROR, and leg 3 exhibits it exactly.
#
# lambda is a MEAN OF log|Z| over the ticks the clock reads, because |Omega_N| = prod_{ticks}|Z|.
# So "per circuit" and "per tick" are not related by a factor of L unless one has already granted
# that the product accrues ONE FACTOR PER EDGE TICK -- which is the disputed stipulation itself.
# The exact relation on K1 (L_F = L_C = 3) is a partition by RESIDUE CLASS:
#     mean_{all n} log|Z_e(n)| = (1/3)[ mean_{n=0(3)} + mean_{n=1(3)} + mean_{n=2(3)} ]
# and  mean_{n=0(3)} log|Z_e(n)|  IS  the circuit rate, exactly, because T^3 = M.
#
# ISOLATION LEDGER (leg 3)
#   HELD FIXED: carrier, connection, states, observable, code path, seed, tick count.
#   MOVED, ONE THING: which residue class of ticks the average is taken over.
import numpy as np
from w11c_lib import K1, B0b, ops, pi_of, states_same_pi, m_jensen, arms_differ, generic_conn

def residue_means(TF, TC, S, N, mod, aa=1, bb=1):
    """mean of log|Z| along the ray (aa,bb), split by t mod `mod`.  Returns array [state, res].
    Vectorised over states: one matmul per tick for the whole family, identical code path
    for every state and every row."""
    UF = np.linalg.matrix_power(TF, aa); UC = np.linalg.matrix_power(TC, bb)
    XF = np.column_stack(S).astype(complex); XC = XF.copy()
    tot = np.zeros((len(S), mod)); cnt = np.zeros(mod)
    for t in range(1, N+1):
        XF = UF @ XF; XC = UC @ XC
        z = np.abs(np.sum(np.conj(XF) * XC, axis=0))
        tot[:, t % mod] += np.where(z > 1e-300, np.log(np.maximum(z, 1e-300)), -690.0)
        cnt[t % mod] += 1
    return tot / cnt[None, :], cnt

# ------------------------------------------------------------------ K1, registrar's own setup
print("== 3A  REPRODUCING THE REGISTRAR'S D1 SETUP EXACTLY (K1, its connection, its 3 states) ==")
K = K1()
aK = np.array([1.0, 0.37, 0.91, 2**0.5, 0.23, 1.77])
TF, TC, MF, MC, WF, WC = ops(K, aK)
sA = np.sqrt(np.array([0.40, 0.15, 0.15, 0.15, 0.15])) + 0j
sB = np.sqrt(np.array([0.40, 0.30, 0.00, 0.05, 0.25])) + 0j
sC = sA * np.exp(1j * np.array([0.0, 1.3, -0.7, 2.2, 0.4]))
S3 = [sA, sB, sC]
arms_differ("registrar states A,B,C", sA, sB, sC)
PI = pi_of(K, sA)
assert np.allclose(pi_of(K, sB), PI) and np.allclose(pi_of(K, sC), PI)
mP = m_jensen(PI)
print(f"   pi = {np.round(PI,6)}   m(P) = {mP:.12f}")
print(f"   genericity of (f,c) = ({np.angle(WF)%(2*np.pi):.9f}, {np.angle(WC)%(2*np.pi):.9f}):"
      f" c = 2+sqrt(2) is algebraic irrational, f = 2.28 rational -> no relation m f + n c = 2 pi j.")

for N in (30000, 300000):
    R, cnt = residue_means(TF, TC, S3, N, 3)
    allmean = R.mean(axis=1)                       # equal counts up to +-1
    print(f"\n   N = {N} edge ticks")
    print(f"     {'':<22}{'state A':>15}{'state B':>15}{'state C':>15}{'spread':>11}")
    for r in range(3):
        lab = f"mean over n = {r} mod 3"
        print(f"     {lab:<22}{R[0,r]:>15.9f}{R[1,r]:>15.9f}{R[2,r]:>15.9f}"
              f"{R[:,r].max()-R[:,r].min():>11.2e}"
              + ("   <- CIRCUIT CLOCK, = m(P)" if r == 0 else "   <- PARTIAL CIRCUIT"))
    print(f"     {'mean over ALL n':<22}{allmean[0]:>15.9f}{allmean[1]:>15.9f}{allmean[2]:>15.9f}"
          f"{allmean.max()-allmean.min():>11.2e}   <- registrar's EDGE rate/tick")
    print(f"     {'the same x3':<22}{3*allmean[0]:>15.9f}{3*allmean[1]:>15.9f}{3*allmean[2]:>15.9f}"
          f"{'':>11}   <- registrar's 'per circuit'")
    print(f"     IDENTITY CHECK  3 x (mean over all n)  ==  sum of the three residue means:"
          f"  max |diff| = {np.abs(3*allmean - R.sum(axis=1)).max():.2e}")
    print(f"     residue 0 minus m(P): {np.abs(R[:,0]-mP).max():.2e}   (finite-N ergodic error)")
print("""
   READING.  The registrar's 'EDGE rescaled x3' number is
        (circuit rate)  +  (mean over n=1 mod 3)  +  (mean over n=2 mod 3),
   i.e. the circuit rate PLUS TWO EXTRA FACTORS PER CIRCUIT in the product |Omega|.  It is not
   m(P) rescaled and was never going to be.  The x3 is not a unit conversion: it PRESUPPOSES that
   the record accrues one factor per edge tick, which is the very stipulation under test.
   WHAT SURVIVES OF D1: the state-dependence of the residue-1 and residue-2 means, which is real.
   WHAT FALLS: 'the EDGE rate is NOT m(P)/3' as evidence -- it is a comparison of a mean over one
   sample set with three times a mean over a larger one.""")

# ------------------------------------------------------------------ does the RATE spread survive?
print("\n== 3B  IS THE *RATE* STATE-DEPENDENT UNDER THE EDGE CLOCK?  N1 IS A RATE STATEMENT. ==")
print("   ONE VARIABLE: the carrier (K1 vs B0b).  Same clock (1,1), same generic connection")
print("   construction, same state family size, same seed, same code path.")
B = B0b()
wB = np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); wB /= wB.sum()
PIB = pi_of(B, np.sqrt(wB)+0j)
for carrier, pi in ((K1(), np.array([0.0,0.30,0.30,0.40])), (B, PIB)):
    a = generic_conn(carrier, np.random.default_rng(7+carrier.nv))
    tF, tC, mF, mC, wf, wc = ops(carrier, a)
    S = states_same_pi(carrier, pi, 16, np.random.default_rng(20260817))
    mp = m_jensen(pi)
    print(f"\n   --- {carrier.name}  L=({carrier.LF},{carrier.LC})  m(P) = {mp:.12f}")
    print(f"       {'N':>9} {'EDGE rate spread':>18} {'CIRCUIT rate spread':>21} {'edge mean':>14}")
    for N in (5000, 50000, 500000):
        re_, _ = residue_means(tF, tC, S, N, 1, 1, 1)
        rc_, _ = residue_means(tF, tC, S, N // max(carrier.LF, carrier.LC), 1,
                               carrier.LF, carrier.LC)
        print(f"       {N:>9} {re_[:,0].max()-re_[:,0].min():>18.3e}"
              f" {rc_[:,0].max()-rc_[:,0].min():>21.3e} {re_[:,0].mean():>14.9f}")
    print(f"       circuit-clock rate - m(P) = {abs(rc_[:,0].mean()-mp):.2e}  [converging to m(P)]")
