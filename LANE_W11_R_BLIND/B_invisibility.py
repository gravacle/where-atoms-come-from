"""LEG B -- three ready states with IDENTICAL pi.  Is the incidence invisible?
   ARMS ARE DIFFED EXPLICITLY BEFORE ANY NUMBER IS REPORTED (W-08 isolation-audit rule:
   the commonest fatal defect is ZERO variables moved)."""
import numpy as np, sys
sys.path.insert(0,'.')
from wcore import *
np.set_printoptions(precision=12, linewidth=200)

f, c = 1.0, np.sqrt(2.0)        # the ONLY generic connection the corpus publishes (S4:603, W-10 N-4)
k = K1(f, c)
print("LEG B -- INVISIBILITY OF INCIDENCE UNDER THE TWO CLOCKS.  K1, f=1.0, c=sqrt(2) (generic).")
print("   edge-phase split gamma_F:", k.phF, " sum =", k.phF.sum())
print("   edge-phase split gamma_C:", k.phC, " sum =", k.phC.sum())

# ---------------- the three ready states -------------------------------------------
pA = np.array([0.40,0.15,0.15,0.15,0.15])                    # S2/S3's published ready state
sA = np.sqrt(pA).astype(complex)
pB = np.array([0.40,0.25,0.05,0.02,0.28])                    # SAME pi, moved WITHIN classes
sB = np.sqrt(pB).astype(complex)
ph = np.array([0.0, 0.7, -1.9, 2.3, 0.4])                    # SAME |s|, phases only
sC = np.sqrt(pA)*np.exp(1j*ph)

states = {'sA (published)':sA, 'sB (within-class moved)':sB, 'sC (phases only)':sC}

print("\nB.0  *** ARMS DIFF -- proving the three arms are not the same bytes ***")
for n,s in states.items():
    print("   %-24s |s|^2 = %s" % (n, np.array2string(np.abs(s)**2, precision=6)))
    print("   %-24s arg s = %s" % ("", np.array2string(np.angle(s), precision=6)))
print("   ||sA - sB|| = %.6f   ||sA - sC|| = %.6f   ||sB - sC|| = %.6f" %
      (np.linalg.norm(sA-sB), np.linalg.norm(sA-sC), np.linalg.norm(sB-sC)))
print("   pi(sA) = %s" % np.array2string(k.pi(sA), precision=15))
print("   pi(sB) = %s" % np.array2string(k.pi(sB), precision=15))
print("   pi(sC) = %s" % np.array2string(k.pi(sC), precision=15))
print("   max |pi(sA)-pi(sB)| = %.2e   max |pi(sA)-pi(sC)| = %.2e   -> pi IDENTICAL, states NOT"
      % (np.abs(k.pi(sA)-k.pi(sB)).max(), np.abs(k.pi(sA)-k.pi(sC)).max()))

# ---------------- circuit clock -----------------------------------------------------
K = 2000
Zc = {n: Z_circuit(k, s, K) for n,s in states.items()}
Ac = np.array([np.abs(Zc[n]) for n in states])
print("\nB.1  CIRCUIT clock  (tick = one whole circuit in each branch; branch op = M = T^3)")
print("   |Z_k| spread across the three states, max over k<=%d : %.3e" % (K, np.ptp(Ac,axis=0).max()))
print("   k=1..6 for each state:")
for n in states: print("      %-24s %s" % (n, np.array2string(np.abs(Zc[n])[:6], precision=12)))

# ---------------- edge clock --------------------------------------------------------
N = 2000
Ze = {n: Z_edge(k, s, N) for n,s in states.items()}
Ae = np.array([np.abs(Ze[n]) for n in states])
sp = np.ptp(Ae, axis=0)
print("\nB.2  EDGE clock  (tick = one edge in each branch; branch op = T)")
print("   |Z^T_n| spread across the three states, max over n<=%d : %.3e" % (N, sp.max()))
print("   spread at n = 1..12 : %s" % np.array2string(sp[:12], precision=3))
print("   spread on n in 3Z (n=3,6,...,%d): max = %.3e" % (N - N%3, sp[2::3].max()))
print("   spread off 3Z                    : min = %.3e" % np.min(np.delete(sp, np.s_[2::3])))
print("   n=1..6 for each state:")
for n in states: print("      %-24s %s" % (n, np.array2string(np.abs(Ze[n])[:6], precision=12)))

# ---------------- the isolation that matters ---------------------------------------
print("\nB.3  *** THE ISOLATION.  Hold the CLOCK fixed at whole circuits and swap the")
print("     transport operator M <-> T.  This moves the variable the question NAMES. ***")
worst = 0.0
for n,s in states.items():
    for kk in range(1, 51):
        zM = Z_circuit(k, s, kk)[-1]
        zT = Z_pair(k, s, 3*kk, 3*kk)      # same clock, edge operator iterated
        worst = max(worst, abs(zM - zT))
print("   max |<M_F^k s, M_C^k s> - <T_F^{3k} s, T_C^{3k} s}>| over 3 states, k<=50 : %.2e" % worst)
print("   ZERO VARIABLES MOVED.  At every time BOTH conventions name, they are the same number.")
print("   The reported difference lives entirely at times the circuit clock does not name.")

# ---------------- rates -------------------------------------------------------------
print("\nB.4  RATES")
mP = mahler4(*k.pi(sA))
print("   m(P), P = %s   (N1's polynomial)   = %.12f" % (np.array2string(k.pi(sA), precision=3), mP))
for NN in (2000, 20000, 200000):
    r_c = [rate(Z_circuit(k,s,NN)) for s in states.values()]
    r_e = [rate(Z_edge(k,s,NN))    for s in states.values()]
    print("   N=%7d  circuit rate/tick: %s  spread %.2e" % (NN, ["%.9f"%x for x in r_c], np.ptp(r_c)))
    print("             edge    rate/tick: %s  spread %.2e" % (["%.9f"%x for x in r_e], np.ptp(r_e)))
    print("             3 x edge rate    : %s   (compare m(P)=%.9f)" % (["%.9f"%(3*x) for x in r_e], mP))
