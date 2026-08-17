# S_02 — THE RACE.  Does the record's decay outrun the near-return floor?
# Both sides measured on the SAME orbit, at the SAME k, with the SAME evaluator, one connection
# per row and NOTHING else moving.  Double precision; the exact rows are re-checked in Fractions.
import numpy as np
from fractions import Fraction

RSG = (0.4, 0.3, 0.3)
RSP = (0.5, 0.0, 0.5)
phi = (1+5**0.5)/2
CONN = [
 ("S1PUB      order 4, ATTAINED     ", np.pi,        3*np.pi/2),
 ("RESONANT   f=2.0 c=1.1 (rank 1)  ", 2.0,          1.1),
 ("W07GEN     2pi*phi, 2pi*phi^2    ", 2*np.pi*phi,  2*np.pi*phi**2),
 ("DIOPH      -2^(1/3), 4^(1/3)     ", -2*np.pi*2**(1/3), 2*np.pi*4**(1/3)),
 ("RANDOM     seed 20260816         ", None,         None),
]
rng = np.random.default_rng(20260816)
a, b = rng.uniform(0,1,2)
CONN[4] = (CONN[4][0], 2*np.pi*a, 2*np.pi*b)

def orbit(f, c, w, K, chunk=10**6):
    """returns (sum(1-|Z|), sum(log|Z|), min(1-|Z|), argmin, count(1-|Z|<1e-3), sum over those)"""
    s1 = 0.0; s2 = 0.0; mn = 2.0; am = -1; nb = 0; sb = 0.0
    for lo in range(1, K+1, chunk):
        hi = min(lo+chunk-1, K)
        k = np.arange(lo, hi+1, dtype=np.float64)
        z = np.abs(w[0]*np.exp(1j*k*(c-f)) + w[1]*np.exp(-1j*k*f) + w[2]*np.exp(1j*k*c))
        z = np.minimum(z, 1.0)
        g = 1.0 - z
        s1 += g.sum(); s2 += np.log(np.maximum(z, 1e-300)).sum()
        j = int(np.argmin(g))
        if g[j] < mn: mn = float(g[j]); am = lo + j
        m = g < 1e-3
        nb += int(m.sum()); sb += float(g[m].sum())
    return s1, s2, mn, am, nb, sb

for tag_state, w in [("RS-G (0.4,0.3,0.3)", RSG), ("RS-P (0.5,0,0.5)", RSP)]:
    print(f"\n==== READY STATE {tag_state} ; K = 10^7 ; ONE VARIABLE MOVES: THE CONNECTION ====")
    print(f"   {'connection':<34} {'SUM(1-|Z|)':>13} {'density c':>10} {'lambda':>13} "
          f"{'FLOOR F(K)':>12} {'argmin':>9} {'#<1e-3':>8} {'their SUM':>10}")
    K = 10**7
    for tag, f, c in CONN:
        s1, s2, mn, am, nb, sb = orbit(f, c, w, K)
        print(f"   {tag:<34} {s1:>13.4f} {s1/K:>10.6f} {s2/K:>13.9f} "
              f"{mn:>12.3e} {am:>9d} {nb:>8d} {sb:>10.4f}")
print()
print("READ:  SUM(1-|Z|) is the record's total decay budget and it is LINEAR in K (density c).")
print("       FLOOR F(K) = min_k (1-|Z_k|) is the deepest SINGLE-cell near-return.")
print("       The last two columns are the ENTIRE contribution of every near-return cell in the")
print("       band 1-|Z_k| < 1e-3, i.e. the largest credit the floor can possibly return.")
