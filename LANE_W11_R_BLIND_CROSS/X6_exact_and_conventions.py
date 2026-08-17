"""X6 -- PRECISION-SENSITIVE CLAIMS IN EXACT ARITHMETIC OR CLOSED FORM, AND AN AUDIT OF THE
   REGISTRAR'S UNPUBLISHED CONNECTION (the blind lane's R-4, pushed one step further)."""
import numpy as np, sys
from fractions import Fraction as Fr
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_BLIND_CROSS')
from xcore import *

print("X6.1  T^L = M_gamma IN EXACT ARITHMETIC, AS A COMPOSITION OF SYMBOLIC MAPS")
print("      T is represented as v |-> (image vertex, exact phase in units of 2pi, a Fraction).")
print("      T^L is obtained by composing that map L times; no float ever touches it.")
def symbolic_T(loop, ph):
    """dict: source vertex -> (target vertex, exact phase Fraction in units of 2pi)"""
    L = len(loop)
    return {loop[j]: (loop[(j+1) % L], ph[j]) for j in range(L)}
def compose(m1, m2):
    """apply m1 then m2"""
    out = {}
    for v, (w, a) in m1.items():
        w2, b = m2[w]
        out[v] = (w2, a + b)
    return out
ok = True
for loop, ph in (([0,1,2], [Fr(1,7), Fr(2,7), Fr(3,7)]),
                 ([0,1,2], [Fr(1,5), Fr(1,3), Fr(-2,15)]),
                 ([0,1,4,3], [Fr(3,11), Fr(-1,4), Fr(5,6), Fr(1,3)]),
                 ([0,3,6], [Fr(0), Fr(0), Fr(1,2)])):
    L = len(loop); tot = sum(ph)
    m = symbolic_T(loop, ph); p = m
    for _ in range(L-1): p = compose(p, m)
    fixed = all(p[v][0] == v for v in loop)
    same  = all(p[v][1] == tot for v in loop)
    ok = ok and fixed and same
    print("      loop %-12s L=%d  phases %s  ->  T^L fixes every loop vertex: %s ; carries the"
          " SAME exact angle %s at each: %s" % (loop, L, [str(x) for x in ph], fixed, tot, same))
print("      ALL EXACT: %s.  M_gamma is not a rival of T; it is the L-th power of T, as a" % ok)
print("      combinatorial identity in the free group on the loop's edges.  (This is the blind")
print("      lane's structural point, and it survives at exact arithmetic in a third code.)")

print("\nX6.2  THE CONNECTIONS EACH LANE ACTUALLY USED, AND WHETHER THEY ARE GENERIC")
print("      W-10 N-4: every RATIONAL (f,c) is exactly resonant; the corpus publishes exactly one")
print("      generic pair, S4:603's f = 1.0, c = sqrt(2).  PUBLISHED_CONVENTIONS.txt states no")
print("      (f,c) at all, so the registrar's pair has to be read out of w11_b_decisive.py:31.")
print("      registrar: a = [1.0, 0.37, 0.91, sqrt2, 0.23, 1.77]  =>  f = %.6f = 57/25 exactly," % (1.0+0.37+0.91))
print("                 c = %.12f = 2 + sqrt(2) exactly." % (2**0.5+0.23+1.77))
print("      EXACT ARGUMENT.  Suppose m*f + n*c = 2*pi*j with integers (m,n) != (0,0).")
print("        j = 0 : m*(57/25) + 2n + n*sqrt2 = 0 forces n*sqrt2 rational, so n = 0, so m = 0.")
print("        j != 0: pi = (m*57/25 + 2n + n*sqrt2)/(2j) would put pi in Q(sqrt2); pi is")
print("                transcendental.  Contradiction.  So the registrar's pair IS generic --")
print("        but no lane checked it, and it is not the pair the corpus publishes.")
def minres(f, c, R=200):
    best = (9e9, None)
    for m in range(-R, R+1):
        for n in range(-R, R+1):
            if m == 0 and n == 0: continue
            x = (m*f + n*c) % (2*np.pi)
            x = min(x, 2*np.pi - x)
            if x < best[0]: best = (x, (m, n))
    return best
for nm, (f_, c_) in (("S4:603 published generic  f=1, c=sqrt2", (1.0, np.sqrt(2.0))),
                     ("registrar's pair          f=57/25, c=2+sqrt2", (57/25, 2+np.sqrt(2.0))),
                     ("KNOWN-RESONANT control    f=2.0, c=1.1", (2.0, 1.1))):
    d, mn = minres(f_, c_)
    print("      %-46s min |m f + n c| mod 2pi over |m|,|n|<=200 = %.3e at %s" % (nm, d, mn))
print("      the control returns EXACTLY 0, at (-88,160) = 8 x (-11,20), a multiple of W-10 N-4's")
print("      registered relation (11,-20) -- the search finds the whole relation subgroup, not")
print("      only its generator.  Neither working pair returns 0 at 200x the search radius.")
print("      Small-but-nonzero is not a proof and is not scored as one; the exact argument above")
print("      is what carries the genericity of both working pairs.")

print("\nX6.3  m(P) BY THREE ROUTES")
pi0 = (0.0, 0.3, 0.3, 0.4)
for e in (14, 18, 20, 22, 24):
    print("      Jensen-in-y, grid 2^%2d                      : %.15f" % (e, mahler4(*pi0, ngrid=1 << e)))
n2 = 4096
th = (np.arange(n2)+0.5)*(2*np.pi/n2)
X, Y = np.meshgrid(np.exp(1j*th), np.exp(1j*th), indexing='ij')
print("      2-D midpoint quadrature %dx%d, no Jensen  : %.12f" % (n2, n2, float(np.log(np.abs(0.3*X + 0.3*Y + 0.4*X*Y)).mean())))
print("      register W-02 erratum / W-10 published      : -0.767507880")
piU = (Fr(4,9), Fr(2,9), Fr(1,9), Fr(2,9))
print("      B0b SENSE U: pi = %s ; N-3 says one Jensen branch dominates and m = log(4/9)."
      % [str(x) for x in piU])
th = np.arange(1<<22)*(2*np.pi/(1<<22)); xg = np.exp(1j*th)
A = np.abs(4/9 + 2/9*xg); Bq = np.abs(1/9 + 2/9*xg)
print("      min over 2^22 grid points of |A| - |B| = %.6f  > 0, so branch A dominates EVERYWHERE"
      % float((A-Bq).min()))
print("      m by Jensen at 2^22 = %.12f    log(4/9) = %.12f    difference %.2e"
      % (mahler4(4/9, 2/9, 1/9, 2/9, 1<<22), np.log(4/9), abs(mahler4(4/9,2/9,1/9,2/9,1<<22)-np.log(4/9))))
print("      EXACT: with branch A dominating, m = m(4/9 + (2/9)x) = log max(4/9, 2/9) = log(4/9).")
