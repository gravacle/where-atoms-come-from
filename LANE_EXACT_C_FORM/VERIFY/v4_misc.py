"""V4 -- three smaller adversarial checks.

(a) THE SIGN FLIPS IN t8_weakfield.txt.  The lane's D column goes -2.9e-5 (lam=0.1131),
    +2.0e-3 (0.16), +5.4e-3 (0.2263), -3.97e-2 (0.32).  Two sign changes.  Recomputed at 60
    digits in the exactly-equivalent one-bath-qubit model.

(b) TEST 6 IS ARITHMETICALLY VACUOUS.  defect(m) = Qchi(m) - m*Qchi(1), and Qchi(m) = m*chi >= 0
    because Holevo chi is non-negative.  Hence |defect| <= m*Qchi(1) ALWAYS.  The "pairwise
    field" prediction C(m,2)*d2 exceeds that bound for every m >= 3, so it could never have been
    observed, whatever the physics.  Numbers below.

(c) THE F_2 / INTEGER LAYER, recomputed independently in exact arithmetic.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from fractions import Fraction
from mpmath import mp, mpf, cos, sin, sqrt, log, tanh
mp.dps = 50
from lane_utils import stab_nn2, stab_blocks, derived_logical_span, sp, f2_rank, pauli_vec, in_span

# ---------------------------------------------------------------- (a)
BETA = mpf(2); Z0 = -tanh(BETA)
TIMES = [mpf(1) + (mpf(12)*i)/24 for i in range(25)]
def bloch(K, lam, t):
    r = sqrt(1 + (lam*K)**2); a = lam*K/r; c = mpf(1)/r; th = 2*t*r
    C = cos(th); S = sin(th)
    return (a*c*Z0*(1-C), -a*Z0*S, Z0*C + c*c*Z0*(1-C))
def H2(x):
    if x <= 0 or x >= 1: return mpf(0)
    return -(x*log(x) + (1-x)*log(1-x))/log(2)
def Sv(v): return H2((1+sqrt(v[0]**2+v[1]**2+v[2]**2))/2)
def mix(pairs, lam, t):
    v = [mpf(0)]*3
    for wgt, K in pairs:
        b = bloch(K, lam, t)
        for i in range(3): v[i] += wgt*b[i]
    return v
def chi(plus, minus, lam):
    acc = mpf(0)
    for t in TIMES:
        vp = mix(plus, lam, t); vm = mix(minus, lam, t)
        av = [(vp[i]+vm[i])/2 for i in range(3)]
        d = Sv(av) - (Sv(vp)+Sv(vm))/2
        acc += d if d > 0 else mpf(0)
    return acc/len(TIMES)
h = mpf(1)/2; q = mpf(1)/4
ALONE  = ([(1,1)], [(1,-1)])
WITH_A = ([(h,2),(h,0)], [(h,0),(h,-2)])
WITH_AB= ([(q,3),(h,1),(q,-1)], [(q,1),(h,-1),(q,-3)])
print("="*104)
print("V4(a)  D RECOMPUTED AT 50 DIGITS AT THE COUPLINGS WHERE THE LANE'S D CHANGES SIGN")
print("="*104)
lanecol = {'0.08':'-7.7906e-05','0.1131':'-2.9051e-05','0.16':'+2.0090e-03',
           '0.2263':'+5.4024e-03','0.32':'-3.9694e-02','0.8':'-0.320493142660'}
print(f"  {'lam':>8}{'D exact (50 dp)':>26}{'D lane (float64)':>20}{'agree?':>10}")
for s, v in lanecol.items():
    lam = mpf(s)
    c0 = chi(*ALONE, lam); cA = chi(*WITH_A, lam); cAB = chi(*WITH_AB, lam)
    D = (c0-cAB) - 2*(c0-cA)
    ok = abs(D - mpf(v)) < abs(mpf(v))*mpf('1e-4') + mpf('1e-12')
    print(f"  {s:>8}{mp.nstr(D,12):>26}{v:>20}{str(ok):>10}")
print("  The lane's D is NOT monotone and changes sign twice inside the scanned range; the")
print("  'clean two-power gap' is read off a 6-point window at the bottom of that range only.")

# ---------------------------------------------------------------- (b)
print()
print("="*104)
print("V4(b)  TEST 6: THE HYPOTHESIS IT REJECTS IS FORBIDDEN BY chi >= 0, NOT BY THE PHYSICS")
print("="*104)
Q = {1:0.521527300760, 2:0.272817377944, 3:0.215349659532, 4:0.180650299961,
     5:0.167652293211, 6:0.145842259357, 7:0.141549402657, 8:0.124984198966}
d2 = Q[2] - 2*Q[1]
print(f"  {'m':>3}{'defect':>14}{'|defect| MAX POSSIBLE':>24}{'pairwise pred |.|':>20}"
      f"{'pred/max':>10}{'|defect|/max':>14}")
for m in range(2, 9):
    defect = Q[m] - m*Q[1]
    mx = m*Q[1]                       # because Q(m) >= 0
    pred = abs((m*(m-1)/2)*d2)
    print(f"  {m:>3}{defect:>14.6f}{mx:>24.6f}{pred:>20.6f}{pred/mx:>10.3f}{abs(defect)/mx:>14.3f}")
print("  For every m >= 3 the pairwise prediction EXCEEDS the largest defect chi >= 0 permits,")
print("  by up to 5.2x at m = 8.  And the measured defect sits at 97% of that ceiling.  A test")
print("  whose alternative hypothesis is arithmetically impossible does not discriminate forms.")
print("  ALSO: defect(m) = -m*Q(1) + Q(m) with 0 <= Q(m) <= Q(1) is LINEAR IN m by construction;")
print("  the fitted slope 0.5439 is just Q(1) = 0.5215 plus the slow drift of Q(m).")

# ---------------------------------------------------------------- (c)
print()
print("="*104)
print("V4(c)  INDEPENDENT EXACT RECHECK OF THE F_2 / INTEGER LAYER")
print("="*104)
print(f"  {'m':>4}{'n=4m':>6}{'#records (2m?)':>16}{'F_2 rank of Gram':>18}{'additivity defect':>19}")
ok = True
for m in (1,2,3,4,6,8,10,12):
    n = 4*m
    st = stab_blocks(m, 4)
    S, L, pairs = derived_logical_span(st, n)
    k = len(pairs)
    G = [[sp(a, b, n) for (b, _) in pairs] + [sp(a, d, n) for (_, d) in pairs] for (a, _) in pairs]
    G2 = [[sp(x, y, n) for y in [p for pr in pairs for p in pr]] for x in [p for pr in pairs for p in pr]]
    r = f2_rank([row[:] for row in G2], 2*k)
    print(f"  {m:>4}{n:>6}{2*k:>16}{r:>18}{(2*k - m*2*1) - 0:>19}")
    ok &= (2*k == 2*m) and (r == 2*m)
print(f"  every m: #records = 2m and Gram rank = 2m (full)?  {ok}")
print("  -> the additivity of the topological counts is CONFIRMED exactly.  It is also forced:")
print("     the stabiliser group of a direct sum of blocks IS the direct sum, so k, rank and")
print("     log2(dim) are additive by definition of the construction.  A quantity that is")
print("     additive by definition cannot be evidence that a SOURCE is additive.")
print()
n = 12
S12 = stab_nn2(n); Sb, Lb, _ = derived_logical_span(S12, n)
bad = 0
for i in range(n):
    for j in range(i+1, n):
        for p in range(n):
            for qq in range(p+1, n):
                a = pauli_vec(n,(i,j),()); b = pauli_vec(n,(),(p,qq))
                if len({i,j} & {p,qq}) == 0 and sp(a,b,n) != 0: bad += 1
print(f"  n=12 exhaustive over ALL weight-2 X/Z pairs with disjoint support: "
      f"anticommuting cases = {bad}  (must be 0; it is the qubit-wise sum, i.e. a definition)")
