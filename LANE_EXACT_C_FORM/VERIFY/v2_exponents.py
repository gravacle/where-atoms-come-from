"""V2 -- THE SUPERPOSITION EXPONENTS, RECOMPUTED IN 60-DIGIT ARITHMETIC.

The lane fits chi ~ lam^1.9974, Delta_A ~ lam^3.9811, D ~ lam^5.9159 on SIX float64 points over
lam = 0.005..0.02828, and says the exponents "land on the integers 2, 4, 6".  Two objections:
  (1) two free parameters against six points, over 0.75 of a decade;
  (2) the float64 window is bounded BELOW by their own 1e-15 floor (at lam = 0.005, |D| is
      1.7e-11 and falls as lam^6, so lam = 0.0016 already sits on the floor).  The exponent can
      therefore never be pinned in float64.

Their venue reduces EXACTLY (V1) to one bath qubit driven by h_K = Z + lam*K*X with
K = s_C + s_A + s_B.  That model is analytic, so I redo it at 60 digits and take LOCAL log-log
slopes down to lam = 1e-8, four orders below anything float64 can reach.
"""
from mpmath import mp, mpf, cos, sin, sqrt, log, tanh, matrix
mp.dps = 60

BETA = mpf(2); Z0 = -tanh(BETA)          # Bloch z of the thermal bath qubit (H_b = Z)
TIMES = [mpf(1) + (mpf(12)*i)/24 for i in range(25)]

def bloch(K, lam, t):
    """Bloch vector of the bath qubit after time t in sector with integer K."""
    r = sqrt(1 + (lam*K)**2); a = lam*K/r; c = mpf(1)/r; th = 2*t*r
    C = cos(th); S = sin(th)
    return (a*c*Z0*(1-C), -a*Z0*S, Z0*C + c*c*Z0*(1-C))

def H2(x):
    if x <= 0 or x >= 1: return mpf(0)
    return -(x*log(x) + (1-x)*log(1-x))/log(2)

def Sv(v):
    n = sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return H2((1+n)/2)

def mix(pairs, lam, t):
    """weighted average of Bloch vectors: pairs = [(weight, K), ...]"""
    v = [mpf(0)]*3
    for w, K in pairs:
        b = bloch(K, lam, t)
        for i in range(3): v[i] += w*b[i]
    return v

def chi(plus, minus, lam):
    """time-averaged Holevo chi for the two readout branches given as sector mixtures"""
    acc = mpf(0)
    for t in TIMES:
        vp = mix(plus, lam, t); vm = mix(minus, lam, t)
        av = [(vp[i]+vm[i])/2 for i in range(3)]
        d = Sv(av) - (Sv(vp)+Sv(vm))/2
        acc += d if d > 0 else mpf(0)
    return acc/len(TIMES)

h = mpf(1)/2; q = mpf(1)/4
ALONE = ([(1, 1)], [(1, -1)])                                  # C only
WITH_A = ([(h, 2), (h, 0)], [(h, 0), (h, -2)])                 # C and A
WITH_AB = ([(q, 3), (h, 1), (q, -1)], [(q, 1), (h, -1), (q, -3)])  # C, A and B

def quantities(lam):
    c0 = chi(*ALONE, lam)
    cA = chi(*WITH_A, lam)
    cAB = chi(*WITH_AB, lam)
    dA = c0 - cA; dAB = c0 - cAB
    D = dAB - 2*dA
    return c0, dA, D

print("="*112)
print("V2(a)  CROSS-CHECK: 60-DIGIT ANALYTIC MODEL vs THE LANE'S FLOAT64 TABLE (t8_weakfield.txt)")
print("="*112)
lane = {mpf('0.005'): ('0.00006958165690', '2.4916e-08', '-1.7259e-11'),
        mpf('0.01'):  ('0.00027819662584', '3.9728e-07', '-1.0877e-09'),
        mpf('0.02'):  ('0.00111071783552', '6.2695e-06', '-6.5485e-08'),
        mpf('0.32'):  ('0.19933326983049', '8.8421e-02', '-3.9694e-02')}
print(f"  {'lam':>8}{'chi_C alone (60 dp)':>26}{'lane':>20}{'Delta_A (60 dp)':>22}{'lane':>14}"
      f"{'D (60 dp)':>16}{'lane':>14}")
for lam in sorted(lane):
    c0, dA, D = quantities(lam)
    print(f"  {float(lam):>8}{mp.nstr(c0,14):>26}{lane[lam][0]:>20}{mp.nstr(dA,8):>22}"
          f"{lane[lam][1]:>14}{mp.nstr(D,8):>16}{lane[lam][2]:>14}")

print()
print("="*112)
print("V2(b)  LOCAL LOG-LOG SLOPES, TAKEN FOUR ORDERS BELOW THE FLOAT64 FLOOR")
print("="*112)
lams = [mpf(10)**(-k) for k in (1, 2, 3, 4, 5, 6, 7, 8)]
vals = []
for lam in lams:
    vals.append((lam,) + quantities(lam))
print(f"  {'lam':>10}{'chi_C alone':>22}{'Delta_A':>22}{'D':>22}{'sign(D)':>9}")
for lam, c0, dA, D in vals:
    print(f"  {mp.nstr(lam,3):>10}{mp.nstr(c0,12):>22}{mp.nstr(dA,12):>22}{mp.nstr(D,12):>22}"
          f"{('+' if D>0 else '-'):>9}")
print()
print(f"  {'window':>20}{'slope chi':>14}{'slope Delta_A':>16}{'slope D':>14}")
for i in range(len(vals)-1):
    l0, c0, d0, D0 = vals[i]; l1, c1, d1, D1 = vals[i+1]
    dl = log(l1/l0)
    print(f"  {mp.nstr(l0,2)+' -> '+mp.nstr(l1,2):>20}"
          f"{mp.nstr(log(c1/c0)/dl,8):>14}{mp.nstr(log(d1/d0)/dl,8):>16}"
          f"{mp.nstr(log(D1/D0)/dl,8):>14}")
print()
print("  VERDICT ON THE EXPONENTS: read the last rows.  If they converge to 2, 4, 6 the lane's")
print("  numbers are right but the *interpretation* (a weak field superposing) is the ordinary")
print("  perturbative statement: chi is O(lam^2) per source, the first cross term is O(lam^4),")
print("  the first genuinely-three-body term is O(lam^6).  Any weakly coupled channel does this.")
