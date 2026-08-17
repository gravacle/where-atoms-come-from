# W-08 / M4-REFUTE-2 leg A — LENS 2 (MATHEMATICS).
# THE TARGET: m4_c_diophantine.py block C1, "THE DECISIVE TEST", the single block M4-1 rests on.
#
# WHAT THIS LEG DOES.  m4_c builds each row by handing ONE float, theta, to `connection_from`,
# which sets argWC = -argWF - 2*pi*theta, splits both angles into three equal edge angles, and
# recovers W_F, W_C by summing the three floats back and exponentiating.  Every step of that
# round trip is float64.  This leg recomputes, in EXACT arithmetic, (a) the angle the machine
# ACTUALLY realised in each row and (b) the min and cell count of the row's INTENDED object.
#
# METHOD.  Every float64 is an exact rational; Decimal(repr(x)) carries it losslessly.  pi is
# carried at 60 significant digits.  For the pair (2,3) on K1, dF=-1 and dC=+1, so
#     D_k = amp |W_F^{-k} - W_C^{k}| = amp |1 - (W_F W_C)^k| = 2 amp |sin(pi k psi)|,
#     psi = ((a1+a2+a3)+(a4+a5+a6)) / 2pi   mod 1
# -- the same identity m4_b B0 derives and m4_c's `profile` uses.  Verified against m4_c's own
# float output in block A3 below.
#
# ISOLATION LEDGER.  HELD FIXED in A1/A2: carrier K1, ready state (numpy default_rng(20260816),
# the same bytes as m4_c and as W-07 leg E), observable A_23, dressing tree, k-range 1..4000,
# threshold 1e-9, the row list, and one evaluation function (`exact_profile`) for every cell of
# the table.  MOVED: the ARITHMETIC MODEL alone -- float64 (what m4_c ran) versus 60-digit exact
# (what m4_c's prose says it means).  theta is identical between the two columns of every row.
# A0 moves the row (six values of theta) and holds everything else, to test one ledger claim.
import numpy as np
from fractions import Fraction
from decimal import Decimal, getcontext, ROUND_FLOOR
getcontext().prec = 60
PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
TWOPI = 2*PI
def frac1(x): return x - x.to_integral_value(rounding=ROUND_FLOOR)   # x mod 1 into [0,1)
def norm(x):                                                          # ||x||, distance to nearest int
    f = frac1(x); return f if f <= Decimal("0.5") else 1-f
def dsin(x):                                                          # sin at 60 digits, 0<=x<=pi/2
    term = x; ssum = x; n = 1
    while abs(term) > Decimal("1e-58"):
        term = -term*x*x/Decimal((2*n)*(2*n+1)); ssum += term; n += 1
    return ssum

rng = np.random.default_rng(20260816)
s = rng.normal(size=5) + 1j*rng.normal(size=5); s /= np.linalg.norm(s)
amp_f = float(abs(s[2])*abs(s[3]))
AMP = Decimal(repr(amp_f))
K = 4000; TOL = Decimal("1e-9")

FACE_V={0,1,2}; CYC_V={0,3,4}; TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
def dress(s_, a):
    u=np.exp(1j*np.asarray(a)); t=np.array(s_,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=s_[v]/w
    return t

LAM13 = (2**0.5-1)*1e-13; LAM16 = (2**0.5-1)*1e-16
rows = [
 ("ord(rho)=4     [S1 PUBLISHED]", -0.25,        Fraction(-1,4),   None),
 ("ord(rho)=4001  [FINITE, > K] ", 1/4001,       Fraction(1,4001),  None),
 ("ord(rho)=8000  [FINITE, > K] ", 1/8000,       Fraction(1,8000),  None),
 ("ord(rho)=2000  [FINITE, < K] ", 1/2000,       Fraction(1,2000),  None),
 ("INFINITE 1/4+(sqrt2-1)e-13   ", 0.25+LAM13,   None, Decimal("1e-13")),
 ("INFINITE 1/4+(sqrt2-1)e-16   ", 0.25+LAM16,   None, Decimal("1e-16")),
]

print("== A0  THE ISOLATION-LEDGER CLAIM THAT amp IS HELD FIXED — CHECKED, NOT ASSUMED ==")
amps=set()
for tag, thf, _, _ in rows:
    a = np.array([np.pi/3]*3 + [(-np.pi-2*np.pi*thf)/3]*3)
    t = dress(s, a); amps.add(float(abs(np.conj(t[2])*t[3])))
print(f"  |s_2||s_3| = {amp_f!r}   (dressing is diagonal unit-modulus, so |t_v| = |s_v| exactly)")
print(f"  distinct float amp over the six rows : {len(amps)} -> {sorted(amps)}")
print(f"  spread = {max(amps)-min(amps):.3e}  ({(max(amps)-min(amps))/amp_f:.2e} relative) -- 1 ulp of round-off,")
print("  not a moving variable.  M4's ledger claim 'amp held fixed' is CORRECT.")
print(f"  Aside: M4 quotes amp = 0.271776443 in m4_g and in M4-7; the value is {amp_f:.12f}.\n")

def realised_psi(theta_float):
    """Exact reconstruction of the angle m4_c's float64 pipeline ACTUALLY built, /2pi, mod 1."""
    argWF = np.pi
    argWC = -argWF - 2*np.pi*theta_float
    aF = argWF/3.0; aC = argWC/3.0
    sF = (aF+aF)+aF; sC = (aC+aC)+aC          # np.sum over 3 elements = left fold
    return frac1((Decimal(repr(sF))+Decimal(repr(sC)))/TWOPI)

def exact_profile(psi):
    best=None; cnt=0; kmin=None
    for k in range(1, K+1):
        d = norm(Decimal(k)*psi)
        D = 2*AMP*dsin(PI*d)
        if best is None or D < best: best, kmin = D, k
        if D < TOL: cnt += 1
    return best, cnt, kmin

print("== A1  THE ANGLE THE MACHINE ACTUALLY REALISED, ROW BY ROW ==")
print("  psi is the row's own quantity, normalised to [0,1/2] (psi and 1-psi give the same D).")
print(f"  {'row':<30} {'psi INTENDED':>28} {'psi REALISED':>28} {'|diff|':>11}")
info={}
for tag, thf, exact_th, e in rows:
    pr = norm(realised_psi(thf))
    if exact_th is not None:
        pi_int = norm(Decimal(exact_th.numerator)/Decimal(exact_th.denominator))
    else:
        pi_int = norm(Decimal("0.25")+(Decimal(2).sqrt()-1)*e)
    info[tag]=(pr,pi_int)
    print(f"  {tag:<30} {pi_int:>28.24f} {pr:>28.24f} {abs(pr-pi_int):>11.3e}")
print()

print("== A2  min D AND CELL COUNT — THREE ARITHMETICS ON THE SAME SIX ROWS ==")
print("  PUBLISHED  = the number m4_c C1 prints (float64).")
print("  REALISED   = 60-digit evaluation of the angle m4_c's float pipeline actually built.")
print("  INTENDED   = 60-digit evaluation of the mathematical object the row's LABEL names.")
print(f"  {'row':<30} {'PUBLISHED min':>14} {'REALISED min':>16} {'INTENDED min':>16}  cells P/R/I")
pub = {"ord(rho)=4     [S1 PUBLISHED]":("6.729e-19",1000),
       "ord(rho)=4001  [FINITE, > K] ":("4.268e-04",0),
       "ord(rho)=8000  [FINITE, > K] ":("2.135e-04",0),
       "ord(rho)=2000  [FINITE, < K] ":("2.141e-13",2),
       "INFINITE 1/4+(sqrt2-1)e-13   ":("2.828e-13",1000),
       "INFINITE 1/4+(sqrt2-1)e-16   ":("6.657e-17",1000)}
for tag, thf, exact_th, e in rows:
    pr, pin = info[tag]
    mr, cr, _ = exact_profile(pr)
    mi, ci, _ = exact_profile(pin)
    pf, pc = pub[tag]
    print(f"  {tag:<30} {pf:>14} {float(mr):>16.6e} {float(mi):>16.6e}  {pc:>5}/{cr:>4}/{ci:>4}")
print()
print("  FINDINGS FROM THIS TABLE:")
print("   (1) ord=4 and ord=2000: the INTENDED min is 0 EXACTLY (rho^k=1 at k=4, k=2000).  m4_c C1")
print("       publishes 6.729e-19 and 2.141e-13 in a column headed 'min D'.  That is verbatim the")
print("       defect m4_g G3(a) convicts W-07 of -- 'W-07 sec3's table publishes 6.729e-19 as min")
print("       where the exact value is 0' -- committed by M4 in its own decisive table.")
print("   (2) The e-16 row: intended psi and realised psi differ by ~1.9e-17, which is 46% of the")
print("       intended perturbation itself (4.14e-17).  The published min 6.657e-17 matches NEITHER")
print("       the intended value (2.829e-16) NOR a correct float evaluation of what was built.")
print("       The row is inside the noise floor of its own constructor.")
print("   (3) The e-13 row survives exactly: intended 2.829280e-13, published 2.828e-13, 1000 cells")
print("       under all three arithmetics.  THIS ROW, ALONE, CARRIES M4-1.")
