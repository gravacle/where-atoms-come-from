# S_05 — IS ord(rho) THE OPERATIVE VARIABLE OF W-07's OWN TABLE?
# W-07 leg E's observable is D_k = amp*|rho^k - 1|, K = 4000, threshold 1e-9, amp from its own
# seeded state.  I rebuild amp from W-07's own code path, then move ONE coordinate at a time.
# EXACT arithmetic (Fraction) for every rational theta; the sub-threshold counts are integers.
import numpy as np
from fractions import Fraction
from decimal import Decimal, getcontext
getcontext().prec = 60

# --- amp, rebuilt exactly as w07_e_isolation.py builds it (seed 20260816, pair (u,v)=(2,3)) ---
rng = np.random.default_rng(20260816)
FACE_V={0,1,2}; CYC_V={0,3,4}; TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
s = rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)
a = np.array([np.pi/3]*3+[np.pi/2]*3)
u_ = np.exp(1j*a); t = np.array(s,dtype=complex)
for v,p in TREE.items():
    w=1.0+0j
    for e in p: w*=u_[e]
    t[v]=s[v]/w
amp = abs(np.conj(t[2])*t[3])
print(f"amp (W-07 leg E, pair (2,3), seed 20260816) = {amp:.15f}")
K, TOL = 4000, 1e-9

def count_exact(theta_frac):
    """theta an exact Fraction.  D_k = amp*2|sin(pi k theta)| in 60-digit decimal."""
    n = 0; mn = None
    pi = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
    for k in range(1, K+1):
        r = theta_frac*k
        fr = r - int(r)
        if fr > Fraction(1,2): fr = 1-fr
        if fr < 0: fr = -fr
        d = Decimal(2)*Decimal(amp)*( (pi*Decimal(fr.numerator)/Decimal(fr.denominator)).__abs__() ).exp().ln()*0  # placeholder
        # use sin directly at 60 digits via the series-free identity sin(x) ~ mpmath-free:
        x = pi*Decimal(fr.numerator)/Decimal(fr.denominator)
        # Taylor for sin on [0, pi/2]; 25 terms is far more than enough at 60 dp
        term = x; ssum = x
        for j in range(1, 25):
            term = -term*x*x/Decimal((2*j)*(2*j+1)); ssum += term
        d = Decimal(2)*Decimal(amp)*ssum
        if mn is None or d < mn: mn = d
        if d < Decimal(TOL): n += 1
    return n, mn

def order_of(theta_frac):
    return theta_frac.denominator   # ord(e^{2 pi i theta}) = denominator in lowest terms

print()
print("== S5a  W-07's OWN ROW, AND A FINITE-ORDER COUNTEREXAMPLE TO ITS OWN SHARP FORM ==")
print("   W-07 sec3: '1000 = 4000/4 is ord(rho)' -> the sharp form is count = floor(K/ord).")
rows = [
  ("S1 PUBLISHED  theta = -1/4          ", Fraction(-1,4)),
  ("theta = 1/4 + 1e-13   (EXACT rational)", Fraction(1,4)+Fraction(1,10**13)),
  ("theta = 1/1000 + 4e-16 (EXACT rational)", Fraction(1,1000)+Fraction(4,10**16)),
  ("theta = 1/2                         ", Fraction(1,2)),
  ("theta = 1/2000                      ", Fraction(1,2000)),
]
print(f"   {'theta':<40} {'ord(rho)':>16} {'floor(K/ord)':>13} {'ACTUAL count<1e-9':>18} {'min D':>12}")
for tag, th in rows:
    o = order_of(th); n, mn = count_exact(th)
    print(f"   {tag:<40} {o:>16d} {o<=K and K//o or 0:>13d} {n:>18d} {float(mn):>12.3e}")
print("   ROW 2 IS THE KILL:  ord(rho) = 10^13 > K, so W-07's own formula predicts 0 exact zeros")
print("   -- correct, there are none -- while the observable W-07 TABULATES returns 1000.")
print("   Irrationality was never needed.  'finite versus infinite' was never the axis.")

print()
print("== S5b  ROWS 2 AND 3 HAVE IDENTICAL min D AND COUNTS 1000 vs 4 ==")
print("   min_k ||k theta|| is held to the last digit (4*1e-13 = 1000*4e-16 = 4e-13) and the")
print("   count moves by 250x.  M4's OWN named variable, min_k ||k theta||, does not carry it.")
print("   The count is carried by the DENOMINATOR q of the nearby rational: count = floor(K/q).")
