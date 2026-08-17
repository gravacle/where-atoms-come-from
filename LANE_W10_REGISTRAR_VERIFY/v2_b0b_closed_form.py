# Registrar's verification 2 — S4:599 says B0b's U "is genuinely 4-term and does not factor.
# QUADRATURE ONLY, -0.810930216216". Both halves of that need checking.
from fractions import Fraction as F
import numpy as np
a, b, c, d = F(4,9), F(2,9), F(1,9), F(2,9)      # B0b SENSE U pushforward (p00,p10,p01,p11)
print("== DOES B0b's P FACTOR? ==")
print(f"  P = {a} + {b}x + {c}y + {d}xy      p00*p11 = {a*d},  p10*p01 = {b*c}")
print(f"  factors iff p00*p11 == p10*p01 :  {a*d == b*c}   -> S4 IS RIGHT, it does not factor.\n")
print("== BUT DOES IT HAVE A CLOSED FORM ANYWAY? ==")
print("  Jensen in y:  m(P) = (1/2pi) INT log max(|a+b e^{it}|, |c+d e^{it}|) dt.")
print(f"  |a+b e^it|^2 = a^2+b^2+2ab cos t = ({a*a+b*b}) + ({2*a*b}) cos t")
print(f"  |c+d e^it|^2 = c^2+d^2+2cd cos t = ({c*c+d*d}) + ({2*c*d}) cos t")
diffA, diffB = (a*a+b*b)-(c*c+d*d), 2*a*b - 2*c*d
print(f"  difference   = {diffA} + {diffB} cos t   >= 0  iff  cos t >= {-diffA/diffB}")
print(f"  and -diffA/diffB = {float(-diffA/diffB):.6f} <= -1, so THE FIRST BRANCH DOMINATES EVERYWHERE.")
print(f"  Hence m(P) = m(a + b x) = log(max(a,b)) = log({max(a,b)}) by Jensen again.\n")
lam = np.log(float(max(a,b)))
print(f"  CLOSED FORM  lambda = log(4/9) = {lam:.15f}")
print(f"  S4:599's quadrature value      = -0.810930216216")
print(f"  agreement to {abs(lam - -0.810930216216):.1e}")
# independent high-resolution Jensen check
n = 1<<22; t = 2*np.pi*np.arange(n)/n; ct = np.cos(t)
A = np.sqrt(float(a*a+b*b) + float(2*a*b)*ct); B = np.sqrt(float(c*c+d*d) + float(2*c*d)*ct)
print(f"  independent Jensen quadrature at n=2^22 = {np.log(np.maximum(A,B)).mean():.15f}")
print(f"  branches ever cross?  {bool((B > A).any())}")
print()
print("  ==> S4:599's 'QUADRATURE ONLY' IS WRONG. lambda = log(4/9) EXACTLY.")
print("      Non-factoring does not imply no closed form: one Jensen branch dominating suffices.")
