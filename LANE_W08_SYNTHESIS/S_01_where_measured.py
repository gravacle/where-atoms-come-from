# S_01 — WHERE WERE THE REGISTER'S RECURRENCE FIGURES MEASURED?
# W-07 sec1/sec4 and M4's headline both assert every recurrence figure in the corpus was measured
# at S1's published order-4 connection.  M4 refuter 1 says that is false.  Deciding it here.
# Double precision; the S1PUB values are re-checked exactly at the end.
import numpy as np
from fractions import Fraction

def Zabs(f, c, w, K):
    k = np.arange(1, K+1, dtype=np.float64)
    z = w[0]*np.exp(1j*k*(c-f)) + w[1]*np.exp(-1j*k*f) + w[2]*np.exp(1j*k*c)
    return np.abs(z)

RSG = (0.4, 0.3, 0.3)          # p=(0.4,.15,.15,.15,.15) -- S3 audit's own test-point state
RSP = (0.5, 0.0, 0.5)          # S1's published ready state

print("== S1a  THE RESONANT POINT f=2.0, c=1.1, RS-G, k<=200000 ==")
K = 200000
d = Zabs(2.0, 1.1, RSG, K)
i = int(np.argmin(d)); print(f"   min |Z_k| = {d.min():.6f}   at k = {i+1}")
# running record-breakers of the max
best = -1.0; recs = []
for j in range(K):
    if d[j] > best: best = d[j]; recs.append((j+1, d[j]))
print(f"   record-breakers of max|Z_k| to k<=200000 : {len(recs)}")
for kk, vv in recs: print(f"       k = {kk:>7}   |Z_k| = {vv:.9f}")
print(f"   W-01 registered '0.0247 at n=42' and 'recurs to 0.99994'.")
print(f"   S3 audit COR-E registered sup_(k>=42)|Z_k| = 0.999941 and, over k<=200000, 0.999999981 at k=106123.")

print()
print("== S1b  THE SAME OBSERVABLE AT S1'S PUBLISHED ORDER-4 CONNECTION, RS-G ==")
d4 = Zabs(np.pi, 3*np.pi/2, RSG, 10**6)
vals = np.unique(np.round(d4, 12))
print(f"   distinct values of |Z_k| over k<=1e6 : {len(vals)}  ->  {vals}")
print(f"   min = {d4.min():.9f}   max = {d4.max():.15f}")
print(f"   is 0.024654 attained at any k<=1e6 ?  min |  |Z_k| - 0.024654 | = {np.abs(d4-0.024654).min():.6e}")
print(f"   is 0.999941 attained at any k<=1e6 ?  min |  |Z_k| - 0.999941 | = {np.abs(d4-0.999941).min():.6e}")
# exact check of the four order-4 values in Gaussian rationals: x=conj(-1)=-1, y=-i, xy=i
print("   EXACT (Gaussian rationals, x=-1, y=-i, xy=i, w=(2/5,3/10,3/10)):")
p11, p10, p01 = Fraction(2,5), Fraction(3,10), Fraction(3,10)
for k in range(1, 5):
    xk  = (-1)**k
    yk  = [1, -1j, -1, 1j][k % 4]     # (-i)^k
    xyk = [1, 1j, -1, -1j][k % 4]     # (i)^k
    re = p11*int(np.real(xyk)) + p10*xk + p01*int(np.real(yk))
    im = p11*int(np.imag(xyk)) + p01*int(np.imag(yk))
    print(f"      k={k}:  Z_k = {re} + {im} i     |Z_k|^2 = {re*re+im*im}")

print()
print("== S1c  RS-P (S1's OWN PUBLISHED READY STATE) AT THE RESONANT POINT AND AT S1PUB ==")
for tag, f, c in [("S1PUB order-4  ", np.pi, 3*np.pi/2), ("RESONANT f=2,c=1.1", 2.0, 1.1)]:
    dd = Zabs(f, c, RSP, 10**6)
    print(f"   {tag}: min={dd.min():.9f}  max={dd.max():.15f}  cells |Z|>1-1e-12 : {int((dd>1-1e-12).sum())}")
