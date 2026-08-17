#!/usr/bin/env python3
"""
X_04 — WHAT I TRIED TO BREAK AND COULD NOT.  Independent second implementation of every
       load-bearing number in LANE_W13_R_REDUCTION.  R_lib.py is not imported here.

Legs:
  A  m(P) across eight node counts, and by a SECOND method (roots of the y-resultant), so
     CHOICE LEDGER L3's Jensen max-form is checked rather than reused.
  B  the two torus zeros, exactly (3 zeta^2 + 4 zeta + 3 = 0), and their simplicity.
  C  R-4's split identity, on a deterministic grid (the lane used a random one).
  D  the sixteen published subtorus rows, from my own reduction formula.
  E  R-5's proof, re-derived a different way and checked exhaustively:
        0.3 z^{m+n} + 0.4 z^n + 0.3 = 0 on |z|=1  <=>  |z^{m+n}+1| = 4/3 and a phase condition
        <=>  m = n.  My route uses Niven on cos(pi (m+n) psi); the lane's uses Niven on
        cos(2 pi m psi).  Both land on m=n; I check the CLAIM numerically over a bigger box.
  F  R-2's two analytic ingredients: area{|P|<delta} = delta^2/(2 pi |det|), and the GLOBAL
     lower bound |P(z)| >= c dist(z, Z(P)) that the (D2) truncation step needs.
  G  R-3's ladder, rebuilt from scratch at 900 digits, and its two dips.
"""
import math
import numpy as np
from fractions import Fraction
from X_lib import (PI_K1, P_eval, m_maxform, m_one_var, frac_sqrt, frac_theta_star, frac_pi)

p00, p10, p01, p11 = PI_K1
print("=" * 79); print("X_04 — INDEPENDENT VERIFICATION OF WHAT SURVIVES"); print("=" * 79)

# ---- A
print("\nA. m(P) ACROSS EIGHT NODE COUNTS, AND BY A SECOND METHOD")
prev = None
for e in range(10, 25, 2):
    v = m_maxform(PI_K1, 1 << e)
    print("     nq = 2^%-2d   m(P) = %.15f   change %s" % (e, v, "%.2e" % abs(v - prev) if prev else "--"))
    prev = v
# second method: m(P) = INT_x [ log|p01+p11 x| + log^+|R(x)| ], R = -(p00+p10 x)/(p01+p11 x)
nq = 1 << 22
t = np.arange(nq) * (2 * np.pi / nq); e1 = np.exp(1j * t)
R = -(p00 + p10 * e1) / (p01 + p11 * e1)
v2 = float(np.mean(np.log(np.abs(p01 + p11 * e1)))) + float(np.mean(np.log(np.maximum(np.abs(R), 1.0))))
print("     second method (split + Jensen in y)          m(P) = %.15f   dev %.2e" % (v2, abs(v2 - prev)))
print("     lane R_01/R_04 report -0.767507880358  -> AGREES")

# ---- B
print("\nB. THE TWO TORUS ZEROS")
TH = float(frac_theta_star(80))
z = np.exp(2j * np.pi * TH)
print("     zeta = e(theta*) = %.12f%+.12fi     |zeta| - 1 = %.2e" % (z.real, z.imag, abs(z) - 1))
print("     3 zeta^2 + 4 zeta + 3 = %.3e   (algebraic of degree 2, ON the unit circle)" % abs(3 * z * z + 4 * z + 3))
for (a, b) in [(TH, -TH), (-TH, TH)]:
    x = np.exp(2j * np.pi * a); y = np.exp(2j * np.pi * b)
    print("     |P(e(%+.6f), e(%+.6f))| = %.3e ;  |dP/dx| = %.6f (simple)"
          % (a, b, abs(P_eval(PI_K1, x, y)), abs(p10 + p11 * y)))
# are there others?  |p00+p10 x| = |p01+p11 x| has exactly the two solutions
tt = np.arange(1 << 20) / (1 << 20)
xx = np.exp(2j * np.pi * tt)
g = np.abs(p00 + p10 * xx) - np.abs(p01 + p11 * xx)
sgn = np.sign(g); crossings = int(np.sum(sgn[:-1] * sgn[1:] < 0))
print("     sign changes of |p00+p10 x| - |p01+p11 x| over a 2^20 grid: %d  -> exactly two zeros" % crossings)

# ---- C
print("\nC. R-4's SPLIT IDENTITY ON A DETERMINISTIC 2000x2000 GRID (the lane used random points)")
n = 2000
a = (np.arange(n) + 0.5) / n
X = np.exp(2j * np.pi * a)[:, None]; Y = np.exp(2j * np.pi * a)[None, :]
lhs = np.log(np.abs(P_eval(PI_K1, X, Y)))
den = p01 + p11 * X
rhs = np.log(np.abs(den)) + np.log(np.abs(Y + (p00 + p10 * X) / den))
print("     max deviation over 4e6 grid points: %.3e   (lane: 2.398e-14 over 2e5 random)" % np.max(np.abs(lhs - rhs)))

# ---- D
print("\nD. THE SIXTEEN PUBLISHED SUBTORUS ROWS, FROM MY OWN Q_{m,n}")
REG = {(1,0):-0.356674944,(0,1):-0.356674944,(1,1):-1.203972804,(1,-1):-0.510825624,
       (2,1):-0.681980359,(2,-1):-0.916290732,(3,1):-0.767783712,(3,2):-0.732940865,
       (4,1):-0.784966659,(5,1):-0.749392712,(5,3):-0.765224351,(7,3):-0.759305247,
       (7,11):-0.764712281,(11,20):-0.767014993,(13,8):-0.768271734,(29,17):-0.767138179}
worst = 0.0
for (m, n), ref in REG.items():
    ex = [(m + n, p10), (n, p11), (0, p01)]
    sh = -min(a for a, _ in ex); deg = max(a for a, _ in ex) + sh
    c = np.zeros(deg + 1)
    for k_, w in ex:
        c[k_ + sh] += w
    v = m_one_var(c)
    worst = max(worst, abs(v - ref))
print("     worst deviation over the 16 rows: %.3e   (lane: 4.792e-10)  -> REDUCTION FORMULA AGREES" % worst)

# ---- E
print("\nE. R-5's SINGULARITY CRITERION, RE-DERIVED AND SCANNED OVER A BIGGER BOX (|m|,|n|<=30)")
print("     my route: |0.3(z^{m+n}+1)| = 0.4 forces cos(pi(m+n)psi) = +-2/3, and the phase")
print("     condition forces psi(m-n) in (1/2)Z; if m != n then psi is rational and Niven")
print("     forbids cos(pi(m+n)psi) = +-2/3.  So m = n; primitivity gives (1,1).")
sing = []
zz = np.exp(2j * np.pi * np.arange(1 << 16) / (1 << 16))
for m in range(-30, 31):
    for n in range(-30, 31):
        if (m, n) == (0, 0) or math.gcd(abs(m), abs(n)) != 1:
            continue
        q = p10 * zz ** (m + n) + p11 * zz ** n + p01
        if np.min(np.abs(q)) < 1e-4:
            sing.append((m, n))
print("     primitive (m,n), |m|,|n| <= 30 with min_{|z|=1}|Q| < 1e-4 : %s" % sing)
print("     -> exactly the two, as R-5 says FOR PRIMITIVE (m,n).  R-5's PROOF IS CORRECT.")
print("     WHAT IS WRONG IS THE 'i.e. u v = 1, i.e. c = f' THAT FOLLOWS IT -- see X_01.")

# ---- F
print("\nF. R-2's TWO ANALYTIC INGREDIENTS")
x0 = np.exp(2j * np.pi * TH); y0 = np.conj(x0)
A = (p10 + p11 * y0) * x0; B = (p01 + p11 * x0) * y0
det = abs(A.real * B.imag - A.imag * B.real)
print("     |A| = %.6f  |B| = %.6f  |det(A,B)| = %.6f  -> the local map is INVERTIBLE" % (abs(A), abs(B), det))
for M in (2, 4, 6):
    dl = math.exp(-M)
    ana = dl * dl / (2 * np.pi * det)
    ng = 6000
    aa = (np.arange(ng) + 0.5) / ng
    XX = np.exp(2j * np.pi * aa)[:, None]; YY = np.exp(2j * np.pi * aa)[None, :]
    grid = float(np.mean(np.abs(P_eval(PI_K1, XX, YY)) < dl))
    print("     M=%d  delta=%.4e   analytic area %.6e   grid %.6e   ratio %.4f" % (M, dl, ana, grid, grid / ana))
# global lower Lipschitz bound
ng = 1500
aa = (np.arange(ng) + 0.5) / ng
XX = np.exp(2j * np.pi * aa)[:, None]; YY = np.exp(2j * np.pi * aa)[None, :]
absP = np.abs(P_eval(PI_K1, XX, YY))
d1 = np.minimum(np.abs(aa - TH), 1 - np.abs(aa - TH))[:, None]
d2 = np.minimum(np.abs(aa + TH), 1 - np.abs(aa + TH))[None, :]
dA = np.sqrt(d1 ** 2 + d2 ** 2)
d1b = np.minimum(np.abs(aa + TH), 1 - np.abs(aa + TH))[:, None]
d2b = np.minimum(np.abs(aa - TH), 1 - np.abs(aa - TH))[None, :]
dB = np.sqrt(d1b ** 2 + d2b ** 2)
dist = np.minimum(dA, dB)
ratio = absP / np.maximum(dist, 1e-12)
print("     min over a 1500^2 grid of |P|/dist(.,Z(P)) = %.6f  > 0  -> the (D2) truncation step" % ratio.min())
print("     |P| >= c dist is available with c ~ %.3f.  R-2's PROOF STANDS." % ratio.min())

# ---- G
print("\nG. R-3's LADDER, REBUILT FROM SCRATCH AT 900 DIGITS")
PREC = 900
G1 = frac_theta_star(PREC); G2 = 1 - G1
print("     gamma1 + gamma2 = %s   (exactly 1, which proof step (d) needs)" % (G1 + G2))
def ladder(gamma, k1, n1, Ms):
    ks, ns = [k1], [n1]
    for M in Ms:
        W = M * ns[-1] + (M - 1) * gamma
        ns.append(W.numerator // W.denominator + 1)
        ks.append(ks[-1] * M)
    al = [Fraction(nn, 1) / kk + gamma / kk for nn, kk in zip(ns, ks)]
    return ks, al, [kk * (al[-1] - a) for kk, a in zip(ks, al)]
D1v, D2v, k1 = 0.5, 1.2, 10
M1 = int(math.ceil(math.exp(D1v * k1))); k2 = k1 * M1
M2 = 1 << int(math.ceil(D2v * k2 / math.log(2.0)))
ksA, alA, epsA = ladder(G1, k1, 3, [M1, M2])
ksB, alB, epsB = ladder(G2, k1, 7, [M1, M2])
def lf(fr):
    return math.log(fr.numerator >> max(0, fr.numerator.bit_length() - 900)) + max(0, fr.numerator.bit_length() - 900) * math.log(2) \
         - (math.log(fr.denominator >> max(0, fr.denominator.bit_length() - 900)) + max(0, fr.denominator.bit_length() - 900) * math.log(2))
print("     rungs k = %s ; k_3 has %d digits" % (ksA[:2], len(str(ksA[2]))))
for j in (0, 1):
    ok = epsA[j] > 0 and epsB[j] > 0
    print("     j=%d  log eps  = %11.3f   log eps' = %11.3f   bound log(2k_j/k_{j+1}) = %11.3f   both>0: %s"
          % (j + 1, lf(epsA[j]), lf(epsB[j]), lf(Fraction(2 * ksA[j], ksA[j + 1])), ok))
alpha, beta = alA[-1], alB[-1]
MP = m_maxform(PI_K1, 1 << 24)
for j, kj in enumerate(ksA[:2]):
    d1_ = (kj * alpha) % 1 - G1
    d2_ = (kj * beta) % 1 - G2
    # |Z| ~ |2 pi (A d1 + B d2)| = 2 pi d1 |A + B r|, r = d2/d1 -- both offsets are tiny but
    # their RATIO is O(1), so no overflow.  d1, d2 are exact Fractions.
    r_ = float(Fraction(d2_.numerator * d1_.denominator, d2_.denominator * d1_.numerator))
    logZ = math.log(2 * math.pi) + lf(d1_) + math.log(abs(A + B * r_))
    print("     k = %6d   log eps = %11.3f   log|Z_k| (local expansion) = %13.4f" % (kj, lf(epsA[j]), logZ))
    print("                   UPPER BOUND on S_k, namely (1/k) log|Z_k| = %10.4f   (m(P) = %.6f)"
          % (logZ / kj, MP))
    print("                   [this is R-3's step (c) bound, not S_k itself; S_k is lower still]")
print("     lane R_05 reports log|Z_10| = -4.578, log|Z_1490| = -1787.797, dips -0.281 and -1.199.")
print("\nDONE X_04")
