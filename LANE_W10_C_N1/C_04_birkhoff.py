#!/usr/bin/env python3
"""
LANE W-10 / C — STEP 4.  DOES (1/N) sum log|Z_k| CONVERGE TO m(P)?
THREE CONNECTIONS, ON THE TWO FOUR-CLASS CARRIERS, WITH K1 AS THE THREE-CLASS CONTROL.

THE ONE VARIABLE across the three rows of each block: THE ARITHMETIC OF THE CONNECTION
(Diophantine / exactly resonant / finite order).  Ready state, carrier, evaluator, K, chunking
and checkpoints are byte-identical across the three.  The arms' (alpha,beta) are PRINTED so the
reader can see they differ.

PHASES.  k*alpha mod 1 is accumulated by EXACT int64 modular reduction against a common
denominator D, plus a float64 correction for the residue — naive k*alpha in float64 loses
~1e-9 by k = 1e7.  This is M1_03's device, reimplemented here.

Precision: float64 for the averages; the resonant and finite-order LIMITS are given in closed
form and printed beside them.
"""
import numpy as np
import mpmath as mp
import sys
from fractions import Fraction as Fr

mp.mp.dps = 40
K = 10 ** 7
CHECK = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
CHUNK = 10 ** 6
D = 2 ** 40

def orbit_average(p, A_num, B_num, Dn, dA=0.0, dB=0.0, K=K, checkpoints=CHECK):
    """alpha = A_num/Dn + dA, beta = B_num/Dn + dB.  Returns {N: average}, min|Z|, #zeros."""
    p00, p10, p01, p11 = p
    tot = 0.0; out = {}; minabs = np.inf; nzero = 0; done = 0
    while done < K:
        n = min(CHUNK, K - done)
        k = np.arange(done + 1, done + n + 1, dtype=np.int64)
        fa = np.mod(((k * A_num) % Dn).astype(np.float64) / Dn + k.astype(np.float64) * dA, 1.0)
        fb = np.mod(((k * B_num) % Dn).astype(np.float64) / Dn + k.astype(np.float64) * dB, 1.0)
        x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb)
        a = np.abs(p00 + p10 * x + p01 * y + p11 * x * y)
        minabs = min(minabs, float(a.min())); nzero += int(np.sum(a == 0.0))
        cs = np.cumsum(np.log(np.maximum(a, 1e-323)))
        for cp in checkpoints:
            if done < cp <= done + n:
                out[cp] = (tot + cs[cp - done - 1]) / cp
        tot += float(cs[-1]); done += n
    return out, minabs, nzero

def m_jensen(p, n=1 << 20):
    p00, p10, p01, p11 = [float(q) for q in p]
    t = 2 * np.pi * (np.arange(n) + 0.5) / n; c = np.cos(t)
    A2 = p00 * p00 + p10 * p10 + 2 * p00 * p10 * c
    B2 = p01 * p01 + p11 * p11 + 2 * p01 * p11 * c
    return float(np.mean(0.5 * np.log(np.maximum(A2, B2))))

def finite_order_limit(p, alphas, betas):
    """Exact average of log|P| over a finite orbit, in mpmath."""
    s = mp.mpf(0)
    for a, b in zip(alphas, betas):
        x = mp.expjpi(2 * a); y = mp.expjpi(2 * b)
        v = mp.mpf(p[0]) + mp.mpf(p[1]) * x + mp.mpf(p[2]) * y + mp.mpf(p[3]) * x * y
        s += mp.log(abs(v))
    return s / len(alphas)

def subtorus_limit(p, m, n, N=1 << 22):
    """H = closure of {(k a, k b)} when m*alpha + n*beta = 0, gcd(m,n)=1:
       H = {(n s, -m s)}, so the limit is a ONE-VARIABLE Mahler measure."""
    p00, p10, p01, p11 = [float(q) for q in p]
    s = 2 * np.pi * (np.arange(N) + 0.5) / N
    x = np.exp(1j * n * s); y = np.exp(-1j * m * s)
    return float(np.mean(np.log(np.abs(p00 + p10 * x + p01 * y + p11 * x * y))))

# ---------------------------------------------------------------- the ready states of record
CASES = [
    ("B0b  SENSE U  (4,2,1,2)/9   FOUR-CLASS", (4/9, 2/9, 1/9, 2/9), float(mp.log(mp.mpf(4)/9))),
    ("B4   SENSE U  (1,1,1,3)/6   FOUR-CLASS", (1/6, 1/6, 1/6, 3/6), float(-mp.log(2))),
    ("K1   SENSE U  (0,2,2,1)/5   three-class CONTROL", (0.0, 0.4, 0.4, 0.2), None),
    ("SENSE C 4-class (1,1,1,1)/4  FOUR-CLASS, HAS TORUS ZEROS", (.25, .25, .25, .25),
     float(-mp.log(4))),
]

# ---------------------------------------------------------------- the three connections
alphaA = -(2.0 ** (1.0 / 3.0)) % 1.0        # u = e^{2 pi i alpha};  alpha = -f/2pi
betaA = (4.0 ** (1.0 / 3.0)) % 1.0
A_numA = int(np.floor(alphaA * D)); B_numA = int(np.floor(betaA * D))
dAA = alphaA - A_numA / D; dBA = betaA - B_numA / D

# RESONANT.  The relation 11*alpha + 20*beta = 0 is made EXACT in the integer representation
# (M1_03's device, reimplemented): alpha = -1/pi, beta = 11/(20 pi).
S = 2 ** 35
A0 = int(round(S / np.pi))
D_C = 20 * S
A_numC = (-20 * A0) % D_C
B_numC = (11 * A0) % D_C
RELATION_RESIDUE = 11 * (-20 * A0) + 20 * (11 * A0)     # must be exactly 0
dAC = dBC = 0.0

CONNS = [
    ("DIOPHANTINE  alpha=-2^(1/3), beta=4^(1/3) mod 1  (Schmidt: badly approximable)",
     A_numA, B_numA, D, dAA, dBA, "m(P)"),
    ("RESONANT     f=2.0, c=1.1   ->  -11 f + 20 c = 0, primitive relation (11,20)",
     A_numC, B_numC, D_C, dAC, dBC, "subtorus"),
    ("FINITE ORDER W_F=-1, W_C=-i  ->  alpha=1/2, beta=3/4, ord(rho)=4  [S1 sec6]",
     2, 3, 4, 0.0, 0.0, "finite"),
]

if __name__ == "__main__":
    print("=" * 104)
    print("C_04 — THE THREE LIMITS.  K = 1e7, float64, exact int64 phase reduction.")
    print("ONE VARIABLE PER BLOCK: the connection's arithmetic.  The arms are printed; DIFF THEM.")
    print("=" * 104)
    print(f"\n  ARM PARAMETERS, so that a byte-identical control is impossible to hide:")
    for nm, An, Bn, Dn, dA, dB, kind in CONNS:
        print(f"    {nm[:60]:60s}  alpha = {An/Dn + dA:.15f}   beta = {Bn/Dn + dB:.15f}")
    a = [An / Dn + dA for nm, An, Bn, Dn, dA, dB, k in CONNS]
    b = [Bn / Dn + dB for nm, An, Bn, Dn, dA, dB, k in CONNS]
    assert len(set(np.round(a, 12))) == 3 and len(set(np.round(b, 12))) == 3, "ARMS COINCIDE"
    print("    all three (alpha,beta) pairs are distinct -> the control is not vacuous")
    print(f"    resonance exactness in the integer representation: 11*A + 20*B = "
          f"{RELATION_RESIDUE}  (must be 0)")

    for label, p, closed in CASES:
        mP = m_jensen(p)
        print("\n" + "-" * 104)
        print(f"  {label}")
        print(f"     m(P) by Jensen (n=2^20) = {mP:.12f}" +
              (f"   closed form {closed:.12f}   |diff| = {abs(mP-closed):.2e}"
               if closed is not None else "   [no closed form quoted]"))
        hi = (p[0] + p[1]) - (p[2] + p[3]); lo = abs(p[0] - p[1]) - abs(p[2] - p[3])
        print(f"     zeros on T^2 ?  {hi*lo <= 1e-15}   (hi={hi:+.4f}, lo={lo:+.4f})")
        for nm, An, Bn, Dn, dA, dB, kind in CONNS:
            out, minabs, nz = orbit_average(p, An, Bn, Dn, dA, dB)
            if kind == "m(P)":
                target, tname = mP, "m(P)"
            elif kind == "subtorus":
                target, tname = subtorus_limit(p, 11, 20), "subtorus m over H"
            else:
                target = float(finite_order_limit(p, [0, .5, 0, .5], [0, .75, .5, .25]))
                tname = "finite-orbit average"
            print(f"     {nm}")
            row = "   ".join(f"N=1e{int(np.log10(cp))}: {out[cp]:.9f}" for cp in CHECK)
            print(f"        {row}")
            print(f"        target ({tname}) = {target:.12f}   |avg(1e7) - target| = "
                  f"{abs(out[K]-target):.3e}   |avg(1e7) - m(P)| = {abs(out[K]-mP):.3e}"
                  f"   min|Z_k| = {minabs:.3e}  zeros hit: {nz}")
    print("\n" + "=" * 104)
    print("S4's OWN CROSS-CHECK CONNECTION, REPRODUCED.  S4:604-612 verifies its closed forms")
    print("against 'direct schedule-B simulation at f = 1.0, c = sqrt(2), N = 2e6' and reports")
    print("B0b direct -0.810929681 (dev 5.4e-07) and B4 direct -0.693146936 (dev 2.4e-07).")
    print("=" * 104)
    alS = (-1.0 / (2 * np.pi)) % 1.0            # f = 1.0  -> alpha = -f/2pi
    beS = (np.sqrt(2.0) / (2 * np.pi)) % 1.0    # c = sqrt2 -> beta = c/2pi
    DS = 2 ** 40
    AnS = int(np.floor(alS * DS)); BnS = int(np.floor(beS * DS))
    dAS = alS - AnS / DS; dBS = beS - BnS / DS
    for label, p, closed in CASES[:2]:
        out, mn, nz = orbit_average(p, AnS, BnS, DS, dAS, dBS, K=2 * 10 ** 6,
                                    checkpoints=[2 * 10 ** 6])
        v = out[2 * 10 ** 6]
        print(f"  {label[:38]:38s} N=2e6 direct = {v:.9f}   m(P) = {m_jensen(p):.9f}"
              f"   dev = {abs(v-m_jensen(p)):.1e}")
    print("  S4's B0b direct -0.810929681 and B4 direct -0.693146936 sit at the same distance")
    print("  from the closed form as this independent run; the finite-N wobble is the whole gap.")

    print("\n" + "=" * 104)
    print("""READ THIS BLOCK BEFORE THE NUMBERS ARE USED.

  ON B0b AND B4 THE THREE-CONNECTION COMPARISON IS A CONTROL THAT COULD NOT HAVE FAILED IN THE
  DIRECTION THE BRIEF EXPECTS.  Both carriers' P are ZERO-FREE on T^2 (C_03 section 4), so
  log|P| is CONTINUOUS there, and Weyl's theorem on the orbit closure H gives convergence for
  EVERY connection with no Diophantine input whatever.  The three limits therefore exist in all
  three rows; what separates them is only WHICH limit -- m(P) when H = T^2, the subtorus Mahler
  measure when H is a circle, the finite-orbit average when H is finite.  That is a statement
  about H, not about approximation quality, and it was already M1's T4.

  THE ROW WHERE SOMETHING COULD FAIL is SENSE C -- the corpus's OTHER published four-class
  column (lambda = -1.386294361120 for both B0b and B4, S4:582).  There P = (1+x)(1+y)/4 has a
  ONE-DIMENSIONAL zero set, log|P| is unbounded, and convergence is NOT free.  It is tested in
  C_05.""")
    sys.exit(0)
