#!/usr/bin/env python3
"""
LANE W-10 / C — REFUTER 2 — STEP 3.
DRIVE K UNTIL THE RULER IS FINER THAN THE EFFECT, AND SEE WHICH TARGET THE RESONANT ARM
ACTUALLY CONVERGES TO ON THE TWO FOUR-CLASS CARRIERS.

R2_02 showed the effect the resonant arm is claimed to detect is 1.17e-08 (B0b), 4.84e-10 (B4)
and EXACTLY 0 (SENSE C), against a K = 1e7 ruler of 5.4e-08.  This step re-runs the arm with a
finer ruler.  It also runs the connection S4 ITSELF published as its verification connection --
f = 1.0, c = sqrt(2), S4:603 and S4:1054 -- which lane C's C-4 headline asserts does not exist
("lambda = m(P) is FALSE at every connection this corpus has ever published").

PHASE ARITHMETIC -- INDEPENDENT OF LANE C's, AND STRICTLY FINER.
  Phases are carried as EXACT integers mod 2^80 in two int64 limbs of 40 bits, advanced by
  CHUNKED modular accumulation, so
    * there is no K ceiling (lane C's direct k*A_num wraps int64 at K = 1.13e7 -- R2_01 sec 3),
    * the represented rotation number is within 2^-81 of the named real, so the phase drift
      over the whole run is < 1e-15 turns even at K = 1e9 (lane C's float64 alpha drifts
      1.1e-09 turns by K = 1e7 -- R2_01 sec 4).

THE RESONANT ARM IS BUILT ON ITS SUBTORUS BY CONSTRUCTION, NOT BY CANCELLATION.  With the
primitive relation (m,n) = (11,20), H = {(z^20, z^-11)} and the orbit is z_k = e^{2 pi i k s},
s = -1/(20 pi):  20 s = -1/pi = alpha and -11 s = 11/(20 pi) = beta.  So 11 alpha + 20 beta = 0
holds identically in the representation, at every K, with no rounding to cancel.

ONE VARIABLE PER BLOCK: the connection's arithmetic.  Ready states, evaluator, chunking and
checkpoints are the same object across every row (one pass over k serves all four states).
Precision: float64 for the accumulated logs; phases exact as above; targets from R2_02 (mpmath).
"""
import numpy as np
import mpmath as mp
import time

mp.mp.dps = 50
B = 40
MASK = (1 << B) - 1
SCALE = 1 << B
CHUNK = 1 << 20

CASES = [
    ("B0b  (4,2,1,2)/9", np.array([4, 2, 1, 2], dtype=float) / 9),
    ("B4   (1,1,1,3)/6", np.array([1, 1, 1, 3], dtype=float) / 6),
    ("K1   (0,2,2,1)/5", np.array([0, 2, 2, 1], dtype=float) / 5),
    ("SENSE C (1,1,1,1)/4", np.array([1, 1, 1, 1], dtype=float) / 4),
]


def to_limbs(x):
    """x in [0,1) as an exact 80-bit fixed-point integer, two 40-bit limbs."""
    n = int(mp.nint(mp.mpf(x) * mp.mpf(2) ** 80)) % (1 << 80)
    return (n >> B) & MASK, n & MASK


def scale_limbs(hi, lo, c):
    """(hi,lo) * c mod 2^80 for a small integer c (|c| <= 64).  Vectorised, exact in int64."""
    neg = c < 0
    c = abs(c)
    lo2 = lo * c
    carry = lo2 >> B
    lo2 = lo2 & MASK
    hi2 = (hi * c + carry) & MASK
    if neg:
        # two's complement mod 2^80
        lo2 = (-lo2) & MASK
        borrow = np.where(lo2 != 0, 1, 0)
        hi2 = (-hi2 - borrow) & MASK
    return hi2, lo2


def frac(hi, lo):
    return hi.astype(np.float64) / SCALE + lo.astype(np.float64) / (float(SCALE) * SCALE)


def run(kind, params, K, checkpoints, cases=CASES):
    """kind = 'pair'  -> params = (alpha, beta) as mpmath reals; phases k*alpha, k*beta
       kind = 'sub'   -> params = (s, m, n); phases k*s*n and -k*s*m  (exactly on H)"""
    if kind == "pair":
        Ahi, Alo = to_limbs(params[0]); Bhi, Blo = to_limbs(params[1])
        gens = [(Ahi, Alo), (Bhi, Blo)]
        mults = None
    else:
        s, m, n = params
        Shi, Slo = to_limbs(s)
        gens = [(Shi, Slo)]
        mults = (n, -m)

    bases = [[0, 0] for _ in gens]
    tot = np.zeros(len(cases)); out = {c: [None] * len(cases) for c in checkpoints}
    done = 0
    while done < K:
        nn = min(CHUNK, K - done)
        j = np.arange(1, nn + 1, dtype=np.int64)
        phases = []
        for gi, (Ahi, Alo) in enumerate(gens):
            bh, bl = bases[gi]
            lo = bl + j * Alo
            carry = lo >> B
            lo = lo & MASK
            hi = (bh + j * Ahi + carry) & MASK
            phases.append((hi, lo))
            l2 = bl + nn * Alo
            c2 = l2 >> B
            bases[gi] = [(bh + nn * Ahi + c2) & MASK, l2 & MASK]
        if mults is None:
            fa = frac(*phases[0]); fb = frac(*phases[1])
        else:
            h1, l1 = scale_limbs(phases[0][0], phases[0][1], mults[0])
            h2, l2 = scale_limbs(phases[0][0], phases[0][1], mults[1])
            fa = frac(h1, l1); fb = frac(h2, l2)
        x = np.exp(2j * np.pi * fa); y = np.exp(2j * np.pi * fb); xy = x * y
        for ci, (lab, p) in enumerate(cases):
            a = np.abs(p[0] + p[1] * x + p[2] * y + p[3] * xy)
            cs = np.cumsum(np.log(np.maximum(a, 1e-323)))
            for cp in checkpoints:
                if done < cp <= done + nn:
                    out[cp][ci] = (tot[ci] + cs[cp - done - 1]) / cp
            tot[ci] += float(cs[-1])
        done += nn
    return out


# ------------------------------------------------------------------ targets, from R2_02
mP = {
    "B0b  (4,2,1,2)/9": mp.log(mp.mpf(4) / 9),
    "B4   (1,1,1,3)/6": -mp.log(2),
    "K1   (0,2,2,1)/5": mp.mpf('-0.75657358572374996479'),
    "SENSE C (1,1,1,1)/4": -mp.log(4),
}
mQ = {
    "B0b  (4,2,1,2)/9": mp.mpf('-0.810930204535048912'),
    "B4   (1,1,1,3)/6": mp.mpf('-0.693147181044366374'),
    "K1   (0,2,2,1)/5": mp.mpf('-0.756337009107433728'),
    "SENSE C (1,1,1,1)/4": -mp.log(4),          # EXACT: Q = (1+z^11)(1+z^20)/4
}

if __name__ == "__main__":
    K_BIG = 1_000_000_000
    K_MED = 200_000_000
    CPS_BIG = [10 ** 7, 10 ** 8, 3 * 10 ** 8, K_BIG]
    CPS_MED = [10 ** 7, 5 * 10 ** 7, K_MED]

    print("=" * 112)
    print("R2_03 — LONGER RULER.  Exact 80-bit phases, chunked modular accumulation, no K ceiling.")
    print("=" * 112)

    # ---------------------------------------------------------------- 1. the resonant arm
    print("\n" + "-" * 112)
    print("1.  THE RESONANT ARM f = 2.0, c = 1.1, BUILT ON ITS SUBTORUS BY CONSTRUCTION.")
    print("    s = -1/(20 pi); the orbit is (z^20, z^-11), z = e^{2 pi i k s}.  K up to 1e9.")
    print("    QUESTION: does the average leave m(P) for the subtorus value m(Q), as C-4 claims?")
    print("-" * 112)
    t0 = time.time()
    s_res = -1 / (20 * mp.pi)
    res = run("sub", (s_res % 1, 11, 20), K_BIG, CPS_BIG)
    print(f"    [{time.time()-t0:.1f}s]")
    for ci, (lab, p) in enumerate(CASES):
        print(f"\n    {lab}   m(P) = {mp.nstr(mP[lab],15)}   m(Q_(11,20)) = {mp.nstr(mQ[lab],15)}"
              f"   effect = {float(abs(mQ[lab]-mP[lab])):.3e}")
        for cp in CPS_BIG:
            v = res[cp][ci]
            d1 = abs(mp.mpf(v) - mP[lab]); d2 = abs(mp.mpf(v) - mQ[lab])
            near = "m(P)" if d1 < d2 else "m(Q)"
            if float(abs(mQ[lab] - mP[lab])) == 0.0:
                near = "TIE (the two targets are the SAME number)"
            print(f"       K={cp:>12,d}  avg = {v:.12f}   |avg-m(P)| = {float(d1):.3e}"
                  f"   |avg-m(Q)| = {float(d2):.3e}   nearer: {near}")

    # ---------------------------------------------------------------- 2. the Diophantine arm
    print("\n" + "-" * 112)
    print("2.  THE DIOPHANTINE ARM alpha = -2^(1/3), beta = 4^(1/3), at 80-bit exact phases.")
    print("    QUESTION: is C_04's 5.4e-08 at K = 1e7 finite-N discrepancy, or bias?")
    print("-" * 112)
    t0 = time.time()
    aD = (2 - mp.cbrt(2)) % 1
    bD = (mp.cbrt(4) - 1) % 1
    dio = run("pair", (aD, bD), K_MED, CPS_MED)
    print(f"    [{time.time()-t0:.1f}s]")
    for ci, (lab, p) in enumerate(CASES):
        row = "   ".join(f"K={cp:.0e}: {dio[cp][ci]:.10f} (err {float(abs(mp.mpf(dio[cp][ci])-mP[lab])):.1e})"
                         for cp in CPS_MED)
        print(f"    {lab:22s} {row}")

    # ------------------------------------------- 3. S4's OWN published verification connection
    print("\n" + "-" * 112)
    print("3.  f = 1.0, c = sqrt(2) -- S4's OWN PUBLISHED VERIFICATION CONNECTION [S4:603, S4:1054].")
    print("    L = {0}: m alpha + n beta in Z  =>  (-m + n sqrt2)/(2 pi) in Z  =>  an algebraic")
    print("    number equals 2 pi j, so j = 0, then n sqrt2 = m forces m = n = 0.  GENERIC.")
    print("    C-4 claims lambda = m(P) is false at EVERY connection the corpus has published.")
    print("-" * 112)
    t0 = time.time()
    aV = (-1 / (2 * mp.pi)) % 1
    bV = (mp.sqrt(2) / (2 * mp.pi)) % 1
    ver = run("pair", (aV, bV), K_MED, CPS_MED)
    print(f"    [{time.time()-t0:.1f}s]")
    for ci, (lab, p) in enumerate(CASES):
        row = "   ".join(f"K={cp:.0e}: {ver[cp][ci]:.10f} (err {float(abs(mp.mpf(ver[cp][ci])-mP[lab])):.1e})"
                         for cp in CPS_MED)
        print(f"    {lab:22s} {row}")

    # ------------------------------------------- 4. a second published generic connection
    print("\n" + "-" * 112)
    print("4.  f = 1.0, c = 3.0 -- a second published connection [S4:330].  Also generic.")
    print("-" * 112)
    t0 = time.time()
    aW = (-1 / (2 * mp.pi)) % 1
    bW = (mp.mpf(3) / (2 * mp.pi)) % 1
    w = run("pair", (aW, bW), 50_000_000, [10 ** 7, 50_000_000])
    print(f"    [{time.time()-t0:.1f}s]")
    for ci, (lab, p) in enumerate(CASES):
        row = "   ".join(f"K={cp:.0e}: {w[cp][ci]:.10f} (err {float(abs(mp.mpf(w[cp][ci])-mP[lab])):.1e})"
                         for cp in [10 ** 7, 50_000_000])
        print(f"    {lab:22s} {row}")

    print("\nDONE.")
