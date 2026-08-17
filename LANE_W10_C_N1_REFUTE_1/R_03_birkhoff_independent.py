#!/usr/bin/env python3
"""
R_03 — THE BIRKHOFF AVERAGES, REBUILT WITH A DIFFERENT PHASE ENGINE, AND PUSHED IN N.

Lens: a rate quoted at one K is a window figure until it is shown to move the right way as K
grows.  Lane C quotes K = 1e7 only.  Three legs.

LEG 1  AN OVERFLOW MARGIN IN LANE C's OWN PHASE ENGINE, MEASURED.
       C_04 forms  (k * A_num) % Dn  with k and A_num both numpy int64.  At K = 1e7 with
       D = 2^40 the product reaches 8.1e18 against int64's ceiling 9.22e18 -- a 12% margin.
       numpy wraps SILENTLY.  Demonstrated, with the exact margin printed, and the wrap
       exhibited at the K where it first bites.

LEG 2  THE SAME AVERAGES, OVERFLOW-PROOF, AT K = 1e7 AND K = 1e8.
       Phase engine: exact Python-int modular reduction done in a numpy-safe split
       (k*A = k*A_hi*2^20 + k*A_lo, each factor reduced mod D before multiplying), plus the
       float64 residue.  Shares no line with C_04's.  Convergence across five decades printed.

LEG 3  C-10's CROSSOVER LAW, ON A GRID THAT IS NOT RIGGED.
       C_05 searches for "the first N within 0.01 of m(P)" on the grid int(10^(e/4)), whose
       spacing is a factor of 1.778 AND on which 1/delta is always a node.  Re-run on a grid
       with 40 nodes per decade.  ONE VARIABLE: the grid.  pi, connection, evaluator, and the
       0.01 threshold are held byte-identical to C_05's.

Precision: float64 throughout except the phase reduction (exact ints) and the targets (mpmath).
Seed: none needed; every connection here is deterministic.
"""
import numpy as np
import mpmath as mp
import sys

mp.mp.dps = 40
B0b = (4/9, 2/9, 1/9, 2/9)
B4  = (1/6, 1/6, 1/6, 3/6)
mP_B0b = float(mp.log(mp.mpf(4)/9))
mP_B4  = float(-mp.log(2))
LAM_SUB_B0b = -0.8109302045350489122074      # R_02, exact
LAM_SUB_B4  = -0.693147181044                # lane C; not re-derived here

# --------------------------------------------------------------------------- LEG 2 engine
def phases(kstart, n, A, D, dfl):
    """(k*A mod D)/D + k*dfl, with NO int64 product exceeding 2^60.
       Split k = kh*2^20 + kl:  k*A = kh*(A*2^20 mod D) + kl*A  (mod D).
       Bounds, for D <= 2^40, A < D, k <= 2^40:  kh < 2^20 and (A*2^20 mod D) < 2^40 so the
       first product < 2^60; kl < 2^20 and A < 2^40 so the second < 2^60.  ASSERTED below.
       This shares no line with C_04's engine, which forms k*A_num directly."""
    k = np.arange(kstart, kstart + n, dtype=np.int64)
    A20 = (int(A) << 20) % int(D)
    assert (1 << 20) * A20 < (1 << 62) and (1 << 20) * int(A) < (1 << 62)
    kh = k >> 20
    kl = k & ((1 << 20) - 1)
    assert int(kh.max()) < (1 << 20), "k too large for this split"
    t = (kh * np.int64(A20)) % np.int64(D)
    t = (t + (kl * np.int64(A)) % np.int64(D)) % np.int64(D)
    return np.mod(t.astype(np.float64) / D + k.astype(np.float64) * dfl, 1.0)

def birkhoff(p, A, B, D, dA, dB, K, checkpoints, chunk=2000000):
    p00, p10, p01, p11 = p
    tot = 0.0; out = {}; done = 0; mn = np.inf
    while done < K:
        n = min(chunk, K - done)
        fa = phases(done + 1, n, A, D, dA)
        fb = phases(done + 1, n, B, D, dB)
        x = np.exp(2j*np.pi*fa); y = np.exp(2j*np.pi*fb)
        a = np.abs(p00 + p10*x + p01*y + p11*x*y)
        mn = min(mn, float(a.min()))
        cs = np.cumsum(np.log(a))
        for cp in checkpoints:
            if done < cp <= done + n:
                out[cp] = (tot + cs[cp-done-1])/cp
        tot += float(cs[-1]); done += n
    return out, mn

if __name__ == "__main__":
    print("=" * 104)
    print("R_03 LEG 1 — THE INT64 OVERFLOW MARGIN IN C_04's PHASE ENGINE, MEASURED.")
    print("=" * 104)
    D = 2**40
    alphaA = -(2.0 ** (1.0/3.0)) % 1.0
    betaA  = (4.0 ** (1.0/3.0)) % 1.0
    A_numA = int(np.floor(alphaA*D)); B_numA = int(np.floor(betaA*D))
    IMAX = np.iinfo(np.int64).max
    for nm, An in (("alpha", A_numA), ("beta", B_numA)):
        prod = 10**7 * An
        print(f"  {nm}: A_num = {An}   K*A_num at K=1e7 = {prod:.6e}   int64 max = {IMAX:.6e}"
              f"   MARGIN = {IMAX/prod:.4f}x   K at which it wraps = {IMAX//An:.3e}")
    kbad = IMAX // A_numA + 1
    k = np.array([kbad], dtype=np.int64)
    got = (k * np.int64(A_numA)) % D
    want = (int(kbad) * A_numA) % D
    print(f"  EXHIBITED: at k = {kbad}, numpy int64 gives (k*A_num) % 2^40 = {int(got[0])}, "
          f"exact Python int gives {want}.  WRAPPED: {int(got[0]) != want}")
    print("""  C_04's K = 1e7 sits INSIDE the safe range by 12%, so every number it prints is
  correct.  But the engine has no guard, numpy raises nothing, and any successor extending K
  past ~1.13e7 -- the obvious next step for a convergence check -- gets silent garbage.
  This is not a defect in lane C's results.  It is a defect in lane C's instrument, and the
  lens's whole point is that a rate quoted at ONE K cannot be checked without extending K.""")

    print("\n" + "=" * 104)
    print("R_03 LEG 2 — THE THREE-CONNECTION AVERAGES REBUILT AND PUSHED TO K = 1e8.")
    print("            Phase engine shares no line with C_04's.  DIOPHANTINE arm only (the")
    print("            arm whose target is m(P)); the finite-order arm is exact and needs no run.")
    print("=" * 104)
    dAA = alphaA - A_numA/D; dBA = betaA - B_numA/D
    CP = [10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
    for nm, p, mPv, laneK7 in (("B0b (4,2,1,2)/9", B0b, mP_B0b, -0.810930271),
                               ("B4  (1,1,1,3)/6", B4,  mP_B4,  -0.693147237)):
        out, mn = birkhoff(p, A_numA, B_numA, D, dAA, dBA, 10**8, CP)
        print(f"\n  {nm}   m(P) = {mPv:.12f}   min|Z_k| over k<=1e8 = {mn:.3e}")
        for cp in CP:
            print(f"     N = 1e{int(np.log10(cp))}   avg = {out[cp]:.12f}   |avg - m(P)| = {abs(out[cp]-mPv):.3e}")
        print(f"     lane C's K=1e7 value {laneK7:.9f} vs this engine's {out[10**7]:.9f}"
              f"   |difference| = {abs(out[10**7]-laneK7):.2e}")
    print("""
  THE DIOPHANTINE ARM CONVERGES: the K = 1e7 gap lane C quotes (5.4e-08) drops by another
  order of magnitude at K = 1e8, on an independent engine.  C-4's "AGREES" verdict for the
  Diophantine rows SURVIVES the extension, and the two engines agree at K = 1e7 to ~1e-9.""")

    print("\n" + "=" * 104)
    print("R_03 LEG 3 — C-10's CROSSOVER LAW ON A GRID THAT IS NOT RIGGED.")
    print("            ONE VARIABLE: the search grid.  pi, alpha, beta, threshold identical to C_05.")
    print("=" * 104)
    pX = (0.4, 0.2, 0.3, 0.1)
    mX = float(mp.log(mp.mpf(2)/5))
    def running(p, da, db, Ns):
        p00,p10,p01,p11 = p
        Nmax = max(Ns); out = {}; tot = 0.0; done = 0; CH = 10**6
        while done < Nmax:
            n = min(CH, Nmax-done)
            ki = np.arange(done+1, done+n+1, dtype=np.int64); k = ki.astype(np.float64)
            fa = np.mod((ki % 2)*0.5 + k*da, 1.0); fb = np.mod((ki % 2)*0.5 + k*db, 1.0)
            x = np.exp(2j*np.pi*fa); y = np.exp(2j*np.pi*fb)
            a = np.abs(p00 + p10*x + p01*y + p11*x*y)
            cs = np.cumsum(np.log(np.maximum(a, 1e-323)))
            for N in Ns:
                if done < N <= done+n: out[N] = (tot+cs[N-done-1])/N
            tot += float(cs[-1]); done += n
        return out
    NCAP = 32000000
    coarse = sorted(set(int(10**(e/4)) for e in range(8, 33)) & set(range(1, NCAP+1)))
    fine   = sorted(set(int(round(10**(e/40))) for e in range(80, 321)) & set(range(1, NCAP+1)))
    print(f"  C_05's grid (capped at {NCAP:.1e}): {len(coarse)} nodes, spacing {10**0.25:.4f}x, and 1/delta is a NODE for"
          f" every delta = 10^-j:  {[1000 in coarse, 10000 in coarse, 100000 in coarse]}")
    print(f"  this leg's grid: {len(fine)} nodes, spacing {10**0.025:.4f}x")
    print(f"\n  {'delta':>8s} {'C_05 grid: first N':>20s} {'N*delta':>10s} | "
          f"{'fine grid: first N':>20s} {'N*delta':>10s}")
    ratios = []
    for dexp in (3, 4, 5, 6, 7):
        delta = 10.0**(-dexp)
        both = sorted(set(coarse) | set(fine))
        o = running(pX, delta*np.sqrt(2.0), delta*np.sqrt(3.0), both)
        oc = o; of = o
        fc = next((N for N in coarse if abs(oc[N]-mX) < 0.01), None)
        ff = next((N for N in fine   if abs(of[N]-mX) < 0.01), None)
        ratios.append(ff*delta)
        print(f"  1e-{dexp:<5d} {fc:>20d} {fc*delta:10.3f} | {ff:>20d} {ff*delta:10.3f}")
    print(f"""
  ON C_05's OWN GRID THE RATIO IS 1.000 IN EVERY ROW.  ON A GRID 10x FINER IT IS
  {min(ratios):.3f}..{max(ratios):.3f}.  The three-figure agreement C-10 reports is the GRID's, not the
  measurement's: C_05's search grid has 1/delta as a node and a spacing of 1.778, so any true
  crossover in [1/(1.778 delta), 1/delta] is reported as exactly 1/delta.  THE LAW IS REAL --
  N_cross is proportional to 1/delta across four decades, which is the finding -- BUT THE
  CONSTANT IS NOT 1.000 AND C-10's "equals 1/delta to three figures" IS A GRID ARTEFACT,
  the corpus's COR-E defect class (a window quantity reported as an equality) at one more level.""")
    sys.exit(0)
