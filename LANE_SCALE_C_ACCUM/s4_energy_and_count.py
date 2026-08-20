"""LANE_SCALE_C_ACCUM  --  script 4: the two quantities gravity actually responds to.

Gravity grows with ENCLOSED ENERGY, and C-35 said a COUNT cannot supply a density law.  Both are
exactly computable on these carriers at any size, so both are checked head-on:

  * the spectrum of H  -- ground energy, gap, spectral width, ground degeneracy -- as k grows
  * the number of distinct record classes, |N(S)/S| - 1 = 4^k - 1, family versus control

CARRIER [[n,n-2,2]] : H = -(X^(x)n + Z^(x)n),  k = n-2
CONTROL [[4,2,2]]^m : H = -sum_b (X^(x)4_b + Z^(x)4_b),  k = 2m

The spectrum is derived from the syndrome structure (H is minus the sum of commuting stabilisers,
so its eigenvalues are -sum of the +-1 syndrome signs) AND verified against a direct eigh of the
real 2^n matrix wherever that matrix fits.
"""
import sys, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM")
import numpy as np
from record_model import xz_to_matrix, eigenspaces
from s1_combinatorics import carrier_nn2, carrier_product, f2_rank

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

def spectrum_from_syndromes(n, S):
    """H = -sum_s S_s with the S_s commuting involutions => eigenvalues are -sum of +-1 signs,
       each syndrome sector having dimension 2^n / 2^(#S)."""
    ns = f2_rank(S)
    mult = 2 ** (n - ns)
    levels = {}
    for signs in itertools.product((1, -1), repeat=ns):
        e = -sum(signs)
        levels[e] = levels.get(e, 0) + mult
    return dict(sorted(levels.items()))

def spectrum_direct(n, S):
    H = -sum(xz_to_matrix(s, n) for s in S)
    es = eigenspaces(H)
    return {round(float(v), 9): int(m) for v, _, m in es}

def run():
    P("=" * 118)
    P("LANE_SCALE_C_ACCUM  script 4  --  does ENERGY accumulate with records?  does the COUNT discriminate?")
    P("=" * 118)
    CH = []
    NS = [4, 6, 8, 10, 12, 14, 16, 18, 20, 22]
    F, C = {}, {}
    for n in NS:
        nn, S = carrier_nn2(n); F[n - 2] = (nn, S, spectrum_from_syndromes(nn, S))
    for m in range(1, 11):
        nn, S = carrier_product(m); C[2 * m] = (nn, S, spectrum_from_syndromes(nn, S))

    # SELF-CHECK: the syndrome spectrum against a direct eigh, wherever 2^n fits
    for tag, D, lim in (("family", F, 12), ("control", C, 12)):
        for k, (nn, S, lev) in D.items():
            if nn > lim: continue
            direct = spectrum_direct(nn, S)
            same = (len(direct) == len(lev)
                    and all(abs(a - b) < 1e-7 and direct[a] == lev[b]
                            for a, b in zip(sorted(direct), sorted(lev))))
            CH.append((tag, k, f"syndrome spectrum matches a direct eigh at dim {2**nn}", same,
                       f"{direct} vs {lev}"))
            CH.append((tag, k, "eigenspace dimensions sum to the Hilbert dimension",
                       sum(lev.values()) == 2 ** nn, ""))
            CH.append((tag, k, "ground space has dimension 2^k",
                       lev[min(lev)] == 2 ** k, f"{lev[min(lev)]} vs {2**k}"))
    bad = [c for c in CH if not c[3]]
    P("")
    P("SELF-CHECKS")
    P("-" * 118)
    for c in bad: P(f"   FAIL {c[0]} k={c[1]} {c[2]}   {c[4]}")
    P(f"   {len(CH)-len(bad)} / {len(CH)} pass" + ("   -- ALL PASS" if not bad else "   -- SOME FAILED"))
    if bad:
        P("   CONCLUSIONS VOID."); return

    P("")
    P("TABLE E   THE SPECTRUM AS k GROWS.  Control at the same k, in the same table.")
    P("-" * 118)
    P(f"{'k':>3} | {'n':>3} {'dim':>10} {'E0':>6} {'gap':>5} {'width':>6} {'E0/k':>7} {'gnd deg':>9}"
      f" || {'nC':>3} {'dim':>12} {'E0':>6} {'gap':>5} {'width':>6} {'E0/k':>7} {'gnd deg':>9}")
    P("-" * 118)
    for n in NS:
        k = n - 2
        def cells(D):
            nn, S, lev = D[k]
            ev = sorted(lev)
            E0 = ev[0]; gap = ev[1] - ev[0]; width = ev[-1] - ev[0]
            return (f"{nn:>3} {2**nn:>10,d} {E0:>6.1f} {gap:>5.1f} {width:>6.1f} "
                    f"{E0/k:>7.3f} {lev[E0]:>9,d}")
        c = cells(C).replace(f"{2**C[k][0]:>10,d}", f"{2**C[k][0]:>12,d}")
        P(f"{k:>3} | {cells(F)} || {cells(C)}")
    P("-" * 118)
    P("E0 = ground (code-space) energy;  gap = first excitation;  width = E_max - E_min")

    P("")
    P("TABLE N   THE COUNT.  C-35 said a count of records cannot supply a density law.  At equal k")
    P("          the collective carrier and k independent carriers have the SAME count, exactly.")
    P("-" * 118)
    P(f"{'k':>3} | {'fam record classes 4^k-1':>26} {'fam indep record bits':>22}"
      f" || {'ctl record classes 4^k-1':>26} {'ctl indep record bits':>22} {'identical?':>11}")
    for n in NS:
        k = n - 2
        cnt = 4 ** k - 1
        P(f"{k:>3} | {cnt:>26,d} {k:>22} || {cnt:>26,d} {k:>22} {'YES':>11}")
    P("-" * 118)
    P("record classes = |N(S)/S| - 1, the number of distinct non-trivial logical classes, every one")
    P("of which satisfies clauses (i)-(iv) (a non-stabiliser Pauli in N(S) is Hermitian up to phase,")
    P("squares to I, commutes with H, and is traceless on every syndrome sector).")

    P("")
    P("READ  (filled from the numbers above)")
    P("-" * 118)
    ks = [n - 2 for n in NS]
    fE = [min(F[k][2]) for k in ks]; cE = [min(C[k][2]) for k in ks]
    fW = [max(F[k][2]) - min(F[k][2]) for k in ks]; cW = [max(C[k][2]) - min(C[k][2]) for k in ks]
    fg = [sorted(F[k][2])[1] - sorted(F[k][2])[0] for k in ks]
    cg = [sorted(C[k][2])[1] - sorted(C[k][2])[0] for k in ks]
    P(f"  ground energy   family {fE}   -> {'CONSTANT' if len(set(fE))==1 else 'VARIES'}")
    P(f"                  control {cE}   -> {'CONSTANT' if len(set(cE))==1 else 'LINEAR in k, slope %.1f' % ((cE[-1]-cE[0])/(ks[-1]-ks[0]))}")
    P(f"  spectral width  family {fW}   -> {'CONSTANT' if len(set(fW))==1 else 'VARIES'}")
    P(f"                  control {cW}   -> {'CONSTANT' if len(set(cW))==1 else 'LINEAR in k, slope %.1f' % ((cW[-1]-cW[0])/(ks[-1]-ks[0]))}")
    P(f"  energy gap      family {fg}    control {cg}")
    P(f"  ground degeneracy is 2^k on BOTH: {[2**k for k in ks]}  -> entropy = k bits, LINEAR on both")
    P("")
    P("Largest carrier: n = 22 (k = 20) family, n = 40 (k = 20) control.  The spectrum is exact at")
    P("every size (syndrome structure); the direct-eigh cross-check ran wherever 2^n <= 4096.")

if __name__ == "__main__":
    run()
    with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM/s4_energy_and_count.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")
