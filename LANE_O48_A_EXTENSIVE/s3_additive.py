"""S3 -- ADDITIVITY OVER DISJOINT REGIONS, with a LIVE CONTROL that registers non-additivity.

Regions:  A = sites 1..m,  B = sites m+1..n,  disjoint and covering.
The only thing that can join them is the boundary bond J_m Z_m Z_{m+1}.

  DECOUPLED arm (the disjoint-regions test) : J_m := 0.  A and B share no term.
  COUPLED arm  (the LIVE CONTROL)           : J_m as drawn.  A and B share one term.

Both arms are measured by the SAME code path on the SAME whole-system enumeration, so a zero
defect in the decoupled arm cannot be the instrument being blind: the control arm reports a
non-zero defect from the identical routine.  All arithmetic is exact integer.

Quantities tested, on both arms:
  SPREAD            S = max_s E(s) - min_s E(s)                configuration-independent
  VARIANCE          Var_s E(s)                                 configuration-independent
  ENERGY of a named configuration  E(s)                        configuration-dependent
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, configs, energies_int

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 118)
p("S3  ADDITIVITY OVER DISJOINT REGIONS")
p("=" * 118)

def measure(n, a, sfix=None):
    """Exact spread, variance and E(sfix) by full 2^n enumeration. a has length n-1."""
    s = configs(n); E = energies_int(s, a)
    spread = int(E.max() - E.min())
    sq = sum(int(v) * int(v) for v in E)
    var = sq // (1 << n)                    # mean is exactly 0
    ev = None
    if sfix is not None:
        idx = 0
        for b in sfix: idx = (idx << 1) | (0 if b > 0 else 1)
        ev = int(E[idx])
    return spread, var, ev

p("")
p("-" * 118)
p("EXHAUSTIVE ENUMERATION OF BOTH ARMS.  n = m + k, region A = 1..m, region B = m+1..n.")
p("Q(A) and Q(B) are measured on the isolated sub-chains; Q(AuB) on the whole chain.")
p("DEFECT := Q(AuB) - Q(A) - Q(B), exact integer.  A zero defect means EXACT additivity.")
p("-" * 118)
p(f"{'n':>3} {'m':>3} {'k':>3} | {'DECOUPLED defect: spread':>25} {'variance':>26} {'E(s*)':>18}"
  f" | {'CONTROL defect: spread':>23} {'variance':>26} {'E(s*)':>18}")

rows = []
for (m, k) in [(2, 2), (3, 2), (3, 3), (4, 3), (4, 4), (5, 4), (5, 5), (6, 5), (6, 6),
               (7, 6), (7, 7), (8, 7), (8, 8), (9, 8), (10, 8)]:
    n = m + k
    a_full = couplings(n - 1)
    a_cut = list(a_full); a_cut[m - 1] = 0          # boundary bond removed -> disjoint regions
    aA, aB = a_full[:m - 1], a_full[m:]
    # a named configuration, the same one on the whole and on each part
    sfix = [1 if ((i * 5 + 3) % 7) < 4 else -1 for i in range(n)]
    for arm, aa in (("DEC", a_cut), ("CTL", a_full)):
        Sw, Vw, Ew = measure(n, aa, sfix)
        SA, VA, EA = measure(m, aA, sfix[:m])
        SB, VB, EB = measure(k, aB, sfix[m:])
        d = (Sw - SA - SB, Vw - VA - VB, Ew - EA - EB)
        if arm == "DEC": dd = d
        else: dc = d
    rows.append((n, m, k, dd, dc, a_full[m - 1]))
    p(f"{n:>3} {m:>3} {k:>3} | {dd[0]:>25} {dd[1]:>26} {dd[2]:>18}"
      f" | {dc[0]:>23} {dc[1]:>26} {dc[2]:>18}")

p("")
p("READ: every DECOUPLED defect is EXACTLY ZERO -- integer zero, no tolerance -- for the spread,")
p("      the variance and the named configuration's energy alike. Every CONTROL defect is")
p("      NON-ZERO in the same three columns from the same routine. The instrument registers")
p("      non-additivity when there is some, and reports none for disjoint regions.")
p("      ALL THREE QUANTITIES ARE EXACTLY ADDITIVE OVER DISJOINT REGIONS.")

p("")
p("-" * 118)
p("WHAT THE CONTROL DEFECT IS.  If the boundary defect were a VOLUME term it would grow with the")
p("region sizes; if a SURFACE term it would not. Exact values against the boundary coupling a_m.")
p("-" * 118)
p(f"{'n':>3} {'m':>3} {'k':>3} {'a_m (the one shared bond)':>26} {'CTL spread defect':>19} "
  f"{'/ 2a_m':>9} {'CTL var defect':>26} {'/ a_m^2':>9}")
for (n, m, k, dd, dc, am) in rows:
    p(f"{n:>3} {m:>3} {k:>3} {am:>26} {dc[0]:>19} {dc[0]/(2*am):>9.4f} {dc[1]:>26} "
      f"{dc[1]/(am*am):>9.4f}")
p("READ: the control defect is EXACTLY 2|J_m| for the spread and EXACTLY J_m^2 for the variance,")
p("      at every region size. IT DOES NOT GROW WITH THE REGIONS. The failure of additivity is a")
p("      PURE SURFACE TERM, set by the single bond crossing the cut and by nothing else.")
p("      This is the contact-or-nothing shape C-47 already recorded, appearing again here.")

p("")
p("-" * 118)
p("MANY DISJOINT REGIONS AT ONCE. Cut the chain into r equal blocks by deleting r-1 bonds, and")
p("compare the whole against the SUM over blocks. Exact, and the same control is run alongside.")
p("-" * 118)
p(f"{'n':>3} {'r':>3} {'DEC spread defect':>19} {'DEC var defect':>26} {'CTL spread defect':>19} "
  f"{'CTL var defect':>26}")
for n, r in [(8, 2), (8, 4), (12, 2), (12, 3), (12, 4), (12, 6), (16, 2), (16, 4), (16, 8)]:
    a_full = couplings(n - 1); w = n // r
    cuts = [j * w - 1 for j in range(1, r)]
    a_cut = list(a_full)
    for c in cuts: a_cut[c] = 0
    for arm, aa in (("DEC", a_cut), ("CTL", a_full)):
        Sw, Vw, _ = measure(n, aa)
        sS = sV = 0
        for j in range(r):
            seg = a_full[j * w:(j + 1) * w - 1]
            s_, v_, _ = measure(w, seg)
            sS += s_; sV += v_
        if arm == "DEC": d1 = (Sw - sS, Vw - sV)
        else: d2 = (Sw - sS, Vw - sV)
    p(f"{n:>3} {r:>3} {d1[0]:>19} {d1[1]:>26} {d2[0]:>19} {d2[1]:>26}")
p("READ: additivity over r disjoint blocks is exact for every r tested; the control is non-zero")
p("      throughout. Q(A_1 u ... u A_r) = sum_j Q(A_j) with no correction of any kind.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/s3_additive.txt", "w").write("\n".join(OUT) + "\n")
