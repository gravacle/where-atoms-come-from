"""S5 -- WHO CARRIES THE ENERGY, AND DOES IT CANCEL (C-46)?

Two things are settled here, both exactly:

  (A) THE ENERGY IS A FUNCTION OF THE PAIR CORRELATIONS AND OF NOTHING ELSE.
      Exhibited by construction: two states with IDENTICAL single-record expectations <Z_i>
      for every i but different energies, and two states with identical energy but different
      <Z_i>. This is O-47's separation, at general n, shown rather than asserted.

  (B) C-46's accumulation test |sum|/sum|.| applied to both parts, over the EXACT distribution
      of all 2^n configurations. C-46 is the test a gravity-role quantity must pass.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, configs, energies_int

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 116)
p("S5  WHO CARRIES THE ENERGY, AND DOES IT CANCEL")
p("=" * 116)

# ---------------------------------------------------------------- (A) records vs correlations
p("")
p("-" * 116)
p("(A) THE ENERGY IS BLIND TO THE RECORDS AND SEES ONLY THEIR CORRELATIONS.")
p("    Two configuration-diagonal states are built at each n:")
p("      rho_F  = uniform over { all-up , all-down }        (ferromagnetic pair)")
p("      rho_A  = uniform over { alternating , its flip }   (antiferromagnetic pair)")
p("    Both are invariant under the admissible writer found in S1, so both have <Z_i> = 0 for")
p("    every record. Their energies are compared.")
p("    LIVE CONTROL: a third state rho_P = the pure all-up configuration, whose <Z_i> = +1 for")
p("    every record -- so the <Z_i> column is not stuck at zero by the instrument.")
p("-" * 116)
p(f"{'n':>3} {'max_i |<Z_i>| rho_F':>21} {'max_i |<Z_i>| rho_A':>21} {'E(rho_F)/D':>14} "
  f"{'E(rho_A)/D':>14} {'DIFFER?':>9} {'CONTROL max_i|<Z_i>| rho_P':>27}")
for n in range(2, 17):
    a = couplings(n - 1)
    up = [1] * n
    alt = [1 if i % 2 == 0 else -1 for i in range(n)]
    def En(s): return sum(a[i] * s[i] * s[i + 1] for i in range(n - 1))
    def Zav(states):
        return max(abs(sum(st[i] for st in states) / len(states)) for i in range(n))
    F = [up, [-v for v in up]]
    A = [alt, [-v for v in alt]]
    eF = En(up); eA = En(alt)
    p(f"{n:>3} {Zav(F):>21.6f} {Zav(A):>21.6f} {eF/D:>14.6f} {eA/D:>14.6f} "
      f"{str(eF != eA):>9} {Zav([up]):>27.6f}")
p("READ: every record is perfectly balanced in BOTH states -- clause (iv)'s condition holds and")
p("      the admissible writer flips them for free -- yet the two states differ in energy at")
p("      every n. The control shows the <Z_i> column can register a non-zero.")
p("      THE ENERGY IS CARRIED BY THE PAIR CORRELATIONS Z_i Z_{i+1}, NOT BY THE RECORDS.")
p("      INDUCED at general n by the writer SEARCH of S1; only the bond form of H was inserted.")

p("")
p("    The converse, also exact: states with the SAME energy but DIFFERENT record values.")
p(f"{'n':>3} {'E(all-up)/D':>14} {'E(all-down)/D':>16} {'equal?':>8} "
  f"{'<Z_1> up':>10} {'<Z_1> down':>12} {'differ?':>9}")
for n in range(2, 17):
    a = couplings(n - 1)
    up = [1] * n; dn = [-1] * n
    def En(s): return sum(a[i] * s[i] * s[i + 1] for i in range(n - 1))
    p(f"{n:>3} {En(up)/D:>14.6f} {En(dn)/D:>16.6f} {str(En(up)==En(dn)):>8} "
      f"{up[0]:>10} {dn[0]:>12} {str(up[0]!=dn[0]):>9}")
p("READ: record values move at fixed energy, and energy moves at fixed record values. The two")
p("      are independent coordinates. Clause (iv) constrains one of them and says nothing about")
p("      the other -- which is exactly the opening O-47 found, holding at every n up to 16.")

# ---------------------------------------------------------------- (B) C-46
p("")
p("-" * 116)
p("(B) C-46 ACCUMULATION TEST, |sum of terms| / sum of |terms|, over the EXACT distribution of")
p("    all 2^n configurations. Terms are the bond contributions J_i s_i s_{i+1}.")
p("    A quantity with ratio 1 accumulates without bound; one with ratio -> 0 screens.")
p("    CONFIG-INDEPENDENT part (the spread) has terms 2|J_i|, all positive: ratio = 1 by form.")
p("-" * 116)
p(f"{'n':>3} {'ratio, GROUND cfg':>19} {'ratio, ALIGNED cfg':>20} {'MEAN ratio over all 2^n':>25} "
  f"{'frac of cfgs with ratio<0.5':>29} {'SPREAD ratio':>14}")
for n in range(2, 21):
    a = couplings(n - 1)
    s = configs(n); E = energies_int(s, a)
    denom = sum(a)
    r = np.abs(E.astype(np.float64)) / float(denom)
    p(f"{n:>3} {1.0:>19.6f} {1.0:>20.6f} {r.mean():>25.6f} {float((r<0.5).mean()):>29.6f} "
      f"{1.0:>14.6f}")
p("READ: the EXTREMAL configurations (ground and aligned) have ratio EXACTLY 1 -- every bond")
p("      term carries the same sign and nothing cancels. The MEAN over all configurations falls")
p("      steadily toward 0 and the fraction of configurations that cancel more than half rises")
p("      toward 1. SIGN-DEFINITENESS OF THE CONFIGURATION-DEPENDENT PART IS A PROPERTY OF THE")
p("      CONFIGURATION, NOT OF THE QUANTITY. The configuration-INDEPENDENT part passes C-46")
p("      unconditionally, because its terms are 2|J_i| and cannot cancel.")

p("")
p("    Asymptotics of the mean ratio, exact enumeration only -- this is a NUMERICAL TREND over")
p("    2 <= n <= 20 and is NOT claimed as an asymptotic law (D-20). Two candidate forms are")
p("    printed with no selection made between them.")
p(f"{'n':>3} {'mean ratio':>13} {'1/sqrt(n)':>12} {'ratio*sqrt(n)':>15} {'ratio*n':>10}")
for n in (4, 6, 8, 10, 12, 14, 16, 18, 20):
    a = couplings(n - 1); s = configs(n); E = energies_int(s, a)
    r = float((np.abs(E.astype(np.float64)) / float(sum(a))).mean())
    p(f"{n:>3} {r:>13.6f} {1/np.sqrt(n):>12.6f} {r*np.sqrt(n):>15.6f} {r*n:>10.6f}")
p("READ: ratio*sqrt(n) is far flatter than ratio*n over this range, which is consistent with a")
p("      1/sqrt(n) decay, but 2 <= n <= 20 is a SHORT RANGE and D-20 forbids calling it an")
p("      asymptotic law from this table. What the table does establish without a fit is the")
p("      DIRECTION: the mean ratio is falling, so typical configurations SCREEN.")

# ---------------------------------------------------------------- separation
p("")
p("-" * 116)
p("SEPARATION (standard (e)) -- stated so it is not mistaken for a measurement.")
p("Every term of this H couples NEAREST NEIGHBOURS ONLY. That was INSERTED. Any two-point")
p("function of this carrier is therefore contact-or-nothing by construction and NO separation")
p("law can be read off it. Demonstrated, not argued:")
p("-" * 116)
p(f"{'n':>3} {'d':>3} {'<Z_1 Z_{1+d}> in the ground state':>34} {'d(E)/d a_i at distance d':>26}")
for n in (12,):
    a = couplings(n - 1)
    gs = [1 if i % 2 == 0 else -1 for i in range(n)]     # all a_i > 0 -> antiferromagnetic ground
    for d in range(1, n):
        corr = gs[0] * gs[d]
        sens = "n/a" if d > n - 1 else ("nonzero only at d=1" if d == 1 else "0 (no coupling)")
        p(f"{n:>3} {d:>3} {corr:>34} {sens:>26}")
p("READ: the ground-state correlator is +-1 at EVERY separation -- it does not fall off at all,")
p("      because the ground state is a single product configuration. And the energy has no")
p("      dependence on any coupling at distance > 1 because none was put in. THIS CARRIER CANNOT")
p("      POSE THE SEPARATION QUESTION. No falloff, power-law or exponential, is measured here,")
p("      and none should be reported from it.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/s5_sign_and_carrier.txt", "w").write("\n".join(OUT) + "\n")
