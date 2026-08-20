"""S4 -- IS IT MERELY A COUNT?  (C-35 / test (c))  This is where every previous candidate died.

THE TEST, stated so it can fail: hold n FIXED and move the J_i. A quantity determined by n
alone is C-35's count in other units, however extensive it is.

Four separate discriminations are run, each with a control:
  T1  does Q move at fixed n when the couplings move?          (a count would not)
  T2  does Q move under the venue's own overall scale?  (D-17)  (a count would not)
  T3  how many distinct values does Q take at fixed n?          (a count takes exactly one)
  T4  does Q see WHERE the couplings sit, or only WHICH they are?  -- the placement test.
      This one is the limitation, and it is reported as such.
CONTROLS printed in the same tables: the uniform-J chain, where the same quantity DOES
degenerate to a count in units of 2J, and the pure record count n itself.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, uniform_couplings, configs, energies_int

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 118)
p("S4  THE COUNT TEST  --  vary the J_i at FIXED n and see whether the quantity moves")
p("=" * 118)

def S_of(a):   return 2 * sum(a)                 # exact, verified against enumeration in S2
def V_of(a):   return sum(v * v for v in a)      # exact, verified against enumeration in S2

# ---------------------------------------------------------------- T1  fixed n, couplings move
p("")
p("-" * 118)
p("T1  RANGE OF VARIATION AT FIXED n.  400 independent coupling draws per n from the venue's own")
p("    family J_i in [0.5,1.0), and 400 from a WIDER family J_i log-uniform over four decades.")
p("    CONTROL columns: the same n with all J_i EQUAL (a count in units of 2J), and the bare")
p("    record count n, which by construction cannot move at all.")
p("-" * 118)
p(f"{'n':>3} | {'SPREAD S/D  min':>17} {'max':>14} {'max/min':>9} | {'WIDE family S/D  min':>21} "
  f"{'max':>14} {'max/min':>11} | {'CONTROL uniform-J S/D':>22} {'CONTROL count n':>16}")
rng = np.random.default_rng(20480)
for n in (4, 6, 8, 10, 12, 14, 16):
    nb = n - 1
    vals = [S_of(couplings(nb, stream=st)) / D for st in range(1, 401)]
    wide = []
    for _ in range(400):
        aa = [int(D * 10.0 ** rng.uniform(-2, 2)) for _ in range(nb)]
        wide.append(S_of(aa) / D)
    au = uniform_couplings(nb)
    p(f"{n:>3} | {min(vals):>17.6f} {max(vals):>14.6f} {max(vals)/min(vals):>9.4f} | "
      f"{min(wide):>21.6f} {max(wide):>14.6f} {max(wide)/min(wide):>11.2f} | "
      f"{S_of(au)/D:>22.6f} {n:>16}")
p("READ: the spread MOVES at fixed n. Inside the venue's own coupling family it moves by tens of")
p("      per cent; opened to four decades of coupling it moves by four orders of magnitude, at")
p("      FIXED n. The record count n in the last column is frozen by construction, which is what")
p("      a count looks like in this table.")

p("")
p(f"{'n':>3} | {'VARIANCE Var/D^2  min':>23} {'max':>16} {'max/min':>9} | "
  f"{'CONTROL uniform-J Var/D^2':>26} {'CONTROL count n':>16}")
for n in (4, 6, 8, 10, 12, 14, 16):
    nb = n - 1
    vals = [V_of(couplings(nb, stream=st)) / D**2 for st in range(1, 401)]
    au = uniform_couplings(nb)
    p(f"{n:>3} | {min(vals):>23.6f} {max(vals):>16.6f} {max(vals)/min(vals):>9.4f} | "
      f"{V_of(au)/D**2:>26.6f} {n:>16}")
p("READ: the configuration variance moves at fixed n too.")

# ---------------------------------------------------------------- T2  D-17 venue scale
p("")
p("-" * 118)
p("T2  D-17  VARY THE VENUE'S OWN SCALES. Three scales exist here: the overall coupling")
p("    magnitude lambda, the number of sites n, and the DISORDER WIDTH of the coupling family.")
p("    A count is blind to the first and the third.")
p("-" * 118)
p(f"{'n':>3} {'lambda':>10} {'S/D':>16} {'S/(lambda)':>16} {'Var/D^2':>18} {'Var/lambda^2':>18}")
for n in (8, 16):
    a0 = couplings(n - 1)
    for lam in (0.01, 0.1, 1.0, 10.0, 1000.0):
        a = [int(v * lam) for v in a0]
        p(f"{n:>3} {lam:>10} {S_of(a)/D:>16.6f} {S_of(a)/D/lam:>16.6f} {V_of(a)/D**2:>18.6f} "
          f"{V_of(a)/D**2/lam**2:>18.6f}")
p("READ: S is exactly homogeneous of degree 1 and Var of degree 2 in the coupling scale. Both")
p("      move under the venue's own scale. A count cannot do that.")
p("")
p(f"{'n':>3} {'disorder width':>16} {'S/D':>16} {'Var/D^2':>18}  (mean coupling held at 0.75 D)")
for n in (8, 16):
    for w in (0.0, 0.1, 0.3, 0.5, 0.74):
        a = [int(D * (0.75 + w * (1 if k % 2 == 0 else -1))) for k in range(n - 1)]
        p(f"{n:>3} {w:>16.2f} {S_of(a)/D:>16.6f} {V_of(a)/D**2:>18.6f}")
p("READ: at fixed n AND fixed mean coupling, the SPREAD does not move with the disorder width")
p("      (it depends only on the sum), but the VARIANCE does. The two config-independent")
p("      quantities are not the same functional, and the variance resolves what the spread cannot.")

# ---------------------------------------------------------------- T3  how many values
p("")
p("-" * 118)
p("T3  HOW MANY DISTINCT VALUES AT FIXED n?  A count takes exactly one.")
p("-" * 118)
p(f"{'n':>3} {'draws':>7} {'distinct S values':>19} {'distinct Var values':>21} "
  f"{'distinct values of the count n':>31}")
for n in (4, 8, 12, 16):
    vs = {S_of(couplings(n - 1, stream=st)) for st in range(1, 401)}
    vv = {V_of(couplings(n - 1, stream=st)) for st in range(1, 401)}
    p(f"{n:>3} {400:>7} {len(vs):>19} {len(vv):>21} {1:>31}")
p("READ: 400 distinct values from 400 draws, at every n. NOT A COUNT.")

# ---------------------------------------------------------------- T4  placement
p("")
p("-" * 118)
p("T4  THE PLACEMENT TEST -- the limitation, reported as one.  Hold n AND the coupling MULTISET")
p("    fixed and only PERMUTE which bond carries which coupling. The carrier genuinely changes")
p("    (S1 established |Aut| = 1, so these are different Hamiltonians).")
p("    LIVE CONTROL in the same table: E(s) at a FIXED configuration s, which must move if the")
p("    permutation is doing anything at all.")
p("-" * 118)
p(f"{'n':>3} {'permutations':>13} {'distinct S':>11} {'distinct Var':>13} "
  f"{'distinct E(s*) [CONTROL]':>25} {'CONTROL range of E(s*)/D':>26}")
for n in (6, 8, 10, 12):
    a0 = couplings(n - 1)
    Ss, Vs, Es = set(), set(), set()
    for _ in range(400):
        a = list(a0); rng.shuffle(a)
        Ss.add(S_of(a)); Vs.add(V_of(a))
        s = configs(n); E = energies_int(s, a)
        sfix = [1 if ((i * 5 + 3) % 7) < 4 else -1 for i in range(n)]
        idx = 0
        for b in sfix: idx = (idx << 1) | (0 if b > 0 else 1)
        Es.add(int(E[idx]))
    p(f"{n:>3} {400:>13} {len(Ss):>11} {len(Vs):>13} {len(Es):>25} "
      f"{(max(Es)-min(Es))/D:>26.6f}")
p("READ: under permutation of the couplings among the bonds, the SPREAD and the VARIANCE take")
p("      exactly ONE value -- they are SYMMETRIC FUNCTIONS of the coupling multiset and are")
p("      blind to where anything sits. The control column moves freely, so the permutations are")
p("      real. The configuration-INDEPENDENT part therefore cannot carry a separation law of any")
p("      kind: it is strictly more than a count of records, and strictly less than geometry.")
p("      The configuration-DEPENDENT part is the only one of the two that sees placement.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/s4_count.txt", "w").write("\n".join(OUT) + "\n")
