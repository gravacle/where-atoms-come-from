"""V3 -- (A) is the passing quantity an ordinary spectral norm of H?
        (B) is its sign-definiteness UNFALSIFIABLE inside this family?
        (C) do the finding's reported S4 numbers match the lane's own output file?
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, uniform_couplings, configs, energies_int, dense_H

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 112)
p("V3-A  IS S JUST 2*||H||_op AND Var JUST Tr(H^2)/dim?  (both non-negative BY CONSTRUCTION)")
p("=" * 112)
p(f"{'n':>3} {'S (enum)/D':>13} {'2*||H||_op':>13} {'match':>7} {'Var/D^2':>13} "
  f"{'Tr(H^2)/2^n':>14} {'match':>7} {'||H||_op == sum|J| ?':>21}")
for n in (2, 3, 4, 5, 6, 7, 8):
    a = couplings(n - 1)
    s = configs(n); E = energies_int(s, a)
    S = int(E.max() - E.min())
    sq = sum(int(v) * int(v) for v in E)
    var = sq // (1 << n)
    H = dense_H(n, a)
    ev = np.linalg.eigvalsh(H)
    op = float(np.abs(ev).max())
    trh2 = float(np.trace(H @ H).real) / (1 << n)
    p(f"{n:>3} {S/D:>13.6f} {2*op:>13.6f} {str(abs(S/D - 2*op) < 1e-9):>7} {var/D**2:>13.6f} "
      f"{trh2:>14.6f} {str(abs(var/D**2 - trh2) < 1e-9):>7} {str(abs(op - sum(a)/D) < 1e-9):>21}")
p("READ: S is exactly twice the operator norm of H and Var is exactly Tr(H^2)/dim. Both are")
p("      standard spectral norms of the Hamiltonian, defined without reference to any record.")

p("")
p("=" * 112)
p("V3-B  CAN THE C-46 SIGN TEST ON S EVER FAIL IN THIS FAMILY?  A test that cannot fail is")
p("      not evidence. Try to break it: mixed-sign couplings, antiferromagnetic, alternating,")
p("      and a coupling set summing to zero.")
p("=" * 112)
p(f"{'coupling family':>34} {'sum_i J_i / D (can cancel)':>28} {'S/D':>12} {'S == 2*sum|J| ?':>16} "
  f"{'C-46 ratio on S terms':>22}")
rng = np.random.default_rng(5)
n = 8
base = couplings(n - 1)
fams = {
    "lane's own (all +)": base,
    "all antiferromagnetic (all -)": [-v for v in base],
    "alternating signs": [v if i % 2 == 0 else -v for i, v in enumerate(base)],
    "random signs": [v * int(rng.choice([-1, 1])) for v in base],
    "signs chosen so sum_i J_i = 0": None,
}
b = list(base)
b[0] = -sum(b[1:])          # forces sum to zero exactly
fams["signs chosen so sum_i J_i = 0"] = b
for name, a in fams.items():
    s = configs(n); E = energies_int(s, a)
    S = int(E.max() - E.min())
    terms = [2 * abs(v) for v in a]
    ratio = abs(sum(terms)) / sum(abs(t) for t in terms)
    p(f"{name:>34} {sum(a)/D:>28.6f} {S/D:>12.6f} {str(S == 2*sum(abs(v) for v in a)):>16} "
      f"{ratio:>22.6f}")
p("READ: the C-46 ratio on S is 1 for EVERY member of the family, including one whose couplings")
p("      sum exactly to zero. The terms are |J_i| because an absolute value was taken when the")
p("      spread was decomposed. No carrier in this family can make the test fail.")
p("      CONTRAST, in the same table's units: the CONFIGURATION-dependent bond terms J_i t_i")
p("      (no absolute value) do cancel -- that is the lane's own falling ratio.")
p(f"{'':>34} {'sum_i J_i t_i for a random s':>28}")
for k in range(3):
    a = base
    t = rng.choice([-1, 1], size=n - 1)
    tt = [int(a[i]) * int(t[i]) for i in range(n - 1)]
    p(f"{'random configuration ' + str(k):>34} {abs(sum(tt))/sum(abs(v) for v in tt):>28.6f}"
      f"   <- unsigned decomposition, ratio < 1")

p("")
p("=" * 112)
p("V3-C  TRANSCRIPTION CHECK: the finding's S4 table vs the lane's own s4_count.txt")
p("=" * 112)
def S_of(a): return 2 * sum(a)
for n in (8, 16):
    vals = [S_of(couplings(n - 1, stream=st)) / D for st in range(1, 401)]
    p(f"  n={n}: RECOMPUTED min {min(vals):.6f}  max {max(vals):.6f}  max/min {max(vals)/min(vals):.4f}")
p("  lane file s4_count.txt   n=8 : min 8.648892  max 12.857377  max/min 1.4866")
p("  lane file s4_count.txt   n=16: min 19.584455 max 25.179459  max/min 1.2857")
p("  FINDING's reported table n=8 : min 9.8688    max 12.0313    max/min 1.2191")
p("  FINDING's reported table n=16: min 22.0093   max 24.9861    max/min 1.1352")
p("READ: filled from the numbers above.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/VERIFY/v3_trivial_and_transcription.txt","w").write("\n".join(OUT)+"\n")
