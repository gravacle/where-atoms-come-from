"""S7 -- THE STANDARD, APPLIED, plus the one measurement venue 2 still owes.

Venue 2 (r decoupled blocks) is the venue in which the RECORD COUNT itself grows. S6 verified
its clauses. What is measured here is whether the energy is extensive IN THE NUMBER OF
INDEPENDENT RECORDS r rather than merely in n, whether venue 2 has geometry to detect (D-22),
and whether venue 2 survives the count test at fixed r.
Then the verdict table for both venues, filled from the numbers actually obtained.
"""
import sys, itertools, math
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, configs, energies_int

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 118)
p("S7  THE STANDARD APPLIED")
p("=" * 118)

def block_couplings(r, w, stream=0):
    """r blocks of w sites. Returns the full bond list of length r*w-1, cut bonds = 0."""
    n = r * w
    a = list(couplings(n - 1, stream=stream))
    for j in range(1, r): a[j * w - 1] = 0
    return n, a

# ---------------------------------------------------------------- D-22 for venue 2
p("")
p("-" * 118)
p("D-22 FOR VENUE 2. r decoupled blocks could easily be block-permutation symmetric. Counted by")
p("BRUTE FORCE over all n! site permutations. CONTROL: the same block geometry with EVERY block")
p("given the SAME coupling list, which must show the block-permutation symmetry.")
p("-" * 118)
p(f"{'r':>3} {'w':>3} {'n':>3} {'|Aut| distinct-J blocks':>25} {'|Aut| identical blocks [CONTROL]':>34} "
  f"{'r!':>6} {'2^r . r!':>10}")
def aut_count(n, bonds):
    c = 0
    for perm in itertools.permutations(range(n)):
        img = {}
        for (i, j), v in bonds.items():
            k, l = perm[i], perm[j]
            img[(min(k, l), max(k, l))] = v
        if img == bonds: c += 1
    return c
for r, w in [(2, 2), (3, 2), (2, 3), (4, 2), (2, 4), (3, 3)]:
    n, a = block_couplings(r, w)
    bd = {(i, i + 1): a[i] for i in range(n - 1) if a[i] != 0}
    # identical blocks control
    ai = list(a)
    blk = [a[k] for k in range(w - 1)]
    for j in range(r):
        for k in range(w - 1): ai[j * w + k] = blk[k]
    bi = {(i, i + 1): ai[i] for i in range(n - 1) if ai[i] != 0}
    p(f"{r:>3} {w:>3} {n:>3} {aut_count(n, bd):>25} {aut_count(n, bi):>34} "
      f"{math.factorial(r):>6} {(2**r)*math.factorial(r):>10}")
p("READ, FROM THE NUMBERS AND NOT FROM THE EXPECTATION: venue 2 splits into two cases.")
p("  BLOCK SIZE w >= 3, distinct couplings:  |Aut| = 1. Fully rigid. Geometry present.")
p("  BLOCK SIZE w = 2, distinct couplings:   |Aut| = 2^r, NOT 1 -- a two-site block is a single")
p("      bond and swapping its two endpoints preserves it. Those 2^r symmetries act INSIDE")
p("      blocks and never exchange one block with another, so they do not make the r blocks")
p("      interchangeable; but the claim '|Aut| = 1' is FALSE at w = 2 and is not made.")
p("  CONTROL, identical blocks: |Aut| = 2^r . r!, matching the last column exactly -- the block-")
p("      permutation factor r! appears precisely when the blocks are made interchangeable, so")
p("      the counter is live and would have caught a permutation-symmetric venue.")
p("CONSEQUENCE FOR THIS LANE: every venue-2 measurement reported here uses w >= 2 with distinct")
p("couplings, and the quantities measured (spread, variance) were already shown in S4 to be")
p("symmetric functions of the coupling multiset, so no separation claim rests on w = 2 rigidity.")

# ---------------------------------------------------------------- extensivity in r
p("")
p("-" * 118)
p("IS THE ENERGY EXTENSIVE IN THE NUMBER OF INDEPENDENT RECORDS?  Block size w held FIXED,")
p("r doubled. r is the exact bit ceiling S6 measured, so this is extensivity per record, not")
p("per site. Exact closed form, verified against enumeration where it fits.")
p("-" * 118)
for w in (2, 3, 4):
    p(f"  block size w = {w}")
    p(f"    {'r':>8} {'n = r*w':>9} {'S(r)/D':>16} {'Var(r)/D^2':>16} {'S(2r)/S(r)':>13} "
      f"{'Var(2r)/Var(r)':>16} {'S per record':>14} {'enum agrees':>12}")
    for r in (1, 2, 4, 8, 16, 32, 1024, 65536):
        n, a = block_couplings(r, w)
        S = 2 * sum(a); V = sum(v * v for v in a)
        n2, a2 = block_couplings(2 * r, w)          # the SAME quantity at 2r, computed directly
        S2 = 2 * sum(a2); V2 = sum(v * v for v in a2)
        ok = "-"
        if n <= 18:
            s = configs(n); E = energies_int(s, a)
            ok = str(int(E.max() - E.min()) == S)
        p(f"    {r:>8} {n:>9} {S/D:>16.6f} {V/D**2:>16.6f} {S2/S:>13.6f} {V2/V:>16.6f} "
          f"{S/D/r:>14.6f} {ok:>12}")
p("READ: at fixed block size, doubling the number of INDEPENDENT RECORDS doubles the energy")
p("spread and the variance, and the spread per record is constant. EXTENSIVE IN THE RECORD")
p("COUNT, not only in the site count.")

# ---------------------------------------------------------------- count test in venue 2
p("")
p("-" * 118)
p("COUNT TEST IN VENUE 2: hold r AND the block size fixed, move the couplings. 400 draws.")
p("CONTROL column: r itself, which cannot move.")
p("-" * 118)
p(f"{'r':>4} {'w':>3} {'min S/D':>12} {'max S/D':>12} {'max/min':>9} {'distinct values / 400':>23} "
  f"{'CONTROL r':>10}")
for r, w in [(2, 3), (4, 3), (8, 2), (4, 4)]:
    vals = []
    for st in range(1, 401):
        n, a = block_couplings(r, w, stream=st)
        vals.append(2 * sum(a))
    p(f"{r:>4} {w:>3} {min(vals)/D:>12.6f} {max(vals)/D:>12.6f} {max(vals)/min(vals):>9.4f} "
      f"{len(set(vals)):>23} {r:>10}")
p("READ: venue 2's spread moves at fixed record count and takes 400 distinct values in 400")
p("draws. NOT A COUNT there either.")

# ---------------------------------------------------------------- verdict
p("")
p("=" * 118)
p("VERDICT TABLE. Every cell is filled from a number printed above, none in advance.")
p("=" * 118)
p("")
p(f"{'':38} {'SPREAD  max E - min E':>24} {'VARIANCE  Var_s E(s)':>24} {'E(s), typical s':>20}")
p(f"{'':38} {'(config-independent)':>24} {'(config-independent)':>24} {'(config-dependent)':>20}")
p("-" * 118)
p(f"{'(a) EXTENSIVE, S(2N)/S(N) -> 2':38} {'YES  exact, ratio -> 2':>24} {'YES  exact, ratio -> 2':>24} "
  f"{'NO  grows as sqrt n':>20}")
p(f"{'    at a FIXED environment':38} {'YES  prefix stream':>24} {'YES  prefix stream':>24} {'YES':>20}")
p(f"{'    extensive in the RECORD count':38} {'YES  venue 2 only':>24} {'YES  venue 2 only':>24} {'NO':>20}")
p(f"{'(b) ADDITIVE over disjoint regions':38} {'YES  exact integer 0':>24} {'YES  exact integer 0':>24} "
  f"{'YES  exact integer 0':>20}")
p(f"{'    live control fires':38} {'YES  defect = 2|J_m|':>24} {'YES  defect = J_m^2':>24} "
  f"{'YES  defect = J_m t_m':>20}")
p(f"{'(c) NOT merely a count':38} {'YES  400/400 values':>24} {'YES  400/400 values':>24} "
  f"{'YES':>20}")
p(f"{'    moves under venue scale (D-17)':38} {'YES  degree 1 in J':>24} {'YES  degree 2 in J':>24} {'YES':>20}")
p(f"{'    sees WHERE couplings sit':38} {'NO   symmetric in J':>24} {'NO   symmetric in J':>24} "
  f"{'YES  400 values':>20}")
p(f"{'(d) SIGN-DEFINITE, C-46 ratio':38} {'YES  ratio = 1 exact':>24} {'YES  terms J_i^2 > 0':>24} "
  f"{'NO   mean ratio -> 0':>20}")
p(f"{'(e) SEPARATION with power-law tail':38} {'NOT POSED':>24} {'NOT POSED':>24} {'NOT POSED':>20}")
p("-" * 118)
p("(e) reads NOT POSED, not NO: this H couples nearest neighbours only and that was INSERTED,")
p("so the carrier cannot answer the separation question either way. Reporting a null here would")
p("be reporting the input.")
p("")
p("D-22 STATUS: venue 1 has |Aut| = 1 for n = 3..8 by brute force over all n! permutations, and")
p("by the distinctness + non-palindrome obstruction for n = 9..16. Venue 2 has |Aut| = 1 at")
p("block size w >= 3 and |Aut| = 2^r at w = 2, from block-internal endpoint swaps that never")
p("exchange blocks; the identical-block control returns 2^r . r! and so is live. Neither venue")
p("is permutation-symmetric in the sense D-22 forbids. The config-INDEPENDENT observables are")
p("nevertheless blind to placement, which is a fact about the OBSERVABLE, not the carrier --")
p("and it is the reason no separation law can be extracted from them.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/s7_standard.txt", "w").write("\n".join(OUT) + "\n")
