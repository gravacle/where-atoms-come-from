"""S6 -- HOW MANY INDEPENDENT RECORD BITS DOES THE CHAIN ACTUALLY HOLD?

D-18 is not finished at "each Z_i passes (i)-(iv)". n operators can each be a record and still
be the SAME BIT. This script asks the question the extensivity result depends on, and answers
it with the program's own multi-record machinery (RecordModel.commuting_family and
.independently_writable) wherever the dense object is affordable, and by exact level
combinatorics everywhere.

The answer is not the flattering one, and a SECOND VENUE is then built where it is.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import D, couplings, configs, energies_int, dense_H, dense_Z
from record_model import RecordModel, eigenspaces, clause_iii, clause_iv

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 118)
p("S6  INDEPENDENT RECORD BITS  --  n records is not n bits")
p("=" * 118)

def level_stats(n, a):
    s = configs(n); E = energies_int(s, a)
    u, inv = np.unique(E, return_inverse=True)
    mult = np.bincount(inv)
    return len(u), int(mult.min()), int(mult.max())

# ---------------------------------------------------------------- venue 1: the coupled chain
p("")
p("-" * 118)
p("VENUE 1  THE FULLY COUPLED CHAIN, the venue everything above was measured in.")
p("An eigenspace of dimension m can hold at most log2(m) independent record bits, because k")
p("independent bits demand 2^k joint eigenspaces inside the level.")
p("Measured with the program's own machinery where the dense object fits, and by exact level")
p("combinatorics at every n.")
p("-" * 118)
p(f"{'n':>3} {'#records passing (i)-(iv)':>26} {'max level dim':>14} {'log2 -> bit CEILING':>21} "
  f"{'MODEL commuting_family':>23} {'MODEL independently_writable':>29}")
for n in range(2, 17):
    a = couplings(n - 1)
    nl, mn, mx = level_stats(n, a)
    ceil = int(np.log2(mx))
    if n <= 8:
        H = dense_H(n, a)
        m = RecordModel(H)
        Rs = [dense_Z(n, i) for i in range(n)]
        fam = m.commuting_family(Rs)
        iw = m.independently_writable(fam)
        p(f"{n:>3} {n:>26} {mx:>14} {ceil:>21} {len(fam):>23} {len(iw):>29}")
    else:
        p(f"{n:>3} {n:>26} {mx:>14} {ceil:>21} {'(dense too large)':>23} {'(dense too large)':>29}")
p("")
p("READ, AND IT IS THE UNFLATTERING READ: every eigenspace of the fully coupled chain has")
p("dimension EXACTLY 2 at every n. All n of the Z_i satisfy clauses (i)-(iv), but inside any")
p("level they are the SAME BIT up to sign -- the chain holds ONE independent record bit however")
p("long it is. The program's own commuting_family returns a family of size 1, independently of n.")
p("The energy spread is extensive; the number of independent records is NOT.")

p("")
p("WHY: with distinct couplings, the map s -> (s_i s_{i+1})_i is exactly 2-to-1 onto the bond")
p("strings, so each energy level is precisely a pair {s, -s}. The bond correlations LABEL the")
p("levels; the records live in the two-fold degeneracy inside a level and nowhere else.")
p(f"{'n':>3} {'#levels':>9} {'2^(n-1)':>10} {'equal?':>8} {'every level dim 2?':>20} "
  f"{'is Z_iZ_{i+1} non-trivial (clause iii)?':>41}")
for n in range(2, 9):
    a = couplings(n - 1)
    nl, mn, mx = level_stats(n, a)
    H = dense_H(n, a); es = eigenspaces(H)
    ct = any(clause_iii(dense_Z(n, i) @ dense_Z(n, i + 1), es) for i in range(n - 1))
    p(f"{n:>3} {nl:>9} {1<<(n-1):>10} {str(nl==(1<<(n-1))):>8} {str(mn==mx==2):>20} {str(ct):>41}")
p("READ: the number of levels is exactly 2^(n-1) and every level has dimension 2. The pair")
p("correlations FAIL clause (iii) -- they are constant on every eigenspace. So they are")
p("conserved charges, not records; and no admissible unitary can change any of them, since an")
p("admissible U preserves each level and the correlation is a scalar there.")
p("THAT IS AN EXACT ARGUMENT, and it is stronger than a search: it is not that no cheap writer")
p("of the correlation was found, it is that NONE EXISTS.")

# ---------------------------------------------------------------- venue 2: r blocks
p("")
p("-" * 118)
p("VENUE 2  THE SAME CHAIN CUT INTO r DISJOINT BLOCKS (the additivity geometry of S3), which is")
p("where the record count itself becomes extensive. Bonds inside blocks are the same distinct")
p("couplings; the r-1 bonds crossing the cuts are absent.")
p("Clauses (i)-(iv) are re-verified from scratch here (D-18): a new carrier is a new check.")
p("-" * 118)
p(f"{'n':>3} {'r':>3} {'block size':>11} {'max level dim':>14} {'bit CEILING':>12} "
  f"{'MODEL fam':>10} {'MODEL indep-writable':>21} {'(iii) all i':>12} {'(iv) max|Tr|':>13}")
for n, r in [(4, 2), (6, 2), (6, 3), (8, 2), (8, 4), (10, 5), (12, 3), (12, 4), (12, 6),
             (14, 7), (16, 4), (16, 8)]:
    w = n // r
    a_full = couplings(n - 1); a = list(a_full)
    for j in range(1, r): a[j * w - 1] = 0
    nl, mn, mx = level_stats(n, a)
    ceil = int(np.log2(mx))
    s = configs(n); E = energies_int(s, a)
    u, inv = np.unique(E, return_inverse=True)
    c3 = True; mt = 0
    for i in range(n):
        si = s[:, i].astype(np.int64)
        tri = np.zeros(len(u), dtype=np.int64); np.add.at(tri, inv, si)
        mt = max(mt, int(np.abs(tri).max()))
        pos = np.zeros(len(u), dtype=np.int64); neg = np.zeros(len(u), dtype=np.int64)
        np.add.at(pos, inv, (si > 0).astype(np.int64)); np.add.at(neg, inv, (si < 0).astype(np.int64))
        if not np.any((pos > 0) & (neg > 0)): c3 = False
    if n <= 8:
        H = dense_H(n, a); m = RecordModel(H)
        Rs = [dense_Z(n, i) for i in range(n)]
        fam = m.commuting_family(Rs); iw = m.independently_writable(fam)
        fs, iws = str(len(fam)), str(len(iw))
    else:
        fs = iws = "(too large)"
    p(f"{n:>3} {r:>3} {w:>11} {mx:>14} {ceil:>12} {fs:>10} {iws:>21} {str(c3):>12} {mt:>13}")
p("")
p("READ: cutting the chain into r blocks gives every level dimension 2^r and a bit ceiling of r.")
p("The program's own machinery confirms r independent, INDEPENDENTLY WRITABLE records wherever")
p("the dense object fits. Clauses (iii) and (iv) hold for every Z_i on this carrier too --")
p("clause (iv)'s trace is exact integer zero for all n and all i.")
p("So the number of independent records IS extensive in this venue, at r = n/w for fixed block")
p("size w, and S3 already showed the energy spread is EXACTLY additive over exactly these blocks.")
p("")
p("THE PRICE, STATED PLAINLY: in venue 2 the blocks are DECOUPLED, so the extensive spread is")
p("a sum of independent block contributions with no interaction between records in different")
p("blocks. Extensive and additive, yes; interacting, no. That is the same contact-or-nothing")
p("boundary C-47 recorded, met again from the other side.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/s6_independence.txt", "w").write("\n".join(OUT) + "\n")
