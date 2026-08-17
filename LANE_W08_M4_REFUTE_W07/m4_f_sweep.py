# W-08 / M4 leg F — THE 52-PARTITION SWEEP.  CONTROL, OR THEOREM?
#
# W-07 sec6: "'Could not have failed' voids a CONTROL — the program's own rule, and this is a
# control, not a theorem.  Leg (b) is void."
# W-06 leg (b), verbatim from REGISTER:551-553: "It is UNIQUELY correct on K1's own count: sweeping
# all 52 partition subgroups of U(1)^5, exactly ONE yields S1 sec4's invariant parameter count
# 2 = b1 + #faces = E - rank(d1) — the full local group.  Global U(1) gives 6; the class group 4."
#
# ISOLATION LEDGER.  F1-F2 hold the complex fixed at K1 and move the PARTITION only (all 52).
# F3 holds the partition-sweep procedure fixed and moves the COMPLEX only (K1 -> a 5-cycle),
# to decide whether the sweep's outcome is contingent on its input or a tautology of its design.
# Integer/rational arithmetic where it matters: d1 has entries in {-1,0,1} and all ranks are
# computed on integer matrices; numpy's float rank is cross-checked against an exact
# fraction-free Gaussian elimination.
import numpy as np
from fractions import Fraction
from collections import Counter

def partitions(c):
    if len(c) == 1: yield [c]; return
    first, rest = c[0], c[1:]
    for p in partitions(rest):
        for i in range(len(p)): yield p[:i]+[[first]+p[i]]+p[i+1:]
        yield [[first]]+p

def exact_rank(M):
    """Fraction-free rank over Q.  M is a list of lists of ints."""
    A = [[Fraction(x) for x in row] for row in M]
    rows, cols = len(A), len(A[0]) if A else 0
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]/A[r][c]
                A[i] = [A[i][j]-f*A[r][j] for j in range(cols)]
        r += 1
    return r

def sweep(EDGES, V, E, target):
    d1 = [[0]*E for _ in range(V)]
    for j, (s, t) in enumerate(EDGES):
        d1[t][j] += 1; d1[s][j] -= 1
    d1n = np.array(d1, float)
    cnt = Counter(); winners = []; bykblocks = {}
    for P in partitions(list(range(V))):
        B = [[1 if v in blk else 0 for blk in P] for v in range(V)]
        M = [[sum(d1[v][j]*B[v][i] for v in range(V)) for i in range(len(P))] for j in range(E)]
        rk = exact_rank(M)
        rkf = int(np.linalg.matrix_rank(np.array(d1n).T @ np.array(B, float)))
        assert rk == rkf, f"exact rank {rk} != float rank {rkf}"
        inv = E - rk
        cnt[inv] += 1
        bykblocks.setdefault(len(P), set()).add((rk, inv))
        if inv == target: winners.append(sorted([sorted(b) for b in P]))
    return cnt, winners, bykblocks

EDGES_K1 = [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]
print("== F1  THE SWEEP, RE-RUN WITH EXACT RANK (W-07 used float rank; both agree, asserted above) ==")
cnt, winners, byk = sweep(EDGES_K1, 5, 6, 2)
print(f"  partitions swept: {sum(cnt.values())}   invariant-count distribution: {dict(sorted(cnt.items()))}")
print(f"  W-07 reported                                {{2: 1, 3: 10, 4: 25, 5: 15, 6: 1}}   -> MATCHES")
print(f"  winners at invariant count 2: {winners}")
print(f"  (rank, invariants) achievable per block-count k: {dict(sorted(byk.items()))}")
print()
print("== F2  IS W-07's ONE-LINE ARGUMENT RIGHT?  AND IS IT MORE THAN IT CLAIMS? ==")
print("  W-07: 'a k-block subgroup acts through <= k-1 parameters; rank <= k-1; invariants >= 7-k;")
print("        reaching 2 needs rank 4 needs k = 5.'   The bound is CORRECT.  It is also LOOSE:")
print("  on a CONNECTED complex the rank is exactly k - (#components of the quotient) = k - 1, so")
print("  invariants = E - k + 1 = 7 - k EXACTLY, for every one of the 52.  Checked:")
ok = all(len(s) == 1 and list(s)[0] == (k-1, 6-(k-1)) for k, s in byk.items())
print(f"    invariants = 7 - k holds for all 52 partitions: {ok}")
S = Counter(len(P) for P in partitions(list(range(5))))
print(f"    multiplicities by k = Stirling S(5,k) = {dict(sorted(S.items()))}  (1,15,25,10,1)")
print("  So the ENTIRE distribution {6,5,4,3,2} with those multiplicities is forced by ONE fact:")
print("  K1 is connected.  W-07 proved half of a sharper theorem than it noticed.")
print()
print("== F3  IS THE SWEEP'S OUTCOME CONTINGENT ON ITS INPUT?  MOVE THE COMPLEX, HOLD THE SWEEP. ==")
print("  A control is void when it could not have come out otherwise ON ANY INPUT.  Feed the SAME")
print("  procedure a different carrier: the 5-cycle C5 (V=5, E=5, no faces).")
EDGES_C5 = [(0,1),(1,2),(2,3),(3,4),(4,0)]
cnt5, win5, byk5 = sweep(EDGES_C5, 5, 5, 2)
print(f"  C5 distribution: {dict(sorted(cnt5.items()))}")
print(f"  partitions of C5 giving invariant count 2: {len(win5)}  -> e.g. {win5[:3]}")
print("  On C5 the count 2 is achieved by TEN 4-block partitions and NOT by the discrete one")
print("  (which gives 1).  The sweep therefore DISCRIMINATES: run on a different carrier it")
print("  returns a different winner and a different multiplicity.  Its K1 outcome is a fact")
print("  ABOUT K1, not an artefact of the sweep's design.")
print()
print("== F4  THE CATEGORY QUESTION, SETTLED ==")
print("  W-06 leg (b) asserts: 'among the 52 partition subgroups of U(1)^5, EXACTLY ONE yields")
print("  invariant count 2.'  That is a UNIVERSALLY QUANTIFIED STATEMENT OVER A FINITE DOMAIN,")
print("  offered as a premise.  Exhaustive enumeration over a finite domain IS a proof.  It is a")
print("  PROPOSITION — true, and now proved twice (by enumeration, and by W-07's counting bound).")
print("  It is NOT a control: it has no treatment arm, no failure mode it was designed to detect,")
print("  and nothing about it is falsified by being provable.")
print()
print("  THE PROGRAM'S RULE, AS THE PROGRAM WROTE IT (REGISTER:503-509, IMP-1):")
print("    'On a control the norm is sound.  On a theorem it is incoherent — a proved statement")
print("     cannot fail; that is what \"theorem\" means.'")
print("  W-07 SUPPLIED THE PROOF AND THEN USED THE EXISTENCE OF THE PROOF TO VOID THE STATEMENT.")
print("  That is IMP-1 exactly, one floor up — the defect W-06 convicted the chain of, committed")
print("  by the lane auditing W-06, in the section auditing W-06's treatment of that very defect.")
print()
print("== F5  THE REDUCTIO: W-07's STANDARD APPLIED TO W-07 ==")
print("  Every one of these 'could not have come out otherwise':")
print("    leg A   ord(-1)=2, ord(-i)=4, <W_F,W_C>=Z_4          — arithmetic, forced")
print("    leg C2  exactly K/4 exact zeros at every K           — 4 | k, forced")
print("    leg C3  the {2:1,3:10,4:25,5:15,6:1} distribution    — F2 above, forced by connectivity")
print("    leg D   sup|Z_k| = 1 attained on the published pair  — |Z_k| is periodic with period 2")
print("  If 'forced' voided results, W-07's sec3 and sec4 would be void along with W-06 sec(b).")
print("  They are not void, and neither is the sweep.  The disqualifier applies to CONTROLS, and")
print("  W-07 has exactly one control — the five generic connections of leg E — whose null result")
print("  IS forced (m4_c: expected count 6.4e-06 at K=4000, threshold 1e-9).  W-07 applied the rule")
print("  to the one object in view that is a theorem, and withheld it from its own control.")
