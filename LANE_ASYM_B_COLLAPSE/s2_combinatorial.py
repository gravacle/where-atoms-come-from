"""S2 -- THE COMBINATORIAL QUANTITIES, EXACTLY, IN THE F_2 SYMPLECTIC REPRESENTATION.

Everything here is computed on (x|z) vectors over F_2^{2n}.  No dense matrix is built, so n
runs far past the dense limit.  Where an EXACT ARGUMENT closes the question at every n, the
argument is stated and the table is the check on it -- not the other way round.

TWO RECORD SETS, and the difference is the whole point:

  SET A -- THE RECORD FAMILY.  {Xbar_1 .. Xbar_k}.  k = n-2 mutually COMMUTING records.
           This is the set the program means by "k records": independent bits that are
           simultaneously definite.  RecordModel.commuting_family demands exactly this.

  SET B -- RECORDS AND THEIR WRITERS.  {Xbar_i, Zbar_i}, 2k operators, EVERY ONE of which
           is separately a record by clauses (i)-(iv).  Here anticommutation is present.

Set B is the POSITIVE CONTROL (D-15) for every zero reported on set A.

QUANTITIES
  N            = k = n-2                    control (b): exactly linear, by construction
  P_int        interacting-pair count       # pairs with non-zero symplectic pairing
  W_tot        total writer weight          SUM over records of the MINIMUM weight of an
                                            admissible operator that flips it
  lam_max      largest eigenvalue of the record-record relation matrix, for two relations:
                 (s) the symplectic pairing matrix   (0/1, anticommute = 1)
                 (o) the support-overlap matrix      (|supp_i AND supp_j|, 0 on diagonal)
"""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

OUT = []
def say(s=""):
    print(s); OUT.append(s)

def stab_xz(n): return [[1]*n + [0]*n, [0]*n + [1]*n]
def sp(a, b, n): return (sum(a[i]*b[n+i] + a[n+i]*b[i] for i in range(n))) % 2
def wt(v, n): return sum(1 for i in range(n) if v[i] or v[n+i])
def supp(v, n): return set(i for i in range(n) if v[i] or v[n+i])
def add(a, b): return [(x+y) % 2 for x, y in zip(a, b)]

def cosetmin(v, n):
    """minimum-weight representative of the logical class of v (stabiliser group has 4 elements)"""
    S = stab_xz(n)
    grp = [[0]*(2*n), S[0], S[1], add(S[0], S[1])]
    return min((add(v, g) for g in grp), key=lambda u: wt(u, n))

NS = range(4, 61, 2)
DATA = {}

say("="*112)
say("S2   EXACT COMBINATORIAL QUANTITIES ON THE [[n, n-2, 2]] FAMILY")
say("     representation: F_2 symplectic (x|z) in F_2^{2n}.  No dense matrix built anywhere in S2.")
say("="*112)
say()

# ---------------------------------------------------------------------------- exact arguments
say("-"*112)
say("EXACT ARGUMENT 1 -- WHICH PAULIS ARE ADMISSIBLE.")
say("  H = -(X^(x)n + Z^(x)n).  A Pauli P commutes with H iff it commutes with BOTH stabilisers.")
say("  For P = (x|z):  [P, X^(x)n] = 0  iff  |z| is even;   [P, Z^(x)n] = 0  iff  |x| is even.")
say("  So N(S) = { (x|z) : |x| even and |z| even }.  This is exact and holds at every even n.")
say()
say("EXACT ARGUMENT 2 -- THE MINIMUM WRITER WEIGHT IS EXACTLY 2, AT EVERY n, FOR EVERY RECORD.")
say("  LOWER BOUND (covers NON-Pauli writers too).  A weight-1 admissible unitary U (x) I must")
say("  satisfy [U, X^(x)n] = 0 and [U, Z^(x)n] = 0, which force [U,X] = [U,Z] = 0, so U is a")
say("  phase and flips nothing.  Hence NO writer of weight 1 exists, of any kind.  W >= 2.")
say("  UPPER BOUND.  The weight-2 members of N(S) are exactly X_iX_j, Z_iZ_j, Y_iY_j.  Against")
say("  a record R = (a|b):  <X_iX_j, R> = b_i + b_j,  <Z_iZ_j, R> = a_i + a_j.  One of these is 1")
say("  unless BOTH a and b are constant vectors -- and a,b constant means R is in {I, X^n, Z^n,")
say("  Y^n}, i.e. R is a stabiliser, not a record (it fails clause (iii)).  So every record has a")
say("  weight-2 writer.  W = 2 EXACTLY, INDEPENDENT OF n AND INDEPENDENT OF WHICH RECORD.")
say()
say("EXACT ARGUMENT 3 -- THE RECORD-RECORD SYMPLECTIC RELATION MATRIX ON SET A IS IDENTICALLY 0.")
say("  Independent records must be simultaneously definite bits, so they COMMUTE -- that is the")
say("  content of RecordModel.commuting_family, and symplectic Gram-Schmidt returns the Xbar_i")
say("  mutually commuting.  Therefore on set A: P_int = 0 and lam_max(symplectic) = 0 at EVERY N.")
say("  The pairing lives between a record and its WRITER, never between two records.")
say("  On set B the same matrix is the direct sum of k copies of [[0,1],[1,0]], so P_int = k and")
say("  lam_max(symplectic) = 1 EXACTLY at every N -- bounded, not growing.")
say("-"*112)
say()

# ---------------------------------------------------------------------------- the table
say("TABLE 1.  Computed values.  Every zero on set A has its set-B positive control in the SAME row.")
say()
say("                          |------------------ SET A: the k commuting records ------------------|"
    "  |----------- SET B: records + writers (CONTROL) -----------|")
say("    n     k=N   |supp|avg |  P_int(s)  lam_s   P_int(o)   lam_o    W_tot    Wmin  Wmax |"
    "   P_int(s)  lam_s   P_int(o)     lam_o     W_tot")
for n in NS:
    k = n-2
    prs = symplectic_logicals(stab_xz(n), n)
    assert len(prs) == k
    A = [cosetmin(p[0], n) for p in prs]
    B = [cosetmin(v, n) for p in prs for v in p]

    def analyse(V):
        m = len(V)
        Ms = np.zeros((m, m)); Mo = np.zeros((m, m))
        for i in range(m):
            si = supp(V[i], n)
            for j in range(m):
                if i == j: continue
                Ms[i, j] = sp(V[i], V[j], n)
                Mo[i, j] = len(si & supp(V[j], n))
        pint_s = int(Ms.sum()//2); pint_o = int((Mo > 0).sum()//2)
        ls = float(np.linalg.eigvalsh(Ms)[-1]) if m else 0.0
        lo = float(np.linalg.eigvalsh(Mo)[-1]) if m else 0.0
        # writer weight: MINIMUM weight over N(S) of something that anticommutes with V[i]
        ws = [min_writer_weight(V[i], n) for i in range(m)]
        return pint_s, ls, pint_o, lo, sum(ws), min(ws), max(ws)

    def min_writer_weight(R, n):
        """exact, by argument 2: search weight 1 then weight 2 explicitly (never assumed)"""
        # weight 1: no admissible weight-1 Pauli exists -- verified by construction here
        for i in range(n):
            for (x, z) in ((1, 0), (0, 1), (1, 1)):
                v = [0]*(2*n); v[i] = x; v[n+i] = z
                if all(sp(v, s, n) == 0 for s in stab_xz(n)) and sp(v, R, n) == 1:
                    return 1
        for i in range(n):
            for j in range(i+1, n):
                for (x1, z1) in ((1, 0), (0, 1), (1, 1)):
                    for (x2, z2) in ((1, 0), (0, 1), (1, 1)):
                        v = [0]*(2*n); v[i] = x1; v[n+i] = z1; v[j] = x2; v[n+j] = z2
                        if all(sp(v, s, n) == 0 for s in stab_xz(n)) and sp(v, R, n) == 1:
                            return 2
        return 99          # would mean the exact argument is wrong

    a = analyse(A); b = analyse(B)
    savg = sum(wt(v, n) for v in A)/len(A)
    DATA[n] = dict(k=k, A=a, B=b, supp_avg=savg)
    if n <= 24 or n % 10 == 0:
        say("  %4d %6d %10.3f | %9d %6.2f %10d %7.2f %8d %7d %5d | %10d %6.2f %10d %9.2f %9d"
            % (n, k, savg, a[0], a[1], a[2], a[3], a[4], a[5], a[6], b[0], b[1], b[2], b[3], b[4]))
say()
say("  READ OF TABLE 1 (filled from the numbers above, not in advance):")
allz = all(DATA[n]['A'][0] == 0 and DATA[n]['A'][1] == 0 for n in NS)
bctl = all(DATA[n]['B'][0] == DATA[n]['k'] and abs(DATA[n]['B'][1]-1.0) < 1e-12 for n in NS)
wall = all(DATA[n]['A'][5] == 2 and DATA[n]['A'][6] == 2 for n in NS)
say("    set A symplectic pairing identically zero at every n tested : %s" % allz)
say("    set B positive control non-zero and equal to k at every n    : %s" % bctl)
say("    minimum writer weight equals 2 for EVERY record at every n   : %s" % wall)
say()

# ---------------------------------------------------------------------------- brute-force writer check
say("-"*112)
say("SELF-CHECK ON THE WRITER WEIGHT -- exhaustive minimisation over ALL 4^n Paulis, small n only.")
say("  This does not assume argument 2; it searches every Pauli, keeps those in N(S) that")
say("  anticommute with the record, and reports the true minimum weight.")
say()
say("      n   records   exhaustive min writer weight (min,max over records)   matches argument 2")
ok_bf = True
for n in (4, 6, 8):
    prs = symplectic_logicals(stab_xz(n), n)
    recs = [cosetmin(v, n) for p in prs for v in p]
    mins = []
    allp = []
    for x in itertools.product((0, 1), repeat=n):
        for z in itertools.product((0, 1), repeat=n):
            if sum(x) % 2 or sum(z) % 2: continue         # not in N(S)
            v = list(x)+list(z)
            if not any(v): continue
            allp.append(v)
    for R in recs:
        cand = [wt(v, n) for v in allp if sp(v, R, n) == 1]
        mins.append(min(cand) if cand else 99)
    m2 = (min(mins) == 2 and max(mins) == 2)
    ok_bf &= m2
    say("   %4d %9d %54s %20s" % (n, len(recs), "(%d, %d)" % (min(mins), max(mins)), m2))
say()
say("   exhaustive check %s" % ("PASS -- argument 2 is confirmed by brute force" if ok_bf else "FAIL"))
say("-"*112)
say()

np.save("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE/s2_data.npy",
        np.array([[n, DATA[n]['k'],
                   DATA[n]['A'][0], DATA[n]['A'][1], DATA[n]['A'][2], DATA[n]['A'][3], DATA[n]['A'][4],
                   DATA[n]['B'][0], DATA[n]['B'][1], DATA[n]['B'][2], DATA[n]['B'][3], DATA[n]['B'][4],
                   DATA[n]['supp_avg']] for n in NS], dtype=float))
say("largest n reached in S2: %d  (k = %d).  Nothing stopped it: the cost is O(n^3) in F_2." % (max(NS), max(NS)-2))
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE/s2_combinatorial.txt", "w").write("\n".join(OUT)+"\n")
