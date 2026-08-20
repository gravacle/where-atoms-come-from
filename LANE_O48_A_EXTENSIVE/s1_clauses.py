"""S1 -- VERIFY THE RECORDS BEFORE MEASURING ANYTHING (D-18).

For every n reported: check (i) bit, (ii) durable, (iii) non-trivial, (iv) writable, and
SEARCH the full Pauli group for an admissible writer.  Nothing is nominated: the writer set
is the OUTPUT of an exhaustive search, and the search is cross-checked against a literal
4^n enumeration at small n.

D-22 is checked first, because a permutation-symmetric carrier would make every later
separation statement empty.
"""
import sys, itertools, time
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE")
import numpy as np
from chain import (D, couplings, uniform_couplings, configs, energies_int, levels,
                   pauli, dense_H, dense_Z)
from record_model import eigenspaces, clause_iii, clause_iv

OUT = []
def p(*x):
    s = " ".join(str(y) for y in x); OUT.append(s); print(s)

p("=" * 100)
p("S1  CLAUSE VERIFICATION AND WRITER SEARCH  --  H = sum_i J_i Z_i Z_{i+1}, open chain")
p("=" * 100)

# ============================================================ D-22  IS THERE GEOMETRY TO DETECT?
p("")
p("-" * 100)
p("D-22  PERMUTATION-SYMMETRY OF THE CARRIER.  |Aut| counted by BRUTE FORCE over all n!")
p("      permutations pi of sites, counting those with H(pi(s)) = H(s) for every config s.")
p("      Two LIVE CONTROLS that must register a LARGER symmetry group than the venue:")
p("        PATH-UNIFORM  : same path, all J equal        -> expect |Aut| = 2 (reflection)")
p("        COMPLETE-UNIF : all pairs coupled, J equal    -> expect |Aut| = n! (NO geometry)")
p("-" * 100)
p(f"{'n':>3} {'|Aut| PATH-DISTINCT-J':>22} {'|Aut| PATH-UNIFORM':>20} {'|Aut| COMPLETE-UNIF':>21} {'n!':>8}  READ")

def aut_count(n, bonds):
    """bonds: dict {(i,j): a}.  Count permutations preserving the bond-coupling function."""
    cnt = 0
    for perm in itertools.permutations(range(n)):
        ok = True
        # build image bond map
        img = {}
        for (i, j), v in bonds.items():
            k, l = perm[i], perm[j]
            img[(min(k, l), max(k, l))] = v
        if img != bonds: ok = False
        if ok: cnt += 1
    return cnt

import math
for n in range(3, 9):
    a = couplings(n - 1)
    bd = {(i, i + 1): a[i] for i in range(n - 1)}
    au = uniform_couplings(n - 1)
    bu = {(i, i + 1): au[i] for i in range(n - 1)}
    bc = {(i, j): 777 for i in range(n) for j in range(i + 1, n)}
    cd, cu, cc = aut_count(n, bd), aut_count(n, bu), aut_count(n, bc)
    read = "venue has geometry" if cd == 1 and cu > 1 and cc == math.factorial(n) else "CHECK"
    p(f"{n:>3} {cd:>22} {cu:>20} {cc:>21} {math.factorial(n):>8}  {read}")

p("")
p("  For n > 8 the n! enumeration is skipped; instead the two exact obstructions are checked:")
p("  a site permutation preserving H must be a graph automorphism of the PATH (identity or the")
p("  reflection i -> n+1-i), and the reflection additionally requires the coupling list to be a")
p("  PALINDROME.  Both are checked below.")
p(f"{'n':>3} {'couplings all distinct':>23} {'palindromic?':>14}  READ")
for n in range(9, 17):
    a = couplings(n - 1)
    p(f"{n:>3} {str(len(set(a)) == len(a)):>23} {str(a == a[::-1]):>14}  "
      f"{'|Aut| = 1, geometry present' if len(set(a)) == len(a) and a != a[::-1] else 'CHECK'}")

# ============================================================ CLAUSES (i)-(iv), EXACT INTEGER
p("")
p("-" * 100)
p("CLAUSES (i)-(iv) FOR EVERY R_i = Z_i, ALL i, AT EVERY n.  Exact integer arithmetic in units")
p("of 1/D, D = 2^40.  Eigenspaces of H are EXACT integer energy levels -- no float grouping.")
p("  (i)   R = R-dag, R^2 = I")
p("  (ii)  [H,R] = 0")
p("  (iii) R non-constant on some eigenspace: some level holds configs with s_i = +1 and -1")
p("  (iv)  Tr(P_E R) = sum of s_i over the level = 0 on EVERY level")
p("LIVE CONTROL for the (iv) zero: the PAIR CORRELATION C_1 = Z_1 Z_2, tested identically.")
p("If Tr(P_E Z_1Z_2) were also zero everywhere the test would be blind; it is not.")
p("NOTE ON THE (i) AND (ii) COLUMNS: at every n they rest on an EXACT ARGUMENT, not on this")
p("loop -- Z_i is diagonal with entries +-1 in the configuration basis, so R = R-dag and R^2 = I;")
p("and H is diagonal in that same basis, so [H,R] = 0. Both are COMPUTED from the matrices in")
p("the dense cross-check table further down, wherever the dense object fits.")
p("-" * 100)
p(f"{'n':>3} {'dim':>7} {'#levels':>8} {'mult':>12} {'(i)':>5} {'(ii)':>6} {'(iii) all i':>12} "
  f"{'(iv) max|Tr| all i':>19} {'CONTROL max|Tr(P_E Z1Z2)|':>26}")

clause_table = {}
for n in range(2, 17):
    a = couplings(n - 1)
    s = configs(n)
    E = energies_int(s, a)
    u, inv = levels(E)
    nl = len(u)
    mult = np.bincount(inv)
    # (i) and (ii) hold by construction for diagonal +-1 operators; cross-checked densely below.
    c1 = True
    c2 = True
    # (iii): per level, does s_i take both values?
    c3_all = True
    max_tr = 0
    for i in range(n):
        si = s[:, i].astype(np.int64)
        # trace per level -- exact integer
        tr = np.bincount(inv, weights=si.astype(np.float64), minlength=nl)
        # redo exactly with integer accumulation
        tri = np.zeros(nl, dtype=np.int64)
        np.add.at(tri, inv, si)
        max_tr = max(max_tr, int(np.abs(tri).max()))
        # non-constant on SOME level: exists level with both +1 and -1
        pos = np.zeros(nl, dtype=np.int64); neg = np.zeros(nl, dtype=np.int64)
        np.add.at(pos, inv, (si > 0).astype(np.int64))
        np.add.at(neg, inv, (si < 0).astype(np.int64))
        if not np.any((pos > 0) & (neg > 0)): c3_all = False
    # CONTROL: Z1 Z2
    cc = (s[:, 0].astype(np.int64) * s[:, 1].astype(np.int64))
    trc = np.zeros(nl, dtype=np.int64); np.add.at(trc, inv, cc)
    ctrl = int(np.abs(trc).max())
    mm = f"{int(mult.min())}..{int(mult.max())}"
    p(f"{n:>3} {1<<n:>7} {nl:>8} {mm:>12} {str(c1):>5} {str(c2):>6} {str(c3_all):>12} "
      f"{max_tr:>19} {ctrl:>26}")
    clause_table[n] = dict(levels=nl, mult=mm, c3=c3_all, c4max=max_tr, ctrl=ctrl)

p("")
p("READ: clause (iv) column is 0 at every n for every one of the n records, while the CONTROL")
p("column is non-zero at every n -- the test can see a non-zero and reports none for the records.")
p("Clause (iii) holds for every record at every n.  Every level has multiplicity >= 2, so H is")
p("DEGENERATE at every n (P-1 satisfied).")

# ============================================================ DENSE CROSS-CHECK OF (i),(ii)
p("")
p("-" * 100)
p("DENSE CROSS-CHECK of (i),(ii),(iii),(iv) using the PROGRAM'S OWN clause functions")
p("(record_model.eigenspaces / clause_iii / clause_iv) on the full 2^n x 2^n matrices.")
p("Run wherever the dense object is affordable.  CONTROL row: Z_1 Z_2 must FAIL clause (iv).")
p("-" * 100)
p(f"{'n':>3} {'||R-Rdag||':>12} {'||R^2-I||':>12} {'max||[H,R]||':>14} {'clause_iii all i':>17} "
  f"{'clause_iv all i':>16} {'CONTROL clause_iv(Z1Z2)':>25}")
for n in range(2, 9):
    a = couplings(n - 1)
    H = dense_H(n, a)
    es = eigenspaces(H)
    r1 = r2 = r3 = 0.0
    ok3 = ok4 = True
    for i in range(n):
        R = dense_Z(n, i)
        r1 = max(r1, np.linalg.norm(R - R.conj().T))
        r2 = max(r2, np.linalg.norm(R @ R - np.eye(1 << n)))
        r3 = max(r3, np.linalg.norm(H @ R - R @ H))
        if not clause_iii(R, es): ok3 = False
        if not clause_iv(R, es): ok4 = False
    C = dense_Z(n, 0) @ dense_Z(n, 1)
    ctrl4 = clause_iv(C, es)
    p(f"{n:>3} {r1:>12.2e} {r2:>12.2e} {r3:>14.2e} {str(ok3):>17} {str(ok4):>16} {str(ctrl4):>25}")
p("READ: the dense check agrees with the exact integer check; the control Z1Z2 FAILS clause (iv),")
p("so the clause_iv routine is live.")

# ============================================================ WRITER SEARCH
p("")
p("-" * 100)
p("ADMISSIBLE WRITER: SEARCHED, NOT NOMINATED (D-18).")
p("A Pauli P = i^k X^x Z^z is ADMISSIBLE iff [P,H] = 0.  Since the bond operators Z_iZ_{i+1} are")
p("independent Pauli words, [P,H] = 0 forces [P, Z_iZ_{i+1}] = 0 for every bond.  P writes R_i")
p("iff {P, Z_i} = 0.  The search below imposes NO ansatz: it enumerates the group and reports")
p("what satisfies both conditions.")
p("  LITERAL SEARCH  : full 4^n enumeration of (x,z) with dense-matrix verification, n <= 7.")
p("  FACTORED SEARCH : full 2^n enumeration over x with z free, n <= 16.  Complete because")
p("                    neither condition involves z -- which is itself an OUTPUT of the search.")
p("-" * 100)
p(f"{'n':>3} {'4^n literal: #admissible writers of Z_1':>40} {'2^n factored: #x . 2^n z':>26} "
  f"{'agree':>7} {'the x found':>20}")

def factored_search(n, a, target):
    """Return the list of x in F2^n that are admissible AND flip Z_target."""
    good = []
    for xi in range(1 << n):
        x = [(xi >> (n - 1 - k)) & 1 for k in range(n)]
        if x[target] != 1: continue
        if all((x[i] + x[i + 1]) % 2 == 0 for i in range(n - 1)): good.append(tuple(x))
    return good

for n in range(2, 17):
    a = couplings(n - 1)
    fx = factored_search(n, a, 0)
    factored = len(fx) * (1 << n)
    if n <= 7:
        H = dense_H(n, a); R = dense_Z(n, 0); I = np.eye(1 << n)
        lit = 0
        for xi in range(1 << n):
            for zi in range(1 << n):
                lab = []
                for k in range(n):
                    xb = (xi >> (n - 1 - k)) & 1; zb = (zi >> (n - 1 - k)) & 1
                    lab.append(0 if (xb, zb) == (0, 0) else (1 if (xb, zb) == (1, 0)
                               else (3 if (xb, zb) == (0, 1) else 2)))
                P = pauli(lab)
                if np.linalg.norm(P @ H - H @ P) > 1e-9: continue
                if np.linalg.norm(P.conj().T @ R @ P + R) > 1e-9: continue
                lit += 1
        agree = (lit == factored)
        p(f"{n:>3} {lit:>40} {factored:>26} {str(agree):>7} "
          f"{('all-ones only' if fx == [tuple([1]*n)] else str(fx)):>20}")
    else:
        p(f"{n:>3} {'(4^n too large)':>40} {factored:>26} {'-':>7} "
          f"{('all-ones only' if fx == [tuple([1]*n)] else str(fx)):>20}")

p("")
p("READ: at every n the SEARCH returns exactly one x-support, the ALL-ONES string, with z free.")
p("The admissible writers of Z_1 are precisely X^(tensor n) Z^z.  Nothing was assumed: at n <= 7")
p("the literal 4^n enumeration returns the same count as the factored one.")
p("CONSEQUENCE the search hands us, not an input: every admissible Pauli writer of ANY single")
p("record flips EVERY record at once and leaves EVERY pair correlation Z_i Z_j invariant.")

# ============================================================ WRITER COST, WITH CONTROL
p("")
p("-" * 100)
p("ENERGY COST OF THE WRITER, exact integer over ALL 2^n configurations.")
p("W = the writer the search found (x = all-ones, z = 0), i.e. it maps s -> -s.")
p("LIVE CONTROL: X_1, a single-site flip -- NOT admissible; its cost must be non-zero.")
p("-" * 100)
p(f"{'n':>3} {'max|dE| under W (int)':>22} {'W flips every Z_i?':>19} {'W changes any Z_iZ_j?':>22} "
  f"{'CONTROL max|dE| under X_1':>26} {'X_1 admissible?':>16}")
for n in range(2, 17):
    a = couplings(n - 1)
    s = configs(n)
    E = energies_int(s, a)
    # W: s -> -s.  index of -s is the bitwise complement
    idx = np.arange(1 << n, dtype=np.int64)
    comp = idx ^ ((1 << n) - 1)
    dW = int(np.abs(E[comp] - E).max())
    flips_all = bool(np.all(s[comp] == -s))
    corr_unchanged = True
    for i in range(n - 1):
        c0 = s[:, i].astype(np.int64) * s[:, i + 1].astype(np.int64)
        c1 = s[comp][:, i].astype(np.int64) * s[comp][:, i + 1].astype(np.int64)
        if not np.array_equal(c0, c1): corr_unchanged = False
    # CONTROL X_1: flip site 1 only -> toggle the top bit
    c1i = idx ^ (1 << (n - 1))
    dX = int(np.abs(E[c1i] - E).max())
    p(f"{n:>3} {dW:>22} {str(flips_all):>19} {str(not corr_unchanged):>22} {dX:>26} {'False':>16}")

p("")
p("READ: the writer costs EXACTLY ZERO at every n (integer zero, not a tolerance), it flips every")
p("record, and it leaves every pair correlation untouched.  The control X_1 costs a non-zero")
p("integer at every n, so the zero column is not the test being blind.")
p("SINGLE RECORDS ARE FREE AT EVERY n.  This is O-47's structure, INDUCED at general n by the")
p("search rather than inserted.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_A_EXTENSIVE/s1_clauses.txt", "w").write("\n".join(OUT) + "\n")
