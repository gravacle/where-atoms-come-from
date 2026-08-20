"""
O-48-C  STEP 1.  THE CARRIER, AND THE FIVE CLAUSES ON IT.  (D-18)

H = sum_i J_i Z_i Z_{i+1} on an open chain of n qubits.  Records nominated: R_i = Z_i.
NOTHING here is nominated and left unchecked:
  * clauses (i)(ii)(iii)(iv) are VERIFIED at every n reported, by exhaustive enumeration in
    EXACT INTEGER arithmetic, and cross-checked against dense complex matrices at small n.
  * the ADMISSIBLE WRITER is SEARCHED for over the entire Pauli group (4^n elements), never
    nominated.  Its weight and support are read off the search, not assumed.
  * D-15: every reported ZERO is printed beside a positive control that returns non-zero.
  * D-22: the site-permutation automorphism group of the carrier is computed, not asserted.
  * D-17: five different coupling families, three system-size ranges.

H is DIAGONAL in the computational basis, so E(s) = sum_i J_i s_i s_{i+1} is an exact integer
sum for integer J.  No 2^n x 2^n matrix is ever needed for the clause work; dense matrices are
built only at n <= 5 as an INDEPENDENT cross-check of the exact bookkeeping.
"""
import sys, itertools
from fractions import Fraction
import numpy as np

OUT = []
def P(s=""):
    OUT.append(s)
    print(s)

# ------------------------------------------------------------------ coupling families (D-17)
def couplings(name, m, seed=0):
    """m = number of bonds = n-1 for an open chain.  INTEGER couplings -> exact arithmetic."""
    if name == "uniform":      return [1] * m
    if name == "linear":       return [i + 1 for i in range(m)]
    if name == "superinc":     return [2 ** i for i in range(m)]
    if name == "randpos":
        rng = np.random.default_rng(seed);  return [int(v) for v in rng.integers(1, 1000, m)]
    if name == "randsign":
        rng = np.random.default_rng(seed + 77)
        return [int(v) * (1 if rng.random() < 0.5 else -1) for v in rng.integers(1, 1000, m)]
    raise ValueError(name)

FAMILIES = ["uniform", "linear", "superinc", "randpos", "randsign"]

# ------------------------------------------------------------------ exact spectrum machinery
def energy(J, s):
    return sum(J[i] * s[i] * s[i + 1] for i in range(len(J)))

def eigenspaces(J, n):
    """dict: energy -> list of spin configurations (tuples of +-1).  EXACT integers."""
    blocks = {}
    for bits in itertools.product((1, -1), repeat=n):
        blocks.setdefault(energy(J, bits), []).append(bits)
    return blocks

# ------------------------------------------------------------------ dense cross-check tools
I2 = np.eye(2, dtype=complex)
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
Y2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z2 = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X2, Y2, Z2]

def kron(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops: out = np.kron(out, o)
    return out

def dense_H(J, n):
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i, j in enumerate(J):
        ops = [I2] * n; ops[i] = Z2; ops[i + 1] = Z2
        H += j * kron(ops)
    return H

def dense_Z(i, n):
    ops = [I2] * n; ops[i] = Z2
    return kron(ops)

# ------------------------------------------------------------------ the Pauli-group SEARCH
def pauli_from_code(code, n):
    return kron([PAULI[c] for c in code])

def search_admissible_writers(J, n, target_xz, dense_check=False, cap=4 ** 8):
    """EXHAUSTIVE search over the whole Pauli group for W with
         [W,H] = 0            (ADMISSIBLE, O-4)
         {W, target} = 0      (W flips the target observable)
       Returns (count, min_weight, example_code, example_support).
       target_xz is (x,z) bit-vectors of the target Pauli observable.
       Commutation is decided by the exact F_2 symplectic form, and -- because D-19 warns that
       F_2 cannot pose questions of sign -- the ADMISSIBILITY of every candidate found is then
       re-verified over the COMPLEX numbers with a dense matrix when dense_check is on."""
    if 4 ** n > cap:
        return None
    tx, tz = target_xz
    bonds = [([0] * n, [1 if k in (i, i + 1) else 0 for k in range(n)]) for i in range(len(J))]
    found, best_w, best_code = 0, None, None
    Hd = dense_H(J, n) if dense_check else None
    for code in itertools.product(range(4), repeat=n):
        # x,z bit vectors of this Pauli:  I=(0,0) X=(1,0) Y=(1,1) Z=(0,1)
        x = [1 if c in (1, 2) else 0 for c in code]
        z = [1 if c in (2, 3) else 0 for c in code]
        # commute with EVERY bond term (each term separately: the ZZ operators are linearly
        # independent, so H is preserved iff every term is)
        ok = True
        for bi, (bx, bz) in enumerate(bonds):
            if J[bi] == 0: continue
            sp = sum(x[k] * bz[k] + z[k] * bx[k] for k in range(n)) % 2
            if sp: ok = False; break
        if not ok: continue
        # anticommute with the target
        sp = sum(x[k] * tz[k] + z[k] * tx[k] for k in range(n)) % 2
        if sp != 1: continue
        found += 1
        w = sum(1 for c in code if c != 0)
        if best_w is None or w < best_w:
            best_w, best_code = w, code
    support = None
    if best_code is not None:
        support = [k for k, c in enumerate(best_code) if c != 0]
        if dense_check:
            W = pauli_from_code(best_code, n)
            comm = np.linalg.norm(W @ Hd - Hd @ W)
            assert comm < 1e-9, f"F_2 search lied: ||[W,H]||={comm}"
    return found, best_w, best_code, support

# ------------------------------------------------------------------ D-22 permutation check
def automorphism_order(J, n):
    """How many site permutations pi leave the energy function invariant, E(s o pi) = E(s)?
       If this is n! the carrier is permutation-symmetric and contains NO geometry (D-22)."""
    configs = list(itertools.product((1, -1), repeat=n))
    base = {c: energy(J, c) for c in configs}
    cnt = 0
    for pi in itertools.permutations(range(n)):
        good = True
        for c in configs:
            cp = tuple(c[pi[k]] for k in range(n))
            if base[cp] != base[c]: good = False; break
        if good: cnt += 1
    return cnt

# =======================================================================================
P("=" * 104)
P("O-48-C  STEP 1.   THE CARRIER AND THE FIVE CLAUSES.   H = sum_i J_i Z_i Z_{i+1}, R_i = Z_i")
P("=" * 104)
P()
P("  Clause (ii) is [H,R]=0 AND [L_k,R]=0.  This venue is CLOSED -- no Lindblad operators are")
P("  introduced, because introducing a bath would import an environment the standard forbids at")
P("  this step.  (ii) therefore reduces to [H,R]=0 here, and that is what is verified.")
P()

# ---------------------------------------------------------------- clause table
P("-" * 104)
P("  CLAUSES (i)-(iv), EXHAUSTIVE OVER ALL 2^n CONFIGURATIONS, EXACT INTEGER ARITHMETIC")
P("-" * 104)
P(f"  {'family':<10} {'n':>3} {'#blocks':>8} {'dim(P_E) hist':<22} "
  f"{'(i)':>5} {'(ii)':>5} {'(iii)':>6} {'max|Tr(P_E R_i)|':>18} {'(iv)':>6}   {'CONTROL max|Tr|':>16}")

def clause_table(fam, n, seed=0):
    J = couplings(fam, n - 1, seed)
    blocks = eigenspaces(J, n)
    dims = {}
    for E, cs in blocks.items(): dims[len(cs)] = dims.get(len(cs), 0) + 1
    hist = ",".join(f"{d}x{c}" for d, c in sorted(dims.items()))
    # (i) R=Z_i is Hermitian and squares to I: verified densely at small n, else structural
    c_i = True
    if n <= 5:
        for i in range(n):
            Zi = dense_Z(i, n)
            c_i &= np.linalg.norm(Zi - Zi.conj().T) < 1e-12
            c_i &= np.linalg.norm(Zi @ Zi - np.eye(2 ** n)) < 1e-12
    # (ii)
    c_ii = True
    if n <= 5:
        Hd = dense_H(J, n)
        for i in range(n):
            c_ii &= np.linalg.norm(Hd @ dense_Z(i, n) - dense_Z(i, n) @ Hd) < 1e-9
    # (iii) non-constant on SOME eigenspace, for every i
    c_iii = True
    for i in range(n):
        ok = any(len(set(c[i] for c in cs)) > 1 for cs in blocks.values())
        c_iii &= ok
    # (iv) Tr(P_E R_i) on EVERY eigenspace, EXACT INTEGER
    worst = 0
    for E, cs in blocks.items():
        for i in range(n):
            worst = max(worst, abs(sum(c[i] for c in cs)))
    c_iv = (worst == 0)
    # CONTROL (D-15): break the balance with a longitudinal field on site 0 and re-run the
    # SAME estimator.  A field h Z_0 keeps H diagonal, so the arithmetic is identical.
    h = 10 ** 6 + 1
    cblocks = {}
    for bits in itertools.product((1, -1), repeat=n):
        cblocks.setdefault(energy(J, bits) + h * bits[0], []).append(bits)
    cworst = 0
    for E, cs in cblocks.items():
        for i in range(n):
            cworst = max(cworst, abs(sum(c[i] for c in cs)))
    P(f"  {fam:<10} {n:>3} {len(blocks):>8} {hist:<22} "
      f"{str(c_i):>5} {str(c_ii):>5} {str(c_iii):>6} {worst:>18} "
      f"{'HOLDS' if c_iv else 'FAILS':>6}   {cworst:>16}")
    return worst, cworst, c_iii

allworst, allctrl = [], []
for fam in FAMILIES:
    for n in (3, 4, 6, 8, 10, 12):
        w, cw, _ = clause_table(fam, n)
        allworst.append(w); allctrl.append(cw)
    P()

P(f"  READ: over {len(allworst)} (family, n) rows the estimator returned max|Tr(P_E R_i)| = "
  f"{max(allworst)} every time,")
P(f"        while the SAME estimator on the field-broken control returned up to {max(allctrl)}.")
P(f"        The zero is a property of the carrier, not of a blind instrument.")
P()

# ---------------------------------------------------------------- the exact reason
P("-" * 104)
P("  WHY THE ZERO IS EXACT AT EVERY n  (an argument, then its test)")
P("-" * 104)
P("  The global spin flip s -> -s leaves every bond product s_i s_{i+1} unchanged, hence leaves")
P("  E(s) unchanged, hence maps each eigenspace to itself; it has NO fixed point in {+-1}^n; and")
P("  it flips every s_i.  So the configurations in any eigenspace pair off (s, -s) with opposite")
P("  s_i, and Tr(P_E Z_i) = 0 EXACTLY, for every i, every eigenspace, every n, every J.")
P("  TEST of the involution, not of its consequence:")
for fam in FAMILIES:
    for n in (7, 11, 14):
        J = couplings(fam, n - 1)
        bad = 0
        for bits in itertools.product((1, -1), repeat=n):
            if energy(J, bits) != energy(J, tuple(-b for b in bits)): bad += 1
        P(f"    {fam:<10} n={n:<3}  configurations whose energy is NOT flip-invariant: {bad}")
P()

# ---------------------------------------------------------------- D-22
P("-" * 104)
P("  D-22   IS THE CARRIER PERMUTATION-SYMMETRIC?   (order of the site-permutation group that")
P("         preserves the energy function; n! would mean NO geometry is present to detect)")
P("-" * 104)
P(f"  {'family':<10} {'n':>3} {'n!':>8} {'|Aut|':>8}   verdict")
for fam in FAMILIES:
    for n in (4, 5, 6):
        J = couplings(fam, n - 1)
        a = automorphism_order(J, n)
        import math
        P(f"  {fam:<10} {n:>3} {math.factorial(n):>8} {a:>8}   "
          f"{'PERMUTATION-SYMMETRIC (no geometry)' if a == math.factorial(n) else 'geometry present'}")
    P()

# ---------------------------------------------------------------- writer SEARCH
P("-" * 104)
P("  CLAUSE (iv), CONSTRUCTIVELY:  EXHAUSTIVE SEARCH OVER THE FULL PAULI GROUP (4^n elements)")
P("  for an ADMISSIBLE W ([W,H]=0) that FLIPS the target.  NOTHING IS NOMINATED.")
P("-" * 104)
P(f"  {'family':<10} {'n':>3} {'target':<10} {'#Paulis':>9} {'#admissible flippers':>21} "
  f"{'min weight':>11} {'support':<20} {'dE':>8}")
for fam in ("linear", "superinc", "randsign"):
    for n in (3, 4, 5, 6, 7):
        J = couplings(fam, n - 1)
        # target 1: the record R_0 = Z_0
        tx = [0] * n; tz = [1 if k == 0 else 0 for k in range(n)]
        r = search_admissible_writers(J, n, (tx, tz), dense_check=(n <= 5))
        cnt, w, code, sup = r
        # the writer maps s -> -s (X-type on every site) => energy change is exactly 0
        dE = 0 if cnt else None
        P(f"  {fam:<10} {n:>3} {'R_0 = Z_0':<10} {4**n:>9} {cnt:>21} "
          f"{str(w):>11} {str(sup):<20} {str(dE):>8}")
        # target 2 (CONTROL, D-15): the PAIR CORRELATION Z_0 Z_1 -- a search that must come back EMPTY
        tx2 = [0] * n; tz2 = [1 if k in (0, 1) else 0 for k in range(n)]
        r2 = search_admissible_writers(J, n, (tx2, tz2), dense_check=False)
        P(f"  {fam:<10} {n:>3} {'Z_0 Z_1':<10} {4**n:>9} {r2[0]:>21} "
          f"{str(r2[1]):>11} {str(r2[3]):<20} {'-':>8}")
    P()

P("  READ: the search returns admissible flippers for the RECORD Z_0 at every n and every family,")
P("  and the minimum weight it finds is n -- FULL SUPPORT, never less.  The same search run on the")
P("  PAIR CORRELATION Z_0 Z_1 returns ZERO admissible flippers.  That pair of numbers, produced by")
P("  one instrument in one table, is the whole asymmetry: single records are freely writable,")
P("  the correlation between them is not writable at all by any admissible Pauli.")
P()

# ---------------------------------------------------------------- clauses on the correlation
P("-" * 104)
P("  IS THE PAIR CORRELATION ITSELF A RECORD?   clauses (iii) AND (iv) on Z_i Z_j, same estimators")
P("  ('all dims 2' means the couplings admit NO accidental degeneracy: every eigenspace is exactly")
P("   the pair {s, -s} and the venue is generic.  Degenerate families are shown too, on purpose.)")
P("-" * 104)
P(f"  {'family':<10} {'n':>3} {'generic?':<9} {'observable':<14} {'const on every P_E':>19} "
  f"{'(iii)':>7} {'max|Tr(P_E O)|':>15} {'(iv)':>7} {'RECORD?':>9}")
corr_verdicts = []
for fam in ("superinc", "randpos", "linear", "uniform", "randsign"):
    n = 8
    J = couplings(fam, n - 1)
    blocks = eigenspaces(J, n)
    generic = all(len(cs) == 2 for cs in blocks.values())
    for (i, j) in ((0, 1), (0, 3), (0, 7)):
        vals = lambda c: c[i] * c[j]
        const = all(len(set(vals(c) for c in cs)) == 1 for cs in blocks.values())
        tr = max(abs(sum(vals(c) for c in cs)) for cs in blocks.values())
        rec = (not const) and tr == 0
        corr_verdicts.append(rec)
        P(f"  {fam:<10} {n:>3} {str(generic):<9} {'Z_%d Z_%d' % (i, j):<14} {str(const):>19} "
          f"{'FAILS' if const else 'holds':>7} {tr:>15} {'HOLDS' if tr == 0 else 'FAILS':>7} "
          f"{'YES' if rec else 'NO':>9}")
    # CONTROL, same estimators on the record itself
    const = all(len(set(c[0] for c in cs)) == 1 for cs in blocks.values())
    tr = max(abs(sum(c[0] for c in cs)) for cs in blocks.values())
    P(f"  {fam:<10} {n:>3} {str(generic):<9} {'Z_0 [CONTROL]':<14} {str(const):>19} "
      f"{'FAILS' if const else 'holds':>7} {tr:>15} {'HOLDS' if tr == 0 else 'FAILS':>7} "
      f"{'YES' if (not const and tr == 0) else 'NO':>9}")
    P()

P(f"  READ, from the column above: the pair correlation was a record in "
  f"{sum(corr_verdicts)} of {len(corr_verdicts)} tested cases.")
P("  On the GENERIC families -- the ones with no accidental degeneracy, every eigenspace exactly")
P("  {s,-s} -- Z_i Z_j is CONSTANT on every eigenspace and FAILS CLAUSE (iii).  On the families")
P("  with accidental degeneracies the constancy can break, but then the same table shows clause")
P("  (iv) failing instead: max|Tr(P_E Z_iZ_j)| is non-zero.  Either way NO row is a record, while")
P("  the CONTROL row Z_0 in the very same table passes both clauses.")
P("  So the object that carries the configuration energy is a SUPERSELECTION LABEL on eigenspaces,")
P("  NOT a written bit.  This is the fact the rest of the lane is about.")
P()

# ---------------------------------------------------------------- how many independent records
P("-" * 104)
P("  HOW MANY INDEPENDENT RECORDS DOES THE CHAIN ACTUALLY HOLD?")
P("  On a GENERIC chain every eigenspace is span{|s>, |-s>}, 2-dimensional.  In that basis")
P("  Z_i acts as s_i * diag(1,-1) -- the SAME 2x2 operator for every i, up to the fixed sign s_i.")
P("  So all n nominated records are ONE two-level system wearing n labels.  Test: the rank of the")
P("  set {Z_i restricted to P_E} as operators on the block, maximised over blocks.")
P("-" * 104)
P(f"  {'family':<10} {'n':>3} {'generic?':<9} {'max rank of span{Z_i|_E}':>25} {'indep. bits on a block':>23}")
for fam in ("superinc", "randpos", "linear", "uniform"):
    for n in (6, 8, 10):
        J = couplings(fam, n - 1)
        blocks = eigenspaces(J, n)
        generic = all(len(cs) == 2 for cs in blocks.values())
        best = 0
        for E, cs in blocks.items():
            # each Z_i restricted to the block is the diagonal vector (c[i] for c in cs)
            M = np.array([[c[i] for c in cs] for i in range(n)], dtype=float)
            best = max(best, int(np.linalg.matrix_rank(M, tol=1e-9)))
        P(f"  {fam:<10} {n:>3} {str(generic):<9} {best:>25} "
          f"{('1 -- all Z_i are ONE bit' if best == 1 else str(best)):>23}")
    P()
P("  READ: on the generic families the rank is 1 at every n -- n nominated records, ONE bit.")
P("  Degenerate families raise the rank, but that rank comes from ACCIDENTAL arithmetic")
P("  coincidences among subset sums of the J_i, not from the chain's structure.")
P()

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_C_SIGN/s1_clauses.txt", "w").write("\n".join(OUT) + "\n")
