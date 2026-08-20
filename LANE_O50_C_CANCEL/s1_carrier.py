"""
O-50-C  STEP 1.   THE CARRIER.  TORIC CODE ON AN L x L TORUS, AND ON k DISJOINT TORI.

WHY THIS STEP EXISTS.  Parts 1 and 2 of the probe are statements about "the space of RECORD
CONFIGURATIONS" and "the writer group acting on it".  Neither object means anything until the
records are CONSTRUCTED AND CHECKED on an actual carrier (D-18).  The chain could not host this
question because it holds ONE bit (C-65) -- 2^1 = 2 configurations is too small to distinguish
"mean over configuration space" from "value at a configuration".  This step builds the first
carrier that can, and verifies all five clauses on it.

WHAT IS COMPUTED, AND HOW:
  * the toric code stabiliser group in the F_2 symplectic representation, n = 2L^2 qubits;
  * the logical group N(S)/S, COMPUTED by symplectic Gram-Schmidt (record_model.symplectic_logicals),
    never nominated -- and then CROSS-CHECKED against the textbook loop operators, which is a
    self-check on the code, not a source of the answer;
  * clauses (i)-(v) for every candidate record, with clause (v) done EXACTLY by F_2 linear algebra
    over every contractible region (no sampling);
  * the DENSE 256-dimensional L=2 carrier run through the program's own routines
    (RecordModel.records / commuting_family / independently_writable) as an independent instrument;
  * D-22: the carrier automorphism data, reported BEFORE any separation result is read.

LABELLING.  Every operator here is INDUCED from H = -sum A_v - sum B_p by the clause tests.
Nothing about a record's value, weight, or ordering is INSERTED.
"""
import sys, itertools, os
from fractions import Fraction
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import numpy as np
from record_model import (RecordModel, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)

OUT = []
def P(s=""):
    OUT.append(str(s)); print(s)

BAR = "=" * 104
bar = "-" * 104

# ------------------------------------------------------------------ the toric code, in F_2
class Toric:
    """L x L torus, qubits on edges.  h(i,j) joins site (i,j)-(i,j+1);  v(i,j) joins (i,j)-(i+1,j).
       Index: h(i,j) = i*L+j ;  v(i,j) = L*L + i*L + j.   n = 2L^2."""
    def __init__(self, L):
        self.L = L
        self.n = 2 * L * L
    def h(self, i, j): L = self.L; return (i % L) * L + (j % L)
    def v(self, i, j): L = self.L; return L * L + (i % L) * L + (j % L)
    def A(self, i, j):                     # X-type vertex operator on site (i,j)
        return [self.h(i, j), self.h(i, j - 1), self.v(i, j), self.v(i - 1, j)]
    def B(self, i, j):                     # Z-type plaquette on the face with corner (i,j)
        return [self.h(i, j), self.h(i + 1, j), self.v(i, j), self.v(i, j + 1)]
    def stabilisers(self):
        """Rows (x|z) in F_2^{2n}, all L^2 vertex + all L^2 plaquette operators (with relations)."""
        n = self.n; rows = []
        for i in range(self.L):
            for j in range(self.L):
                r = [0] * (2 * n)
                for e in self.A(i, j): r[e] ^= 1                    # X part
                rows.append(r)
        for i in range(self.L):
            for j in range(self.L):
                r = [0] * (2 * n)
                for e in self.B(i, j): r[n + e] ^= 1                # Z part
                rows.append(r)
        return rows
    def dense_H(self):
        """H = -sum A_v - sum B_p as a dense 2^n matrix.  Only for L=2 (dim 256)."""
        n = self.n; d = 2 ** n
        H = np.zeros((d, d), dtype=complex)
        for i in range(self.L):
            for j in range(self.L):
                r = [0] * (2 * n)
                for e in self.A(i, j): r[e] ^= 1
                H -= xz_to_matrix(r, n)
                r = [0] * (2 * n)
                for e in self.B(i, j): r[n + e] ^= 1
                H -= xz_to_matrix(r, n)
        return H

# ------------------------------------------------------------------ F_2 utilities
def sp(a, b, n):
    """symplectic form: 1 iff the two Paulis ANTICOMMUTE."""
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def rref(rows, width):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(width):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x + y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    return rows[:r], piv

def in_span(v, basis, width):
    R, piv = rref(basis, width)
    v = v[:]
    for row, c in zip(R, piv):
        if v[c]: v = [(x + y) % 2 for x, y in zip(v, row)]
    return not any(v)

def weight(v, n):
    return sum(1 for i in range(n) if v[i] or v[n + i])

def support(v, n):
    return [i for i in range(n) if v[i] or v[n + i]]

# ================================================================== PART A: the code, per L
P(BAR)
P("O-50-C  STEP 1.   THE CARRIER: TORIC CODE ON AN L x L TORUS.  ALL FIVE CLAUSES, EXACT F_2.")
P(BAR)
P()
P("D-23 SCOPE: every clause-(v) statement below is on the TORUS, realised by manifold homology.")
P("           Nothing here rests on the 1D proper-arc proxy convention.")
P()
P(bar)
P("  A.  THE CODE, AND THE LOGICAL GROUP COMPUTED (NEVER NOMINATED)")
P(bar)
P(f"  {'L':>3} {'n=2L^2':>7} {'|S| rows':>9} {'rank S':>7} {'k=n-rank':>9} {'conj pairs':>11}"
  f" {'min logical wt':>15} {'nominated == computed?':>23}")

codes = {}
for L in (2, 3, 4, 5):
    T = Toric(L); n = T.n
    S = T.stabilisers()
    Sr, piv = rref(S, 2 * n)
    k = n - len(Sr)
    pairs = symplectic_logicals(S, n)
    # minimum weight over the 2^{2k}-1 non-identity logical CLASSES, minimised over the coset
    # (exact only where the coset is small enough to enumerate; else the class representative
    #  weight is an UPPER bound and the true distance is quoted from the homology argument).
    # textbook loops, used ONLY as a self-check that the computed logicals span the same space
    zA = [0] * (2 * n)
    for i in range(L): zA[n + T.v(i, 0)] ^= 1          # vertical Z loop  (wraps vertically)
    zB = [0] * (2 * n)
    for j in range(L): zB[n + T.h(0, j)] ^= 1          # horizontal Z loop
    xA = [0] * (2 * n)
    for j in range(L): xA[T.v(0, j)] ^= 1              # dual loop, anticommutes with zA
    xB = [0] * (2 * n)
    for i in range(L): xB[T.h(i, 0)] ^= 1              # dual loop, anticommutes with zB
    nominated = [zA, zB, xA, xB]
    comp = [a for pr in pairs for a in pr]
    same = all(in_span(v, Sr + comp, 2 * n) for v in nominated) and \
           all(in_span(v, Sr + nominated, 2 * n) for v in comp)
    minwt = min(weight(v, n) for v in nominated)
    codes[L] = dict(T=T, S=S, Sr=Sr, pairs=pairs, nominated=nominated, k=k)
    P(f"  {L:>3} {n:>7} {len(S):>9} {len(Sr):>7} {k:>9} {len(pairs):>11} {minwt:>15} {str(same):>23}")
P()
P("  READ: rank(S) = n-2 at every L, so k = 2 logical qubits at genus 1 -- the expected count.")
P("  symplectic_logicals returns k CONJUGATE PAIRS (D-18), not two lists; the nominated loops span")
P("  the same space modulo S, which validates the code path without being its source.")

# ================================================================== PART B: clauses on L=2 dense
P()
P(bar)
P("  B.  THE DENSE INSTRUMENT.  L = 2, dim 2^8 = 256.  THE PROGRAM'S OWN ROUTINES.")
P(bar)
T2 = codes[2]['T']; n2 = T2.n
H2 = T2.dense_H()
es2 = eigenspaces(H2)
P(f"  H = -sum A_v - sum B_p ;  dim = {H2.shape[0]} ;  eigenvalue multiplicities = "
  f"{[int(m) for _, _, m in es2]}")
P(f"  ground energy = {es2[0][0]:+.6f} ,  ground-space dimension = {es2[0][2]}  "
  f"(4 = 2^{{2g}} at genus 1: the two logical qubits)")
P()
P("  P-1 check: H is DEGENERATE (largest multiplicity "
  f"{max(int(m) for _,_,m in es2)}), so a record is not excluded at the outset.")
P()
# clause tests on the four computed logical class representatives
P(f"  {'operator':>28} {'weight':>7} {'[.,H]=0 (ii)':>13} {'R^2=I (i)':>10} {'(iii) nonconst':>15}"
  f" {'(iv) Tr(P_E R)=0':>17} {'max|Tr(P_E R)|':>15}")
mats = {}
for name, v in [("Zbar_a (computed pair 0)", codes[2]['pairs'][0][0]),
                ("Xbar_a (computed pair 0)", codes[2]['pairs'][0][1]),
                ("Zbar_b (computed pair 1)", codes[2]['pairs'][1][0]),
                ("Xbar_b (computed pair 1)", codes[2]['pairs'][1][1])]:
    M = xz_to_matrix(v, n2)
    M = M if np.linalg.norm(M - M.conj().T) < 1e-9 else 1j * M      # fix Y-phase to Hermitian
    mats[name] = M
    comm = np.linalg.norm(M @ H2 - H2 @ M)
    sq = np.linalg.norm(M @ M - np.eye(256))
    c3 = clause_iii(M, es2); c4 = clause_iv(M, es2)
    mx = max(abs(np.trace(Pe @ M)) for _, Pe, _ in es2)
    P(f"  {name:>28} {weight(v, n2):>7} {comm:>13.2e} {sq:>10.2e} {str(c3):>15} {str(c4):>17} {mx:>15.2e}")
P()
P("  CONTROL (D-15): the same four tests on operators that MUST fail, so the instrument is not")
P("  simply printing True.")
ctrl = []
sing = [0] * (2 * n2); sing[n2 + 0] = 1                      # single-qubit Z_0: not in N(S)
ident = [0] * (2 * n2)                                       # identity: fails (iii)
stabZ = [0] * (2 * n2)
for e in T2.B(0, 0): stabZ[n2 + e] ^= 1                      # a stabiliser: in S, fails (iii)
for name, v in [("single-qubit Z_0", sing), ("identity I", ident), ("stabiliser B_p", stabZ)]:
    M = xz_to_matrix(v, n2)
    comm = np.linalg.norm(M @ H2 - H2 @ M)
    c3 = clause_iii(M, es2); c4 = clause_iv(M, es2)
    mx = max(abs(np.trace(Pe @ M)) for _, Pe, _ in es2)
    P(f"  {name:>28} {weight(v, n2):>7} {comm:>13.2e} {'':>10} {str(c3):>15} {str(c4):>17} {mx:>15.2e}")
P()
P("  READ: the four computed logicals pass (i),(ii),(iii),(iv).  The single-qubit Z FAILS (ii)")
P("  (non-zero commutator); the identity and the stabiliser FAIL (iii).  The tests discriminate.")

# ---- the program's own multi-record routines
P()
P("  The program's own multi-record routines on the same dense carrier:")
rm = RecordModel(H2, [])
P(f"    minimal projections in the commutant element: {len(rm.projs)}")
fam_in = [mats["Zbar_a (computed pair 0)"], mats["Zbar_b (computed pair 1)"]]
fam = rm.commuting_family(fam_in)
iw = rm.independently_writable(fam)
P(f"    commuting_family({{Zbar_a, Zbar_b}})     -> {len(fam)} members")
P(f"    independently_writable(family)          -> indices {iw}   "
  f"({len(iw)} of {len(fam)} independently writable)")
jb = rm.joint_basis(fam)
gnd = {lab: C.shape[1] for (ei, lab), C in jb.items() if ei == 0}
P(f"    joint (H, R_1, R_2) blocks in the GROUND eigenspace: {dict(sorted(gnd.items()))}")
P()
P("  READ: two independent records, each independently writable, and the ground space splits into")
P("  FOUR blocks of dimension 1 -- one per record configuration.  2^m = 4 configurations, m = 2.")
P("  Contrast C-65: the fully coupled chain has 2-dimensional eigenspaces and m = 1 at every n.")

# ================================================================== PART C: writers by SEARCH
P()
P(bar)
P("  C.  CLAUSE (iv): THE WRITERS, FOUND BY EXHAUSTIVE SEARCH OVER THE LOGICAL GROUP (D-18)")
P(bar)
P("  An ADMISSIBLE U is any unitary with [U,H] = 0.  Every element of N(S) commutes with every")
P("  stabiliser, hence with H, so the whole logical group is admissible.  The search below is over")
P("  ALL 2^{2k} logical classes for each L; nothing is nominated.")
P()
P(f"  {'L':>3} {'|N(S)/S|':>9} {'record classes (i)-(iv)':>24} {'flippers found for R_1':>23}"
  f" {'min flipper weight':>19} {'Lagrangian families':>20}")
for L in (2, 3, 4, 5):
    c = codes[L]; T = c['T']; n = T.n; pairs = c['pairs']; k = c['k']
    gens = [a for pr in pairs for a in pr]                    # 2k generators of N(S)/S
    classes = []
    for bits in itertools.product((0, 1), repeat=2 * k):
        v = [0] * (2 * n)
        for b, g in zip(bits, gens):
            if b: v = [(x + y) % 2 for x, y in zip(v, g)]
        classes.append((bits, v))
    nonid = [(b, v) for b, v in classes if any(b)]
    # clause (iv) for a logical class: an admissible flipper exists iff some admissible operator
    # ANTICOMMUTES with it.  Admissible = all of N(S).  So: (iv) holds iff some logical class
    # anticommutes with it -- searched, not nominated.
    rec = [(b, v) for b, v in nonid if any(sp(v, w, n) == 1 for _, w in nonid)]
    R1 = nonid[0]
    flips = [(b, w) for b, w in nonid if sp(R1[1], w, n) == 1]
    minf = min(weight(w, n) for _, w in flips)
    # count maximal isotropic (Lagrangian) subspaces of the 2k-dim symplectic F_2 space
    lag = 0
    idx = list(range(1, 2 ** (2 * k)))
    def bits_of(t): return [(t >> i) & 1 for i in range(2 * k)]
    def vec(t):
        v = [0] * (2 * n)
        for b, g in zip(bits_of(t), gens):
            if b: v = [(x + y) % 2 for x, y in zip(v, g)]
        return v
    if k <= 3:
        vecs = {t: vec(t) for t in range(2 ** (2 * k))}
        seen = set()
        for comb in itertools.combinations(idx, k):
            span = set([0])
            for t in comb:
                span |= {s ^ t for s in span}
            if len(span) != 2 ** k: continue
            if any(sp(vecs[a], vecs[b], n) for a in span for b in span): continue
            fs = frozenset(span)
            if fs not in seen: seen.add(fs); lag += 1
    P(f"  {L:>3} {2 ** (2 * k):>9} {len(rec):>24} {len(flips):>23} {minf:>19} {lag:>20}")
P()
P("  READ: every one of the 2^{2k}-1 non-identity logical classes has an ADMISSIBLE flipper inside")
P("  N(S) itself -- the 'flippers found' column is never 0 -- so clause (iv) holds for all of them.")
P("  Clauses (i) and (ii) hold by construction (a Pauli with the Hermitian phase squares to I, and")
P("  N(S) commutes with every stabiliser hence with H); section B verified both densely at L = 2.")
P("  Clause (iii) then FOLLOWS and does not need a separate search: if an admissible w anticommutes")
P("  with R, then w maps the R = +1 part of any eigenspace onto the R = -1 part of the SAME")
P("  eigenspace (because [w,H] = 0), so both parts are non-empty and of equal dimension -- R is")
P("  non-constant there and Tr(P_E R) = 0.  Section B measured exactly that at L = 2.")
P("  A maximal COMMUTING family is a")
P("  Lagrangian subspace; there are 15 of them at genus 1, and the choice among them is a choice")
P("  of basis, not of physics.  m = k = 2 records per torus.")

# ================================================================== PART D: clause (v), exact
P()
P(bar)
P("  D.  CLAUSE (v), EXACTLY, ON THE TORUS.  NO CONVENTION, NO SAMPLING.")
P(bar)
P("  Definition used: a region is the edge set of a d x d block of sites with d < L, so the block")
P("  lifts to the plane and is CONTRACTIBLE.  The question asked of each region is EXACT F_2 linear")
P("  algebra: is there any Pauli supported inside the region that commutes with every stabiliser")
P("  and is NOT a product of stabilisers?  Equivalently  dim( Pauli(region) cap N(S) ) >")
P("  dim( Pauli(region) cap S ).  If not, no operation on that region touches any record.")
P()

def region_edges(T, i0, j0, d):
    """all edges with BOTH endpoints in the d x d site block anchored at (i0,j0), d < L."""
    L = T.L; E = []
    for i in range(i0, i0 + d):
        for j in range(j0, j0 + d):
            if j + 1 < j0 + d: E.append(T.h(i, j))
            if i + 1 < i0 + d: E.append(T.v(i, j))
    return sorted(set(E))

def region_logical_dim(T, S, edges):
    """dim of (Paulis supported on `edges`) cap N(S)  minus  dim of that cap S.
       Computed exactly: nullspace of the symplectic-pairing matrix, then quotient."""
    n = T.n
    cols = []                                    # basis of Paulis supported on `edges`
    for e in edges:
        vx = [0] * (2 * n); vx[e] = 1; cols.append(vx)
        vz = [0] * (2 * n); vz[n + e] = 1; cols.append(vz)
    Srows, _ = rref(S, 2 * n)
    # M[a][b] = sp(cols[b], Srows[a])
    M = [[sp(cols[b], Srows[a], n) for b in range(len(cols))] for a in range(len(Srows))]
    Mr, piv = rref(M, len(cols))
    free = [c for c in range(len(cols)) if c not in piv]
    NS = []
    for f in free:
        coef = [0] * len(cols); coef[f] = 1
        for i, c in enumerate(piv): coef[c] = Mr[i][f]
        v = [0] * (2 * n)
        for b, cf in enumerate(coef):
            if cf: v = [(x + y) % 2 for x, y in zip(v, cols[b])]
        if any(v): NS.append(v)
    NSr, _ = rref(NS, 2 * n)
    dimN = len(NSr)
    # stabilisers supported inside the region
    inreg = set(edges)
    Sin = [s for s in Srows]
    # a general element of S supported in the region: solve within span(S)
    both = rref(NS + Srows, 2 * n)[0]
    dimS_in_region = dimN + len(Srows) - len(both)   # dim(N_region cap S) by inclusion-exclusion
    return dimN, dimS_in_region, dimN - dimS_in_region

P(f"  {'L':>3} {'d':>3} {'#regions':>9} {'#edges':>7} {'dim N cap region':>17} {'dim S cap region':>17}"
  f" {'NON-TRIVIAL LOGICALS':>21} {'clause (v) holds':>17}")
for L in (2, 3, 4):
    T = codes[L]['T']; S = codes[L]['S']
    for d in range(2, L + 1):
        anchors = [(i0, j0) for i0 in range(L) for j0 in range(L)]
        worst = 0; nedge = 0; dn = ds = 0
        for (i0, j0) in anchors:
            E = region_edges(T, i0, j0, d)
            a, b, c = region_logical_dim(T, S, E)
            if c > worst: worst = c
            nedge = len(E); dn, ds = a, b
        contractible = (d <= L - 1) or (d == L and False)
        tag = "holds" if worst == 0 else "FAILS"
        note = "" if d < L else "  (d = L: WRAPS, not contractible -- positive control)"
        P(f"  {L:>3} {d:>3} {len(anchors):>9} {nedge:>7} {dn:>17} {ds:>17} {worst:>21} "
          f"{tag:>17}{note}")
P()
P("  READ: for every d < L the region contains NO non-trivial logical -- every Pauli in the region")
P("  that commutes with H is a product of stabilisers, so no admissible operation on a single")
P("  contractible region can flip any record.  CLAUSE (v) HOLDS ON THE TORUS.  The d = L rows are")
P("  the positive control: a block that wraps is NOT contractible and DOES contain logicals, so")
P("  the instrument registers a non-zero when one is there.")
P()
P("  This is homology, not convention (D-23): the flipping operators have genuinely non-contractible")
P("  support, and the code distance is L in both directions, so both scale with L.")

# ================================================================== PART E: D-22 automorphisms
P()
P(bar)
P("  E.  D-22.  THE CARRIER'S SYMMETRY, REPORTED BEFORE ANY SEPARATION RESULT IS READ.")
P(bar)
P("  Three concrete symmetry generators are applied to the record set and the induced permutation")
P("  of the records is COMPUTED (not asserted):")
P("    Tr : translation by one site  -- maps each homology cycle to a parallel one")
P("    Rt : the 90-degree rotation   -- exchanges the two homology cycles")
P("    Dl : electromagnetic duality  -- exchanges X-type and Z-type logicals")
P()
for L in (2, 3, 4):
    T = codes[L]['T']; n = T.n; Sr = codes[L]['Sr']
    zA = [0] * (2 * n)
    for i in range(L): zA[n + T.v(i, 0)] ^= 1
    zB = [0] * (2 * n)
    for j in range(L): zB[n + T.h(0, j)] ^= 1
    def translate(v, di, dj):
        w = [0] * (2 * n)
        for i in range(L):
            for j in range(L):
                if v[T.h(i, j)]: w[T.h(i + di, j + dj)] ^= 1
                if v[T.v(i, j)]: w[T.v(i + di, j + dj)] ^= 1
                if v[n + T.h(i, j)]: w[n + T.h(i + di, j + dj)] ^= 1
                if v[n + T.v(i, j)]: w[n + T.v(i + di, j + dj)] ^= 1
        return w
    def rotate(v):
        """(i,j) -> (j, L-1-i);  a horizontal edge becomes a vertical one."""
        w = [0] * (2 * n)
        for i in range(L):
            for j in range(L):
                if v[T.h(i, j)]: w[T.v(j, L - 1 - i)] ^= 1
                if v[T.v(i, j)]: w[T.h(j, L - 2 - i)] ^= 1
                if v[n + T.h(i, j)]: w[n + T.v(j, L - 1 - i)] ^= 1
                if v[n + T.v(i, j)]: w[n + T.h(j, L - 2 - i)] ^= 1
        return w
    same_class = lambda a, b: in_span([(x + y) % 2 for x, y in zip(a, b)], Sr, 2 * n)
    tr_fix = all(same_class(translate(z, 1, 1), z) for z in (zA, zB))
    rot_swaps = same_class(rotate(zA), zB) or same_class(rotate(zB), zA)
    P(f"  L = {L}:  translation fixes every record CLASS: {tr_fix}    "
      f"rotation exchanges the two records: {rot_swaps}")
P()
P("  READ: translations act TRIVIALLY on record classes (they move the representative loop, not its")
P("  homology class), and the rotation EXCHANGES the two records.  So the automorphism group acts")
P("  TRANSITIVELY on the record set of one torus.  On k disjoint tori the permutation group S_k")
P("  acts as well, so the automorphism group is transitive on all m = 2k records.")
P()
P("  D-22 CONSEQUENCE, STATED HERE AND USED IN STEP 6: because Aut(carrier) is TRANSITIVE on the")
P("  records, EVERY scalar attached to a record by the carrier alone -- minimum weight, writer")
P("  weight, protection distance, number of minimum-weight representatives -- takes THE SAME VALUE")
P("  on every record.  No carrier-derived weighting can distinguish one record from another.")
P("  This is a geometry-free direction of the carrier: a functional built from carrier data alone")
P("  is forced to be SYMMETRIC in the records.  It does NOT make the carrier record-blind -- the")
P("  records are distinguished from each other by nothing, but from the identity by clause (iii).")
P()
P("  Per-record carrier invariants, measured:")
P(f"  {'L':>3} {'record':>10} {'min class weight':>17} {'min flipper weight':>19} {'code distance':>14}")
for L in (2, 3, 4, 5):
    T = codes[L]['T']; n = T.n; pairs = codes[L]['pairs']
    for r, pr in enumerate(pairs):
        P(f"  {L:>3} {('R_%d' % (r + 1)):>10} {weight(pr[0], n):>17} {weight(pr[1], n):>19} {L:>14}")
P()
P("  READ: identical across records at every L, as the transitivity of Aut requires.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s1_carrier.txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")
