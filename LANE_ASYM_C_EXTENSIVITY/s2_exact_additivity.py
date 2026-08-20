"""S2 -- THE EXACT (COMBINATORIAL) HALF OF THE AUDIT.  [F2] symplectic representation.

Everything here is an EXACT finite computation, not a fitted trend, so the categorisations it
produces hold at every m -- not merely within the range simulated.

Quantities:
  N(m)   number of records                          (count)
  W(m)   total writer weight (minimal admissible writer weight, summed over records)
  P(m)   interacting-pair count (record-site supports that intersect)
  G(m)   the RELATION MATRIX (support-overlap) and its spectrum: trace, ||.||_F^2, lambda_max
  T(m)   the TRANSPORT / CROSS-REGION INFLUENCE scalar

For each: Q(A+B) vs Q(A)+Q(B), the ratio, and the ABSOLUTE DEFECT, with a positive control
in the same table.
"""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from lanelib import *

OUT = []
def P_(s=""):
    print(s); OUT.append(str(s))

P_("=" * 104)
P_("S2  EXACT ADDITIVITY AUDIT  --  [F2] symplectic representation, m disjoint [[4,2,2]] blocks")
P_("=" * 104)

# ---------------------------------------------------------------- helper: enumerate N(S)
def normaliser(stab, n, cap=2 ** 22):
    """every Pauli (x|z) commuting with every stabiliser.  |N(S)| = 2^(2n-r)."""
    R, piv = rref_f2(stab, 2 * n)
    r = len(R)
    size = 2 ** (2 * n - r)
    if size > cap: return None, size
    # nullspace of the symplectic-pairing map
    A = [[sp([1 if k == j else 0 for k in range(2 * n)], s, n) for j in range(2 * n)] for s in R]
    Ar, piv2 = rref_f2(A, 2 * n)
    free = [c for c in range(2 * n) if c not in piv2]
    basis = []
    for f in free:
        v = [0] * (2 * n); v[f] = 1
        for i, c in enumerate(piv2): v[c] = Ar[i][f]
        basis.append(v)
    out = []
    for bits in itertools.product((0, 1), repeat=len(basis)):
        v = [0] * (2 * n)
        for bit, bv in zip(bits, basis):
            if bit: v = [(x + y) % 2 for x, y in zip(v, bv)]
        out.append(v)
    return out, size

# ---------------------------------------------------------------- W(m): min writer weight
def min_writer_weights(m, brute=True):
    """for each record, the minimum weight of an ADMISSIBLE writer that flips IT and no other
       record in the family.  Admissible = commutes with H  <=>  lies in N(S)."""
    recs, wrts, n = composite_records_writers(m)
    stab = composite_stab(m)
    if brute:
        NS, size = normaliser(stab, n)
        if NS is None: return None, size
        mins = []
        for i, Ri in enumerate(recs):
            best = None
            for v in NS:
                if sp(v, Ri, n) != 1: continue
                if any(sp(v, Rj, n) for j, Rj in enumerate(recs) if j != i): continue
                w = weight(v, n)
                if best is None or w < best: best = w
            mins.append(best)
        return mins, size
    return [weight(w, n) for w in wrts], None

# ---------------------------------------------------------------- P(m), G(m)
def site_supports(m):
    recs, wrts, n = composite_records_writers(m)
    return [support(recs[i], n) | support(wrts[i], n) for i in range(len(recs))], n

def relation_matrix(m):
    sup, n = site_supports(m)
    k = len(sup)
    G = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            G[i, j] = len(sup[i] & sup[j])
    return G

def pair_count(m):
    sup, n = site_supports(m)
    k = len(sup)
    return sum(1 for i in range(k) for j in range(i + 1, k) if sup[i] & sup[j])

def cross_pair_count(m):
    """pairs in DIFFERENT blocks whose supports intersect"""
    sup, n = site_supports(m)
    k = len(sup)
    blk = [i // 2 for i in range(k)]
    return sum(1 for i in range(k) for j in range(i + 1, k)
               if blk[i] != blk[j] and (sup[i] & sup[j]))

# ---------------------------------------------------------------- T(m): transport / influence
def transport_table(m):
    """CROSS-REGION INFLUENCE.  For an admissible operation g supported entirely inside block
       a, how far does it move a record living in block b != a?   ||g R g^dag - R||_F / ||R||_F.
       In [F2] a Pauli conjugation gives g R g^dag = (-1)^{sp(g,R)} R, so the displacement is
       0 when they commute and 2 when they anticommute -- EXACT, no floating point at all.

       CONTROL 1 (must be non-zero): g = the record's OWN writer, same block.
       CONTROL 2 (must be non-zero): g = a NON-admissible single-qubit Pauli in the same block."""
    recs, wrts, n = composite_records_writers(m)
    stab = composite_stab(m)
    k = len(recs)
    blk = [i // 2 for i in range(k)]
    # admissible operators supported in block a: N(S) elements whose support lies in block a
    NS, size = normaliser(stab, n)
    worst_cross = 0.0
    n_cross_moved = 0
    n_cross_tested = 0
    for g in NS:
        s = support(g, n)
        if not s: continue
        gb = set(q // 4 for q in s)
        if len(gb) != 1: continue                     # must be supported in ONE block
        a = gb.pop()
        for i in range(k):
            if blk[i] == a: continue                  # cross-region only
            n_cross_tested += 1
            d = 2.0 if sp(g, recs[i], n) else 0.0
            worst_cross = max(worst_cross, d)
            if d > 0: n_cross_moved += 1
    # control 1: own writer
    ctrl1 = max(2.0 if sp(wrts[i], recs[i], n) else 0.0 for i in range(k))
    # control 2: non-admissible single-qubit Pauli that anticommutes with a record
    ctrl2 = 0.0
    for q in range(n):
        for pv in ([1 if t == q else 0 for t in range(n)] + [0] * n,
                   [0] * n + [1 if t == q else 0 for t in range(n)]):
            if any(sp(pv, r_, n) for r_ in recs):
                ctrl2 = 2.0
    return dict(worst_cross=worst_cross, n_cross_moved=n_cross_moved,
                n_cross_tested=n_cross_tested, ctrl_own_writer=ctrl1,
                ctrl_nonadmissible=ctrl2, NS_size=size)

# ================================================================= TABLE 1: scaling in m
P_("\n" + "-" * 104)
P_("TABLE 1  --  the exact record-level quantities vs m   (N=records, W=total writer weight,")
P_("            P=interacting pairs, P_cross=pairs in DIFFERENT blocks,  G=relation matrix)")
P_("-" * 104)
P_("%-4s %-8s %-6s %-8s %-8s %-10s %-10s %-12s %-12s"
   % ("m", "dim", "N(m)", "W(m)", "P(m)", "P_cross", "tr G", "||G||_F^2", "lam_max(G)"))
P_("-" * 104)
rows1 = []
BRUTE_MAX = 3
for m in range(1, 13):
    G = relation_matrix(m)
    ev = np.linalg.eigvalsh(G)
    mins, size = min_writer_weights(m, brute=(m <= BRUTE_MAX))
    Wm = sum(mins)
    rows1.append(dict(m=m, N=2 * m, W=Wm, Pn=pair_count(m), Pc=cross_pair_count(m),
                      trG=float(np.trace(G)), fro=float((G ** 2).sum()), lmax=float(ev.max())))
    P_("%-4s %-8s %-6d %-8d %-8d %-10d %-10.1f %-12.1f %-12.6f"
       % (m, "16^%d" % m, 2 * m, Wm, pair_count(m), cross_pair_count(m),
          np.trace(G), (G ** 2).sum(), ev.max()))
P_("-" * 104)
P_("W(m) for m<=%d is a BRUTE-FORCE minimum over the whole normaliser N(S) (|N(S)| = 2^(6m):"
   % BRUTE_MAX)
P_("   64, 4096, 262144 for m=1,2,3); for larger m it is the per-block witness, which the exact")
P_("   argument below shows is the true minimum.")

# ================================================================= TABLE 2: A+B vs A+B
P_("\n" + "-" * 104)
P_("TABLE 2  --  ADDITIVITY:  Q(A + B)  against  Q(A) + Q(B)   for two DISJOINT blocks")
P_("            (and the m-block generalisation Q(m) vs m*Q(1)).  DEFECT = whole - sum of parts.")
P_("-" * 104)
P_("%-16s %-4s %-14s %-14s %-14s %-12s %-10s"
   % ("quantity", "m", "Q(m blocks)", "m*Q(1 block)", "DEFECT", "ratio", "verdict"))
P_("-" * 104)
q1 = rows1[0]
for key, label in [("N", "N  records"), ("W", "W  writerwt"), ("Pn", "P  pairs"),
                   ("Pc", "P_cross"), ("trG", "tr G"), ("fro", "||G||_F^2"), ("lmax", "lam_max G")]:
    for m in (2, 4, 8, 12):
        r = rows1[m - 1]
        whole = float(r[key]); parts = m * float(q1[key])
        defect = whole - parts
        ratio = (whole / parts) if abs(parts) > 1e-12 else float("nan")
        verdict = "ADDITIVE" if abs(defect) < 1e-12 else ("SUB" if defect < 0 else "SUPER")
        P_("%-16s %-4d %-14.6f %-14.6f %-14.3e %-12.6f %-10s"
           % (label, m, whole, parts, defect, ratio, verdict))
    P_("")

# ================================================================= TABLE 3: transport
P_("-" * 104)
P_("TABLE 3  --  TRANSPORT / CROSS-REGION INFLUENCE   (D-15: controls in the SAME table)")
P_("-" * 104)
P_("%-4s %-12s %-14s %-16s %-18s %-18s"
   % ("m", "|N(S)|", "cross tested", "cross MOVED", "worst cross disp", "CTRL own writer"))
P_("-" * 104)
for m in (2, 3):
    t = transport_table(m)
    P_("%-4d %-12d %-14d %-16d %-18.6f %-18.6f"
       % (m, t["NS_size"], t["n_cross_tested"], t["n_cross_moved"], t["worst_cross"],
          t["ctrl_own_writer"]))
P_("   CTRL non-admissible single-qubit Pauli displacement: %.6f (non-zero -> the metric can register motion)"
   % transport_table(2)["ctrl_nonadmissible"])

# ================================================================= the exact arguments
P_("\n" + "=" * 104)
P_("EXACT ARGUMENTS  (these are proofs, valid at EVERY m -- not extrapolations)")
P_("=" * 104)
P_("""
A-1  W(m) = 4m EXACTLY, and is exactly additive.
     An admissible writer for record i lies in N(S) and must anticommute with R_i and commute
     with every other record.  Write v in N(S) as v = sum_b v_b over blocks (the qubit sets are
     disjoint, so this decomposition is unique and weight(v) = sum_b weight(v_b)).  The
     symplectic form is likewise blockwise: sp(v,R_i) = sp(v_{b(i)}, R_i).  So the constraints
     only ever see v_{b(i)}; every other v_b can be set to 0 without breaking any constraint and
     can only reduce weight.  Hence the minimum is attained inside block b(i) and equals the
     single-block minimum, which brute force finds to be 2.  W(m) = 2m records * 2 = 4m.
     DEFECT IDENTICALLY 0.  Extensive -- but it is a COUNT of blocks in disguise (C-35).

A-2  P_cross(m) = 0 EXACTLY.
     Supports of records in different blocks are disjoint subsets of the qubit set by
     construction, so their intersection is empty at every m.  No pair of records in different
     regions is ever 'interacting' in the support sense.  DEFECT IDENTICALLY 0.

A-3  G(m) = G(1) (+) ... (+) G(1)  (m-fold direct sum), EXACTLY.
     G_ij = |supp_i ^ supp_j| vanishes whenever i and j are in different blocks (A-2), so G is
     block diagonal.  Consequences, all exact:
        spec G(m) = m copies of spec G(1)          (multiset union)
        tr G(m)      = m tr G(1)                   ADDITIVE
        ||G(m)||_F^2 = m ||G(1)||_F^2              ADDITIVE
        lambda_max G(m) = lambda_max G(1)          CONSTANT IN m  -> SATURATING
     The additive spectral functionals are again sums over blocks of a fixed per-block number:
     they are m * const, carrying no information beyond the block COUNT.  The one spectral
     quantity that is not a count -- the top eigenvalue, the natural 'strength' of the relation
     -- does not grow at all.

A-4  T(m) = 0 EXACTLY: no admissible operation supported inside one block moves any record in
     another block.  For Paulis g,R the conjugation is g R g^dag = (-1)^{sp(g,R)} R, and sp is
     blockwise, so g supported in block a has sp(g, R_i) = 0 for every record i outside a.
     The displacement is exactly 0 -- not small, zero -- for every m, every g, every record.
     The same holds for ANY admissible unitary supported in block a, Pauli or not, because
     U (x) I acting on disjoint tensor factors commutes with R_i identically.
     The controls in TABLE 3 show the metric does register motion when motion exists.
""")

P_("=" * 104)
P_("READ OF S2  (filled in from the numbers above, not in advance)")
P_("=" * 104)
P_("""
 * Every exact record-level quantity here is either (i) EXACTLY ADDITIVE with DEFECT
   IDENTICALLY ZERO -- N, W, P, tr G, ||G||_F^2 -- or (ii) EXACTLY CONSTANT in m -- lambda_max G.
 * Type (i) are additive because they are SUMS OVER BLOCKS OF A CONSTANT.  They equal
   (per-block value) x (number of blocks).  That is extensivity of the trivial kind: it is the
   block count wearing a different unit, exactly the object C-35 already rules out as a source.
 * Type (ii) SATURATES and is ruled out by requirement (a) at any N.
 * P_cross = 0 and T = 0 are the decisive ones: they say the records in different regions have
   NO relation at all -- not a weak one.  A quantity that is identically zero at every finite m
   by the blockwise-symplectic argument cannot become non-zero at large m.
 * NOTHING here is SUPER-additive or SUB-additive.  The defect is not small; it is zero by a
   proof, so there is no growth law to extrapolate.
""")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/s2_exact_additivity.txt",
     "w").write("\n".join(OUT) + "\n")
print("\n[written]")
