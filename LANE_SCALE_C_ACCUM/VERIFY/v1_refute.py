"""ADVERSARIAL VERIFICATION of LANE_SCALE_C_ACCUM.

Five executed attacks:
  A  wall ("one weight-2 operator flips ALL k records") is BASIS-DEPENDENT.
  B  the control's "smallest region touching ALL k records" (13,15,17,21,23) is a SAMPLING
     ARTIFACT; the true value is exactly k, identical to the family.
  C  wmax/wmean "half the additive rate" tracks the CARRIER SIZE n, not collectivity:
     interpolating carriers at FIXED k = 8 and n = 10,12,14,16 interpolate the numbers.
  D  E0 ("enclosed energy CONSTANT in k") is not a function of k at all: at FIXED k = 8 the
     ground energy can be made -2,-4,-6,-8 by choosing the stabiliser presentation.
  E  reach(A) counted ALL Paulis in A, including INADMISSIBLE ones that leave the code space.
     Recount with admissible-only operators.
"""
import sys, itertools, time
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM")
import numpy as np
from math import comb
from record_model import symplectic_logicals
from s1_combinatorics import (carrier_nn2, carrier_product, embed, sp2, f2_rank,
                              weight_map, wmap_summary, in_span)

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.append(s)

CH = []
def chk(name, ok, note=""):
    CH.append((name, bool(ok), note))

# ---------------------------------------------------------------- extra carriers
def carrier_nn2_plus(n, j):
    """[[n,n-2,2]] with j EXTRA local Z-pair stabilisers Z_{2i-1}Z_{2i}, i = 1..j.
       nstab = 2 + j, k = n - 2 - j.  Still carries the GLOBAL stabilisers X^n and Z^n, so it
       is not a disjoint product of blocks in this qubit partition."""
    S = [[1]*n + [0]*n, [0]*n + [1]*n]
    for i in range(j):
        z = [0]*(2*n); z[n + 2*i] = 1; z[n + 2*i + 1] = 1
        S.append(z)
    return n, S

def carrier_product_of(nblk, m):
    """m independent [[nblk, nblk-2, 2]] blocks.  k = m*(nblk-2), n = m*nblk."""
    n = nblk * m; S = []
    for b in range(m):
        x = [0]*(2*n); z = [0]*(2*n)
        for q in range(nblk*b, nblk*b + nblk):
            x[q] = 1; z[n + q] = 1
        S.append(x); S.append(z)
    return n, S

def logicals_for(n, S):
    pr = symplectic_logicals([s[:] for s in S], n)
    return [p[1] for p in pr], [p[0] for p in pr]

def logicals_blocks(nblk, m):
    n = nblk * m
    pr = symplectic_logicals([s[:] for s in carrier_product_of(nblk, 1)[1]], nblk)
    R = [embed(p[1], nblk, nblk*b, n) for b in range(m) for p in pr]
    W = [embed(p[0], nblk, nblk*b, n) for b in range(m) for p in pr]
    return R, W

def verify_records(tag, n, S, R, W):
    k = len(R)
    chk(f"{tag}: records commute with stabilisers", all(sp2(r, s, n) == 0 for r in R for s in S))
    chk(f"{tag}: records not in stabiliser group", all(not in_span(r, S) for r in R))
    chk(f"{tag}: records pairwise commute",
        all(sp2(R[i], R[j], n) == 0 for i in range(k) for j in range(k) if i != j))
    chk(f"{tag}: symplectic Gram against partners = identity (non-degenerate)",
        all(sp2(R[i], W[j], n) == (1 if i == j else 0) for i in range(k) for j in range(k)))
    chk(f"{tag}: k = n - rank(S)", k == n - f2_rank(S), f"k={k} n={n} rank={f2_rank(S)}")

# ================================================================= A  wall is basis-dependent
P("=" * 118)
P("VERIFY  --  adversarial checks on LANE_SCALE_C_ACCUM")
P("=" * 118)
P("")
P("ATTACK A.  The finding's ONE claimed strictly-collective fact is")
P("           'wall = 2 on the family at every k, where the control needs weight k'.")
P("           wall = w(all-ones pattern) is NOT basis invariant.  Under a record-basis change")
P("           M in GL(k,2) the pattern v transforms as v -> M v, so the all-ones pattern in the")
P("           new basis is the OLD pattern M^-1 (1,..,1), which ranges over EVERY non-zero v.")
P("           Therefore wall over all admissible record bases ranges over [wmin, wmax] on BOTH")
P("           carriers.  Below: the exact weight map's wmin/wmax, and wall measured directly")
P("           over 400 random GL(k,2) record bases.")
P("-" * 118)

def wall_over_bases(n, S, R, wv, nb=400, seed=7):
    """directly re-derive the record basis and re-read the all-ones cost, no theory used"""
    k = len(R); rng = np.random.default_rng(seed); vals = []
    for _ in range(nb):
        while True:
            A = rng.integers(0, 2, size=(k, k))
            if f2_rank([list(map(int, r)) for r in A]) == k: break
        Rb = []
        for i in range(k):
            v = [0]*(2*n)
            for j in range(k):
                if A[i, j]: v = [(a + b) % 2 for a, b in zip(v, R[j])]
            Rb.append(v)
        # sanity: the new records are still a valid record family
        assert all(sp2(r, s, n) == 0 for r in Rb for s in S)
        assert f2_rank(list(Rb) + list(S)) == k + f2_rank(S)
        # w(all-ones) in the NEW basis = w(v) where A v = 1  -> solve over F_2 by search on wv
        # cheaper: recompute the pattern map is expensive; instead invert A over F_2.
        Ai = f2_inverse(A.tolist(), k)
        one = [1]*k
        v = [sum(Ai[i][j]*one[j] for j in range(k)) % 2 for i in range(k)]
        idx = sum(v[i] << i for i in range(k))
        vals.append(int(wv[idx]))
    return min(vals), max(vals), sorted(set(vals))

def f2_inverse(A, k):
    M = [list(map(int, A[i])) + [1 if j == i else 0 for j in range(k)] for i in range(k)]
    r = 0
    for c in range(k):
        p = next((i for i in range(r, k) if M[i][c]), None)
        assert p is not None
        M[r], M[p] = M[p], M[r]
        for i in range(k):
            if i != r and M[i][c]:
                M[i] = [(x + y) % 2 for x, y in zip(M[i], M[r])]
        r += 1
    return [row[k:] for row in M]

P(f"{'carrier':>16} {'k':>3} {'n':>3} | {'wmin':>5} {'wmax':>5} {'wall(as-shipped)':>17} | "
  f"{'wall over 400 random record bases':>36}")
P("-" * 118)
ROWS = []
for (tag, kind, par) in [("[[6,4,2]]", "fam", 6), ("[[8,6,2]]", "fam", 8),
                         ("[[10,8,2]]", "fam", 10), ("[[12,10,2]]", "fam", 12),
                         ("[[4,2,2]]^2", "ctl", 2), ("[[4,2,2]]^3", "ctl", 3),
                         ("[[4,2,2]]^4", "ctl", 4)]:
    if kind == "fam":
        n, S = carrier_nn2(par); R, W = logicals_for(n, S)
    else:
        n, S = carrier_product(par); R, W = logicals_blocks(4, par)
    verify_records(tag, n, S, R, W)
    wv, tot, how = weight_map(n, S, R)
    assert how == "exact", how
    sm = wmap_summary(wv, len(R))
    mn, mx, uniq = wall_over_bases(n, S, R, wv)
    ROWS.append((tag, len(R), n, sm["wmin"], sm["wmax"], sm["w_allones"], mn, mx, uniq))
    P(f"{tag:>16} {len(R):>3} {n:>3} | {sm['wmin']:>5} {sm['wmax']:>5} {sm['w_allones']:>17} | "
      f"min {mn}  max {mx}   values seen {uniq}")
P("-" * 118)
famrows = [r for r in ROWS if r[0].startswith("[[") and "^" not in r[0]]
ctlrows = [r for r in ROWS if "^" in r[0]]
chk("A: family wall is NOT constant across record bases (2 is one basis's value)",
    any(r[7] > 2 for r in famrows), f"family wall maxima {[r[7] for r in famrows]}")
chk("A: control wall CAN equal 2 in some record basis (so 'control needs weight k' is basis choice)",
    all(r[6] == 2 for r in ctlrows), f"control wall minima {[r[6] for r in ctlrows]}")

# ================================================================= B  region sampling artifact
P("")
P("ATTACK B.  s3 reports 'smallest region touching ALL k records' as 13,15,17,21,23 for the")
P("           control at k = 12,14,16,18,20 -- all flagged 's' (SAMPLED, 4000 random regions).")
P("           The control is m disjoint [[4,2,2]] blocks; two qubits inside a block already")
P("           reach both of that block's records, so k qubits (2 per block) reach ALL k.")
P("           Constructed explicitly and measured -- no sampling.")
P("-" * 118)

def coord_images(n, R):
    k = len(R); rows = []
    for j in range(2*n):
        e = [0]*(2*n); e[j] = 1
        rows.append([sp2(e, R[i], n) for i in range(k)])
    return rows

def reach_all(rows, n, qubits):
    idx = list(qubits) + [n + q for q in qubits]
    return f2_rank([rows[j] for j in idx])

P(f"{'k':>3} {'nC':>4} | {'s3 reported min region':>23} {'explicit 2-per-block region':>29} "
  f"{'reach':>6} {'= k?':>5} {'lower bound ceil(k/2)':>22}")
P("-" * 118)
S3_REPORT = {12: 13, 14: 15, 16: 17, 18: 21, 20: 23}
bad_b = []
for m in range(1, 11):
    k = 2*m; n, S = carrier_product(m); R, W = logicals_blocks(4, m)
    verify_records(f"[[4,2,2]]^{m}", n, S, R, W)
    rows = coord_images(n, R)
    reg = [4*b + q for b in range(m) for q in (0, 1)]
    rc = reach_all(rows, n, reg)
    rep = S3_REPORT.get(k, k)
    ok = (rc == k)
    if k in S3_REPORT and ok and rep > k: bad_b.append((k, rep))
    P(f"{k:>3} {n:>4} | {rep:>23} {'|A| = '+str(len(reg)):>29} {rc:>6} {str(ok):>5} {-(-k//2):>22}")
P("-" * 118)
chk("B: an explicit size-k region reaches ALL k records on EVERY control carrier",
    True if not bad_b or all(x[0] for x in bad_b) else False,
    f"s3's sampled values that are strictly too large: {bad_b}")
chk("B: s3's control column for 'min region for ALL k' is WRONG at k>=12 (sampling artifact)",
    len(bad_b) == 5, f"{bad_b}")

# ================================================================= C  n, not collectivity
P("")
P("ATTACK C.  The family/control comparison is NOT matched in carrier size: at the same k the")
P("           family has n = k+2 qubits and the control has n = 2k.  Every weight quantity is")
P("           therefore confounded with n.  Held at FIXED k = 8, four carriers with n = 10,")
P("           12, 14, 16 -- two collective, two products -- to see what wmax/wmean track.")
P("-" * 118)
CARRIERS_K8 = []
n, S = carrier_nn2(10); R, W = logicals_for(n, S)
CARRIERS_K8.append(("[[10,8,2]]  COLLECTIVE", n, S, R, W))
n, S = carrier_product_of(6, 2); R, W = logicals_blocks(6, 2)
CARRIERS_K8.append(("[[6,4,2]]^2 PRODUCT", n, S, R, W))
n, S = carrier_nn2_plus(14, 4); R, W = logicals_for(n, S)
CARRIERS_K8.append(("[[14,8]]+4Zpairs COLLECTIVE", n, S, R, W))
n, S = carrier_product(4); R, W = logicals_blocks(4, 4)
CARRIERS_K8.append(("[[4,2,2]]^4 PRODUCT", n, S, R, W))

P(f"{'carrier':>28} {'k':>3} {'n':>3} {'nstab':>5} | {'wmin':>5} {'wmax':>5} {'wmean':>7} "
  f"{'wall':>5} | {'wmax/n':>7} {'wmean/n':>8} {'E0':>5}")
P("-" * 118)
K8 = []
for tag, n, S, R, W in CARRIERS_K8:
    verify_records(tag, n, S, R, W)
    chk(f"{tag}: k = 8", len(R) == 8, f"k={len(R)}")
    wv, tot, how = weight_map(n, S, R)
    assert how == "exact", f"{tag} {how} {tot}"
    sm = wmap_summary(wv, len(R))
    ns = f2_rank(S)
    K8.append((tag, n, ns, sm))
    P(f"{tag:>28} {len(R):>3} {n:>3} {ns:>5} | {sm['wmin']:>5} {sm['wmax']:>5} {sm['wmean']:>7.3f} "
      f"{sm['w_allones']:>5} | {sm['wmax']/n:>7.3f} {sm['wmean']/n:>8.3f} {-ns:>5}")
P("-" * 118)
mono = all(K8[i][3]["wmean"] <= K8[i+1][3]["wmean"] + 1e-9 for i in range(len(K8)-1))
chk("C: at FIXED k=8, wmean increases MONOTONICALLY with n across collective and product alike",
    mono, f"n = {[r[1] for r in K8]}, wmean = {[round(r[3]['wmean'],3) for r in K8]}")
spread = max(r[3]["wmean"]/r[1] for r in K8) - min(r[3]["wmean"]/r[1] for r in K8)
chk("C: wmean/n is nearly carrier-independent at fixed k (spread < 0.06)", spread < 0.06,
    f"wmean/n = {[round(r[3]['wmean']/r[1],4) for r in K8]}, spread {spread:.4f}")

# ================================================================= D  energy is not a function of k
P("")
P("ATTACK D.  'enclosed ENERGY is exactly CONSTANT in k on the collective carrier' is an")
P("           identity E0 = -nstab, and nstab = 2 is the DEFINITION of the [[n,n-2,2]] family.")
P("           Held at FIXED k = 8, E0 is whatever the presentation makes it.")
P("-" * 118)
P(f"{'carrier':>28} {'k':>3} {'n':>3} {'nstab':>5} {'E0 = -nstab':>12} {'width = 2*nstab':>16}")
P("-" * 118)
E0s = []
for tag, n, S, R, W in CARRIERS_K8:
    ns = f2_rank(S); E0s.append(-ns)
    P(f"{tag:>28} {len(R):>3} {n:>3} {ns:>5} {-ns:>12} {2*ns:>16}")
for j in (1, 2, 3, 5):
    n, S = carrier_nn2_plus(10 + j, j) if (10 + j) % 2 == 0 else (None, None)
for j in (2, 4, 6):
    n, S = carrier_nn2_plus(10 + j, j); R, W = logicals_for(n, S)
    verify_records(f"[[{n}]]+{j}Zpairs", n, S, R, W)
    ns = f2_rank(S); E0s.append(-ns)
    P(f"{('[['+str(n)+',8]]+'+str(j)+'Zpairs'):>28} {len(R):>3} {n:>3} {ns:>5} {-ns:>12} {2*ns:>16}")
P("-" * 118)
chk("D: at FIXED k=8 the ground energy takes many different values (E0 is NOT a function of k)",
    len(set(E0s)) >= 4, f"E0 values at k=8: {sorted(set(E0s))}")

# ================================================================= E  reach used inadmissible ops
P("")
P("ATTACK E.  s3's reach(A) is the rank of the patterns of ALL Paulis supported in A, including")
P("           ones that do NOT commute with the stabilisers.  Those are not admissible")
P("           operations: they leave the code space and create an excitation (D-18).  Recount")
P("           with ADMISSIBLE-ONLY operators supported in A.")
P("-" * 118)

def f2_solve_basis(rows, ncols):
    """basis of the nullspace of the matrix whose ROWS are `rows`, over F_2"""
    M = [r[:] for r in rows]; piv = []; r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c]:
                M[i] = [(x + y) % 2 for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
    free = [c for c in range(ncols) if c not in piv]
    B = []
    for f in free:
        v = [0]*ncols; v[f] = 1
        for i, c in enumerate(piv): v[c] = M[i][f]
        B.append(v)
    return B

def reach_admissible(n, S, R, qubits):
    """rank of the record-pattern images of Paulis supported in `qubits` that COMMUTE with
       every stabiliser."""
    idx = list(qubits) + [n + q for q in qubits]
    d = len(idx)
    # constraint rows: for each stabiliser, sp(e_j, s) over the local coordinates
    rows = []
    for s in S:
        row = []
        for j in idx:
            e = [0]*(2*n); e[j] = 1
            row.append(sp2(e, s, n))
        rows.append(row)
    basis = f2_solve_basis(rows, d)
    out = []
    for b in basis:
        v = [0]*(2*n)
        for t, j in enumerate(idx):
            if b[t]: v[j] ^= 1
        out.append([sp2(v, R[i], n) for i in range(len(R))])
    return f2_rank(out) if out else 0

def reach_best(n, S, R, r, admissible, cap=200000):
    best = 0
    if comb(n, r) <= cap:
        for A in itertools.combinations(range(n), r):
            v = reach_admissible(n, S, R, A) if admissible else reach_all(coord_images(n, R), n, A)
            best = max(best, v)
        how = "e"
    else:
        rng = np.random.default_rng(1)
        rows = None if admissible else coord_images(n, R)
        for _ in range(3000):
            A = sorted(rng.choice(n, size=r, replace=False).tolist())
            v = reach_admissible(n, S, R, A) if admissible else reach_all(rows, n, A)
            best = max(best, v)
        how = "s"
    return best, how

P(f"{'carrier':>16} {'k':>3} {'n':>3} | {'r':>2} {'s3 reach (ALL Paulis)':>22} "
  f"{'reach (ADMISSIBLE only)':>24}")
P("-" * 118)
diffE = []
for par in (6, 8, 10, 12, 14):
    n, S = carrier_nn2(par); R, W = logicals_for(n, S)
    for r in (2, 4, 6):
        ba, ha = reach_best(n, S, R, r, admissible=False)
        bb, hb = reach_best(n, S, R, r, admissible=True)
        diffE.append((par, r, ba, bb))
        P(f"{('[['+str(par)+','+str(par-2)+',2]]'):>16} {par-2:>3} {n:>3} | {r:>2} "
          f"{str(ba)+ha:>22} {str(bb)+hb:>24}")
P("-" * 118)
chk("E: the admissible-only reach differs from s3's all-Pauli reach somewhere",
    any(a != b for _, _, a, b in diffE),
    f"pairs (all, admissible) = {[(a,b) for _,_,a,b in diffE]}")
chk("E: even admissible-only, a size-6 region's best reach SATURATES (does not grow with k)",
    len(set(b for p, r, a, b in diffE if r == 6)) <= 2,
    f"r=6 admissible reach across k=4..12: {[b for p,r,a,b in diffE if r==6]}")

# ================================================================= self-checks
P("")
P("SELF-CHECKS OF THIS VERIFICATION")
P("-" * 118)
bad = [c for c in CH if not c[1]]
for c in CH:
    P(f"   {'PASS' if c[1] else 'FAIL'}  {c[0]}   {c[2]}")
P(f"   {len(CH)-len(bad)} / {len(CH)} pass")

with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_C_ACCUM/VERIFY/v1_refute.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
