#!/usr/bin/env python3
# =============================================================================
# T-51 FOURTH DESIGN -- GEOMETRY RESPONSE: is the earned metric responsive to
# written content, or provably rigid?
#
# This is the owed fourth instrument design of T-51's round (the original failed
# in flight; register debt in C-92's row).  Angle: the field as a response of the
# earned geometry itself -- d_W (C-78), the writer-cost landscape w_min/w_enc
# (C-80), cut-rank (C-81) -- conditioned on content written vs unwritten.
#
# VERDICT COMPUTED HERE: the rigidity branch is real and PROVABLE.  On stabiliser
# venues every earned-geometry integer is exactly content-independent (three
# mechanisms, Parts A/B/D), while the real-valued amplitude functional OVER the
# rigid class structure responds (Parts C/E) -- with onset order equal to the
# CONNECTED w_enc per placement, independently reproducing the T-51 auditor's
# counterexample (old-spec 3 / connected 4 / measured onset 4.0) from separate
# machinery.
#
# Discipline: D-1 absolute (no classical form required or tested anywhere; every
# shape below is an output); D-24 (all separations in earned quantities; lattice
# coordinates appear only as construction labels, never in claims); D-15 (two-way
# controls; certificates labeled as certificates, never counted as controls);
# D-8 (every gated verdict could come out the other way -- gates marked
# CERTIFICATE are algebraic identities and are labeled so).
#
# Sealed machinery: model/geometry.py (C-78 dW_class_matrix, C-80 conventions
# coset_min_np/coupling_cost/Torus, C-81 _subgroup_in), read-only.
# Writes: stdout only (redirect to FOURTH_rigidity_OUT.txt).
# =============================================================================
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
import numpy as np
import geometry as G
from geometry import (toric_stabilizers, sp_pair, supp_mask, vec_to_mask, pc,
                      solve_affine_f2, span_numpy, Torus, _independent_subset,
                      _generator_graph_dist, _subgroup_in, rank_bits)
from record_model import symplectic_logicals

REPORT = []
def gate(name, ok, detail=""):
    REPORT.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")

# ============================================================================
# PART A -- C-78 tier: d_W under content writes.
# Independent per-pair recomputation (no difference-class cache), then the
# translation-invariance identity.  R-I mechanism: an admissible writer's label
# action a(w) = (sp(w, Zbar_1), sp(w, Zbar_2)) is a function of w alone, so
# {w : takes s to s'} = {w : a(w) = s XOR s'} -- the same set for every
# translate (s+c, s'+c).  A content write IS a translation of the base point.
# Hence d_W cannot respond, at any L, algebraically.  The computation exhibits
# the identity at L = 2, 3, 4 (CERTIFICATE: it cannot fail given R-I).
# ============================================================================
print("=" * 78)
print("PART A -- d_W under content writes (C-78 instrument, per-pair recomputation)")
print("=" * 78)

def dW_per_pair(L):
    n, stab, stars, plaqs = toric_stabilizers(L)
    pairs = symplectic_logicals(stab, n)
    logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
    Zl_z = [(l[1] >> n) for l in logi]
    plaq_rows = [p for p in plaqs]
    configs = [(a, b) for a in range(2) for b in range(2)]
    dmat = {}
    for s in configs:
        for sp_ in configs:
            t = (s[0] ^ sp_[0], s[1] ^ sp_[1])
            rows = plaq_rows + [Zl_z[0], Zl_z[1]]
            rhs = [0] * len(plaq_rows) + [t[0], t[1]]
            x0, nb, _r = solve_affine_f2(rows, rhs, n)
            assert x0 is not None
            coset = span_numpy(nb) ^ np.uint64(x0)
            dmat[(s, sp_)] = int(np.bitwise_count(coset).min())
    return dmat

for L in (2, 3, 4):
    dpp = dW_per_pair(L)
    sealed = G.dW_class_matrix(L)
    ok_seal = all(dpp[k] == sealed["dmat"][k] for k in dpp)
    ok_trans = True
    cfgs = [(a, b) for a in range(2) for b in range(2)]
    for s in cfgs:
        for sp_ in cfgs:
            for c in cfgs:
                sc = (s[0] ^ c[0], s[1] ^ c[1]); spc = (sp_[0] ^ c[0], sp_[1] ^ c[1])
                if dpp[(sc, spc)] != dpp[(s, sp_)]:
                    ok_trans = False
    gate(f"A.dW L={L} per-pair recomputation == sealed class matrix", ok_seal,
         f"class_min={sealed['class_min']}")
    gate(f"A.dW L={L} invariant under all content writes (CERTIFICATE, R-I)",
         ok_trans, "d_W(s+c,s'+c)=d_W(s,s'), all 64 triples")

# ============================================================================
# Shared hole-venue machinery (C-80 conventions).
# Probe: star-hole pair (adjacent vertices; direct Z-connector weight 1).
# Source: plaquette-hole pair; its value b is the WRITTEN CONTENT -- written
# (b=-1) and unwritten (b=+1) are two exact sectors of ONE Hamiltonian.
# ============================================================================
def edge_vertices(T, e):
    LxLy = T.Lx * T.Ly
    if e < LxLy:
        y, x = divmod(e, T.Lx)
        return ((x, y), ((x + 1) % T.Lx, y))
    y, x = divmod(e - LxLy, T.Lx)
    return ((x, y), (x, (y + 1) % T.Ly))

def connected_probe_path(T, PV, w):
    """CONNECTED w_enc gate (the judged pre-commission repair): z-support is one
       connected edge set whose odd-degree vertices are exactly the probe pair."""
    from collections import defaultdict, deque
    n = T.n
    es = [e for e in range(n) if (w >> (n + e)) & 1]
    if not es:
        return False
    deg = defaultdict(int); vert_edges = defaultdict(set)
    for e in es:
        a, b = edge_vertices(T, e)
        deg[a] += 1; deg[b] += 1
        vert_edges[a].add(e); vert_edges[b].add(e)
    if {v for v, d in deg.items() if d % 2 == 1} != set(PV):
        return False
    seen = {es[0]}; dq = deque([es[0]])
    while dq:
        e = dq.popleft()
        a, b = edge_vertices(T, e)
        for e2 in vert_edges[a] | vert_edges[b]:
            if e2 not in seen:
                seen.add(e2); dq.append(e2)
    return len(seen) == len(es)

def probe_coset(T, PV, src, span_kind="Z"):
    """The probe's label-exact writer class: sp(w, A_p)=1 for both probe stars,
       sp(w, q)=0 for every remaining star, sp(w, p)=0 for EVERY plaquette
       (source included: no co-write -- this is what makes the class identical
       between carrier-present and carrier-absent venues, R-II).
       Span 'Z' mirrors the sealed coupling_cost convention (independent
       plaquettes + both Z logicals); 'ZX' adds the closed X-loop group
       (remaining stars + X logicals): Y-dressing allowed, a could-fail check
       that dressing does not beat the pure-Z minima."""
    n = T.n
    A_p = [T.star(*v) for v in PV]
    rem_stars = [T.star(x, y) for y in range(T.Ly) for x in range(T.Lx)
                 if (x, y) not in PV]
    plaqs = T.all_plaqs()
    w0 = 1 << (n + T.h(PV[0][0], PV[0][1]))     # direct connector: shared edge
    gens = list(_independent_subset(plaqs)) + [T.zbar1(), T.zbar2()]
    if span_kind == "ZX":
        gens += list(_independent_subset(rem_stars)) + [T.xbar1(), T.xbar2()]
    assert all(sp_pair(w0, a, n) == 1 for a in A_p)
    assert all(sp_pair(w0, q, n) == 0 for q in rem_stars)
    assert all(sp_pair(w0, p, n) == 0 for p in plaqs)
    arr = span_numpy([np.uint64(g) for g in gens]) ^ np.uint64(w0)
    W_src = T.dual_path_x(*src)     # declared crossing detector (winding class
                                    # convention: the sealed dual_path_x route)
    out = []
    for w in map(int, arr):
        assert sp_pair(w, A_p[0], n) == 1 and sp_pair(w, A_p[1], n) == 1
        assert all(sp_pair(w, q, n) == 0 for q in rem_stars)
        assert all(sp_pair(w, p, n) == 0 for p in plaqs)
        out.append((w, pc(supp_mask(w, n)), sp_pair(w, W_src, n)))
    return out

def histo(elems, pred=lambda r: True):
    h = {}
    for r in elems:
        if pred(r):
            h[r[1]] = h.get(r[1], 0) + 1
    return dict(sorted(h.items()))

def static_columns(T, PV, src):
    elems = probe_coset(T, PV, src, "Z")
    w_direct = min(r[1] for r in elems if r[2] == 0)
    w_enc_old = min(r[1] for r in elems if r[2] == 1)
    conn = [r[1] for r in elems if r[2] == 1 and connected_probe_path(T, PV, r[0])]
    return elems, w_direct, w_enc_old, (min(conn) if conn else None)

# ============================================================================
# PART B -- the static writer-cost landscape around written vs unwritten
# content, 3x2 torus (R-II mechanism).  The class is defined by symplectic
# pairings against OPERATORS; the written value b lives in the STATE; no term
# of the enumeration reads b.  Full weight histograms computed; the auditor's
# connectivity counterexample reproduced with independent machinery.
# ============================================================================
print()
print("=" * 78)
print("PART B -- probe writer-cost landscape, hole venue 3x2 (static tier)")
print("=" * 78)
T2 = Torus(3, 2)
PV2 = [(0, 0), (1, 0)]
for src in [((1, 1), (2, 1)), ((0, 1), (2, 1))]:
    elems, w_d, w_old, w_conn = static_columns(T2, PV2, src)
    print(f"  source {src}: histogram(all) {histo(elems)}")
    print(f"    enclosing-class histogram {histo(elems, lambda r: r[2] == 1)}")
    print(f"    w_direct={w_d}  w_enc(old-spec)={w_old}  w_enc(CONNECTED)={w_conn}")
    elems_zx = probe_coset(T2, PV2, src, "ZX")
    gate(f"B.mixed-span (Y-dressed) min == pure-Z min, source {src}",
         min(histo(elems_zx)) == min(histo(elems)),
         f"{min(histo(elems_zx))} vs {min(histo(elems))}")
gate("B.auditor counterexample reproduced independently: "
     "old-spec w_enc=3, CONNECTED w_enc=4 at source ((1,1),(2,1))",
     static_columns(T2, PV2, ((1, 1), (2, 1)))[2:] == (3, 4)
     if True else False,
     f"got {static_columns(T2, PV2, ((1, 1), (2, 1)))[2:]}")
print("  NOTE (R-II): no term of these enumerations takes the written value b as")
print("  input -- histogram equality between written/unwritten sectors, and between")
print("  carrier-present/carrier-absent venues (sp(w,p)=0 for all p either way), is")
print("  an algebraic identity: CERTIFICATE, never counted as a control.  The")
print("  could-fail shadow of this certificate is Part C's onset-order battery.")
loc = T2.all_stars() + T2.all_plaqs()
for src in [((1, 1), (2, 1)), ((0, 1), (2, 1))]:
    dgs = [_generator_graph_dist(loc, T2.n, (pv[1] * 3 + pv[0]),
                                 6 + (sv[1] * 3 + sv[0]))
           for pv in PV2 for sv in src]
    print(f"  d_gen probe-to-source {src}: min {min(dgs)}  (DESCRIPTIVE ONLY -- "
          f"mixed-type separation carries no sealed law; per the T-51 judgment,")
print("  connected w_enc alone carries the earned separation)")

# ============================================================================
# Sector-exact dynamics machinery (one Hamiltonian; written/unwritten are two
# exact sectors: the mediator V = lam*sum_e Z_e conserves every plaquette).
# ============================================================================
def build_venue(Lx, Ly, PV, src_plaqs):
    T = Torus(Lx, Ly)
    rem_stars = [T.star(x, y) for y in range(Ly) for x in range(Lx)
                 if (x, y) not in PV]
    rem_plaqs = [T.plaq(x, y) for y in range(Ly) for x in range(Lx)
                 if (x, y) not in src_plaqs]
    src = [T.plaq(*p) for p in src_plaqs]
    return T, rem_stars, rem_plaqs, src

def sector_states(T, rem_plaqs, src, b, windings=(0, 0)):
    n = T.n
    xs = np.arange(1 << n, dtype=np.uint64)
    sel = np.ones(1 << n, dtype=bool)
    for p in rem_plaqs:
        sel &= (np.bitwise_count(xs & np.uint64(p >> n)) % 2 == 0)
    for p in src:
        sel &= (np.bitwise_count(xs & np.uint64(p >> n)) % 2
                == (0 if b == +1 else 1))
    for wv, m in zip(windings, [T.zbar1(), T.zbar2()]):
        sel &= (np.bitwise_count(xs & np.uint64(m >> n)) % 2 == wv)
    return xs[sel]

def sector_H(T, rem_stars, rem_plaqs, states, lam, probe_bias=0.0, A_pick=None):
    n = T.n
    dim = len(states)
    pos = {int(s): i for i, s in enumerate(states)}
    H = np.zeros((dim, dim))
    for i, s in enumerate(map(int, states)):
        d = 0.0
        for p in rem_plaqs:
            d -= (-1.0) ** pc(s & (p >> n))
        for e in range(n):
            d += lam * ((-1.0) ** ((s >> e) & 1))
        H[i, i] = d
    terms = [(a, 1.0) for a in rem_stars]
    if probe_bias and A_pick is not None:
        terms.append((A_pick, probe_bias))
    for a, c in terms:
        xm = a & ((1 << n) - 1)
        for i, s in enumerate(map(int, states)):
            H[pos[s ^ xm], i] += -c
    return H, pos

def doublet(T, rem_stars, rem_plaqs, src, b, lam, windings, A_probe):
    st = sector_states(T, rem_plaqs, src, b, windings)
    H, pos = sector_H(T, rem_stars, rem_plaqs, st, lam)
    ev, evec = np.linalg.eigh(H)
    xm = A_probe & ((1 << T.n) - 1)
    perm = np.array([pos[int(s) ^ xm] for s in st])
    Aop = np.zeros((len(st), len(st)))
    Aop[perm, np.arange(len(st))] = 1.0
    wit = abs(evec[:, 0] @ Aop @ evec[:, 1])
    return ev[1] - ev[0], wit, st

def onset(vals):
    return [float(np.log(abs(vals[i + 1]) / abs(vals[i])) / np.log(2.0))
            for i in range(len(vals) - 1) if vals[i] != 0 and vals[i + 1] != 0]

# ============================================================================
# PART C -- 3x2 dynamics: the two-tier split, per placement (R-II's could-fail
# shadow).  Integer tier: Delta's onset order (= w_direct) equal in BOTH
# sectors; F's onset order = CONNECTED w_enc OF THAT PLACEMENT (4 at the
# auditor's placement, 3 at the second -- both could come out otherwise).
# Amplitude tier: F != 0 -- the response exists and lives ONLY here.
# ============================================================================
print()
print("=" * 78)
print("PART C -- sector-exact dynamics 3x2: rigid integers, responsive amplitudes")
print("=" * 78)
lams = [0.004, 0.008, 0.016, 0.032, 0.064]
for src, w_conn_expect in [(((1, 1), (2, 1)), 4), (((0, 1), (2, 1)), 3)]:
    T, rs, rp, sr = build_venue(3, 2, PV2, list(src))
    A_probe = T.star(0, 0)
    rows = {}
    for b in (+1, -1):
        for lam in lams + [0.05]:
            D, wit, st = doublet(T, rs, rp, sr, b, lam, (0, 0), A_probe)
            assert len(st) == 32
            rows[(b, lam)] = (D, wit)
    print(f"  source {src} (connected w_enc = {w_conn_expect}), winding (0,0):")
    print(f"  {'lam':>7} {'Delta(b=+1)':>14} {'Delta(b=-1)':>14} {'F':>14} "
          f"{'wit+':>6} {'wit-':>6}")
    for lam in lams + [0.05]:
        Dp, wp = rows[(+1, lam)]; Dm, wm = rows[(-1, lam)]
        print(f"  {lam:7.3f} {Dp:14.6e} {Dm:14.6e} {Dm - Dp:14.6e} "
              f"{wp:6.3f} {wm:6.3f}")
    oDp = onset([rows[(+1, l)][0] for l in lams])
    oDm = onset([rows[(-1, l)][0] for l in lams])
    oF = onset([rows[(-1, l)][0] - rows[(+1, l)][0] for l in lams])
    print(f"    onset ladders: Delta+ {['%.3f' % o for o in oDp]}  "
          f"Delta- {['%.3f' % o for o in oDm]}  F {['%.3f' % o for o in oF]}")
    gate(f"C.{src} Delta onset = 1 = w_direct in BOTH sectors [integer tier]",
         abs(oDp[-1] - 1) < 0.1 and abs(oDm[-1] - 1) < 0.1,
         f"slopes {oDp[-1]:.3f}, {oDm[-1]:.3f}")
    gate(f"C.{src} F onset -> CONNECTED w_enc = {w_conn_expect} of this placement",
         abs(oF[0] - w_conn_expect) < 0.1,
         f"first-rung slope {oF[0]:.3f} (ladder {['%.2f' % o for o in oF]})")
    gate(f"C.{src} F != 0 at lam=0.05 [amplitude tier responds]",
         abs(rows[(-1, 0.05)][0] - rows[(+1, 0.05)][0]) > 1e-10,
         f"F = {rows[(-1, 0.05)][0] - rows[(+1, 0.05)][0]:+.3e}")
# winding annex (the judged mandatory sweep before any sign attribution)
T, rs, rp, sr = build_venue(3, 2, PV2, [(1, 1), (2, 1)])
print("  winding annex, source ((1,1),(2,1)), lam=0.05 "
      "(sign is a winding-sector quantity on small tori):")
for wsec in [(0, 0), (1, 0), (0, 1), (1, 1)]:
    FF = {b: doublet(T, rs, rp, sr, b, 0.05, wsec, T.star(0, 0))[0]
          for b in (+1, -1)}
    print(f"    winding {wsec}: F = {FF[-1] - FF[+1]:+.6e}")

# ============================================================================
# PART D -- region response (C-81 tier, R-III mechanism).  At lam=0: a region's
# reduced state responds to a distant write IFF the region supports an
# enclosure detector (a Z-loop separating the source holes -- equivalently,
# every admissible representative of the content's writer crosses it oddly);
# where it responds the response is maximal (the region READS the bit), and
# even there the region's entropy -- the integer tier -- is exactly unmoved
# (the write is Pauli transport; spectra are conserved).  Cut-rank: state-free.
# At lam>0: the ENTROPY differential becomes nonzero -- the amplitude-tier
# response, concentrated on crossing regions.  All four directions could fail.
# ============================================================================
print()
print("=" * 78)
print("PART D -- region response at lam=0 and lam>0 (3x2, source ((1,1),(2,1)))")
print("=" * 78)
T, rs, rp, sr = build_venue(3, 2, PV2, [(1, 1), (2, 1)])
A_probe = T.star(0, 0)
n2 = T.n

def ground_full(b, lam):
    st = sector_states(T, rp, sr, b, (0, 0))
    H, pos = sector_H(T, rs, rp, st, lam,
                      probe_bias=(0.1 if lam == 0 else 0.0), A_pick=A_probe)
    ev, evec = np.linalg.eigh(H)
    psi = np.zeros(1 << n2)
    psi[np.asarray(st, dtype=np.int64)] = evec[:, 0]
    return psi

def rho_region(psi, edges):
    rest = [e for e in range(n2) if e not in edges]
    idx = np.arange(1 << n2)
    a = np.zeros(1 << n2, dtype=np.int64)
    for i, e in enumerate(edges):
        a |= ((idx >> e) & 1) << i
    bp = np.zeros(1 << n2, dtype=np.int64)
    for i, e in enumerate(rest):
        bp |= ((idx >> e) & 1) << i
    M = np.zeros((1 << len(edges), 1 << len(rest)))
    M[a, bp] = psi
    return M @ M.T

def S_bits(rho):
    p = np.linalg.eigvalsh(rho)
    p = p[p > 1e-14]
    return float(-(p * np.log2(p)).sum())

R_cross = sorted({T.h(2, 1), T.h(2, 0), T.v(2, 1), T.v(0, 1)})  # encloses src hole
R_away = sorted({T.h(0, 0), T.h(1, 0), T.v(0, 0), T.v(1, 0)})   # avoidable by class
print(f"  R_cross (edge labels) = {R_cross}   R_away = {R_away}")
for lam in (0.0, 0.05):
    psis = {b: ground_full(b, lam) for b in (+1, -1)}
    for name, R in [("R_away ", R_away), ("R_cross", R_cross)]:
        r = {b: rho_region(psis[b], R) for b in (+1, -1)}
        dtr = 0.5 * np.abs(np.linalg.eigvalsh(r[+1] - r[-1])).sum()
        S = {b: S_bits(r[b]) for b in (+1, -1)}
        print(f"  lam={lam:4.2f} {name}: trace-dist = {dtr:.3e}   "
              f"S(+1)={S[+1]:.6f}  S(-1)={S[-1]:.6f}  dS={S[-1] - S[+1]:+.3e} bits")
        if lam == 0.0 and name.strip() == "R_away":
            gate("D.lam=0 R_away: reduced state EXACTLY content-independent",
                 dtr < 1e-10, f"trace-dist {dtr:.1e}")
        if lam == 0.0 and name.strip() == "R_cross":
            gate("D.lam=0 R_cross: reduced STATE responds (it reads the bit)",
                 dtr > 1e-3, f"trace-dist {dtr:.1e}")
            gate("D.lam=0 R_cross: ENTROPY (integer tier) exactly unmoved",
                 abs(S[+1] - S[-1]) < 1e-10, f"dS {S[-1] - S[+1]:.1e}")
gens = rs + rp
full = (1 << (2 * n2)) - 1
for name, R in [("R_away", R_away), ("R_cross", R_cross)]:
    inside = 0
    for e in R:
        inside |= (1 << e) | (1 << (n2 + e))
    cr = (rank_bits(gens) - len(_subgroup_in(gens, inside, full))
          - len(_subgroup_in(gens, full & ~inside, full)))
    print(f"  cut-rank({name}) = {cr}  (state-free by construction: CERTIFICATE "
          f"-- C-81's static arm on this venue)")

# ============================================================================
# PART E -- 3x3: the same two-tier split at two earned separations
# (256-dim sectors; the commissioned Second Lump venue).
# ============================================================================
print()
print("=" * 78)
print("PART E -- 3x3 confirmation at two earned separations")
print("=" * 78)
T3 = Torus(3, 3)
PV3 = [(0, 0), (1, 0)]
for src3 in [((1, 1), (2, 1)), ((1, 2), (2, 2))]:
    _, w_d, w_old, w_conn = static_columns(T3, PV3, src3)
    T, rs, rp, sr = build_venue(3, 3, PV3, list(src3))
    A3 = T.star(0, 0)
    rows = {}
    for b in (+1, -1):
        for lam in lams:
            D, wit, st = doublet(T, rs, rp, sr, b, lam, (0, 0), A3)
            assert len(st) == 256
            rows[(b, lam)] = D
    oF = onset([rows[(-1, l)] - rows[(+1, l)] for l in lams])
    oDp = onset([rows[(+1, l)] for l in lams])
    oDm = onset([rows[(-1, l)] for l in lams])
    F05 = rows[(-1, 0.064)] - rows[(+1, 0.064)]
    print(f"  source {src3}: w_direct={w_d} w_enc(old)={w_old} w_enc(CONN)={w_conn}")
    print(f"    Delta+ onset {['%.3f' % o for o in oDp]}  "
          f"Delta- onset {['%.3f' % o for o in oDm]}  "
          f"F onset {['%.3f' % o for o in oF]}  F(lam=0.064)={F05:+.3e}")
    gate(f"E.{src3} F onset first rung within 0.15 of CONNECTED w_enc = {w_conn}",
         abs(oF[0] - w_conn) < 0.15, f"slope {oF[0]:.3f}")
    gate(f"E.{src3} Delta onset = 1 both sectors [integer tier]",
         abs(oDp[0] - 1) < 0.05 and abs(oDm[0] - 1) < 0.05,
         f"slopes {oDp[0]:.3f}, {oDm[0]:.3f}")

print()
print("=" * 78)
print("GATE SUMMARY")
print("=" * 78)
npass = sum(1 for _, ok, _ in REPORT if ok)
for name, ok, det in REPORT:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
print(f"{npass}/{len(REPORT)} gates pass")
