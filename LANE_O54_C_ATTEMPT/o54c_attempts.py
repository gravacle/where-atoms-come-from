"""LANE_O54_C_ATTEMPT -- THE ADVERSARIAL CONSTRUCTION ATTEMPT: try hard to build a
PURE-GAMMA FALLOFF between two separated regions.

THE HYPOTHESIS UNDER ATTACK (O-54, stated so it can fail):
  "the falloff is not Gamma's to give -- pure-Gamma two-region relations are TOPOLOGICAL
   (separation-independent or contact-only), and any genuine falloff enters through the
   MEDIATOR (a coupling term under the openness clause)."
This lane's job is to make that FAIL: construct a separation-DEPENDENT pure-Gamma
two-region quantity and check whether any such dependence has falloff character.

ATTEMPTS:
  (a) homology-lattice coupling of two defects ("handles" realized as plaquette punctures,
      the puncture route to enlarged H_1; owner of defect-based logical qubits:
      Bravyi-Kitaev quant-ph/9811052) -- does the intersection/Gram structure or any
      interface aggregate vary with earned separation?
  (b) region families on the torus -- stabilizer mutual information and interface rank
      between two separated regions vs their earned separation (owners of topological
      mutual information: Kitaev-Preskill 2006, Levin-Wen 2006; stabilizer entropy formula:
      Fattal-Cubitt-Yamamoto-Bravyi-Chuang quant-ph/0406168; stabilizer perimeter law:
      Hamma-Ionicioiu-Zanardi 2005).
  (c) the writer-weight landscape -- minimal weight of an admissible operation coupling
      the record at hole A to the record at hole B, vs earned separation (C-78 instrument:
      d_W = minimal admissible-writer weight; exhaustive coset scans, exact F_2).

DISCIPLINE: D-1 no classical gravitational form in any construction step (the comparison
target is named ONLY in the final comparison section).  D-15 controls beside every zero and
every flat line.  D-24 audit table at the end.  Exact arithmetic throughout; every
PASS/FAIL gated by a computed boolean.  INSERTED vs INDUCED labeled on every effect.
Carrier venue: toric code (Kitaev quant-ph/9707021); venue Aut structure cited from C-74.
"""
import sys, time
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O54_C_ATTEMPT")
from o54c_lib import (pc, rank_f2, supp_mask, weight_xz, sp_pair, Torus,
                      region_entropy, mutual_info, straddle_count, qubit_graph,
                      graph_dist_regions, generator_graph_dist, coset_min_np)
from collections import defaultdict

T0 = time.time()
GATES = []
def gate(name, ok, extra=""):
    GATES.append((name, bool(ok)))
    print(("PASS  " if ok else "FAIL  ") + name + (("  " + extra) if extra else ""))

def independent_subset(gens):
    out, piv = [], {}
    for g in gens:
        m = g
        while m:
            t = m.bit_length() - 1
            if t in piv:
                m ^= piv[t]
            else:
                piv[t] = m
                out.append(g)
                break
    return out

def classify(rows):
    """rows: list of (d, value) with d >= separated threshold. Exact classification."""
    if not rows:
        return "EMPTY"
    grp = defaultdict(list)
    for d, v in rows:
        grp[d].append(v)
    ds = sorted(grp)
    if all(v == rows[0][1] for _, v in rows):
        return "CONSTANT"
    inc = all(max(grp[ds[i]]) <= min(grp[ds[i + 1]]) for i in range(len(ds) - 1))
    dec = all(min(grp[ds[i]]) >= max(grp[ds[i + 1]]) for i in range(len(ds) - 1))
    if inc and not dec:
        return "INCREASING"
    if dec and not inc:
        return "DECREASING"
    return "MIXED"

SEQUENCES = []   # (label, origin INDUCED/INSERTED, rows [(d, val)], cov) for the falloff
                 # detector.  cov = None when the value is expected to be a function of d
                 # alone; otherwise a list of (covariate_key, d, val) rows -- the detector
                 # then requires the value to be resolved WITHIN each covariate class
                 # (constant or non-decreasing in d), and a falloff verdict would need a
                 # decreasing-in-d dependence either globally or inside some class.

def add_seq(label, origin, rows, cov=None):
    SEQUENCES.append((label, origin, rows, cov))

def is_function_of_d(rows):
    grp = defaultdict(set)
    for d, v in rows:
        grp[d].add(v)
    return all(len(s) == 1 for s in grp.values())

print("=" * 100)
print("O54-C  ADVERSARIAL ATTEMPT AT A PURE-GAMMA FALLOFF   (exact F_2; no floats on the measurement path)")
print("=" * 100)

# =====================================================================================
# SECTION 0 -- CARRIER VERIFICATION (every venue used below)
# =====================================================================================
print("\n-- SECTION 0: carrier verification --")
VENUES = {}
for (Lx, Ly) in [(8, 8), (12, 12), (3, 3), (3, 7), (4, 6)]:
    T = Torus(Lx, Ly)
    stars, plaqs = T.all_stars(), T.all_plaqs()
    n = T.n
    comm_ok = all(sp_pair(s, p, n) == 0 for s in stars for p in plaqs)
    rs, rp = rank_f2(stars), rank_f2(plaqs)
    logi = [T.xbar1(), T.xbar2(), T.zbar1(), T.zbar2()]
    # pairing matrix of logicals: xbar1-zbar2 and xbar2-zbar1 conjugate pairs
    P = [[sp_pair(a, b, n) for b in logi] for a in logi]
    want = [[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]]
    gate("carrier (%d,%d): stars x plaqs all commute" % (Lx, Ly), comm_ok)
    gate("carrier (%d,%d): rank(stars)=LxLy-1 and rank(plaqs)=LxLy-1" % (Lx, Ly),
         rs == Lx * Ly - 1 and rp == Lx * Ly - 1, "rs=%d rp=%d" % (rs, rp))
    gate("carrier (%d,%d): logical pairing matrix as computed conjugate structure" % (Lx, Ly),
         P == want)
    VENUES[(Lx, Ly)] = (T, stars, plaqs)

# =====================================================================================
# SECTION A -- ATTEMPT (b): TWO-REGION AGGREGATES ON THE TORUS vs EARNED SEPARATION
# =====================================================================================
print("\n" + "=" * 100)
print("SECTION A -- attempt (b): two-region interface aggregates vs earned separation")
print("  EARNED SEPARATION INSTRUMENT: BFS distance on the interaction graph computed from")
print("  the Hamiltonian generators' supports alone (no coordinates imported); C-78 lineage.")
print("  SEPARATED := no single generator touches both regions (straddle count 0).")
print("=" * 100)

for L in (8, 12):
    T, stars, plaqs = VENUES[(L, L)]
    n = T.n
    local = stars + plaqs
    adj = qubit_graph(local, n)
    states = {
        "Zbar-fixed": local + [T.zbar1(), T.zbar2()],
        "Xbar-fixed": local + [T.xbar1(), T.xbar2()],
    }
    for sname, gens in states.items():
        G = rank_f2(gens)
        gate("A L=%d state %s is pure (independent-generator count == n)" % (L, sname),
             G == n, "G=%d n=%d" % (G, n))

    # ---- patch family (contractible regions), axis and diagonal placements
    print("\nA.%d two 2x2 patches on the (%d,%d) torus  (regions: all edges in the window)" % (L, L, L))
    A = T.patch(0, 0, 2)
    rows_patch = []
    print("  placement       d_G  straddle  S_A S_B S_AB  MI   IR(AuB)")
    placements = [("axis x0=%d" % x0, T.patch(x0, 0, 2)) for x0 in range(2, L - 1)] + \
                 [("diag x0=%d" % x0, T.patch(x0, x0, 2)) for x0 in range(2, L - 2)]
    for label, B in placements:
        if set(A) & set(B):
            continue
        d = graph_dist_regions(adj, A, B)
        st = straddle_count(local, n, A, B)
        comp = [q for q in range(n) if q not in set(A) | set(B)]
        ir_union = straddle_count(local, n, sorted(set(A) | set(B)), comp)
        for sname, gens in states.items():
            SA, SB, SAB, MI = mutual_info(gens, n, n, A, B)
            rows_patch.append((label, sname, d, st, SA, SB, SAB, MI, ir_union))
    for r in rows_patch:
        print("  %-13s %-10s d=%d  str=%d  %2d %2d %3d   %2d   %3d" %
              (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]))
    sep_rows = [r for r in rows_patch if r[3] == 0]
    con_rows = [r for r in rows_patch if r[3] > 0]
    gate("A L=%d patches: MI == 0 on EVERY separated placement, both states" % L,
         all(r[7] == 0 for r in sep_rows), "%d separated rows" % len(sep_rows))
    gate("A L=%d patches: D-15 contact control: MI > 0 at contact (both states)" % L,
         len(con_rows) > 0 and all(r[7] > 0 for r in con_rows),
         "contact MI values %s" % sorted(set(r[7] for r in con_rows)))
    gate("A L=%d patches: interface rank between separated regions == 0 (straddle count)" % L,
         all(r[3] == 0 for r in sep_rows))
    gate("A L=%d patches: IR(AuB) constant on separated placements (contact-only step)" % L,
         len(set(r[8] for r in sep_rows)) == 1 and
         all(r[8] != sep_rows[0][8] for r in con_rows),
         "sep %s contact %s" % (sorted(set(r[8] for r in sep_rows)), sorted(set(r[8] for r in con_rows))))
    add_seq("patch MI (L=%d, both states)" % L, "INDUCED",
                      [(r[2], r[7]) for r in sep_rows])
    add_seq("patch straddle IR (L=%d)" % L, "INDUCED",
                      [(r[2], r[3]) for r in sep_rows])
    add_seq("patch IR(AuB) (L=%d)" % L, "INDUCED",
                      [(r[2], r[8]) for r in sep_rows])
    add_seq("patch S_A (L=%d)" % L, "INDUCED",
                      [(r[2], r[4]) for r in sep_rows])

    # D-22: translation invariance of a mid-separation row
    Bmid = T.patch(4, 0, 2)
    d1 = graph_dist_regions(adj, A, Bmid)
    mi1 = mutual_info(states["Zbar-fixed"], n, n, A, Bmid)[3]
    A2 = T.patch(3, 2, 2); B2 = T.patch(7, 2, 2)
    d2 = graph_dist_regions(adj, A2, B2)
    mi2 = mutual_info(states["Zbar-fixed"], n, n, A2, B2)[3]
    gate("A L=%d D-22 translation: shifted placement gives same (d_G, MI)" % L,
         d1 == d2 and mi1 == mi2, "(%s,%s) vs (%s,%s)" % (d1, mi1, d2, mi2))

    # ---- band family (non-contractible regions)
    print("\nA.%d two width-2 non-contractible bands on the (%d,%d) torus" % (L, L, L))
    Ab = T.band(0, 2)
    rows_band = []
    for x0 in range(2, L - 1):
        Bb = T.band(x0, 2)
        if set(Ab) & set(Bb):
            continue
        d = graph_dist_regions(adj, Ab, Bb)
        st = straddle_count(local, n, Ab, Bb)
        for sname, gens in states.items():
            SA, SB, SAB, MI = mutual_info(gens, n, n, Ab, Bb)
            rows_band.append((x0, sname, d, st, SA, SB, SAB, MI))
    for r in rows_band:
        print("  x0=%-3d %-10s d=%d  str=%d  S_A=%2d S_B=%2d S_AB=%2d  MI=%d" %
              (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]))
    bsep = [r for r in rows_band if r[3] == 0]
    bcon = [r for r in rows_band if r[3] > 0]
    for sname in states:
        vals = sorted(set(r[7] for r in bsep if r[1] == sname))
        gate("A L=%d bands state %s: MI CONSTANT over all separated separations" % (L, sname),
             len(vals) == 1, "constant value = %s (TOPOLOGICAL, exact)" % vals)
        add_seq("band MI (L=%d, %s)" % (L, sname), "INDUCED",
                          [(r[2], r[7]) for r in bsep if r[1] == sname])
    gate("A L=%d bands: D-15 beside the flat line: contact rows exist and differ or equal by computation" % L,
         len(bcon) > 0, "contact MI %s vs separated %s" %
         (sorted(set(r[7] for r in bcon)), sorted(set(r[7] for r in bsep))))

# ---- INSERTED range-r control carrier: the instrument MUST show a falloff where one exists
print("\nA.ctrl  D-15 INSERTED CONTROL: ring carriers with finite-range pair correlations")
print("        (correlations INSERTED by construction; this law is the control's, never Gamma's)")
for (nq, r) in [(24, 3), (20, 5)]:
    starts = [i for i in range(nq) if (i % (2 * r)) < r]
    gens = []
    for i in starts:
        j = (i + r) % nq
        gens.append((1 << i) | (1 << j))                     # XX
        gens.append((1 << (nq + i)) | (1 << (nq + j)))       # ZZ
    G = rank_f2(gens)
    gate("A.ctrl ring nq=%d r=%d: pure (rank == n)" % (nq, r), G == nq)
    wlen = 6 if nq == 24 else 5
    a0 = 4 if nq == 24 else 2
    A = list(range(a0, a0 + wlen))
    seq = []
    for off in range(a0 + wlen, a0 + wlen + 7):
        B = [q % nq for q in range(off, off + wlen)]
        if set(A) & set(B):
            continue
        MI = mutual_info(gens, nq, nq, A, B)[3]
        seq.append((off - (a0 + wlen), MI))
    print("  ring nq=%d r=%d: MI vs window gap: %s" % (nq, r, seq))
    nonconst = len(set(v for _, v in seq)) > 1
    hits0 = any(v == 0 for _, v in seq) and any(v > 0 for _, v in seq)
    dec = classify(seq) == "DECREASING"
    gate("A.ctrl ring nq=%d r=%d: instrument SHOWS a falloff where one exists (decreasing to 0)" % (nq, r),
         nonconst and hits0 and dec)
    add_seq("ring MI (nq=%d r=%d)" % (nq, r), "INSERTED", seq)
ctrl_tails = []
for lbl, org, seq, _cov in SEQUENCES[-2:]:
    ctrl_tails.append(max([d for d, v in seq if v > 0], default=-1))
gate("A.ctrl the two inserted ranges give DIFFERENT falloff tails (instrument resolves scale)",
     ctrl_tails[0] != ctrl_tails[1], "last nonzero gaps %s" % ctrl_tails)

# =====================================================================================
# SECTION B -- ATTEMPT (a): TWO-DEFECT HOMOLOGY LATTICE vs EARNED SEPARATION
# =====================================================================================
print("\n" + "=" * 100)
print("SECTION B -- attempt (a): two punctures (enlarged H_1), homology-lattice coupling vs separation")
print("  Defect record: Z_L = boundary loop of removed plaquette u (= B_u), X_L = dual connector p_uv.")
print("  Owner of the defect construction: Bravyi-Kitaev quant-ph/9811052; carrier Kitaev quant-ph/9707021.")
print("=" * 100)

HOLE_CASES = {
    (4, 6): [(1, 0), (2, 0), (0, 1), (0, 2), (0, 3), (1, 1), (2, 1), (1, 2), (2, 2), (1, 3), (2, 3)],
    (3, 7): [(1, 0), (0, 1), (0, 2), (0, 3), (1, 1), (1, 2), (1, 3)],
}
HOLE_RESULTS = {}   # (venue, v) -> dict with d_gen, wmin, nmin, zmin, nz, gram, MI data

for venue, vlist in HOLE_CASES.items():
    Lx, Ly = venue
    T, stars, plaqs = VENUES[venue]
    n = T.n
    u = (0, 0)
    Bu = T.plaq(*u)
    local = stars + plaqs
    plaq_index = {(x, y): Ly and (y * Lx + x) for y in range(Ly) for x in range(Lx)}
    gen_idx_u = len(stars) + plaq_index[u]
    stars_ind = independent_subset(stars)
    gens_X = stars_ind + [T.xbar1(), T.xbar2()]
    rkX = rank_f2(gens_X)
    gate("B (%d,%d): X-writer span independent (dim %d)" % (Lx, Ly, len(gens_X)), rkX == len(gens_X))
    xg = [g & ((1 << n) - 1) for g in gens_X]   # pure x-part masks for the scan

    print("\nB venue (%d,%d): hole u=%s against %d partner holes" % (Lx, Ly, u, len(vlist)))
    print("  v        d_gen d_q   Gram-equal  w_min(X-conn) N_min  z_min(loop) N_z   MI(ring,ring) Z-state/X-state")
    gram0 = None
    for v in vlist:
        Bv = T.plaq(*v)
        p = T.dual_path_x(u, v)
        rem_plaqs = [T.plaq(x, y) for y in range(Ly) for x in range(Lx) if (x, y) not in (u, v)]
        # B1 admissibility of the connector
        adm = all(sp_pair(p, s, n) == 0 for s in stars) and \
              all(sp_pair(p, q, n) == 0 for q in rem_plaqs) and \
              sp_pair(p, Bu, n) == 1 and sp_pair(p, Bv, n) == 1
        # earned separations (two instruments)
        gen_idx_v = len(stars) + plaq_index[v]
        d_gen = generator_graph_dist(local, n, gen_idx_u, gen_idx_v)
        adjq = qubit_graph(local, n)
        su = [q for q in range(n) if (supp_mask(Bu, n) >> q) & 1]
        sv = [q for q in range(n) if (supp_mask(Bv, n) >> q) & 1]
        d_q = graph_dist_regions(adjq, su, sv)
        # Gram matrix of the class set {Bu, p, zbar1, zbar2, xbar1, xbar2}
        cls = [Bu, p, T.zbar1(), T.zbar2(), T.xbar1(), T.xbar2()]
        gram = tuple(tuple(sp_pair(a, b, n) for b in cls) for a in cls)
        if gram0 is None:
            gram0 = gram
        # exhaustive X-connector coset scan (pure-Gamma writer landscape)
        wmin, nmin, tot = coset_min_np(xg, p & ((1 << n) - 1))
        # exhaustive Z-loop coset scan (flat-capable control on the SAME instrument)
        plq_ind = independent_subset(rem_plaqs)
        gens_Z = [g >> n for g in plq_ind] + [T.zbar1() >> n, T.zbar2() >> n]
        assert rank_f2(gens_Z) == len(gens_Z)
        zmin, nz, totz = coset_min_np(gens_Z, Bu >> n)
        # COMPLETENESS of the two scans (the exhaustiveness claim rests on this):
        #   X side: admissible X-ops = dual chains with boundary in {u,v} (dim LxLy+1
        #   cycles + connector coset); writers flipping the hole-Z record = p * cycles;
        #   the span IS the full dual-cycle group (dim gated) and every generator commutes
        #   with Bu, so coset = ALL such writers.
        #   Z side: writers flipping the hole-X record = Bu * {cycles pairing 0 with p};
        #   that subgroup has dim LxLy (Bu itself pairs 1), the span has dim LxLy and every
        #   generator pairs 0 with p, so coset = ALL such writers.
        complete_X = all(sp_pair(g, Bu, n) == 0 for g in gens_X)
        complete_Z = all(sp_pair(g, p, n) == 0 for g in plq_ind + [T.zbar1(), T.zbar2()])
        # HOMOLOGICAL CROSSING TYPE of the connector class (computed, gauge-invariant on the
        # punctured venue): which torus-wrapping Z-loop classes pair oddly with p.  Any such
        # loop is itself an admissible writer of the hole record, so the local-record cost
        # is predicted to SWITCH: z_pred = min(4, weights of odd-pairing wrap loops).
        odd_wrap_weights = []
        for r_ in range(Ly):          # horizontal wrap loop at row r_, weight Lx
            loop = 0
            for x_ in range(Lx):
                loop |= 1 << (n + T.h(x_, r_))
            if sp_pair(p, loop, n) == 1:
                odd_wrap_weights.append(Lx)
        for c_ in range(Lx):          # vertical wrap loop at column c_, weight Ly
            loop = 0
            for y_ in range(Ly):
                loop |= 1 << (n + T.v(c_, y_))
            if sp_pair(p, loop, n) == 1:
                odd_wrap_weights.append(Ly)
        ctype = tuple(sorted(set(odd_wrap_weights)))   # e.g. () or (3,) or (3,7)
        z_pred = min([4] + odd_wrap_weights)
        # hole-annulus MI in both defect states
        ring_u, ring_v = su, sv
        mi_pair = []
        for extra, tag in [(Bu, "Zhole"), (p, "Xhole")]:
            gens_state = stars + rem_plaqs + [T.zbar1() if tag == "Zhole" else T.xbar1(),
                                              T.zbar2() if tag == "Zhole" else T.xbar2(), extra]
            Gs = rank_f2(gens_state)
            assert Gs == n, (venue, v, tag, Gs)
            if set(ring_u) & set(ring_v):
                mi_pair.append(None)
            else:
                mi_pair.append(mutual_info(gens_state, n, n, ring_u, ring_v)[3])
        st = straddle_count(local, n, su, sv)
        HOLE_RESULTS[(venue, v)] = dict(d_gen=d_gen, d_q=d_q, adm=adm, gram=gram,
                                        wmin=wmin, nmin=nmin, zmin=zmin, nz=nz,
                                        mi=mi_pair, straddle=st, tot=tot,
                                        ctype=ctype, z_pred=z_pred, disp=v,
                                        complete=(complete_X and complete_Z))
        print("  %-8s  %2d   %2s   %-5s        %2d        %5d      %2d      %4d    %s  ctype=%s z_pred=%d"
              % (str(v), d_gen, str(d_q), gram == gram0, wmin, nmin, zmin, nz, mi_pair,
                 ctype, z_pred))
    res = [HOLE_RESULTS[(venue, v)] for v in vlist]
    gate("B (%d,%d): connector admissible and record-conjugate for every pair" % (Lx, Ly),
         all(r["adm"] for r in res))
    gate("B (%d,%d): scan COMPLETENESS: both coset spans equal the full writer subgroups "
         "(every span generator commutes with the record operator; dims match)" % (Lx, Ly),
         all(r["complete"] for r in res))
    gate("B (%d,%d): homology-lattice Gram matrix IDENTICAL across ALL separations (TOPOLOGICAL)" % (Lx, Ly),
         all(r["gram"] == gram0 for r in res))
    gate("B (%d,%d): two earned separation instruments agree: d_gen == d_qubit + 1 on every pair" % (Lx, Ly),
         all(r["d_gen"] == r["d_q"] + 1 for r in res))
    sepr = [r for r in res if r["straddle"] == 0 and r["mi"][0] is not None]
    for k, tag in [(0, "Z-hole-fixed"), (1, "X-hole-fixed")]:
        vals = sorted(set(r["mi"][k] for r in sepr))
        gate("B (%d,%d): hole-annulus MI (%s) CONSTANT over separated pairs" % (Lx, Ly, tag),
             len(vals) == 1, "value %s over d_gen %s" % (vals, sorted(set(r["d_gen"] for r in sepr))))
        add_seq("hole-annulus MI (%s,%s)" % (venue, tag), "INDUCED",
                          [(r["d_gen"], r["mi"][k]) for r in sepr])
    add_seq("hole Gram entry (u-conn pairing) %s" % (venue,), "INDUCED",
                      [(r["d_gen"], r["gram"][0][1]) for r in res])

# =====================================================================================
# SECTION C -- ATTEMPT (c): THE WRITER-WEIGHT LANDSCAPE (exhaustive, exact)
# =====================================================================================
print("\n" + "=" * 100)
print("SECTION C -- attempt (c): minimal admissible-writer weight coupling record A to record B")
print("  vs earned separation.  C-78 instrument (writer weight = boundary-crossing cost), exhaustive scans.")
print("=" * 100)

for venue, vlist in HOLE_CASES.items():
    Lx, Ly = venue
    print("\nC venue (%d,%d):  d_gen | w_min(coupling writer) N_min | z_min(local record) N_z" % (Lx, Ly))
    rows = []
    for v in vlist:
        r = HOLE_RESULTS[(venue, v)]
        rows.append((v, r["d_gen"], r["wmin"], r["nmin"], r["zmin"], r["nz"], r["tot"],
                     r["ctype"], r["z_pred"]))
    for v, d, w, nm, z, nz, tot, ct, zp in sorted(rows, key=lambda t: t[1]):
        print("   v=%-7s  d=%d  |  w_min=%2d  N_min=%5d  |  z_min=%d (pred %d, ctype %s)  N_z=%d"
              "   (coset size 2^%d exhaustive)"
              % (str(v), d, w, nm, z, zp, ct, nz, tot.bit_length() - 1))
    gate("C (%d,%d): EXACT LAW  w_min == d_gen  on EVERY pair (linear, slope 1, intercept 0, earned units)"
         % (Lx, Ly), all(w == d for _, d, w, *_ in rows))
    # THE LOCAL-RECORD COST: an adversarial find that is STILL topological, not a falloff.
    # z_min is NOT always 4: any torus-wrapping loop whose class pairs oddly with the
    # connector also writes the hole record, so z_min = min(4, weight of such a wrap class).
    # That is a function of the connector's HOMOLOGICAL CROSSING TYPE (a discrete topological
    # datum), not of the separation magnitude: gated three ways below.
    gate("C (%d,%d): local-record cost z_min == min(4, odd-crossing wrap-class weight) EXACTLY "
         "on every pair (topological switch formula)" % (Lx, Ly),
         all(z == zp for *_, z, nz, tot, ct, zp in rows))
    byct = defaultdict(set)
    for *_, z, nz, tot, ct, zp in rows:
        byct[ct].add(z)
    gate("C (%d,%d): z_min CONSTANT within each crossing type (d-independent at fixed type: "
         "TOPOLOGICAL, not a magnitude law)" % (Lx, Ly),
         all(len(s) == 1 for s in byct.values()),
         "type->value %s" % dict((k, sorted(v)) for k, v in byct.items()))
    gate("C (%d,%d): D-15 beside the flat/topological line: the SAME instrument returned the "
         "LINEAR law on the coupling class of the SAME venue (flat vs growing distinguished)"
         % (Lx, Ly), all(z == zp for *_, z, nz, tot, ct, zp in rows) and
         all(w == d for _, d, w, *_ in rows))
    bydir = defaultdict(set)
    for v, d, w, *_ in rows:
        bydir[d].add(w)
    gate("C (%d,%d): isotropy in the earned metric: equal d_gen => equal w_min across orientations" % (Lx, Ly),
         all(len(s) == 1 for s in bydir.values()))
    # direction-family covariate for the degeneracy count (descriptive: the venue's
    # translation classes, D-22/C-74; the headline claims do not rest on it)
    def dirfam(v):
        return "axis-x" if v[1] == 0 else ("axis-y" if v[0] == 0 else "dx=%d" % v[0])
    add_seq("coupling-writer w_min %s" % (venue,), "INDUCED",
                      [(d, w) for _, d, w, *_ in rows])
    add_seq("coupling-writer degeneracy N_min %s" % (venue,), "INDUCED",
                      [(d, nm) for _, d, _, nm, *_ in rows],
                      cov=[(dirfam(v), d, nm) for v, d, _, nm, *_ in rows])
    add_seq("local-record z_min %s" % (venue,), "INDUCED",
                      [(d, z) for _, d, _, _, z, *_ in rows],
                      cov=[(ct, d, z) for _, d, _, _, z, nz, tot, ct, zp in rows])
    fams = defaultdict(list)
    for v, d, _, nm, *_ in rows:
        fams[dirfam(v)].append((d, nm))
    fam_ok = all(classify(rw) in ("CONSTANT", "INCREASING") for rw in fams.values())
    nvals = defaultdict(list)
    for _, d, _, nm, *_ in rows:
        nvals[d].append(nm)
    dsorted = sorted(nvals)
    gate("C (%d,%d): writer DEGENERACY N_min: within every direction family CONSTANT or GROWING "
         "with d, and its ceiling grows (entropy-like; NOT a falloff; NOT a function of d alone)"
         % (Lx, Ly), fam_ok and max(nvals[dsorted[-1]]) > max(nvals[dsorted[0]]),
         "by family " + str(dict((k, sorted(v_)) for k, v_ in fams.items())))

# ---- C reduction control: full-Pauli exhaustive on (3,3) -- z-parts never lower the minimum
print("\nC.red  FULL-PAULI REDUCTION CONTROL on (3,3), holes u=(0,0), v=(2,1)")
T, stars, plaqs = VENUES[(3, 3)]
n = T.n
u, v = (0, 0), (2, 1)
Bu, Bv = T.plaq(*u), T.plaq(*v)
p = T.dual_path_x(u, v)   # note: wraps are allowed; path built by +1 steps
rem = [T.plaq(x, y) for y in range(3) for x in range(3) if (x, y) not in (u, v)]
stars_ind = independent_subset(stars)
rem_ind = independent_subset(rem)
# full centralizer elements commuting with Bu (X-record writer coset rep p)
span_full_X = stars_ind + rem_ind + [T.zbar1(), T.zbar2(), T.xbar1(), T.xbar2(), Bu]
assert rank_f2(span_full_X) == len(span_full_X)
gate("C.red span checks: every span element commutes with Bu; rep p anticommutes",
     all(sp_pair(g, Bu, n) == 0 for g in span_full_X) and sp_pair(p, Bu, n) == 1)
wfull, nfull, totf = coset_min_np(span_full_X, p, weight="xz", n=n)
xg33 = [g & ((1 << n) - 1) for g in stars_ind + [T.xbar1(), T.xbar2()]]
wx33, nx33, _ = coset_min_np(xg33, p & ((1 << n) - 1))
local33 = stars + plaqs
pidx = {(x, y): y * 3 + x for y in range(3) for x in range(3)}
d33 = generator_graph_dist(local33, n, len(stars) + pidx[u], len(stars) + pidx[v])
print("  full-Pauli coset 2^%d exhaustive: w_min=%d;  X-only: w_min=%d;  d_gen=%d"
      % (totf.bit_length() - 1, wfull, wx33, d33))
gate("C.red full-Pauli min == X-only min == d_gen (z-parts never lower the coupling cost; "
     "supp(P) = supp_x u supp_z >= supp_x, verified exhaustively)", wfull == wx33 == d33)
span_full_Z = stars_ind + rem_ind + [T.zbar1(), T.zbar2(), T.xbar1(), T.xbar2(), p]
assert rank_f2(span_full_Z) == len(span_full_Z)
gate("C.red span checks (Z side): every span element commutes with p; rep Bu anticommutes",
     all(sp_pair(g, p, n) == 0 for g in span_full_Z) and sp_pair(Bu, p, n) == 1)
zfull, _, totz = coset_min_np(span_full_Z, Bu, weight="xz", n=n)
odd33 = []
for r_ in range(3):
    loop = 0
    for x_ in range(3):
        loop |= 1 << (n + T.h(x_, r_))
    if sp_pair(p, loop, n) == 1:
        odd33.append(3)
for c_ in range(3):
    loop = 0
    for y_ in range(3):
        loop |= 1 << (n + T.v(c_, y_))
    if sp_pair(p, loop, n) == 1:
        odd33.append(3)
zp33 = min([4] + odd33)
print("  full-Pauli Z-record coset 2^%d exhaustive: z_min=%d (topological-switch prediction %d)"
      % (totz.bit_length() - 1, zfull, zp33))
gate("C.red full-Pauli local-record min == topological-switch prediction (full Pauli space "
     "confirms the wrap-class switch; mixed z/x parts never beat it)", zfull == zp33)

# ---- C INSERTED shortcut control: a wormhole flattens the law (and is labeled INSERTED)
print("\nC.worm  D-15 INSERTED SHORTCUT CONTROL on (4,6): remove a THIRD plaquette (1,0), making the")
print("        single edge shared by cells (0,0)-(1,0) an admissible weight-1 coupling writer.")
print("        The flat law below is the INSERTED carrier's, never Gamma's.")
T, stars, plaqs = VENUES[(4, 6)]
n = T.n
u = (0, 0); extra_hole = (1, 0)
Bu = T.plaq(*u)
worm_rows = []
for v in [(2, 0), (0, 3), (2, 3)]:
    rem = [T.plaq(x, y) for y in range(6) for x in range(4) if (x, y) not in (u, v, extra_hole)]
    stars_ind = independent_subset(stars)
    bridge = T.dual_path_x(extra_hole, v)          # boundary {extra_hole, v}
    rep = 1 << T.v(1, 0)                            # single edge: boundary {u, extra_hole}
    adm = all(sp_pair(rep, s, n) == 0 for s in stars) and \
          all(sp_pair(rep, q, n) == 0 for q in rem) and sp_pair(rep, Bu, n) == 1
    gens = [g & ((1 << n) - 1) for g in stars_ind + [T.xbar1(), T.xbar2()]] + \
           [bridge & ((1 << n) - 1)]
    assert rank_f2(gens) == len(gens)
    wmin, nm, tot = coset_min_np(gens, rep)
    d_gen = HOLE_RESULTS[((4, 6), v)]["d_gen"]
    worm_rows.append((v, d_gen, wmin, adm))
    print("  v=%-7s d_gen(u,v)=%d  ->  w_min with wormhole = %d  (coset 2^%d exhaustive)"
          % (str(v), d_gen, wmin, tot.bit_length() - 1))
gate("C.worm wormhole writer admissible on the modified carrier", all(a for *_, a in worm_rows))
gate("C.worm INSERTED shortcut makes the coupling cost separation-INDEPENDENT (constant 1): the "
     "instrument distinguishes flat from linear; the pure-Gamma carrier chose LINEAR",
     len(set(w for _, _, w, _ in worm_rows)) == 1 and worm_rows[0][2] == 1 and
     len(set(d for _, d, _, _ in worm_rows)) > 1)
add_seq("wormhole w_min (INSERTED carrier)", "INSERTED",
                  [(d, w) for _, d, w, _ in worm_rows])

# =====================================================================================
# SECTION D -- THE FALLOFF SEARCH VERDICT (computed over every measured sequence)
# =====================================================================================
print("\n" + "=" * 100)
print("SECTION D -- systematic falloff search over EVERY measured two-region sequence")
print("  (separated regime d >= separation threshold; contact steps reported as contact-only)")
print("=" * 100)
def seq_verdict(rows, cov):
    """Exact per-sequence verdict.
    A FALLOFF requires a dependence that DECREASES with earned separation: either the
    value is a function of d and classifies DECREASING, or -- when the value is not a
    function of d alone -- some computed covariate class shows a decreasing-in-d
    dependence.  A value that is not a function of d and whose covariate classes are each
    constant is TOPOLOGICAL: set by a discrete datum, not by the separation magnitude."""
    if is_function_of_d(rows):
        c = classify(rows)
        return c, (c == "DECREASING")
    if cov is None:
        return "NOT-f(d), UNRESOLVED", True   # conservatively counts toward falloff-found
    grp = defaultdict(list)
    for k, d, v_ in cov:
        grp[k].append((d, v_))
    sub = {k: classify(rw) for k, rw in grp.items()}
    if any(c == "DECREASING" or c == "MIXED" for c in sub.values()):
        return "NOT-f(d), class DEC/MIXED %s" % sub, True
    if all(c == "CONSTANT" for c in sub.values()):
        return "TOPOLOGICAL (const per class)", False
    return "GROWING per class", False

falloff_induced = []
induced_verdicts, inserted_verdicts = [], []
print("  %-46s %-9s %-34s rows" % ("sequence", "origin", "verdict"))
for label, origin, rows, cov in SEQUENCES:
    verd, is_falloff = seq_verdict(rows, cov)
    print("  %-46s %-9s %-34s %s" % (label, origin, verd,
                                     rows if len(rows) <= 8 else rows[:8] + ["..."]))
    if origin == "INDUCED":
        induced_verdicts.append(verd)
        if is_falloff:
            falloff_induced.append(label)
    else:
        inserted_verdicts.append(verd)
gate("D VERDICT: NO pure-Gamma (INDUCED) quantity decreases with earned separation, globally "
     "or within any computed covariate class -- no falloff was constructible",
     len(falloff_induced) == 0, "falloff-shaped induced sequences: %s" % falloff_induced)
gate("D every INDUCED sequence is CONSTANT, TOPOLOGICAL (constant per computed discrete "
     "class), or GROWING -- and the growing ones exist (linear writer cost, degeneracy)",
     all(v in ("CONSTANT", "INCREASING", "TOPOLOGICAL (const per class)", "GROWING per class")
         for v in induced_verdicts) and "INCREASING" in induced_verdicts)
gate("D the INSERTED controls (finite-range carrier, wormhole carrier) are exactly where the "
     "decreasing and shortcut-flat laws live (mediator-shaped; labeled, never Gamma's)",
     "DECREASING" in inserted_verdicts and "CONSTANT" in inserted_verdicts)

print("\nSTRUCTURAL OBSERVATION (venue-scoped, corner tier): on a stabilizer venue every Pauli")
print("expectation is -1, 0 or +1 (P in +/-S or not), so no continuously decaying two-point")
print("profile exists to be found; a pure-Gamma falloff could only appear as a decreasing")
print("INTEGER sequence in an aggregate.  Every induced integer aggregate measured above is")
print("constant, contact-only, or increasing.  This scopes the verdict to the corner tier;")
print("it does not preclude other venues, and it is why the adversarial attempt ranged over")
print("aggregates (MI, IR, Gram, writer cost, writer degeneracy) rather than correlators.")

# =====================================================================================
# GATE SUMMARY, D-24 AUDIT, FINAL COMPARISON
# =====================================================================================
npass = sum(1 for _, ok in GATES if ok)
print("\n" + "=" * 100)
print("GATES: %d/%d passed" % (npass, len(GATES)))
ALL = npass == len(GATES)
print("OVERALL %s  (computed conjunction of %d gates)" % ("PASS" if ALL else "FAIL", len(GATES)))
for name, ok in GATES:
    if not ok:
        print("  FAILED: " + name)

print("""
====================================================================================================
D-24 AUDIT TABLE (every concept used earns its place)
----------------------------------------------------------------------------------------------------
  region (qubit subset)   INSTRUMENT       chosen families (patches, bands, hole annuli); results
                                           gated over ALL placements, not a favored one
  separation              EARNED           BFS on the interaction/generator-overlap graph computed
                                           from generator SUPPORTS alone -- no coordinates imported;
                                           two independent instruments (qubit graph, generator graph)
                                           agree exactly (d_gen = d_q + 1, gated); C-78 lineage
  writer weight           EARNED (C-78)    minimal admissible-writer weight, exhaustive coset scans;
                                           reused instrument, cited, exhaustiveness preserved
  entropy of a region     BORROWED         S(A) = |A| - dim S_A, Fattal-Cubitt-Yamamoto-Bravyi-Chuang
                                           quant-ph/0406168; used as instrument, exact integers
  mutual information      BORROWED         instrument combination; topological MI owners:
                                           Kitaev-Preskill 2006, Levin-Wen 2006; stabilizer perimeter
                                           law: Hamma-Ionicioiu-Zanardi 2005
  interface rank          EARNED (C-79)    straddle count between regions; reused instrument, cited
  defect/hole records     BORROWED         Bravyi-Kitaev quant-ph/9811052 (construction); carrier
                                           Kitaev quant-ph/9707021; OURS: the writer-landscape law
                                           w_min = d_gen gated exhaustively, and its degeneracy count
  homology Gram matrix    INDUCED          symplectic pairings of computed classes; exact
  dimension               NOT NEEDED       no dimension claim made in this lane
  Euclidean coordinates   NOT USED         every distance is BFS on computed graphs
  linear law (w = d)      EARNED RESULT    exhaustive, exact, isotropy and translation gated
  classical gravity       NOT USED         no gravitational form in any construction step (D-1);
                                           the comparison target appears ONLY in the final section
  ground-state choice     DECLARED         both Zbar-fixed and Xbar-fixed states run; both hole
                                           states run; every gate quantified over the choice
  INSERTED controls       LABELED          finite-range ring carrier, wormhole carrier: their laws
                                           are reported as the controls', never as Gamma's
----------------------------------------------------------------------------------------------------
LARGEST EXACT OBJECTS: exhaustive coset scans to 2^25 = 33,554,432 admissible writers per hole pair
on the (4,6) torus (48 qubits), exact popcount minima with full degeneracy counts; full-Pauli
exhaustive 2^20 on (3,3); F_2 rank/entropy computations to n = 288 qubits (L = 12 torus).
No floats on any measurement path.
""")

print("=" * 100)
print("FINAL COMPARISON SECTION -- the only place the comparison target is named")
print("=" * 100)
print("""
The comparison target is the Newton-form falloff (phi ~ 1/r, F ~ 1/r^2).

WHAT THE ADVERSARIAL ATTEMPT FOUND (all exact, all gated):
  (a) HOMOLOGY-LATTICE COUPLING: the Gram/intersection structure of the two-defect class
      lattice is IDENTICAL at every earned separation -- TOPOLOGICAL.  No falloff.
  (b) TWO-REGION AGGREGATES: interface rank between separated regions is EXACTLY ZERO
      (contact-only step); patch MI is EXACTLY ZERO on every separated placement (both
      states, both venues); band MI is an exact separation-independent CONSTANT
      (non-contractible pair) -- TOPOLOGICAL.  No falloff.
  (c) WRITER-WEIGHT LANDSCAPE: the minimal admissible coupling-writer weight obeys the
      EXACT LAW  w_min = d  (earned separation, slope 1, intercept 0, every pair, every
      orientation, exhaustive to 2^25) -- a genuine pure-Gamma separation law with the
      OPPOSITE character of a falloff: the cost GROWS linearly with separation
      (string-tension/confinement character; lattice-gauge string picture owners:
      Wegner 1971, Kitaev quant-ph/9707021).  Its degeneracy N_min GROWS with separation
      (entropy-like).  Nothing decays.
  (d) THE ADVERSARIAL FIND OF THE LANE: the LOCAL-record writer cost is not always the
      hole-boundary cost 4 -- it SWITCHES to the venue's wrap-class weight exactly when
      the connector's homological crossing type contains that class (z_min = min(4,
      odd-crossing wrap weight), gated exactly on every pair and in the full-Pauli
      space).  This is the closest thing to a separation effect pure Gamma produced,
      and it is a TWO-VALUED function of a DISCRETE homological datum -- constant in d
      within each crossing type, with both values already present at d = 1.  Topological
      in precisely the hypothesis's sense; not a falloff.

HYPOTHESIS VERDICT: the O-54 hypothesis SURVIVES this adversarial attempt, on the corner
tier, exhaustively: pure-Gamma two-region relations came out topological, contact-only, or
LINEARLY GROWING -- never a falloff.  The falloff remains the mediator's to give (the world
tier's measured exponent-3 law, T42_B, is the coupling's contribution and stays labeled as
such).  The division of labor C-77 would take as its refined shape -- Gamma gives the
boundary-degree law and a linear (confinement-character) separation cost; the coupling gives
the falloff; alpha gives the costs -- is CONSISTENT with everything measured here.

HONEST LIMIT: this is a search over one carrier family (toric venues to 48/288 qubits),
puncture-realized handles rather than closed higher-genus surfaces, and integer stabilizer
aggregates.  A pure-Gamma falloff on some other venue is not excluded by this lane; what is
excluded is its constructibility from these exhaustively-scanned families.

NEXT STEP (named, per program rule): the linear law w_min = d is Gamma's OWN two-region
relation and it composes with C-79's boundary-degree law.  The named next increment: couple
the linear writer-cost landscape to the mediator (census access/coupling term) and derive
which falloff exponents are REACHABLE through a mediator constrained by Gamma's linear cost
-- i.e., does confinement-grade cost + openness-clause coupling FORCE the mediator's
exponent family, and does the genus-proper (closed higher-genus) venue preserve the
linear law?
""")
print("total runtime %.1fs" % (time.time() - T0))
