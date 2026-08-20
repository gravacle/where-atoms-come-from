"""T43-A  EXTERNALLY CERTIFIABLE CONTENT -- CORNER TIER.  Exact F_2, toric code L up to 12.

THE PRE-REGISTERED QUESTION (REGISTER_V001.md, final section, registered BEFORE this ran):
of the volume-many records a region stores, how many can the outside world independently
CERTIFY (learn the value of) or WRITE, through the interface?  Boundary-bounded in both
tiers -> the surface has answered (certifiability reading).  Volume-bounded -> the surface
declines.  This lane computes the CORNER tier.  The lane returns DATA; the rule decides.

OPERATIONAL DEFINITIONS (D-24: defined ON THE SURFACE, not imported):

  STORED (baseline, T42_C's capP, recomputed):  capP = 2|R| - r_in = log2 of the number of
      operation classes supported in R with distinct effect (syndrome + record action);
      r_in = dim of the stabiliser subgroup supported wholly inside R.

  CERTIFY (learn from outside):  an interior quantity g (a Pauli supported in R) is
      CERTIFIABLE iff there is an operator h supported ENTIRELY OUTSIDE R with g+h in S
      (the stabiliser span): measuring h outside yields g's value, because the code state
      fixes the product g*h.  Certifiable classes = (P_R  intersect  (S + P_out)) modulo
      (S intersect P_R).  CERT = its dimension.  Computed TWO independent ways per region:
        route 1 (sector ranks):  CERT = sum over CSS sectors of (aR + aRc - a),
                                 a / aR / aRc = rank of the generators / restricted to R /
                                 restricted to the complement;
        route 2 (direct span intersection): dim(P_R int (S+P_out)) by explicit rank
                                 arithmetic on the assembled span, minus r_in (with r_in
                                 itself re-derived by an explicit interior-element witness
                                 construction, every element verified supported in R).
      CERT_GI = the record-relevant (gauge-invariant) part of the same count: certifiable
      classes that commute with every stabiliser (logical/wrap classes).

  WRITE (flip from outside):  an interior record bit [g] is FLIPPABLE by writer w iff w is
      ADMISSIBLE (commutes with every stabiliser generator: creates no syndrome -- the
      on-surface admissibility of C-74/C-78) and anticommutes with g (then w flips [g]'s
      value; the pairing is class-well-defined because admissible writers commute with S).
        WRITE0   = flippable ledger classes, writers supported STRICTLY OUTSIDE R.
                   (Locality: disjoint supports commute -- must compute to 0.  This zero is
                   a THEOREM OF LOCALITY, not an instrument verdict; the operational write
                   instrument is the next line.)
        WRITE_B  = flippable ledger classes, writers supported on OUTSIDE + THE BOUNDARY
                   LAYER B = (qubits of R lying in the support of at least one straddling
                   generator) -- the qubits of R that the interface constraints touch, a
                   lattice count, earned not imported.
        WRITE_GI_B = the record-relevant part: how many of the region's LOGICAL/record
                   classes flip through the interface layer.

CONVENTIONS: toric lattice, stars, plaquettes, region-as-induced-edge-set IDENTICAL to
LANE_T42_C_BOUNDCAP/t42_corner.py (= LANE_O50_A_ACTION/f2lib.py, the C-74 lane).  A Pauli
is (x|z); CSS split makes every computation decouple into X and Z sectors.  Exact integer /
F_2 bitmask arithmetic; NO floats anywhere in this file.

FIT NOTHING: polynomial degrees are read by constant k-th finite differences on exact
integer sequences (the earned instrument of T42_C); every verdict line is a computed
boolean.  D-15: controls beside every bound, including constructions where the outside CAN
reach volume-many records (comb, scatter) and one where it can reach NONE (g=0 adjacency).
D-1: no classical gravitational form appears in any construction step; the recovery target
is named only in the READING section at the end of the OUT.
"""
import random

# ---------------------------------------------------------------- exact F_2 on bitmasks
def pc(x):
    return bin(x).count("1")

def rank_bits(vecs):
    piv = {}
    r = 0
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                r += 1
                break
    return r

def basis_bits(vecs):
    piv = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                break
    return list(piv.values())

def kernel_basis(rows, cols):
    """kernel of x -> (row . x)_rows over F_2, x supported on `cols` (bit positions)."""
    colvecs = []
    for idx, c in enumerate(cols):
        v = 0
        for i, rw in enumerate(rows):
            if (rw >> c) & 1:
                v |= (1 << i)
        colvecs.append((v, 1 << idx))
    piv = {}
    ker = []
    for v, tag in colvecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                pv, ptag = piv[h]
                v ^= pv
                tag ^= ptag
            else:
                piv[h] = (v, tag)
                break
        if v == 0:
            ker.append(tag)
    out = []
    for tag in ker:
        mask = 0
        for idx, c in enumerate(cols):
            if (tag >> idx) & 1:
                mask |= (1 << c)
        out.append(mask)
    return out

def subgroup_in(gens, inside, full):
    """basis of {sum of a subset of gens : support inside `inside`} -- explicit witness
       construction: eliminate the OUTSIDE restrictions with combination tags; every kernel
       tag's full sum is an element supported inside (asserted).  Independent of the rank-
       difference formula for r_in."""
    outside = full & ~inside
    piv = {}
    tags = []
    for idx, g in enumerate(gens):
        v = g & outside
        t = 1 << idx
        while v:
            h = v.bit_length() - 1
            if h in piv:
                pv, pt = piv[h]
                v ^= pv
                t ^= pt
            else:
                piv[h] = (v, t)
                break
        if v == 0:
            tags.append(t)
    elems = []
    for t in tags:
        e = 0
        i = 0
        while t:
            if t & 1:
                e ^= gens[i]
            t >>= 1
            i += 1
        assert e & outside == 0
        if e:
            elems.append(e)
    return basis_bits(elems)

# ---------------------------------------------------------------- toric carrier (= T42_C / O50A f2lib)
class Toric:
    def __init__(self, L):
        self.L = L
        self.n = 2 * L * L
        self.h = lambda i, j: (i % L) * L + (j % L)
        self.v = lambda i, j: L * L + (i % L) * L + (j % L)
        self.stars, self.plaqs = [], []
        for i in range(L):
            for j in range(L):
                s = 0
                for e in (self.h(i, j), self.h(i, j - 1), self.v(i, j), self.v(i - 1, j)):
                    s |= 1 << e
                self.stars.append(s)
                t = 0
                for e in (self.h(i, j), self.h(i + 1, j), self.v(i, j), self.v(i, j + 1)):
                    t |= 1 << e
                self.plaqs.append(t)
        self.full = (1 << self.n) - 1

    def edges_of_vertexset(self, VS):
        L = self.L
        R = 0
        for (i, j) in VS:
            if ((i % L), ((j + 1) % L)) in VS:
                R |= 1 << self.h(i, j)
            if (((i + 1) % L), (j % L)) in VS:
                R |= 1 << self.v(i, j)
        return R

    def straddle_edges(self, VS):
        L = self.L
        cnt = 0
        for (i, j) in VS:
            for (di, dj, e) in ((0, 1, self.h(i, j)), (1, 0, self.v(i, j)),
                                (0, -1, self.h(i, (j - 1) % L)), (-1, 0, self.v((i - 1) % L, j))):
                if (((i + di) % L), ((j + dj) % L)) not in VS:
                    cnt += 1
        return cnt

# ---------------------------------------------------------------- region analysis (generic venue)
def analyse(GX, GZ, n, VS, R, tag, cls, T=None, toric=True):
    full = (1 << n) - 1
    Rc = full & ~R
    AQ = pc(R)
    PER = T.straddle_edges(VS) if (T is not None and VS is not None) else -1
    cols_R = [e for e in range(n) if (R >> e) & 1]

    # ---- sector ranks
    aX  = rank_bits(list(GX));  aZ  = rank_bits(list(GZ))
    aXR = rank_bits([g & R  for g in GX]);  aZR = rank_bits([g & R  for g in GZ])
    aXRc= rank_bits([g & Rc for g in GX]);  aZRc= rank_bits([g & Rc for g in GZ])
    r_in  = (aX - aXRc) + (aZ - aZRc)
    r_out = (aX - aXR)  + (aZ - aZR)
    STORED = 2 * AQ - r_in                      # T42_C capP, recomputed
    CERT = (aXR + aXRc - aX) + (aZR + aZRc - aZ)          # route 1

    # ---- r_in by explicit interior-witness construction (independent route)
    SxR = subgroup_in(GX, R, full)
    SzR = subgroup_in(GZ, R, full)
    r_in_w = len(SxR) + len(SzR)

    # ---- CERT route 2: direct span intersection dim(P_R int (S + P_out)) - r_in
    #      Pauli (x|z) encoded as x | (z<<n); stars are (s|0), plaquettes (0|p).
    span = [g for g in GX] + [(g << n) for g in GZ]
    for e in range(n):
        if not ((R >> e) & 1):
            span.append(1 << e)            # X_e outside
            span.append(1 << (n + e))      # Z_e outside
    d_span = rank_bits(span)
    allv = list(span)
    for e in cols_R:
        allv.append(1 << e)
        allv.append(1 << (n + e))
    d_all = rank_bits(allv)                # must be 2n (P_R + S + P_out = everything)
    CERT_direct = (2 * AQ + d_span - d_all) - r_in_w

    # ---- record content (T42_C direct method) + gauge-invariant certifiable part
    kerX_R = kernel_basis(list(GZ), cols_R)   # X-ops in R commuting with all Z-gens
    kerZ_R = kernel_basis(list(GX), cols_R)   # Z-ops in R commuting with all X-gens
    bGX = basis_bits(list(GX))
    bGZ = basis_bits(list(GZ))
    capX = rank_bits(bGX + kerX_R) - len(bGX)
    capZ = rank_bits(bGZ + kerZ_R) - len(bGZ)
    cap_record = capX + capZ
    # certifiable gauge-invariant classes: dim(ker int span(G|_R)) - r_in, per sector
    def gi_cert(ker, gens_res, r_in_sector):
        dA = len(ker)
        dB = rank_bits(list(gens_res))
        dAB = rank_bits(ker + list(gens_res))
        return (dA + dB - dAB) - r_in_sector
    CERT_GI = gi_cert(kerX_R, [g & R for g in GX], aX - aXRc) \
            + gi_cert(kerZ_R, [g & R for g in GZ], aZ - aZRc)

    # ---- boundary layer B: qubits of R touched by straddling generators (lattice count)
    Bmask = 0
    for g in list(GX) + list(GZ):
        if (g & R) and (g & Rc):
            Bmask |= g & R
    O  = Rc | Bmask
    O0 = Rc
    cols_O  = [e for e in range(n) if (O  >> e) & 1]
    cols_O0 = [e for e in range(n) if (O0 >> e) & 1]

    # ---- admissible writers (commute with every generator = zero syndrome), and the
    #      flip-pairing rank against the region's ledger
    WX_B = kernel_basis(list(GZ), cols_O)     # X-type writers on outside+layer
    WZ_B = kernel_basis(list(GX), cols_O)     # Z-type writers on outside+layer
    WRITE_B = rank_bits([w & R for w in WX_B]) + rank_bits([w & R for w in WZ_B])
    WX_0 = kernel_basis(list(GZ), cols_O0)
    WZ_0 = kernel_basis(list(GX), cols_O0)
    WRITE0 = rank_bits([w & R for w in WX_0]) + rank_bits([w & R for w in WZ_0])

    # ---- record-relevant writes: pair record-class representatives with writers
    def reps_mod_span(ker, gen_basis):
        piv = {}
        for g in gen_basis:
            v = g
            while v:
                h = v.bit_length() - 1
                if h in piv:
                    v ^= piv[h]
                else:
                    piv[h] = v
                    break
        reps = []
        for v0 in ker:
            v = v0
            while v:
                h = v.bit_length() - 1
                if h in piv:
                    v ^= piv[h]
                else:
                    piv[h] = v
                    reps.append(v0)
                    break
        return reps
    repsX = reps_mod_span(kerX_R, bGX)        # X-type record reps (flip via Z writers)
    repsZ = reps_mod_span(kerZ_R, bGZ)        # Z-type record reps (flip via X writers)
    def pairing_rank(reps, writers):
        rows = []
        for r0 in reps:
            m = 0
            for j, w in enumerate(writers):
                if pc(r0 & w) & 1:
                    m |= 1 << j
            rows.append(m)
        return rank_bits(rows)
    WRITE_GI_B = pairing_rank(repsX, WZ_B) + pairing_rank(repsZ, WX_B)

    return dict(tag=tag, cls=cls, AQ=AQ, PER=PER, r_in=r_in, r_in_w=r_in_w, r_out=r_out,
                STORED=STORED, CERT=CERT, CERT_direct=CERT_direct, d_all=d_all, n2=2 * n,
                cap_record=cap_record, CERT_GI=CERT_GI, BQ=pc(Bmask),
                WRITE0=WRITE0, WRITE_B=WRITE_B, WRITE_GI_B=WRITE_GI_B,
                HIDDEN=STORED - CERT)

# ---------------------------------------------------------------- degree by finite differences
def degree(seq, maxdeg=4):
    """degree of the integer sequence if some k-th difference is identically zero with the
       (k-1)-th a nonzero constant; None if no difference goes constant in the window."""
    cur = list(seq)
    if len(set(cur)) == 1:
        return 0 if cur[0] != 0 else 0
    for k in range(1, maxdeg + 1):
        cur = [b - a for a, b in zip(cur, cur[1:])]
        if len(cur) < 2:
            return None
        if len(set(cur)) == 1:
            return k if cur[0] != 0 else k - 1
    return None

def diffs_str(seq):
    d1 = [b - a for a, b in zip(seq, seq[1:])]
    d2 = [b - a for a, b in zip(d1, d1[1:])]
    d3 = [b - a for a, b in zip(d2, d2[1:])]
    return f"D1={d1} D2={d2} D3={d3}"

def rect_vs(L, a, b, i0, j0):
    return {(((i0 + s) % L), ((j0 + t) % L)) for s in range(a) for t in range(b)}

def gate(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok

# ---------------------------------------------------------------- main
def main():
    ok = True
    print("T43-A  EXTERNALLY CERTIFIABLE CONTENT OF A REGION -- CORNER TIER, EXACT F_2")
    print("=" * 100)

    # ---------- per-L carrier sanity ----------
    toric = {}
    for L in (3, 4, 5, 6, 8, 10, 12):
        T = Toric(L)
        toric[L] = T
        m = L * L
        allst = 0
        for s in T.stars:
            allst ^= s
        allpl = 0
        for p in T.plaqs:
            allpl ^= p
        okL = (allst == 0 and allpl == 0 and
               rank_bits(list(T.stars)) == m - 1 and rank_bits(list(T.plaqs)) == m - 1)
        print(f"L={L:>2}: n={T.n:>3} qubits; rank(stars)=rank(plaqs)=m-1, full products vanish: {okL}")
        ok &= okL

    ALL = []

    def run(T, VS, R, tag, cls):
        r = analyse(T.stars, T.plaqs, T.n, VS, R, tag, cls, T=T)
        r["L"] = T.L
        ALL.append(r)
        return r

    # ---------- C-74 / T42_C baseline census: all rectangles, all offsets, L=3..6 ----------
    for L in (3, 4, 5, 6):
        T = toric[L]
        for a in range(1, L):
            for b in range(1, L):
                for i0 in range(L):
                    for j0 in range(L):
                        VS = rect_vs(L, a, b, i0, j0)
                        R = T.edges_of_vertexset(VS)
                        cls = ("RECT_THICK" if (a >= 2 and b >= 2) else "STRIP")
                        run(T, VS, R, f"{a}x{b}@({i0},{j0})", cls)

    # ---------- square sweeps at large L (offset (0,0)) ----------
    for L in (8, 10, 12):
        T = toric[L]
        for s in range(2, L):
            VS = rect_vs(L, s, s, 0, 0)
            run(T, VS, T.edges_of_vertexset(VS), f"{s}x{s}@(0,0)", "SQ_SWEEP")
        for b in range(2, L):                     # strips at large L
            VS = rect_vs(L, 1, b, 0, 0)
            run(T, VS, T.edges_of_vertexset(VS), f"1x{b}@(0,0)", "STRIP")

    # offset spot-check, s=5 at L=12
    T = toric[12]
    off_rows = []
    for (i0, j0) in ((0, 0), (3, 7), (11, 2), (6, 6), (9, 10)):
        VS = rect_vs(12, 5, 5, i0, j0)
        off_rows.append(run(T, VS, T.edges_of_vertexset(VS), f"5x5@({i0},{j0})", "OFFSET"))

    # ---------- non-contractible bands (C-74 positive control) ----------
    band_rows = []
    for L in (4, 5, 6, 8, 12):
        T = toric[L]
        for a in range(1, L):
            VS = {(i, j) for i in range(a) for j in range(L)}
            band_rows.append(run(T, VS, T.edges_of_vertexset(VS), f"band{a}xL", "BAND"))
        VS = {(j, i) for i in range(2) for j in range(L)}   # transposed spot
        band_rows.append(run(T, VS, T.edges_of_vertexset(VS), "bandLx2", "BAND"))

    # ---------- punctured control: added interface must RAISE the certifiable count ----------
    T = toric[12]
    punct_rows = []
    holes = [(2, 2), (2, 6), (6, 2), (6, 6)]
    base = rect_vs(12, 9, 9, 0, 0)
    for k in range(5):
        VS = base - set(holes[:k])
        punct_rows.append(run(T, VS, T.edges_of_vertexset(VS), f"9x9-{k}holes", "PUNCT"))

    # comb: fully punctured geometry -- interface made maximal
    VS = {(i, j) for i in range(9) for j in (0, 2, 4, 6, 8)}
    comb = run(T, VS, T.edges_of_vertexset(VS), "comb9x5cols", "COMB")

    # ---------- scatter control (deterministic seed), matched to the 7x7 block at L=8 ----------
    T = toric[8]
    VS = rect_vs(8, 7, 7, 0, 0)
    block77 = run(T, VS, T.edges_of_vertexset(VS), "7x7@(0,0)", "BLOCK_MATCH")
    rnd = random.Random(43)
    scat_rows = []
    for k in range(3):
        edges = rnd.sample(range(T.n), block77["AQ"])
        R = 0
        VSs = set()
        for e in edges:
            R |= 1 << e
            if e < 64:
                i, j = divmod(e, 8)
                VSs |= {(i, j), (i, (j + 1) % 8)}
            else:
                i, j = divmod(e - 64, 8)
                VSs |= {(i, j), ((i + 1) % 8, j)}
        scat_rows.append(run(toric[8], VSs, R, f"scatter{k}", "SCATTER"))

    # ---------- g=0 adjacency control: same lattice, adjacency (coupling) removed ----------
    g0_rows = []
    T = toric[8]
    GZ0 = [1 << e for e in range(T.n)]            # fully local single-site constraints
    for s in range(2, 8):
        VS = rect_vs(8, s, s, 0, 0)
        R = T.edges_of_vertexset(VS)
        r = analyse([], GZ0, T.n, VS, R, f"g0_{s}x{s}", "G0_LOCAL", T=T)
        r["L"] = 8
        g0_rows.append(r)
    free_rows = []
    for s in range(2, 8):
        VS = rect_vs(8, s, s, 0, 0)
        R = T.edges_of_vertexset(VS)
        r = analyse([], [], T.n, VS, R, f"free_{s}x{s}", "G0_FREE", T=T)
        r["L"] = 8
        free_rows.append(r)

    # ================================================================ table
    print()
    print("=" * 100)
    print("TABLE (exact; STORED = capP = 2|R|-r_in beside CERT and WRITE in the SAME table,")
    print("as the rule requires.  BQ = boundary-layer qubit count, HID = STORED - CERT)")
    print("=" * 100)
    hdr = f"{'L':>2} {'class':<12} {'tag':<16} {'AQ':>4} {'PER':>4} {'STORED':>6} {'CERT':>5} " \
          f"{'HID':>5} {'capR':>4} {'cGI':>4} {'W0':>3} {'W_B':>4} {'WGI':>4} {'BQ':>4}"
    print(hdr)
    shown = set()
    for r in ALL + g0_rows + free_rows:
        if r["cls"] in ("RECT_THICK", "STRIP") and r["L"] <= 6 and not r["tag"].endswith("@(0,0)"):
            continue
        key = (r["L"], r["cls"], r["tag"])
        if key in shown:
            continue
        shown.add(key)
        print(f"{r['L']:>2} {r['cls']:<12} {r['tag']:<16} {r['AQ']:>4} {r['PER']:>4} "
              f"{r['STORED']:>6} {r['CERT']:>5} {r['HIDDEN']:>5} {r['cap_record']:>4} "
              f"{r['CERT_GI']:>4} {r['WRITE0']:>3} {r['WRITE_B']:>4} {r['WRITE_GI_B']:>4} {r['BQ']:>4}")

    # ================================================================ gates
    print()
    print("=" * 100)
    print("GATES -- every verdict is a computed boolean over the full region list")
    print("=" * 100)
    rect = [r for r in ALL if r["cls"] in ("RECT_THICK", "STRIP") and r["L"] <= 6]
    thick = [r for r in rect if r["cls"] == "RECT_THICK"]
    strips = [r for r in ALL if r["cls"] == "STRIP"]
    sq = {L: [r for r in ALL if r["cls"] == "SQ_SWEEP" and r["L"] == L] for L in (8, 10, 12)}
    contractible = rect + [r for r in ALL if r["cls"] in ("SQ_SWEEP", "OFFSET", "PUNCT", "COMB",
                                                          "BLOCK_MATCH")]

    print("BASELINE REPLICATION (control: the instrument reproduces the sealed T42_C venue):")
    ok &= gate("C-74/T42_C census: rectangle count == 1480 (L=3..6, all shapes, all offsets)",
               len(rect) == 1480, f"count={len(rect)}")
    ok &= gate("cap_record == 0 on EVERY contractible region; > 0 on EVERY band (C-74)",
               all(r["cap_record"] == 0 for r in contractible) and
               all(r["cap_record"] > 0 for r in band_rows))
    a34 = next(r for r in thick if r["L"] == 4 and r["tag"] == "3x3@(0,0)")
    a35 = next(r for r in thick if r["L"] == 5 and r["tag"] == "3x3@(0,0)")
    a24 = next(r for r in thick if r["L"] == 5 and r["tag"] == "2x4@(0,0)")
    ok &= gate("sealed T42_C anchors reproduced: 3x3 has STORED(capP)=19 (L=4 and L=5), "
               "2x4 has 17; their certifiable counts equal the sealed IR2=14",
               a34["STORED"] == 19 and a35["STORED"] == 19 and a24["STORED"] == 17 and
               a34["CERT"] == a35["CERT"] == a24["CERT"] == 14)

    print()
    print("CERTIFIABLE CONTENT -- TWO INDEPENDENT ROUTES AND THE INTERFACE IDENTITY:")
    every = ALL + g0_rows + free_rows
    ok &= gate("CERT (sector ranks) == CERT (direct span intersection), EVERY region incl. controls",
               all(r["CERT"] == r["CERT_direct"] for r in every))
    ok &= gate("r_in rank-formula == r_in explicit-witness construction, EVERY region",
               all(r["r_in"] == r["r_in_w"] for r in every))
    ok &= gate("span completeness: P_R + S + P_out has full dimension 2n, EVERY region",
               all(r["d_all"] == r["n2"] for r in every))
    ok &= gate("THE IDENTITY: certifiable content == the T-42 interface law 2*PER-10 on every "
               "thick rectangle, every offset (the interface aggregate IS certifiability)",
               all(r["CERT"] == 2 * r["PER"] - 10 for r in thick +
                   [x for x in ALL if x["cls"] in ("SQ_SWEEP", "OFFSET", "BLOCK_MATCH")]))
    ok &= gate("thin regions hide nothing: CERT == STORED on EVERY strip (all interface)",
               all(r["CERT"] == r["STORED"] for r in strips))
    ok &= gate("thick regions hide bulk: HIDDEN = STORED - CERT > 0 on every thick square s>=3",
               all(r["HIDDEN"] > 0 for L in (8, 10, 12) for r in sq[L] if r["AQ"] >= 12))

    print()
    print("THE LOCALITY THEOREM AND THE WRITE INSTRUMENT:")
    ok &= gate("WRITE0 == 0 on EVERY region (strictly-outside admissible writers flip NOTHING:"
               " disjoint supports commute -- computed, not narrated)",
               all(r["WRITE0"] == 0 for r in every))
    ok &= gate("WRITE_B > 0 on every NON-EMPTY toric region: the boundary layer is where "
               "writing happens",
               all(r["WRITE_B"] > 0 for r in ALL if r["AQ"] > 0))

    print()
    print("EXACT DEGREES BY FINITE DIFFERENCES (squares s x s, offset (0,0); the earned instrument):")
    verdicts = {}
    for L in (10, 12):
        rows = sorted(sq[L], key=lambda r: r["AQ"])
        S_ = [r["STORED"] for r in rows]
        C_ = [r["CERT"] for r in rows]
        W_ = [r["WRITE_B"] for r in rows]
        B_ = [r["BQ"] for r in rows]
        H_ = [r["HIDDEN"] for r in rows]
        G_ = [r["CERT_GI"] for r in rows]
        WG = [r["WRITE_GI_B"] for r in rows]
        print(f"  L={L}: s=2..{L-1}")
        print(f"    STORED {S_}")
        print(f"           {diffs_str(S_)}")
        print(f"    CERT   {C_}")
        print(f"           {diffs_str(C_)}")
        print(f"    W_B    {W_}")
        print(f"           {diffs_str(W_)}")
        print(f"    BQ     {B_}   HIDDEN {H_}")
        print(f"    CERT_GI{G_}   WRITE_GI_B {WG}")
        # W_B has an exact ALL-LAYER CROSSOVER: for s <= 3 the whole square lies in the
        # boundary layer (BQ == AQ) and W_B == STORED; the boundary law holds from s >= 4.
        cross = all(rows[i]["BQ"] == rows[i]["AQ"] and W_[i] == S_[i] for i in (0, 1))
        verdicts[L] = dict(dS=degree(S_), dC=degree(C_), dW=degree(W_[2:]), dB=degree(B_),
                           dH=degree(H_), allGI0=(set(G_) == {0} and set(WG) == {0}),
                           cross=cross,
                           w_law=all(w == c + 6 for w, c in zip(W_[2:], C_[2:])))
        print(f"    degrees: STORED={verdicts[L]['dS']} CERT={verdicts[L]['dC']} "
              f"WRITE_B(s>=4)={verdicts[L]['dW']} BQ={verdicts[L]['dB']} HIDDEN={verdicts[L]['dH']}")
    ok &= gate("STORED is VOLUME-degree: deg_s == 2 at L=10 and L=12",
               all(verdicts[L]["dS"] == 2 for L in (10, 12)))
    ok &= gate("CERT is BOUNDARY-degree: deg_s == 1 at L=10 and L=12",
               all(verdicts[L]["dC"] == 1 for L in (10, 12)))
    ok &= gate("WRITE_B crossover (exact): at s=2,3 the WHOLE square is boundary layer "
               "(BQ == AQ) and W_B == STORED -- a small region is entirely writable",
               all(verdicts[L]["cross"] for L in (10, 12)))
    ok &= gate("WRITE_B is BOUNDARY-degree past the crossover: deg_s == 1 on s>=4 "
               "at L=10 and L=12",
               all(verdicts[L]["dW"] == 1 for L in (10, 12)))
    ok &= gate("exact write law past the crossover: WRITE_B == CERT + 6 (= 8s - 4) on every "
               "square s >= 4 -- the writable and certifiable counts obey the SAME boundary "
               "law, offset by a constant",
               all(verdicts[L]["w_law"] for L in (10, 12)))
    ok &= gate("boundary layer BQ is BOUNDARY-degree: deg_s == 1",
               all(verdicts[L]["dB"] == 1 for L in (10, 12)))
    ok &= gate("HIDDEN (stored minus certifiable) is VOLUME-degree: deg_s == 2",
               all(verdicts[L]["dH"] == 2 for L in (10, 12)))
    ok &= gate("record-relevant counts vanish on contractible squares: CERT_GI == WRITE_GI_B == 0",
               all(verdicts[L]["allGI0"] for L in (10, 12)))
    s4 = {L: next(r for r in sq[L] if r["tag"] == "4x4@(0,0)") for L in (8, 10, 12)}
    ok &= gate("venue-size independence: s=4 square identical (STORED, CERT, WRITE_B, BQ) "
               "at L=8, 10, 12",
               len({(s4[L]["STORED"], s4[L]["CERT"], s4[L]["WRITE_B"], s4[L]["BQ"])
                    for L in (8, 10, 12)}) == 1)
    ok &= gate("translation invariance: 5x5 at 5 offsets, L=12, identical full row",
               len({(r["STORED"], r["CERT"], r["WRITE_B"], r["BQ"], r["WRITE_GI_B"])
                    for r in off_rows}) == 1)

    print()
    print("CONTROL (i) -- D-15: A REGION WHOSE INTERFACE IS MADE ARTIFICIALLY LARGE.")
    print("The boundary-bound must be a MEASUREMENT, not blindness: adding interface must RAISE")
    print("the certifiable count, and a fully punctured geometry must reach volume-many records.")
    pc_ = [r["CERT"] for r in punct_rows]
    ps_ = [r["STORED"] for r in punct_rows]
    print(f"  9x9 block at L=12, k=0..4 deep punctures: CERT={pc_}  STORED={ps_}")
    ok &= gate("CERT rises STRICTLY with every added puncture",
               all(b > a for a, b in zip(pc_, pc_[1:])))
    ok &= gate("puncture increment is the exact interface quantum: +6 per deep hole, every step",
               all(b - a == 6 for a, b in zip(pc_, pc_[1:])))
    ok &= gate("comb (fully punctured): CERT == STORED == 2|R| -- the outside reaches "
               "VOLUME-MANY records when the geometry is all interface",
               comb["CERT"] == comb["STORED"] == 2 * comb["AQ"],
               f"CERT={comb['CERT']} STORED={comb['STORED']} AQ={comb['AQ']}")
    ok &= gate("scatter (matched |R| to the 7x7 block, L=8): CERT_scatter > CERT_block in "
               "every seed -- dispersed content is more certifiable, not less",
               all(r["CERT"] > block77["CERT"] for r in scat_rows),
               f"block CERT={block77['CERT']}, scatter CERT={[r['CERT'] for r in scat_rows]}")

    print()
    print("CONTROL (ii) -- THE NON-CONTRACTIBLE BAND (C-74: record bits ARE reachable there).")
    print(f"  {'L':>2} {'band':<10} {'capR':>4} {'CERT_GI':>7} {'WRITE_GI_B':>10}")
    for r in band_rows:
        print(f"  {r['L']:>2} {r['tag']:<10} {r['cap_record']:>4} {r['CERT_GI']:>7} {r['WRITE_GI_B']:>10}")
    def band_width(r):
        t = r["tag"]
        if t == "bandLx2":
            return 2
        return int(t[4:t.index("x")])
    ok &= gate("every band's record content is certifiable from outside: CERT_GI > 0 on EVERY band",
               all(r["CERT_GI"] > 0 for r in band_rows))
    ok &= gate("thin bands' record bits are WRITABLE through the interface (C-74 reachability "
               "SHOWN): WRITE_GI_B == cap_record on every band of width 1 or 2",
               all(r["WRITE_GI_B"] == r["cap_record"] for r in band_rows if band_width(r) <= 2))
    ok &= gate("THICKNESS PROTECTION, exact profile (reported as found): WRITE_GI_B == 1 at "
               "width 3 (one of the two record bits still reaches the layer), == 0 at every "
               "width >= 4 -- deep record bits are certifiable but NOT writable through the "
               "interface layer",
               all(r["WRITE_GI_B"] == 1 for r in band_rows if band_width(r) == 3) and
               all(r["WRITE_GI_B"] == 0 for r in band_rows if band_width(r) >= 4))

    print()
    print("CONTROL (iii) -- g=0 ADJACENCY: the SAME lattice with the coupling/adjacency removed.")
    print("Fully local single-site constraints (G0_LOCAL) and the free venue (G0_FREE):")
    for fam, rowsF in (("G0_LOCAL", g0_rows), ("G0_FREE", free_rows)):
        S_ = [r["STORED"] for r in rowsF]
        C_ = [r["CERT"] for r in rowsF]
        W_ = [r["WRITE_B"] for r in rowsF]
        print(f"  {fam}: STORED={S_} (deg {degree(S_)}), CERT={C_}, WRITE_B={W_}")
        ok &= gate(f"{fam}: STORED is VOLUME-degree (deg_s == 2) while CERT == 0 and "
                   f"WRITE_B == 0 on every region -- with no straddling constraints the outside "
                   f"certifies NOTHING; certifiability is carried by adjacency",
                   degree(S_) == 2 and set(C_) == {0} and set(W_) == {0})

    print()
    print("D-22 VENUE CHECK: Aut(carrier) = 8L^2 (C-74, ESTABLISHED, borrowed from the register,")
    print("      not recomputed): the venue is far from permutation-symmetric; for L=8,10,12:")
    print("      512, 800, 1152.  Geometry is present to detect.")

    # ================================================================ the rule's input
    print()
    print("=" * 100)
    print("THE PRE-REGISTERED RULE'S INPUT (corner tier) -- computed booleans, no chosen reading")
    print("=" * 100)
    cert_boundary = (all(verdicts[L]["dC"] == 1 for L in (10, 12)) and
                     all(r["CERT"] == r["CERT_direct"] for r in every))
    writ_boundary = (all(r["WRITE0"] == 0 for r in every) and
                     all(verdicts[L]["dW"] == 1 for L in (10, 12)))
    stored_volume = all(verdicts[L]["dS"] == 2 for L in (10, 12))
    instr_sens = (all(b > a for a, b in zip(pc_, pc_[1:])) and
                  comb["CERT"] == 2 * comb["AQ"] and
                  all(set([r["CERT"], r["WRITE_B"]]) == {0} for r in g0_rows + free_rows))
    print(f"  CERTIFIABLE content boundary-degree (deg 1 = s^(D-1), D=2), two routes agreeing: {cert_boundary}")
    print(f"  WRITABLE content: 0 strictly outside; boundary-degree through the interface layer: {writ_boundary}")
    print(f"  STORED content volume-degree (deg 2 = s^D): {stored_volume}")
    print(f"  instrument sensitive both ways (punctures raise it; comb reaches volume; g=0 reaches zero): {instr_sens}")
    corner_boundary_bounded = cert_boundary and writ_boundary and stored_volume and instr_sens
    print(f"  => CORNER TIER: externally certifiable content is BOUNDARY-BOUNDED "
          f"(computed conjunction): {corner_boundary_bounded}")

    print()
    print("READING (the only place the recovery target is named -- D-1): what the outside world")
    print("can learn of a region -- which is what a horizon entropy counts -- scales with the")
    print("boundary; what the region stores scales with the bulk; and the difference (HIDDEN)")
    print("is exactly the volume-degree remainder.  Ownership: that operators in a correctable")
    print("region can be cleaned to its complement, and cut-rank entropy of stabiliser states,")
    print("are BORROWED territory (Bravyi-Terhal cleaning; Kitaev-Preskill / Hamma-Ionicioiu-")
    print("Zanardi).  OURS: the operational certify/write ledger on the record surface, the")
    print("exact identity certifiable==cut-rank by two routes, the locality zero + boundary-")
    print("layer write instrument, the thickness protection of record bits, the puncture/comb/")
    print("g=0 sensitivity controls, and the two-tier test under the pre-registered rule.")

    print()
    n_reg = len(every)
    print(f"REGIONS ANALYSED: {n_reg} (census {len(rect)}, sweeps {sum(len(sq[L]) for L in sq)}, "
          f"strips>{6} {len([r for r in strips if r['L'] > 6])}, bands {len(band_rows)}, "
          f"punctured {len(punct_rows)}+comb, scatter {len(scat_rows)}, offsets {len(off_rows)}, "
          f"g0 {len(g0_rows) + len(free_rows)})")
    print(f"LARGEST EXACT OBJECT: L=12, n=288 qubits, F_2 vectors of width 576; no floats anywhere")
    print()
    print(f"ALL GATES PASS: {ok}")
    return ok

if __name__ == "__main__":
    main()
