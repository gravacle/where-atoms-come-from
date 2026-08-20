"""T42-C  THE BOUNDARY-CAPACITY RELATION -- CORNER TIER.  Exact F_2, toric code L = 3..6.

WHAT IS COMPUTED, per region R (a set of qubits = edges):
  (a) RECORD CONTENT  cap_record = dim( (N(S) cap P_R) / (S cap P_R) )  -- the number of
      independent record-relevant (logical) classes the region supports.  C-74 baseline:
      0 on every contractible region, >0 on every non-contractible one.
  (b) INTERFACE RANK, two exact measures:
        SGEN = COUNT of stabiliser generators whose support straddles the cut (in R and out);
        IR2  = r_S - r_in - r_out  -- the stabiliser CUT-RANK: independent constraints that
               genuinely straddle, modulo products falling wholly inside or outside.
      BORROWED: the cut-rank instrument is Kitaev / Hamma-Ionicioiu-Zanardi territory (region
      entropy of stabiliser states).  OURS: its exact tabulated law on THIS venue's regions,
      the identity LAW-1 below, and the thick/thin crossover.
      Cross-check IR2b = rank of generator restrictions to R minus r_in; must equal IR2.
  (c) DISTINGUISHABLE-CONTENT CAPACITY:
        cap_pauli  = 2|R| - r_in  = log2( # of operation classes in R with distinct effect
                     (syndrome + record action) ),  where r_in = rank of the stabiliser
                     subgroup supported wholly inside R;
        cap_record = the record-relevant part of the same count (see (a)).

CONVENTIONS: toric lattice, stars, plaquettes, and the region-as-induced-edge-set-of-a-
vertex-block convention are IDENTICAL to LANE_O50_A_ACTION/f2lib.py + s4_clause_v.py (C-74).
A Pauli is (x|z) in F_2^{2n}.  Stars are X-type, plaquettes Z-type; the CSS split makes every
computation decouple into an X-sector and a Z-sector, and everything below is exact integer /
F_2 arithmetic on bitmasks.  No floats anywhere in this file.

FIT NOTHING: every 'law' line is an exact equality checked on every region of the named
class by a computed boolean.  Constants in laws are read off one region and then required
to hold exactly on all others INCLUDING the irregular and holed discriminator shapes.
D-15: every law ships with a control class on which the law FAILS (computed, not narrated).
D-1 : no classical gravitational form appears anywhere in the construction; the final
comparison section of the lane report (OUT of t42_world.py) is the only place the recovery
target is named.
"""
import sys, random

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
    """kernel of the F_2 map  x in F_2^{cols} -> (row . x)_rows ; x returned as bitmasks
       supported on `cols` (a list of bit positions)."""
    m = len(rows)
    # augment each column vector with an identity tag to read off kernel combinations
    # we eliminate on the COLUMN vectors: col j is (bits of rows hitting col j | tag j)
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

# ---------------------------------------------------------------- toric carrier (= O50A f2lib)
class Toric:
    """qubits on edges of an L x L periodic lattice.
       h(i,j): edge (i,j)-(i,j+1), index i*L+j ; v(i,j): edge (i,j)-(i+1,j), index L*L+i*L+j.
       star  A(i,j) = X on {h(i,j), h(i,j-1), v(i,j), v(i-1,j)}
       plaq  B(i,j) = Z on {h(i,j), h(i+1,j), v(i,j), v(i,j+1)}"""
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
        """induced edge set: an edge is in R iff BOTH endpoints are in VS."""
        L = self.L
        R = 0
        for (i, j) in VS:
            if ((i % L), ((j + 1) % L)) in VS:
                R |= 1 << self.h(i, j)
            if (((i + 1) % L), (j % L)) in VS:
                R |= 1 << self.v(i, j)
        return R

    def straddle_edges(self, VS):
        """edges with EXACTLY ONE endpoint in VS -- the perimeter COUNT (a lattice count,
           nothing else)."""
        L = self.L
        cnt = 0
        for (i, j) in VS:
            for (di, dj, e) in ((0, 1, self.h(i, j)), (1, 0, self.v(i, j)),
                                (0, -1, self.h(i, (j - 1) % L)), (-1, 0, self.v((i - 1) % L, j))):
                if (((i + di) % L), ((j + dj) % L)) not in VS:
                    cnt += 1
        return cnt

# ---------------------------------------------------------------- region analysis
def analyse(T, VS, R, tag, cls):
    L, n = T.L, T.n
    Rc = T.full ^ R
    AQ = pc(R)
    AV = len(VS)
    PER = T.straddle_edges(VS)
    stars_in = sum(1 for s in T.stars if (s & Rc) == 0)
    plaqs_in = sum(1 for p in T.plaqs if (p & Rc) == 0)
    SGEN = sum(1 for g in T.stars + T.plaqs if (g & R) and (g & Rc))
    # ranks of restrictions
    rk_st_R  = rank_bits([s & R  for s in T.stars])
    rk_st_Rc = rank_bits([s & Rc for s in T.stars])
    rk_pl_R  = rank_bits([p & R  for p in T.plaqs])
    rk_pl_Rc = rank_bits([p & Rc for p in T.plaqs])
    m = L * L
    # subgroup wholly inside / outside (kernel of restriction, minus the 1-dim trivial kernel
    # per sector: the product of ALL stars (resp. plaquettes) has empty support -- verified
    # once per L in main())
    r_in  = (m - rk_st_Rc - 1) + (m - rk_pl_Rc - 1)
    r_out = (m - rk_st_R  - 1) + (m - rk_pl_R  - 1)
    r_S   = rank_bits([s for s in T.stars]) + rank_bits([p for p in T.plaqs])
    IR2   = r_S - r_in - r_out
    IR2b  = (rk_st_R - (m - rk_st_Rc - 1)) + (rk_pl_R - (m - rk_pl_Rc - 1))
    SYN   = rk_st_R + rk_pl_R          # rank of the region-complement coupling map
                                       # P_R -> syndromes (CSS: sectors decouple)
    cap_pauli = 2 * AQ - r_in
    # record content, DIRECT (membership in the full star / plaq span, exactly as C-74 did):
    cols = [e for e in range(n) if (R >> e) & 1]
    kerX = kernel_basis([p for p in T.plaqs], cols)   # X-ops in R commuting with all plaqs
    kerZ = kernel_basis([s for s in T.stars], cols)   # Z-ops in R commuting with all stars
    bst = basis_bits(list(T.stars))
    bpl = basis_bits(list(T.plaqs))
    capX = rank_bits(bst + kerX) - len(bst)
    capZ = rank_bits(bpl + kerZ) - len(bpl)
    cap_record = capX + capZ
    # sector formula for the same quantity (must agree -- computed boolean)
    cap_formula = ((AQ - rk_pl_R) - (m - rk_st_Rc - 1)) + ((AQ - rk_st_R) - (m - rk_pl_Rc - 1))
    # topology counts (pure graph counts on (VS, R))
    comp = components(T, VS, R)
    c1 = AQ - AV + comp
    return dict(L=L, tag=tag, cls=cls, AQ=AQ, AV=AV, PER=PER, SGEN=SGEN,
                stars_in=stars_in, plaqs_in=plaqs_in, r_in=r_in, r_out=r_out, r_S=r_S,
                IR2=IR2, IR2b=IR2b, SYN=SYN, cap_pauli=cap_pauli,
                cap_record=cap_record, cap_formula=cap_formula, comp=comp, c1=c1)

def components(T, VS, R):
    L = T.L
    seen = set()
    comp = 0
    for v0 in VS:
        if v0 in seen:
            continue
        comp += 1
        stack = [v0]
        seen.add(v0)
        while stack:
            (i, j) = stack.pop()
            for (di, dj, e) in ((0, 1, T.h(i, j)), (1, 0, T.v(i, j)),
                                (0, -1, T.h(i, (j - 1) % L)), (-1, 0, T.v((i - 1) % L, j))):
                if (R >> e) & 1:
                    w = (((i + di) % L), ((j + dj) % L))
                    if w in VS and w not in seen:
                        seen.add(w)
                        stack.append(w)
    return comp

# ---------------------------------------------------------------- region generators
def rect_vs(L, a, b, i0, j0):
    return {(((i0 + s) % L), ((j0 + t) % L)) for s in range(a) for t in range(b)}

def main():
    rows = []
    for L in (3, 4, 5, 6):
        T = Toric(L)
        m = L * L
        # once per L: the trivial-kernel assumption behind r_in / r_out
        allst = 0
        for s in T.stars: allst ^= s
        allpl = 0
        for p in T.plaqs: allpl ^= p
        ok_triv = (allst == 0 and allpl == 0 and
                   rank_bits(list(T.stars)) == m - 1 and rank_bits(list(T.plaqs)) == m - 1)
        print(f"L={L}: n={T.n} qubits; rank(stars)=rank(plaqs)=m-1 and full products vanish: {ok_triv}")
        assert ok_triv
        # RECTANGLES, every a,b in 1..L-1, EVERY offset  (the C-74 census, 1480 total)
        for a in range(1, L):
            for b in range(1, L):
                for i0 in range(L):
                    for j0 in range(L):
                        VS = rect_vs(L, a, b, i0, j0)
                        R = T.edges_of_vertexset(VS)
                        cls = ("RECT_THICK" if (a >= 2 and b >= 2) else "STRIP")
                        rows.append(analyse(T, VS, R, f"{a}x{b}@({i0},{j0})", cls))
        # L-SHAPES (offset 0 only; translation acts transitively -- checked on rectangles)
        if L >= 4:
            for (a, b, c, d) in [(3, 3, 1, 1)] + ([(4, 4, 1, 1), (4, 4, 2, 2), (4, 3, 2, 1)] if L >= 5 else []) \
                               + ([(5, 5, 2, 2), (5, 5, 1, 3), (5, 4, 2, 2)] if L >= 6 else []):
                if a <= L - 1 and b <= L - 1:
                    VS = rect_vs(L, a, b, 0, 0) - rect_vs(L, c, d, 0, 0)  # corner removed
                    # remove the corner block at the (0,0) corner
                    VS = rect_vs(L, a, b, 0, 0) - {(s, t) for s in range(c) for t in range(d)}
                    R = T.edges_of_vertexset(VS)
                    rows.append(analyse(T, VS, R, f"L{a}x{b}-{c}x{d}", "LSHAPE"))
        # HOLED: (a x b) block minus interior vertices
        if L >= 4:
            a = b = L - 1
            VS0 = rect_vs(L, a, b, 0, 0)
            mid = (a // 2, b // 2)
            VS = VS0 - {mid}
            R = T.edges_of_vertexset(VS)
            deep = (2 <= mid[0] <= a - 3) and (2 <= mid[1] <= b - 3)
            rows.append(analyse(T, VS, R, f"H{a}x{b}-v{mid}", "HOLED_DEEP" if deep else "HOLED_SHALLOW"))
            if L == 6:
                for extra in [{(2, 2)}, {(1, 2)}, {(2, 2), (2, 4)}, {(1, 1), (3, 3)}]:
                    VS = VS0 - extra
                    R = T.edges_of_vertexset(VS)
                    deep = all((2 <= p <= a - 3) and (2 <= q <= b - 3) for (p, q) in extra)
                    rows.append(analyse(T, VS, R, f"H{a}x{b}-{sorted(extra)}",
                                        "HOLED_DEEP" if deep else "HOLED_SHALLOW"))
        # BANDS (non-contractible; the C-74 positive control)
        for a in range(1, L):
            VS = {(i, j) for i in range(a) for j in range(L)}
            R = T.edges_of_vertexset(VS)
            rows.append(analyse(T, VS, R, f"band{a}xL", "BAND"))
            VS = {(j, i) for i in range(a) for j in range(L)}
            R = T.edges_of_vertexset(VS)
            rows.append(analyse(T, VS, R, f"bandLx{a}", "BAND"))
        # SCATTER control (deterministic seed): random edge sets, |R| matched to the big block
        rnd = random.Random(42 + L)
        target = None
        for r in rows[::-1]:
            if r["L"] == L and r["cls"] == "RECT_THICK":
                target = r["AQ"]; break
        for k in range(3):
            edges = rnd.sample(range(T.n), target)
            R = 0
            for e in edges: R |= 1 << e
            # vertex set = all endpoints (for the counts); region still the edge set R
            VS = set()
            L_ = L
            for e in edges:
                if e < L * L:
                    i, j = divmod(e, L)
                    VS |= {(i, j), (i, (j + 1) % L)}
                else:
                    i, j = divmod(e - L * L, L)
                    VS |= {(i, j), ((i + 1) % L, j)}
            rows.append(analyse(T, None if False else VS, R, f"scatter{k}", "SCATTER"))
    return rows

def gate(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok

if __name__ == "__main__":
    rows = main()
    print()
    print("=" * 100)
    print("FULL TABLE (exact integers; AQ=|R| qubits, AV=vertices, PER=straddling-edge count,")
    print("SGEN=straddling generators, r_in/r_out=stabiliser rank wholly inside/outside,")
    print("IR2=cut-rank, SYN=coupling rank, capP=cap_pauli, capR=cap_record, c1=cycle rank)")
    print("=" * 100)
    hdr = f"{'L':>2} {'class':<13} {'tag':<22} {'AQ':>3} {'AV':>3} {'PER':>4} {'SGEN':>4} " \
          f"{'s_in':>4} {'p_in':>4} {'r_in':>4} {'IR2':>4} {'SYN':>4} {'capP':>5} {'capR':>4} {'c1':>3}"
    print(hdr)
    shown = set()
    for r in rows:
        # print rectangles only at offset (0,0) (translation-duplicates suppressed in the
        # listing, NOT in the checks -- every one of the 1480 is in `rows` and in every gate)
        if r["cls"] in ("RECT_THICK", "STRIP") and not r["tag"].endswith("@(0,0)"):
            continue
        key = (r["L"], r["tag"])
        if key in shown: continue
        shown.add(key)
        print(f"{r['L']:>2} {r['cls']:<13} {r['tag']:<22} {r['AQ']:>3} {r['AV']:>3} {r['PER']:>4} "
              f"{r['SGEN']:>4} {r['stars_in']:>4} {r['plaqs_in']:>4} {r['r_in']:>4} {r['IR2']:>4} "
              f"{r['SYN']:>4} {r['cap_pauli']:>5} {r['cap_record']:>4} {r['c1']:>3}")

    print()
    print("=" * 100)
    print("GATES -- every verdict below is a computed boolean over the full region list")
    print("=" * 100)
    allrows = rows
    rect = [r for r in rows if r["cls"] in ("RECT_THICK", "STRIP")]
    thick = [r for r in rows if r["cls"] == "RECT_THICK"]
    strip = [r for r in rows if r["cls"] == "STRIP"]
    lsh  = [r for r in rows if r["cls"] == "LSHAPE"]
    hd   = [r for r in rows if r["cls"] == "HOLED_DEEP"]
    hs   = [r for r in rows if r["cls"] == "HOLED_SHALLOW"]
    band = [r for r in rows if r["cls"] == "BAND"]
    scat = [r for r in rows if r["cls"] == "SCATTER"]

    ok = True
    ok &= gate("C-74 REPLICATION (baseline control): rectangle census count == 1480",
               len(rect) == 1480, f"count={len(rect)}")
    ok &= gate("C-74 REPLICATION: cap_record == 0 on EVERY contractible region "
               "(rectangles, L-shapes, holed)",
               all(r["cap_record"] == 0 for r in rect + lsh + hd + hs))
    ok &= gate("C-74 POSITIVE CONTROL (D-15): cap_record > 0 on EVERY non-contractible band",
               all(r["cap_record"] > 0 for r in band),
               f"values={sorted(set(r['cap_record'] for r in band))}")
    ok &= gate("cap_record == 2 on every band (both wrap classes of one direction)",
               all(r["cap_record"] == 2 for r in band))
    ok &= gate("cap_record <= 4 everywhere (record space is 4-dimensional, G-9)",
               all(r["cap_record"] <= 4 for r in allrows))
    ok &= gate("cap_record: direct span-membership computation == sector formula, every region",
               all(r["cap_record"] == r["cap_formula"] for r in allrows))
    ok &= gate("CROSS-CHECK: cut-rank two ways, IR2 == IR2b (group quotient == restriction rank), "
               "every region", all(r["IR2"] == r["IR2b"] for r in allrows))
    ok &= gate("COUPLING = CUT-RANK + INTERIOR:  SYN == IR2 + r_in, every region",
               all(r["SYN"] == r["IR2"] + r["r_in"] for r in allrows))
    print()
    print("LAW-1 (OURS; exact, every region of every class incl. scatter):")
    print("      IR2 + 2*r_in + cap_record == 2*AQ")
    print("      (interface rank + twice the interior + the record content exhausts the region's")
    print("       Pauli coordinates -- boundary structure and content share one exact budget)")
    ok &= gate("LAW-1 holds on every region",
               all(r["IR2"] + 2 * r["r_in"] + r["cap_record"] == 2 * r["AQ"] for r in allrows))
    ok &= gate("D-15 CONTROL for LAW-1 instrument: the WRONG identity "
               "IR2 + r_in + cap_record == 2*AQ FAILS somewhere",
               any(r["IR2"] + r["r_in"] + r["cap_record"] != 2 * r["AQ"] for r in allrows))
    print()
    print("LAW-2 (cut-rank law; instrument BORROWED Kitaev/Hamma-Ionicioiu-Zanardi, exact form on")
    print("       these regions OURS):  IR2 == 2*PER - 10   on every THICK contractible shape")
    ok &= gate("LAW-2 on all thick rectangles (a,b>=2), every offset, L=3..6",
               all(r["IR2"] == 2 * r["PER"] - 10 for r in thick))
    ok &= gate("LAW-2 UNCHANGED on L-shapes (irregular discriminator)",
               all(r["IR2"] == 2 * r["PER"] - 10 for r in lsh) if lsh else False)
    ok &= gate("LAW-2 UNCHANGED on deep-holed shapes (holes move PER; law tracks PER exactly)",
               all(r["IR2"] == 2 * r["PER"] - 10 for r in hd) if hd else False)
    ok &= gate("D-15 CONTROL: LAW-2 FAILS on strips (1 x b: no interior)",
               any(r["IR2"] != 2 * r["PER"] - 10 for r in strip))
    ok &= gate("D-15 CONTROL: LAW-2 FAILS on scattered edge-sets of the same size",
               any(r["IR2"] != 2 * r["PER"] - 10 for r in scat))
    ok &= gate("STRIP law (their own exact law, AREA not perimeter): IR2 == 2*AQ on every strip",
               all(r["IR2"] == 2 * r["AQ"] for r in strip))
    ok &= gate("SCATTER law: IR2 == 2*AQ - 2*r_in - cap_record collapses to == 2*AQ "
               "(r_in = 0: nothing is wholly inside)",
               all(r["IR2"] == 2 * r["AQ"] and r["r_in"] == 0 for r in scat))
    print()
    print("INTERIOR ACCOUNTING (what r_in is, earned by count):")
    ok &= gate("r_in == stars_in + plaqs_in on every SOLID contractible shape "
               "(rectangles incl. strips, L-shapes)",
               all(r["r_in"] == r["stars_in"] + r["plaqs_in"] for r in rect + lsh))
    ok &= gate("HOLED shapes: r_in == stars_in + plaqs_in + (#holes) -- each vertex-hole "
               "contributes ONE extra Z-type element (the loop around it), none X-type",
               all(r["r_in"] == r["stars_in"] + r["plaqs_in"] + (r["c1"] - r["plaqs_in"])
                   for r in hd + hs))
    print()
    print("PERIMETER vs AREA -- THE MATCHED-PAIR DISCRIMINATOR (no fit, no exponent, pure pairs):")
    pairs_same_per = []
    pairs_same_aq = []
    cont = thick + lsh + hd
    for i in range(len(cont)):
        for j in range(i + 1, len(cont)):
            a, b = cont[i], cont[j]
            if a["L"] != b["L"]: continue
            if a["PER"] == b["PER"] and a["AQ"] != b["AQ"]:
                pairs_same_per.append((a, b))
            if a["AQ"] == b["AQ"] and a["PER"] != b["PER"]:
                pairs_same_aq.append((a, b))
    ok &= gate(f"SAME perimeter, DIFFERENT area ({len(pairs_same_per)} pairs): IR2 EQUAL in every pair",
               all(a["IR2"] == b["IR2"] for a, b in pairs_same_per) and len(pairs_same_per) > 0)
    ok &= gate(f"SAME perimeter, DIFFERENT area: cap_pauli DIFFERS in every pair "
               "(capacity is NOT a perimeter quantity)",
               all(a["cap_pauli"] != b["cap_pauli"] for a, b in pairs_same_per))
    ok &= gate(f"SAME area, DIFFERENT perimeter ({len(pairs_same_aq)} pairs): IR2 DIFFERS in every pair",
               all(a["IR2"] != b["IR2"] for a, b in pairs_same_aq) and len(pairs_same_aq) > 0)
    ex = pairs_same_per[0] if pairs_same_per else None
    if ex:
        a, b = ex
        print(f"      example: {a['tag']} (AQ={a['AQ']}) vs {b['tag']} (AQ={b['AQ']}), "
              f"PER={a['PER']} both, IR2={a['IR2']}=={b['IR2']}")
    print()
    print("WHAT cap_pauli OBEYS (mixed, bulk-dominant -- stated, not fitted):")
    ok &= gate("cap_pauli == 2*AQ - r_in by construction; on thick rectangles this equals "
               "2*AV + (PER/2) - 5 exactly (checked), i.e. BULK term + boundary term",
               all(2 * r["AV"] + r["PER"] // 2 - 5 == r["cap_pauli"] for r in thick))
    ok &= gate("D-15 CONTROL: the pure-perimeter candidate cap_pauli == 2*PER - c FAILS on thick "
               "rectangles for every constant c (no single c works)",
               len(set(2 * r["PER"] - r["cap_pauli"] for r in thick)) > 1)
    ok &= gate("D-15 CONTROL: the pure-area candidate cap_pauli == 2*AQ - c FAILS "
               "(no single c works)",
               len(set(2 * r["AQ"] - r["cap_pauli"] for r in thick)) > 1)
    print()
    print("D-22 VENUE CHECK: Aut(carrier) = 8L^2 exactly (C-74, ESTABLISHED -- borrowed from the")
    print("      register, not recomputed): the carrier is far from permutation-symmetric; the")
    print("      venue contains geometry to detect.  8L^2 for L=3..6: 72, 128, 200, 288.")
    print()
    n_reg = len(allrows)
    print(f"REGIONS ANALYSED: {n_reg} total "
          f"(rect census {len(rect)}, L-shapes {len(lsh)}, holed {len(hd) + len(hs)}, "
          f"bands {len(band)}, scatter {len(scat)})")
    print(f"LARGEST OBJECT: L=6, n=72 qubits, F_2 vectors of width 144; all arithmetic exact F_2/integer")
    print()
    print(f"ALL GATES PASS: {ok}")
