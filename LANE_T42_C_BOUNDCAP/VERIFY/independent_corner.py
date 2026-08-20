"""VERIFY -- independent reimplementation of the T42-C corner tier.

Written from PUBLISHED_CONVENTIONS.txt, NOT from t42_corner.py: different edge indexing
(interleaved 2*(i*L+j)+dir instead of blocked), different elimination code, different
region-generation code, and a genuinely different assembly of cap_record.
Every verdict is a computed boolean.
"""

# ---------- my own F_2 elimination ----------
def myrank(vs):
    basis = []
    for v in vs:
        for b in basis:
            v = min(v, v ^ b)  # not standard reduction; do proper below
    # proper: rebuild
    piv = {}
    r = 0
    for v in vs:
        w = v
        while w:
            top = w.bit_length() - 1
            if top in piv:
                w ^= piv[top]
            else:
                piv[top] = w
                r += 1
                break
    return r

def myspan_dim_supported(vs, mask):
    """dim of the subspace of span(vs) supported inside `mask` (bits outside mask zero).
       = dim span - rank of the projection to the complement of mask."""
    d = myrank(vs)
    proj = [v & ~mask for v in vs]
    return d - myrank(proj)

def mykernel_dim(rows, cols_mask, ncols_positions):
    """dim of { x supported on cols_mask : for every row r, |x & r| even }.
       computed by my own column elimination."""
    cols = [c for c in range(ncols_positions) if (cols_mask >> c) & 1]
    colvecs = []
    for c in cols:
        v = 0
        for i, rw in enumerate(rows):
            if (rw >> c) & 1:
                v |= 1 << i
        colvecs.append(v)
    return len(cols) - myrank(colvecs)

def mykernel_vectors(rows, cols_mask, ncols_positions):
    """actual kernel vectors (as bitmasks on the original positions)."""
    cols = [c for c in range(ncols_positions) if (cols_mask >> c) & 1]
    colvecs = []
    for idx, c in enumerate(cols):
        v = 0
        for i, rw in enumerate(rows):
            if (rw >> c) & 1:
                v |= 1 << i
        colvecs.append((v, 1 << idx))
    piv = {}
    out = []
    for v, tag in colvecs:
        while v:
            top = v.bit_length() - 1
            if top in piv:
                pv, pt = piv[top]
                v ^= pv
                tag ^= pt
            else:
                piv[top] = (v, tag)
                break
        if v == 0:
            m = 0
            for idx, c in enumerate(cols):
                if (tag >> idx) & 1:
                    m |= 1 << c
            out.append(m)
    return out

# ---------- my own toric carrier: edge (i,j,d) -> index 2*(i*L+j)+d, d=0 horiz, d=1 vert ----------
class MyToric:
    def __init__(self, L):
        self.L = L
        self.n = 2 * L * L
        self.stars = []
        self.plaqs = []
        for i in range(L):
            for j in range(L):
                # star at vertex (i,j): h(i,j), h(i,j-1), v(i,j), v(i-1,j)
                s = 0
                for (a, b, d) in ((i, j, 0), (i, (j - 1) % L, 0), (i, j, 1), ((i - 1) % L, j, 1)):
                    s |= 1 << (2 * (a * L + b) + d)
                self.stars.append(s)
                # plaquette at (i,j): h(i,j), h(i+1,j), v(i,j), v(i,j+1)
                p = 0
                for (a, b, d) in ((i, j, 0), ((i + 1) % L, j, 0), (i, j, 1), (i, ((j + 1) % L), 1)):
                    p |= 1 << (2 * (a * L + b) + d)
                self.plaqs.append(p)

    def e(self, i, j, d):
        L = self.L
        return 2 * ((i % L) * L + (j % L)) + d

    def induced(self, VS):
        L = self.L
        R = 0
        for (i, j) in VS:
            if ((i % L), ((j + 1) % L)) in VS:
                R |= 1 << self.e(i, j, 0)
            if (((i + 1) % L), (j % L)) in VS:
                R |= 1 << self.e(i, j, 1)
        return R

    def straddle(self, VS):
        L = self.L
        c = 0
        for (i, j) in VS:
            if ((i % L), ((j + 1) % L)) not in VS: c += 1
            if ((i % L), ((j - 1) % L)) not in VS: c += 1
            if (((i + 1) % L), (j % L)) not in VS: c += 1
            if (((i - 1) % L), (j % L)) not in VS: c += 1
        return c


def analyse(T, VS, R):
    L, n = T.L, T.n
    m = L * L
    Rc = ((1 << n) - 1) ^ R
    AQ = bin(R).count("1")
    AV = len(VS)
    PER = T.straddle(VS)
    # inside/outside stabiliser subgroup dims, via span-support (my route)
    r_in = myspan_dim_supported(T.stars, R) + myspan_dim_supported(T.plaqs, R)
    r_out = myspan_dim_supported(T.stars, Rc) + myspan_dim_supported(T.plaqs, Rc)
    r_S = myrank(T.stars) + myrank(T.plaqs)
    IR2 = r_S - r_in - r_out
    SYN = myrank([s & R for s in T.stars]) + myrank([p & R for p in T.plaqs])
    cap_pauli = 2 * AQ - r_in
    # cap_record my route: kernel dims minus dim of full-span elements supported in R
    kerX = mykernel_vectors(T.plaqs, R, n)   # X ops in R commuting with all plaquettes
    kerZ = mykernel_vectors(T.stars, R, n)
    inX = myspan_dim_supported(T.stars, R)   # star-span elements supported in R
    inZ = myspan_dim_supported(T.plaqs, R)
    # quotient by FULL-span membership restricted to R: an op in R that lies in span(stars)
    # must lie in the star-span-supported-in-R subspace, so quotient dim:
    capX = myrank(kerX) - inX
    capZ = myrank(kerZ) - inZ
    cap_record = capX + capZ
    # sanity: kernel contains the inside span
    # lattice counts, my own BFS/union-find
    comp = 0
    seen = set()
    for v0 in VS:
        if v0 in seen: continue
        comp += 1
        st = [v0]; seen.add(v0)
        while st:
            (i, j) = st.pop()
            for (di, dj, d, a, b) in ((0, 1, 0, i, j), (0, -1, 0, i, (j - 1) % L),
                                      (1, 0, 1, i, j), (-1, 0, 1, (i - 1) % L, j)):
                if (R >> T.e(a, b, d)) & 1:
                    w = (((i + di) % L), ((j + dj) % L))
                    if w in VS and w not in seen:
                        seen.add(w); st.append(w)
    c1 = AQ - AV + comp
    # dual: plaquettes touching R + components of the dual region graph
    parent = {}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry: parent[rx] = ry
    for idx in range(m):
        pass
    for c in range(n):
        if not ((R >> c) & 1): continue
        owners = [k for k in range(m) if (T.plaqs[k] >> c) & 1]
        assert len(owners) == 2
        for o in owners:
            if o not in parent: parent[o] = o
        union(owners[0], owners[1])
    PD = len(parent)
    comp_D = len({find(x) for x in parent})
    c1_D = AQ - PD + comp_D
    return dict(AQ=AQ, AV=AV, PER=PER, PD=PD, comp=comp, comp_D=comp_D, c1=c1, c1_D=c1_D,
                r_in=r_in, r_out=r_out, IR2=IR2, SYN=SYN, cap_pauli=cap_pauli,
                cap_record=cap_record)


def rect_vs(L, a, b, i0, j0):
    return {(((i0 + s) % L), ((j0 + t) % L)) for s in range(a) for t in range(b)}


def gate(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok


if __name__ == "__main__":
    import random
    rows = []
    for L in (3, 4, 5, 6):
        T = MyToric(L)
        for a in range(1, L):
            for b in range(1, L):
                for i0 in range(L):
                    for j0 in range(L):
                        VS = rect_vs(L, a, b, i0, j0)
                        rows.append((L, f"{a}x{b}@({i0},{j0})",
                                     "RECT_THICK" if (a >= 2 and b >= 2) else "STRIP",
                                     analyse(T, VS, T.induced(VS))))
        # L-shapes, holed, bands: same families as the lane (from its published listing)
        if L >= 4:
            fams = [(3, 3, 1, 1)] + ([(4, 4, 1, 1), (4, 4, 2, 2), (4, 3, 2, 1)] if L >= 5 else []) \
                   + ([(5, 5, 2, 2), (5, 5, 1, 3), (5, 4, 2, 2)] if L >= 6 else [])
            for (a, b, c, d) in fams:
                if a <= L - 1 and b <= L - 1:
                    VS = rect_vs(L, a, b, 0, 0) - {(s, t) for s in range(c) for t in range(d)}
                    rows.append((L, f"L{a}x{b}-{c}x{d}", "LSHAPE", analyse(T, VS, T.induced(VS))))
            a = b = L - 1
            VS0 = rect_vs(L, a, b, 0, 0)
            mid = (a // 2, b // 2)
            VS = VS0 - {mid}
            rows.append((L, f"H{a}x{b}-v{mid}", "HOLED", analyse(T, VS, T.induced(VS))))
            if L == 6:
                for extra in [{(2, 2)}, {(1, 2)}, {(2, 2), (2, 4)}, {(1, 1), (3, 3)}]:
                    VS = VS0 - extra
                    rows.append((L, f"H{a}x{b}-{sorted(extra)}", "HOLED", analyse(T, VS, T.induced(VS))))
        for a in range(1, L):
            VS = {(i, j) for i in range(a) for j in range(L)}
            rows.append((L, f"band{a}xL", "BAND", analyse(T, VS, T.induced(VS))))
            VS = {(j, i) for i in range(a) for j in range(L)}
            rows.append((L, f"bandLx{a}", "BAND", analyse(T, VS, T.induced(VS))))
        # scatter with the SAME seed protocol as the lane (seed 42+L, sample of range(n))
        # but note the lane's indexing differs; identical numbers are NOT expected, only
        # the qualitative gates.  We use our own seed.
        rnd = random.Random(1000 + L)
        target = (L - 1) * (L - 1) * 2 - 2 * (L - 1) + 1  # AQ of the (L-1)x(L-1) block
        # simpler: recompute from the block row
        blk = [r for r in rows if r[0] == L and r[1] == f"{L-1}x{L-1}@(0,0)"][0][3]
        target = blk["AQ"]
        for k in range(3):
            edges = rnd.sample(range(T.n), target)
            R = 0
            for e0 in edges: R |= 1 << e0
            VS = set()
            for c in range(T.n):
                if (R >> c) & 1:
                    q, d = divmod(c, 2)
                    i, j = divmod(q, L)
                    if d == 0:
                        VS |= {(i, j), (i, (j + 1) % L)}
                    else:
                        VS |= {(i, j), ((i + 1) % L, j)}
            rows.append((L, f"myscatter{k}", "SCATTER", analyse(T, VS, R)))

    allr = [r[3] for r in rows]
    rect = [r[3] for r in rows if r[2] in ("RECT_THICK", "STRIP")]
    thick = [r[3] for r in rows if r[2] == "RECT_THICK"]
    strip = [r[3] for r in rows if r[2] == "STRIP"]
    irr = [r[3] for r in rows if r[2] in ("LSHAPE", "HOLED")]
    band = [r[3] for r in rows if r[2] == "BAND"]
    scat = [r[3] for r in rows if r[2] == "SCATTER"]

    def get(L, tag):
        return next(r[3] for r in rows if r[0] == L and r[1] == tag)

    ok = True
    ok &= gate("census count 1480", len(rect) == 1480, f"{len(rect)}")
    ok &= gate("cap_record == 0 on every contractible region (rect+L+holed)",
               all(r["cap_record"] == 0 for r in rect + irr))
    ok &= gate("cap_record in {1,2} on bands, >0 all",
               all(r["cap_record"] in (1, 2) for r in band))
    ok &= gate("LAW-1 IR2 + 2 r_in + cap_record == 2 AQ, every region",
               all(r["IR2"] + 2 * r["r_in"] + r["cap_record"] == 2 * r["AQ"] for r in allr))
    ok &= gate("LAW-2 IR2 == 2 PER - 10 on THICK rectangles",
               all(r["IR2"] == 2 * r["PER"] - 10 for r in thick))
    n_strip_viol = sum(1 for r in strip if r["IR2"] != 2 * r["PER"] - 10)
    ok &= gate("LAW-2 fails on strips (head-count > 0), so the finding headline "
               "'IR2 = 2*PER - 10 exactly on all 1480 rectangles' is FALSE as stated "
               "(the census includes 606 strips; the law's true scope is the 874 thick "
               "rectangles)", n_strip_viol > 0,
               f"{n_strip_viol}/{len(strip)} strips violate (the rest coincide, e.g. 1x2)")
    ok &= gate("LAW-3 IR2 == 2(AV+PD-AQ) - 2(comp+comp_D) + cap_record, every region",
               all(r["IR2"] == 2 * (r["AV"] + r["PD"] - r["AQ"]) - 2 * (r["comp"] + r["comp_D"])
                   + r["cap_record"] for r in allr))
    ok &= gate("r_in == c1 + c1_D - cap_record, every region",
               all(r["r_in"] == r["c1"] + r["c1_D"] - r["cap_record"] for r in allr))
    ok &= gate("bridge AV+PD-AQ == PER-3 on thick rectangles",
               all(r["AV"] + r["PD"] - r["AQ"] == r["PER"] - 3 for r in thick))
    ok &= gate("bridge breaks somewhere on irregulars",
               any(r["AV"] + r["PD"] - r["AQ"] != r["PER"] - 3 for r in irr))
    ok &= gate("strip law IR2 == 2 AQ",
               all(r["IR2"] == 2 * r["AQ"] for r in strip))
    ok &= gate("cap_pauli == 2AV + PER/2 - 5 on thick rectangles",
               all(r["cap_pauli"] == 2 * r["AV"] + r["PER"] // 2 - 5 for r in thick))
    n_strip_cp = sum(1 for r in strip if r["AQ"] > 0 and r["cap_pauli"] != 2 * r["AV"] + r["PER"] // 2 - 5)
    ok &= gate("cap_pauli formula FAILS on strips too (so 'on rectangles' also over-scoped)",
               n_strip_cp > 0, f"{n_strip_cp} strip violations")
    ok &= gate("scatter qualitative (own seed, own indexing): scatter r_in < matched-block "
               "r_in and scatter IR2 > matched-block IR2, at L=5,6",
               all(rs["r_in"] < get(Lv, f"{Lv-1}x{Lv-1}@(0,0)")["r_in"]
                   and rs["IR2"] > get(Lv, f"{Lv-1}x{Lv-1}@(0,0)")["IR2"]
                   for Lv in (5, 6)
                   for rs in [r[3] for r in rows if r[0] == Lv and r[2] == "SCATTER"]))
    # spot-compare specific rows against the lane's sealed table values
    checks = [
        (get(6, "5x5@(0,0)"), dict(AQ=40, PER=20, IR2=30, cap_pauli=55, cap_record=0)),
        (get(6, "4x4@(0,0)"), dict(AQ=24, PER=16, IR2=22, cap_pauli=35)),
        (get(6, "3x5@(0,0)"), dict(AQ=22, PER=16, IR2=22, cap_pauli=33)),
        (get(6, "L5x5-2x2"), dict(AQ=32, PER=20, IR2=28, cap_pauli=46)),
        (get(6, "H5x5-v(2, 2)"), dict(AQ=36, PER=24, IR2=36, cap_pauli=54)),
        (get(6, "1x5@(0,0)"), dict(AQ=4, PER=12, IR2=8, cap_pauli=8)),
        (get(6, "band5xL"), dict(AQ=54, PER=12, IR2=22, cap_pauli=66, cap_record=2)),
        (get(4, "3x3@(0,0)"), dict(PER=12, IR2=14)),
        (get(4, "L3x3-1x1"), dict(PER=12, IR2=12)),
    ]
    ok &= gate("spot rows match the sealed OUT table (independent recomputation)",
               all(all(row[k] == v for k, v in want.items()) for row, want in checks))
    # the finding JSON's table excerpt numbers that CONTRADICT the sealed OUT:
    r55 = get(6, "5x5@(0,0)"); rb5 = get(6, "band5xL")
    r44 = get(6, "4x4@(0,0)"); r35 = get(6, "3x5@(0,0)")
    ok &= gate("finding-excerpt errors confirmed: 5x5 capP is 55 not 65; band5xL capP is 66 "
               "not 64; 4x4/3x5 capP are 35/33 not 39/38",
               r55["cap_pauli"] == 55 and rb5["cap_pauli"] == 66
               and r44["cap_pauli"] == 35 and r35["cap_pauli"] == 33)
    print()
    print(f"regions: {len(allr)}  (thick {len(thick)}, strips {len(strip)})")
    print(f"ALL INDEPENDENT GATES PASS: {ok}")
