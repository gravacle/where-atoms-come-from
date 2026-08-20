"""ADVERSARIAL VERIFY for T43-A (corner tier).  Fresh implementation, not copied.

Probes:
  P1  independent recompute of STORED / CERT / HIDDEN on L=12 squares s=2..11,
      strips, punctures, comb, g0 -- fresh F_2 code, projection-rank route.
  P2  EXPLICIT CERTIFICATES: for every certifiable class on sample regions, construct
      h supported strictly outside R with g+h in S; verify supp(h) & R == 0 and
      membership by elimination.  (Attacks: is "measuring outside yields g" real?)
  P3  DEFINITION VARIANTS (attack: is the certifiability definition rigged?):
      CERT_B    = certifiers allowed on outside + boundary layer B
      CERT_CENT = certify against the full centralizer (S + logical wraps; i.e. a
                  specific ground state with definite wrap values)
      WRITE_ANY = arbitrary (non-admissible, syndrome-creating) writers through the
                  layer: upper bound 2*BQ, computed
      Each swept over s = 2..11 at L = 12; degrees by finite differences.
  P4  volume-control sensitivity re-check: comb CERT == STORED == 2|R|; punctures
      +6 each; g0 CERT == 0 with STORED volume-degree; scatter (same seed recipe).
  P5  bands: thickness profile WRITE_GI_B via anticommutation pairing, fresh code.
"""

def rref_rank(rows):
    """rank over F_2 by forward elimination with lowest-set-bit pivots (different
       pivot order from the lane's highest-bit routine)."""
    pivs = []
    r = 0
    for v in rows:
        for p in pivs:
            low = p & (-p)
            if v & low:
                v ^= p
        if v:
            pivs.append(v)
            r += 1
    return r


def rref_basis(rows):
    pivs = []
    for v in rows:
        for p in pivs:
            low = p & (-p)
            if v & low:
                v ^= p
        if v:
            pivs.append(v)
    return pivs


def solve_membership(target, gens):
    """return list of gen-indices summing to target, or None.  Tagged elimination."""
    pivs = []  # (vec, tagset)
    for i, g in enumerate(gens):
        v, t = g, {i}
        for (pv, pt) in pivs:
            low = pv & (-pv)
            if v & low:
                v ^= pv
                t = t ^ pt
        if v:
            pivs.append((v, t))
    v, t = target, set()
    for (pv, pt) in pivs:
        low = pv & (-pv)
        if v & low:
            v ^= pv
            t = t ^ pt
    return sorted(t) if v == 0 else None


class Lat:
    def __init__(self, L):
        self.L = L
        self.n = 2 * L * L
        self.stars, self.plaqs = [], []
        for i in range(L):
            for j in range(L):
                s = (1 << self.h(i, j)) | (1 << self.h(i, (j - 1) % L)) \
                    | (1 << self.v(i, j)) | (1 << self.v((i - 1) % L, j))
                self.stars.append(s)
                p = (1 << self.h(i, j)) | (1 << self.h((i + 1) % L, j)) \
                    | (1 << self.v(i, j)) | (1 << self.v(i, (j + 1) % L))
                self.plaqs.append(p)

    def h(self, i, j):
        return (i % self.L) * self.L + (j % self.L)

    def v(self, i, j):
        return self.L * self.L + (i % self.L) * self.L + (j % self.L)

    def region(self, VS):
        R = 0
        for (i, j) in VS:
            if ((i % self.L), ((j + 1) % self.L)) in VS:
                R |= 1 << self.h(i, j)
            if (((i + 1) % self.L), (j % self.L)) in VS:
                R |= 1 << self.v(i, j)
        return R

    def logicals(self):
        """X wraps on horizontal/vertical cycles and Z wraps (dual)."""
        L = self.L
        lx1 = 0
        lx2 = 0
        lz1 = 0
        lz2 = 0
        for k in range(L):
            lx1 |= 1 << self.v(k, 0)      # X along a vertical non-contractible cycle
            lx2 |= 1 << self.h(0, k)      # X along a horizontal cycle
            lz1 |= 1 << self.h(k, 0)      # Z on the dual cycle
            lz2 |= 1 << self.v(0, k)
        return [lx1, lx2], [lz1, lz2]


def sector_counts(gens, R, n):
    """per-sector: (interior_rank, proj_R_rank, CERT_sector) fresh route:
       CERT = dim proj_R(S) - dim interior = certifiable classes."""
    a = rref_rank(list(gens))
    aR = rref_rank([g & R for g in gens])
    Rc = ((1 << n) - 1) & ~R
    aRc = rref_rank([g & Rc for g in gens])
    interior = a - aRc
    return interior, aR, aR - interior


def analyse(T, VS, gens_x=None, gens_z=None):
    GX = T.stars if gens_x is None else gens_x
    GZ = T.plaqs if gens_z is None else gens_z
    n = T.n
    R = T.region(VS)
    AQ = bin(R).count("1")
    ix, pxR, cx = sector_counts(GX, R, n)
    iz, pzR, cz = sector_counts(GZ, R, n)
    STORED = 2 * AQ - (ix + iz)
    CERT = cx + cz
    # boundary layer
    Rc = ((1 << n) - 1) & ~R
    B = 0
    for g in list(GX) + list(GZ):
        if (g & R) and (g & Rc):
            B |= g & R
    BQ = bin(B).count("1")
    return dict(R=R, AQ=AQ, STORED=STORED, CERT=CERT, HID=STORED - CERT, B=B, BQ=BQ)


def cert_with_span(T, R, span_gens_x, span_gens_z, allowed_outside_mask):
    """dim(P_R ^ (span + P_O)) - dim(span ^ P_R), per CSS sector, direct formula.
       O = allowed_outside_mask.  A+B fills sector => dim = AQ + rank(span|e_O) - n."""
    n = T.n
    AQ = bin(R).count("1")
    total = 0
    for gens in (span_gens_x, span_gens_z):
        rows = list(gens) + [1 << e for e in range(n) if (allowed_outside_mask >> e) & 1]
        d_span = rref_rank(rows)
        # check A+B = full sector
        d_all = rref_rank(rows + [1 << e for e in range(n) if (R >> e) & 1])
        assert d_all == n, "span+P_R does not fill the sector"
        inter = AQ + d_span - n
        a = rref_rank(list(gens))
        Rc = ((1 << n) - 1) & ~R
        aRc = rref_rank([g & Rc for g in gens])
        interior = a - aRc
        total += inter - interior
    return total


def degree(seq, maxk=4):
    cur = list(seq)
    if len(set(cur)) == 1:
        return 0
    for k in range(1, maxk + 1):
        cur = [b - a for a, b in zip(cur, cur[1:])]
        if len(cur) >= 2 and len(set(cur)) == 1:
            return k if cur[0] != 0 else k - 1
    return None


def rect(L, a, b, i0=0, j0=0):
    return {(((i0 + s) % L), ((j0 + t) % L)) for s in range(a) for t in range(b)}


def main():
    ok = True

    def gate(label, cond, detail=""):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
        ok = ok and cond

    print("VERIFY-P1: independent recompute, L=12 squares s=2..11")
    T = Lat(12)
    S_, C_, H_, BQ_ = [], [], [], []
    for s in range(2, 12):
        r = analyse(T, rect(12, s, s))
        S_.append(r["STORED"]); C_.append(r["CERT"]); H_.append(r["HID"]); BQ_.append(r["BQ"])
    print(f"  STORED={S_}")
    print(f"  CERT  ={C_}")
    print(f"  HID   ={H_}")
    print(f"  BQ    ={BQ_}")
    gate("STORED matches lane (2s^2+2s-5)", S_ == [2*s*s + 2*s - 5 for s in range(2, 12)])
    gate("CERT matches lane (8s-10)", C_ == [8*s - 10 for s in range(2, 12)])
    gate("HIDDEN matches lane (2s^2-6s+5)", H_ == [2*s*s - 6*s + 5 for s in range(2, 12)])
    gate("degrees: STORED=2 CERT=1 HIDDEN=2",
         degree(S_) == 2 and degree(C_) == 1 and degree(H_) == 2)
    strip_eq = []
    for b in range(2, 12):
        r = analyse(T, rect(12, 1, b))
        strip_eq.append(r["CERT"] == r["STORED"])
    gate("strips: CERT == STORED (nothing hidden)", all(strip_eq))

    print()
    print("VERIFY-P2: explicit certificates -- construct h strictly outside with g+h in S")
    # sample regions: 4x4 and 7x7 at L=12
    for s in (4, 7):
        VS = rect(12, s, s)
        R = T.region(VS)
        Rc = ((1 << T.n) - 1) & ~R
        n_cert, n_checked = 0, 0
        for gens in (T.stars, T.plaqs):
            projR = rref_basis([g & R for g in gens])
            intb = interior_basis(gens, R, T.n)
            # class reps: projR elements independent modulo interior span
            reps = []
            acc = list(intb)
            rk = rref_rank(acc)
            for v in projR:
                if rref_rank(acc + [v]) > rk:
                    acc.append(v)
                    rk += 1
                    reps.append(v)
            for g in reps:
                # find stabiliser element w with w&R == g  (solve in restriction)
                sol = solve_membership(g, [x & R for x in gens])
                assert sol is not None
                w = 0
                for i in sol:
                    w ^= gens[i]
                h = w & Rc
                assert (w & R) == g and (h & R) == 0
                # membership: g+h == w is a sum of generators by construction
                n_cert += 1
            n_checked += len(reps)
        print(f"  {s}x{s}: {n_cert} certifiable classes, every one given an explicit "
              f"strictly-outside certificate h with g+h in S")
        gate(f"{s}x{s}: certificate count == CERT", n_cert == 8*s - 10, f"count={n_cert}")

    print()
    print("VERIFY-P3: definition variants (the rigging attack), L=12 squares s=2..11")
    lx, lz = T.logicals()
    CB_, CC_, W2B_ = [], [], []
    for s in range(2, 12):
        VS = rect(12, s, s)
        r = analyse(T, VS)
        R = r["R"]
        O_layer = (((1 << T.n) - 1) & ~R) | r["B"]
        cb = cert_with_span(T, R, T.stars, T.plaqs, O_layer)
        cc = cert_with_span(T, R, list(T.stars) + lx, list(T.plaqs) + lz,
                            ((1 << T.n) - 1) & ~R)
        CB_.append(cb); CC_.append(cc); W2B_.append(2 * r["BQ"])
    d1 = lambda q: [b - a for a, b in zip(q, q[1:])]
    print(f"  CERT_B    (certifiers on outside+layer) = {CB_}  D1={d1(CB_)}")
    print(f"  CERT_CENT (fixed ground state, S+wraps) = {CC_}  D1={d1(CC_)}")
    print(f"  WRITE_ANY (arbitrary writers, bound 2BQ)= {W2B_}  degree={degree(W2B_)}")
    gate("CERT_B: all-layer crossover at s=2,3 (CERT_B == STORED there), then "
         "BOUNDARY-degree (deg 1, D1 == 16 const) for s >= 4 -- layer access does "
         "not flip the verdict",
         CB_[0] == S_[0] and CB_[1] == S_[1] and degree(CB_[2:]) == 1)
    gate("CERT_CENT stays BOUNDARY-degree for s >= 4 (D1 == 8 const; = CERT + 4 wraps): "
         "fixing wrap values does not flip it",
         degree(CC_[2:]) == 1 and all(cc == c + 4 for cc, c in zip(CC_[2:], C_[2:])))
    gate("arbitrary writers through the layer are bounded by 2BQ, BOUNDARY-degree",
         degree(W2B_) == 1)
    gate("variants never exceed STORED", all(cb <= st and cc <= st
         for cb, cc, st in zip(CB_, CC_, S_)))

    print()
    print("VERIFY-P4: sensitivity controls, fresh code")
    base = rect(12, 9, 9)
    holes = [(2, 2), (2, 6), (6, 2), (6, 6)]
    pc_ = []
    for k in range(5):
        pc_.append(analyse(T, base - set(holes[:k]))["CERT"])
    print(f"  punctures CERT = {pc_}")
    gate("punctures: CERT = 62,68,74,80,86 (+6 each)", pc_ == [62, 68, 74, 80, 86])
    combVS = {(i, j) for i in range(9) for j in (0, 2, 4, 6, 8)}
    c = analyse(T, combVS)
    gate("comb: CERT == STORED == 2|R|", c["CERT"] == c["STORED"] == 2 * c["AQ"],
         f"CERT={c['CERT']} AQ={c['AQ']}")
    T8 = Lat(8)
    g0S, g0C = [], []
    for s in range(2, 8):
        VS = rect(8, s, s)
        R = T8.region(VS)
        gz = [1 << e for e in range(T8.n)]
        ix, pxR, cx = sector_counts([], R, T8.n) if False else (0, 0, 0)
        # X sector empty gens; Z sector single-site
        iz, pzR, cz = sector_counts(gz, R, T8.n)
        AQ = bin(R).count("1")
        g0S.append(2 * AQ - iz)
        g0C.append(cz)
    print(f"  g0: STORED={g0S} (deg {degree(g0S)}), CERT={g0C}")
    gate("g0: STORED volume-degree, CERT identically 0",
         degree(g0S) == 2 and set(g0C) == {0})
    import random
    rnd = random.Random(43)
    b77 = analyse(T8, rect(8, 7, 7))
    scat_c = []
    for k in range(3):
        edges = rnd.sample(range(T8.n), b77["AQ"])
        R = 0
        for e in edges:
            R |= 1 << e
        ix, _, cx = sector_counts(T8.stars, R, T8.n)
        iz, _, cz = sector_counts(T8.plaqs, R, T8.n)
        scat_c.append(cx + cz)
    print(f"  scatter CERT = {scat_c} vs block {b77['CERT']}")
    gate("scatter beats block in every seed (86/84/86 vs 46)",
         scat_c == [86, 84, 86] and b77["CERT"] == 46)

    print()
    print("VERIFY-P5: band thickness profile, fresh pairing code, L=12 widths 1..5")
    prof = []
    for a in range(1, 6):
        VS = {(i, j) for i in range(a) for j in range(12)}
        R = T.region(VS)
        Rc = ((1 << T.n) - 1) & ~R
        B = 0
        for g in list(T.stars) + list(T.plaqs):
            if (g & R) and (g & Rc):
                B |= g & R
        O = Rc | B
        colsO = [e for e in range(T.n) if (O >> e) & 1]
        # admissible writers on O: X-writers commute with all plaqs, Z-writers with stars
        WX = kernel(T.plaqs, colsO, T.n)
        WZ = kernel(T.stars, colsO, T.n)
        # record reps: interior ops commuting with all opposite gens, modulo stabiliser span
        wgi = 0
        for (gens_same, gens_opp, WR) in ((T.stars, T.plaqs, WZ), (T.plaqs, T.stars, WX)):
            colsR = [e for e in range(T.n) if (R >> e) & 1]
            ker = kernel(gens_opp, colsR, T.n)
            gb = rref_basis(list(gens_same))
            acc = list(gb)
            rk = rref_rank(acc)
            reps = []
            for v in ker:
                if rref_rank(acc + [v]) > rk:
                    acc.append(v)
                    rk += 1
                    reps.append(v)
            rows = []
            for r0 in reps:
                m = 0
                for j, w in enumerate(WR):
                    if bin(r0 & w).count("1") & 1:
                        m |= 1 << j
                rows.append(m)
            wgi += rref_rank(rows)
        prof.append(wgi)
    print(f"  WRITE_GI_B widths 1..5 = {prof}")
    gate("thickness profile 1,2,1,0,0 (protection at width>=4)", prof == [1, 2, 1, 0, 0])

    print()
    print(f"VERIFY ALL: {ok}")
    return ok


def interior_basis(gens, R, n):
    Rc = ((1 << n) - 1) & ~R
    pivs = []  # (outside_part, full_sum)
    out = []
    for g in gens:
        v, f = g & Rc, g
        for (pv, pf) in pivs:
            low = pv & (-pv)
            if v & low:
                v ^= pv
                f ^= pf
        if v:
            pivs.append((v, f))
        else:
            if f:
                assert f & Rc == 0
                out.append(f)
    return rref_basis(out)


def _full_interior(gens, R, n):
    return interior_basis(gens, R, n)


def kernel(rows, cols, n):
    colv = []
    for idx, c in enumerate(cols):
        v = 0
        for i, rw in enumerate(rows):
            if (rw >> c) & 1:
                v |= 1 << i
        colv.append((v, 1 << idx))
    pivs = []
    ker = []
    for v, t in colv:
        for (pv, pt) in pivs:
            low = pv & (-pv)
            if v & low:
                v ^= pv
                t ^= pt
        if v:
            pivs.append((v, t))
        else:
            ker.append(t)
    out = []
    for t in ker:
        m = 0
        for idx, c in enumerate(cols):
            if (t >> idx) & 1:
                m |= 1 << c
        out.append(m)
    return out


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
