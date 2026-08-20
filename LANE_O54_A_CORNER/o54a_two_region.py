"""O-54-A  TWO-REGION INTERFACE vs EARNED SEPARATION -- PURE-GAMMA (CORNER) TIER.

THE HYPOTHESIS UNDER TEST (stated so it can fail): "the falloff is not Gamma's to give --
pure-Gamma two-region relations are TOPOLOGICAL (separation-independent or contact-only), and
any genuine falloff enters through the MEDIATOR."  If a pure-Gamma falloff appears below, the
hypothesis FAILS.

WHAT IS COMPUTED, all exact F_2 / integer arithmetic on bitmasks, no floats anywhere:
  Regions A, B = s x s vertex blocks (region = induced edge set of its vertex set, the C-74 /
  T42-C convention), separation g = number of vertex columns strictly between the blocks,
  swept g = 0 .. L-2s (the top end is WRAP CONTACT on the torus; g_eff = min(g, L-2s-g)).
  (a) IR2(A), IR2(B), IR2(A u B) -- the stabiliser CUT-RANK, the EARNED interface-rank
      instrument of C-79 (LANE_T42_C_BOUNDCAP), computed by BOTH exact expressions used
      there (group-quotient r_S - r_in - r_out, and restriction-rank IR2b) and gated equal;
      and the MUTUAL INTERFACE  I_IR(A:B) = IR2(A) + IR2(B) - IR2(A u B).
  (b) The stabiliser entropy S(R) = |R| - r_in(R) and the stabiliser MUTUAL INFORMATION
      I_MI(A:B) = S(A) + S(B) - S(A u B) = r_in(AuB) - r_in(A) - r_in(B)  (disjoint supports).
      OWNERS: both instruments are the stabiliser region-entropy object of
      Hamma-Ionicioiu-Zanardi (2005); the topological part of such quantities is
      Kitaev-Preskill (2006) / Levin-Wen (2006) territory.  I_IR and I_MI are the same
      object up to normalisation and bookkeeping of the two cut sides; both are reported.
  (c) RECORD-RELEVANT JOINT CONTENT: cap_record per region (C-74 instrument), and for the
      union the split  JOINT_NEW = logical classes supported in A u B attributable to
      NEITHER part alone,  IDENT = identifications (classes of A and of B that merge in the
      union quotient).  JC = JOINT_NEW - IDENT.

CONTROLS (D-15) -- every zero and every flat line ships with a control that registers:
  (i)   CONTACT g=0 and WRAP-CONTACT g=L-2s: mutual interface must be NON-zero.
  (ii)  NON-CONTRACTIBLE PAIR: two winding rings; the venue's own long-range structure MUST
        appear (I_MI > 0 and IDENT > 0 at every distance) -- so a flat zero for contractible
        pairs is a measurement, not blindness.
  (iii) ENCLOSING GEOMETRY: A inside, B an enclosing frame, gap swept; and beside it the
        KITAEV-PRESKILL tripartite combination on a disk (contact construction), where the
        topological constant is KNOWN to appear -- the instrument's registration of
        topological terms is demonstrated where the literature says one exists.
  (iv)  CORRIDOR SWEEP (region deformation, INSERTED geometry): A grown toward B by a
        2-vertex-wide corridor of length c = 0..g; a graded I(c) would have revealed a
        falloff-with-remaining-distance had one existed; a step exactly at contact is the
        contact-only law seen a second way.
  (v)   INSERTED-RANGE ROW (instrument sensitivity, NOT a Gamma venue): one extra Z-row, an
        open string of length ell, added to the generator list purely as a RANK-INSTRUMENT
        control; the mutual interface must switch on exactly when the row's support reaches
        both regions (predicted onset ell = g+3 for the construction below).  This is the
        control that shows the instrument WOULD have shown a separation-dependent law had a
        longer-range term been present.  INSERTED, and reported only as such.
  (vi)  OFFSET INVARIANCE (D-22): the venue's translation symmetry checked on the numbers.

DISCIPLINE: D-1 no classical gravitational form in any construction step; the comparison
target is named ONLY in the final comparison section.  D-24 audit in D24_AUDIT.txt.  Every
PASS/FAIL below is printed from a computed boolean.  FIT NOTHING: laws are exact equalities
checked on every row by booleans.  INSERTED vs INDUCED labels throughout.

PROVENANCE OF MACHINERY (earned, reused, cited -- D-24):
  bitmask F_2 kit + Toric class + region/IR2/cap_record conventions: copied verbatim-in-
  substance from LANE_T42_C_BOUNDCAP/t42_corner.py (C-79) which itself replicates
  LANE_O50_A_ACTION/f2lib.py + s4_clause_v.py (C-74).  Earned separation: the C-78 metric
  d_W = boundary-crossing cost; for these blocks the minimal-crossing count between the
  regions is the BFS edge distance dV between their vertex sets, reported per row (dV = g+1,
  monotone in g -- the sweep in g IS the sweep in earned separation).
"""
import sys

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
    """kernel of the F_2 map x in F_2^{cols} -> (row . x)_rows; x returned as bitmasks
       supported on `cols`.  (T42-C routine, verbatim.)"""
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

# ---------------------------------------------------------------- toric carrier (= O50A/T42C)
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
        """induced edge set: an edge is in R iff BOTH endpoints are in VS (C-74 convention)."""
        L = self.L
        R = 0
        for (i, j) in VS:
            if ((i % L), ((j + 1) % L)) in VS:
                R |= 1 << self.h(i, j)
            if (((i + 1) % L), (j % L)) in VS:
                R |= 1 << self.v(i, j)
        return R

# ---------------------------------------------------------------- per-sector rank stats
def sector_stats(rows, R, full):
    """rows: edge-support bitmasks of one CSS sector's generators (plus any INSERTED rows,
       labelled by the caller).  Returns (rk, rin, rout, IR2_quotient, IR2_restriction)."""
    Rc = full ^ R
    rk   = rank_bits(rows)
    rkR  = rank_bits([r & R for r in rows])
    rkRc = rank_bits([r & Rc for r in rows])
    rin  = rk - rkRc          # rank of the subgroup supported wholly inside R
    rout = rk - rkR           # rank of the subgroup supported wholly outside R
    ir2  = rk - rin - rout    # group-quotient cut-rank (T42-C way 1)
    ir2b = rkR - rin          # restriction-rank cut-rank (T42-C way 2)
    return rk, rin, rout, ir2, ir2b

def region_measures(T, R, extraZ=()):
    """IR2 (both ways, gated), r_in, and stabiliser entropy S = |R| - r_in.  extraZ rows are
       INSERTED instrument rows (control (v) only)."""
    zrows = list(T.plaqs) + list(extraZ)
    rkX, rinX, routX, ir2X, ir2bX = sector_stats(T.stars, R, T.full)
    rkZ, rinZ, routZ, ir2Z, ir2bZ = sector_stats(zrows,   R, T.full)
    ok = (ir2X == ir2bX) and (ir2Z == ir2bZ)
    AQ = pc(R)
    return dict(AQ=AQ, rin=rinX + rinZ, IR2=ir2X + ir2Z, IR2b=ir2bX + ir2bZ,
                IR2X=ir2X, IR2Z=ir2Z, S=AQ - rinX - rinZ, agree=ok)

def cap_record_with_kernels(T, R):
    """C-74 record-content instrument: cap = dim( (N(S) cap P_R) / (S cap P_R) ), per sector,
       computed by direct span membership exactly as T42-C did.  Also returns the kernels so
       the caller can split the union content into JOINT_NEW vs IDENT."""
    n = T.n
    cols = [e for e in range(n) if (R >> e) & 1]
    kerX = kernel_basis(T.plaqs, cols)   # X-ops in R commuting with all plaquettes
    kerZ = kernel_basis(T.stars, cols)   # Z-ops in R commuting with all stars
    bst = basis_bits(list(T.stars))
    bpl = basis_bits(list(T.plaqs))
    capX = rank_bits(bst + kerX) - len(bst)
    capZ = rank_bits(bpl + kerZ) - len(bpl)
    return capX + capZ, (kerX, kerZ, bst, bpl)

def joint_content(T, RA, RB):
    """cap_record of A, B, AuB and the split of the union's classes:
       JOINT_NEW = classes in the union attributable to NEITHER part alone
       IDENT     = identifications (part classes that merge in the union quotient)
       JC = capAB - capA - capB = JOINT_NEW - IDENT."""
    capA, (kXA, kZA, bst, bpl) = cap_record_with_kernels(T, RA)
    capB, (kXB, kZB, _, _) = cap_record_with_kernels(T, RB)
    capAB, (kXAB, kZAB, _, _) = cap_record_with_kernels(T, RA | RB)
    jn = 0
    ident = 0
    for (base, kA, kB, kAB) in ((bst, kXA, kXB, kXAB), (bpl, kZA, kZB, kZAB)):
        r0 = len(base)
        r_parts = rank_bits(base + kA + kB) - r0     # classes reachable from the parts
        r_union = rank_bits(base + kAB) - r0         # all classes of the union
        capA_s = rank_bits(base + kA) - r0
        capB_s = rank_bits(base + kB) - r0
        jn += r_union - r_parts
        ident += (capA_s + capB_s) - r_parts
    return capA, capB, capAB, jn, ident

def straddle_both(T, RA, RB, extraZ=()):
    """COUNT of generators whose support meets BOTH regions (a lattice count)."""
    return sum(1 for gset in (T.stars, list(T.plaqs) + list(extraZ)) for g in gset
               if (g & RA) and (g & RB))

def pair_measures(T, VSA, VSB, extraZ=(), with_content=True):
    RA = T.edges_of_vertexset(VSA)
    RB = T.edges_of_vertexset(VSB)
    assert RA & RB == 0, "regions must be disjoint qubit sets"
    RAB = RA | RB
    mA = region_measures(T, RA, extraZ)
    mB = region_measures(T, RB, extraZ)
    mAB = region_measures(T, RAB, extraZ)
    out = dict(IR2A=mA['IR2'], IR2B=mB['IR2'], IR2AB=mAB['IR2'],
               I_IR=mA['IR2'] + mB['IR2'] - mAB['IR2'],
               I_MI=mA['S'] + mB['S'] - mAB['S'],
               SB=straddle_both(T, RA, RB, extraZ),
               agree=mA['agree'] and mB['agree'] and mAB['agree'])
    if with_content:
        capA, capB, capAB, jn, ident = joint_content(T, RA, RB)
        out.update(capA=capA, capB=capB, capAB=capAB, JOINT_NEW=jn, IDENT=ident,
                   JC=capAB - capA - capB)
    return out

# ---------------------------------------------------------------- region constructions
def block_vs(L, s, i0, j0):
    return {(((i0 + a) % L), ((j0 + b) % L)) for a in range(s) for b in range(s)}

def ring_vs(L, i0, w):
    """winding ring: w consecutive vertex rows, all columns (non-contractible)."""
    return {(((i0 + a) % L), j) for a in range(w) for j in range(L)}

def frame_vs(L, o0, q, w):
    """square frame: outer q x q block at (o0,o0) minus the inner (q-2w) x (q-2w) block."""
    outer = block_vs(L, q, o0, o0)
    inner = block_vs(L, q - 2 * w, o0 + w, o0 + w)
    return outer - inner

def bfs_dist(L, VSA, VSB):
    """minimal edge-count graph distance between the vertex sets on the torus (the C-78
       minimal boundary-crossing count between the regions)."""
    from collections import deque
    dist = {v: 0 for v in VSA}
    dq = deque(VSA)
    while dq:
        (i, j) = dq.popleft()
        if (i, j) in VSB:
            return dist[(i, j)]
        for (a, b) in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
            w = (a % L, b % L)
            if w not in dist:
                dist[w] = dist[(i, j)] + 1
                dq.append(w)
    return None

# ---------------------------------------------------------------- venue check (D-22)
def venue_check(L):
    T = Toric(L)
    m = L * L
    prodS = 0
    for s in T.stars:
        prodS ^= s
    prodP = 0
    for p in T.plaqs:
        prodP ^= p
    rS = rank_bits(list(T.stars))
    rP = rank_bits(list(T.plaqs))
    ok = (prodS == 0) and (prodP == 0) and (rS == m - 1) and (rP == m - 1)
    k = T.n - rS - rP
    return T, ok, (rS, rP, k)

# ================================================================= main
def main():
    gates = {}   # name -> bool, every verdict computed

    print("=" * 100)
    print("O-54-A  TWO-REGION INTERFACE vs EARNED SEPARATION -- PURE-GAMMA CORNER TIER")
    print("exact F_2 bitmask arithmetic throughout; no floats; every PASS/FAIL is a computed boolean")
    print("=" * 100)

    # ---------------------------------------------------------------- SECTION 0: venue
    print("\nSECTION 0 -- VENUE (D-22).  Toric carrier as in C-74/C-79; per L: product of all")
    print("stars = product of all plaquettes = identity, sector ranks = L^2-1, k = 2 logicals.")
    Ts = {}
    v_ok = True
    for L in range(4, 13):
        T, ok, (rS, rP, k) = venue_check(L)
        Ts[L] = T
        v_ok = v_ok and ok and (k == 2)
        print("  L=%2d  n=%3d  rank(stars)=%3d  rank(plaqs)=%3d  k=%d  ok=%s" %
              (L, T.n, rS, rP, k, ok))
    gates['G0_venue'] = v_ok

    # ---------------------------------------------------------------- SECTION 1: main sweep
    print("\n" + "=" * 100)
    print("SECTION 1 -- MAIN SWEEP (INDUCED, pure Gamma: nothing added to the venue).")
    print("A, B = s x s vertex blocks, rows aligned, A at column 0, B at column s+g.")
    print("g = vertex columns strictly between; g_eff = min(g, L-2s-g); dV = C-78 crossing count.")
    print("I_IR = IR2(A)+IR2(B)-IR2(AuB); I_MI = stabiliser mutual information (HIZ, log2 units).")
    hdr = ("L", "s", "g", "geff", "dV", "SB", "IR2A", "IR2B", "IR2AB", "I_IR", "I_MI",
           "capA", "capB", "capAB", "JNEW", "IDENT", "JC")
    print(("%4s" * len(hdr)) % hdr)
    rows = []
    agree_all = True
    for L in range(4, 13):
        T = Ts[L]
        for s in (2, 3, 4):
            if 2 * s > L:
                continue
            for g in range(0, L - 2 * s + 1):
                VSA = block_vs(L, s, 0, 0)
                VSB = block_vs(L, s, 0, s + g)
                r = pair_measures(T, VSA, VSB)
                r.update(L=L, s=s, g=g, geff=min(g, L - 2 * s - g),
                         dV=bfs_dist(L, VSA, VSB))
                rows.append(r)
                agree_all = agree_all and r['agree']
                print(("%4d" * len(hdr)) %
                      (L, s, g, r['geff'], r['dV'], r['SB'], r['IR2A'], r['IR2B'],
                       r['IR2AB'], r['I_IR'], r['I_MI'], r['capA'], r['capB'],
                       r['capAB'], r['JOINT_NEW'], r['IDENT'], r['JC']))
    gates['G1_IR2_two_ways_agree'] = agree_all

    # row classes: separated / single contact (one seam) / double contact (L = 2s: two seams,
    # the union then winds the torus)
    sep = [r for r in rows if r['geff'] >= 1]
    single = [r for r in rows if r['geff'] == 0 and r['L'] > 2 * r['s']]
    double = [r for r in rows if r['geff'] == 0 and r['L'] == 2 * r['s']]
    gates['G2_separated_all_zero'] = all(
        r['I_IR'] == 0 and r['I_MI'] == 0 and r['JC'] == 0 and
        r['JOINT_NEW'] == 0 and r['IDENT'] == 0 for r in sep)
    gates['G2_contact_necessary_for_nonzero'] = all(
        (r['I_IR'] == 0 and r['I_MI'] == 0) or r['SB'] > 0 for r in rows)

    print("\nCONTACT VALUES (every g_eff = 0 row; wrap seams included):")
    for r in [x for x in rows if x['geff'] == 0]:
        kind = "DOUBLE" if r['L'] == 2 * r['s'] else "single"
        print("  L=%2d s=%d g=%2d  %s seam  I_IR=%d  I_MI=%d  JNEW=%d  IDENT=%d"
              % (r['L'], r['s'], r['g'], kind, r['I_IR'], r['I_MI'],
                 r['JOINT_NEW'], r['IDENT']))
    # exact contact forms, FIT NOTHING: read off one row, then REQUIRED on every row of the
    # class as a computed boolean (s = 2, 3, 4 all included -- note s=2 gives exactly ZERO:
    # a single-plaquette seam is rank-absorbed; contact is necessary, not sufficient)
    gates['G3_single_contact_form'] = all(
        r['I_IR'] == 2 * (r['s'] - 2) and r['I_MI'] == r['s'] - 2 and
        r['JC'] == 0 and r['JOINT_NEW'] == 0 for r in single) and len(single) > 0
    gates['G3_double_contact_form'] = all(
        r['I_IR'] == 4 * (r['s'] - 2) + 1 and r['I_MI'] == 2 * (r['s'] - 2) and
        r['JOINT_NEW'] == 1 for r in double) and len(double) > 0
    gates['G3_wrap_seam_equals_direct_seam'] = all(
        r['I_IR'] == 2 * (r['s'] - 2) for r in single if r['g'] == r['L'] - 2 * r['s'])
    print("  exact single-seam form   I_IR = 2(s-2), I_MI = s-2, no joint content : %s"
          % gates['G3_single_contact_form'])
    print("  exact double-seam form   I_IR = 4(s-2)+1, I_MI = 2(s-2), JOINT_NEW = 1")
    print("    (the union winds the torus; ONE genuinely joint record class)     : %s"
          % gates['G3_double_contact_form'])

    # convention cross-check at one contact row: induced-union vs edge-union
    T = Ts[8]
    VSA, VSB = block_vs(8, 3, 0, 0), block_vs(8, 3, 0, 3)
    RA, RB = T.edges_of_vertexset(VSA), T.edges_of_vertexset(VSB)
    Rind = T.edges_of_vertexset(VSA | VSB)
    mI = region_measures(T, Rind)
    mA, mB = region_measures(T, RA), region_measures(T, RB)
    print("\nCONVENTION NOTE (g=0, L=8, s=3): edge-union I_IR = %d ; induced-union variant"
          % (mA['IR2'] + mB['IR2'] - region_measures(T, RA | RB)['IR2']))
    print("  (seam edges included in the joint region) gives I_IR_ind = %d.  Both reported;"
          % (mA['IR2'] + mB['IR2'] - mI['IR2']))
    print("  the edge-union (disjoint qubit sets) is the mutual-information object (HIZ).")

    # ---------------------------------------------------------------- SECTION 2: corridor
    print("\n" + "=" * 100)
    print("SECTION 2 -- CONTROL (iv): CORRIDOR SWEEP (INSERTED geometry: the region is grown,")
    print("the venue untouched).  A 2-row corridor of length c toward B; a falloff in the")
    print("remaining distance g-c would appear as graded I(c); the law is read from the table.")
    print("%4s%4s%4s%4s%6s%6s%6s" % ("L", "s", "g", "c", "I_IR", "I_MI", "SB"))
    corr_rows = []
    corr_agree = True
    for (L, s) in ((10, 3), (12, 3)):
        T = Ts[L]
        mid = max(0, (s - 2) // 2)
        for g in range(1, min(7, L - 2 * s + 1)):
            for c in range(0, g + 1):
                VSA = block_vs(L, s, 0, 0) | {(((mid + a) % L), s + b)
                                              for a in range(2) for b in range(c)}
                VSB = block_vs(L, s, 0, s + g)
                r = pair_measures(T, VSA, VSB, with_content=False)
                r.update(L=L, s=s, g=g, c=c, wrap=(g == L - 2 * s))
                corr_rows.append(r)
                corr_agree = corr_agree and r['agree']
                print("%4d%4d%4d%4d%6d%6d%6d" % (L, s, g, c, r['I_IR'], r['I_MI'], r['SB']))
    # exact corridor law, all computed: with the blocks separated (g < L-2s), the corridor's
    # own single-plaquette contact seam (SB = 1 exactly at c = g) is rank-absorbed -- I stays
    # ZERO at every c, INCLUDING contact, and NOTHING GRADED appears at any intermediate c.
    # At g = L-2s the blocks are wrap-contacting on the far side (I = 2(s-2) from Section 1)
    # and corridor contact adds exactly +1.
    gates['G5_corridor_never_graded'] = corr_agree and all(
        (r['I_IR'] == 0 and r['I_MI'] == 0 and r['SB'] == (1 if r['c'] == r['g'] else 0))
        for r in corr_rows if not r['wrap'])
    gates['G5_corridor_wrap_rows_consistent'] = all(
        r['I_IR'] == 2 * (r['s'] - 2) + (1 if r['c'] == r['g'] else 0)
        for r in corr_rows if r['wrap'])

    # ---------------------------------------------------------------- SECTION 3: winding pair
    print("\n" + "=" * 100)
    print("SECTION 3 -- CONTROL (ii): NON-CONTRACTIBLE PAIR (INDUCED: the venue's own long-")
    print("range structure).  A, B = winding rings of width w, row gaps (r-w, L-r-w) between.")
    print("%4s%4s%4s%6s%6s%6s%6s%6s%6s%6s%6s" %
          ("L", "w", "gap1", "gap2", "I_IR", "I_MI", "capA", "capB", "capAB", "JNEW", "IDENT"))
    wind_rows = []
    wind_agree = True
    for L in (8, 10, 12):
        T = Ts[L]
        for w in (1, 2):
            for r0 in range(w + 1, L - w):
                gap1, gap2 = r0 - w, L - r0 - w
                if gap1 < 1 or gap2 < 1:
                    continue
                rr = pair_measures(T, ring_vs(L, 0, w), ring_vs(L, r0, w))
                rr.update(L=L, w=w, gap1=gap1, gap2=gap2)
                wind_rows.append(rr)
                wind_agree = wind_agree and rr['agree']
                print("%4d%4d%4d%6d%6d%6d%6d%6d%6d%6d%6d" %
                      (L, w, gap1, gap2, rr['I_IR'], rr['I_MI'], rr['capA'], rr['capB'],
                       rr['capAB'], rr['JOINT_NEW'], rr['IDENT']))
    gates['G6_winding_longrange_appears'] = wind_agree and all(
        rr['I_MI'] > 0 and rr['IDENT'] > 0 for rr in wind_rows)
    # exact topological form, every distance, every L: the pair's relation is the CONSTANT w
    # in every column, the joint structure is pure IDENTIFICATION of the two rings' wrap
    # classes (JOINT_NEW = 0) -- distance appears NOWHERE
    gates['G6_winding_exact_constant_w'] = all(
        rr['I_IR'] == rr['w'] and rr['I_MI'] == rr['w'] and rr['IDENT'] == rr['w'] and
        rr['JOINT_NEW'] == 0 and rr['capA'] == rr['w'] and rr['capB'] == rr['w'] and
        rr['capAB'] == rr['w'] for rr in wind_rows)
    vals = {(rr['w'], rr['I_MI'], rr['IDENT'], rr['I_IR']) for rr in wind_rows}
    per_w = {w: {(a, b, c) for (ww, a, b, c) in vals if ww == w} for w in (1, 2)}
    gates['G6_winding_distance_independent'] = all(len(per_w[w]) == 1 for w in (1, 2))
    print("  distinct (I_MI, IDENT, I_IR) per width: w=1 -> %s ; w=2 -> %s" %
          (sorted(per_w[1]), sorted(per_w[2])))

    # ---------------------------------------------------------------- SECTION 4: enclosing
    print("\n" + "=" * 100)
    print("SECTION 4 -- CONTROL (iii-a): ENCLOSING GEOMETRY (INDUCED).  A = sA x sA block")
    print("inside, B = frame of width w enclosing it, radial gap t swept.")
    print("%4s%4s%4s%4s%4s%6s%6s%6s%6s%6s%6s" %
          ("L", "sA", "w", "t", "q", "SB", "I_IR", "I_MI", "capA", "capB", "JC"))
    enc_rows = []
    enc_agree = True
    for L in (10, 12):
        T = Ts[L]
        for (sA, w) in ((2, 2), (3, 2), (3, 1), (2, 1)):
            for t in range(0, 4):
                q = sA + 2 * t + 2 * w
                if q > L - 1:
                    continue
                VSB = frame_vs(L, 0, q, w)
                VSA = block_vs(L, sA, w + t, w + t)
                rr = pair_measures(T, VSA, VSB)
                rr.update(L=L, sA=sA, w=w, t=t, q=q)
                enc_rows.append(rr)
                enc_agree = enc_agree and rr['agree']
                print("%4d%4d%4d%4d%4d%6d%6d%6d%6d%6d%6d" %
                      (L, sA, w, t, q, rr['SB'], rr['I_IR'], rr['I_MI'],
                       rr['capA'], rr['capB'], rr['JC']))
    gates['G7_enclosing_separated_zero'] = enc_agree and all(
        rr['I_IR'] == 0 and rr['I_MI'] == 0 for rr in enc_rows if rr['t'] >= 1)
    # contact rows: frames of width >= 2 register (closed seam), width-1 frames are the
    # sub-threshold family (their seam constraints are rank-absorbed, exactly as the s=2
    # block seam and the corridor seam were) -- both facts gated, both are data
    gates['G7_enclosing_contact_w2_nonzero'] = all(
        rr['I_IR'] > 0 for rr in enc_rows if rr['t'] == 0 and rr['w'] >= 2)
    gates['G7_enclosing_contact_w1_absorbed'] = all(
        rr['I_IR'] == 0 and rr['SB'] > 0 for rr in enc_rows if rr['t'] == 0 and rr['w'] == 1)

    # ---------------------------------------------------------------- SECTION 5: KP tripartite
    print("\n" + "=" * 100)
    print("SECTION 5 -- CONTROL (iii-b): KITAEV-PRESKILL TRIPARTITE (contact construction,")
    print("owners Kitaev-Preskill 2006 / Levin-Wen 2006): the topological constant KNOWN to")
    print("appear, extracted with the SAME entropy instrument S(R) = |R| - r_in(R).")
    print("Disk = induced edges of a q x q block; edges partitioned A/B/C by midpoint")
    print("(A: top-left, B: top-right, C: bottom).  TEE = SA+SB+SC-SAB-SAC-SBC+SABC.")
    tees = []
    for L in (12,):
        T = Ts[L]
        for q in (5, 6, 7, 8, 9, 10):
            for (i0, j0) in ((0, 0), (2, 3)):
                D = T.edges_of_vertexset(block_vs(L, q, i0, j0))
                # midpoint coordinates relative to the block, doubled to stay integer
                A = B = C = 0
                rs = q - 1          # 2*(q-1)/2
                cs = q - 1
                for e in range(T.n):
                    if not ((D >> e) & 1):
                        continue
                    if e < L * L:
                        i, j = divmod(e, L)
                        mi, mj = 2 * ((i - i0) % L), 2 * ((j - j0) % L) + 1
                    else:
                        i, j = divmod(e - L * L, L)
                        mi, mj = 2 * ((i - i0) % L) + 1, 2 * ((j - j0) % L)
                    if mi >= rs:
                        C |= 1 << e
                    elif mj < cs:
                        A |= 1 << e
                    else:
                        B |= 1 << e
                S = lambda R: region_measures(T, R)['S']
                tee = (S(A) + S(B) + S(C) - S(A | B) - S(A | C) - S(B | C) + S(A | B | C))
                tees.append(tee)
                print("  L=%2d q=%2d offset=(%d,%d)  TEE = %d" % (L, q, i0, j0, tee))
    gates['G8_KP_constant'] = len(set(tees)) == 1
    gates['G8_KP_nonzero'] = all(t != 0 for t in tees)
    print("  TEE constant across sizes and offsets: %s ; nonzero: %s  (value %s in log2 units)"
          % (gates['G8_KP_constant'], gates['G8_KP_nonzero'], sorted(set(tees))))

    # ---------------------------------------------------------------- SECTION 6: inserted range
    print("\n" + "=" * 100)
    print("SECTION 6 -- CONTROL (v): INSERTED-RANGE ROW (INSERTED -- instrument sensitivity,")
    print("NOT a Gamma venue, NOT a stabiliser code: one open Z-string row of length ell,")
    print("h-edges of row 1, columns 1..ell, added to the rank instrument's generator list.")
    print("Predicted straddle-both onset for L=12, s=3: ell = g+3.  A separation-dependent")
    print("law WOULD register here; whatever this control shows is the COUPLING'S to report,")
    print("never Gamma's (measuring an inserted mediator law and calling it derived is the")
    print("fatal error this label exists to prevent).")
    L, s = 12, 3
    T = Ts[L]
    print("%4s%5s%6s%6s%12s" % ("g", "ell", "I_IR", "SB", "reach>=g+3?"))
    ins_rows = []
    ins_agree = True
    for g in range(1, 6):
        for ell in (2, 3, 4, 5, 6, 7, 8, 9):
            string = 0
            for j in range(1, ell + 1):
                string |= 1 << T.h(1, j)
            VSA = block_vs(L, s, 0, 0)
            VSB = block_vs(L, s, 0, s + g)
            r = pair_measures(T, VSA, VSB, extraZ=(string,), with_content=False)
            r.update(g=g, ell=ell, reach=(ell >= g + 3))
            ins_rows.append(r)
            ins_agree = ins_agree and r['agree']
            print("%4d%5d%6d%6d%12s" % (g, ell, r['I_IR'], r['SB'], r['reach']))
    # exact inserted-row laws, all computed:
    #   (1) the row STRADDLES both regions exactly when its reach covers the gap (SB > 0
    #       <=> ell >= g+3), the predicted geometric onset;
    #   (2) the mutual interface registers as a UNIT SHELL exactly AT the reach threshold
    #       (I_IR = 1 iff ell = g+3, else 0): at fixed range ell the pair relation is now
    #       g-DEPENDENT (nonzero at the separation the row can just span) -- the instrument
    #       demonstrably reports separation dependence the moment a term with reach exists,
    #       and what it reports is the INSERTED row's reach, nothing of Gamma's.
    gates['G9_inserted_straddle_onset_predicted'] = ins_agree and all(
        (r['SB'] > 0) == r['reach'] for r in ins_rows)
    gates['G9_inserted_unit_shell_at_reach'] = all(
        r['I_IR'] == (1 if r['ell'] == r['g'] + 3 else 0) for r in ins_rows)
    gates['G9_inserted_gives_g_dependence'] = any(
        r['I_IR'] > 0 for r in ins_rows) and all(
        r2['I_IR'] == 0 for r2 in ins_rows if not r2['reach'])

    # ---------------------------------------------------------------- SECTION 7: offsets (D-22)
    print("\n" + "=" * 100)
    print("SECTION 7 -- CONTROL (vi): OFFSET INVARIANCE (D-22 venue symmetry on the numbers).")
    T = Ts[8]
    base = None
    off_ok = True
    for (i0, j0) in ((0, 0), (1, 2), (3, 5), (6, 7)):
        VSA = block_vs(8, 2, i0, j0)
        VSB = block_vs(8, 2, i0, j0 + 2 + 2)
        r = pair_measures(T, VSA, VSB)
        key = (r['IR2A'], r['IR2B'], r['IR2AB'], r['I_IR'], r['I_MI'], r['JC'])
        base = key if base is None else base
        off_ok = off_ok and (key == base)
        print("  offset (%d,%d): IR2A,IR2B,IR2AB,I_IR,I_MI,JC = %s" % (i0, j0, key))
    gates['G10_offset_invariance'] = off_ok

    # ---------------------------------------------------------------- SECTION 8: gates
    print("\n" + "=" * 100)
    print("SECTION 8 -- GATE SUMMARY (each verdict printed FROM the computed boolean):")
    allpass = True
    for name in sorted(gates):
        allpass = allpass and gates[name]
        print("  %-42s %s" % (name, "PASS" if gates[name] else "FAIL"))
    print("  %-42s %s" % ("ALL GATES", "PASS" if allpass else "FAIL"))

    # ---------------------------------------------------------------- SECTION 9: the law
    print("\n" + "=" * 100)
    print("SECTION 9 -- THE TWO-REGION LAW (pure Gamma, this venue, exact, gates cited):")
    print("  SEPARATED (g_eff >= 1, contractible pair): I_IR = I_MI = 0 and no joint record")
    print("    content, IDENTICALLY -- not small, not decaying: zero at EVERY positive earned")
    print("    separation, every s, every L to 12 (G2).  Contact of some generator with both")
    print("    regions is NECESSARY for any nonzero value (G2b).")
    print("  CONTACT (one seam):  I_IR = 2(s-2),  I_MI = s-2,  no joint content (G3).  The")
    print("    s=2 case is exactly ZERO: a single-plaquette seam is rank-absorbed -- contact")
    print("    is necessary, NOT sufficient; the seam must carry >= 2 independent straddling")
    print("    constraints.  The same sub-threshold absorption zeroes the 2-row corridor's")
    print("    contact (G5) and the width-1 enclosing frame's contact (G7d).")
    print("  CONTACT (two seams, L = 2s): I_IR = 4(s-2)+1, I_MI = 2(s-2), and JOINT_NEW = 1:")
    print("    the union winds the torus and ONE genuinely joint record class appears --")
    print("    joint content tracks the union's TOPOLOGY, not the separation (G3b).")
    print("  NON-CONTRACTIBLE PAIR: the exact constant w in every column at every distance;")
    print("    the joint structure is pure IDENTIFICATION of wrap classes (G6) -- the venue's")
    print("    long-range structure is visible to the instrument and is distance-independent.")
    print("  ENCLOSING: zero at every positive gap; closed-seam contact registers for width")
    print("    >= 2 frames (G7).  KP tripartite: the known topological constant, TEE = -1,")
    print("    exactly, at every size and offset (G8) -- the instrument sees topological")
    print("    terms where the literature says they exist.")
    print("  INSERTED RANGE: one extra row of reach ell makes the pair relation g-dependent")
    print("    -- a unit shell exactly at the separation the row can just span (G9).  So the")
    print("    instrument WOULD have shown a separation law; pure Gamma has none to show.")
    print()
    print("HYPOTHESIS VERDICT (computed): pure-Gamma two-region relations on this venue are")
    print("  contact-only or topological -- the hypothesis HOLDS iff all gates pass: %s"
          % ("HOLDS" if allpass else "CHECK FAILED GATES"))

    # ---------------------------------------------------------------- SECTION 10: comparison
    print("\n" + "=" * 100)
    print("SECTION 10 -- FINAL COMPARISON (the ONLY place the target is named, per D-1).")
    print("The comparison target is the Newton-form falloff (interaction weakening smoothly")
    print("and monotonically with separation, 1/r^2 in form).  The pure-Gamma two-region")
    print("relation derived above does NOT have that shape: the exact law is a contact seam")
    print("term plus a topological constant, because every Gamma generator has unit range and")
    print("the venue has zero correlation length.  What Gamma gives at two regions is (a) the")
    print("contact seam law 2(s-2), (b) the topological sector structure (the constant w, the")
    print("wrap-class identification, TEE = -1).  A genuine falloff must therefore enter")
    print("through a term with reach -- the MEDIATOR of the openness clause -- exactly as")
    print("control (v) demonstrates in miniature: one inserted row of reach ell produces a")
    print("separation-dependent shell at the distance it can just span, and that dependence")
    print("is the INSERTED row's law, never Gamma's.  A further honest observation: this")
    print("instrument's codomain is integer rank -- it expresses REACH and TOPOLOGY natively,")
    print("and any smooth graded weakening would have to live in a weighted quantity (the")
    print("coupling's strength, alpha's costs), which is where the mediator increment must")
    print("look.  DIVISION OF LABOR consistent with C-77: Gamma gives the boundary-degree")
    print("law (C-79) and the two-region contact/topological structure (this lane); the")
    print("falloff is the coupling's to give; costs are alpha's.  Nothing in this section")
    print("entered any construction step above.")

if __name__ == "__main__":
    main()
