"""THE GEOMETRY LAYER -- the C-77 increment machinery, folded into the model (T-46).

Everything here is PORTED from the sealed lanes, with fidelity to sealed behavior over
elegance: the exact F_2 bitmask machinery, the writer searches, and the cut-ranks are the
lanes' own algorithms, adapted only as far as packaging requires.  Exact arithmetic
throughout the measurement path: python-int bitmasks over F_2, integer counts, numpy
uint64 popcount scans for exhaustive coset minima.  No floats anywhere except in the
formation section, which reproduces the sealed encoding-level tables (their floats).

WHAT EACH SECTION IMPLEMENTS, WITH ITS CLAIM ROW AND SEALED SOURCE:
  1  exact F_2 kit                              (shared; LANE_T42_A_DISTANCE/t42_lib.py,
                                                 LANE_T42_C_BOUNDCAP/t42_corner.py)
  2  EARNED DISTANCE d_W                 C-78   (LANE_T42_A_DISTANCE/t42_a_toric.py)
  3  INTERFACE RANK + BOUNDARY LAW      C-79   (LANE_T42_C_BOUNDCAP/t42_corner.py,
     corner tier: IR2 = 2*PER-10, capP,          LANE_T42_D_DERIVE/t42d_derive.py)
     32*capP = (IR2+10)^2+8*(IR2+10)-160
  4  WORLD TIER: IFACE = 6n^2,           C-79   (LANE_T42_C_BOUNDCAP/t42_world.py,
     CAP_d1 = n^3-(n-2)^3 = 6n^2-12n+8,          LANE_T42_D_DERIVE/t42d_derive.py)
     IFACE^3 = 216*CAP_total^2
  5  TWO-REGION TAXONOMY                 C-80   (LANE_O54_A_CORNER/o54a_two_region.py,
     I(A:B)=0 separated, seam law                LANE_O54_C_ATTEMPT/o54c_attempts.py +
     I_IR = 2(s-2), winding constant w,          o54c_lib.py)
     coupling-writer law w_min = d
  6  CERTIFIABILITY                      C-81   (LANE_T43_A_CORNER/t43_corner.py)
     CERT = cut-rank identity = 8s-10,
     WRITE0 = 0
  7  THE CERTIFICATION WINDOW LAW        C-83   (LANE_T45_CLOCK/t45_clock.py; per-epoch
     CERT_W(n) = min(n^3, W*(6n^2-12n+8))        min-cut L1 verified in C-82 /
                                                 LANE_T43_B_WORLD)
  8  ENCODING-LEVEL FORMATION        C-71/C-72  (LANE_T32_REALSOURCE/t32_source.py,
     occupancy vs orientation accumulation,      LANE_T34_NAND/t34_nand.py)
     unwritten nulls with declared tolerance

BORROWED INSTRUMENTS, OWNERS NAMED WHERE THE LANES NAMED THEM: the toric carrier is
Kitaev (quant-ph/9707021); cut-rank / stabiliser region entropy is Hamma-Ionicioiu-Zanardi
(2005) and Fattal-Cubitt-Yamamoto-Bravyi-Chuang (quant-ph/0406168) territory; topological
mutual information Kitaev-Preskill (2006) / Levin-Wen (2006); defect records Bravyi-Kitaev
(quant-ph/9811052).  OURS: the exact laws, constants, scopes, and the writer-search
instruments, per the register rows above.

The sealed lanes remain the source of truth; validate_geometry.py gates every number this
module reproduces against its SEALED value.  Returns are DATA."""
import sys as _sys
import os as _os

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

import numpy as _np
from record_model import symplectic_logicals


# =====================================================================================
# 1  EXACT F_2 KIT  (LANE_T42_A_DISTANCE/t42_lib.py and LANE_T42_C_BOUNDCAP/t42_corner.py,
#    verbatim-in-substance)
# =====================================================================================
def pc(x):
    """popcount of a python int."""
    return bin(x).count("1")


def rank_bits(vecs):
    """F_2 rank of a list of bitmask rows.  (t42_corner.py, verbatim.)"""
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
    """A pivot basis of the span.  (t42_corner.py, verbatim.)"""
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
    """Kernel of the F_2 map x -> (row . x)_rows for x supported on `cols` (bit positions);
       kernel elements returned as bitmasks on those columns.  (t42_corner.py, verbatim.)"""
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


def solve_affine_f2(rows, rhs, nvars):
    """Solve M x = b over F_2 on bitmask rows; returns (x0, nullbasis, rank) or
       (None, None, rank) if inconsistent.  (t42_lib.py, verbatim.)"""
    aug = []
    for r, b in zip(rows, rhs):
        aug.append([r, b & 1])
    red = []
    for r, b in aug:
        for pcol, prow, pb in red:
            if (r >> pcol) & 1:
                r ^= prow
                b ^= pb
        if r == 0:
            if b:
                return None, None, len(red)
            continue
        col = r.bit_length() - 1
        nred = []
        for pcol2, prow2, pb2 in red:
            if (prow2 >> col) & 1:
                prow2 ^= r
                pb2 ^= b
            nred.append((pcol2, prow2, pb2))
        red = nred
        red.append((col, r, b))
    rank = len(red)
    x0 = 0
    for col, r, b in red:
        if b:
            x0 |= (1 << col)
    for r, b in zip(rows, rhs):
        assert pc(r & x0) % 2 == (b & 1)
    pivcols = set(c for c, _, _ in red)
    nullbasis = []
    for f in range(nvars):
        if f in pivcols:
            continue
        v = 1 << f
        for col, r, b in red:
            if (r >> f) & 1:
                v |= (1 << col)
        nullbasis.append(v)
    for v in nullbasis:
        for r in rows:
            assert pc(r & v) % 2 == 0
    return x0, nullbasis, rank


def span_numpy(basis):
    """All elements of the span as a numpy uint64 array (masks must fit 64 bits).
       (t42_lib.py, verbatim.)"""
    arr = _np.zeros(1, dtype=_np.uint64)
    for g in basis:
        arr = _np.concatenate([arr, arr ^ _np.uint64(g)])
    return arr


def sp_pair(u, v, n):
    """Symplectic F_2 pairing of (x|z) Paulis given as 2n-bit masks (low n = x part).
       (t42_lib.py / o54c_lib.py, verbatim.)"""
    mask = (1 << n) - 1
    xu, zu = u & mask, u >> n
    xv, zv = v & mask, v >> n
    return (pc(xu & zv) + pc(zu & xv)) % 2


def supp_mask(g, n):
    """Qubit-support mask: union of x and z supports.  (o54c_lib.py, verbatim.)"""
    return (g & ((1 << n) - 1)) | (g >> n)


def vec_to_mask(vec):
    """list of 0/1 -> int mask.  (t42_lib.py, verbatim.)"""
    m = 0
    for i, b in enumerate(vec):
        if b:
            m |= (1 << i)
    return m


# =====================================================================================
# 2  EARNED DISTANCE d_W  --  C-78  (LANE_T42_A_DISTANCE/t42_a_toric.py)
# =====================================================================================
def toric_stabilizers(L):
    """Toric code on the L x L torus, qubits on edges, n = 2 L^2; stars X-type,
       plaquettes Z-type.  (t42_lib.py, verbatim.)"""
    n = 2 * L * L

    def h(i, j):
        return (i % L) * L + (j % L)

    def v(i, j):
        return L * L + (i % L) * L + (j % L)

    stars, plaqs = [], []
    for i in range(L):
        for j in range(L):
            stars.append(vec_to_mask([1 if e in (h(i, j), h(i, j - 1), v(i, j), v(i - 1, j))
                                      else 0 for e in range(n)]))
            plaqs.append(vec_to_mask([1 if e in (h(i, j), h(i + 1, j), v(i, j), v(i, j + 1))
                                      else 0 for e in range(n)]))
    stab = []
    for s in stars:
        stab.append([(s >> e) & 1 for e in range(n)] + [0] * n)
    for p in plaqs:
        stab.append([0] * n + [(p >> e) & 1 for e in range(n)])
    return n, stab, stars, plaqs


_DW_CACHE = {}


def dW_class_matrix(L):
    """C-78: the earned-distance class matrix on the toric torus at side L.

    d_W(s,s') = minimal admissible-writer weight taking configuration s to s'; writers are
    SEARCHED per configuration pair over the FULL affine space of admissible x-type
    operations with that label action (the z=0 reduction verified exhaustively in the
    sealed lane at L=2,3), the minimum weight taken by exact popcount over the whole coset.
    Configurations are the 4 ground labels against the COMPUTED conjugate logical pairs
    (record_model.symplectic_logicals -- computed, never nominated).

    CLAIM ROW: C-78 (LANE_T42_A_DISTANCE, sealed).  Headline: class minima {L, L, 2L};
    d_W = minimal boundary-crossing cost; the metric axioms hold; d_code = L.
    OWNERS: carrier Kitaev quant-ph/9707021; the distance instrument is OURS (C-78).

    Returns dict(dmat={(s,s'): d}, class_min={t: d}, dcode=int, order=[(0,0),(1,0),(0,1),(1,1)]).
    """
    if L in _DW_CACHE:
        return _DW_CACHE[L]
    n, stab, stars, plaqs = toric_stabilizers(L)
    pairs = symplectic_logicals(stab, n)
    assert len(pairs) == 2, "toric torus must yield exactly 2 conjugate logical pairs"
    logi = [(vec_to_mask(X), vec_to_mask(Z)) for X, Z in pairs]
    # verify computed pairs really are conjugate & mutually commuting (the lane's own gate)
    assert all(sp_pair(logi[i][0], logi[i][1], n) == 1 for i in range(2))
    assert all(sp_pair(logi[i][a], logi[j][b], n) == 0
               for i in range(2) for j in range(2) for a in range(2) for b in range(2) if i != j)
    Zlogs = [logi[0][1], logi[1][1]]
    Zl_z = [zl >> n for zl in Zlogs]        # z-parts pair against x-type writers
    plaq_rows = [p for p in plaqs]          # z-part rows of the Z-type stabilizers
    configs = [(a, b) for a in range(2) for b in range(2)]
    dmat, class_min = {}, {}
    for s in configs:
        for sp_ in configs:
            t = (s[0] ^ sp_[0], s[1] ^ sp_[1])
            if t in class_min:
                dmat[(s, sp_)] = class_min[t]
                continue
            rows = plaq_rows + [Zl_z[0], Zl_z[1]]
            rhs = [0] * len(plaq_rows) + [t[0], t[1]]
            x0, nb, _rank = solve_affine_f2(rows, rhs, n)
            assert x0 is not None
            coset = span_numpy(nb) ^ _np.uint64(x0)
            w = int(_np.bitwise_count(coset).min())
            class_min[t] = w
            dmat[(s, sp_)] = w
    dcode = min(class_min[t] for t in class_min if t != (0, 0))
    out = dict(dmat=dmat, class_min=class_min, dcode=dcode,
               order=[(0, 0), (1, 0), (0, 1), (1, 1)])
    _DW_CACHE[L] = out
    return out


def dW(L, s, s2):
    """C-78: earned distance between two configurations s, s2 (each a label pair in
       {0,1}^2) on the L x L toric carrier.  See dW_class_matrix."""
    return dW_class_matrix(L)["dmat"][(tuple(s), tuple(s2))]


# =====================================================================================
# 3  INTERFACE RANK AND THE BOUNDARY LAW (CORNER)  --  C-79
#    (LANE_T42_C_BOUNDCAP/t42_corner.py; LANE_T42_D_DERIVE/t42d_derive.py)
# =====================================================================================
class Toric:
    """Toric carrier, qubits on the edges of an L x L periodic lattice; the C-74/C-79
       conventions (LANE_T42_C_BOUNDCAP/t42_corner.py = LANE_O50_A_ACTION/f2lib.py,
       verbatim-in-substance).  h(i,j): edge (i,j)-(i,j+1), index i*L+j; v(i,j): edge
       (i,j)-(i+1,j), index L*L+i*L+j.  Stars X-type, plaquettes Z-type."""

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
        """Induced edge set: an edge is in R iff BOTH endpoints are in VS (C-74)."""
        L = self.L
        R = 0
        for (i, j) in VS:
            if ((i % L), ((j + 1) % L)) in VS:
                R |= 1 << self.h(i, j)
            if (((i + 1) % L), (j % L)) in VS:
                R |= 1 << self.v(i, j)
        return R

    def straddle_edges(self, VS):
        """Edges with EXACTLY ONE endpoint in VS: the perimeter count."""
        L = self.L
        cnt = 0
        for (i, j) in VS:
            for (di, dj, e) in ((0, 1, self.h(i, j)), (1, 0, self.v(i, j)),
                                (0, -1, self.h(i, (j - 1) % L)), (-1, 0, self.v((i - 1) % L, j))):
                if (((i + di) % L), ((j + dj) % L)) not in VS:
                    cnt += 1
        return cnt


_TORIC_CACHE = {}


def toric(L):
    """Cached Toric carrier."""
    if L not in _TORIC_CACHE:
        _TORIC_CACHE[L] = Toric(L)
    return _TORIC_CACHE[L]


def block_vs(L, s, i0=0, j0=0):
    """s x s vertex block at offset (i0, j0) (the C-74 region family)."""
    return {(((i0 + a) % L), ((j0 + b) % L)) for a in range(s) for b in range(s)}


def interface_rank(L, s, i0=0, j0=0):
    """C-79 corner tier: the cut-rank IR2 and content capacity capP of the s x s block
       region on the L x L toric carrier.

    IR2 = r_S - r_in - r_out, the stabiliser CUT-RANK (independent constraints genuinely
    straddling the cut), computed BOTH ways of the sealed lane (group quotient and
    restriction rank) and gated equal here by assertion.  capP = 2|R| - r_in, the
    distinguishable-content capacity.  On thick blocks (2 <= s <= L-2): PER = 4s,
    IR2 = 2*PER - 10 = 8s - 10, capP = 2s^2 + 2s - 5, and the parameter-free corner form
        32*capP == (IR2+10)^2 + 8*(IR2+10) - 160.

    CLAIM ROW: C-79 (LANE_T42_C_BOUNDCAP LAW-2 + LANE_T42_D_DERIVE R1, sealed).
    BORROWED: the cut-rank instrument is Kitaev / Hamma-Ionicioiu-Zanardi territory
    (stabiliser region entropy).  OURS: the exact laws, constants, and scope.

    Returns dict(AQ, PER, r_in, r_out, IR2, capP)."""
    T = toric(L)
    VS = block_vs(L, s, i0, j0)
    R = T.edges_of_vertexset(VS)
    Rc = T.full ^ R
    AQ = pc(R)
    PER = T.straddle_edges(VS)
    rM = rank_bits(list(T.stars))
    rP = rank_bits(list(T.plaqs))
    r_inX = rM - rank_bits([m & Rc for m in T.stars])
    r_inZ = rP - rank_bits([m & Rc for m in T.plaqs])
    r_outX = rM - rank_bits([m & R for m in T.stars])
    r_outZ = rP - rank_bits([m & R for m in T.plaqs])
    r_in, r_out = r_inX + r_inZ, r_outX + r_outZ
    IR2 = (rM + rP) - r_in - r_out
    # the sealed lane's cross-check: restriction-rank route must agree
    IR2b = (rank_bits([m & R for m in T.stars]) - r_inX) + \
           (rank_bits([m & R for m in T.plaqs]) - r_inZ)
    assert IR2 == IR2b, "cut-rank two ways disagree -- port broken"
    return dict(AQ=AQ, PER=PER, r_in=r_in, r_out=r_out, IR2=IR2, capP=2 * AQ - r_in)


def corner_form_residual(IR2, capP):
    """C-79: residual of the parameter-free corner form; 0 iff the sealed relation
       32*capP == (IR2+10)^2 + 8*(IR2+10) - 160 holds.  (LANE_T42_D_DERIVE, sealed.)"""
    return 32 * capP - ((IR2 + 10) ** 2 + 8 * (IR2 + 10) - 160)


# =====================================================================================
# 4  WORLD TIER  --  C-79  (LANE_T42_C_BOUNDCAP/t42_world.py; LANE_T42_D_DERIVE)
# =====================================================================================
def world_counts(n):
    """C-79 world tier: exact access counts for the n x n x n block of barrier records
       (census model; one record per grain, nearest-neighbour adjacency, access only
       through the block's own boundary).

    IFACE     = block-complement adjacencies (independent access channels) = 6n^2;
    CAP_d1    = grains reachable without traversing another grain = n^3 - (n-2)^3
              = 6n^2 - 12n + 8  (the C-82 per-epoch min-cut law L1(n));
    CAP_total = n^3 (the census volume capacity);
    and the parameter-free world form  IFACE^3 == 216 * CAP_total^2.

    All three are COUNTED here (loops over grains, as the sealed lane did), never assumed
    from the closed forms; the closed forms are what the validator gates.

    CLAIM ROW: C-79 (LANE_T42_C_BOUNDCAP world tier + LANE_T42_D_DERIVE R1, sealed);
    the CAP_d1 closed form is C-82's per-epoch law (LANE_T43_B_WORLD).

    Returns dict(IFACE, CAP_d1, CAP_total)."""
    IFACE = 0
    CAP_d1 = 0
    for x in range(n):
        for y in range(n):
            for z in range(n):
                out = 0
                for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                                   (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                    X, Y, Z = x + dx, y + dy, z + dz
                    if not (0 <= X < n and 0 <= Y < n and 0 <= Z < n):
                        out += 1
                IFACE += out
                if out:
                    CAP_d1 += 1
    return dict(IFACE=IFACE, CAP_d1=CAP_d1, CAP_total=n ** 3)


# =====================================================================================
# 5  TWO-REGION TAXONOMY  --  C-80
#    (LANE_O54_A_CORNER/o54a_two_region.py; LANE_O54_C_ATTEMPT/o54c_attempts.py + o54c_lib.py)
# =====================================================================================
def _sector_stats(rows, R, full):
    """(rk, r_in, r_out, IR2_quotient, IR2_restriction) for one CSS sector.
       (o54a_two_region.py, verbatim.)"""
    Rc = full ^ R
    rk = rank_bits(rows)
    rkR = rank_bits([r & R for r in rows])
    rkRc = rank_bits([r & Rc for r in rows])
    rin = rk - rkRc
    rout = rk - rkR
    return rk, rin, rout, rk - rin - rout, rkR - rin


def _region_measures(T, R):
    """IR2 (both ways, gated) and stabiliser entropy S = |R| - r_in.
       (o54a_two_region.py region_measures, verbatim-in-substance.)"""
    rkX, rinX, routX, ir2X, ir2bX = _sector_stats(T.stars, R, T.full)
    rkZ, rinZ, routZ, ir2Z, ir2bZ = _sector_stats(T.plaqs, R, T.full)
    assert ir2X == ir2bX and ir2Z == ir2bZ
    AQ = pc(R)
    return dict(AQ=AQ, rin=rinX + rinZ, IR2=ir2X + ir2Z, S=AQ - rinX - rinZ)


def _cap_record_with_kernels(T, R):
    """C-74 record-content instrument, with kernels returned for the joint split.
       (o54a_two_region.py, verbatim.)"""
    n = T.n
    cols = [e for e in range(n) if (R >> e) & 1]
    kerX = kernel_basis(T.plaqs, cols)
    kerZ = kernel_basis(T.stars, cols)
    bst = basis_bits(list(T.stars))
    bpl = basis_bits(list(T.plaqs))
    capX = rank_bits(bst + kerX) - len(bst)
    capZ = rank_bits(bpl + kerZ) - len(bpl)
    return capX + capZ, (kerX, kerZ, bst, bpl)


def _joint_content(T, RA, RB):
    """cap_record of A, B, AuB and the union-class split JOINT_NEW / IDENT.
       (o54a_two_region.py, verbatim.)"""
    capA, (kXA, kZA, bst, bpl) = _cap_record_with_kernels(T, RA)
    capB, (kXB, kZB, _, _) = _cap_record_with_kernels(T, RB)
    capAB, (kXAB, kZAB, _, _) = _cap_record_with_kernels(T, RA | RB)
    jn = 0
    ident = 0
    for (base, kA, kB, kAB) in ((bst, kXA, kXB, kXAB), (bpl, kZA, kZB, kZAB)):
        r0 = len(base)
        r_parts = rank_bits(base + kA + kB) - r0
        r_union = rank_bits(base + kAB) - r0
        capA_s = rank_bits(base + kA) - r0
        capB_s = rank_bits(base + kB) - r0
        jn += r_union - r_parts
        ident += (capA_s + capB_s) - r_parts
    return capA, capB, capAB, jn, ident


def _straddle_both(T, RA, RB):
    """COUNT of generators whose support meets BOTH regions."""
    return sum(1 for gset in (T.stars, T.plaqs) for g in gset
               if (g & RA) and (g & RB))


def _pair_measures(T, VSA, VSB):
    """(o54a_two_region.py pair_measures, verbatim-in-substance.)"""
    RA = T.edges_of_vertexset(VSA)
    RB = T.edges_of_vertexset(VSB)
    assert RA & RB == 0, "regions must be disjoint qubit sets"
    mA, mB, mAB = _region_measures(T, RA), _region_measures(T, RB), _region_measures(T, RA | RB)
    capA, capB, capAB, jn, ident = _joint_content(T, RA, RB)
    return dict(IR2A=mA['IR2'], IR2B=mB['IR2'], IR2AB=mAB['IR2'],
                I_IR=mA['IR2'] + mB['IR2'] - mAB['IR2'],
                I_MI=mA['S'] + mB['S'] - mAB['S'],
                SB=_straddle_both(T, RA, RB),
                capA=capA, capB=capB, capAB=capAB,
                JOINT_NEW=jn, IDENT=ident, JC=capAB - capA - capB)


def two_region_blocks(L, s, g):
    """C-80: mutual interface of two s x s blocks at separation g (vertex columns strictly
       between; blocks row-aligned, A at column 0, B at column s+g) on the L x L torus.

    THE SEALED LAW (LANE_O54_A_CORNER): SEPARATED (g_eff >= 1): I_IR = I_MI = 0 and no
    joint record content, identically -- zero at EVERY positive earned separation.
    SINGLE CONTACT (g_eff = 0, L > 2s): I_IR = 2(s-2), I_MI = s-2 (s=2 rank-absorbed to
    exactly zero).  DOUBLE CONTACT (L = 2s): I_IR = 4(s-2)+1, JOINT_NEW = 1.

    CLAIM ROW: C-80 (LANE_O54_A_CORNER, sealed).  BORROWED: stabiliser entropy
    Fattal et al.; topological MI Kitaev-Preskill / Levin-Wen; cut-rank HIZ.
    OURS: the exact seam laws and the contact-only/topological taxonomy.

    Returns the pair-measure dict plus geff."""
    T = toric(L)
    VSA = block_vs(L, s, 0, 0)
    VSB = block_vs(L, s, 0, s + g)
    r = _pair_measures(T, VSA, VSB)
    r['geff'] = min(g, L - 2 * s - g)
    return r


def ring_vs(L, i0, w):
    """Winding ring: w consecutive vertex rows, all columns (non-contractible)."""
    return {(((i0 + a) % L), j) for a in range(w) for j in range(L)}


def winding_pair(L, w, r0):
    """C-80: mutual interface of two non-contractible winding rings of width w with first
       rows 0 and r0 on the L x L torus.

    THE SEALED LAW (LANE_O54_A_CORNER, control (ii)): the venue's own long-range structure
    is the exact CONSTANT w at every distance: I_IR = I_MI = IDENT = w, JOINT_NEW = 0,
    capA = capB = capAB = w -- pure identification of wrap classes, distance-independent.

    CLAIM ROW: C-80 (LANE_O54_A_CORNER, sealed).  Returns the pair-measure dict."""
    T = toric(L)
    return _pair_measures(T, ring_vs(L, 0, w), ring_vs(L, r0, w))


# ---- the O54-C venue (Lx x Ly torus, x,y convention) for the coupling-writer law -------
class Torus:
    """L_x x L_y toric code, o54c_lib.py conventions, verbatim-in-substance.
       h(x,y) = edge (x,y)-(x+1,y), index y*Lx+x; v(x,y) = edge (x,y)-(x,y+1),
       index Lx*Ly + y*Lx + x.  Paulis are 2n-bit ints (low n = x part)."""

    def __init__(self, Lx, Ly):
        self.Lx, self.Ly = Lx, Ly
        self.n = 2 * Lx * Ly

    def h(self, x, y):
        return (y % self.Ly) * self.Lx + (x % self.Lx)

    def v(self, x, y):
        return self.Lx * self.Ly + (y % self.Ly) * self.Lx + (x % self.Lx)

    def star(self, x, y):
        m = 0
        for e in (self.h(x, y), self.h(x - 1, y), self.v(x, y), self.v(x, y - 1)):
            m |= 1 << e
        return m

    def plaq(self, x, y):
        m = 0
        for e in (self.h(x, y), self.h(x, y + 1), self.v(x, y), self.v(x + 1, y)):
            m |= 1 << (self.n + e)
        return m

    def all_stars(self):
        return [self.star(x, y) for y in range(self.Ly) for x in range(self.Lx)]

    def all_plaqs(self):
        return [self.plaq(x, y) for y in range(self.Ly) for x in range(self.Lx)]

    def zbar1(self):
        m = 0
        for y in range(self.Ly):
            m |= 1 << (self.n + self.v(0, y))
        return m

    def zbar2(self):
        m = 0
        for x in range(self.Lx):
            m |= 1 << (self.n + self.h(x, 0))
        return m

    def xbar1(self):
        m = 0
        for y in range(self.Ly):
            m |= 1 << self.h(0, y)
        return m

    def xbar2(self):
        m = 0
        for x in range(self.Lx):
            m |= 1 << self.v(x, 0)
        return m

    def dual_path_x(self, u, vv):
        """X-type connector for removed plaquettes u -> vv: dual path, x first then y.
           (o54c_lib.py, verbatim.)"""
        (ux, uy), (vx, vy) = u, vv
        m = 0
        x, y = ux, uy
        while x != vx:
            nx = (x + 1) % self.Lx
            m ^= 1 << self.v(x + 1, y)
            x = nx
        while y != vy:
            ny = (y + 1) % self.Ly
            m ^= 1 << self.h(x, y + 1)
            y = ny
        return m


def _independent_subset(gens):
    """(o54c_attempts.py, verbatim.)"""
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


def _generator_graph_dist(local_gens, n, gi, gj):
    """BFS distance in the generator-overlap graph.  (o54c_lib.py, verbatim.)"""
    from collections import deque
    supps = [supp_mask(g, n) for g in local_gens]
    N = len(local_gens)
    adj = [[] for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            if supps[i] & supps[j]:
                adj[i].append(j)
                adj[j].append(i)
    dist = {gi: 0}
    dq = deque([gi])
    while dq:
        i = dq.popleft()
        if i == gj:
            return dist[i]
        for j in adj[i]:
            if j not in dist:
                dist[j] = dist[i] + 1
                dq.append(j)
    return None


def coset_min_np(gens, rep, base_bits=22):
    """Exhaustive min weight (and degeneracy, and coset size) over rep * span(gens),
       single-part masks fitting uint64.  (o54c_lib.py, verbatim-in-substance.)"""
    base = gens[:base_bits]
    rest = gens[base_bits:]
    arr = _np.zeros(1, dtype=_np.uint64)
    for g in base:
        arr = _np.concatenate([arr, arr ^ _np.uint64(g)])
    best = None
    cnt = 0
    total = (1 << len(gens))
    for i in range(1 << len(rest)):
        off = rep
        k = i
        j = 0
        while k:
            if k & 1:
                off ^= rest[j]
            k >>= 1
            j += 1
        w = _np.bitwise_count(arr ^ _np.uint64(off))
        m = int(w.min())
        c = int((w == m).sum())
        if best is None or m < best:
            best, cnt = m, c
        elif m == best:
            cnt += c
    return best, cnt, total


def coupling_cost(Lx, Ly, u, v):
    """C-80: the coupling-writer law w_min = d between two hole records.

    Two plaquettes u, v are removed (defect records; owner of the construction
    Bravyi-Kitaev quant-ph/9811052).  The minimal admissible X-type writer coupling the
    two hole records (coset representative: the dual connector path) is scanned
    EXHAUSTIVELY over the full admissible span (independent stars + the two X logicals),
    exactly as the sealed lane did; the earned separation d_gen is BFS distance on the
    generator-overlap graph, computed from generator supports alone.

    THE SEALED LAW (LANE_O54_C_ATTEMPT section C): w_min == d_gen on every pair --
    linear, slope 1, intercept 0, in earned units (confinement character, no falloff).

    CLAIM ROW: C-80 (LANE_O54_C_ATTEMPT, sealed).
    Returns dict(d_gen, w_min, N_min, coset_bits)."""
    T = Torus(Lx, Ly)
    n = T.n
    stars, plaqs = T.all_stars(), T.all_plaqs()
    local = stars + plaqs
    p = T.dual_path_x(u, v)
    rem_plaqs = [T.plaq(x, y) for y in range(Ly) for x in range(Lx) if (x, y) not in (u, v)]
    Bu, Bv = T.plaq(*u), T.plaq(*v)
    # the lane's admissibility gate on the connector
    assert all(sp_pair(p, s0, n) == 0 for s0 in stars)
    assert all(sp_pair(p, q, n) == 0 for q in rem_plaqs)
    assert sp_pair(p, Bu, n) == 1 and sp_pair(p, Bv, n) == 1
    stars_ind = _independent_subset(stars)
    gens_X = stars_ind + [T.xbar1(), T.xbar2()]
    # the lane's independence gate: the X-writer span is independent
    assert rank_bits(list(gens_X)) == len(gens_X)
    # completeness (the lane's gate): every span generator commutes with Bu
    assert all(sp_pair(g, Bu, n) == 0 for g in gens_X)
    xg = [g & ((1 << n) - 1) for g in gens_X]
    wmin, nmin, tot = coset_min_np(xg, p & ((1 << n) - 1))
    plaq_index = {(x, y): y * Lx + x for y in range(Ly) for x in range(Lx)}
    d_gen = _generator_graph_dist(local, n, len(stars) + plaq_index[u],
                                  len(stars) + plaq_index[v])
    return dict(d_gen=d_gen, w_min=wmin, N_min=nmin, coset_bits=tot.bit_length() - 1)


# =====================================================================================
# 6  CERTIFIABILITY  --  C-81  (LANE_T43_A_CORNER/t43_corner.py)
# =====================================================================================
def _subgroup_in(gens, inside, full):
    """Basis of the subgroup of span(gens) supported inside `inside`, by explicit witness
       construction.  (t43_corner.py, verbatim.)"""
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


def certifiability(L, s, i0=0, j0=0):
    """C-81: the certify/write ledger of the s x s block on the L x L toric carrier.

    STORED = capP = 2|R| - r_in (the T42_C capacity, recomputed).
    CERT   = certifiable-from-outside classes = dim( P_R int (S + P_out) / (S int P_R) ),
             computed BOTH sealed routes (sector ranks; direct span intersection with the
             witness-constructed r_in) and gated equal by assertion.
    WRITE0 = ledger classes flippable by admissible writers supported STRICTLY OUTSIDE R
             -- the locality theorem: identically 0 (disjoint supports commute).
    WRITE_B = the same through the boundary layer (qubits of R touched by straddling
             generators).

    THE SEALED IDENTITY (LANE_T43_A_CORNER): CERT == the C-79 interface law 2*PER - 10
    (= 8s - 10 on squares) on every thick rectangle -- certifiable content IS the
    cut-rank.  WRITE0 == 0 on every region.

    CLAIM ROW: C-81 (LANE_T43_A_CORNER, sealed; the boundary-bounded ruling is C-83).
    BORROWED: cleaning of correctable regions is Bravyi-Terhal territory; cut-rank
    entropy HIZ/Kitaev-Preskill.  OURS: the operational certify/write ledger, the exact
    identity by two routes, the locality zero.

    Returns dict(AQ, PER, STORED, CERT, WRITE0, WRITE_B, HIDDEN, r_in)."""
    T = toric(L)
    GX, GZ, n = T.stars, T.plaqs, T.n
    VS = block_vs(L, s, i0, j0)
    R = T.edges_of_vertexset(VS)
    full = (1 << n) - 1
    Rc = full & ~R
    AQ = pc(R)
    PER = T.straddle_edges(VS)
    cols_R = [e for e in range(n) if (R >> e) & 1]
    # sector ranks (route 1)
    aX = rank_bits(list(GX))
    aZ = rank_bits(list(GZ))
    aXR = rank_bits([g & R for g in GX])
    aZR = rank_bits([g & R for g in GZ])
    aXRc = rank_bits([g & Rc for g in GX])
    aZRc = rank_bits([g & Rc for g in GZ])
    r_in = (aX - aXRc) + (aZ - aZRc)
    STORED = 2 * AQ - r_in
    CERT = (aXR + aXRc - aX) + (aZR + aZRc - aZ)
    # route 2: direct span intersection with witness-constructed r_in
    SxR = _subgroup_in(list(GX), R, full)
    SzR = _subgroup_in(list(GZ), R, full)
    r_in_w = len(SxR) + len(SzR)
    assert r_in == r_in_w, "r_in two ways disagree -- port broken"
    span = [g for g in GX] + [(g << n) for g in GZ]
    for e in range(n):
        if not ((R >> e) & 1):
            span.append(1 << e)
            span.append(1 << (n + e))
    d_span = rank_bits(span)
    allv = list(span)
    for e in cols_R:
        allv.append(1 << e)
        allv.append(1 << (n + e))
    d_all = rank_bits(allv)
    assert d_all == 2 * n, "span completeness failed -- port broken"
    CERT_direct = (2 * AQ + d_span - d_all) - r_in_w
    assert CERT == CERT_direct, "CERT two ways disagree -- port broken"
    # the write instrument
    Bmask = 0
    for g in list(GX) + list(GZ):
        if (g & R) and (g & Rc):
            Bmask |= g & R
    cols_O = [e for e in range(n) if ((Rc | Bmask) >> e) & 1]
    cols_O0 = [e for e in range(n) if (Rc >> e) & 1]
    WX_B = kernel_basis(list(GZ), cols_O)
    WZ_B = kernel_basis(list(GX), cols_O)
    WRITE_B = rank_bits([w & R for w in WX_B]) + rank_bits([w & R for w in WZ_B])
    WX_0 = kernel_basis(list(GZ), cols_O0)
    WZ_0 = kernel_basis(list(GX), cols_O0)
    WRITE0 = rank_bits([w & R for w in WX_0]) + rank_bits([w & R for w in WZ_0])
    return dict(AQ=AQ, PER=PER, STORED=STORED, CERT=CERT, WRITE0=WRITE0,
                WRITE_B=WRITE_B, HIDDEN=STORED - CERT, r_in=r_in)


# =====================================================================================
# 7  THE CERTIFICATION WINDOW LAW  --  C-83  (LANE_T45_CLOCK; per-epoch law C-82)
# =====================================================================================
def L1(n):
    """C-82's verified per-epoch min-cut: L1(n) = 6n^2 - 12n + 8 = n^3 - (n-2)^3
       (LANE_T43_B_WORLD, sealed; equals the world CAP_d1 of section 4)."""
    return 6 * n * n - 12 * n + 8


def cert_window(n, W):
    """C-83: the certification-window law CERT_W(n) = min(n^3, W * L1(n)).

    W = tau / t_epoch is the surface's own certification window (clause (ii')'s
    durability tolerance supplies the clock).  The two limits are the sealed headline:
      W = 1         reproduces C-82's per-epoch law L1(n) exactly (whenever L1(n) <= n^3);
      W = infinity  reproduces the volume n^3 -- DEF-A's immortal-record corner
                    (t_m -> infinity), the 'anytime' reading.

    CLAIM ROW: C-83 (LANE_T45_CLOCK, sealed; per-epoch law C-82 / LANE_T43_B_WORLD).
    Pass W = float('inf') for the corner."""
    if W == float('inf'):
        return n ** 3
    return min(n ** 3, int(W) * L1(n))


# =====================================================================================
# 8  ENCODING-LEVEL FORMATION  --  C-71/C-72  (LANE_T32_REALSOURCE; LANE_T34_NAND)
# =====================================================================================
# the sealed lanes' physical constants
E_CHARGE = 1.602176634e-19        # C, per electron (T-34)
N_E = 100                          # electrons per programmed floating gate (T-34)
M_GRAIN = 4.0e5 * 1.26e-24         # A m^2, CoCrPt grain moment M_s * V (T-32)


def orientation_patterns(seed=7, Ns=(10, 100, 1000, 10000, 100000)):
    """C-71/C-72, orientation encoding (T-32, magnetic media): the sealed pattern table.

    Reproduces LANE_T32_REALSOURCE/t32_source.py's EXACT draw order (rng seed 7; per N:
    one draw for random data, one for the AC-erased control), so the sealed screening
    rows reproduce bit-for-bit.  Patterns are +-1 orientations:
      DC-SATURATED   all +1              -> accumulates, ratio 1 (the (d) scope);
      random data    fair +-1 draws      -> SCREENS (~1/sqrt(N));
      DC-free coded  alternating +1,-1   -> screens exactly (sum 0 at even N);
      AC-erased ctrl fair +-1 draws      -> the demagnetised null (screens).

    CLAIM ROW: C-71/C-72 as rescoped by the solidity review (register: one-signedness
    belongs to DC-SATURATED writing; real data tracks screen).  Standard physics
    (saturation remanence) honestly owned; OURS is the encoding-level statement.

    Returns dict[(name, N)] -> numpy array of +-1 (or 0) values."""
    rng = _np.random.default_rng(seed)
    out = {}
    for N in Ns:
        out[("DC-SATURATED", N)] = _np.ones(N)
        out[("random data", N)] = rng.integers(0, 2, N) * 2 - 1
        out[("DC-free coded", N)] = _np.tile([1, -1], N // 2 + 1)[:N]
        out[("AC-erased ctrl", N)] = rng.integers(0, 2, N) * 2 - 1
    return out


def occupancy_patterns(seed=11, N=1000, draws=1000):
    """C-71/C-72, occupancy encoding (T-34, NAND trapped charge): the sealed pattern set.

    Reproduces LANE_T34_NAND/t34_nand.py's EXACT draw order (rng seed 11: `draws` random
    0/1 pages of N cells, THEN the +-5-electron over-erase error page).  Occupancy values
    are 0/1 (a cell is programmed or it is not); the write injects ONE carrier sign
    (electrons), so the accumulated charge of any written page is one-signed: ratio
    |sum|/sum|.| == 1 for EVERY pattern -- mechanism-constitutive.  The unwritten page's
    accumulated quantity is null within the declared +-5e/cell over-erase tolerance.

    CLAIM ROW: C-71/C-72 (LANE_T34_NAND, sealed; PROVED bar's second mechanism).

    Returns dict(written=[draws arrays of 0/1], all_programmed=array of 1s,
                 unwritten_e=array of signed electron counts, tol_e=int)."""
    rng = _np.random.default_rng(seed)
    written = [rng.integers(0, 2, N) for _ in range(draws)]
    err = rng.integers(-5, 6, N)
    return dict(written=written, all_programmed=_np.ones(N),
                unwritten_e=err, tol_e=5 * N)


def accumulation(per_record, pattern):
    """The accumulation instrument shared by both encodings: signed sum, absolute sum,
       and the accumulation-vs-screening discriminator |sum|/sum|.| (C-46 lineage; the
       same object as ProjectModel.configuration).  Returns dict(sum, abs_sum, ratio)."""
    p = _np.asarray(pattern, dtype=float)
    ssum, sabs = per_record * p.sum(), per_record * _np.abs(p).sum()
    return dict(sum=ssum, abs_sum=sabs,
                ratio=(abs(ssum) / sabs if sabs > 0 else None))
