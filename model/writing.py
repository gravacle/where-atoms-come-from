"""THE WRITING TIER -- the C-91 machinery, folded into the model (T-54/T-55, family: writing).

THE PRINCIPAL'S DIRECTIVE (2026-08-21): the URM is the world model -- the framework new
observations are added INTO.  This module homes the writing tier so that a NEW observation
of writing-tier kind (a real surface's dE shifting under adjacency -- C-93's named next
step; a new writing process's conservation bit -- O-58 N3) has a place to ENTER, not just
a museum of the sealed T-48 numbers.

Everything here is PORTED from the sealed lanes, with fidelity to sealed behavior over
elegance, exactly as model/geometry.py did for T-46: exact python ints and Fractions on
every measurement path; floats appear ONLY in the labeled world-entry / comparison
functions at the bottom, run over results already computed exactly.

WHAT EACH SECTION IMPLEMENTS, WITH ITS CLAIM ROW AND SEALED SOURCE (all C-91, FORMAL):
  1  exact Fraction matrix kit                  (LANE_T48_A_DERIVATION/t48a_derivation.py,
                                                 verbatim-in-substance)
  2  KERNEL-TIER VENUES + ENERGY CENSUS  C-91   (LANE_T48_A S1: the D=1 Ising-ring bond
     ring / toric-plaquette / Z_3^3 grain        graph, the D=2 plaquette venue from
     venues, built from carrier supports;        model/geometry.py's toric_stabilizers --
     the moving-sector census dE==0 <=>          the C-78 lineage instrument -- and the
     one-adjacent, zero WITH positive count      D=3 grain venue; sealed census 1024/2048
                                                 and 2304/4608)
  3  KERNEL-TIER ENSEMBLES + VERDICTS    C-91   (LANE_T48_A S0/S3/S4/S5: unitarity forces
     unitarity => conservation for EVERY         conservation for every weighting;
     weighting; conserving <=> critical,         conserving <=> critical as det(I-K) == 0
     det(I-K) == 0 beside exact nonzero          with the exact nonzero CTRL-LEAK dets
     controls; CTRL-BIAS-LINK: conserving        beside it; the biased-but-conserving
     without uniformity, still critical,         counterexample; the polarized-bath
     invariance broken; CTRL-LEAK;               control: bath polarization is bias, and
     CTRL-BATH                                   bias is mass)
  4  CORNER TIER: UNIFORMITY EARNED      C-91   (LANE_T48_B_CORNER: the plaquette venue
     FROM THE WRITER ALGEBRA; mu_c               from carrier supports via geometry.Torus
     LOCATED IN-LANE                             -- o54c lineage; the elementary writer's
     invariant tuple (1,2,2,2,1) identical       invariant tuple computed by symplectic
     on every link => every algebra-             pairing on every link; mu_c = 1/deg by
     measurable ensemble is link-uniform;        Perron row sums + exact sector sandwich +
     the conserving member is unique and         exact resolvent singular AT 1/deg,
     its amplitude is 1/deg = mu_c               nonsingular beside; the weight-1 coset
                                                 stratum (w_min, N_min) == (1, 1) by
                                                 geometry.coset_min_np, exhaustive)
  5  WORLD TIER: THE THREE ENSEMBLES     C-91   (LANE_T48_C_WORLD/t48c_lib.py + driver,
     E1 TRANSPORT: conserving and critical       verbatim-in-substance: torus3 venue,
     (mu = 1/deg = mu_c) at EVERY dE and         kernel_pos / kernel_edge / extract /
     barrier; E2 TRAIL-WITH-RETREAT:             op-identity / drift instruments; the
     conserving at every dE, uniform             closed-form gap
     exactly at dE = 0; E3 TRAIL-WITH-           ln(mu_c/mu) = ln(1 + e^{dE/kT}/l)
     DECAY (the model's own erase                with f0 and E_b dropping out exactly;
     channel): NEVER critical, the mass          sealed anchors: E1 mu = {1/6}; E2 split
     gap in closed form, f0 and E_b drop         (1/7, 2/7) at b = 1/2; E3a mu table and
     out                                         the T44-B comparison row mu(1/2) = 1/8)
  6  THE OBSERVATION GATE (D-25)         C-91   (this module, T-54/T-55: how a NEW
     surface_boltzmann / surface_gap:            surface's number enters -- through
     a new world surface enters through          provenance, then rational brackets, then
     provenance; its gap is bracketed by         TWO FULL EXACT KERNEL COMPUTATIONS, the
     exact kernel computations, never read       closed form CHECKED against the computed
     off the closed form alone;                  ratio at both brackets, never sourced
     kernel_pos_field: the C-93 entry            from it)
     point for site-dependent dE

VERDICTS ARE COMPUTED BOOLEANS; RETURNS ARE DATA.  deg is COMPUTED from each venue's own
rows, never declared; 1/deg is formed from the computed deg; mu_c is located in-lane
(Perron + sector sandwich + resolvent), never imported as a number (D-8).  Every zero
carries a positive control beside it in checks_writing.py (D-15).  World-tier surfaces
enter ONLY through the D-25 provenance gate (surface_boltzmann REFUSES without it).

BORROWED INSTRUMENTS, OWNERS NAMED WHERE THE LANES NAMED THEM: unital-channel /
doubly-stochastic correspondence and Birkhoff's theorem (standard quantum information /
linear algebra; re-verified on explicit operators in the sealed lane); Perron 1907 /
Frobenius 1912 and Gershgorin 1931 (applied only through exact row sums and resolvents on
the venue); stochastic kernels and detailed balance, standard Markov-chain theory (Feller
I; Kolmogorov criterion); Stinespring 1955 (why row sums 1 is the structural property of
an energy-conserving dilation writer -- the row sums themselves are computed, never
cited); Hashimoto 1989 (the non-backtracking operator, used only as the criticality
REFERENCE for the extension-only counting, earned by the same row-sum instrument);
Goldstein 1951 / Kac 1974 (persistent walks, comparison-only remark); the toric carrier
is Kitaev quant-ph/9707021.  OURS: the writer-kernel instrument, the venues, the earned
uniformity, the closed-form gap, per register row C-91.

The sealed lanes remain the source of truth; checks_writing.py gates every number this
module reproduces against its SEALED value, plus API-fidelity probes beyond the gated
range (definition-not-shortcut)."""
import os as _os
import sys as _sys

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

from fractions import Fraction as Fr
from math import comb as _comb

import numpy as _np  # noqa: F401  (geometry's coset scan uses it; kept visible)

from geometry import (toric_stabilizers, Torus, sp_pair, coset_min_np,
                      _independent_subset, rank_bits)

F0, F1 = Fr(0), Fr(1)


# =====================================================================================
# 1  EXACT FRACTION MATRIX KIT  (LANE_T48_A_DERIVATION/t48a_derivation.py,
#    verbatim-in-substance)
# =====================================================================================
def zeros(n):
    return [[F0] * n for _ in range(n)]


def eye(n):
    m = zeros(n)
    for i in range(n):
        m[i][i] = F1
    return m


def msca(c, A):
    return [[c * x for x in row] for row in A]


def msub(A, B):
    n = len(A)
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]


def meq(A, B):
    n = len(A)
    return all(A[i][j] == B[i][j] for i in range(n) for j in range(n))


def mvec(A, v):
    n = len(A)
    return [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]


def mat_row_sums(K):
    return [sum(row) for row in K]


def mat_col_sums(K):
    n = len(K)
    return [sum(K[i][j] for i in range(n)) for j in range(n)]


def is_stochastic(K):
    return all(s == F1 for s in mat_row_sums(K))


def is_doubly_stochastic(K):
    return is_stochastic(K) and all(s == F1 for s in mat_col_sums(K))


def det_exact(M):
    """Exact determinant by Fraction Gaussian elimination (partial pivot on nonzero).
       (t48a_derivation.py, verbatim.)"""
    n = len(M)
    A = [row[:] for row in M]
    det = F1
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return F0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = F1 / A[c][c]
        for r in range(c + 1, n):
            if A[r][c] != 0:
                f = A[r][c] * inv
                for k in range(c, n):
                    A[r][k] -= f * A[c][k]
    return det


# =====================================================================================
# 2  KERNEL-TIER VENUES + THE ENERGY CENSUS  --  C-91  (LANE_T48_A S1)
#    Venues are built FROM CARRIER SUPPORTS; deg is read off the rows downstream.
# =====================================================================================
def _maskset(fs):
    m = 0
    for b in fs:
        m |= 1 << b
    return m


def ring_venue(N):
    """The D=1 Ising-ring venue: N sites, bonds b=(b, b+1); H = -sum Z_b Z_{b+1}; walls
       (defects) live on bonds; the venue is the bond graph through shared sites.
       CLAIM ROW: C-91 kernel tier (LANE_T48_A_DERIVATION S1, sealed; N=8 in-lane).
       Returns dict(adj=list[set], n, deltas=[int mask del(X_i)], sites_per_bond_ok)."""
    bond_sites = [(b, (b + 1) % N) for b in range(N)]
    site_bonds = [[] for _ in range(N)]
    for b, (s1, s2) in enumerate(bond_sites):
        site_bonds[s1].append(b)
        site_bonds[s2].append(b)
    well_formed = all(len(bs) == 2 for bs in site_bonds)
    deltas = [_maskset(site_bonds[i]) for i in range(N)]
    adj = [set() for _ in range(N)]
    for i in range(N):
        a, b = sorted(site_bonds[i])
        adj[a].add(b)
        adj[b].add(a)
    return dict(adj=adj, n=N, deltas=deltas, well_formed=well_formed)


def ring_syndrome(z, N):
    """Domain-wall syndrome of Z-basis state z on the N-ring.  (t48a, verbatim.)"""
    return sum(((z >> b) & 1) ^ ((z >> ((b + 1) % N)) & 1) and (1 << b) for b in range(N))


def plaquette_venue(L):
    """The D=2 kernel-tier venue: plaquette adjacency of the L x L toric carrier, rebuilt
       from plaquette SUPPORTS (one venue link = one shared carrier edge = one
       grain-boundary crossing, the Gamma price).  Carrier instrument:
       geometry.toric_stabilizers (C-78 lineage; Kitaev quant-ph/9707021).
       CLAIM ROW: C-91 kernel tier (LANE_T48_A_DERIVATION S1, sealed; L=3,4 in-lane).
       Returns dict(adj, n, deltas=[int, per elementary writer], n_writers, stars, plaqs,
       edge_in_two, shared_single)."""
    n_e, _stab, stars, plaqs = toric_stabilizers(L)
    edge_plaqs = [[p for p in range(len(plaqs)) if (plaqs[p] >> e) & 1] for e in range(n_e)]
    edge_in_two = all(len(ep) == 2 for ep in edge_plaqs)
    npl = len(plaqs)
    adj = [set() for _ in range(npl)]
    shared = {}
    for e, pq in enumerate(edge_plaqs):
        if len(pq) == 2:
            p, q = pq
            adj[p].add(q)
            adj[q].add(p)
            shared[(p, q)] = shared.get((p, q), 0) + 1
    deltas = []
    for e in range(n_e):
        if len(edge_plaqs[e]) == 2:
            p, q = edge_plaqs[e]
            deltas.append((1 << p) | (1 << q))
    return dict(adj=adj, n=npl, deltas=deltas, n_writers=len(deltas), stars=stars,
                plaqs=plaqs, edge_in_two=edge_in_two,
                shared_single=all(v == 1 for v in shared.values()))


def grain_venue(L=3):
    """The D=3 kernel-tier venue: Z_L^3 grains, face adjacency (T44-B world lineage).
       CLAIM ROW: C-91 kernel tier (LANE_T48_A_DERIVATION S1, sealed; L=3 in-lane).
       Returns dict(adj, n, transpose=the axis-transposing generator (i,j,k)->(j,i,k),
       used by the CTRL-BIAS-LINK invariance witness)."""
    g_idx = {(i, j, k): i * L * L + j * L + k
             for i in range(L) for j in range(L) for k in range(L)}
    n = L ** 3
    adj = [set() for _ in range(n)]
    for (i, j, k), gi in g_idx.items():
        for ax, dd in ((0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1)):
            c = [i, j, k]
            c[ax] = (c[ax] + dd) % L
            adj[gi].add(g_idx[tuple(c)])
    g_lst = sorted(g_idx, key=g_idx.get)
    transpose = [g_idx[(j, i, k)] for (i, j, k) in g_lst]
    return dict(adj=adj, n=n, transpose=transpose, g_idx=g_idx)


def plaquette_transpose(L):
    """The axis-transposing generator of the plaquette venue: cell (i,j) -> (j,i)
       (toric_stabilizers plaquette index i*L+j).  Used by the CTRL-BIAS-LINK
       invariance witness (LANE_T48_A S4.7)."""
    return [j * L + i for i in range(L) for j in range(L)]


def check_aut(adj, gmap):
    """Does the vertex map gmap preserve adjacency exactly?  (t48a, verbatim.)"""
    n = len(adj)
    return all(sorted(gmap[y] for y in adj[x]) == sorted(adj[gmap[x]]) for x in range(n))


def venue_degree(adj):
    """The venue's degree, COMPUTED from its own rows (D-8: never declared).
       Returns the constant if the rows agree, else None."""
    degs = {len(x) for x in adj}
    return degs.pop() if len(degs) == 1 else None


def moving_census(deltas, nbits, even_only):
    """The energy-sector census of the induced writer action on syndromes:
       dE = 2*(popcount(s XOR del) - popcount(s)) under H = -sum(checks); the sealed
       identity is dE == 0 EXACTLY on the defect-moving sector (one check adjacent) and
       dE in {-4, 0, +4}.  CLAIM ROW: C-91 kernel tier (LANE_T48_A S1.4/S1.11, sealed:
       D=1 census 1024 of 2048; D=2 L=3 census 2304 of 4608).
       Returns dict(sectors_ok, moving, total)."""
    ok = True
    moving = 0
    total = 0
    for s in range(1 << nbits):
        if even_only and bin(s).count("1") % 2:
            continue
        ps = bin(s).count("1")
        for d in deltas:
            dE = 2 * (bin(s ^ d).count("1") - ps)
            one_adj = bin(s & d).count("1") == 1
            ok &= (dE == 0) == one_adj and dE in (-4, 0, 4)
            moving += 1 if dE == 0 else 0
            total += 1
    return dict(sectors_ok=ok, moving=moving, total=total)


# =====================================================================================
# 3  KERNEL-TIER ENSEMBLES AND VERDICTS  --  C-91  (LANE_T48_A S0/S3/S4/S5)
#    Unitarity forces conservation for EVERY weighting; conserving <=> critical.
# =====================================================================================
def kernel_uniform(adj, c):
    """The invariant lazy family c*I + ((1-c)/deg)*A; c = 0 is E-LOC.  deg computed
       from the rows.  (t48a S3/S4, verbatim-in-substance.)"""
    n = len(adj)
    deg = venue_degree(adj)
    t = (F1 - Fr(c)) / deg
    K = zeros(n)
    for x in range(n):
        K[x][x] = Fr(c)
        for y in adj[x]:
            K[x][y] = t
    return K


def link_amplitudes(adj, K):
    """The set of per-link kernel entries over the venue's own links."""
    return {K[x][y] for x in range(len(adj)) for y in adj[x]}


def crit_det(K):
    """The exact criticality witness det(I - K); the sealed identity is 0 for every
       measure-conserving kernel, with exact nonzero controls beside (CTRL-LEAK)."""
    return det_exact(msub(eye(len(K)), K))


def leak_kernel(adj, survival):
    """CTRL-LEAK: the declared-survival leaky ensemble survival * E-LOC.
       CLAIM ROW: C-91 kernel tier (LANE_T48_A S5, sealed: survival 9/10; the exact
       nonzero dets are the sealed anchors gated in checks_writing.py)."""
    return msca(Fr(survival), kernel_uniform(adj, F0))


def biased_kernel_2d(L, a, b):
    """CTRL-BIAS-LINK on the L x L plaquette venue: amplitude a on x-links, b on
       y-links; conserving iff 2a+2b == 1 (computed downstream, never assumed).
       CLAIM ROW: C-91 kernel tier (LANE_T48_A S4.4-S4.7, sealed at L=3, a=1/3, b=1/6:
       doubly stochastic, NON-uniform, still critical, invariance broken exactly at the
       axis-transposing generator)."""
    a, b = Fr(a), Fr(b)
    n = L * L
    K = zeros(n)
    for i in range(L):
        for j in range(L):
            x = i * L + j
            for di, w in (((1, 0), a), ((-1, 0), a), ((0, 1), b), ((0, -1), b)):
                y = ((i + di[0]) % L) * L + ((j + di[1]) % L)
                K[x][y] += w
    return K


def biased_kernel_3d(L, ws):
    """CTRL-BIAS-LINK on the Z_L^3 grain venue: per-axis amplitudes ws (3 of them).
       CLAIM ROW: C-91 kernel tier (LANE_T48_A S4.4-S4.7, sealed at L=3,
       ws = (1/4, 1/6, 1/12))."""
    ws = [Fr(w) for w in ws]
    V = grain_venue(L)
    K = zeros(V["n"])
    for (i, j, k), x in V["g_idx"].items():
        for ax in range(3):
            for dd in (1, -1):
                c = [i, j, k]
                c[ax] = (c[ax] + dd) % L
                K[x][V["g_idx"][tuple(c)]] += ws[ax]
    return K


def invariance_violations(adj, K, gmap):
    """Count of kernel entries NOT invariant under the venue map gmap: the computed
       witness that ensemble-weight symmetry is extra data (LANE_T48_A S4.7).  The
       biased-but-conserving kernel registers a POSITIVE count here while the uniform
       kernel registers the exact zero beside it (D-15 pairing in checks_writing.py)."""
    return sum(1 for x in range(len(adj)) for y in adj[x]
               if K[x][y] != K[gmap[x]][gmap[y]])


def glob_kernel(nsynd_bits, deltas, weights, extra_trivial=0):
    """The global writer ensemble on the FULL (even-parity) syndrome space:
       T[s][s XOR delta] += w; trivial-action writers add weight to the diagonal.
       THE SEALED IDENTITY (C-91 kernel tier, LANE_T48_A S3.8-S3.12): the induced kernel
       is doubly stochastic for EVERY weight row summing to 1 -- unitarity (involutive
       permutation action), never uniformity, conserves measure.
       Returns (synds, T) with T dict-of-dict."""
    synds = [s for s in range(1 << nsynd_bits) if bin(s).count("1") % 2 == 0]
    T = {s: {} for s in synds}
    for s in synds:
        if extra_trivial:
            T[s][s] = T[s].get(s, F0) + Fr(extra_trivial)
        for d, w in zip(deltas, weights):
            t = s ^ d
            T[s][t] = T[s].get(t, F0) + Fr(w)
    return synds, T


def dict_doubly_stochastic(synds, T):
    """Exact double stochasticity of a dict-of-dict kernel.  (t48a, verbatim.)"""
    rows_ok = all(sum(T[s].values()) == F1 for s in synds)
    cols = {s: F0 for s in synds}
    for s in synds:
        for t, w in T[s].items():
            cols[t] += w
    return rows_ok and all(v == F1 for v in cols.values())


def pair_sector(T, adj, p_track, q_hold):
    """The pair-sector reading of the global ensemble (the T-44 connecting-string
       model's per-crossing mu): starting from the 2-defect configuration
       {p_track, q_hold}, decompose one step into tracked-end motion / origin-end
       motion / pair creation.  CLAIM ROW: C-91 kernel tier (LANE_T48_A S3.13-S3.16,
       sealed on T3: split 4/18 + 4/18 + 10/18; conditional per-link 1/4).
       Returns dict(tracked={cell: amp}, origin_sum, creation, conditional=set,
       tracked_on_neighbors, nonadjacent)."""
    s_pair = (1 << p_track) | (1 << q_hold)
    row = T[s_pair]
    tracked, origin_side, creation = {}, {}, F0
    for t, w in row.items():
        if bin(t).count("1") == 2 and (t >> q_hold) & 1:
            pt = (t & ~(1 << q_hold)).bit_length() - 1
            if pt != p_track:
                tracked[pt] = w
        elif bin(t).count("1") == 2 and (t >> p_track) & 1:
            qt = (t & ~(1 << p_track)).bit_length() - 1
            origin_side[qt] = w
        else:
            creation += w
    act = sum(tracked.values())
    cond = {v / act for v in tracked.values()} if act else set()
    return dict(tracked=tracked, origin_sum=sum(origin_side.values()), creation=creation,
                conditional=cond, tracked_on_neighbors=set(tracked) == adj[p_track],
                nonadjacent=p_track not in adj[q_hold])


def bath_dilation(p):
    """CTRL-BATH: the iv' dilation route at its smallest honest size -- system qubit +
       bath qubit, U = exchange (|1,0> <-> |0,1|, else fix), bath traced with weights
       (p, 1-p).  CLAIM ROW: C-91 kernel tier (LANE_T48_A S5.5-S5.8, sealed: column
       sums [1,1] at p=1/2, [3/2,1/2] at 3/4, [2,0] at 1 -- a polarized bath is a
       non-measure-conserving writer ensemble; bath polarization is bias, and bias is
       mass).  Returns dict(unitary, conserves_excitation, trace_preserving,
       doubly_stochastic, col_sums)."""
    p = Fr(p)
    U = [[F0] * 4 for _ in range(4)]  # basis |s b>: 00, 01, 10, 11
    U[0][0] = F1
    U[3][3] = F1
    U[1][2] = F1
    U[2][1] = F1
    UtU = [[sum(U[k][i] * U[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    exc = [0, 1, 1, 2]
    img = [0, 2, 1, 3]
    T = [[F0, F0], [F0, F0]]
    for s in (0, 1):
        for bq, wb in ((0, p), (1, F1 - p)):
            out = img[2 * s + bq]
            T[s][out >> 1] += wb
    return dict(unitary=meq(UtU, eye(4)),
                conserves_excitation=all(exc[img[i]] == exc[i] for i in range(4)),
                trace_preserving=all(x == F1 for x in mat_row_sums(T)),
                doubly_stochastic=is_doubly_stochastic(T),
                col_sums=mat_col_sums(T))


# =====================================================================================
# 4  CORNER TIER: UNIFORMITY EARNED FROM THE WRITER ALGEBRA; mu_c IN-LANE  --  C-91
#    (LANE_T48_B_CORNER; carrier + coset instruments from geometry.py, o54c lineage)
# =====================================================================================
def corner_venue(Lx, Ly):
    """The corner venue REBUILT FROM CARRIER SUPPORTS: the plaquette dual lattice of the
       Lx x Ly toric carrier (geometry.Torus, o54c conventions), one link per shared
       carrier edge.  CLAIM ROW: C-91 corner tier (LANE_T48_B_CORNER S0, sealed at
       (4,6), (3,7), (5,5): #links == #carrier edges, all multiplicities 1, deg = 4
       computed, no self-links).
       Returns dict(T, n, cells, idx, rows=[(j, mult) adjacency rows], links,
       edge_of_link, xg=the admissible X-writer span for the coset instrument,
       stars, edge_in_two, mult_all_one, no_self)."""
    T = Torus(Lx, Ly)
    n = T.n
    cells = [(x, y) for y in range(Ly) for x in range(Lx)]
    idx = {c: i for i, c in enumerate(cells)}
    nc = len(cells)
    pe = [set() for _ in range(nc)]
    for ci, (x, y) in enumerate(cells):
        m = T.plaq(x, y) >> n
        e = 0
        while m:
            if m & 1:
                pe[ci].add(e)
            m >>= 1
            e += 1
    edge_in = [0] * n
    for es in pe:
        for e in es:
            edge_in[e] += 1
    links = []
    edge_of_link = {}
    for i in range(nc):
        for j in range(i + 1, nc):
            common = pe[i] & pe[j]
            if common:
                links.append((i, j, len(common)))
                edge_of_link[(i, j)] = sorted(common)
    rows = [[] for _ in range(nc)]
    for (i, j, m) in links:
        rows[i].append((j, m))
        rows[j].append((i, m))
    stars = T.all_stars()
    stars_ind = _independent_subset(stars)
    gens_X = stars_ind + [T.xbar1(), T.xbar2()]
    assert rank_bits(list(gens_X)) == len(gens_X)
    xg = [g & ((1 << n) - 1) for g in gens_X]
    return dict(T=T, n=n, cells=cells, idx=idx, rows=rows, links=links,
                edge_of_link=edge_of_link, xg=xg, stars=stars, plaq_edges=pe,
                edge_in_two=all(c == 2 for c in edge_in),
                mult_all_one=all(m == 1 for (_i, _j, m) in links),
                no_self=all(all(j != i for (j, _m) in rows[i]) for i in range(nc)))


def chain_venue(L):
    """The D=1 chain venue C_L (LANE_T48_B's own D-15 discriminator).  Rows of
       (neighbor, multiplicity)."""
    return [[((i - 1) % L, 1), ((i + 1) % L, 1)] for i in range(L)]


def writer_invariants(cv):
    """THE WRITER ALGEBRA'S OWN INVARIANTS of the elementary writer on EVERY link:
       (Gamma price, |syndrome|, #stars containing its edge, #plaquettes containing its
       edge, link multiplicity), each computed by popcount / symplectic pairing /
       membership -- never declared.  THE SEALED IDENTITY (C-91 corner tier,
       LANE_T48_B_CORNER S3): the tuple is (1, 2, 2, 2, 1) IDENTICALLY on every link,
       so every algebra-measurable ensemble is link-uniform -- uniformity EARNED.
       Returns list of tuples, one per link."""
    T, n, cells = cv["T"], cv["n"], cv["cells"]
    out = []
    for (i, j, m) in cv["links"]:
        e = cv["edge_of_link"][(i, j)][0]
        wmask = 1 << e
        price = bin(wmask).count("1")
        syn = sum(1 for (x, y) in cells if sp_pair(wmask, T.plaq(x, y), n) == 1)
        in_stars = sum(1 for s in cv["stars"] if (s >> e) & 1)
        in_plaqs = sum(1 for es in cv["plaq_edges"] if e in es)
        out.append((price, syn, in_stars, in_plaqs, m))
    return out


def elementary_coset(cv, link):
    """(w_min, N_min, coset bits) of the elementary writer's admissible coset, scanned
       EXHAUSTIVELY (geometry.coset_min_np, o54c lineage).  THE SEALED IDENTITY
       (C-91 corner tier, LANE_T48_B S0): (w_min, N_min) == (1, 1) -- the weight-1
       stratum of the admissible-writer coset IS the link, in weight AND count."""
    e = cv["edge_of_link"][tuple(sorted(link))][0]
    wmin, nmin, tot = coset_min_np(cv["xg"], 1 << e)
    return wmin, nmin, tot.bit_length() - 1


def resolvent_rows(rows, mu, src):
    """EXACT rational solve of (I - mu*A) x = e_src on (neighbor, mult) adjacency rows;
       None iff singular (zero pivot column).  (t48b resolvent, verbatim-in-substance.)"""
    n = len(rows)
    mu = Fr(mu)
    M = [[F1 if i == j else F0 for j in range(n)] for i in range(n)]
    for i in range(n):
        for (j, m) in rows[i]:
            M[i][j] -= mu * m
    rhs = [F0] * n
    rhs[src] = F1
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] != 0), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        rhs[col], rhs[piv] = rhs[piv], rhs[col]
        inv = F1 / M[col][col]
        M[col] = [x * inv for x in M[col]]
        rhs[col] *= inv
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [xr - f * xc for xr, xc in zip(M[r], M[col])]
                rhs[r] -= f * rhs[col]
    return rhs


def mu_c_locate(rows, src, beside):
    """mu_c LOCATED IN-LANE, never imported: Perron row-sum sandwich (constant row sums
       deg => spectral radius deg exactly, Gershgorin closing the bound from above) plus
       the exact resolvent SINGULAR at 1/deg and NONSINGULAR at every declared point
       beside it.  CLAIM ROW: C-91 corner tier (LANE_T48_B S1, sealed: 1/4 on the
       plaquette venues, 1/2 on the chain).
       Returns dict(deg, rows_constant, mu_c, singular_at_mu_c, nonsingular_beside)."""
    sums = {sum(m for (_j, m) in r) for r in rows}
    rows_constant = len(sums) == 1
    deg = sums.pop() if rows_constant else None
    mu_c = Fr(1, deg) if deg else None
    return dict(deg=deg, rows_constant=rows_constant, mu_c=mu_c,
                singular_at_mu_c=(resolvent_rows(rows, mu_c, src) is None) if deg else None,
                nonsingular_beside=[resolvent_rows(rows, mb, src) is not None
                                    for mb in beside])


def sector_sandwich(M=300, IND=10000):
    """The exact sector-sandwich lemmas locating mu_c = 1/4 in the venue limit (Z^2) and
       1/2 on the D=1 chain: 16^m/(2m+1)^2 <= C(2m,m)^2 <= 16^m and C(2m,m)(2m+1) >= 4^m,
       with exact induction steps to IND.  (t48b S1, verbatim.)  Returns dict of booleans."""
    up_ind = all((2 * m + 1) * (2 * m + 3) <= (2 * m + 2) ** 2 for m in range(1, IND + 1))
    upper = all(_comb(2 * m, m) ** 2 * (2 * m + 1) <= 16 ** m for m in range(1, M + 1))
    lower = all(_comb(2 * m, m) * (2 * m + 1) >= 4 ** m for m in range(1, M + 1))
    mono = all((2 * m + 1) ** 2 * (m + 1) >= 4 * m * (m + 1) ** 2 for m in range(1, IND + 1))
    return dict(upper=upper, upper_induction=up_ind, lower=lower, lower_induction=mono,
                all_ok=upper and up_ind and lower and mono)


def conserving_member(rows):
    """The unique measure-conserving member of the link-uniform (algebra-measurable)
       ensemble class: deg * t == 1 with deg COMPUTED from the rows.  THE SEALED
       IDENTITY (C-91 corner tier, LANE_T48_B S4 K1): t* = 1/deg == mu_c -- the induced
       per-link amplitude of the conserving ensemble IS the critical one.
       Returns dict(deg, t_star, unique=deg*t_star == 1)."""
    loc = mu_c_locate(rows, 0, beside=())
    t_star = Fr(1, loc["deg"]) if loc["deg"] else None
    return dict(deg=loc["deg"], t_star=t_star,
                unique=(t_star is not None and t_star * loc["deg"] == 1))


# =====================================================================================
# 5  WORLD TIER: THE THREE ENSEMBLES  --  C-91  (LANE_T48_C_WORLD/t48c_lib.py, ported)
# =====================================================================================
DIRS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
OPP = (1, 0, 3, 2, 5, 4)


def torus3(n):
    """The census access geometry's grain lattice closed up: periodic n x n x n cells,
       adjacent <=> share a face.  (t48c_lib.py, verbatim.)
       Returns cells, index, nbr[i][d]."""
    cells = [(x, y, z) for z in range(n) for y in range(n) for x in range(n)]
    idx = {c: i for i, c in enumerate(cells)}
    nbr = []
    for (x, y, z) in cells:
        row = []
        for (dx, dy, dz) in DIRS:
            row.append(idx[((x + dx) % n, (y + dy) % n, (z + dz) % n)])
        nbr.append(row)
    return cells, idx, nbr


def bfs_dist(nbr, src):
    """BFS distance on the venue.  (t48c_lib.py, verbatim.)"""
    from collections import deque
    dist = {src: 0}
    dq = deque([src])
    while dq:
        i = dq.popleft()
        for j in nbr[i]:
            if j not in dist:
                dist[j] = dist[i] + 1
                dq.append(j)
    return dist


def l1_wrap(c, n):
    return sum(min(v, n - v) for v in c)


def kernel_pos(nbr, diag, link):
    """Position-space one-writer-step kernel: W[i][i] = diag, W[i][j] = link per
       face-neighbor.  Dict rows, exact Fractions.  (t48c_lib.py, verbatim.)"""
    W = []
    for i, row in enumerate(nbr):
        r = {i: Fr(diag)}
        for j in row:
            r[j] = r.get(j, 0) + Fr(link)
        W.append(r)
    return W


def kernel_pos_field(nbr, diag_of, link_of):
    """THE C-93 ENTRY POINT (T-54/T-55; named, NOT run here as physics): the
       position-space writer kernel with SITE-DEPENDENT amplitudes -- diag_of(i) and
       link_of(i, d) exact Fractions -- where 'does one written record shift dE for an
       adjacent write?' (C-93's named next step, C-91's own currency: bias is mass) will
       build its kernel from a written pattern's dE field.  At constant rates this IS
       kernel_pos (gated in checks_writing.py); conservation and criticality of any
       field kernel are read by the SAME instruments (srow_sums, extract_pos,
       apply_forward, the exact resolvent/determinant kit above)."""
    W = []
    for i, row in enumerate(nbr):
        r = {i: Fr(diag_of(i))}
        for d, j in enumerate(row):
            r[j] = r.get(j, 0) + Fr(link_of(i, d))
        W.append(r)
    return W


def kernel_edge(nbr, stay, fresh, back):
    """Directed-edge (tip-state) kernel: state (i, d) = tip at grain i, last step along
       DIRS[d]; continue fresh (e != OPP[d]) with `fresh`, retreat with `back`
       (back=None DROPS the retreat channel: that measure is lost), stay with `stay`.
       (t48c_lib.py, verbatim.)"""
    S = []
    for i, row in enumerate(nbr):
        for d in range(6):
            s = i * 6 + d
            r = {s: Fr(stay)}
            for e in range(6):
                t = row[e] * 6 + e
                if e == OPP[d]:
                    if back is not None:
                        r[t] = r.get(t, 0) + Fr(back)
                else:
                    r[t] = r.get(t, 0) + Fr(fresh)
            S.append(r)
    return S


def nb_edge_adjacency(nbr):
    """The venue's non-backtracking directed-edge adjacency (Hashimoto 1989; the
       criticality REFERENCE for extension-only counting, earned by row sums).
       (t48c_lib.py, verbatim.)"""
    S = []
    for i, row in enumerate(nbr):
        for d in range(6):
            r = {}
            for e in range(6):
                if e != OPP[d]:
                    t = row[e] * 6 + e
                    r[t] = r.get(t, 0) + 1
            S.append(r)
    return S


def srow_sums(W):
    return [sum(r.values()) for r in W]


def scol_sums(W, n):
    cs = [F0] * n
    for r in W:
        for j, v in r.items():
            cs[j] += v
    return cs


def apply_forward(W, vec):
    """Measure evolution new[j] = sum_i vec[i] W[i][j], exact.  (t48c_lib.py, verbatim.)"""
    new = [F0] * len(W)
    for i, v in enumerate(vec):
        if v:
            for j, w in W[i].items():
                new[j] += v * w
    return new


def extract_pos(W, nbr):
    """Induced per-link amplitude, position space: I - W = c (I - sum m S); returns
       (set of c, set of m).  Uniformity is GATED downstream, never inserted.
       (t48c_lib.py, verbatim.)"""
    cs, ms = set(), set()
    for i, r in enumerate(W):
        c = 1 - r[i]
        cs.add(c)
        for j in nbr[i]:
            ms.add(r[j] / c)
    return cs, ms


def extract_edge(W, nbr):
    """Induced per-link amplitudes, directed-edge space: (set of c, set of m_fresh,
       set of m_back, set of per-state redistribution totals).  (t48c_lib.py, verbatim.)"""
    cs, mf, mb, tots = set(), set(), set(), set()
    for i, row in enumerate(nbr):
        for d in range(6):
            s = i * 6 + d
            r = W[s]
            c = 1 - r.get(s, F0)
            cs.add(c)
            tot = F0
            for e in range(6):
                t = row[e] * 6 + e
                a = r.get(t, F0)
                if e == OPP[d]:
                    if a:
                        mb.add(a / c)
                    tot += a / c
                else:
                    mf.add(a / c)
                    tot += a / c
            tots.add(tot)
    return cs, mf, mb, tots


def op_identity_pos(W, nbr, c, m):
    """Entrywise gate: I - W == c (I - m A) on the venue.  (t48c_lib.py, verbatim.)"""
    for i, r in enumerate(W):
        if 1 - r[i] != c:
            return False
        row_counts = {}
        for j in nbr[i]:
            row_counts[j] = row_counts.get(j, 0) + 1
        for j, mult in row_counts.items():
            if -r[j] != -c * m * mult:
                return False
        if set(r.keys()) - {i} - set(row_counts.keys()):
            return False
    return True


def drift_per_state(W, nbr):
    """Exact one-step expected displacement of each edge state.  (t48c_lib.py, verbatim.)"""
    out = {}
    for i, row in enumerate(nbr):
        for d in range(6):
            s = i * 6 + d
            v = [F0] * 3
            for t, a in W[s].items():
                if t != s:
                    e = t % 6
                    for ax in range(3):
                        v[ax] += a * DIRS[e][ax]
            out[s] = tuple(v)
    return out


def trail_energetics(nbr):
    """The configuration-energy grounding of the world ensembles: E_config = dE *
       (#written grains), delta-counts computed as exact integers on an EXPLICIT written
       trail, never symbolic.  THE SEALED IDENTITY (C-91 world tier, LANE_T48_C_WORLD
       G07/G08): a TRANSPORT step changes the count by 0 in EVERY direction -- the
       computed zero that (with venue symmetry) earns E1's direction-uniform amplitude
       -- while a trail WRITE changes it by +1 (every fresh direction) and the ERASE by
       -1 (the D-15 nonzero beside the zero).
       Returns dict(transport=[dN x 6], write=[dN per fresh direction], erase=[dN])."""
    x0 = 0
    a = nbr[x0][0]
    b = nbr[a][0]
    c = nbr[b][2]
    trail = [x0, a, b, c]
    x_tip = trail[-1]
    dN_transport, dN_write, dN_erase = [], [], []
    for d in range(len(nbr[0])):
        y = nbr[x0][d]
        dN_transport.append(len({y}) - len({x0}))
        yt = nbr[x_tip][d]
        if yt not in trail:
            dN_write.append(len(set(trail) | {yt}) - len(set(trail)))
    dN_erase.append(len(set(trail) - {x_tip}) - len(set(trail)))
    return dict(transport=dN_transport, write=dN_write, erase=dN_erase)


def detailed_balance(u, b):
    """Detailed balance COMPUTED from the model's own two-state kernel (the per-attempt
       form of project_model's corrected activation convention): K2 = [[1-u, u],
       [v, 1-v]], stationary pi solved exactly.  THE SEALED IDENTITY (C-91 world tier,
       LANE_T48_C_WORLD G09): pi K2 == pi and pi_meta/pi_stable == v/u == b exactly.
       Returns dict(stationary, ratio_is_b)."""
    u, b = Fr(u), Fr(b)
    v = u * b
    pi = (v / (u + v), u / (u + v))
    piK = (pi[0] * (1 - u) + pi[1] * v, pi[0] * u + pi[1] * (1 - v))
    return dict(stationary=piK == pi, ratio_is_b=(pi[0] / pi[1] == b and v / u == b))


# ---------------------------------------------------------------- the three ensembles
def ensemble_transport(nbr, a):
    """E1 TRANSPORT (iv' read literally): one writer step = ONE energy-conserving
       dilation event moving the written unit across one grain boundary; erase-behind
       releases exactly what write-ahead costs (dN == 0 in every direction, the sealed
       G07 zero).  Per-link amplitude = the single event's activation factor a (both
       saddle readings a = u and a = u*b are honest; the sealed result is independent of
       the choice, gated in checks_writing.py).
       THE SEALED IDENTITY (C-91 world tier, LANE_T48_C_WORLD V1): conserving AND
       critical -- mu = 1/deg = mu_c -- at EVERY dE and EVERY barrier."""
    deg = len(nbr[0])
    a = Fr(a)
    return kernel_pos(nbr, 1 - deg * a, a)


def ensemble_trail_retreat(nbr, u, b):
    """E2 TRAIL WITH RETREAT (the raw two rates, strings counted WITH backtracking):
       per attempt the tip writes a fresh neighbor (v = u*b each, 5 directions) or the
       tip grain erases and the string retreats (u, the backtracking link), or nothing.
       THE SEALED IDENTITY (C-91 world tier, LANE_T48_C_WORLD V2): conserving at every
       dE; uniform (m_back == m_fresh == 1/deg) EXACTLY at dE = 0 (b = 1) and split
       (b, 1)/(5b + 1) beside; zero stationary spatial drift beside nonzero per-state
       persistence at dE != 0."""
    u, v = Fr(u), Fr(u) * Fr(b)
    return kernel_edge(nbr, 1 - 5 * v - u, v, u)


def ensemble_trail_decay(nbr, u, b, counting="H1"):
    """E3 TRAIL WITH DECAY: the model's OWN erase channel as absorption -- extend (v
       each) or the tip decays and the string leaves the ensemble (probability u, the
       model's own erase probability).  counting='H1' is the walk ensemble of the sealed
       T-44 coupling (extend to all deg neighbors); counting='NB' is the
       non-backtracking variant (its criticality reference is the venue's own
       directed-edge operator, earned by the same row-sum instrument).
       THE SEALED IDENTITY (C-91 world tier, LANE_T48_C_WORLD V3): NEVER critical --
       mu = b/(deg*b + 1) (H1) or b/(5b + 1) (NB), f0 and E_b dropping out exactly; the
       mass gap in closed form ln(mu_c/mu) = ln(1 + e^{dE/kT}/l)."""
    u, v = Fr(u), Fr(u) * Fr(b)
    if counting == "H1":
        deg = len(nbr[0])
        return kernel_pos(nbr, 1 - deg * v - u, v)
    return kernel_edge(nbr, 1 - 5 * v - u, v, None)


# ---------------------------------------------------------------- verdict instruments
def transport_verdict(nbr, W):
    """The computed verdicts of a position-space writer kernel: conservation (every row
       sum exactly 1), earned uniformity, the induced mu, criticality (mu == 1/deg with
       deg COMPUTED from the venue), the entrywise operator identity, and the massless
       signature (I - W) 1 == 0.  Returns dict; every entry a computed value."""
    deg = len(nbr[0]) if all(len(r) == len(nbr[0]) for r in nbr) else None
    mu_c_ref = Fr(1, deg) if deg else None
    rs = set(srow_sums(W))
    conserving = rs == {F1}
    cs, ms = extract_pos(W, nbr)
    uniform = len(cs) == 1 and len(ms) == 1
    mu = next(iter(ms)) if uniform else None
    c = next(iter(cs)) if uniform else None
    return dict(deg=deg, mu_c_ref=mu_c_ref, row_sums=rs, conserving=conserving,
                uniform=uniform, mu=mu,
                at_criticality=bool(uniform and mu == mu_c_ref),
                op_identity=bool(uniform and op_identity_pos(W, nbr, c, mu)),
                massless_signature=all(1 - sum(r.values()) == 0 for r in W))


def retreat_verdict(nbr, W):
    """The computed verdicts of the directed-edge trail kernel: conservation, double
       stochasticity, the fresh/back amplitudes, uniformity (m_back == m_fresh), and
       whether the per-state redistribution totals are exactly 1.  Returns dict."""
    deg = len(nbr[0])
    rs = set(srow_sums(W))
    cols = set(scol_sums(W, len(nbr) * 6))
    cs, mf, mb, tots = extract_edge(W, nbr)
    m_fresh = next(iter(mf)) if len(mf) == 1 else None
    m_back = next(iter(mb)) if len(mb) == 1 else None
    return dict(deg=deg, row_sums=rs, conserving=rs == {F1},
                doubly_stochastic=rs == {F1} and cols == {F1},
                m_fresh=m_fresh, m_back=m_back,
                uniform=bool(m_fresh is not None and m_fresh == m_back),
                redistribution_exact=tots == {F1},
                mu_c_ref=Fr(1, deg))


def decay_verdict(nbr, W):
    """The computed verdicts of the H1-counting decay kernel: the constant row sum
       (1 - loss), the loss, the induced mu, and THE COMPUTED MASS RATIO mu_c/mu formed
       from the kernel's own entries (the closed form 1 + e^{dE/kT}/deg is CHECKED
       against this in checks_writing.py, never sourced from it).  Returns dict."""
    deg = len(nbr[0])
    mu_c_ref = Fr(1, deg)
    rs = set(srow_sums(W))
    beta = next(iter(rs)) if len(rs) == 1 else None
    cs, ms = extract_pos(W, nbr)
    uniform = len(cs) == 1 and len(ms) == 1
    mu = next(iter(ms)) if uniform else None
    return dict(deg=deg, mu_c_ref=mu_c_ref, row_sum=beta,
                loss=(1 - beta) if beta is not None else None,
                conserving=rs == {F1}, uniform=uniform, mu=mu,
                mass_ratio=(mu_c_ref / mu) if mu else None,
                below_criticality=bool(mu is not None and mu < mu_c_ref))


def decay_verdict_nb(nbr, W):
    """The computed verdicts of the NB-counting decay kernel; the criticality reference
       deg_NB is EARNED from the venue's own directed-edge operator row sums (Hashimoto
       1989 as reference only), never declared.  Returns dict."""
    B = nb_edge_adjacency(nbr)
    rsB = set(srow_sums(B))
    deg_nb = next(iter(rsB)) if len(rsB) == 1 else None
    mu_c_nb = Fr(1, deg_nb) if deg_nb else None
    rs = set(srow_sums(W))
    beta = next(iter(rs)) if len(rs) == 1 else None
    cs, mf, mb, tots = extract_edge(W, nbr)
    mu = next(iter(mf)) if len(mf) == 1 else None
    return dict(deg_nb=deg_nb, mu_c_nb=mu_c_nb, row_sum=beta,
                loss=(1 - beta) if beta is not None else None,
                back_channel_empty=mb == set(), mu=mu,
                redistribution_total=next(iter(tots)) if len(tots) == 1 else None,
                mass_ratio=(mu_c_nb / mu) if mu else None,
                below_criticality=bool(mu is not None and mu_c_nb is not None
                                       and mu < mu_c_nb))


def closed_form_gap_ratio(b, l):
    """THE CLOSED FORM of the C-91 mass ratio, mu_c/mu = 1 + e^{dE/kT}/l with
       e^{dE/kT} = 1/b, as an exact rational.  This is the CHECK side: checks_writing.py
       gates it EQUAL to decay_verdict's computed mass_ratio -- the closed form is never
       the source of the computed number (D-8)."""
    return 1 + F1 / (Fr(l) * Fr(b))


# =====================================================================================
# 6  THE OBSERVATION GATE (D-25)  --  how a NEW surface's number enters the model
# =====================================================================================
def surface_boltzmann(s):
    """The Boltzmann dials of a REAL record surface: b = exp(-dE/kT) and
       u = exp(-E_b/kT), floats (the world entry; everything downstream that decides is
       exact).  D-25 ENFORCED: REFUSES a surface that does not carry provenance (build
       it through URM.surface).  Declines (None) where the surface is not thermally
       activated, exactly as the model's laws do."""
    if not getattr(s, "provenance", None):
        raise ValueError(
            "writing tier REFUSES: a world-tier surface must enter through the D-25 "
            "provenance gate (URM.surface) -- the model is grounded in real record "
            "data, never the toy category (the principal, 2026-08-20).")
    if not s.thermal:
        return None
    import grounded as G
    from math import exp
    kT = G.KB * s.T
    return dict(b=exp(-s.dE / kT), u=exp(-s.E_b / kT), dE_over_kT=s.dE / kT)


def rational_bracket(x, den):
    """The declared-denominator rational bracket lo <= x <= hi, lo and hi exact
       Fractions with hi - lo == 1/den."""
    from math import floor
    lo = Fr(floor(x * den), den)
    return lo, lo + Fr(1, den)


def surface_gap(s, n=4, den=10 ** 9, u_samples=(Fr(1, 20), Fr(1, 100))):
    """A NEW OBSERVATION ENTERS: the written-trail mass gap of a real record surface
       (C-91 world tier), computed -- not read off the closed form.  The surface's
       b = exp(-dE/kT) is bracketed by exact rationals (declared denominator), the E3
       decay kernel is BUILT AND MEASURED exactly at both brackets on the declared n^3
       venue at every declared u sample (E_b/f0 independence is thereby re-computed for
       THIS surface's entry, not remembered), and the closed form is CHECKED against
       every computed ratio.  The float gap ln(mu_c/mu) is then certified INSIDE the
       computed bracket.  D-25: refuses without provenance (surface_boltzmann).
       Returns dict(b, bracket, computed_ratios, closed_form_agrees, u_independent,
       gap_ln, gap_bracket_ln, contained) or None where the surface declines."""
    from math import log
    dial = surface_boltzmann(s)
    if dial is None:
        return None
    _cells, _idx, nbr = torus3(n)
    lo, hi = rational_bracket(dial["b"], den)
    ratios = {}
    agree = True
    for b_r in (lo, hi):
        per_u = set()
        for u_r in u_samples:
            W = ensemble_trail_decay(nbr, u_r, b_r, "H1")
            v = decay_verdict(nbr, W)
            per_u.add(v["mass_ratio"])
            agree &= (v["mass_ratio"] == closed_form_gap_ratio(b_r, v["deg"]))
        ratios[b_r] = per_u
    u_independent = all(len(per_u) == 1 for per_u in ratios.values())
    r_lo = next(iter(ratios[lo]))   # ratio is decreasing in b: r_lo >= true >= r_hi
    r_hi = next(iter(ratios[hi]))
    gap_ln = log(1 + (1 / dial["b"]) / len(nbr[0])) if dial["b"] > 0 else None
    bracket_ln = (log(r_hi), log(r_lo))
    return dict(b=dial["b"], dE_over_kT=dial["dE_over_kT"], bracket=(lo, hi),
                computed_ratios={str(k): sorted(str(x) for x in v)
                                 for k, v in ratios.items()},
                closed_form_agrees=agree, u_independent=u_independent,
                gap_ln=gap_ln, gap_bracket_ln=bracket_ln,
                contained=bool(gap_ln is not None
                               and bracket_ln[0] <= gap_ln <= bracket_ln[1]))
