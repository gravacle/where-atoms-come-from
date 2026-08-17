"""
REFUTATION LANE for W19-B.  LENS: MATHEMATICS.
This library IMPORTS the target lane's own instrument (lib_b.py) unchanged, so that every number
below is computed by THE SAME code the target used.  Nothing here re-implements their entropies.
What is added here is only what is needed to check their CRITERION STATEMENTS against the actual
definitions:

  1. MixE            -- Pauli expectations for a MIXED state (their PauliExpect is pure-state only).
                        Needed because the whole "threshold is tight" claim rests on global purity.
  2. all_subspaces   -- COMPLETE enumeration of Pauli *-subalgebras on a region, over GF(2).
                        Needed to ask whether their five-item choice list is complete or convenient.
  3. rdelta_zurek    -- the PUBLISHED redundancy R_delta = N/f_delta (Blume-Kohout & Zurek), which
                        is NOT the count-over-a-fixed-partition the target lane used.
  4. sbs_defects_in_basis -- the target's SBS test with the pointer basis of S made a FREE VARIABLE,
                        because Korbicz's SBS is an EXISTENTIAL over bases and their code fixes Z.

CEILING DECLARED UP FRONT: dense numpy, L <= 12 for constructed states (dim 4096), L <= 10 for
dense diagonalisation, and algebra_entropy costs 2^c * 4^k so k <= 6 in practice.  Every place I
hit one of these is named in the output.
"""
import sys, os, itertools
import numpy as np

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_B_CRITERION")
from lib_b import *          # noqa: F401,F403  -- THE TARGET'S OWN INSTRUMENT, UNMODIFIED
import lib_b


# ------------------------------------------------------------------ 1. mixed states
class MixE:
    """<P>_rho for rho = sum_j w_j |psi_j><psi_j|.  Duck-types lib_b.PauliExpect: algebra_entropy
       and mutual_information call E(x, z) and nothing else, so this drops straight in."""
    def __init__(self, psis, weights, L, mode=None):
        self.parts = [pauli_table(p, L, mode) for p in psis]
        self.w = list(weights)
        self.L = L
    def __call__(self, x, z):
        return sum(w * E(x, z) for w, E in zip(self.w, self.parts))


def rho_from_mixture(psis, weights):
    D = len(psis[0])
    r = np.zeros((D, D), dtype=complex)
    for w, p in zip(weights, psis):
        r += w * np.outer(p, p.conj())
    return r


def reduce_rho(rho, L, keep):
    """Partial trace of a full density matrix onto `keep`.  Same bit convention as lib_b:
       link l sits in bit l of the state index, hence numpy axis (L-1-l)."""
    keep = sorted(keep)
    other = [i for i in range(L) if i not in keep]
    pk = [L - 1 - l for l in keep] + [L - 1 - l for l in other]
    perm = pk + [L + p for p in pk]
    dk = 1 << len(keep); do = 1 << len(other)
    t = rho.reshape([2] * (2 * L)).transpose(perm).reshape(dk, do, dk, do)
    return np.einsum("iaja->ij", t)


def mi_ext_rho(rho, L, A, B):
    return (vn_entropy(reduce_rho(rho, L, A)) + vn_entropy(reduce_rho(rho, L, B))
            - vn_entropy(reduce_rho(rho, L, sorted(list(A) + list(B)))))


# ------------------------------------------------------------------ 2. complete subalgebra lattice
def span(vecs, n):
    """All elements of the GF(2) span of `vecs` (packed ints, n bits)."""
    s = {0}
    for v in vecs:
        s |= {u ^ v for u in s}
    return frozenset(s)


def all_subspaces(n):
    """EVERY GF(2)-subspace of F_2^n, as frozensets of packed ints.  COMPLETE, NOT SAMPLED:
       every subspace has exactly one reduced row echelon basis, so enumerating (pivot set, free
       entries) enumerates the subspaces exactly once each.  The count is checked against the
       Galois number G_n in r0 for n = 1..6.
       (A first version enumerated all subsets of nonzero vectors and de-duplicated by span; that
       is also exhaustive but costs C(2^n-1, n) and was unusable at n = 6.  Same output, verified
       by the Galois-number check and by the fact that the returned list is canonically sorted.)"""
    out = set()
    for k in range(n + 1):
        for piv in itertools.combinations(range(n), k):
            free = [(i, c) for i, p in enumerate(piv) for c in range(p + 1, n) if c not in piv]
            for bits in itertools.product((0, 1), repeat=len(free)):
                rows = [1 << p for p in piv]
                for (i, c), b in zip(free, bits):
                    if b: rows[i] |= (1 << c)
                out.add(span(rows, n))
    return sorted(out, key=lambda s: (len(s), sorted(s)))


GALOIS = {0: 1, 1: 2, 2: 5, 3: 16, 4: 67, 5: 374, 6: 2825}   # number of subspaces of F_2^n


def vec_to_sp(v, R, n):
    """v in F_2^{2n}: low n bits = x-part on R (in sorted order), high n bits = z-part."""
    pos = sorted(R)
    x = z = 0
    for i, l in enumerate(pos):
        if (v >> i) & 1:        x |= (1 << l)
        if (v >> (n + i)) & 1:  z |= (1 << l)
    return hermitize(SP(x, z, 1.0))


def gauge_invariant_mask(car, R):
    """Which v in F_2^{2n} are GLOBALLY gauge-invariant Paulis supported in R.
       Condition (lib_b's own): z . G_v = 0 mod 2 for every Gauss string, x unconstrained."""
    n = len(R)
    ok = []
    for v in range(1 << (2 * n)):
        sp = vec_to_sp(v, R, n)
        if all(not anticomm(sp, SP(g, 0)) for g in car["gauss"]):
            ok.append(v)
    return ok


def center_kind(car, gens):
    """Classify the CENTRE of the algebra: 'trivial', 'electric' (every central generator is a
       pure-X string), 'magnetic' (every central generator carries a Z), or 'mixed'."""
    if not gens:
        return "trivial(unit)", 0
    pairs, cen = algebra_structure(gens)
    if not cen:
        return "trivial(factor)", 0
    haveZ = [c.z != 0 for c in cen]
    if not any(haveZ):  return "ELECTRIC", len(cen)
    if all(haveZ):      return "MAGNETIC", len(cen)
    return "mixed", len(cen)


def support_links(gens, L):
    s = set()
    for g in gens:
        for l in range(L):
            if ((g.x >> l) & 1) or ((g.z >> l) & 1):
                s.add(l)
    return sorted(s)


def is_gauge_invariant(car, gens):
    return all(all(not anticomm(g, SP(gv, 0)) for gv in car["gauss"]) for g in gens)


# ------------------------------------------------------------------ 3. the two R_delta's
def rdelta_count(I_singles, HS, delta=0.1, tol=1e-9):
    """THE TARGET LANE'S R_delta: count the members of ONE FIXED partition that clear the bar.
       (b1_discrimination.py line 'Rdelta = int(sum(1 for v in single if ...))')"""
    return int(sum(1 for v in I_singles if v >= (1 - delta) * HS - tol))


def rdelta_zurek(Ibysize, N, HS, delta=0.1, tol=1e-9):
    """THE PUBLISHED R_delta (Blume-Kohout & Zurek 2006; Zurek 2009):  R_delta = N / m*, where m*
       is the SMALLEST fragment size whose AVERAGE mutual information reaches (1-delta) H(S).
       Ibysize[m] = mean over all fragments of size m of I(S:F).  Returns (R, m*)."""
    for m in range(1, N + 1):
        if Ibysize[m] >= (1 - delta) * HS - tol:
            return N / m, m
    return 0.0, None


def rdelta_packing(qual_masks, N):
    """THE READING THE BRIEF ACTUALLY STATES, computed EXACTLY: 'partition E into DISJOINT
       fragments ... R_delta counts how many disjoint fragments independently reach it'.
       Given the set of fragment bitmasks that clear the bar, return the MAXIMUM number of
       PAIRWISE DISJOINT ones.  This is a set-packing problem; solved exactly by DP over the
       2^N availability masks (N <= 10 here, so this is cheap and not a heuristic).
       NOTE this fixes an over-report in the naive 'floor(N/m*)' proxy: when only ONE fragment
       of the winning size qualifies (the Bell-localised arm), floor(N/m*) says N and the exact
       packing says 1."""
    qual = sorted(qual_masks)
    bybit = {}
    for q in qual:
        low = (q & -q).bit_length() - 1
        bybit.setdefault(low, []).append(q)
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def f(avail):
        if avail == 0:
            return 0
        low = (avail & -avail).bit_length() - 1
        best = f(avail & ~(1 << low))          # leave this link unused
        for q in bybit.get(low, ()):
            if q & ~avail == 0:
                best = max(best, 1 + f(avail & ~q))
        return best
    r = f((1 << N) - 1)
    f.cache_clear()
    return r


# ------------------------------------------------------------------ 4. SBS with the basis freed
def rot_link(psi, L, l, V):
    """Apply the 2x2 unitary V to link l of a state vector (bit l of the index)."""
    t = psi.reshape([2] * L)
    t = np.moveaxis(t, L - 1 - l, 0).reshape(2, -1)
    t = V @ t
    return np.moveaxis(t.reshape([2] * L), 0, L - 1 - l).reshape(-1)


def bloch_basis(theta, phi):
    """The unitary whose columns are the Bloch-sphere basis {|n>, |-n>}; V^dagger rotates that
       basis onto the computational one."""
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -np.exp(-1j * phi) * s],
                     [np.exp(1j * phi) * s, c]], dtype=complex)


def sbs_defects(psi, L, S, frags, extra=(), tol=1e-12):
    """EXACTLY the target's sbs_report arithmetic (b2_sbs_vs_rdelta.py), re-expressed so the
       pointer basis is whatever basis the CALLER has already rotated into.  D_coh dephases S in
       the computational basis; D_prod is strong independence; D_orth is distinguishability."""
    order = list(S) + sum([list(f) for f in frags], []) + list(extra)
    rho = rdm_ordered(psi, L, order)
    dims = [2] + [1 << len(f) for f in frags] + ([1 << len(extra)] if extra else [])
    d = int(np.prod(dims)); step = d // 2
    dephased = np.zeros_like(rho)
    for i in range(2):
        dephased[i * step:(i + 1) * step, i * step:(i + 1) * step] = \
            rho[i * step:(i + 1) * step, i * step:(i + 1) * step]
    D_coh = 0.5 * trace_norm(rho - dephased)
    ps, sig = [], []
    for i in range(2):
        blk = rho[i * step:(i + 1) * step, i * step:(i + 1) * step]
        p = float(np.trace(blk).real); ps.append(p)
        sig.append(blk / p if p > tol else blk)
    sub = dims[1:]
    D_prod = 0.0; marg = []
    for i in range(2):
        m = [ptrace_blocks(sig[i], sub, [k]) for k in range(len(sub))]
        marg.append(m)
        prod = np.array([[1.0 + 0j]])
        for mm in m: prod = np.kron(prod, mm)
        D_prod = max(D_prod, 0.5 * trace_norm(sig[i] - prod))
    D_orth = 0.0
    for k in range(len(sub)):
        D_orth = max(D_orth, trace_norm(msqrt(marg[0][k]) @ msqrt(marg[1][k])))
    return dict(p=ps, D_coh=D_coh, D_prod=D_prod, D_orth=D_orth,
                worst=max(D_coh, D_prod, D_orth))


def rot_block(psi, L, links, V):
    """Apply the 2^|links| x 2^|links| unitary V to the tensor factor formed by `links`, in the
       order given (links[0] most significant)."""
    n = len(links)
    other = [i for i in range(L) if i not in links]
    perm = [L - 1 - l for l in links] + [L - 1 - l for l in other]
    t = psi.reshape([2] * L).transpose(perm).reshape(1 << n, -1)
    t = V @ t
    t = t.reshape([2] * L)
    inv = np.argsort(perm)
    return t.transpose(inv).reshape(-1)


def sbs_defects_multi(psi, L, S, frags, extra=(), tol=1e-12):
    """The same three SBS defects for a system of ANY dimension 2^|S|, dephasing in whatever basis
       the caller has already rotated the system block into.  (b2's sbs_report hard-codes a
       2-dimensional system AND the Z basis; both restrictions are lifted here.)"""
    order = list(S) + sum([list(f) for f in frags], []) + list(extra)
    rho = rdm_ordered(psi, L, order)
    nS = 1 << len(S)
    dims = [nS] + [1 << len(f) for f in frags] + ([1 << len(extra)] if extra else [])
    d = int(np.prod(dims)); step = d // nS
    dephased = np.zeros_like(rho)
    for i in range(nS):
        dephased[i * step:(i + 1) * step, i * step:(i + 1) * step] = \
            rho[i * step:(i + 1) * step, i * step:(i + 1) * step]
    D_coh = 0.5 * trace_norm(rho - dephased)
    ps, sig = [], []
    for i in range(nS):
        blk = rho[i * step:(i + 1) * step, i * step:(i + 1) * step]
        p = float(np.trace(blk).real); ps.append(p)
        sig.append(blk / p if p > tol else blk)
    sub = dims[1:]
    D_prod = 0.0; marg = []
    for i in range(nS):
        m = [ptrace_blocks(sig[i], sub, [k]) for k in range(len(sub))]
        marg.append(m)
        prod = np.array([[1.0 + 0j]])
        for mm in m: prod = np.kron(prod, mm)
        if ps[i] > tol:
            D_prod = max(D_prod, 0.5 * trace_norm(sig[i] - prod))
    D_orth = 0.0
    for k in range(len(sub)):
        for i in range(nS):
            for j in range(i + 1, nS):
                if ps[i] > tol and ps[j] > tol:
                    D_orth = max(D_orth, trace_norm(msqrt(marg[i][k]) @ msqrt(marg[j][k])))
    return dict(p=ps, D_coh=D_coh, D_prod=D_prod, D_orth=D_orth,
                worst=max(D_coh, D_prod, D_orth))


def sbs_best_over_bases(psi, L, S, frags, extra=(), ngrid=25, refine=3):
    """KORBICZ'S SBS IS AN EXISTENTIAL OVER POINTER BASES.  Minimise max(D_coh, D_prod, D_orth)
       over the orthonormal bases of the single-qubit system S.  Returns (best_defects, angles).
       CONFOUND RECORDED: this is a GRID plus local refinement, not a proof of the minimum.
       ngrid must be ODD so that theta = pi/2 -- the equatorial (X/Y) bases -- is actually sampled;
       a first version used ngrid = 24, missed pi/2 entirely, and reported 1.000 for a state whose
       X basis gives 2e-17.  That bug is exactly the failure mode this file is about, so it is
       recorded here rather than silently fixed."""
    assert len(S) == 1
    assert ngrid % 2 == 1, "ngrid must be odd so theta = pi/2 is sampled"
    best = None
    grid = [(t, p) for t in np.linspace(0, np.pi, ngrid)
            for p in np.linspace(0, 2 * np.pi, ngrid - 1, endpoint=False)]
    for (t, p) in grid:
        V = bloch_basis(t, p)
        d = sbs_defects(rot_link(psi, L, S[0], V.conj().T), L, S, frags, extra)
        if best is None or d["worst"] < best[0]["worst"]:
            best = (d, (t, p))
    t0, p0 = best[1]; w = np.pi / ngrid
    for _ in range(refine):
        for t in np.linspace(max(0, t0 - w), min(np.pi, t0 + w), 9):
            for p in np.linspace(p0 - w, p0 + w, 9):
                V = bloch_basis(t, p)
                d = sbs_defects(rot_link(psi, L, S[0], V.conj().T), L, S, frags, extra)
                if d["worst"] < best[0]["worst"]:
                    best = (d, (t, p))
        t0, p0 = best[1]; w /= 4
    return best


# ------------------------------------------------------------------ misc
def elec_ghz(car, a=None):
    L = car["L"]; popc = popc_table(1 << L)
    plus = np.ones(1 << L, dtype=complex)
    minus = ((-1.0) ** popc[np.arange(1 << L)]).astype(complex)
    plus /= np.linalg.norm(plus); minus /= np.linalg.norm(minus)
    if a is None: a = 1 / np.sqrt(2)
    psi = a * plus + np.sqrt(1 - a * a) * minus
    return project_physical(psi / np.linalg.norm(psi), car)


def hr(title, ch="="):
    print("\n" + ch * 104)
    print(title)
    print(ch * 104)
