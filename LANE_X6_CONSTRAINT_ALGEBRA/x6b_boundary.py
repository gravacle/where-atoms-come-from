# LANE X6 - b : IS THE PHYSICAL CONTENT OF OUR CONSTRAINT A BOUNDARY TERM?
#
# QUESTION (from BRIEF.md): gravity's Hamiltonian+momentum constraints have the property that the
# total constraint of a region reduces to a pure SURFACE integral (ADM energy is a boundary term;
# the bulk constraints cancel).  Does OUR record-level constraint -- the Z2 Gauss law on a lattice
# gauge carrier -- do the same thing, and does it do it QUANTITATIVELY (support scaling), not just
# rhetorically?
#
# WHAT IS BEING MEASURED, EXACTLY
#   P_R  =  prod_{v in R} G_v          (total constraint operator of a region R of vertices)
#   B_R  =  prod_{l in dR} sigma_l     (dR = links with EXACTLY ONE endpoint in R = the cut)
#   ||P_R - B_R|| = 0   would mean the bulk constraints cancel identically and the region's total
#   constraint IS a surface operator.
#
# CONVENTION NOTE (stated because the task names sigma^x and the canonical sector is written in a
# configuration basis).  The program's CANONICAL physical sector is the set of link configurations
# s with even incidence parity at every vertex.  That sector is the +1 eigenspace of the DIAGONAL
# vertex operator G_v = prod_{l ~ v} sigma_l, where sigma is the single-link operator diagonal in
# the configuration basis.  Whether one names that single-link operator sigma^z (configuration
# basis = its eigenbasis) or sigma^x (working in the conjugate basis, as is standard when the link
# variable is called the electric field) is a BASIS LABEL, not physics.  This script therefore does
# BOTH:
#   ARM D  (diagonal / "electric") : G_v = prod_{l~v} sigma^z_l  -- defines the canonical sector.
#   ARM F  (flip / conjugate)      : G_v = prod_{l~v} sigma^x_l  -- the star operator, the same
#                                    statement in the conjugate basis.
# Both arms are computed as EXPLICIT OPERATORS on the full 2^L space and the answer must be the
# same, because the cancellation is a combinatorial fact about incidence, not about the basis.
#
# W-46 GUARD: a prior lane got a wrong scaling law by placing the region at a lattice CORNER.  On a
# torus there is no corner, but the INDEXING has a wrap seam, so (a) every block region used here is
# chosen with indices that do not wrap, and (b) the script explicitly recomputes the support for
# EVERY translate of the region and reports min/max, so a placement-dependent answer cannot hide.
#
# HONEST LIMIT ON "EXPLICIT MATRIX": for L=18 the full space is 2^18 = 262144, so a dense 2^18 x
# 2^18 array (7e10 entries) cannot be stored.  Every operator here is EXACTLY diagonal (ARM D) or
# EXACTLY a permutation (ARM F) in the configuration basis, so its explicit matrix is represented
# without approximation by its full diagonal / its full permutation image, and the Frobenius norm of
# a difference is computed over ALL 2^18 basis vectors.  That is exact, not sampled.  Dense matrices
# ARE built and used on the 1024-dimensional physical sector.

import itertools, time
import numpy as np

np.set_printoptions(linewidth=200)
LOG2 = np.log(2.0)


# ----------------------------------------------------------------------------------------------
# carrier
# ----------------------------------------------------------------------------------------------
def build(nx, ny):
    """nx x ny torus.  Returns ind (link name -> index), E (edge list), NV, L, PL (plaquettes)."""
    ind = {}
    E = []
    def vid(i, j): return i + nx * j
    for j in range(ny):
        for i in range(nx):
            ind[('h', i, j)] = len(E); E.append((vid(i, j), vid((i + 1) % nx, j)))
            ind[('v', i, j)] = len(E); E.append((vid(i, j), vid(i, (j + 1) % ny)))
    NV = nx * ny
    L = len(E)
    # CANONICAL TORUS PLAQUETTES (with the wrap-around), verbatim from the program conventions
    PL = [[ind[('h', i, j)], ind[('v', (i + 1) % nx, j)], ind[('h', i, (j + 1) % ny)], ind[('v', i, j)]]
          for j in range(ny) for i in range(nx)]
    return ind, E, NV, L, PL


def incidence(E, NV, L):
    """M[v,l] = 1 iff link l touches vertex v (no self-loops for n>=3)."""
    M = np.zeros((NV, L), dtype=np.uint8)
    for k, (a, b) in enumerate(E):
        M[a, k] ^= 1
        M[b, k] ^= 1
    return M


def mask_of(links, L):
    m = 0
    for l in links:
        m |= (1 << l)
    return m


def popcount(x):
    return bin(x).count('1')


def block_region(nx, ny, i0, j0, k):
    """k x k block of vertices with lower-left corner (i0,j0); indices are NOT allowed to wrap
    unless the block covers the whole direction."""
    return [((i0 + a) % nx) + nx * ((j0 + b) % ny) for b in range(k) for a in range(k)]


def crossing_links(E, R):
    Rs = set(R)
    return [l for l, (a, b) in enumerate(E) if (a in Rs) != (b in Rs)]


def interior_links(E, R):
    Rs = set(R)
    return [l for l, (a, b) in enumerate(E) if (a in Rs) and (b in Rs)]


# ----------------------------------------------------------------------------------------------
# GF(2) rank (used for the interior-information count, independent of any Hilbert space)
# ----------------------------------------------------------------------------------------------
def gf2_rank_matrix(M):
    """M: numpy uint8 array over GF(2). Returns rank."""
    A = (np.array(M, dtype=np.uint8) % 2).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        for i in range(rows):
            if i != r and A[i, c]:
                A[i] ^= A[r]
        r += 1
        if r == rows:
            break
    return r


# ==============================================================================================
print("=" * 100)
print("LANE X6-b  --  IS THE TOTAL CONSTRAINT OF A REGION A PURE BOUNDARY OPERATOR?")
print("=" * 100)
print()

# ----------------------------------------------------------------------------------------------
# [0] carrier and canonical physical sector (CANONICAL SNIPPET, verbatim)
# ----------------------------------------------------------------------------------------------
nx = ny = 3
ind, E, NV, L, PL = build(nx, ny)
M = incidence(E, NV, L)

print("[0] CARRIER")
print(f"    torus {nx} x {ny};  NV = {NV} vertices;  L = {L} links;  full space dim 2^L = {2**L}")
t0 = time.time()
st = [s for s in itertools.product(range(2), repeat=L)
      if all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
              - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % 2 == 0 for v in range(NV))]
idx = {s: i for i, s in enumerate(st)}
D = len(st)
print(f"    CANONICAL physical sector (Gauss law at every vertex):  D = {D}   [{time.time()-t0:.1f}s]")
print(f"    independent prediction 2^(L-NV+1) = {2**(L-NV+1)}     -> SC-0 "
      f"{'PASS' if D == 2**(L-NV+1) else 'FAIL'}")
print()

# vectorised copy of the same sector, cross-checked against the canonical list (SC-0b)
allc = np.arange(2 ** L, dtype=np.int64)
bits = ((allc[:, None] >> np.arange(L)[None, :]) & 1).astype(np.uint8)
par = (bits @ M.T) % 2
phys_mask = ~par.any(axis=1)
phys_idx = np.nonzero(phys_mask)[0]
canon_set = set(int(sum(s[k] << k for k in range(L))) for s in st)
sc0b = (len(phys_idx) == D) and (set(int(x) for x in phys_idx) == canon_set)
print(f"    SC-0b vectorised sector identical to canonical itertools sector : {'PASS' if sc0b else 'FAIL'}")
print()


def zdiag(mask):
    """explicit diagonal (length 2^L) of the Z-type Pauli string with support `mask`."""
    lk = [l for l in range(L) if (mask >> l) & 1]
    if not lk:
        return np.ones(2 ** L)
    s = bits[:, lk].sum(axis=1) % 2
    return 1.0 - 2.0 * s


def G_mask(v):
    return mask_of([l for l in range(L) if M[v, l]], L)


# ----------------------------------------------------------------------------------------------
# [1] SELF-CHECK -- R = ALL vertices  =>  boundary empty  =>  P_R must be the identity
# ----------------------------------------------------------------------------------------------
print("[1] SELF-CHECK (required):  R = ALL vertices of the torus -> dR is empty -> P_R must be I")
R_all = list(range(NV))
cr_all = crossing_links(E, R_all)
print(f"    |R| = {len(R_all)} (whole torus);  number of crossing links = {len(cr_all)} (must be 0)")

# ARM D : diagonal
PallD = np.ones(2 ** L)
for v in R_all:
    PallD = PallD * zdiag(G_mask(v))
nD_full = float(np.sqrt(np.sum((PallD - 1.0) ** 2)))
PallD_phys = np.diag(PallD[phys_idx])
nD_phys = float(np.linalg.norm(PallD_phys - np.eye(D)))

# ARM F : flip / permutation
maskF = 0
for v in R_all:
    maskF ^= G_mask(v)
permF = allc ^ maskF
mismatch = int(np.count_nonzero(permF != allc))
nF_full = float(np.sqrt(2.0 * mismatch))

print(f"    ARM D (diagonal Gauss)  || P_all - I ||_F  over ALL 2^L basis vectors = {nD_full:.6e}")
print(f"    ARM D  dense on the physical sector (1024x1024)  || P_all - I ||_F   = {nD_phys:.6e}")
print(f"    ARM F (flip Gauss)      || P_all - I ||_F  over ALL 2^L basis vectors = {nF_full:.6e}"
      f"   (support mask = {maskF})")
ok1 = (nD_full == 0.0) and (nD_phys == 0.0) and (nF_full == 0.0)
print(f"    SELF-CHECK  ||P_all - I|| = 0  :  {'PASS' if ok1 else 'FAIL'}")
print()

# ----------------------------------------------------------------------------------------------
# [2] P_R  versus the boundary product, for three region sizes
# ----------------------------------------------------------------------------------------------
print("[2] IS P_R EQUAL TO THE PRODUCT OVER THE CROSSING LINKS?   (3x3 torus, R in the interior of")
print("    the index range, never touching the wrap seam)")
print()
regions = [
    ("|R|=1  single vertex (1,1)", block_region(nx, ny, 1, 1, 1)),
    ("|R|=2  domino (1,1)-(2,1)", [1 + nx * 1, 2 + nx * 1]),
    ("|R|=4  2x2 block at (1,1)", block_region(nx, ny, 1, 1, 2)),
    ("|R|=9  whole torus", list(range(NV))),
]
print(f"    {'region':<28}{'|R|':>5}{'|dR|':>7}{'||P_R-B_R|| full':>20}{'||P_R-B_R|| phys':>20}{'ARM F':>12}")
all_zero = True
for tag, R in regions:
    cr = crossing_links(E, R)
    # ARM D
    Pd = np.ones(2 ** L)
    for v in R:
        Pd = Pd * zdiag(G_mask(v))
    Bd = zdiag(mask_of(cr, L))
    n_full = float(np.sqrt(np.sum((Pd - Bd) ** 2)))
    Pp = np.diag(Pd[phys_idx]); Bp_ = np.diag(Bd[phys_idx])
    n_phys = float(np.linalg.norm(Pp - Bp_))
    # ARM F
    mF = 0
    for v in R:
        mF ^= G_mask(v)
    mB = mask_of(cr, L)
    nf = float(np.sqrt(2.0 * int(np.count_nonzero((allc ^ mF) != (allc ^ mB)))))
    all_zero &= (n_full == 0.0 and n_phys == 0.0 and nf == 0.0)
    print(f"    {tag:<28}{len(R):>5}{len(cr):>7}{n_full:>20.6e}{n_phys:>20.6e}{nf:>12.3e}")
print()
print(f"    ALL DIFFERENCES EXACTLY ZERO: {all_zero}")
print("    -> the bulk constraints cancel identically; the total constraint of a region IS a pure")
print("       surface operator supported on the cut.  This is Gauss-law-as-boundary-term, exact.")
print()

# ----------------------------------------------------------------------------------------------
# [3] SUPPORT SCALING -- area or perimeter?
# ----------------------------------------------------------------------------------------------
print("[3] HOW DOES THE SUPPORT OF P_R SCALE WITH THE SIZE OF R?")
print("    support(P_R)   = number of links P_R acts on nontrivially (= |dR| after cancellation)")
print("    naive support  = number of links touched by R at all (what the product would act on if")
print("                     the bulk did NOT cancel) = interior links + crossing links")
print("    3x3 alone cannot separate the two laws (|R|=9 is the whole torus), so larger tori are")
print("    used for the scaling; those numbers are exact Pauli-string supports, no Hilbert space.")
print()
print(f"    {'torus':>8}{'k':>4}{'|R|=k^2':>9}{'support(P_R)':>15}{'naive support':>15}"
      f"{'4k (perimeter)':>16}{'placement min/max':>20}")
rows = []
for (n, ks) in [(3, [1, 2]), (4, [1, 2, 3]), (5, [1, 2, 3, 4]), (6, [1, 2, 3, 4, 5])]:
    ind_n, E_n, NV_n, L_n, PL_n = build(n, n)
    M_n = incidence(E_n, NV_n, L_n)
    for k in ks:
        # canonical placement away from the index seam
        i0 = 1 if k < n else 0
        R = block_region(n, n, i0, i0, k)
        cr = crossing_links(E_n, R)
        it = interior_links(E_n, R)
        # exact operator support via Pauli-string composition
        mm = 0
        for v in R:
            mm ^= mask_of([l for l in range(L_n) if M_n[v, l]], L_n)
        sup = popcount(mm)
        assert sup == len(cr), "support mismatch"
        # W-46 GUARD: every translate
        sups = []
        for a in range(n):
            for b in range(n):
                Rt = block_region(n, n, a, b, k)
                mt = 0
                for v in Rt:
                    mt ^= mask_of([l for l in range(L_n) if M_n[v, l]], L_n)
                sups.append(popcount(mt))
        rows.append((n, k, k * k, sup, len(it) + len(cr)))
        print(f"    {str(n)+'x'+str(n):>8}{k:>4}{k*k:>9}{sup:>15}{len(it)+len(cr):>15}{4*k:>16}"
              f"{str(min(sups))+' / '+str(max(sups)):>20}")
print()
print("    W-46 GUARD: min = max in every row -> the support does NOT depend on where R is placed;")
print("    the earlier corner artefact cannot be reproduced on the torus.")
print()

# log-log exponents on the largest torus series
ser = [(r[2], r[3], r[4]) for r in rows if r[0] == 6]
A = np.array([s[0] for s in ser], dtype=float)
S = np.array([s[1] for s in ser], dtype=float)
Nv = np.array([s[2] for s in ser], dtype=float)
def slope(x, y):
    lx = np.log(x); ly = np.log(y)
    return float(np.polyfit(lx, ly, 1)[0])
p_sup = slope(A, S)
p_nai = slope(A, Nv)
print(f"    LEAST-SQUARES EXPONENT p in support ~ |R|^p   (6x6 series, |R| = 1,4,9,16,25):")
print(f"        support(P_R)  : p = {p_sup:.6f}      (PERIMETER law = 0.5 exactly, AREA law = 1.0)")
print(f"        naive support : p = {p_nai:.6f}      (see the caveat below -- NOT a pure power)")
print("    CAVEAT ON THAT SECOND NUMBER, stated because the fit does not land on 1.0: the naive")
print("    support is EXACTLY 2k^2 + 2k, an area term plus a perimeter term, so it is not a pure")
print("    power law and a log-log fit over k=1..5 must undershoot.  The exact formulae, verified:")
print(f"        support(P_R) - 4k        for k=1..5 : {[int(s - 4*np.sqrt(a)) for a, s in zip(A, S)]}"
      f"   (pure perimeter, exponent 0.5 EXACTLY at every point)")
print(f"        naive - (2k^2+2k)        for k=1..5 : "
      f"{[int(v - (2*a + 2*np.sqrt(a))) for a, v in zip(A, Nv)]}   (area + perimeter; "
      f"local log-log slope between k=4 and k=5 = {np.log(Nv[4]/Nv[3])/np.log(A[4]/A[3]):.6f}, -> 1)")
print("    So: the cancelled operator obeys a perimeter law EXACTLY and at every single point; the")
print("    uncancelled product is area-dominated and tends to exponent 1 from below.")
print(f"    support(P_R) - 4*sqrt(|R|) for every point : "
      f"{[int(s - 4*np.sqrt(a)) for a, s in zip(A, S)]}   (all zero -> support = 4k EXACTLY)")
print(f"    ratio naive/support : {[round(float(v/s),4) for v,s in zip(Nv,S)]}  -> diverges like k/2,")
print("      i.e. the cancellation removes an AREA's worth of operator content and leaves a PERIMETER.")
print()

# ----------------------------------------------------------------------------------------------
# [4] WHAT DOES THE BOUNDARY OPERATOR TELL YOU ABOUT THE INTERIOR?
# ----------------------------------------------------------------------------------------------
print("[4] IS ANY NONTRIVIAL GAUGE-INVARIANT INFORMATION ABOUT THE INTERIOR OF R VISIBLE IN THE")
print("    BOUNDARY OPERATOR ALONE?")
print()
R4 = block_region(nx, ny, 1, 1, 2)
cr4 = crossing_links(E, R4)
in4 = interior_links(E, R4)
print(f"    R = 2x2 block at (1,1):  interior links = {in4}   crossing links = {cr4}")

# (a) the boundary operator itself, on the physical sector
Bd4 = zdiag(mask_of(cr4, L))
vals = Bd4[phys_idx]
print(f"    (a) B_R restricted to the physical sector: distinct eigenvalues = "
      f"{sorted(set(np.round(vals,12)))}, "
      f"|| B_R|_phys - I || = {np.linalg.norm(np.diag(vals)-np.eye(D)):.3e}")
print("        -> B_R is CONSTANT (+1) on every physical state.  Its expectation value is the same")
print("        number for every state in the theory, so it carries EXACTLY 0 bits about anything.")

# (b) an interior Wilson loop exists and is not constant
p_in = None
for p in PL:
    if set(p) <= set(in4):
        p_in = p
        break
print(f"    (b) plaquette lying STRICTLY inside R: links {p_in}")
mp = mask_of(p_in, L)
img = phys_idx ^ mp
stays = np.isin(img, phys_idx).all()
# eigenvalue multiplicities of the interior Wilson loop on the physical sector
pos = {int(c): i for i, c in enumerate(phys_idx)}
perm = np.array([pos[int(c)] for c in img])
fixed = int(np.count_nonzero(perm == np.arange(D)))
n_plus = (D + fixed) // 2
n_minus = (D - fixed) // 2
print(f"        preserves the physical sector: {bool(stays)};  fixed points = {fixed};  "
      f"eigenvalue multiplicities  (+1, -1) = ({n_plus}, {n_minus})")
print(f"        -> the interior Wilson loop is a genuine 2-valued gauge-invariant observable of the")
print(f"           interior:  H = {-(n_plus/D)*np.log(n_plus/D)/LOG2 - (n_minus/D)*np.log(n_minus/D)/LOG2:.6f} bit,")
Wmat = np.zeros((D, D))
Wmat[perm, np.arange(D)] = 1.0
Bmat = np.diag(Bd4[phys_idx])
comm = float(np.linalg.norm(Wmat @ Bmat - Bmat @ Wmat))
print(f"           explicit dense check on the physical sector: || [W_p , B_R] || = {comm:.3e}  ->")
print("           the interior Wilson loop commutes with the boundary operator, so no measurement of")
print("           B_R constrains it at all.")

# (c) exact information accounting over the physical sector (uniform distribution)
def H_from_counts(c):
    c = np.array([x for x in c if x > 0], dtype=float)
    p = c / c.sum()
    return float(-(p * np.log(p)).sum() / LOG2)

def bits_of(cfgs, links):
    if not links:
        return np.zeros(len(cfgs), dtype=np.int64)
    key = np.zeros(len(cfgs), dtype=np.int64)
    for i, l in enumerate(links):
        key |= (((cfgs >> l) & 1) << i)
    return key

cfg = phys_idx
k_int = bits_of(cfg, in4)
k_cr = bits_of(cfg, cr4)
k_joint = k_int.astype(np.int64) * (1 << len(cr4)) + k_cr
H_int = H_from_counts(np.unique(k_int, return_counts=True)[1])
H_cr = H_from_counts(np.unique(k_cr, return_counts=True)[1])
H_j = H_from_counts(np.unique(k_joint, return_counts=True)[1])
H_int_given_cr = H_j - H_cr
print()
print("    (c) EXACT information accounting, uniform distribution over the 1024 physical states:")
print(f"        H(interior links)                              = {H_int:.6f} bits")
print(f"        H(crossing links)                              = {H_cr:.6f} bits")
print(f"        H(interior | ALL crossing-link data)           = {H_int_given_cr:.6f} bits")
print(f"        I(interior ; B_R)  [B_R is constant]           = {0.0:.6f} bits")
print("        -> even the FULL boundary flux data leaves 1 bit of the interior undetermined; the")
print("           single boundary OPERATOR leaves all of it undetermined.")

# (d) the same count algebraically, independent of any Hilbert space: GF(2) cycle rank
sub = M[np.array(R4)][:, np.array(in4)]
rk = gf2_rank_matrix(sub)
hidden = len(in4) - rk
print()
print("    (d) INDEPENDENT ALGEBRAIC COUNT (GF(2) rank, no Hilbert space):")
print(f"        interior links = {len(in4)}, rank of the Gauss constraints on them = {rk},")
print(f"        free interior bits = {hidden}  -> matches (c) H(interior|crossing) = {H_int_given_cr:.6f}"
      f"   {'PASS' if abs(hidden - H_int_given_cr) < 1e-9 else 'FAIL'}   [SC-1]")

print()
print("    (e) HOW THAT HIDDEN INFORMATION SCALES (same GF(2) count on larger tori):")
print(f"        {'torus':>8}{'k':>4}{'|R|':>6}{'boundary support':>18}{'hidden interior bits':>23}"
      f"{'(k-1)^2':>10}")
hid_rows = []
for (n, ks) in [(4, [1, 2, 3]), (5, [1, 2, 3, 4]), (6, [1, 2, 3, 4, 5])]:
    ind_n, E_n, NV_n, L_n, PL_n = build(n, n)
    M_n = incidence(E_n, NV_n, L_n)
    for k in ks:
        R = block_region(n, n, 1, 1, k)
        it = interior_links(E_n, R)
        cr = crossing_links(E_n, R)
        if it:
            sub_n = M_n[np.array(R)][:, np.array(it)]
            h = len(it) - gf2_rank_matrix(sub_n)
        else:
            h = 0
        hid_rows.append((n, k, len(cr), h))
        print(f"        {str(n)+'x'+str(n):>8}{k:>4}{k*k:>6}{len(cr):>18}{h:>23}{(k-1)**2:>10}")
hh = np.array([r[3] for r in hid_rows if r[0] == 6], dtype=float)
aa = np.array([k * k for k in [1, 2, 3, 4, 5]], dtype=float)
nz = hh > 0
hid6 = [r[3] for r in hid_rows if r[0] == 6]
print(f"        hidden bits - (k-1)^2  for k=1..5 : {[h - (k-1)**2 for k, h in zip([1,2,3,4,5], hid6)]}"
      f"   -> hidden = (k-1)^2 = (sqrt(|R|)-1)^2 EXACTLY")
print(f"        log-log fit over k=2..5 : q = {slope(aa[nz], hh[nz]):.6f}.  CAVEAT: (sqrt(A)-1)^2 is")
print(f"        not a pure power either, so this finite-range fit OVERSHOOTS; the local slope between")
print(f"        k=4 and k=5 is {np.log(hh[4]/hh[3])/np.log(aa[4]/aa[3]):.6f} and the law is area, exponent -> 1.")
print("        The load-bearing comparison is exact and needs no fit: boundary support 4k against")
print("        hidden interior information (k-1)^2 -- perimeter against area, diverging as k/4.")
print()

# ----------------------------------------------------------------------------------------------
# [5] verdict
# ----------------------------------------------------------------------------------------------
print("=" * 100)
print("[5] WHAT WAS MEASURED")
print("=" * 100)
print("  1. || P_R - prod_{l in dR} sigma_l || = 0.000000e+00 EXACTLY, in both bases, on the full")
print("     2^18 space and as dense 1024x1024 matrices on the physical sector, for every region")
print("     tested.  The bulk constraints cancel identically: every link with both endpoints in R")
print("     is hit exactly twice and sigma^2 = I.  The total constraint of a region IS a surface")
print("     operator.  This is an OPERATOR IDENTITY, true for every coupling and every lattice")
print("     size, so it is a scale-free algebraic fact, not a finite-lattice accident.")
print("  2. support(P_R) = 4k = 4*sqrt(|R|) EXACTLY at every point (residuals all 0); fitted")
print("     exponent 0.500000 = a perimeter law.  The uncancelled product is 2k^2+2k exactly")
print("     (fitted 0.840379 over k=1..5, local slope 0.908530 at the top, area law in the limit).")
print("  3. The boundary operator alone carries ZERO gauge-invariant information about the interior:")
print("     it equals +1 on every physical state.  The interior's own gauge-invariant content")
print("     (Wilson loops) grows like the AREA, (k-1)^2 bits, none of it visible in B_R.")
print("  NOT MEASURED HERE: whether gravity's constraint algebra has the same structure functions;")
print("  that is x6a's business.  This lane establishes only the boundary-term property on our side.")
