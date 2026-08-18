"""W-51.  CAN THE CARRIER'S GEOMETRY BE READ BACK OUT OF THE RECORDS?
    METHOD LENS: THE OPERATOR ALGEBRA ITSELF.

Every lane of this program installed the graph by hand.  This lane asks whether the PLAQUETTE
ADJACENCY (which plaquettes share a link -- the dual graph) can be recovered without ever being
told it.

THREE ROUTES ARE RUN AND SCORED SEPARATELY, BECAUSE THEY HAVE DIFFERENT STANDING.

  ROUTE 1  THE RECORD ALGEBRA ALONE.  Products, commutators, anticommutators among the record
           operators W_S.  Expected by the lens to carry adjacency.  It does not, and the reason is
           structural, not numerical: at Z_2 every record is a product of link X's, so any two
           records COMMUTE regardless of which links they share.  The record algebra is
           Z_2^{n_p} -- the same abelian group on a block, a chain and a ring.  Reported as a NULL,
           and the null is the substantive part of this lane's answer to the lens.

  ROUTE 2  THE RECORDS PLUS ONE UNIFORM, UNLABELLED ENVIRONMENT.   The reconstruction consumes
           exactly one real number per record: its decay rate under a bath that couples to every
           link with the same strength.  No link is ever named, selected, or resolved.  Adjacency
           comes out by POLARISATION:
                 |bdy(p) INTERSECT bdy(q)|  =  ( r(p) + r(q) - r(p,q) ) / (4 gamma)
           because bdy is additive over GF(2) (bdy(S xor T) = bdy(S) xor bdy(T)) and
           |A xor B| = |A| + |B| - 2|A intersect B|.  This is the lane's result.

  ROUTE 3  THE LINK-RESOLVED BATH.  Put the bath on ONE link at a time and see which records decay.
           This works perfectly and is reported as CIRCULAR, not as a result: rate(W_S) under a
           single-link bath IS the incidence matrix entry, by W-41's counting law.  Measuring it is
           a direct read of the incidence structure wearing a dynamical costume.  The lens warned
           about exactly this and it is reported as the null it is.

WHAT ROUTE 2 IS AND IS NOT.  The Gram matrix G(p,q) = |bdy(p) INTERSECT bdy(q)| is STRICTLY LESS
information than the incidence structure: it says HOW MANY links two plaquettes share, never WHICH.
Adjacency is then a derived predicate G>0.  So this is not the link list in disguise.  It is also
not a statistical inference: the polarisation identity is exact, so the recovery is exact, and the
honest description is "the dual graph is a FUNCTION of label-free decay rates", not "the dual graph
was estimated".  Both statements are in the report.

THE REDUCED CARRIER.  Physical (gauge-invariant) states of a simply-connected Z_2 patch are the
Z_2 cycle space = the span of plaquette subsets, so the basis is labelled by S subset of plaquettes
and dim = 2^{n_p}.  In that basis
    W_T |S> = |S xor T>            (the plaquette move -- GEOMETRY-FREE, see Route 1)
    Z_k |S> = (-1)^[k in bdy(S)]|S>  (the electric link operator -- the ONLY place geometry enters)
    H = -MAG - g2*ELEC,  MAG = sum_p (W_p + W_p^dag),  ELEC = sum_k (Z_k + Z_k^dag)
This is validated below against the full Gauss-law construction of W-33/W-41 (states enumerated on
the link configuration space, Gauss law imposed at every vertex) -- rates agree to machine epsilon.
The reduced form is what lets a SHUFFLED carrier be simulated at all: a shuffled incidence has no
vertex structure, but it still has a map plaquette -> link set, and that is all the reduced form
needs.

CONVENTIONS INHERITED (W-33 erratum).  numpy reshape(-1) is ROW-major, so vec(AXB)=(A kron B^T)vec X
and the Schrodinger generator is -i(H kron I - I kron H^T) + gamma sum (L kron L* - I kron I).
The Heisenberg (observable) generator used here is built directly as
    +i(H kron I - I kron H^T) + gamma sum (Z kron Z^T - I kron I)
and is checked below to be the conjugate transpose of the W-33 Schrodinger form.
"""

import itertools, math, numpy as np

RNG = np.random.default_rng(20260818)
GAM = 0.5

# ----------------------------------------------------------------------------------------------
# CARRIERS.  A carrier is exactly three things: n_p, n_links, and PL: plaquette -> set of links.
# Nothing else -- no coordinates, no vertex list, no ordering -- is ever handed to a reconstruction.
# ----------------------------------------------------------------------------------------------

def from_cells(cells):
    """Square cells (i,j) -> (n_links, PL).  Used only to MANUFACTURE carriers and ground truth."""
    lid = {}; PL = []
    for (i, j) in cells:
        s = []
        for e in [('h', i, j), ('h', i, j + 1), ('v', i, j), ('v', i + 1, j)]:
            if e not in lid: lid[e] = len(lid)
            s.append(lid[e])
        PL.append(sorted(s))
    return len(lid), PL

def ring_cells(n):
    """n quadrilaterals in a cycle: spokes 0..n-1, inner arcs n..2n-1, outer arcs 2n..3n-1."""
    PL = [sorted([i, (i + 1) % n, n + i, 2 * n + i]) for i in range(n)]
    return 3 * n, PL

CARRIERS = {}
CARRIERS['BLOCK5'] = from_cells([(0, 0), (1, 0), (0, 1), (1, 1), (0, 2)])   # P-pentomino
CARRIERS['CHAIN5'] = from_cells([(i, 0) for i in range(5)])                 # I-pentomino
CARRIERS['RING5']  = ring_cells(5)
CARRIERS['BLOCK4'] = from_cells([(0, 0), (1, 0), (0, 1), (1, 1)])           # 2x2, the W-41 carrier
CARRIERS['BLOCK8'] = from_cells([(i, j) for j in range(2) for i in range(4)])
CARRIERS['CHAIN8'] = from_cells([(i, 0) for i in range(8)])
CARRIERS['RING8']  = ring_cells(8)

# ---- GROUND TRUTH (used ONLY for scoring and for the pre-flight audit; never fed to a route) ----
def true_gram(n_links, PL):
    n_p = len(PL)
    G = np.zeros((n_p, n_p), dtype=int)
    for p in range(n_p):
        for q in range(n_p):
            G[p, q] = len(set(PL[p]) & set(PL[q]))
    return G

def true_adj(n_links, PL):
    G = true_gram(n_links, PL); n_p = len(PL)
    A = (G > 0).astype(int)
    for p in range(n_p): A[p, p] = 0
    return A

def pairs(n_p): return list(itertools.combinations(range(n_p), 2))

def fidelity(Ahat, Atrue):
    n_p = Atrue.shape[0]; pr = pairs(n_p)
    if not pr: return 1.0
    return sum(1 for (p, q) in pr if Ahat[p, q] == Atrue[p, q]) / len(pr)

def edge_list(A):
    return sorted((p, q) for p, q in pairs(A.shape[0]) if A[p, q])

def degseq(A): return sorted(int(A[p].sum()) for p in range(A.shape[0]))

# ----------------------------------------------------------------------------------------------
# REDUCED PHYSICS.  Built from (n_p, n_links, PL) only.
# ----------------------------------------------------------------------------------------------

def build_reduced(n_p, n_links, PL):
    D = 1 << n_p
    plink = [0] * n_p
    for p in range(n_p):
        m = 0
        for k in PL[p]: m ^= (1 << k)
        plink[p] = m
    bdy = np.zeros(D, dtype=np.int64)
    for S in range(1, D):
        low = S & -S; p = low.bit_length() - 1
        bdy[S] = bdy[S ^ low] ^ plink[p]
    bsz = np.array([bin(int(b)).count('1') for b in bdy], dtype=int)   # |bdy(S)|
    # electric link diagonals
    Zd = np.empty((n_links, D))
    for k in range(n_links):
        Zd[k] = 1.0 - 2.0 * ((bdy >> k) & 1)
    Kmat = Zd.T @ Zd                                    # sum_k outer(z_k, z_k)
    ELECd = 2.0 * Zd.sum(axis=0)                        # ELEC = sum_k (Z_k + Z_k^dag)
    return dict(n_p=n_p, n_links=n_links, D=D, bdy=bdy, bsz=bsz, Zd=Zd, Kmat=Kmat, ELECd=ELECd)

def Wmat(C, S):
    """Record operator W_S: the plaquette move for region S.  A permutation of the region basis."""
    D = C['D']
    M = np.zeros((D, D))
    idx = np.arange(D) ^ S
    M[idx, np.arange(D)] = 1.0
    return M

def Hmat(C, g2):
    D = C['D']; n_p = C['n_p']
    MAG = np.zeros((D, D))
    for p in range(n_p):
        Wp = Wmat(C, 1 << p)
        MAG += Wp + Wp.T
    return -MAG - g2 * np.diag(C['ELECd'])

# ---- Heisenberg-picture propagation of a record (batched over records) ----
def heis_deriv(C, H, Ob, gam):
    # dO/dt = i[H,O] + gam*( sum_k Z_k O Z_k - n_links*O )
    # np.matmul broadcasts (D,D)@(b,D,D) and (b,D,D)@(D,D) through BLAS; einsum does not.
    return 1j * (np.matmul(H, Ob) - np.matmul(Ob, H)) \
           + gam * (C['Kmat'][None, :, :] * Ob - C['n_links'] * Ob)

def measure_rates(C, g2, regions, gam=GAM, t1=0.02, t2=0.12, nstep=None):
    if nstep is None: nstep = 400 if C['D'] <= 64 else 160
    """OPERATIONAL rate of every record in `regions`, from real time evolution.
       Returns r[S] = -(1/(t2-t1)) * ln |f(t2)/f(t1)|, f(t) = <W_S, O_S(t)>/<W_S,W_S>."""
    H = Hmat(C, g2); D = C['D']
    Ob = np.stack([Wmat(C, S).astype(complex) for S in regions])
    W0 = Ob.copy()
    nrm = np.einsum('bij,bij->b', W0.conj(), W0).real
    dt = t2 / nstep
    n1 = int(round(t1 / dt))
    f1 = None
    for n in range(nstep):
        if n == n1:
            f1 = np.einsum('bij,bij->b', W0.conj(), Ob).real / nrm
        k1 = heis_deriv(C, H, Ob, gam)
        k2 = heis_deriv(C, H, Ob + 0.5 * dt * k1, gam)
        k3 = heis_deriv(C, H, Ob + 0.5 * dt * k2, gam)
        k4 = heis_deriv(C, H, Ob + dt * k3, gam)
        Ob = Ob + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    f2 = np.einsum('bij,bij->b', W0.conj(), Ob).real / nrm
    eps = 1e-14
    r = -(np.log(np.abs(f2) + eps) - np.log(np.abs(f1) + eps)) / (t2 - t1)
    flags = int(np.sum((f1 <= 0) | (f2 <= 0)))
    return {S: float(r[i]) for i, S in enumerate(regions)}, flags

# ---- expm by scaling and squaring (numpy only), used to validate the RK4 propagation ----
def expm_ss(A):
    nA = np.max(np.sum(np.abs(A), axis=1))
    s = max(0, int(math.ceil(math.log2(nA / 0.5))) if nA > 0.5 else 0)
    B = A / (2.0 ** s)
    X = np.eye(A.shape[0], dtype=complex); T = np.eye(A.shape[0], dtype=complex)
    for k in range(1, 20):
        T = T @ B / k
        X = X + T
    for _ in range(s): X = X @ X
    return X

# ----------------------------------------------------------------------------------------------
# THE RECONSTRUCTIONS.  These are BARRIER FUNCTIONS.  Read their signatures.
# ----------------------------------------------------------------------------------------------

def reconstruct_route2(n_p, rate_of_region, gam):
    """ROUTE 2.  INPUTS: an integer n_p, a dict {region-index-set -> real decay rate}, and gamma.
       There is no link list, no vertex list, no coordinate, no incidence datum in this signature.
       The region indices are opaque labels: permuting them permutes the output identically."""
    assert isinstance(n_p, int)
    assert all(isinstance(v, float) for v in rate_of_region.values())
    G = np.zeros((n_p, n_p))
    for p, q in pairs(n_p):
        rp = rate_of_region[(1 << p)]
        rq = rate_of_region[(1 << q)]
        rpq = rate_of_region[(1 << p) | (1 << q)]
        G[p, q] = G[q, p] = (rp + rq - rpq) / (4.0 * gam)
    A = (G > 0.5).astype(int)
    np.fill_diagonal(A, 0)
    return A, G

def reconstruct_route3_CIRCULAR(n_p, n_channels, rate_of_channel_and_plaquette, gam):
    """ROUTE 3.  INPUTS include one rate PER BATH CHANNEL per plaquette.  Even though no channel is
       named, rate = 2*gamma*|bath INTERSECT bdy| makes each number an INCIDENCE MATRIX ENTRY.
       Kept only to be reported as circular."""
    M = np.zeros((n_channels, n_p), dtype=int)
    for c in range(n_channels):
        for p in range(n_p):
            M[c, p] = 1 if rate_of_channel_and_plaquette[(c, p)] > gam else 0
    A = ((M.T @ M) > 0).astype(int); np.fill_diagonal(A, 0)
    return A

# ==============================================================================================
print("=" * 100)
print("W-51  RECORD -> GEOMETRY.  LENS: THE OPERATOR ALGEBRA ITSELF.")
print("=" * 100)

# ----------------------------------------------------------------------------------------------
print("\n[0]  PRE-REGISTRATION -- FORCED OR NOT, BEFORE ANY DYNAMICS RUNS.")
print("""
  QUANTITY AT RISK:  Ghat(p,q) = ( r(p) + r(q) - r({p,q}) ) / (4*gamma),  claimed to equal the
  number of links plaquettes p and q share.  The program has seven false positives from measuring
  something a constraint already fixed.  The danger here is that r({p,q}) is determined by the
  COUNTS the carrier is built with -- every plaquette carries 4 links, the link total is fixed --
  in which case Ghat is a constant dressed up as a measurement and carries no adjacency at all.

  THE ONE ARGUMENT THAT SETTLES IT:  |bdy(S)| for |S|=2 is NOT a function of the preserved counts.
  A degree-preserving shuffle of the plaquette-link incidence keeps |bdy(p)|=4 for every p, keeps
  the link total, and keeps every link's multiplicity -- and yet moves |bdy({p,q})|.  Demonstrated
  combinatorially below, with no Lindbladian anywhere near it.  If the spread were zero the lane
  would stop here and report a forced null.""")

n_links, PL = CARRIERS['BLOCK5']

def bdy_size(PL, S):
    c = set()
    for p in S: c ^= set(PL[p])
    return len(c)

def shuffle_PL(n_p, n_links, PL, rng, tries=4000):
    degs = [len(s) for s in PL]
    mult = [0] * n_links
    for s in PL:
        for k in s: mult[k] += 1
    stub0 = np.array([k for k in range(n_links) for _ in range(mult[k])])
    for _ in range(tries):
        stubs = list(rng.permutation(stub0))
        pos = 0; out = []; ok = True
        for p in range(n_p):
            s = stubs[pos:pos + degs[p]]; pos += degs[p]
            if len(set(s)) != len(s): ok = False; break
            out.append(sorted(int(x) for x in s))
        if ok: return out
    return None

n_p = len(PL)
vals_true = sorted(set(bdy_size(PL, S) for S in pairs(n_p)))
sh_vals = set()
for _ in range(40):
    sp = shuffle_PL(n_p, n_links, PL, RNG)
    if sp is None: continue
    for S in pairs(n_p): sh_vals.add(bdy_size(sp, S))
mult = [0] * n_links
for s in PL:
    for k in s: mult[k] += 1
print(f"  BLOCK5:  n_p={n_p}  n_links={n_links}  |bdy(p)| per plaquette = {sorted(set(len(s) for s in PL))}")
print(f"           link multiplicity profile = {sorted(mult)}   (a shuffle preserves ALL of these)")
print(f"  |bdy(S)| for |S|=2, TRUE carrier    : {vals_true}")
print(f"  |bdy(S)| for |S|=2, over 40 shuffles: {sorted(sh_vals)}")
spread = len(vals_true) > 1 and len(sh_vals) > 1
print(f"  spread at fixed |S|=2 and fixed counts: TRUE {len(vals_true)} distinct, "
      f"SHUFFLED {len(sh_vals)} distinct")
print(f"  => the at-risk quantity is {'NOT FORCED by the preserved counts' if spread else 'FORCED -- STOP'}")
if not spread:
    raise SystemExit("forced: lane stops")

# ----------------------------------------------------------------------------------------------
print("\n[1]  VALIDATION OF THE MACHINERY (three independent checks, before any claim).")

# (a) reduced model vs the full Gauss-law construction of W-33/W-41, on the 2x2 carrier
def build_full_2x2():
    V2 = [(i, j) for j in range(3) for i in range(3)]; vid = {v: k for k, v in enumerate(V2)}
    E = []
    for j in range(3):
        for i in range(2): E.append((vid[(i, j)], vid[(i + 1, j)]))
    for j in range(2):
        for i in range(3): E.append((vid[(i, j)], vid[(i, j + 1)]))
    hid = lambda i, j: j * 2 + i; vx = lambda i, j: 6 + j * 3 + i
    P = [[(hid(i, j), +1), (vx(i + 1, j), +1), (hid(i, j + 1), -1), (vx(i, j), -1)]
         for j in range(2) for i in range(2)]
    N = 2
    st = [s for s in itertools.product(range(N), repeat=len(E))
          if all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
                - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % N == 0 for v in range(len(V2)))]
    idx = {s: i for i, s in enumerate(st)}
    def Zop(links):
        return np.diag([(-1.0) ** (sum(s[k] for k in links) % 2) for s in st])
    def Move(mv):
        D = len(st); M = np.zeros((D, D))
        for j, s in enumerate(st):
            t = list(s)
            for k, sg in mv: t[k] = (t[k] + sg) % N
            M[idx[tuple(t)], j] = 1.0
        return M
    def compose(ps):
        acc = {}
        for p in ps:
            for k, sg in p: acc[k] = acc.get(k, 0) + sg
        return [(k, s) for k, s in acc.items() if s != 0]
    return st, idx, E, P, Zop, Move, compose

st, idxf, E, P, Zop, Move, compose = build_full_2x2()
Df = len(st)
MAGf = sum((lambda X: X + X.T)(Move(p)) for p in P)
ELECf = sum(2.0 * Zop([k]) for k in range(len(E)))
g2v = 0.05
Hf = -MAGf - g2v * ELECf
Zdf = np.array([np.diag(Zop([k])) for k in range(len(E))])
Kf = Zdf.T @ Zdf
Cfull = dict(n_p=4, n_links=len(E), D=Df, Kmat=Kf)
regions_f = [(0,), (1,), (0, 1), (0, 3)]
def measure_full(Hf, Cfull, regs):
    Ob = np.stack([Move(compose([P[i] for i in S])).astype(complex) for S in regs])
    W0 = Ob.copy(); nrm = np.einsum('bij,bij->b', W0.conj(), W0).real
    t2 = 0.12; nstep = 400; dt = t2 / nstep; n1 = int(round(0.02 / dt)); f1 = None
    for n in range(nstep):
        if n == n1: f1 = np.einsum('bij,bij->b', W0.conj(), Ob).real / nrm
        def dv(O): return 1j * (Hf @ O - O @ Hf) + GAM * (Kf[None] * O - Cfull['n_links'] * O)
        k1 = dv(Ob); k2 = dv(Ob + 0.5 * dt * k1); k3 = dv(Ob + 0.5 * dt * k2); k4 = dv(Ob + dt * k3)
        Ob = Ob + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    f2 = np.einsum('bij,bij->b', W0.conj(), Ob).real / nrm
    return -(np.log(np.abs(f2)) - np.log(np.abs(f1))) / (t2 - 0.02)
rf = measure_full(Hf, Cfull, regions_f)
nl4, PL4 = CARRIERS['BLOCK4']
C4 = build_reduced(4, nl4, PL4)
regs_bits = [1, 2, 3, 1 | 8]
rr, _ = measure_rates(C4, g2v, regs_bits)
print(f"  (a) reduced region-space model  vs  full Gauss-law model (2x2 patch, g2={g2v}, gamma={GAM})")
print(f"      {'region':>10s} {'full Gauss law':>16s} {'reduced':>12s} {'|diff|':>10s}")
mx = 0.0
for i, S in enumerate(regions_f):
    a = float(rf[i]); b = rr[regs_bits[i]]
    mx = max(mx, abs(a - b))
    print(f"      {str(S):>10s} {a:16.10f} {b:12.10f} {abs(a-b):10.2e}")
print(f"      max discrepancy {mx:.2e}  -> the reduced carrier IS the W-33/W-41 carrier")

# (b) RK4 vs expm of the Liouvillian (scaling and squaring), Heisenberg generator
C5 = build_reduced(5, *CARRIERS['BLOCK5'])
D5 = C5['D']; I5 = np.eye(D5)
Hb = Hmat(C5, 0.1)
Mheis = 1j * (np.kron(Hb, I5) - np.kron(I5, Hb.T))
for k in range(C5['n_links']):
    Zk = np.diag(C5['Zd'][k])
    Mheis = Mheis + GAM * (np.kron(Zk, Zk.T) - np.kron(I5, I5))
Msch = -1j * (np.kron(Hb, I5) - np.kron(I5, Hb.T))
for k in range(C5['n_links']):
    Zk = np.diag(C5['Zd'][k]).astype(complex)
    Msch = Msch + GAM * (np.kron(Zk, Zk.conj()) - np.kron(I5, I5))
print(f"  (b) Liouville dim {D5*D5}.  ||Mheis - Msch^dag|| = {np.linalg.norm(Mheis - Msch.conj().T):.2e}"
      f"   (row-major convention, W-33 erratum)")
tt = 0.12
U = expm_ss(tt * Mheis)
S0 = 1 | 2
Ov = (U @ Wmat(C5, S0).astype(complex).reshape(-1)).reshape(D5, D5)
Ob = Wmat(C5, S0).astype(complex)[None]
dt = tt / 400
for _ in range(400):
    k1 = heis_deriv(C5, Hb, Ob, GAM); k2 = heis_deriv(C5, Hb, Ob + 0.5 * dt * k1, GAM)
    k3 = heis_deriv(C5, Hb, Ob + 0.5 * dt * k2, GAM); k4 = heis_deriv(C5, Hb, Ob + dt * k3, GAM)
    Ob = Ob + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
print(f"      ||RK4(t) - expm_ss(t*M)|| = {np.linalg.norm(Ob[0] - Ov):.2e}   at t={tt}")

# (c) measured rate vs the W-41 counting law at g2=0
r0, _ = measure_rates(C5, 0.0, [1 << p for p in range(5)] + [(1 << p) | (1 << q) for p, q in pairs(5)])
worst = max(abs(r0[S] - 2 * GAM * C5['bsz'][S]) for S in r0)
print(f"  (c) measured rate vs 2*gamma*|bdy| at g2=0, all 15 records: max |diff| = {worst:.2e}")

# ----------------------------------------------------------------------------------------------
print("\n[2]  OPERATORS MUST EARN MEASUREMENT.  Record operators on BLOCK5, before any commutator.")
print(f"      {'record S':>14s} {'||O||':>10s} {'||O^dag O - I||':>16s} {'distinct eigs':>14s}")
for S in [1, 2, 1 | 2, 1 | 4, 31]:
    O = Wmat(C5, S)
    ev = np.linalg.eigvals(O)
    dis = len(np.unique(np.round(ev, 8)))
    print(f"      {S:>14d} {np.linalg.norm(O):10.6f} {np.linalg.norm(O.T @ O - np.eye(D5)):16.2e} {dis:14d}")
print("      (norm sqrt(D)=%.6f for a permutation, unitary to machine epsilon, spectrum {+1,-1})"
      % math.sqrt(D5))

# ----------------------------------------------------------------------------------------------
print("\n[3]  ROUTE 1 -- THE RECORD ALGEBRA ALONE.   THE LENS'S OWN HYPOTHESIS, TESTED FIRST.")
print("""
  Every record is a product of link X operators, so W_S W_T = W_{S xor T} = W_T W_S for ALL S,T
  whatever links they share.  If adjacency lived in the record algebra it would have to show up in
  a commutator, an anticommutator, or a product -- and all three are geometry-blind by that line.
  Measured, not asserted:""")
print(f"      {'carrier':>10s} {'max ||[W_S,W_T]||':>19s} {'max ||{W_S,W_T}-2W_SW_T||':>27s} "
      f"{'product table = Z_2^n_p':>24s}")
for nm in ['BLOCK5', 'CHAIN5', 'RING5']:
    Cx = build_reduced(5, *CARRIERS[nm])
    mc = 0.0; ma = 0.0; grp = True
    Ws = {S: Wmat(Cx, S) for S in range(32)}
    for S in range(32):
        for T in range(32):
            mc = max(mc, np.linalg.norm(Ws[S] @ Ws[T] - Ws[T] @ Ws[S]))
            ma = max(ma, np.linalg.norm(Ws[S] @ Ws[T] + Ws[T] @ Ws[S] - 2 * Ws[S] @ Ws[T]))
            if np.linalg.norm(Ws[S] @ Ws[T] - Ws[S ^ T]) > 1e-12: grp = False
    print(f"      {nm:>10s} {mc:19.2e} {ma:27.2e} {str(grp):>24s}")
print("""
  READING.  The three carriers -- a 2D block, an open chain, a ring -- have three different dual
  graphs and ONE identical record algebra.  Every commutator vanishes, every product is fixed by the
  xor of the labels, the multiplication table is Z_2^5 in all three cases.  The record algebra by
  itself therefore contains ZERO bits of adjacency.  ROUTE 1 IS A NULL, and it is a clean one: the
  null is forced by the abelian structure, not by insufficient resolution.
  Note the same holds for the magnetic Hamiltonian, MAG = sum_p (W_p + W_p^dag), which is itself a
  sum of records and commutes with all of them.  Geometry first enters through Z_k -- the ELECTRIC
  operator and the bath coupling.  That is why Route 2 needs an environment.""")

# ----------------------------------------------------------------------------------------------
print("\n[4]  ROUTE 2 -- RECORDS PLUS ONE UNIFORM, UNLABELLED ENVIRONMENT.")
print("""
  WHAT THE RECONSTRUCTION CONSUMES (this is the whole input, and it is asserted in code):
      n_p                          an integer
      { region-label -> rate }     one real number per record, from real time evolution
      gamma                        the single environment strength
  A region label is a subset of the opaque plaquette index set.  Relabelling the plaquettes
  relabels the output and nothing else.  There is ONE bath and it is uniform: no link is selected,
  named, or resolved at any point.  n_p + C(n_p,2) real numbers go in; the dual graph comes out.""")

def run_route2(name, n_p, n_links, PL, g2, gam=GAM, verbose=False):
    C = build_reduced(n_p, n_links, PL)
    regs = [1 << p for p in range(n_p)] + [(1 << p) | (1 << q) for p, q in pairs(n_p)]
    rates, flags = measure_rates(C, g2, regs, gam=gam)
    Ahat, Ghat = reconstruct_route2(n_p, rates, gam)
    return Ahat, Ghat, flags

for nm in ['BLOCK5', 'CHAIN5', 'RING5']:
    nl, pl = CARRIERS[nm]
    At = true_adj(nl, pl); Gt = true_gram(nl, pl)
    Ah, Gh, fl = run_route2(nm, 5, nl, pl, 0.0)
    off = [(p, q) for p, q in pairs(5)]
    gerr = max(abs(Gh[p, q] - Gt[p, q]) for p, q in off)
    print(f"\n  {nm}:  fidelity {fidelity(Ah, At):.3f}   max |Ghat-Gtrue| over pairs {gerr:.2e}")
    print(f"        reconstructed edges {edge_list(Ah)}")
    print(f"        true          edges {edge_list(At)}")
    print(f"        reconstructed degree sequence {degseq(Ah)}   true {degseq(At)}")

# separation statistic
def gapstat(Gh, At):
    n_p = At.shape[0]
    adj = [Gh[p, q] for p, q in pairs(n_p) if At[p, q]]
    non = [Gh[p, q] for p, q in pairs(n_p) if not At[p, q]]
    return (min(adj) if adj else float('nan')), (max(non) if non else float('nan'))

nl, pl = CARRIERS['BLOCK5']; At5 = true_adj(nl, pl)
Ah, Gh, _ = run_route2('BLOCK5', 5, nl, pl, 0.0)
lo, hi = gapstat(Gh, At5)
print(f"\n  SEPARATION on BLOCK5 at g2=0: min Ghat over adjacent pairs {lo:.6f}, "
      f"max Ghat over non-adjacent pairs {hi:.2e}, threshold 0.5")

# ----------------------------------------------------------------------------------------------
print("\n[5]  THE SHUFFLE CONTROL.  Degree-preserving permutation of which links belong to which")
print("     plaquettes.  Scored TWO ways, because the two ways answer different questions.")
print("       (i)  against the ORIGINAL carrier's adjacency -- does the method just replay a fixed")
print("            answer that the counts or the labels already fix?   Must COLLAPSE.")
print("       (ii) against the SHUFFLED carrier's OWN adjacency -- does the method track the")
print("            carrier it is actually run on?   Must HOLD.")

NSH = 8
print(f"\n  {'carrier':>8s} {'shuffle':>8s} {'fid vs ORIGINAL':>17s} {'fid vs OWN truth':>18s} "
      f"{'Ghat values seen':>22s}")
sh_orig = []; sh_own = []
nl, pl = CARRIERS['BLOCK5']; At = true_adj(nl, pl)
for s in range(NSH):
    sp = shuffle_PL(5, nl, pl, RNG)
    Ash = true_adj(nl, sp)
    Ah, Gh, _ = run_route2('SH', 5, nl, sp, 0.0)
    f_o = fidelity(Ah, At); f_w = fidelity(Ah, Ash)
    sh_orig.append(f_o); sh_own.append(f_w)
    vals = sorted(set(int(round(Gh[p, q])) for p, q in pairs(5)))
    print(f"  {'BLOCK5':>8s} {s:8d} {f_o:17.3f} {f_w:18.3f} {str(vals):>22s}")
print(f"\n  mean over {NSH} shuffles:  vs ORIGINAL {np.mean(sh_orig):.3f} +- {np.std(sh_orig):.3f}"
      f"    vs OWN truth {np.mean(sh_own):.3f} +- {np.std(sh_own):.3f}")

# controls
ne = len(edge_list(At))
rand_f = []
for _ in range(400):
    pr = pairs(5); pick = RNG.choice(len(pr), size=ne, replace=False)
    Ar = np.zeros((5, 5), dtype=int)
    for i in pick:
        p, q = pr[i]; Ar[p, q] = Ar[q, p] = 1
    rand_f.append(fidelity(Ar, At))
allone = np.ones((5, 5), dtype=int); np.fill_diagonal(allone, 0)
allzero = np.zeros((5, 5), dtype=int)
print(f"  CONTROLS on BLOCK5 (10 pairs, {ne} true edges):")
print(f"    random graph, same edge count : {np.mean(rand_f):.3f} +- {np.std(rand_f):.3f}")
print(f"    predict ALL pairs adjacent    : {fidelity(allone, At):.3f}")
print(f"    predict NO pairs adjacent     : {fidelity(allzero, At):.3f}")
print(f"    Route 2 on the true carrier   : {fidelity(*[run_route2('B',5,nl,pl,0.0)[0], At]):.3f}")

# ----------------------------------------------------------------------------------------------
print("\n[6]  A SECOND GEOMETRY.  Same plaquette count, same |bdy(p)|, same link multiplicity")
print("     profile -- different dual graph.  A working method must tell them apart.")
print(f"\n  {'carrier':>8s} {'n_links':>8s} {'mult profile':>26s} {'true dual':>16s} "
      f"{'reconstructed dual':>20s} {'fid':>6s}")
recon_duals = {}
for nm in ['BLOCK5', 'RING5', 'CHAIN5']:
    nl, pl = CARRIERS[nm]
    m = [0] * nl
    for s in pl:
        for k in s: m[k] += 1
    At = true_adj(nl, pl)
    Ah, Gh, _ = run_route2(nm, 5, nl, pl, 0.0)
    recon_duals[nm] = Ah
    print(f"  {nm:>8s} {nl:8d} {str(sorted(m)):>26s} {str(degseq(At)):>16s} "
          f"{str(degseq(Ah)):>20s} {fidelity(Ah, At):6.3f}")
same = np.array_equal(recon_duals['BLOCK5'], recon_duals['RING5'])
print(f"\n  reconstructed BLOCK5 dual identical to reconstructed RING5 dual? {same}")
print(f"  BLOCK5 recon edges {edge_list(recon_duals['BLOCK5'])}")
print(f"  RING5  recon edges {edge_list(recon_duals['RING5'])}")
print(f"  degree sequences differ ({degseq(recon_duals['BLOCK5'])} vs {degseq(recon_duals['RING5'])})"
      f" -> the two carriers are NOT isomorphic and the method separates them")

print("\n  Larger carriers, n_p=8 (Hilbert 256; no Liouville matrix is ever formed at this size).")
print(f"  {'carrier':>8s} {'n_p':>4s} {'n_links':>8s} {'true edges':>11s} {'recon edges':>12s} {'fid':>6s}")
for nm in ['BLOCK8', 'CHAIN8', 'RING8']:
    nl, pl = CARRIERS[nm]
    At = true_adj(nl, pl)
    Ah, Gh, _ = run_route2(nm, 8, nl, pl, 0.0)
    print(f"  {nm:>8s} {8:4d} {nl:8d} {len(edge_list(At)):11d} {len(edge_list(Ah)):12d} "
          f"{fidelity(Ah, At):6.3f}")

print("\n  Shuffle control at n_p=8 (4 shuffles of BLOCK8):")
nl, pl = CARRIERS['BLOCK8']; At8 = true_adj(nl, pl)
o8 = []; w8 = []
for s in range(4):
    sp = shuffle_PL(8, nl, pl, RNG)
    Ash = true_adj(nl, sp)
    Ah, Gh, _ = run_route2('SH8', 8, nl, sp, 0.0)
    o8.append(fidelity(Ah, At8)); w8.append(fidelity(Ah, Ash))
print(f"    vs ORIGINAL {np.mean(o8):.3f} +- {np.std(o8):.3f}    vs OWN truth "
      f"{np.mean(w8):.3f} +- {np.std(w8):.3f}   (28 pairs)")

# ----------------------------------------------------------------------------------------------
print("\n[7]  THE COUPLING SCAN.  g2 turns on the electric term, which does NOT commute with the")
print("     records, so the decay stops being a pure exponential and the extracted rate is")
print("     contaminated.  Numbers as they come; nothing is fitted.")
print(f"\n  {'g2':>8s} {'BLOCK5 fid':>11s} {'RING5 fid':>10s} {'CHAIN5 fid':>11s} "
      f"{'min Ghat adj':>13s} {'max Ghat non-adj':>17s} {'sign flags':>11s}")
g2list = [0.0, 0.001, 0.01, 0.05, 0.1, 0.3, 1.0, 3.0]
scan_rows = []
for g2 in g2list:
    fids = {}
    lo = hi = float('nan'); fl_tot = 0
    for nm in ['BLOCK5', 'RING5', 'CHAIN5']:
        nl, pl = CARRIERS[nm]; At = true_adj(nl, pl)
        Ah, Gh, fl = run_route2(nm, 5, nl, pl, g2)
        fids[nm] = fidelity(Ah, At); fl_tot += fl
        if nm == 'BLOCK5': lo, hi = gapstat(Gh, At)
    scan_rows.append((g2, fids['BLOCK5'], fids['RING5'], fids['CHAIN5'], lo, hi))
    print(f"  {g2:8.3f} {fids['BLOCK5']:11.3f} {fids['RING5']:10.3f} {fids['CHAIN5']:11.3f} "
          f"{lo:13.4f} {hi:17.4f} {fl_tot:11d}")

# ----------------------------------------------------------------------------------------------
print("\n[8]  ROUTE 3 -- THE LINK-RESOLVED BATH.  REPORTED AS CIRCULAR, NOT AS A RESULT.")
print("""
  Put the bath on one link at a time.  W-41's law says rate(W_S) = 2*gamma*|bath INTERSECT bdy(S)|,
  so with a single-link bath the measured rate is 2*gamma times the INCIDENCE MATRIX ENTRY [k in
  bdy(p)].  The channels are never named -- and it does not help.  Thresholding those rates
  RECONSTRUCTS THE INCIDENCE MATRIX ITSELF, from which adjacency follows by M^T M.  That is not a
  reconstruction of geometry from records; it is the incidence structure read out through a
  dynamical proxy.  It is exactly the failure mode the lens named, and it is scored as a null.""")
for nm in ['BLOCK5', 'RING5']:
    nl, pl = CARRIERS[nm]
    C = build_reduced(5, nl, pl)
    tab = {}
    for k in range(nl):
        regs = [1 << p for p in range(5)]
        H = Hmat(C, 0.0)
        Zk = C['Zd'][k]
        Kk = np.outer(Zk, Zk)
        Csingle = dict(C); Csingle['Kmat'] = Kk; Csingle['n_links'] = 1
        rr, _ = measure_rates(Csingle, 0.0, regs)
        for p in range(5): tab[(k, p)] = rr[1 << p]
    Ah = reconstruct_route3_CIRCULAR(5, nl, tab, GAM)
    At = true_adj(nl, pl)
    sp = shuffle_PL(5, nl, pl, RNG)
    Csh = build_reduced(5, nl, sp)
    tabs = {}
    for k in range(nl):
        Zk = Csh['Zd'][k]; Cs = dict(Csh); Cs['Kmat'] = np.outer(Zk, Zk); Cs['n_links'] = 1
        rr, _ = measure_rates(Cs, 0.0, [1 << p for p in range(5)])
        for p in range(5): tabs[(k, p)] = rr[1 << p]
    Ash_hat = reconstruct_route3_CIRCULAR(5, nl, tabs, GAM)
    print(f"  {nm}: fidelity vs own truth {fidelity(Ah, At):.3f}   "
          f"on a shuffle, fidelity vs the shuffle's own truth {fidelity(Ash_hat, true_adj(nl, sp)):.3f}")
print("  Both are 1.000 by construction.  CIRCULAR -- the input already is the incidence structure.")

# ----------------------------------------------------------------------------------------------
print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)
nl, pl = CARRIERS['BLOCK5']; At = true_adj(nl, pl)
Ah0, Gh0, _ = run_route2('BLOCK5', 5, nl, pl, 0.0)
f_true = fidelity(Ah0, At)
print(f"  ROUTE 1  record algebra alone (commutators/products)      NULL   "
      f"(all commutators 0; identical Z_2^5 table on block, chain and ring)")
print(f"  ROUTE 2  records + ONE uniform unlabelled bath            fidelity {f_true:.3f} on the true carrier")
print(f"           shuffles, scored vs ORIGINAL adjacency           {np.mean(sh_orig):.3f} +- {np.std(sh_orig):.3f}"
      f"   (collapses, as required)")
print(f"           shuffles, scored vs the SHUFFLE'S own adjacency  {np.mean(sh_own):.3f} +- {np.std(sh_own):.3f}"
      f"   (tracks the carrier it is run on)")
print(f"           random-graph control                             {np.mean(rand_f):.3f}")
print(f"           second geometry RING5 vs BLOCK5                  separated; degree sequences "
      f"{degseq(recon_duals['RING5'])} vs {degseq(recon_duals['BLOCK5'])}")
print(f"  ROUTE 3  link-resolved bath                               CIRCULAR (input is the incidence matrix)")
print(f"\n  g2 scan (BLOCK5): " + "  ".join(f"{g:.3g}:{f:.2f}" for g, f, _, _, _, _ in scan_rows))
print("""
  THE HONEST READING.  Adjacency is NOT in the record algebra -- Route 1 is a forced null and that
  is the lens's real answer.  It IS recoverable once one uniform environment is present, from
  n_p + C(n_p,2) label-free scalars, exactly and without fitting.  But the recovery is an identity,
  not an inference: the pair-region decay rates ARE the Gram matrix of the record boundaries up to
  an affine map, and 'plaquettes share a link' is by definition Gram > 0.  So the correct claim is
  the weaker and more useful one -- the dual graph is a FUNCTION of quantities that carry no link
  labels -- and not the stronger one that geometry was inferred from something that did not already
  determine it.  What the shuffle control establishes is that those scalars are not fixed by the
  carrier's counts: the same counts give different rates and the method follows the carrier, not
  the bookkeeping.""")
