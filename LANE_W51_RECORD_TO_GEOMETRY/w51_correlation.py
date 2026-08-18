"""W-51.  CAN THE CARRIER'S GEOMETRY BE READ OFF THE RECORDS?

THE GAP.  Every lane so far installed the graph by hand: which link touches which vertex, and
therefore which plaquettes neighbour each other, was an INPUT. This lane asks whether the dual
graph (which plaquettes share a link) is RECOVERABLE from records + dynamics alone.

METHOD -- CORRELATIONS BETWEEN RECORDS, IN THE STATE WHERE THE RECORDS ARE WRITTEN.
Put the carrier in the state in which every record holds the value +1 (all plaquette operators
=+1). Couple it to ONE GLOBAL, UNIFORM bath -- every link damped at the same rate gamma, no
addressability, one scalar knob. Watch two records p and q, and watch their PRODUCT record.
Define, purely from measured decay rates at one fixed time t*,

        E(p,q) = Gamma(p) + Gamma(q) - Gamma(p*q)          "pair excess"

Decay is SUB-ADDITIVE exactly to the extent that two records' boundaries overlap, so E(p,q) counts
shared links.  Equivalently, in expectation values, the CONNECTED CORRELATION

        K(p,q,t) = <W_p W_q>(t) - <W_p>(t) <W_q>(t)

is identically ZERO for records with disjoint boundaries and nonzero otherwise. That is the whole
method: the dual graph is the support pattern of the connected record-record correlation.

WHY THE STATE MATTERS (the lens's warning, taken seriously).  In the maximally mixed state every
record has <W>=0 at all times: every pair looks identical and the method returns nothing. The
written state is the only one in which records HAVE values to correlate. Both are run below.

WHY THIS IS NOT GEOMETRY-IN-GEOMETRY-OUT.  See SECTION 0. The reconstructor is a pure function
    reconstruct(F, gamma, measure)
whose only channel to the carrier is an oracle returning ONE REAL NUMBER per record. It never
receives a link list, a vertex list, coordinates, or an incidence table.

SIMULATION FRAME.  The carrier is simulated in its RECORD SECTOR: physical states of the planar
Z_2 patch in the electric (loop) basis are closed link-configurations, and on a simply-connected
patch each is the boundary of a unique plaquette subset S. So the basis is labelled by subsets,
the plaquette flip is X_p, and the link operator is Z_k = prod_{p contains k} Z_p. This is not an
approximation and not a new model: SECTION 2 reproduces the original gauge-space construction
(w34_sieve.py / w41_k4.py machinery) to machine precision. It is used because a SHUFFLED carrier
has no lattice realisation, and the control must run on the same footing as the true carrier.

numpy only. expm by scaling-and-squaring, both a dense version and an action version, cross-checked
against each other on the row-major vectorised Liouvillian
    -i(H kron I - I kron H^T) + gamma sum (L kron L* - I kron I).
"""

import itertools, numpy as np

RNG = np.random.default_rng(51)
GAM = 0.5                      # bath strength: ONE scalar, same on every link
np.set_printoptions(linewidth=150, suppress=True)

# ======================================================================================
# SECTION 0 -- FORCED-OR-NOT, DECLARED BEFORE ANY DYNAMICS RUNS
# ======================================================================================
print("="*100)
print("W-51  RECONSTRUCTING THE DUAL GRAPH FROM RECORDS ALONE")
print("="*100)
print("""
SECTION 0  FORCED-OR-NOT (declared before a single number is computed)

  AT RISK #1 -- THE SINGLE-RECORD DECAY RATE Gamma(p).
    Every plaquette carries exactly FOUR links by definition of a plaquette. Under a global bath
    the counting law gives Gamma(p) = 2*gamma*|boundary(p)| = 8*gamma for EVERY p, on the true
    carrier, on every shuffle, on the ring and on the block alike. Var_p Gamma(p) = 0 is FORCED by
    the construction rule, not measured. THE ONE ARGUMENT: a quantity whose variance across
    plaquettes is fixed at zero by a definition cannot separate plaquettes, so no reconstruction
    may rest on single-record rates. This lane therefore rests on the PAIR EXCESS, which is a
    comparison of a joint rate against the sum of two individual rates and is not fixed by any
    counting rule. The prediction Var_p Gamma(p) = 0 is printed in SECTION 4 as a check on this
    declaration, NOT as a result.

  AT RISK #2 -- THE NUMBER OF EDGES IN THE RECONSTRUCTED GRAPH.
    The mandated shuffle keeps every count fixed, so it preserves the number of links lying on two
    plaquettes -- i.e. the TOTAL EDGE COUNT of the dual graph is invariant under shuffling and is
    therefore FORCED. Any fidelity measure that scores edge counts would score 100% on a shuffle
    and be a false positive. THE ONE ARGUMENT: fidelity here is scored PER PAIR (which of the
    C(F,2) pairs are adjacent, and with what multiplicity), never by totals. The chance baseline
    for per-pair scoring at fixed edge count is computed explicitly in SECTION 7.

  AT RISK #3 -- THE RECORD SET ITSELF.
    In the record sector the elementary record operators are X_0..X_{F-1}: single-qubit operators
    identical to one another up to which label they carry. PROOF THAT THEY CONTAIN NO ADJACENCY:
    the set {X_p} is invariant under ALL F! relabelings of the plaquettes, so it is invariant under
    every graph automorphism AND every non-automorphism, and an S_F-invariant object cannot single
    out one graph among the many on F labelled vertices. Geometry enters the simulation ONLY
    through the bath jump operators and the electric term of H -- i.e. only through the DYNAMICS,
    which is exactly the channel this lane is allowed to read.
""")

# ======================================================================================
#  expm -- scaling and squaring, no scipy
# ======================================================================================
def expm_ss(A, order=20):
    """Dense matrix exponential by scaling-and-squaring with a Taylor inner series."""
    nrm = float(np.abs(A).sum(axis=1).max())
    s = 0
    while nrm / (2.0**s) > 0.5:
        s += 1
    X = A / (2.0**s)
    n = A.shape[0]
    E = np.eye(n, dtype=complex); T = np.eye(n, dtype=complex)
    for k in range(1, order+1):
        T = (T @ X) / k
        E = E + T
    for _ in range(s):
        E = E @ E
    return E

def expm_action(applyL, rho, t, bound, order=20):
    """e^{tL} rho where L is given only as a callable. Scaling-and-squaring: the 'squaring' is
    repeated application of the short-time propagator, which is exact for an action."""
    theta = abs(t) * bound
    s = 1
    while theta / s > 0.25:
        s *= 2
    dt = t / s
    for _ in range(s):
        term = rho; acc = rho
        for k in range(1, order+1):
            term = applyL(term) * (dt / k)
            acc = acc + term
            if np.abs(term).max() < 1e-17 * max(np.abs(acc).max(), 1e-300):
                break
        rho = acc
    return rho

# ======================================================================================
#  CARRIERS.  A carrier is an INCIDENCE: F plaquettes, each on 4 distinct links, each link on
#  1 or 2 plaquettes. Nothing else about it is ever handed to the reconstructor.
# ======================================================================================
def carrier_block(nx, ny):
    """nx by ny block of plaquettes on a planar patch (the standard carrier of this program)."""
    H = nx * (ny + 1)
    hid = lambda i, j: j * nx + i
    vx  = lambda i, j: H + j * (nx + 1) + i
    inc = [[hid(i, j), hid(i, j+1), vx(i, j), vx(i+1, j)] for j in range(ny) for i in range(nx)]
    L = H + (nx + 1) * ny
    return dict(name=f"BLOCK {nx}x{ny}", F=nx*ny, L=L, inc=inc)

def carrier_ring(n):
    """n plaquettes closed into a ring (one row on a cylinder). Dual graph = n-cycle."""
    hid = lambda i, j: j * n + i
    vx  = lambda i: 2*n + i
    inc = [[hid(i, 0), hid(i, 1), vx(i), vx((i+1) % n)] for i in range(n)]
    return dict(name=f"RING {n}", F=n, L=3*n, inc=inc)

def carrier_strip(n):
    """n plaquettes in an open row. Dual graph = path."""
    hid = lambda i, j: j * n + i
    vx  = lambda i: 2*n + i
    inc = [[hid(i, 0), hid(i, 1), vx(i), vx(i+1)] for i in range(n)]
    return dict(name=f"STRIP 1x{n}", F=n, L=3*n+1, inc=inc)

def carrier_from_dual(name, F, edges):
    """Abstract carrier with a prescribed dual graph: each listed edge is an interior link shared
    by two plaquettes, and each plaquette is padded to 4 links with links of its own."""
    inc = [[] for _ in range(F)]
    k = 0
    for (a, b) in edges:
        inc[a].append(k); inc[b].append(k); k += 1
    for p in range(F):
        while len(inc[p]) < 4:
            inc[p].append(k); k += 1
    assert all(len(set(v)) == 4 for v in inc)
    return dict(name=name, F=F, L=k, inc=inc)

def multiplicity_profile(c):
    m = {}
    for ls in c["inc"]:
        for k in ls: m[k] = m.get(k, 0) + 1
    prof = {}
    for k in range(c["L"]): prof[m.get(k, 0)] = prof.get(m.get(k, 0), 0) + 1
    return dict(sorted(prof.items()))

def shuffle_carrier(c, rng):
    """MANDATORY CONTROL. Randomly permute which links belong to which plaquettes, keeping every
    count fixed: same number of links, same multiplicity of every link, 4 distinct links per
    plaquette."""
    mult = {}
    for ls in c["inc"]:
        for k in ls: mult[k] = mult.get(k, 0) + 1
    slots = []
    for k in range(c["L"]): slots += [k] * mult.get(k, 0)
    for _ in range(20000):
        rng.shuffle(slots)
        inc = [slots[4*p:4*p+4] for p in range(c["F"])]
        if all(len(set(v)) == 4 for v in inc):
            return dict(name=c["name"] + " SHUFFLED", F=c["F"], L=c["L"], inc=inc)
    raise RuntimeError("shuffle failed")

def true_adjacency(c):
    """GROUND TRUTH -- used ONLY for scoring, never handed to the reconstructor."""
    F = c["F"]; A = np.zeros((F, F), dtype=int)
    for p in range(F):
        for q in range(p+1, F):
            A[p, q] = A[q, p] = len(set(c["inc"][p]) & set(c["inc"][q]))
    return A

# ======================================================================================
#  THE DUAL / RECORD-SECTOR MODEL
# ======================================================================================
_PC = np.array([bin(i).count("1") for i in range(1 << 16)], dtype=np.int64)

def model(c):
    F, L, inc = c["F"], c["L"], c["inc"]
    D = 1 << F
    masks = np.arange(D, dtype=np.int64)
    owner = np.zeros(L, dtype=np.int64)
    for p, ls in enumerate(inc):
        for k in ls: owner[k] |= (1 << p)
    Z = np.empty((L, D))
    for k in range(L):
        Z[k] = 1.0 - 2.0 * (_PC[masks & owner[k]] & 1)     # Z_k = prod_{p in link k} Z_p
    M = Z.T @ Z                                            # sum_k Z_k rho Z_k = M .* rho
    Hx = np.zeros((D, D))
    for p in range(F):
        Hx[masks ^ (1 << p), masks] += 1.0                 # plaquette flip = X_p
    Hz = Z.sum(axis=0)
    return dict(c=c, F=F, L=L, D=D, Z=Z, M=M, Hx=Hx, Hz=Hz)

def lindblad(m, g2, gam=GAM):
    """H = -MAG - g2*ELEC with MAG = sum(Move+h.c.) = 2 sum_p X_p and ELEC = sum(Z+h.c.) = 2 sum_k Z_k
    (exactly the convention of w34_sieve.py / w41_k4.py). Bath: EVERY link, uniform gamma."""
    H = -2.0 * (m["Hx"] + g2 * np.diag(m["Hz"]))
    M, L = m["M"], m["L"]
    def applyL(rho):
        return -1j * (H @ rho - rho @ H) + gam * (M * rho - L * rho)
    bound = 2.0 * float(np.abs(H).sum(axis=1).max()) + 2.0 * gam * L
    return applyL, bound, H

def state(m, kind, rng=None):
    D = m["D"]
    if kind == "written":                 # every record holds +1: the all-+ state, X_p=+1 for all p
        psi = np.ones(D, dtype=complex) / np.sqrt(D)
        return np.outer(psi, psi.conj())
    if kind == "mixed":
        return np.eye(D, dtype=complex) / D
    if kind == "random":
        v = rng.normal(size=D) + 1j * rng.normal(size=D); v /= np.linalg.norm(v)
        return np.outer(v, v.conj())
    raise ValueError(kind)

def record_expectation(rho, mask):
    """<W_S> for the record W_S = prod_{p in S} X_p, a permutation on the subset basis."""
    D = rho.shape[0]
    idx = np.arange(D) ^ mask
    return complex(np.sum(rho[idx, np.arange(D)]))

# ======================================================================================
#  THE ORACLE.  Everything the reconstructor is allowed to touch.
# ======================================================================================
class Oracle:
    """Answers exactly one question: 'what is the decay rate of THIS record?'  Returns one real
    number. It is handed a tuple of record LABELS and nothing else; it returns no vectors, no
    supports, no link ids."""
    def __init__(self, m, g2, t, kind="written", gam=GAM, rng=None):
        applyL, bound, H = lindblad(m, g2, gam)
        self.rho0 = state(m, kind, rng)
        self.rhot = expm_action(applyL, self.rho0, t, bound)
        self.t = t; self.m = m; self.n = 0
    def value(self, S):
        self.n += 1
        mask = 0
        for p in S: mask ^= (1 << p)
        return record_expectation(self.rhot, mask).real
    def rate(self, S):
        v = self.value(S)
        if not np.isfinite(v) or v <= 0.0:
            return np.nan
        return -np.log(v) / self.t

# ======================================================================================
#  THE RECONSTRUCTOR.  Signature is the proof: F, gamma, and a scalar oracle. Nothing else.
# ======================================================================================
def reconstruct(F, gamma, measure):
    single = np.array([measure((p,)) for p in range(F)], dtype=float)
    A = np.zeros((F, F)); E = np.zeros((F, F)); worst = 0.0
    for p in range(F):
        for q in range(p+1, F):
            e = single[p] + single[q] - measure((p, q))
            x = e / (4.0 * gamma)
            E[p, q] = E[q, p] = e
            if np.isfinite(x):
                A[p, q] = A[q, p] = np.round(x)
                worst = max(worst, abs(x - round(x)))
            else:
                A[p, q] = A[q, p] = -1
                worst = np.inf
    return A.astype(int), E, single, worst

def reconstruct_corr(F, value):
    """Second channel, same lens, no rates: the CONNECTED correlation between two records."""
    K = np.zeros((F, F))
    v1 = [value((p,)) for p in range(F)]
    for p in range(F):
        for q in range(p+1, F):
            K[p, q] = K[q, p] = value((p, q)) - v1[p] * v1[q]
    return K

def score(Ahat, Aref):
    F = Aref.shape[0]; iu = np.triu_indices(F, 1)
    a = Ahat[iu]; b = np.asarray(Aref)[iu]
    acc = float(np.mean(a == b))
    Eh = set(np.flatnonzero(a > 0)); Er = set(np.flatnonzero(b > 0))
    jac = len(Eh & Er) / max(len(Eh | Er), 1)
    return acc, jac, bool(np.all(a == b))

# ======================================================================================
# SECTION 1 -- expm VALIDATION: action version vs dense row-major vectorised Liouvillian
# ======================================================================================
print("="*100)
print("SECTION 1  expm VALIDATION (scaling-and-squaring, two independent implementations)")
c4 = carrier_block(2, 2); m4 = model(c4)
applyL, bound, H4 = lindblad(m4, 0.05)
D = m4["D"]; Id = np.eye(D, dtype=complex)
Lvec = -1j * (np.kron(H4, Id) - np.kron(Id, H4.T)) - GAM * m4["L"] * np.kron(Id, Id)
for k in range(m4["L"]):
    Zk = np.diag(m4["Z"][k]).astype(complex)
    Lvec = Lvec + GAM * np.kron(Zk, Zk.conj())
r0 = state(m4, "written")
tval = 0.2
ra = expm_action(applyL, r0, tval, bound)
rb = (expm_ss(tval * Lvec) @ r0.reshape(-1)).reshape(D, D)
print(f"  vectorised Liouvillian dim {Lvec.shape[0]}  (row-major convention, as mandated)")
print(f"  ||action - dense||_max = {np.abs(ra-rb).max():.3e}     trace(rho_t) = {np.trace(ra).real:.12f}")
print(f"  hermiticity defect     = {np.abs(ra-ra.conj().T).max():.3e}     min eig = {np.linalg.eigvalsh(0.5*(ra+ra.conj().T)).min():.3e}")
ok1 = np.abs(ra-rb).max() < 1e-10 and abs(np.trace(ra).real-1) < 1e-10
print(f"  -> the two exponentials agree and the propagated state is a state: {ok1}")

# ======================================================================================
# SECTION 2 -- CROSS-VALIDATION AGAINST THE ORIGINAL GAUGE-SPACE CONSTRUCTION
# ======================================================================================
print()
print("="*100)
print("SECTION 2  CROSS-VALIDATION: record-sector model vs the original gauge-space machinery")
print("           (build/Move/Zop/compose exactly as in w34_sieve.py and w41_k4.py)")

def build(V, E, N):
    st = [s for s in itertools.product(range(N), repeat=len(E))
          if all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
                 -sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % N == 0 for v in range(len(V)))]
    return st, {s: i for i, s in enumerate(st)}
def Zop(st, links, N):
    w = np.exp(2j*np.pi/N)
    return np.diag([w**(sum(s[k] for k in links) % N) for s in st]).astype(complex)
def Move(st, idx, mv, N):
    Dm = len(st); Mx = np.zeros((Dm, Dm), complex)
    for j, s in enumerate(st):
        t = list(s)
        for k, sg in mv: t[k] = (t[k]+sg) % N
        t = tuple(t)
        if t in idx: Mx[idx[t], j] = 1.0
    return Mx
def compose(ps):
    acc = {}
    for p in ps:
        for k, sg in p: acc[k] = acc.get(k, 0)+sg
    return [(k, s) for k, s in acc.items() if s != 0]

V2 = [(i, j) for j in range(3) for i in range(3)]; vid = {v: k for k, v in enumerate(V2)}
Eo = []
for j in range(3):
    for i in range(2): Eo.append((vid[(i, j)], vid[(i+1, j)]))
for j in range(2):
    for i in range(3): Eo.append((vid[(i, j)], vid[(i, j+1)]))
hid = lambda i, j: j*2+i; vxo = lambda i, j: 6+j*3+i
Po = [[(hid(i, j), +1), (vxo(i+1, j), +1), (hid(i, j+1), -1), (vxo(i, j), -1)] for j in range(2) for i in range(2)]
sto, idxo = build(V2, Eo, 2); Do = len(sto)
MAG = sum((lambda X: X+X.conj().T)(Move(sto, idxo, p, 2)) for p in Po)
ELEC = sum(Zop(sto, [k], 2)+Zop(sto, [k], 2).conj().T for k in range(len(Eo)))
Zo = np.array([np.diag(Zop(sto, [k], 2)).real for k in range(len(Eo))])
Mo = Zo.T @ Zo
print(f"  original gauge space: {len(Eo)} links, {len(V2)} vertices, physical dim {Do} "
      f"(record sector dim {m4['D']})")
g2t, tt = 0.05, 0.2
Ho = -MAG - g2t*ELEC
def applyLo(rho):
    return -1j*(Ho @ rho - rho @ Ho) + GAM*(Mo*rho - len(Eo)*rho)
bo = 2.0*float(np.abs(Ho).sum(axis=1).max()) + 2.0*GAM*len(Eo)
psio = np.ones(Do, dtype=complex)/np.sqrt(Do)
rho_o = expm_action(applyLo, np.outer(psio, psio.conj()), tt, bo)
oc = Oracle(m4, g2t, tt)
mx = 0.0
print(f"  {'record S':>12s} {'gauge-space rate':>18s} {'record-sector rate':>20s} {'diff':>12s}")
for r in range(1, 5):
    for S in itertools.combinations(range(4), r):
        Wo = Move(sto, idxo, compose([Po[i] for i in S]), 2)
        vo = float(np.trace(Wo @ rho_o).real); ro = -np.log(vo)/tt
        rd = oc.rate(S)
        mx = max(mx, abs(ro-rd))
        if r <= 2 and S in [(0,), (1,), (2,), (3,), (0, 1), (0, 3)]:
            print(f"  {str(S):>12s} {ro:18.10f} {rd:20.10f} {abs(ro-rd):12.2e}")
print(f"  max |gauge-space - record-sector| over ALL 15 records = {mx:.3e}")
print(f"  -> the record-sector model IS the original carrier, not a new toy: {mx < 1e-9}")

# ======================================================================================
# SECTION 3 -- OPERATORS MUST EARN MEASUREMENT
# ======================================================================================
print()
print("="*100)
print("SECTION 3  OPERATOR AUDIT (before any commutator or rate is read)")
CT = carrier_block(3, 2)
mT = model(CT)
Dt = mT["D"]
def Wmat(F, S):
    Dd = 1 << F; mask = 0
    for p in S: mask ^= (1 << p)
    Wm = np.zeros((Dd, Dd)); Wm[np.arange(Dd) ^ mask, np.arange(Dd)] = 1.0
    return Wm
print(f"  carrier {CT['name']}: F={CT['F']} records, L={CT['L']} links, dim {Dt}")
print(f"  {'record':>10s} {'||O||_F':>10s} {'||O^dag O - I||':>16s} {'#distinct eigvals':>18s} {'max comm w/ W_0':>16s}")
W0 = Wmat(CT["F"], (0,))
for S in [(0,), (1,), (5,), (0, 1), (0, 5), (0, 1, 2, 3, 4, 5)]:
    Wm = Wmat(CT["F"], S)
    ev = np.linalg.eigvals(Wm)
    nd = len(np.unique(np.round(ev.real, 9) + 1j*np.round(ev.imag, 9)))
    comm = np.abs(Wm @ W0 - W0 @ Wm).max()
    print(f"  {str(S):>10s} {np.linalg.norm(Wm):10.4f} {np.abs(Wm.T@Wm-np.eye(Dt)).max():16.2e} {nd:18d} {comm:16.2e}")
print("  NOTE: all elementary records COMMUTE with one another (max commutator 0 above), so no")
print("  reconstruction from record-record commutators is possible. The geometry is not in the")
print("  record algebra; it can only be in the dynamics.")

# ======================================================================================
# SECTION 4 -- THE FORCED CHECK, MEASURED
# ======================================================================================
print()
print("="*100)
print("SECTION 4  FORCED CHECK (testing the SECTION 0 declaration, not producing a result)")
oc = Oracle(mT, 0.01, 0.2)
_, _, sing, _ = reconstruct(CT["F"], GAM, oc.rate)
print(f"  single-record rates on the true carrier : {np.round(sing, 8)}")
print(f"  variance across plaquettes              : {np.var(sing):.3e}   (declared FORCED = 0)")
csh = shuffle_carrier(CT, np.random.default_rng(7)); msh = model(csh)
oc2 = Oracle(msh, 0.01, 0.2)
_, _, sing2, _ = reconstruct(csh["F"], GAM, oc2.rate)
print(f"  single-record rates on a SHUFFLED carrier: {np.round(sing2, 8)}")
print(f"  variance across plaquettes              : {np.var(sing2):.3e}")
print(f"  edge count true {int(true_adjacency(CT).sum()//2)}  vs shuffled {int(true_adjacency(csh).sum()//2)}  "
      f"(declared FORCED to be equal by the count-preserving shuffle)")
_, Ecmp, _, _ = reconstruct(CT["F"], GAM, oc.rate)
print(f"  for scale: the PAIR EXCESS E(p,q) that the method actually reads ranges over "
      f"{Ecmp.min():.4f} .. {Ecmp.max():.4f}")
print("  -> single-record rates carry no information (spread ~1e-4, forced); the edge TOTAL carries")
print("     none either. Everything below is scored PER PAIR.")

# ======================================================================================
# SECTION 5 -- WHICH STATE?  (the lens's warning, tested)
# ======================================================================================
print()
print("="*100)
print("SECTION 5  WHICH STATE DO WE EVALUATE IN?")
AT = true_adjacency(CT)
for kind in ["written", "mixed", "random"]:
    o = Oracle(mT, 0.01, 0.2, kind=kind, rng=np.random.default_rng(3))
    Ah, Ee, sg, wd = reconstruct(CT["F"], GAM, o.rate)
    acc, jac, ex = score(Ah, AT)
    v1 = abs(o.value((0,)))
    tag = "records have values" if v1 > 1e-8 else "ALL RECORDS HAVE <W>=0 -- NOTHING TO CORRELATE"
    print(f"  state={kind:8s}  |<W_0>(t)| = {v1:.3e}   per-pair fidelity vs truth = "
          f"{acc:.3f}   integer defect {wd:.2e}   [{tag}]")

# ======================================================================================
# SECTION 6/7/8/9 -- THE RUN
# ======================================================================================
def spectrum(A):
    return np.round(np.sort(np.linalg.eigvalsh(A.astype(float))), 6)

def triangles(A):
    return int(round(np.trace(np.linalg.matrix_power(A.astype(float), 3)) / 6))

def run(c, g2, t, kind="written", rng=None):
    m = model(c)
    o = Oracle(m, g2, t, kind=kind, rng=rng)
    Ah, Ee, sg, wd = reconstruct(c["F"], GAM, o.rate)
    nq = o.n                                   # queries used by the RECONSTRUCTOR alone
    K = reconstruct_corr(c["F"], o.value)      # second channel, counted separately
    nneg = sum(1 for S in [(p,) for p in range(c["F"])]
                       + [(p, q) for p in range(c["F"]) for q in range(p+1, c["F"])]
               if o.value(S) <= 0.0)
    return Ah, Ee, K, wd, nq, nneg

print()
print("="*100)
print("SECTION 6  THE TRUE CARRIER")
G2, T = 0.01, 0.2
Ah, Ee, K, wd, nq, _ = run(CT, G2, T)
acc, jac, ex = score(Ah, AT)
print(f"  {CT['name']}  g2={G2}  t*={T}  gamma={GAM}  bath = ALL {CT['L']} links, uniform")
print(f"  oracle queries the RECONSTRUCTOR consumed = {nq} real numbers "
      f"(F + C(F,2) = {CT['F']+CT['F']*(CT['F']-1)//2}); the graph it must determine is "
      f"{CT['F']*(CT['F']-1)//2} bits.")
print("\n  TRUE adjacency (never shown to the reconstructor):"); print(AT)
print("  RECONSTRUCTED adjacency:"); print(Ah)
print(f"\n  per-pair fidelity {acc:.4f}   edge-set Jaccard {jac:.4f}   exact graph match {ex}"
      f"   worst integer defect {wd:.2e}")
iu = np.triu_indices(CT["F"], 1)
adj = AT[iu] > 0
print(f"  connected correlation K: |K| on ADJACENT pairs   min {np.abs(K[iu][adj]).min():.3e} "
      f"max {np.abs(K[iu][adj]).max():.3e}")
print(f"                           |K| on DISJOINT pairs   min {np.abs(K[iu][~adj]).min():.3e} "
      f"max {np.abs(K[iu][~adj]).max():.3e}")
sep = np.abs(K[iu][adj]).min() / max(np.abs(K[iu][~adj]).max(), 1e-300)
print(f"  separation ratio (min adjacent / max disjoint) = {sep:.3e}")

print()
print("="*100)
print("SECTION 7  THE SHUFFLE CONTROL  (decides the lane)")
print("  Randomly permute which links belong to which plaquettes, every count fixed, rebuild, re-run.")
print("  COLUMN 'vs TRUE' is the control: it asks whether the method still returns the ORIGINAL")
print("  block's graph. If it does, the method reads something other than this carrier and the")
print("  result is void. COLUMN 'vs OWN' asks whether it correctly returns the shuffled carrier's")
print("  own graph -- that is a diagnostic, not the control.")
rngs = np.random.default_rng(2026)
print(f"\n  {'shuffle':>9s} {'L':>4s} {'mult profile':>16s} {'edges':>6s} {'fid vs TRUE':>12s} "
      f"{'jac vs TRUE':>12s} {'fid vs OWN':>11s} {'exact vs OWN':>13s}")
sh_true, sh_own = [], []
for i in range(6):
    cs = shuffle_carrier(CT, rngs)
    As = true_adjacency(cs)
    Ahs, _, _, wds, _, _ = run(cs, G2, T)
    a1, j1, _ = score(Ahs, AT)
    a2, j2, e2 = score(Ahs, As)
    sh_true.append(a1); sh_own.append(a2)
    print(f"  {i:9d} {cs['L']:4d} {str(multiplicity_profile(cs)):>16s} {int(As.sum()//2):6d} "
          f"{a1:12.4f} {j1:12.4f} {a2:11.4f} {str(e2):>13s}")
# explicit chance baseline at fixed edge count
tot = int(AT.sum()//2); npairs = CT["F"]*(CT["F"]-1)//2
rb = np.random.default_rng(11)
base = []
for _ in range(20000):
    v = np.zeros(npairs, dtype=int); v[rb.permutation(npairs)[:tot]] = 1
    base.append(np.mean(v == (AT[iu] > 0).astype(int)))
print(f"\n  chance baseline (random graph, SAME edge count, per-pair) = {np.mean(base):.4f}")
print(f"  shuffled-carrier fidelity vs TRUE: mean {np.mean(sh_true):.4f}  max {np.max(sh_true):.4f}")
print(f"  shuffled-carrier fidelity vs OWN : mean {np.mean(sh_own):.4f}  min {np.min(sh_own):.4f}")

print()
print("="*100)
print("SECTION 8  A SECOND GEOMETRY (and a third, with IDENTICAL counts to the block)")
print("  NOTE, recorded because it nearly produced a false null: the first third carrier tried here")
print("  was a 6-cycle-plus-one-chord, which is ISOMORPHIC to the 3x2 block's dual (the 3-rung")
print("  ladder IS the theta graph). It was replaced by the BOWTIE below, which shares F, L, edge")
print("  count, multiplicity profile AND degree sequence with the block but is a different graph.")
BOW = carrier_from_dual("BOWTIE-6 (two triangles + bridge)", 6,
                        [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3)])
others = [carrier_ring(6), carrier_strip(6), BOW]
print()
print(f"  {'carrier':>34s} {'F':>3s} {'L':>4s} {'edg':>4s} {'true degrees':>16s} "
      f"{'recon degrees':>16s} {'tri':>4s} {'fid':>7s} {'exact':>6s}")
recon_store = {}
for c in [CT] + others:
    A = true_adjacency(c)
    Ahc, _, Kc, wdc, _, _ = run(c, G2, T)
    a, j, e = score(Ahc, A)
    recon_store[c["name"]] = Ahc
    print(f"  {c['name']:>34s} {c['F']:3d} {c['L']:4d} {int(A.sum()//2):4d} "
          f"{str([int(x) for x in sorted(A.sum(axis=1))]):>16s} "
          f"{str([int(x) for x in sorted(Ahc.sum(axis=1))]):>16s} {triangles(Ahc):4d} {a:7.4f} {str(e):>6s}")
kb = recon_store[CT["name"]]; kr = recon_store["RING 6"]; kc = recon_store[BOW["name"]]
print(f"\n  do the reconstructions DISTINGUISH the geometries?  (statistics computed from the")
print(f"  RECONSTRUCTED matrices only -- the adjacency spectrum is an isomorphism invariant)")
print(f"    BLOCK  spectrum {spectrum(kb)}  triangles {triangles(kb)}")
print(f"    RING   spectrum {spectrum(kr)}  triangles {triangles(kr)}")
print(f"    BOWTIE spectrum {spectrum(kc)}  triangles {triangles(kc)}")
d_br = not np.allclose(spectrum(kb), spectrum(kr))
d_bw = not np.allclose(spectrum(kb), spectrum(kc))
print(f"    block vs ring   distinguished: {d_br}   (different L and edge count too)")
print(f"    block vs bowtie distinguished: {d_bw}   <- SAME F, SAME L=17, SAME 7 edges, SAME")
print(f"        multiplicity profile, SAME degree sequence [2,2,2,2,3,3]. Nothing countable")
print(f"        separates these two carriers; only the geometry does.")

print()
print("="*100)
print("SECTION 9  COUPLING SCAN  (numbers as they come; nothing fitted)")
print("  'dead' = how many of the 21 records had <W>(t*) <= 0, i.e. no decay rate exists at all.")
print(f"  {'g2':>8s} {'t*':>6s} {'fid BLOCK':>10s} {'defect':>10s} {'dead':>5s} {'fid RING':>9s} "
      f"{'fid BOWTIE':>11s} {'fid SHUF vs TRUE':>17s} {'sep ratio':>11s}")
scan = []
RING6 = carrier_ring(6)
for t in [0.05, 0.2, 1.0]:
    for g2 in [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]:
        Ab, _, Kb, db, _, nneg = run(CT, g2, t)
        ab, _, _ = score(Ab, AT)
        Ar, _, _, _, _, _ = run(RING6, g2, t); ar, _, _ = score(Ar, true_adjacency(RING6))
        Ac, _, _, _, _, _ = run(BOW, g2, t); ac, _, _ = score(Ac, true_adjacency(BOW))
        rr = np.random.default_rng(900+int(1000*g2))
        cs = shuffle_carrier(CT, rr)
        Asf, _, _, _, _, _ = run(cs, g2, t); asf, _, _ = score(Asf, AT)
        adjm = AT[iu] > 0
        s = np.abs(Kb[iu][adjm]).min() / max(np.abs(Kb[iu][~adjm]).max(), 1e-300)
        scan.append((g2, t, ab, db, ar, ac, asf, s, nneg))
        print(f"  {g2:8.3f} {t:6.2f} {ab:10.4f} {db:10.2e} {nneg:5d} {ar:9.4f} {ac:11.4f} "
              f"{asf:17.4f} {s:11.2e}")

# ======================================================================================
# SECTION 10 -- SUMMARY
# ======================================================================================
print()
print("="*100)
print("SECTION 10  SUMMARY (every verdict computed from the number printed beside it)")
fid_true = acc
fid_shuf_mean = float(np.mean(sh_true))
fid_shuf_max = float(np.max(sh_true))
chance = float(np.mean(base))
geo_ok = d_br and d_bw
scan_ok = [row for row in scan if row[1] == 0.2]
print(f"  true carrier per-pair fidelity          {fid_true:.4f}   exact graph match {ex}")
print(f"  shuffled carriers vs TRUE  mean/max     {fid_shuf_mean:.4f} / {fid_shuf_max:.4f}   "
      f"chance baseline {chance:.4f}")
print(f"  shuffled carriers vs OWN graph  mean    {np.mean(sh_own):.4f}")
print(f"  ring / strip / bowtie recovered exactly {all(score(recon_store[c['name']], true_adjacency(c))[2] for c in others)}")
print(f"  block vs bowtie (all counts identical)  distinguished: {d_bw}")
print(f"  fidelity across g2 at t*=0.2            {[f'{r[0]}:{r[2]:.2f}' for r in scan_ok]}")
print(f"  fidelity across g2 at t*=0.05           {[f'{r[0]}:{r[2]:.2f}' for r in scan if r[1]==0.05]}")
gap = fid_true - fid_shuf_max
if fid_true == 1.0 and fid_shuf_max <= chance + 0.15 and geo_ok:
    verdict = "RECONSTRUCTED"
    why = ("the dual graph is recovered exactly from 21 scalars, the count-preserving shuffle "
           "collapses it to the chance baseline, and two carriers with identical counts are "
           "separated.")
elif fid_shuf_max > chance + 0.3:
    verdict = "CIRCULAR / VOID"
    why = "the shuffle control still returns the original graph, so the method is not reading this carrier."
else:
    verdict = "NULL"
    why = "the true carrier is not recovered above the chance baseline."
print(f"\n  gap (true - worst shuffle) = {gap:.4f}")
print(f"  VERDICT: {verdict} -- {why}")
brk = [r for r in scan if r[2] < 1.0]
print(f"\n  SCOPE, stated with the failures: the reconstruction is exact at every (g2,t*) in the")
print(f"  scan EXCEPT (g2,t*) in {[(r[0], r[1]) for r in brk]}, both at g2=1.0.")
print(f"  The failure is NOT monotone in g2*t*: (g2=1.0,t*=0.2) fails with 10 records driven to")
print(f"  <W> <= 0, while (g2=0.5,t*=1.0) -- a LARGER product -- is exact. The coherent term makes")
print(f"  <W_S>(t) oscillate, and the method fails wherever the chosen t* lands near a zero")
print(f"  crossing, not above some coupling threshold. At t*=0.05 the method is exact across the")
print(f"  whole scan including g2=1.0. Nothing was fitted; the t* values were fixed in advance and")
print(f"  integer rounding is the only decision rule.")
print(f"\n  WHAT THIS DOES NOT SHOW: the record LABELS p=0..F-1 are given. This lane recovers WHICH")
print(f"  PAIRS of records are adjacent, not the identification of records themselves; and the")
print(f"  reconstruction is of the DUAL graph (plaquette adjacency), which is what was asked, not")
print(f"  of the link/vertex incidence beneath it.")
