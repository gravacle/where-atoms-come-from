"""W-44.  A METRIC THAT DOES WORK.  LINK LENGTH AS A QUANTUM DEGREE OF FREEDOM.

THE GAP.  W-27..W-42 measured a carrier with NO metric: links with no length, no proper time, no
clock.  W-43b tried to add one and failed in a specific, instructive way: it gave a site a two-level
"mass", set the system Hamiltonian to ZERO, and coupled the environment DIRECTLY to that mass with
sz (x) sz.  The potential the mass would source never entered any evolution, and the environment
coupling was a projector onto the very operator being read.  Its redundancy number was therefore
guaranteed before the code ran.  This lane does not repeat either mistake:

  * the metric SOURCES the couplings of the gauge dynamics and the traversal of a probe;
  * NO environment operator, and NO bath jump operator, ever touches the metric register.
    Everything the environment learns about the metric has to arrive through the dynamics.

THE LENS: LINK LENGTH.  Each of a chosen set of links carries a two-level length variable
N_k in {0,1} ("short"/"long").  Length does work in three places:

  (M1) TRAVERSAL.  A probe hopping link k has amplitude  tau * (1 - beta*N_k).  A long link takes
       longer to cross.  This factor is REAL and identical in both directions.
  (M2) DWELL PHASE (redshift).  While the probe sits at site k it accrues phase at rate V*N_{link k},
       so the phase along a path is sum over the path of (local length) x (time spent there).
  (M3) CURVED WILSON ACTION.  On the patch, the plaquette coupling of the plaquette containing the
       length link is J -> J*(1 - beta*N), and the electric weight of that link is 1 -> (1 + beta*N).
       This is the metric entering the gauge Hamiltonian the way sqrt(g) enters F^2.

WHY LENGTH IS NOT A SECOND GAUGE FIELD -- and why it is implemented as (M1)+(M2) and not as a
hopping phase.  A length must accumulate with the SAME SIGN in both directions (a path traversed
backwards is just as long).  A Hermitian hop cannot do that: if the forward hop carries exp(i*phi*N)
then its adjoint necessarily carries exp(-i*phi*N), and the round-trip is a signed product, i.e. a
U(1) FLUX, not a length.  Direction-independent accumulation in a Hermitian theory has to come from
a real amplitude (M1) or from a potential accrued in TIME (M2).  That is exactly the structural
difference the lens asks for:  the gauge holonomy is a signed PRODUCT round a CLOSED loop and is
blind to path length; the length is an unsigned SUM along whatever path is actually walked, and each
link's length is by itself a local observable.

ROUTES (requirement 2).  Two different routes, stated up front:
  Parts 2 and 3 propagate an EXACT JOINT PURE STATE under expm(-iHT) -- no superoperator.
  Part 4 builds the FULL LINDBLAD SUPEROPERATOR on a dim-32 system (1024 x 1024) and takes its
  eigendecomposition.  Row-major convention:  -i(H (x) I - I (x) H^T) + gamma*sum(L (x) L* - I (x) I).
  expm is scaling-and-squaring, written here, numpy only, no scipy.
"""

import itertools, numpy as np

# ================================ numerics =====================================================
def expm(A):
    nr = np.linalg.norm(A, np.inf)
    k = max(0, int(np.ceil(np.log2(nr))) + 1) if nr > 0 else 0
    B = A / (2.0 ** k); X = np.eye(A.shape[0], dtype=complex); T = X.copy()
    for m in range(1, 80):
        T = T @ B / m; X = X + T
        if np.linalg.norm(T, np.inf) < 1e-18 * max(1.0, np.linalg.norm(X, np.inf)): break
    for _ in range(k): X = X @ X
    return X

def vn(r):
    ev = np.linalg.eigvalsh((r + r.conj().T) / 2); ev = ev[ev > 1e-12]
    return float(-(ev * np.log2(ev)).sum())

def shannon(p):
    p = np.asarray(p, float); p = p[p > 1e-12]
    return float(-(p * np.log2(p)).sum())

def audit(name, O):
    d = O.shape[0]
    fro = np.linalg.norm(O)
    sp  = float(np.linalg.norm(O, 2))
    ud  = np.linalg.norm(O @ O.conj().T - np.eye(d))
    ne  = len(np.unique(np.round(np.linalg.eigvals(O), 6)))
    print(f"    {name:<34s} ||O||_2={sp:8.4f}  ||O||_F={fro:9.4f}  "
          f"unitarity defect={ud:8.2e}  distinct eigenvalues={ne}")
    return ne

def bit(e, b): return (e >> b) & 1

def kron_list(ops):
    o = np.array([[1]], complex)
    for x in ops: o = np.kron(o, x)
    return o

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], complex)
sz = np.array([[1, 0], [0, -1]], complex)

# ================================ the carrier ==================================================
V2 = [(i, j) for j in range(3) for i in range(3)]; vid = {v: k for k, v in enumerate(V2)}
E = []
for j in range(3):
    for i in range(2): E.append((vid[(i, j)], vid[(i + 1, j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i, j)], vid[(i, j + 1)]))
NL = len(E); NV = 9; CENTER = vid[(1, 1)]
CUT   = [k for k, (a, b) in enumerate(E) if a == CENTER or b == CENTER]
PERIM = [k for k in range(NL) if k not in CUT]
hid = lambda i, j: j * 2 + i; vx = lambda i, j: 6 + j * 3 + i
PLQ = [[(hid(i, j), +1), (vx(i + 1, j), +1), (hid(i, j + 1), -1), (vx(i, j), -1)]
       for j in range(2) for i in range(2)]

RIMV = [vid[(0,0)], vid[(1,0)], vid[(2,0)], vid[(2,1)], vid[(2,2)], vid[(1,2)], vid[(0,2)], vid[(0,1)]]
def link_between(u, v):
    for k, (a, b) in enumerate(E):
        if (a, b) == (u, v) or (b, a) == (u, v): return k
    raise ValueError((u, v))
RIML = [link_between(RIMV[k], RIMV[(k + 1) % 8]) for k in range(8)]
assert sorted(RIML) == sorted(PERIM)

def div(s, v):
    return (sum(s[k] for k, (a, b) in enumerate(E) if a == v)
          - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % 2

# --- pure gauge patch (no matter): the W-33..W-36 / W-40 carrier, physical dim 16 --------------
ST = [s for s in itertools.product(range(2), repeat=NL) if all(div(s, v) == 0 for v in range(NV))]
IX = {s: i for i, s in enumerate(ST)}; DS = len(ST)

def Move(mv, st=None, ix=None):
    st = ST if st is None else st; ix = IX if ix is None else ix
    D = len(st); M = np.zeros((D, D), complex)
    for j, s in enumerate(st):
        t = list(s)
        for k, sg in mv: t[k] = (t[k] + sg) % 2
        t = tuple(t)
        if t in ix: M[ix[t], j] = 1.0
    return M
def Zop(links, st=None):
    st = ST if st is None else st
    return np.diag([(-1.0) ** (sum(s[k] for k in links) % 2) for s in st]).astype(complex)
def compose(ps):
    acc = {}
    for p in ps:
        for k, sg in p: acc[k] = acc.get(k, 0) + sg
    return [(k, s) for k, s in acc.items() if s % 2]

BP  = [Move(p) for p in PLQ]
RG  = Move(compose(PLQ))              # rim Wilson loop = product of all four plaquettes
ZL  = [Zop([k]) for k in range(NL)]
IdS = np.eye(DS, dtype=complex)

# --- ring with explicit matter: probe on 8 rim vertices + static anti-charge at the centre -----
RSTATES = []
for p in range(8):
    for s in itertools.product(range(2), repeat=NL):
        q = lambda v: (1 if v == RIMV[p] else 0) ^ (1 if v == CENTER else 0)
        if all(div(s, v) == q(v) for v in range(NV)): RSTATES.append((p, s))
RIDX = {x: i for i, x in enumerate(RSTATES)}; DR = len(RSTATES)

print("=" * 100)
print("W-44   LINK LENGTH: A METRIC DEGREE OF FREEDOM THAT ENTERS THE DYNAMICS")
print("=" * 100)
print(f"  pure-gauge patch      physical dim {DS}")
print(f"  ring with matter      physical dim {DR}   ({DR//8} gauge states per probe position)")
print(f"  cut links  {CUT}   rim links {RIML}")

# ================================ PART 0: FORCED-OR-NOT, BEFORE ANY DYNAMICS ===================
print()
print("-" * 100)
print("PART 0.  FORCED-OR-NOT.  Named BEFORE any dynamics runs.")
print("-" * 100)
print("""
  AT RISK:  the predictability-sieve result of Part 4, i.e. "what does the dynamics protect once
  the metric is switched on".  Two of the three things it could report are already determined by
  commutators, and only the third is a measurement.

  (i) rate(R) IS FORCED.  R = B1 B2 B3 B4 is a product of plaquette flips supported on PERIM.
      [R, B_p] = 0 for every p, and R is diagonal-free but commutes with every LENGTH operator
      because the length register is a different tensor factor.  Therefore ANY metric-dependence of
      the magnetic couplings J_p(N) leaves [R, H_mag] = 0 exactly.  The bath is L_k = Z_k on the
      four CUT links, and support(R) INTERSECT CUT = empty, so [L_k, R] = 0 for every bath link:
      boundary count = 0.  R decays only through the electric term on PERIM, at O(g^4).
      Conclusion: R staying slow is NOT evidence about the metric.  It cannot come out otherwise.

  (ii) rate(LENGTH) AT eps = 0 IS FORCED.  N (and hence Lambda = sigma_z on the length qubit) is
      diagonal in the length basis and commutes with every B_p, every Z_k and every bath L_k.  If
      the length register has no transverse field then [Lambda, H] = 0 EXACTLY and its decay rate is
      identically zero, so it wins the sieve by a conservation law, not by selection.  Any sieve run
      at eps = 0 is therefore uninformative and is reported only as a baseline.

  (iii) WHAT IS ACTUALLY MEASURED: with eps > 0 (the length genuinely fluctuates) and beta > 0 (the
      length genuinely sources the gauge couplings), the length sector acquires a decay rate through
      its entanglement with the gauge field, which the bath then damps.  Whether that rate lands
      ABOVE or BELOW rate(R) is not fixed by any commutator and is the measurement.

  Same argument for Part 3:  Lambda commutes with every environment coupling operator
  (kappa * Z_cut (x) sigma_z), so at beta = 0 the joint state factorises and I(Lambda:F) = 0
  identically.  Every bit the environment holds about the metric has to be sourced by beta.
""")

# ================================ PART 1: OPERATOR AUDIT =======================================
print("-" * 100)
print("PART 1.  OPERATORS MUST EARN MEASUREMENT.")
print("-" * 100)
LENP = PERIM[0]                    # the patch's length link: a RIM link (so it sits in R's support)
LENC = CUT[0]                      # variant: length on a cut link, next to the environment
Rp = np.kron(RG, I2); Lp = np.kron(IdS, sz)
audit("R  (rim Wilson loop) (x) I_len", Rp)
audit("Lambda (link length) I_g (x) sz", Lp)
audit("R * Lambda (product)", Rp @ Lp)
print(f"    [R, Lambda] = {np.linalg.norm(Rp@Lp - Lp@Rp):.2e}   (they are independent records)")
MAGH = sum(B + B.conj().T for B in BP)
print(f"    ||[R, H_mag]|| = {np.linalg.norm(RG@MAGH - MAGH@RG):.2e}    "
      f"||[R, Z_cut0]|| = {np.linalg.norm(RG@ZL[CUT[0]] - ZL[CUT[0]]@RG):.2e}    "
      f"||[R, Z_perim0]|| = {np.linalg.norm(RG@ZL[PERIM[0]] - ZL[PERIM[0]]@RG):.2e}")

# ================================ PART 2: TRANSPORT ============================================
# length qubits live on RIML[0] (always traversed) and RIML[4] (the hop removed by the cut control)
RLEN = [RIML[0], RIML[4]]; nLr = len(RLEN); ER = 2 ** nLr; DTR = DR * ER
LIVE, DEAD = 0, 1

def hop_block(k, tau):
    M = np.zeros((DR, DR), complex); kn = (k + 1) % 8
    for j, (p, s) in enumerate(RSTATES):
        if p != k: continue
        t = list(s); t[RIML[k]] ^= 1; t = tuple(t)
        i = RIDX.get((kn, t))
        if i is not None: M[i, j] -= tau
    return M + M.conj().T

def site_proj(k):
    return np.diag([1.0 if p == k else 0.0 for p, _ in RSTATES]).astype(complex)

def len_factor(link, beta):
    """diag over the 2^nLr length configs of (1 - beta*N_link); identity if link carries no length"""
    if link not in RLEN: return np.eye(ER, dtype=complex)
    b = RLEN.index(link)
    return np.diag([1.0 - beta * bit(e, b) for e in range(ER)]).astype(complex)

def len_number(link):
    if link not in RLEN: return np.zeros((ER, ER), dtype=complex)
    b = RLEN.index(link)
    return np.diag([float(bit(e, b)) for e in range(ER)]).astype(complex)

def H_ring(tau=1.0, beta=0.0, V=0.0, skip=None):
    H = np.zeros((DTR, DTR), complex)
    for k in range(8):
        if skip is not None and k == skip: continue
        H = H + np.kron(hop_block(k, tau), len_factor(RIML[k], beta))     # M1 traversal amplitude
        if V != 0.0:
            H = H + V * np.kron(site_proj(k), len_number(RIML[k]))        # M2 dwell / redshift phase
    return H

Rring = np.zeros((DR, DR), complex)
for j, (p, s) in enumerate(RSTATES):
    t = list(s)
    for k in PERIM: t[k] ^= 1
    Rring[RIDX[(p, tuple(t))], j] = 1.0
Rfull = np.kron(Rring, np.eye(ER, dtype=complex))
def Lam_ring(b):
    return np.kron(np.eye(DR, dtype=complex),
                   np.diag([1.0 - 2.0 * bit(e, b) for e in range(ER)]).astype(complex))
Llive, Ldead = Lam_ring(LIVE), Lam_ring(DEAD)

print()
print("-" * 100)
print("PART 2.  TRANSPORT.  Does a probe read the LENGTHS it crosses, and does it read them")
print("         differently from the way it reads the gauge holonomy?")
print("-" * 100)
audit("R  (ring rim loop) (x) I_len", Rfull)
audit("Lambda_live  (length of RIML[0])", Llive)
audit("Lambda_dead  (length of RIML[4])", Ldead)
print(f"    total transport dimension = {DTR}")

POS = np.array([p for p, _ in RSTATES])
def probe_dist(v):
    """marginal over the 8 probe sites; the length register and the gauge field are traced out"""
    d = np.zeros(8)
    np.add.at(d, POS, (np.abs(v.reshape(DR, ER)) ** 2).sum(axis=1))
    return d

# structural check: probe-position coherences vanish because the Gauss sector moves with the probe
_chk = np.zeros((8, 8), complex)
_v = np.ones(DR, complex) / np.sqrt(DR)
for i, (p, s) in enumerate(RSTATES):
    for j2, (q, t) in enumerate(RSTATES):
        if s == t: _chk[p, q] += _v[i] * np.conj(_v[j2])
print(f"    probe RDM off-diagonal norm on a uniform state = "
      f"{np.linalg.norm(_chk - np.diag(np.diag(_chk))):.2e}   "
      "(zero: different probe sites live in different Gauss sectors, so the RDM is diagonal)")

def info_ring(psi, Op):
    br = []
    for sgn in (+1, -1):
        Pr = (np.eye(DTR) + sgn * Op) / 2
        v = Pr @ psi; pr = float(np.vdot(v, v).real)
        br.append((pr, probe_dist(v / np.sqrt(pr)) if pr > 1e-14 else None))
    avg = sum(p * d for p, d in br if d is not None)
    return shannon(avg) - sum(p * shannon(d) for p, d in br if d is not None)

def init_ring(seed=5):
    g = np.random.default_rng(seed)
    mask = np.array([1.0 if p == 0 else 0.0 for p, _ in RSTATES])
    w = (g.normal(size=DR) + 1j * g.normal(size=DR)) * mask
    Pp = (np.eye(DR) + Rring) / 2; Pm = (np.eye(DR) - Rring) / 2
    a = Pp @ w; b = Pm @ w; a /= np.linalg.norm(a); b /= np.linalg.norm(b)
    v = a + b; v /= np.linalg.norm(v)
    plus = np.ones(ER, complex) / np.sqrt(ER)          # every length qubit in |+>: <Lambda> = 0
    psi = np.kron(v, plus)
    return psi / np.linalg.norm(psi)

psi0 = init_ring()
print(f"    initial state: probe at site 0,  <R> = {np.vdot(psi0, Rfull@psi0).real:+.2e},  "
      f"<Lambda_live> = {np.vdot(psi0, Llive@psi0).real:+.2e},  "
      f"<Lambda_dead> = {np.vdot(psi0, Ldead@psi0).real:+.2e}")

TAU, BETA_R, VDW = 1.0, 0.6, 0.5
TIMES = [0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0]
DENSE = np.arange(0.0, 24.01, 0.25)

def propagator_set(H):
    """H is Hermitian; diagonalise once and propagate psi0 on a dense grid.
       Verified against the hand-written expm below."""
    Ev, W = np.linalg.eigh((H + H.conj().T) / 2)
    c = W.conj().T @ psi0
    return lambda T: W @ (np.exp(-1j * Ev * T) * c)

def transport_table(tag, beta, V, skip=None):
    H = H_ring(TAU, beta, V, skip)
    prop = propagator_set(H)
    chk = np.linalg.norm(prop(3.0) - expm(-1j * H * 3.0) @ psi0)
    print(f"\n  {tag}")
    print(f"    eigh-route vs hand-written expm at T=3: ||difference|| = {chk:.2e}")
    print(f"    {'T':>6s} {'I(R:probe)':>12s} {'I(Lam_live:probe)':>18s} {'I(Lam_dead:probe)':>18s} {'P(back at 0)':>13s}")
    print("    " + "-" * 72)
    out = {}
    for T in TIMES:
        psi = prop(T)
        d = probe_dist(psi)
        iR, iL, iD = info_ring(psi, Rfull), info_ring(psi, Llive), info_ring(psi, Ldead)
        out[T] = (iR, iL, iD)
        print(f"    {T:6.2f} {iR:12.6f} {iL:18.6f} {iD:18.6f} {d[0]:13.4f}")
    dense = np.array([[info_ring(prop(T), O) for O in (Rfull, Llive, Ldead)] for T in DENSE])
    print(f"    dense grid, {len(DENSE)} times in [0,24]:   "
          f"max  I(R)={dense[:,0].max():.6f}  I(Lam_live)={dense[:,1].max():.6f}  I(Lam_dead)={dense[:,2].max():.6f}")
    print(f"    {'':>34s}   mean I(R)={dense[:,0].mean():.6f}  I(Lam_live)={dense[:,1].mean():.6f}  "
          f"I(Lam_dead)={dense[:,2].mean():.6f}")
    out['dense'] = dense
    return out

full_on  = transport_table(f"METRIC ON   beta={BETA_R} (traversal), V={VDW} (dwell phase), tau={TAU}, full ring",
                           BETA_R, VDW)
full_off = transport_table("CONTROL  beta=0, V=0  -- the no-metric carrier.  Lambda columns must be exactly 0",
                           0.0, 0.0)

print("\n  EXACTNESS CHECK (requirement 6): with the metric coupling off, the gauge column must")
print("  reproduce the metric-free ring bit for bit.")
Hno = np.zeros((DR, DR), complex)
for k in range(8): Hno = Hno + hop_block(k, TAU)
v0 = init_ring()[::ER] * np.sqrt(ER)   # the ring factor of psi0, renormalised
v0 = v0 / np.linalg.norm(v0)
def info_bare(v):
    br = []
    for sgn in (+1, -1):
        Pr = (np.eye(DR) + sgn * Rring) / 2
        w = Pr @ v; pr = float(np.vdot(w, w).real)
        d = np.zeros(8)
        if pr > 1e-14: np.add.at(d, POS, np.abs(w / np.sqrt(pr)) ** 2)
        br.append((pr, d if pr > 1e-14 else None))
    avg = sum(p * d for p, d in br if d is not None)
    return shannon(avg) - sum(p * shannon(d) for p, d in br if d is not None)
print(f"    {'T':>6s} {'I(R:probe) no-metric carrier':>30s} {'I(R:probe) beta=V=0, dim 512':>30s} {'|diff|':>10s}")
print("    " + "-" * 80)
for T in [2.0, 4.0, 8.0, 20.0]:
    a = info_bare(expm(-1j * Hno * T) @ v0); b = full_off[T][0]
    print(f"    {T:6.2f} {a:30.9f} {b:30.9f} {abs(a-b):10.2e}")

print("\n  CONTROL  tau = 0 (no transport at all).  Everything must be 0.")
H0 = H_ring(0.0, BETA_R, VDW)
for T in [4.0, 20.0]:
    psi = expm(-1j * H0 * T) @ psi0
    print(f"    T={T:5.1f}   I(R)={info_ring(psi,Rfull):.3e}   I(Lam_live)={info_ring(psi,Llive):.3e}"
          f"   I(Lam_dead)={info_ring(psi,Ldead):.3e}")

print("\n  THE DECISIVE STRUCTURAL CONTROL -- CUT THE RING (remove hop 4).  Every local interaction")
print("  is kept; only the closed path is removed.  The gauge holonomy is a product round a CLOSED")
print("  loop; a length is a sum along whatever path is walked.  Dwell phase off (V=0) so that the")
print("  length of the untraversed link RIML[4] can enter ONLY by being crossed.")
cut_on = transport_table("CUT RING, metric ON (beta=%.2f, V=0)" % BETA_R, BETA_R, 0.0, skip=4)
print("\n    -> I(Lam_dead : probe) on the cut ring is the mechanism control: RIML[4] is never")
print("       traversed, so if the length acts only through traversal this must be exactly 0.")

print("\n  SUMMARY OF PART 2 (dense grid, 97 times in [0,24]):")
def mx(tab, i): return float(tab['dense'][:, i].max())
def mn(tab, i): return float(tab['dense'][:, i].mean())
print(f"    {'run':>34s} {'max I(R)':>10s} {'max I(Lam_live)':>16s} {'max I(Lam_dead)':>16s}"
      f" {'mean I(R)':>11s} {'mean I(Lam_live)':>17s}")
print("    " + "-" * 112)
for nm, tab in [("full ring, metric ON", full_on), ("full ring, metric OFF", full_off),
                ("CUT ring, metric ON", cut_on)]:
    print(f"    {nm:>34s} {mx(tab,0):10.6f} {mx(tab,1):16.6f} {mx(tab,2):16.6f}"
          f" {mn(tab,0):11.6f} {mn(tab,1):17.6f}")

# ================================ PART 3: REDUNDANCY ===========================================
print()
print("-" * 100)
print("PART 3.  REDUNDANCY.  How much does ONE fragment of a local environment learn about the")
print("         METRIC, versus about the GAUGE LOOP?  Same machinery, same state, same run.")
print("-" * 100)
print("  System = 3x3 patch (dim 16) (x) ONE link-length qubit  ->  dim 32.")
print("  Environment = 5 qubits, each coupled to  Z on a CUT link  --  NEVER to the length register.")
print("  The metric enters through the CURVED WILSON ACTION only:  J_p -> J_p*(1-beta*N) for the")
print("  plaquette holding the length link, and its electric weight 1 -> (1+beta*N).")

NQ = 5; DE = 2 ** NQ
def H_patch(beta, eps, g2, lenlink):
    Nn = np.diag([0.0, 1.0]).astype(complex)
    H = np.zeros((DS * 2, DS * 2), complex)
    for ip, p in enumerate(PLQ):
        fac = (I2 - beta * Nn) if lenlink in [k for k, _ in p] else I2
        Bh = BP[ip] + BP[ip].conj().T
        H = H - np.kron(Bh, fac)
    for k in range(NL):
        w = (I2 + beta * Nn) if k == lenlink else I2
        H = H - g2 * np.kron(ZL[k] + ZL[k].conj().T, w)
    H = H - eps * np.kron(IdS, sx)
    return H

def run_patch(beta, eps, g2, lenlink, kappa, T, seed=3):
    Hs = H_patch(beta, eps, g2, lenlink)
    DSY = DS * 2
    H = np.kron(Hs, np.eye(DE, dtype=complex))
    for k in range(NQ):
        Zk = np.kron(Zop([CUT[k % len(CUT)]]), I2)
        ops = [I2] * NQ; ops[k] = sz
        H = H + kappa * np.kron(Zk, kron_list(ops))
    U = expm(-1j * H * T)
    g = np.random.default_rng(seed)
    w = g.normal(size=DS) + 1j * g.normal(size=DS)
    a = (IdS + RG) / 2 @ w; b = (IdS - RG) / 2 @ w
    a /= np.linalg.norm(a); b /= np.linalg.norm(b)
    psiG = (a + b) / np.sqrt(2.0); psiG /= np.linalg.norm(psiG)      # <R> = 0
    psiS = np.kron(psiG, np.ones(2, complex) / np.sqrt(2.0))         # <Lambda> = 0
    plus = np.ones(2, complex) / np.sqrt(2.0)
    psiE = kron_list([plus.reshape(2, 1)] * NQ).reshape(-1)
    psi = U @ np.kron(psiS, psiE)
    return psi.reshape((DSY,) + (2,) * NQ), DSY

def holevo(psiT, DSY, Op, frag):
    br = []
    for sgn in (+1, -1):
        Pr = (np.eye(DSY) + sgn * Op) / 2
        v = np.tensordot(Pr, psiT, axes=([1], [0])); p = float(np.vdot(v, v).real)
        if p < 1e-14: br.append((0.0, None)); continue
        v = v / np.sqrt(p)
        keep = [0] + [1 + i for i in frag]; tr = [ax for ax in range(1 + NQ) if ax not in keep]
        M = np.transpose(v, keep + tr).reshape(DSY * 2 ** len(frag), -1)
        rho = (M @ M.conj().T).reshape(DSY, 2 ** len(frag), DSY, 2 ** len(frag))
        br.append((p, np.einsum('ijik->jk', rho)))
    avg = sum(p * m for p, m in br if m is not None)
    return vn(avg) - sum(p * vn(m) for p, m in br if m is not None)

def profile(psiT, DSY, Op):
    vals = []
    for f in range(NQ + 1):
        combos = list(itertools.combinations(range(NQ), f))[:20]
        vals.append(float(np.mean([holevo(psiT, DSY, Op, c) for c in combos])))
    return vals

Rsys = np.kron(RG, I2); Lsys = np.kron(IdS, sz)
G2, EPS = 0.01, 0.35

print("\n  SATURATION FIRST.  A fragment plot means nothing until the WHOLE environment holds")
print("  something.  Scan the coupling; read I(.:all 5) for BOTH records.  Operating-point rule,")
print("  fixed before the scan: take the grid point maximising min(I(R:all), I(Lambda:all)), so")
print("  neither record is read at a setting that starves it.  The max-I(R:all) point is also shown.")
print(f"    {'kappa':>7s} {'T':>6s} {'I(R:all)':>11s} {'I(Lambda:all)':>15s} {'min of the two':>16s}")
print("    " + "-" * 62)
bestmin = None; bestR = None
for kap, T in [(0.5, 6.0), (1.0, 6.0), (1.0, 12.0), (2.0, 8.0), (3.0, 12.0),
               (5.0, 12.0), (8.0, 16.0), (12.0, 16.0)]:
    pT, DSY = run_patch(0.8, EPS, G2, LENP, kap, T)
    iR = holevo(pT, DSY, Rsys, tuple(range(NQ))); iL = holevo(pT, DSY, Lsys, tuple(range(NQ)))
    print(f"    {kap:7.2f} {T:6.1f} {iR:11.6f} {iL:15.6f} {min(iR,iL):16.6f}")
    if bestmin is None or min(iR, iL) > bestmin[2]: bestmin = (kap, T, min(iR, iL))
    if bestR is None or iR > bestR[2]: bestR = (kap, T, iR)
KAP, TT = bestmin[0], bestmin[1]
KAPR, TTR = bestR[0], bestR[1]
print(f"    -> operating point kappa={KAP}, T={TT}   (max-I(R:all) point was kappa={KAPR}, T={TTR})")

def redundancy_block(beta, lenlink, tag, eps=EPS, kap=None, T=None):
    kap = KAP if kap is None else kap; T = TT if T is None else T
    pT, DSY = run_patch(beta, eps, G2, lenlink, kap, T)
    vR = profile(pT, DSY, Rsys); vL = profile(pT, DSY, Lsys)
    print(f"\n  {tag}")
    print("    |F|          : " + "  ".join(f"{i:8d}" for i in range(NQ + 1)))
    print("    I(R:F)       : " + "  ".join(f"{v:8.4f}" for v in vR))
    print("    I(Lambda:F)  : " + "  ".join(f"{v:8.4f}" for v in vL))
    rR = vR[1] / vR[-1] if vR[-1] > 1e-9 else float('nan')
    rL = vL[1] / vL[-1] if vL[-1] > 1e-9 else float('nan')
    sR = f"{rR:.4f}" if vR[-1] > 1e-9 else "n/a (I(all)=0)"
    sL = f"{rL:.4f}" if vL[-1] > 1e-9 else "n/a (I(all)=0)"
    print(f"    ratio I(|F|=1)/I(all):   gauge {sR}     metric {sL}")
    return vR, vL, rR, rL

vR_on, vL_on, rR_on, rL_on = redundancy_block(
    0.8, LENP, f"METRIC ON  beta=0.8, eps={EPS}, length on RIM link {LENP} (inside R's support)")
vR_c, vL_c, rR_c, rL_c = redundancy_block(
    0.8, LENC, f"METRIC ON  beta=0.8, eps={EPS}, length on CUT link {LENC} (right next to the bath)")
vR_off, vL_off, rR_off, rL_off = redundancy_block(
    0.0, LENP, "CONTROL beta=0 -- the metric sources nothing.  I(Lambda:F) must be exactly 0 at every |F|;")
print(f"    CONTROL FIRES: max_f |I(Lambda:F)| at beta=0 = {max(abs(x) for x in vL_off):.3e}")
print("\n  EXACTNESS (requirement 6): at beta=0 the length register decouples, so the GAUGE profile")
print("  must equal the metric-free dim-16 carrier's profile under the identical environment.")

def run_bare(kappa, T, seed=3):
    H = np.kron(-MAGH - G2 * sum(ZL[k] + ZL[k].conj().T for k in range(NL)), np.eye(DE, dtype=complex))
    for k in range(NQ):
        ops = [I2] * NQ; ops[k] = sz
        H = H + kappa * np.kron(Zop([CUT[k % len(CUT)]]), kron_list(ops))
    U = expm(-1j * H * T)
    g = np.random.default_rng(seed)
    w = g.normal(size=DS) + 1j * g.normal(size=DS)
    a = (IdS + RG) / 2 @ w; b = (IdS - RG) / 2 @ w
    a /= np.linalg.norm(a); b /= np.linalg.norm(b)
    psiG = (a + b) / np.sqrt(2.0); psiG /= np.linalg.norm(psiG)
    plus = np.ones(2, complex) / np.sqrt(2.0)
    psi = U @ np.kron(psiG, kron_list([plus.reshape(2, 1)] * NQ).reshape(-1))
    return psi.reshape((DS,) + (2,) * NQ)

pB = run_bare(KAP, TT)
_hold = NQ
vB = [float(np.mean([holevo(pB, DS, RG, c) for c in list(itertools.combinations(range(NQ), f))[:20]]))
      for f in range(NQ + 1)]
print("    I(R:F) metric-free dim-16 : " + "  ".join(f"{v:8.4f}" for v in vB))
print("    I(R:F) beta=0    dim-32   : " + "  ".join(f"{v:8.4f}" for v in vR_off))
print(f"    max |difference| = {max(abs(a-b) for a, b in zip(vB, vR_off)):.3e}")
print("\n  THE METRIC DOES WORK, in this measurement:  switching beta on CHANGES the gauge numbers.")
print(f"    max_f |I(R:F)[beta=0.8, len on RIM] - I(R:F)[beta=0]| = "
      f"{max(abs(a-b) for a, b in zip(vR_on, vR_off)):.3e}")
print(f"    max_f |I(R:F)[beta=0.8, len on CUT] - I(R:F)[beta=0]| = "
      f"{max(abs(a-b) for a, b in zip(vR_c, vR_off)):.3e}")

print("\n  CONTROL kappa = 0 -- no environment coupling at all.  Both records must be 0 everywhere.")
pT0, DSY0 = run_patch(0.8, EPS, G2, LENP, 0.0, TT)
v0R = profile(pT0, DSY0, Rsys); v0L = profile(pT0, DSY0, Lsys)
print(f"    kappa=0:  max_f |I(R:F)| = {max(abs(x) for x in v0R):.3e}   "
      f"max_f |I(Lambda:F)| = {max(abs(x) for x in v0L):.3e}")

vR_hi, vL_hi, rR_hi, rL_hi = redundancy_block(
    0.8, LENP, f"SAME MEASUREMENT at the max-I(R:all) point kappa={KAPR}, T={TTR} (length on RIM link {LENP})",
    kap=KAPR, T=TTR)

print("\n  DOES beta CHANGE THE GAUGE NUMBERS TOO?  Sweep beta at the operating point.")
print(f"    {'beta':>6s} {'I(R:all)':>11s} {'I(R:|F|=1)':>12s} {'I(Lambda:all)':>15s} {'I(Lambda:|F|=1)':>17s}")
print("    " + "-" * 66)
BSWEEP = []
for beta in [0.0, 0.2, 0.4, 0.8, 1.2]:
    pT, DSY = run_patch(beta, EPS, G2, LENP, KAP, TT)
    a = profile(pT, DSY, Rsys); b = profile(pT, DSY, Lsys)
    BSWEEP.append((beta, a[-1], a[1], b[-1], b[1]))
    print(f"    {beta:6.2f} {a[-1]:11.6f} {a[1]:12.6f} {b[-1]:15.6f} {b[1]:17.6f}")

# ================================ PART 4: THE SIEVE ============================================
print()
print("-" * 100)
print("PART 4.  THE PREDICTABILITY SIEVE WITH THE METRIC ON.  Nominating nothing.")
print("-" * 100)
print("  Full Lindblad superoperator on the dim-32 system  ->  1024 x 1024, exact eigendecomposition.")
print("  Row-major generator  -i(H (x) I - I (x) H^T) + gamma*sum_k (L_k (x) L_k* - I (x) I),  L_k = Z_k on CUT.")
print("  Identification dictionary: the 16 magnetic subset-operators (x) {I, sx, sy, sz} on the length")
print("  qubit = 64 named operators.  The identity mode is skipped (it is the trace, rate 0 always).")

sy = np.array([[0, -1j], [1j, 0]], complex)
PAULI = [("I_len", I2), ("Lambda_x", sx), ("Lambda_y", sy), ("Lambda (sigma_z)", sz)]
MAGOPS = {}
for r in range(5):
    for S in itertools.combinations(range(4), r):
        O = Move(compose([PLQ[i] for i in S])) if S else IdS.copy()
        nm = "1" if not S else ("RIM LOOP" if r == 4 else f"{r}-plaquette {S}")
        MAGOPS[S] = (nm, O)
DICT = {}
for S, (nm, O) in MAGOPS.items():
    for pn, pm in PAULI:
        F = np.kron(O, pm)
        DICT[(S, pn)] = (f"{nm} (x) {pn}", F / np.linalg.norm(F))

DSY = DS * 2
IdY = np.eye(DSY, dtype=complex)
def sieve(beta, eps, g2=G2, gam=0.5, lenlink=None, topk=8, quiet=False):
    lenlink = LENP if lenlink is None else lenlink
    H = H_patch(beta, eps, g2, lenlink)
    M = -1j * (np.kron(H, IdY) - np.kron(IdY, H.T))
    for k in CUT:
        Lk = np.kron(ZL[k], I2)
        M = M + gam * (np.kron(Lk, Lk.conj()) - np.kron(IdY, IdY))
    w, U = np.linalg.eig(M.conj().T)
    rate = -np.conj(w).real
    order = np.argsort(rate)
    rows = []
    for i in order:
        O = U[:, i].reshape(DSY, DSY); n = np.linalg.norm(O)
        if n < 1e-12: continue
        O = O / n
        best = max(DICT.items(), key=lambda kv: abs(np.vdot(kv[1][1].reshape(-1), O.reshape(-1))))
        ov = abs(np.vdot(best[1][1].reshape(-1), O.reshape(-1)))
        nm = best[1][0] if ov > 0.30 else "mixed / no clean match"
        if nm.startswith("1 (x) I_len"): continue           # the identity mode: trace, rate 0 always
        rows.append((rate[i], nm, ov))
        if len(rows) >= topk: break
    if not quiet:
        print(f"\n  beta={beta}  eps={eps}  g2={g2}  gamma={gam}  length on link {lenlink}")
        print(f"    {'decay rate':>14s}  {'best operator match':<34s} {'overlap':>8s}")
        print("    " + "-" * 62)
        for r, nm, ov in rows:
            print(f"    {r:14.6e}  {nm:<34s} {ov:8.3f}")
    return rows

def named_rate(rows, key):
    for r, nm, ov in rows:
        if key in nm: return r
    return None

# ---- robust route: each named operator's OWN SURVIVAL, from the adjoint propagator ------------
# The eigenmode listing above labels eigenVECTORS, and where the Liouvillian is degenerate the
# eigenvectors are arbitrary mixtures, so the label is not trustworthy there.  The unambiguous
# quantity is the Heisenberg-picture survival of a FIXED operator:
#     vec(O(t)) = exp(M^dagger t) vec(O),    S(t) = |<O, O(t)>| / ||O||^2,
#     rate_eff(O) = -ln S(TS) / TS.
# This is exactly Zurek's sieve criterion -- HOW MUCH OF ITSELF the observable retains -- applied to
# observables instead of states.  Two notes, because both matter:
#   * the LONG-TIME log-slope is NOT usable here: past the fast timescales every operator with any
#     overlap on the globally slowest eigenmode reports that one eigenvalue, so it cannot separate
#     operators.  Checked: at beta=0.8 twelve different operators return 2.1352e-2 to six figures.
#   * the t->0 slope is not usable either: Re<O, L^dag O> has no Hamiltonian contribution, so it
#     returns exactly 0 for every observable that commutes with the four bath operators -- which is
#     the forced answer named in Part 0, not a measurement.
# TS is chosen between the two: long compared with the fast rates (~0.5-1), short compared with
# 1/rate of the slow sector.
TS = 10.0
def own_rates(beta, eps, g2=G2, gam=0.5, lenlink=None):
    lenlink = LENP if lenlink is None else lenlink
    H = H_patch(beta, eps, g2, lenlink)
    M = -1j * (np.kron(H, IdY) - np.kron(IdY, H.T))
    for k in CUT:
        Lk = np.kron(ZL[k], I2)
        M = M + gam * (np.kron(Lk, Lk.conj()) - np.kron(IdY, IdY))
    mu, U = np.linalg.eig(M.conj().T)
    Uinv = np.linalg.inv(U)
    ev = np.exp(mu * TS)
    out = {}
    for key, (nm, F) in DICT.items():
        if nm.startswith("1 (x) I_len"): continue        # the identity: exactly conserved, always
        v = F.reshape(-1); c = Uinv @ v
        s = abs(np.vdot(v, U @ (ev * c)))
        out[nm] = float(-np.log(max(s, 1e-15)) / TS)
    return out

print("\n  BASELINE (i): metric register present but DECOUPLED, beta=0, eps=0.")
print("  Forced, as declared in Part 0: Lambda is exactly conserved and must show rate 0.")
r00 = sieve(0.0, 0.0)
print("\n  BASELINE (ii): beta=0, eps>0.  Lambda precesses but nothing damps it: still forced.")
r0e = sieve(0.0, EPS)
print("\n  THE MEASUREMENT: beta>0 AND eps>0.  The length fluctuates AND sources the gauge couplings,")
print("  so the bath on the cut can reach it through the gauge field.  Nothing is nominated.")
rON = sieve(0.8, EPS)

print("\n  EXACTNESS (requirement 6): with beta=0 the 32-dim Liouvillian must factorise, so its")
print("  gauge-sector rates must equal the metric-free 16-dim carrier's rates exactly.")
Hb = -MAGH - G2 * sum(ZL[k] + ZL[k].conj().T for k in range(NL))
Mb = -1j * (np.kron(Hb, IdS) - np.kron(IdS, Hb.T))
for k in CUT: Mb = Mb + 0.5 * (np.kron(ZL[k], ZL[k].conj()) - np.kron(IdS, IdS))
wb = np.linalg.eig(Mb.conj().T)[0]; rb = np.sort(-np.conj(wb).real)
Hz = H_patch(0.0, 0.0, G2, LENP)
Mz = -1j * (np.kron(Hz, IdY) - np.kron(IdY, Hz.T))
for k in CUT:
    Lk = np.kron(ZL[k], I2); Mz = Mz + 0.5 * (np.kron(Lk, Lk.conj()) - np.kron(IdY, IdY))
wz = np.linalg.eig(Mz.conj().T)[0]; rz = np.sort(-np.conj(wz).real)
print(f"    slowest 6 rates, metric-free dim-16 carrier : {np.array2string(rb[:6], precision=6)}")
print(f"    slowest 6 rates, beta=0 dim-32 carrier      : {np.array2string(rz[:6], precision=6)}")
print(f"    max |difference| over the 256 metric-free rates matched into the 1024 : "
      f"{np.max(np.abs(rb - np.array([rz[np.argmin(np.abs(rz-x))] for x in rb]))):.3e}")

print("\n  RANKING ALL 63 NON-IDENTITY DICTIONARY OPERATORS BY THEIR OWN SURVIVAL rate_eff = -lnS(%.0f)/%.0f," % (TS, TS))
print("  at the measurement point beta=0.8, eps=%.2f.  Slowest first.  Nothing nominated: every" % EPS)
print("  operator in the dictionary is ranked and the ordering is read off afterwards.")
ratesON = own_rates(0.8, EPS)
ordON = sorted(ratesON.items(), key=lambda kv: kv[1])
print(f"    {'rank':>4s} {'decay rate':>16s}  operator")
print("    " + "-" * 56)
for i, (nm, r) in enumerate(ordON[:12]):
    print(f"    {i+1:4d} {r:16.6e}  {nm}")
gap = ordON[1][1] / max(ordON[0][1], 1e-18)
print(f"    slowest / second-slowest  =  {gap:.4g}x")

print("\n  THE COMPETITION, swept.  rate(RIM LOOP (x) I_len)  versus  the slowest operator that")
print("  carries the LENGTH register (any Lambda_x / Lambda_y / Lambda_z factor).")
print(f"    {'beta':>6s} {'eps':>6s} {'rate(RIM LOOP)':>16s} {'slowest length-carrying op':>27s} {'slower':>10s} {'ratio':>11s}  which length op")
print("    " + "-" * 116)
SWEEP = []
for beta in [0.0, 0.4, 0.8, 1.2]:
    for eps in [0.0, 0.35, 1.0]:
        rr = own_rates(beta, eps)
        rR = rr["RIM LOOP (x) I_len"]
        lam = sorted(((v, k) for k, v in rr.items() if "Lambda" in k))
        rL, nL_ = lam[0]
        win = "RIM LOOP" if rR < rL else ("LENGTH" if rL < rR else "tie")
        marg = max(rR, rL) / max(min(rR, rL), 1e-18)
        SWEEP.append((beta, eps, rR, rL, win, marg, nL_))
        print(f"    {beta:6.2f} {eps:6.2f} {rR:16.6e} {rL:27.6e} {win:>10s} {marg:11.4g}x  {nL_}")

print("\n  DOES THE METRIC CHANGE WHAT THE GAUGE RECORD COSTS?  rate(RIM LOOP (x) I_len) versus beta,")
print("  with the length on a RIM link (inside R's support) and on a CUT link (outside it).")
print(f"    {'beta':>6s} {'rate(R), length on RIM':>24s} {'rate(R), length on CUT':>24s}")
print("    " + "-" * 58)
RSWEEP = []
for beta in [0.0, 0.4, 0.8, 1.2]:
    a = own_rates(beta, EPS, lenlink=LENP)["RIM LOOP (x) I_len"]
    b = own_rates(beta, EPS, lenlink=LENC)["RIM LOOP (x) I_len"]
    RSWEEP.append((beta, a, b))
    print(f"    {beta:6.2f} {a:24.6e} {b:24.6e}")

print()
print("=" * 100)
print("NUMBERS ONLY.  No verdicts here.")
print("=" * 100)
print(f"  transport, full ring, metric ON      max I(R:probe)={mx(full_on,0):.6f}  "
      f"max I(Lam_live:probe)={mx(full_on,1):.6f}  max I(Lam_dead:probe)={mx(full_on,2):.6f}")
print(f"  transport, full ring, metric OFF     max I(R:probe)={mx(full_off,0):.6f}  "
      f"max I(Lam_live:probe)={mx(full_off,1):.3e}  max I(Lam_dead:probe)={mx(full_off,2):.3e}")
print(f"  transport, CUT ring,  metric ON      max I(R:probe)={mx(cut_on,0):.3e}  "
      f"max I(Lam_live:probe)={mx(cut_on,1):.6f}  max I(Lam_dead:probe)={mx(cut_on,2):.3e}")
print(f"  redundancy (len on RIM)  gauge: I(1)={vR_on[1]:.4f} I(all)={vR_on[-1]:.4f} ratio={rR_on:.4f}")
print(f"  redundancy (len on RIM)  metric: I(1)={vL_on[1]:.4f} I(all)={vL_on[-1]:.4f} ratio={rL_on:.4f}")
print(f"  redundancy (len on CUT)  gauge: I(1)={vR_c[1]:.4f} I(all)={vR_c[-1]:.4f} ratio={rR_c:.4f}")
print(f"  redundancy (len on CUT)  metric: I(1)={vL_c[1]:.4f} I(all)={vL_c[-1]:.4f} ratio={rL_c:.4f}")
print(f"  redundancy control beta=0        metric: max_f |I| = {max(abs(x) for x in vL_off):.3e}")
for beta, a1, a2, b1, b2 in BSWEEP:
    print(f"  redundancy beta={beta:.2f}  I(R:all)={a1:.6f} I(R:1)={a2:.6f} "
          f"I(Lam:all)={b1:.6f} I(Lam:1)={b2:.6f}")
for beta, eps, rR, rL, win, marg, nmL in SWEEP:
    print(f"  sieve beta={beta:.2f} eps={eps:.2f}   rate(RIM LOOP)={rR:.6e}   "
          f"rate(slowest length op)={rL:.6e}   slower={win}  ratio={marg:.4g}   [{nmL}]")
for beta, a, b in RSWEEP:
    print(f"  rate(RIM LOOP) beta={beta:.2f}   length on RIM={a:.6e}   length on CUT={b:.6e}")
print(f"  sieve slowest overall at beta=0.8, eps={EPS}:  {ordON[0][0]}  rate={ordON[0][1]:.6e}; "
      f"second {ordON[1][0]} rate={ordON[1][1]:.6e}; ratio {gap:.4g}x")
print("=" * 100)
