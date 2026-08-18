"""W-44.  A METRIC THAT DOES WORK.  A CLOCK CARRIED BY THE TRANSPORTED PROBE.

THE GAP.  Everything W-27..W-42 measured is pure gauge structure: the carrier's links have no
length, no proper time, no clock.  W-43b tried to close this and failed in a specific way: it put a
two-level "mass" at a site, set that subsystem's Hamiltonian to ZERO, and coupled the environment
STRAIGHT to the mass.  The potential such a mass would source never entered any evolution, and the
mass was conserved by construction, so its redundancy (1.000) was forced by a commutator, not
measured.  This lane must not repeat that.

WHAT IS INSTALLED HERE.  A single two-level METRIC degree of freedom M (occupation N_M = 0/1: "no
mass"/"mass at the corner", or a superposition of the two).  Its value sources a LAPSE field

        N(x)  =  1 - lam * phi(x) * N_M ,      phi(x) = 1/(1+|x - x_mass|)  (normalised to max 1)

and the lapse multiplies the local Hamiltonian density everywhere:

  (i)   GAUGE CARRIER (the 3x3 patch).  Each plaquette term and each electric link term is weighted
        by the lapse at that plaquette / that link.  Curvature therefore changes how fast each
        region of the gauge field evolves.  The metric is not a spectator: it is a factor in H.
  (ii)  TRANSPORT (the 8-site ring with an explicit probe charge and a static anti-charge).  The
        hop amplitude on each link is weighted by the lapse on that link -- a LINK LENGTH /
        TRAVERSAL TIME.  The probe also feels the potential itself.
  (iii) THE CLOCK.  The probe carries an internal dc-level clock whose RATE is the local lapse:
              H_clock = omega * sum_v |v><v|_probe (x) N(v) (x) C,   C = traceless diag ladder.
        A clock deep in the potential ticks slow.  Two arms of the ring that pass through different
        potential accumulate different proper time, so the clock state becomes correlated with the
        path.  That is the configuration in which proper time is physically load-bearing.
  (iv)  The metric is NOT conserved: H_M = -eps * sx_M.  Whether N_M survives is a question for the
        dynamics, not a construction.  It also back-reacts: N_M multiplies operators that do not
        commute with sx_M, so the gauge field pushes on the geometry as well.

lam = 0 removes every one of (i)-(iv)'s couplings at once and must return exactly the no-metric
answers of W-34 / W-36 / W-37.

SIMULATION ROUTES (stated, per requirement 2):
  * REDUNDANCY  -- exact UNITARY evolution of a joint pure state, system(32) (x) 6 env qubits
    (dim 2048).  No superoperator.  The env coupling is diagonal in the env computational basis, so
    the propagator is block diagonal: 64 exact expm's of 32x32 blocks.  expm is scaling-and-squaring
    written out below; no scipy.
  * SIEVE       -- FULL Lindblad superoperator, exactly exponentiated.  System dim 32, so the
    generator is 1024 x 1024.  Left modes from eig(M^dagger); operator decay rates from
    ||exp(M^dagger T) vec(O)||.  Row-major convention:  M = -i(H (x) I - I (x) H^T) + gamma sum
    (L (x) L* - I (x) I).
  * CLOCK       -- exact unitary evolution, dim 128 * 2 * 4 = 1024.  Hermitian eigendecomposition
    for the time scans; verified against the hand-written expm at one time.
  A clock-carrying SIEVE (patch (x) metric (x) 2-level clock, dim 64) would need a 4096x4096
    generator; it is not run, and that is stated rather than hidden.

NO SMUGGLING.  The bath is W-34's bath verbatim: Z on the four cut links, extended by (x) I on the
metric.  Nothing in the bath, and nothing in the environment coupling, ever mentions the rim loop,
the metric occupation, or the clock.  Every metric-to-environment channel has to be routed through
[N_M, H].
"""
import itertools, numpy as np

# ============================================================================================
#  generic numerics -- expm by scaling-and-squaring, von Neumann entropy
# ============================================================================================
def expm(A):
    nr = np.linalg.norm(A, np.inf)
    k = max(0, int(np.ceil(np.log2(nr))) + 1) if nr > 0 else 0
    B = A / (2.0 ** k)
    X = np.eye(A.shape[0], dtype=complex); T = X.copy()
    for m in range(1, 60):
        T = T @ B / m; X = X + T
        if np.linalg.norm(T, np.inf) < 1e-18 * max(1.0, np.linalg.norm(X, np.inf)): break
    for _ in range(k): X = X @ X
    return X

def vn(rho):
    ev = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    ev = ev[ev > 1e-12]
    return float(-(ev * np.log2(ev)).sum())

def opcard(name, O):
    ev = np.linalg.eigvals(O)
    nd = len(np.unique(np.round(ev, 6)))
    ud = np.linalg.norm(O @ O.conj().T - np.eye(O.shape[0]))
    print(f"    {name:<34s} dim {O.shape[0]:5d}  ||O||_F {np.linalg.norm(O):9.4f}"
          f"  unitarity defect {ud:9.2e}  distinct eigenvalues {nd:3d}")
    return nd

# ============================================================================================
#  the 3x3 patch (identical to W-34 / W-36) + the metric profile it lives in
# ============================================================================================
V2 = [(i, j) for j in range(3) for i in range(3)]; vid = {v: k for k, v in enumerate(V2)}
E = []
for j in range(3):
    for i in range(2): E.append((vid[(i, j)], vid[(i + 1, j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i, j)], vid[(i, j + 1)]))
NL = len(E)
hid = lambda i, j: j * 2 + i
vx  = lambda i, j: 6 + j * 3 + i
P = [[(hid(i, j), +1), (vx(i + 1, j), +1), (hid(i, j + 1), -1), (vx(i, j), -1)]
     for j in range(2) for i in range(2)]
CENTER = vid[(1, 1)]
CUT   = [k for k, (a, b) in enumerate(E) if a == CENTER or b == CENTER]
PERIM = [k for k in range(NL) if k not in CUT]

def build(N=2):
    st = [s for s in itertools.product(range(N), repeat=NL)
          if all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
                - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % N == 0 for v in range(9))]
    return st, {s: i for i, s in enumerate(st)}
ST, IX = build(); DG = len(ST)

def Zop(links):
    return np.diag([(-1.0) ** (sum(s[k] for k in links) % 2) for s in ST]).astype(complex)
def Move(mv):
    M = np.zeros((DG, DG), complex)
    for j, s in enumerate(ST):
        t = list(s)
        for k, sg in mv: t[k] = (t[k] + sg) % 2
        t = tuple(t)
        if t in IX: M[IX[t], j] = 1.0
    return M
def compose(ps):
    acc = {}
    for p in ps:
        for k, sg in p: acc[k] = acc.get(k, 0) + sg
    return [(k, s) for k, s in acc.items() if s % 2]

B  = [Move(p) for p in P]                      # 4 plaquette flips
RG = Move(compose(P))                          # rim Wilson loop  (discrete Stokes)
Zl = [Zop([k]) for k in range(NL)]
I16 = np.eye(DG, dtype=complex)

# ---- the geometry: a mass at the corner vertex (0,0); phi = 1/(1+r), normalised to max 1 --------
MASS_XY = (0.0, 0.0)
def _phi(x, y):
    return 1.0 / (1.0 + np.hypot(x - MASS_XY[0], y - MASS_XY[1]))
plq_xy  = [((k % 2) + 0.5, (k // 2) + 0.5) for k in range(4)]
lnk_xy  = [None] * NL
for j in range(3):
    for i in range(2): lnk_xy[hid(i, j)] = (i + 0.5, float(j))
for j in range(2):
    for i in range(3): lnk_xy[vx(i, j)] = (float(i), j + 0.5)
_raw = [_phi(*p) for p in plq_xy] + [_phi(*l) for l in lnk_xy]
_nrm = max(_raw)
PHI_P = np.array([_phi(*p) / _nrm for p in plq_xy])
PHI_L = np.array([_phi(*l) / _nrm for l in lnk_xy])

sz = np.array([[1, 0], [0, -1]], complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
I2 = np.eye(2, dtype=complex)
NM = (I2 - sz) / 2                              # occupation: 0 = no mass, 1 = mass

G2, EPS, GAM = 0.01, 0.15, 0.5                  # electric coupling, metric tunnelling, bath rate
LAM = 0.6                                        # THE metric coupling.  lam = 0 -> no metric.

def H_patch_metric(lam, g2=G2, eps=EPS):
    """gauge patch (x) metric.  Lapse N(x) = 1 - lam*phi(x)*N_M multiplies every local term."""
    H = np.zeros((DG * 2, DG * 2), complex)
    for k in range(4):
        H -= np.kron(B[k] + B[k].conj().T, I2 - lam * PHI_P[k] * NM)
    for l in range(NL):
        H -= g2 * np.kron(Zl[l] + Zl[l].conj().T, I2 - lam * PHI_L[l] * NM)
    H -= eps * np.kron(I16, sx)
    return H
def H_patch_nometric(g2=G2):
    MAG = sum(b + b.conj().T for b in B)
    ELE = sum(z + z.conj().T for z in Zl)
    return -MAG - g2 * ELE

R_GAUGE  = np.kron(RG, I2)                       # gauge record on the joint carrier
R_METRIC = np.kron(I16, sz)                      # metric record on the joint carrier

print("=" * 108)
print("W-44  A METRIC THAT DOES WORK.  Lapse N(x)=1-lam*phi(x)*N_M multiplies the local Hamiltonian")
print("      density; the transported probe carries a clock that ticks at the local lapse.")
print("=" * 108)
print(f"\n  patch: {DG} physical gauge states, {NL} links, cut links {CUT}, rim links {PERIM}")
print(f"  mass at vertex (0,0).  lapse profile (lam={LAM}):")
print("      plaquette phi : " + "  ".join(f"{p:.3f}" for p in PHI_P))
print("      link      phi : " + "  ".join(f"{p:.3f}" for p in PHI_L))
print(f"      so N(x) ranges over [{1-LAM*max(PHI_P.max(),PHI_L.max()):.3f}, "
      f"{1-LAM*min(PHI_P.min(),PHI_L.min()):.3f}] -- a genuine spatial gradient, not a global rescaling")

# ============================================================================================
#  1.  FORCED-OR-NOT, BEFORE ANY DYNAMICS RUNS
# ============================================================================================
print("\n" + "=" * 108)
print("1.  FORCED-OR-NOT.  Named in advance, settled before a single propagator is built.")
print("=" * 108)
print("""
  AT RISK: I(N_M : F), the Holevo information an environment fragment holds about the METRIC.
  This is exactly the quantity W-43b got for free.  There it was FORCED twice over: the mass had
  H = 0 so [N_M, H] = 0 exactly (a conserved charge cannot be damaged, so its "durability" measured
  nothing), and the environment was coupled DIRECTLY to sz_M (so each fragment held a full copy by
  construction).  Either alone determines the number without any dynamics.

  THE ONE ARGUMENT THAT SETTLES IT HERE -- a commutator, printed below before anything propagates:
      (a) ||[N_M, H]||  must be NONZERO, or the metric record is conserved and its survival is a
          constraint, not a measurement.  H_M = -eps*sx_M gives [N_M, sx_M] = -sx_M*sz_M ... != 0.
      (b) ||[N_M, L_k]|| must be ZERO for every bath / environment operator, or the environment is
          reading the answer off directly.  L_k = Z_cut(k) (x) I is diagonal in the gauge basis and
          acts as the identity on the metric, so this vanishes identically.
  (a) and (b) together mean: the environment can only learn about the metric through the gauge field
  the metric curves.  Any nonzero I(N_M:F) is then earned by the dynamics.

  SECOND, SMALLER RISK, also settled by a commutator: the SIEVE WINNER.  [RG, B_p] = 0 for every
  plaquette (abelian Z_2 moves commute) and [RG, N_M] = 0, so weighting the plaquettes by the lapse
  does NOT weaken the rim loop's protection: [R_GAUGE, H_magnetic(lam)] = 0 for every lam.  The rim
  loop therefore cannot be dislodged by the metric spoiling its commutator.  It can only be
  dislodged if some metric-carrying operator turns out to be even slower.  That is the open question.
""")
_H0 = H_patch_metric(0.0); _HL = H_patch_metric(LAM)
_Ls = [np.kron(Zl[k], I2) for k in CUT]
print(f"    ||[N_M(x)I , H(lam=0)]||           = {np.linalg.norm(np.kron(I16,NM)@_H0-_H0@np.kron(I16,NM)):.6f}")
print(f"    ||[N_M(x)I , H(lam={LAM})]||         = {np.linalg.norm(np.kron(I16,NM)@_HL-_HL@np.kron(I16,NM)):.6f}")
print(f"    max_k ||[N_M(x)I , L_k]||          = "
      f"{max(np.linalg.norm(np.kron(I16,NM)@L-L@np.kron(I16,NM)) for L in _Ls):.3e}")
print(f"    ||[R_GAUGE , H_magnetic(lam=0)]||  = "
      f"{np.linalg.norm(R_GAUGE@H_patch_metric(0.0,g2=0,eps=0)-H_patch_metric(0.0,g2=0,eps=0)@R_GAUGE):.3e}")
print(f"    ||[R_GAUGE , H_magnetic(lam={LAM})]|| = "
      f"{np.linalg.norm(R_GAUGE@H_patch_metric(LAM,g2=0,eps=0)-H_patch_metric(LAM,g2=0,eps=0)@R_GAUGE):.3e}")
print(f"    ||[R_GAUGE , H_full(lam={LAM})]||     = {np.linalg.norm(R_GAUGE@_HL-_HL@R_GAUGE):.6f}   "
      f"(the electric term is what damages it, exactly as at lam=0)")

# ============================================================================================
#  2.  OPERATORS EARN MEASUREMENT
# ============================================================================================
print("\n" + "=" * 108)
print("2.  OPERATORS EARN MEASUREMENT.  Every operator that will be given a commutator or a")
print("    Holevo read, with its norm, unitarity defect and eigenvalue count.")
print("=" * 108)
opcard("R_GAUGE  = rim loop (x) I", R_GAUGE)
opcard("R_METRIC = I (x) sz_M", R_METRIC)
opcard("N_M      = I (x) (1-sz)/2", np.kron(I16, NM))

# ============================================================================================
#  3.  MEASUREMENT (a):  REDUNDANCY.  Same machinery for gauge and metric.
# ============================================================================================
print("\n" + "=" * 108)
print("3.  MEASUREMENT (a).  REDUNDANCY.  One joint pure state; ONE environment; two records read")
print("    off the SAME state with the SAME fragments, so the two numbers are comparable.")
print("    Environment: 6 qubits, each coupled to Z on a CUT link.  It never mentions either record.")
print("=" * 108)

NQ = 6
def env_state(Hsys, Dsys, cutlinks, kappa, T, psiS):
    """exact joint state.  H_int = kappa sum_k Zcut(k) (x) sz_k is diagonal in the env basis,
       so U is block diagonal: one 32x32 expm per env basis state."""
    out = np.zeros((Dsys, 2 ** NQ), complex)
    for e in range(2 ** NQ):
        Hb = Hsys.copy()
        for k in range(NQ):
            s = 1.0 - 2.0 * ((e >> (NQ - 1 - k)) & 1)
            Hb = Hb + kappa * s * cutlinks[k % len(cutlinks)]
        out[:, e] = (expm(-1j * Hb * T) @ psiS) / np.sqrt(2.0 ** NQ)
    return out.reshape((Dsys,) + (2,) * NQ)

def holevo(psiT, Dsys, Proj, frag):
    br = []
    for Pr in Proj:
        v = np.tensordot(Pr, psiT, axes=([1], [0]))
        p = float(np.vdot(v, v).real)
        if p < 1e-14: br.append((0.0, None)); continue
        v = v / np.sqrt(p)
        keep = [0] + [1 + i for i in frag]
        tr = [ax for ax in range(1 + NQ) if ax not in keep]
        M = np.transpose(v, keep + tr).reshape(Dsys * 2 ** len(frag), -1)
        rho = (M @ M.conj().T).reshape(Dsys, 2 ** len(frag), Dsys, 2 ** len(frag))
        br.append((p, np.einsum('ijik->jk', rho)))
    avg = sum(p * m for p, m in br if m is not None)
    return vn(avg) - sum(p * vn(m) for p, m in br if m is not None)

def profile(psiT, Dsys, Rec):
    Pp = (np.eye(Dsys) + Rec) / 2; Pm = (np.eye(Dsys) - Rec) / 2
    out = []
    for f in range(NQ + 1):
        combos = list(itertools.combinations(range(NQ), f))[:20]
        out.append(float(np.mean([holevo(psiT, Dsys, (Pp, Pm), c) for c in combos])))
    return out

def init_joint(lam, seed=3):
    g = np.random.default_rng(seed)
    Pp = (np.eye(DG) + RG) / 2; Pm = (np.eye(DG) - RG) / 2
    w = g.normal(size=DG) + 1j * g.normal(size=DG)
    a = Pp @ w; b = Pm @ w; a /= np.linalg.norm(a); b /= np.linalg.norm(b)
    gs = (a + b) / np.sqrt(2.0); gs /= np.linalg.norm(gs)          # <RG> = 0 exactly
    ms = np.ones(2, complex) / np.sqrt(2.0)                        # <sz_M> = 0 exactly
    return np.kron(gs, ms), gs

CUT32 = [np.kron(Zl[k], I2) for k in CUT]
CUT16 = [Zl[k] for k in CUT]

print("\n  SATURATION FIRST.  A fragment plot means nothing until the WHOLE environment holds the")
print("  bit.  Scan the coupling; read the shape only where I(gauge : all 6) has saturated.")
print(f"    {'kappa':>7s} {'T':>6s} {'I(gauge:all)':>14s} {'I(metric:all)':>15s}")
print("    " + "-" * 46)
psiS32, _ = init_joint(LAM)
best = None
for kap, T in [(0.6, 6.0), (1.2, 6.0), (2.0, 8.0), (3.0, 12.0), (5.0, 12.0), (8.0, 16.0),
               (12.0, 20.0), (20.0, 30.0)]:
    ps = env_state(H_patch_metric(LAM), DG * 2, CUT32, kap, T, psiS32)
    ig = holevo(ps, DG * 2, ((np.eye(DG*2)+R_GAUGE)/2, (np.eye(DG*2)-R_GAUGE)/2), tuple(range(NQ)))
    im = holevo(ps, DG * 2, ((np.eye(DG*2)+R_METRIC)/2, (np.eye(DG*2)-R_METRIC)/2), tuple(range(NQ)))
    print(f"    {kap:7.2f} {T:6.1f} {ig:14.6f} {im:15.6f}")
    if best is None or ig > best[2]: best = (kap, T, ig)
KAP, TT = best[0], best[1]
print(f"    -> reading the shape at kappa={KAP}, T={TT}")

rows = {}
for tag, lam in [("metric ON  (lam=%.1f)" % LAM, LAM), ("metric OFF (lam=0)", 0.0)]:
    psiS, _ = init_joint(lam)
    ps = env_state(H_patch_metric(lam), DG * 2, CUT32, KAP, TT, psiS)
    pg = profile(ps, DG * 2, R_GAUGE)
    pm = profile(ps, DG * 2, R_METRIC)
    rows[tag] = (pg, pm)
    print(f"\n  {tag}")
    print("      |F|            : " + "  ".join(f"{i:8d}" for i in range(NQ + 1)))
    print("      I(gauge  : F)  : " + "  ".join(f"{v:8.5f}" for v in pg))
    print("      I(metric : F)  : " + "  ".join(f"{v:8.5f}" for v in pm))

# the no-metric reference: the 16-dim carrier with no metric factor at all
psiS16 = init_joint(0.0)[1]
ps16 = env_state(H_patch_nometric(), DG, CUT16, KAP, TT, psiS16)
pg16 = profile(ps16, DG, RG)
print("\n  NO-METRIC REFERENCE (dim 16 carrier, no metric factor present at all)")
print("      I(gauge  : F)  : " + "  ".join(f"{v:8.5f}" for v in pg16))
print(f"      max |lam=0 (dim 32)  -  no-metric (dim 16)| over all |F|  = "
      f"{max(abs(a-b) for a,b in zip(rows['metric OFF (lam=0)'][0], pg16)):.3e}")

print(f"\n  {'case':>24s} {'I(|F|=1)':>10s} {'I(|F|=3)':>10s} {'I(all)':>10s} {'I(1)/I(all)':>12s}")
print("  " + "-" * 70)
RED = {}
for tag in rows:
    for nm, v in (("gauge", rows[tag][0]), ("metric", rows[tag][1])):
        tot = v[-1]
        r = v[1] / tot if abs(tot) > 1e-9 else float('nan')
        RED[(tag, nm)] = (v[1], v[3], tot, r)
        rs = f"{r:12.4f}" if tot > 1e-9 else f"{'n/a':>12s}"
        print(f"  {nm+' | '+tag:>24s} {v[1]:10.5f} {v[3]:10.5f} {tot:10.5f} {rs}")
print("\n  CONTROL that must return zero if the mechanism is absent: at lam=0 the metric is coupled")
print("  to nothing the environment touches, so I(metric : F) must vanish identically.")
print(f"    max_F |I(metric:F)| at lam=0 = {max(abs(x) for x in rows['metric OFF (lam=0)'][1]):.3e}")
print(f"    max_F |I(metric:F)| at lam={LAM} = {max(abs(x) for x in rows['metric ON  (lam=%.1f)'%LAM][1]):.3e}")

# ============================================================================================
#  4.  MEASUREMENT (b):  THE PREDICTABILITY SIEVE WITH THE METRIC ON.  NOMINATE NOTHING.
# ============================================================================================
print("\n" + "=" * 108)
print("4.  MEASUREMENT (b).  THE PREDICTABILITY SIEVE, metric term switched ON.  Nothing nominated:")
print("    the slowest left-eigenmodes of the full Lindbladian are read off and only then matched")
print("    against a dictionary of 64 named operators (16 magnetic subsets x {I, sz, sx, sy} on the")
print("    metric).  Bath = W-34's bath verbatim: Z on the 4 cut links, gamma=0.5, (x) I on metric.")
print("=" * 108)

DICT = {}
for r in range(5):
    for S in itertools.combinations(range(4), r):
        Og = Move(compose([P[i] for i in S])) if S else I16.copy()
        gname = ("I_gauge" if not S else "RIM LOOP(all 4)" if r == 4 else
                 f"plaq{S}" if r == 1 else f"{r}-plaq loop{S}")
        for mn, Om in (("(x)I_M", I2), ("(x)sz_M", sz), ("(x)sx_M", sx), ("(x)sy_M", sy)):
            O = np.kron(Og, Om)
            DICT[(S, mn)] = (gname + mn, O / np.linalg.norm(O))

def liouvillian(H, D, Ls, gam=GAM):
    Id = np.eye(D, dtype=complex)
    M = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for L in Ls: M += gam * (np.kron(L, L.conj()) - np.kron(Id, Id))
    return M

def sieve(lam, Trate=30.0, topk=8, label=""):
    D = DG * 2
    H = H_patch_metric(lam)
    M = liouvillian(H, D, [np.kron(Zl[k], I2) for k in CUT])
    Mdag = M.conj().T
    w, U = np.linalg.eig(Mdag)
    rate = -np.real(w)
    order = np.argsort(rate)
    print(f"\n  --- lam = {lam}   {label}")
    print(f"    slowest left-eigenmodes of the Lindbladian (dim {D}, generator {D*D}x{D*D})")
    print(f"      {'decay rate':>14s} {'|diag wt|':>10s}  best dictionary match")
    print("      " + "-" * 74)
    shown = 0
    for i in order:
        if shown >= topk: break
        O = U[:, i].reshape(D, D); n = np.linalg.norm(O)
        if n < 1e-12: continue
        O = O / n
        dw = np.linalg.norm(np.diag(np.diag(O)))
        best = max(DICT.values(), key=lambda kv: abs(np.vdot(kv[1].reshape(-1), O.reshape(-1))))
        ov = abs(np.vdot(best[1].reshape(-1), O.reshape(-1)))
        tag = best[0] if ov > 0.3 else ("diagonal / electric-type" if dw > 0.7 else "no clean match")
        print(f"      {rate[i]:14.6e} {dw:10.3f}  {tag}   (overlap {ov:.3f})")
        shown += 1
    # exact operator decay rates: -(1/T) ln ||exp(Mdag T) vec(O)|| / ||vec(O)||
    ET = expm(Mdag * Trate)
    rr = []
    for key, (nm, O) in DICT.items():
        if key[0] == () and key[1] == "(x)I_M": continue          # identity: rate 0 by trace
        v = O.reshape(-1); vt = ET @ v
        g = -np.log(max(np.linalg.norm(vt) / np.linalg.norm(v), 1e-300)) / Trate
        rr.append((g, nm))
    rr.sort()
    print(f"    operator decay rates  -(1/T)ln||exp(L^dag T)O||/||O||,  T={Trate}   (10 slowest of 63)")
    for g, nm in rr[:10]:
        print(f"      {g:14.6e}   {nm}")
    return rr

print("\n  no-metric reference first: the SAME machinery on the bare 16-dim carrier (W-34's setting).")
Mref = liouvillian(H_patch_nometric(), DG, CUT16).conj().T
ETref = expm(Mref * 30.0)
refrates = {}
for r in range(5):
    for S in itertools.combinations(range(4), r):
        if not S: continue
        Og = Move(compose([P[i] for i in S])); Og = Og / np.linalg.norm(Og)
        v = Og.reshape(-1); vt = ETref @ v
        refrates[S] = -np.log(np.linalg.norm(vt) / np.linalg.norm(v)) / 30.0
for S in sorted(refrates, key=lambda k: refrates[k]):
    nm = "RIM LOOP(all 4)" if len(S) == 4 else (f"plaq{S}" if len(S) == 1 else f"{len(S)}-plaq loop{S}")
    print(f"      {refrates[S]:14.6e}   {nm}")

rr0 = sieve(0.0,  label="(metric present but DECOUPLED -- the control)")
rrL = sieve(LAM,  label="(metric ON: the lapse weights every plaquette and every link)")

d = {nm: g for g, nm in rr0}
mx = max(abs(d[nm + "(x)I_M"] - refrates[S]) for r in range(1, 5)
         for S in itertools.combinations(range(4), r)
         for nm in [("RIM LOOP(all 4)" if len(S) == 4 else
                     f"plaq{S}" if len(S) == 1 else f"{len(S)}-plaq loop{S}")])
print(f"\n  CONTROL: lam=0 gauge-operator rates vs the no-metric carrier -- max |difference| = {mx:.3e}")

win_g, win_n = rrL[0]
run_g, run_n = next((g, n) for g, n in rrL[1:] if n != win_n)
print(f"\n  metric ON: slowest = {win_n}  at rate {win_g:.6e}")
print(f"             next    = {run_n}  at rate {run_g:.6e}")
print(f"             ratio next/slowest = {(run_g/win_g if win_g>0 else float('inf')):.3f}")
d0 = {nm: g for g, nm in rr0}
dL = {nm: g for g, nm in rrL}
print(f"\n  DOES THE METRIC CHANGE THE SIEVE'S NUMBERS?  (lam=0 -> lam={LAM})")
for nm in ["RIM LOOP(all 4)(x)I_M", "RIM LOOP(all 4)(x)sz_M", "I_gauge(x)sz_M",
           "I_gauge(x)sx_M", "plaq(0,)(x)I_M"]:
    print(f"      {nm:<26s} {d0[nm]:14.6e}  ->  {dL[nm]:14.6e}")

# ============================================================================================
#  5.  THE CLOCK.  Transport with an internal clock whose rate is the local lapse.
# ============================================================================================
print("\n" + "=" * 108)
print("5.  THE CLOCK.  The probe of W-37 now carries a dc-level clock ticking at the LOCAL LAPSE.")
print("    Different arms of the ring pass through different potential, accumulate different proper")
print("    time, and the clock state becomes correlated with the path.")
print("=" * 108)

RIMV = [vid[(0, 0)], vid[(1, 0)], vid[(2, 0)], vid[(2, 1)], vid[(2, 2)], vid[(1, 2)], vid[(0, 2)], vid[(0, 1)]]
def link_between(u, v):
    for k, (a, b) in enumerate(E):
        if (a, b) == (u, v) or (b, a) == (u, v): return k
    raise ValueError
RIML = [link_between(RIMV[k], RIMV[(k + 1) % 8]) for k in range(8)]
assert sorted(RIML) == sorted(PERIM)

def div(s, v):
    return (sum(s[k] for k, (a, b) in enumerate(E) if a == v)
          - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % 2
STA = []
for p in range(8):
    for s in itertools.product(range(2), repeat=NL):
        q = lambda v: (1 if v == RIMV[p] else 0) ^ (1 if v == CENTER else 0)
        if all(div(s, v) == q(v) for v in range(9)): STA.append((p, s))
IDA = {x: i for i, x in enumerate(STA)}; DT = len(STA)
POS = np.array([p for p, _ in STA])
print(f"\n  transport carrier: {DT} physical states ({DT//8} gauge states per probe position),"
      f" static anti-charge at the centre")

def shift(s, links):
    t = list(s)
    for k in links: t[k] ^= 1
    return tuple(t)
RT = np.zeros((DT, DT), complex)
for j, (p, s) in enumerate(STA): RT[IDA[(p, shift(s, PERIM))], j] = 1.0
HOP = []
for k in range(8):
    A = np.zeros((DT, DT), complex)
    for j, (p, s) in enumerate(STA):
        if p != k: continue
        i = IDA.get(((k + 1) % 8, shift(s, [RIML[k]])))
        if i is not None: A[i, j] = 1.0
    HOP.append(A)

MASS_SITE = 2                                     # the mass sits beside ring site 2
def dring(a, b): return min((a - b) % 8, (b - a) % 8)
PHI_V = np.array([1.0 / (1.0 + dring(v, MASS_SITE)) for v in range(8)])
PHI_V = PHI_V / PHI_V.max()
PHI_LK = np.array([(PHI_V[k] + PHI_V[(k + 1) % 8]) / 2 for k in range(8)])
PHI_T = np.diag(PHI_V[POS]).astype(complex)
print("      ring potential phi(v) : " + "  ".join(f"{p:.3f}" for p in PHI_V))
print(f"      arm A (sites 1,2,3) mean phi = {PHI_V[[1,2,3]].mean():.3f} ;"
      f"  arm B (sites 7,6,5) mean phi = {PHI_V[[7,6,5]].mean():.3f}")

DC = 4
# The clock generator is taken TRACELESS.  Shifting C by a constant is a redefinition of the clock's
# energy zero -- it changes no clock reading -- but it removes the constant mean force the clock
# coupling would otherwise exert on the probe.  With <C> = 0 in the running clock state the
# "clock only" configuration below carries NO mean potential at all: whatever it does to the fringe
# it does by CORRELATING the clock with the path, and by nothing else.
Cn = (np.diag(np.arange(DC)) - (DC - 1) / 2.0).astype(complex)
Ic = np.eye(DC, dtype=complex)
IT = np.eye(DT, dtype=complex)
TAU, OMEGA, VG, EPSC = 1.0, 1.2, 1.0, 0.15

opcard("R_TRANSPORT = perimeter loop", RT)
opcard("R_TRANSPORT (x) I_M (x) I_C", np.kron(np.kron(RT, I2), Ic))
print(f"    {'clock generator C (traceless)':<34s} dim {DC:5d}  ||O||_F {np.linalg.norm(Cn):9.4f}"
      f"  Hermitian, trace {np.trace(Cn).real:+.1e}, not a record; distinct eigenvalues {DC:3d}")

def H_clockspace(lam, sw_hop=1.0, sw_pot=1.0, sw_clk=1.0, omega=OMEGA, cut=False, eps=EPSC):
    H = np.zeros((DT * 2 * DC, DT * 2 * DC), complex)
    for k in range(8):
        if cut and k == 4: continue
        Lap = I2 - lam * sw_hop * PHI_LK[k] * NM
        blk = np.kron(np.kron(HOP[k], Lap), Ic)
        H -= TAU * (blk + blk.conj().T)
    H -= VG * lam * sw_pot * np.kron(np.kron(PHI_T, NM), Ic)
    H += omega * (np.kron(np.kron(IT, I2), Cn)
                  - lam * sw_clk * np.kron(np.kron(PHI_T, NM), Cn))
    H -= eps * np.kron(np.kron(IT, sx), Ic)
    return H

Pp_T = (np.eye(DT) + RT) / 2; Pm_T = (np.eye(DT) - RT) / 2
def branch_vectors(seed=5):
    g = np.random.default_rng(seed)
    mask = np.array([1.0 if p == 0 else 0.0 for p, _ in STA])
    w = (g.normal(size=DT) + 1j * g.normal(size=DT)) * mask
    a = Pp_T @ w; b = Pm_T @ w
    return a / np.linalg.norm(a), b / np.linalg.norm(b)
AVEC, BVEC = branch_vectors()
CLK0 = np.ones(DC, complex) / np.sqrt(DC)                 # a RUNNING clock
def full(v, m):
    return np.kron(np.kron(v, m), CLK0)
MET = {"|0> no mass": np.array([1, 0], complex),
       "|1> mass":    np.array([0, 1], complex),
       "|+> superpos": np.ones(2, complex) / np.sqrt(2)}

PAIRS = []
_bys = {}
for i, (p, s) in enumerate(STA): _bys.setdefault(s, []).append(i)
for s, lst in _bys.items():
    for i in lst:
        for j in lst: PAIRS.append((i, j))
PAIRS = np.array(PAIRS)

def evolve_factory(H):
    ev, V = np.linalg.eigh(H)
    def U(psi, T): return V @ (np.exp(-1j * ev * T) * (V.conj().T @ psi))
    return U

def psite(psi, site):
    M = psi.reshape(DT, 2 * DC)
    return float(np.sum(np.abs(M[POS == site]) ** 2))
def probe_rdm(psi):
    M = psi.reshape(DT, 2 * DC)
    rho = np.zeros((8, 8), complex)
    G = M @ M.conj().T
    np.add.at(rho, (POS[PAIRS[:, 0]], POS[PAIRS[:, 1]]), G[PAIRS[:, 0], PAIRS[:, 1]])
    return rho
def clock_rdm(psi):
    M = psi.reshape(DT * 2, DC)
    return np.einsum('ic,id->cd', M, M.conj())
def _holevo_R(psi, rdm):
    br = []
    for Pr in (Pp_T, Pm_T):
        v = np.tensordot(Pr, psi.reshape(DT, 2 * DC), axes=([1], [0])).reshape(-1)
        p = float(np.vdot(v, v).real)
        br.append((p, rdm(v / np.sqrt(p)) if p > 1e-14 else None))
    avg = sum(p * m for p, m in br if m is not None)
    return vn(avg) - sum(p * vn(m) for p, m in br if m is not None)
def info_R_probe(psi): return _holevo_R(psi, probe_rdm)
def info_R_clock(psi): return _holevo_R(psi, clock_rdm)

# --- verify the eigendecomposition propagator against the hand-written expm, once ---
Hv = H_clockspace(LAM)
Uv = evolve_factory(Hv)
_p = full(AVEC, MET["|+> superpos"])
print(f"\n    propagator check: ||eigh-route psi(T=1.3) - expm-route psi(T=1.3)|| = "
      f"{np.linalg.norm(Uv(_p,1.3) - expm(-1j*Hv*1.3)@_p):.3e}")

CFG = [("lam=0            (no metric at all)", dict(lam=0.0)),
       ("lam=%.1f hop only" % LAM,             dict(lam=LAM, sw_pot=0.0, sw_clk=0.0)),
       ("lam=%.1f clock only" % LAM,           dict(lam=LAM, sw_hop=0.0, sw_pot=0.0)),
       ("lam=%.1f pot only" % LAM,             dict(lam=LAM, sw_hop=0.0, sw_clk=0.0)),
       ("lam=%.1f FULL" % LAM,                 dict(lam=LAM))]
TS = [2.0, 3.0, 4.0, 6.0, 8.0, 12.0]

print("\n  (5a) AHARONOV-BOHM FRINGE CONTRAST at the recombination site 4, antipodal to the source.")
print("       V = |P_+(4) - P_-(4)| / (P_+(4) + P_-(4)), where +/- are the two flux sectors of the")
print("       SAME dynamics.  The bare ring is exactly symmetric, so site 4 is a DARK PORT: at")
print("       lam=0 one flux sector arrives there and the other cancels exactly, V = 1 identically.")
print("       Anything that distinguishes the two arms degrades it.  Metric prepared in |+>.")
VIS = {}
for nm, kw in CFG:
    U = evolve_factory(H_clockspace(**kw))
    row = []
    for T in TS:
        pa = U(full(AVEC, MET["|+> superpos"]), T)
        pb = U(full(BVEC, MET["|+> superpos"]), T)
        p, q = psite(pa, 4), psite(pb, 4)
        row.append(abs(p - q) / (p + q) if p + q > 1e-12 else 0.0)
    VIS[nm] = row
    print(f"    {nm:<34s} " + "  ".join(f"{v:7.4f}" for v in row) + f"   | mean {np.mean(row):7.4f}")
print("    " + " " * 34 + " " + "  ".join(f"T={t:<5.1f}" for t in TS))
print("\n       the two arrival probabilities themselves, so the dark port is visible, not asserted:")
for nm, kw in [CFG[0], CFG[2], CFG[-1]]:
    U = evolve_factory(H_clockspace(**kw))
    pr = [(psite(U(full(AVEC, MET["|+> superpos"]), T), 4),
           psite(U(full(BVEC, MET["|+> superpos"]), T), 4)) for T in TS]
    print(f"    {nm:<34s} P+ " + "  ".join(f"{p:7.4f}" for p, _ in pr))
    print(f"    {'':<34s} P- " + "  ".join(f"{q:7.4f}" for _, q in pr))

print("\n  (5b) CLOCK ENTROPY  S(clock) in bits -- how much the clock has recorded.  Starts at 0.")
for nm, kw in CFG:
    U = evolve_factory(H_clockspace(**kw))
    row = [vn(clock_rdm(U(full((AVEC + BVEC) / np.sqrt(2), MET["|+> superpos"]), T))) for T in TS]
    print(f"    {nm:<34s} " + "  ".join(f"{v:7.4f}" for v in row))

print("\n  (5c) I(R : probe position) in bits -- W-37's measurement, re-run with the clock present.")
print("       W-37 (no metric, no clock) reached 0.83 bits.  Ceiling 1.")
for nm, kw in CFG:
    U = evolve_factory(H_clockspace(**kw))
    row = [info_R_probe(U(full((AVEC + BVEC) / np.sqrt(2), MET["|+> superpos"]), T)) for T in TS]
    print(f"    {nm:<34s} " + "  ".join(f"{v:7.4f}" for v in row))

print("\n  (5f) I(R : clock) in bits -- does the CLOCK itself become a reader of the loop?")
for nm, kw in CFG:
    U = evolve_factory(H_clockspace(**kw))
    row = [info_R_clock(U(full((AVEC + BVEC) / np.sqrt(2), MET["|+> superpos"]), T)) for T in TS]
    print(f"    {nm:<34s} " + "  ".join(f"{v:7.4f}" for v in row))

print("\n  (5d) METRIC STATE DEPENDENCE of the fringe contrast, FULL coupling.")
for mn, mv in MET.items():
    U = evolve_factory(H_clockspace(lam=LAM))
    row = []
    for T in TS:
        p, q = psite(U(full(AVEC, mv), T), 4), psite(U(full(BVEC, mv), T), 4)
        row.append(abs(p - q) / (p + q) if p + q > 1e-12 else 0.0)
    print(f"    metric {mn:<27s} " + "  ".join(f"{v:7.4f}" for v in row))

print("\n  (5e) CLOCK-RATE SWEEP in the CLOCK-ONLY configuration -- hop amplitudes untouched, no")
print("       potential on the probe, traceless clock generator, so the ONLY thing in the problem")
print("       is that the clock ticks at a different rate on the two arms.  omega=0 freezes the")
print("       clock and must return the dark port exactly.")
print(f"       {'omega':>7s}   " + "  ".join(f"T={t:<5.1f}" for t in TS) + "   mean    S(clock)@8  I(R:probe)@12")
CLKV = {}
for om in [0.0, 0.15, 0.3, 0.6, 1.2, 2.4, 4.8]:
    U = evolve_factory(H_clockspace(lam=LAM, sw_hop=0.0, sw_pot=0.0, omega=om))
    row = []
    for T in TS:
        p, q = psite(U(full(AVEC, MET["|+> superpos"]), T), 4), psite(U(full(BVEC, MET["|+> superpos"]), T), 4)
        row.append(abs(p - q) / (p + q) if p + q > 1e-12 else 0.0)
    sc = vn(clock_rdm(U(full((AVEC + BVEC) / np.sqrt(2), MET["|+> superpos"]), 8.0)))
    ip = info_R_probe(U(full((AVEC + BVEC) / np.sqrt(2), MET["|+> superpos"]), 12.0))
    CLKV[om] = (float(np.mean(row)), sc, ip)
    print(f"       {om:7.2f}   " + "  ".join(f"{v:7.4f}" for v in row)
          + f"  {np.mean(row):7.4f}  {sc:9.4f}  {ip:11.4f}")

print("\n  CONTROLS.")
Uc0 = evolve_factory(H_clockspace(lam=0.0, cut=True))
UcL = evolve_factory(H_clockspace(lam=LAM, cut=True))
print("    CUT RING (one hop removed: every local interaction kept, only the closed path gone).")
print("    If the probe's reading of R is topological it must be exactly zero, clock or no clock.")
for T in [4.0, 8.0, 12.0]:
    a = info_R_probe(Uc0(full((AVEC + BVEC) / np.sqrt(2), MET["|+> superpos"]), T))
    b = info_R_probe(UcL(full((AVEC + BVEC) / np.sqrt(2), MET["|+> superpos"]), T))
    pa, pb = psite(Uc0(full(AVEC, MET["|+> superpos"]), T), 4), psite(Uc0(full(BVEC, MET["|+> superpos"]), T), 4)
    va = abs(pa - pb) / (pa + pb) if pa + pb > 1e-12 else 0.0
    print(f"      T={T:5.1f}   I(R:probe) cut,lam=0 = {a:.3e}   cut,lam={LAM} = {b:.3e}   "
          f"fringe contrast cut,lam=0 = {va:.3e}")
U0 = evolve_factory(H_clockspace(lam=0.0))
print("    lam=0 must leave the clock EXACTLY unentangled (H_clock = omega*I(x)I(x)C commutes with all):")
for T in [4.0, 12.0]:
    print(f"      T={T:5.1f}   S(clock) at lam=0 = "
          f"{vn(clock_rdm(U0(full((AVEC+BVEC)/np.sqrt(2), MET['|+> superpos']), T))):.3e}")

print("    W-37 REPLICA: the bare 128-dim transport carrier, no metric factor and no clock factor")
print("    present at all.  lam=0 with the clock attached must reproduce it to machine precision.")
Hw37 = np.zeros((DT, DT), complex)
for k in range(8):
    Hw37 -= TAU * (HOP[k] + HOP[k].conj().T)
evw, Vw = np.linalg.eigh(Hw37)
def _rdm37(v):
    G = np.outer(v, v.conj()); rho = np.zeros((8, 8), complex)
    np.add.at(rho, (POS[PAIRS[:, 0]], POS[PAIRS[:, 1]]), G[PAIRS[:, 0], PAIRS[:, 1]])
    return rho
def _info37(v):
    br = []
    for Pr in (Pp_T, Pm_T):
        u = Pr @ v; p = float(np.vdot(u, u).real)
        br.append((p, _rdm37(u / np.sqrt(p)) if p > 1e-14 else None))
    avg = sum(p * m for p, m in br if m is not None)
    return vn(avg) - sum(p * vn(m) for p, m in br if m is not None)
psi37 = (AVEC + BVEC) / np.sqrt(2)
dmax = 0.0; LAM0_INFO = {}
print(f"      {'T':>6s} {'I(R:probe) W-37 replica':>24s} {'I(R:probe) lam=0 w/ clock':>26s}")
for T in TS:
    a = _info37(Vw @ (np.exp(-1j * evw * T) * (Vw.conj().T @ psi37)))
    b = info_R_probe(U0(full(psi37, MET["|+> superpos"]), T))
    dmax = max(dmax, abs(a - b)); LAM0_INFO[T] = b
    print(f"      {T:6.1f} {a:24.6f} {b:26.6f}")
print(f"      max |difference| = {dmax:.3e}")

print("\n" + "=" * 108)
print("6.  DOES THE METRIC DO WORK?  Every number below is the same quantity with the metric")
print("    coupling off and on.  If a row does not move, the metric did nothing to it.")
print("=" * 108)
gO = RED[("metric ON  (lam=%.1f)" % LAM, "gauge")]; gF = RED[("metric OFF (lam=0)", "gauge")]
mO = RED[("metric ON  (lam=%.1f)" % LAM, "metric")]; mF = RED[("metric OFF (lam=0)", "metric")]
print(f"    {'quantity':<44s} {'lam=0':>16s} {'lam=%.1f'%LAM:>16s}")
print("    " + "-" * 78)
print(f"    {'I(gauge:|F|=1) bits':<44s} {gF[0]:16.6f} {gO[0]:16.6f}")
print(f"    {'I(gauge:all 6) bits':<44s} {gF[2]:16.6f} {gO[2]:16.6f}")
print(f"    {'I(metric:all 6) bits':<44s} {mF[2]:16.6f} {mO[2]:16.6f}")
print(f"    {'sieve rate  RIM LOOP (x) I_M':<44s} {d0['RIM LOOP(all 4)(x)I_M']:16.6e} {dL['RIM LOOP(all 4)(x)I_M']:16.6e}")
print(f"    {'sieve rate  I_gauge (x) sz_M':<44s} {d0['I_gauge(x)sz_M']:16.6e} {dL['I_gauge(x)sz_M']:16.6e}")
print(f"    {'AB fringe contrast, mean over the T grid':<44s} {np.mean(VIS[CFG[0][0]]):16.6f} {np.mean(VIS[CFG[-1][0]]):16.6f}")
print(f"    {'  same, clock-only configuration':<44s} {CLKV[0.0][0]:16.6f} {CLKV[1.2][0]:16.6f}   (omega=0 vs omega=1.2)")
print(f"    {'I(R : probe position) at T=12':<44s} {LAM0_INFO[12.0]:16.6f} "
      f"{info_R_probe(Uv(full((AVEC+BVEC)/np.sqrt(2), MET['|+> superpos']), 12.0)):16.6f}")
print(f"    {'S(clock) at T=8':<44s} "
      f"{vn(clock_rdm(U0(full((AVEC+BVEC)/np.sqrt(2), MET['|+> superpos']), 8.0))):16.6f} "
      f"{vn(clock_rdm(Uv(full((AVEC+BVEC)/np.sqrt(2), MET['|+> superpos']), 8.0))):16.6f}")
print("\ndone.")
