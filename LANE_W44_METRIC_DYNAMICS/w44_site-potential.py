"""W-44.  A METRIC DEGREE OF FREEDOM THAT ACTUALLY DOES WORK.

THE GAP.  The carrier is pure gauge structure: links with no length, no proper time, no clock.
W-43b tried to add one and failed in a specific, nameable way -- it gave a site a two-level "mass",
set the system Hamiltonian to ZERO, and coupled the environment straight to the mass operator. The
potential that mass would source never appeared in H, so nothing it did could be gravitational; and
an environment coupled directly to the record is a projector onto the answer.

WHAT IS BUILT HERE.  A two-level mass at ONE vertex of the 3x3 patch, with OCCUPATION n in {0,1}
("a mass is here" / "no mass is here"), sourcing a NEWTONIAN POTENTIAL Phi(x) = 1/(|x-x_mass|+a)
over the patch. The potential enters as a LAPSE / REDSHIFT factor multiplying the LOCAL Hamiltonian
density:

      H  =  - sum_p  N_p  (B_p + B_p^dag)  -  g2 sum_k N_k (Z_k + Z_k^dag)  -  eps sigma_x
      N_p = 1 + lambda * f_p * n            f_p = Phi(plaquette p) / max Phi
      N_k = 1 + lambda * f_k * n            f_k = Phi(link k)      / max Phi

That is exactly what a metric does: g_00 at a point rescales how fast the local dynamics runs, so
phase accrues with TIME SPENT and is EVEN under path reversal, unlike a gauge phase. lambda is the
metric coupling. lambda = 0 removes the potential from H entirely and the carrier is the W-34/W-36
carrier tensored with a decoupled spin. eps is a mass tunnelling amplitude: it is what makes the
metric variable a DYNAMICAL degree of freedom rather than a conserved label.

In the transport ring the same mass sources a diagonal SITE POTENTIAL V * w_v * n on the probe's
position -- phase per unit dwell time, even under reversal. Section (c) runs the OCCUPATION coupling
and, beside it, the signed-mass coupling V * w_v * sigma_z, which returns an EXACT zero for a reason
that is itself a measurement of what "even under reversal" costs.

NOTHING IN THIS FILE COUPLES AN ENVIRONMENT, A BATH, OR A PROBE TO THE MASS OPERATOR. Every jump
operator is Z on a gauge link; every environment qubit reads a gauge link; the transport probe hops
gauge links. The mass is reachable ONLY through the potential it sources. That is the whole design.

ROUTES USED (no scipy anywhere; expm is scaling-and-squaring, written below):
  SELECTION  : EIGENDECOMPOSITION of the Lindblad generator. system dim 32 -> generator 1024 x 1024.
               No exponentiation of the generator is performed.
  REDUNDANCY : EXACT UNITARY evolution of the joint pure state. 32 system x 2^5 environment = 1024.
  TRANSPORT  : EXACT UNITARY evolution. 128 (probe+gauge) x 2 (mass) = 256.
  Largest object handled: 1024 x 1024.

FORCED-OR-NOT, DECLARED BEFORE ANY DYNAMICS RUNS -- see the block printed at the top of the output.
"""
import itertools, numpy as np

# ----------------------------------------------------------------------------------------------
# linear algebra (numpy only)
# ----------------------------------------------------------------------------------------------
def expm(A):
    """scaling and squaring with a Taylor series"""
    nr = np.linalg.norm(A, np.inf)
    k = max(0, int(np.ceil(np.log2(nr))) + 1) if nr > 0 else 0
    B = A / (2.0 ** k)
    X = np.eye(A.shape[0], dtype=complex); T = X.copy()
    for m in range(1, 80):
        T = T @ B / m; X = X + T
        if np.linalg.norm(T, np.inf) < 1e-18 * max(1.0, np.linalg.norm(X, np.inf)): break
    for _ in range(k): X = X @ X
    return X

def vn(rho):
    ev = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    ev = ev[ev > 1e-12]
    return float(-(ev * np.log2(ev)).sum())

def audit(name, O):
    """REQUIREMENT 5. No operator gets a commutator read until it has earned one."""
    fro = np.linalg.norm(O)
    spec = np.linalg.norm(O, 2)
    ud = np.linalg.norm(O @ O.conj().T - np.eye(O.shape[0]))
    ev = np.unique(np.round(np.linalg.eigvals(O), 6))
    print(f"    {name:<34s} ||O||_F={fro:9.4f}  ||O||_2={spec:7.4f}  "
          f"unitarity defect={ud:9.2e}  distinct eigenvalues={len(ev)}")
    return len(ev), fro

# ----------------------------------------------------------------------------------------------
# the carrier: 3x3 planar patch, Z_2 on links, Gauss law at every vertex   (identical to W-34/W-36)
# ----------------------------------------------------------------------------------------------
V2 = [(i, j) for j in range(3) for i in range(3)]
vid = {v: k for k, v in enumerate(V2)}
E = []
for j in range(3):
    for i in range(2): E.append((vid[(i, j)], vid[(i + 1, j)]))
for j in range(2):
    for i in range(3): E.append((vid[(i, j)], vid[(i, j + 1)]))
L = len(E)
hid = lambda i, j: j * 2 + i
vx  = lambda i, j: 6 + j * 3 + i
P = [[(hid(i, j), +1), (vx(i + 1, j), +1), (hid(i, j + 1), -1), (vx(i, j), -1)]
     for j in range(2) for i in range(2)]
PXY = [(i + 0.5, j + 0.5) for j in range(2) for i in range(2)]          # plaquette centres
CENTER = vid[(1, 1)]
CUT   = [k for k, (a, b) in enumerate(E) if a == CENTER or b == CENTER]
PERIM = [k for k in range(L) if k not in CUT]

def build(N=2):
    st = [s for s in itertools.product(range(N), repeat=L)
          if all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
                - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % N == 0 for v in range(9))]
    return st, {s: i for i, s in enumerate(st)}
st, idx = build()
DS = len(st)

def Zop(links):
    return np.diag([(-1.0) ** (sum(s[k] for k in links) % 2) for s in st]).astype(complex)
def Move(mv):
    M = np.zeros((DS, DS), complex)
    for j, s in enumerate(st):
        t = list(s)
        for k, sg in mv: t[k] = (t[k] + sg) % 2
        t = tuple(t)
        if t in idx: M[idx[t], j] = 1.0
    return M
def compose(ps):
    acc = {}
    for p in ps:
        for k, sg in p: acc[k] = acc.get(k, 0) + sg
    return [(k, s) for k, s in acc.items() if s % 2]

# link midpoints, for the potential
LXY = []
for j in range(3):
    for i in range(2): LXY.append((i + 0.5, float(j)))
for j in range(2):
    for i in range(3): LXY.append((float(i), j + 0.5))

# ----------------------------------------------------------------------------------------------
# the metric degree of freedom
# ----------------------------------------------------------------------------------------------
MASS_XY = (0.0, 0.0)                      # the mass sits at the corner vertex (0,0)
sz = np.array([[1, 0], [0, -1]], complex)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
I2 = np.eye(2, dtype=complex)
IdS = np.eye(DS, dtype=complex)
NHAT = (I2 - sz) / 2          # the mass OCCUPATION, eigenvalues {0,1}. A mass is present or it is
                              # not; there is no negative mass. Using the polarisation sigma_z
                              # instead installs an accidental antiunitary symmetry -- measured and
                              # reported in section (c).

SOFT = 0.5                                # finite mass size; keeps Phi finite at its own vertex
def newton(xy):
    d = np.hypot(xy[0] - MASS_XY[0], xy[1] - MASS_XY[1])
    return 1.0 / (d + SOFT)
FP = np.array([newton(c) for c in PXY]); FP = FP / FP.max()       # lapse profile, plaquettes
FL = np.array([newton(c) for c in LXY]); FL = FL / FL.max()       # lapse profile, links

BP   = [(lambda X: X + X.conj().T)(Move(p)) for p in P]           # magnetic, plaquette by plaquette
ZK   = [(lambda X: X + X.conj().T)(Zop([k])) for k in range(L)]   # electric, link by link
MAG  = sum(BP)
ELEC = sum(ZK)
RG   = Move(compose(P))                                           # rim Wilson loop = the GAUGE record

def Hsys(lam, eps, g2):
    """lam = metric coupling. lam=0 -> H = H_gauge (x) 1 - eps (1 (x) sigma_x), a decoupled spin."""
    H = np.zeros((DS * 2, DS * 2), complex)
    for a, B in enumerate(BP):
        H += -np.kron(B, I2 + lam * FP[a] * NHAT)
    if g2:
        for k, Z in enumerate(ZK):
            H += -g2 * np.kron(Z, I2 + lam * FL[k] * NHAT)
    if eps:
        H += -eps * np.kron(IdS, sx)
    return H

RG32 = np.kron(RG, I2)                    # gauge record on the enlarged carrier
M32  = np.kron(IdS, sz)                   # metric record
Id32 = np.eye(DS * 2, dtype=complex)

# ==============================================================================================
print("=" * 100)
print("W-44   METRIC / PROPER TIME AS A DEGREE OF FREEDOM THAT ENTERS THE DYNAMICS")
print("=" * 100)
print(f"  carrier: 3x3 patch, Z_2 links, Gauss at every vertex -> gauge dim {DS}")
print(f"  metric : one two-level mass, occupation n in [0,1], at vertex {MASS_XY};\n           it sources Phi(x) = 1/(|x-x_mass| + {SOFT}) and that Phi is the LAPSE multiplying\n           the local Hamiltonian density")
print(f"  lapse profile over the 4 plaquettes  f_p = {np.round(FP,4)}")
print(f"  lapse profile over the 12 links      f_k = {np.round(FL,4)}")
print(f"  bath / environment / probe couple ONLY to gauge links. Nothing couples to the mass.")

print()
print("-" * 100)
print("  FORCED-OR-NOT, DECLARED BEFORE ANY DYNAMICS RUNS")
print("-" * 100)
print("""  AT RISK #1 -- the SIEVE naming the metric operator M = 1 (x) sigma_z a protected record.
     FORCED whenever eps = 0, for ANY lambda, by two commutators:
        [M, H] = 0   because every lambda-term is diagonal in the mass basis,
        [M, L_k] = 0 because every jump operator is Z on a gauge link.
     Then Gdag(M) = 0 exactly: M is a zero-rate mode by construction, not by dynamics. Worse, Gdag
     then COMMUTES with left-multiplication by M, so the ENTIRE spectrum is exactly doubled and the
     sieve cannot separate O from M.O at all.
     SETTLED BY: giving the mass a tunnelling term -eps sigma_x. Then [M,H] != 0 and M must earn its
     rate. EVERY reported sieve number below uses eps > 0. The eps = 0 doubling is printed too, as a
     check that the argument is right rather than a claim.
     A SECOND forced case, same family: at lambda = 0 the mass factor evolves UNITARILY and
     independently of both gauge field and bath, so EVERY mass operator has strictly zero decay rate
     there. So the lambda = 0 sieve necessarily crowns a metric operator, and that crowning means
     nothing. Only lambda > 0 AND eps > 0 is an informative sieve. Both degenerate cases are run and
     printed anyway, because a declared forced case that does not come out forced is a broken
     argument, not a discovery.

  AT RISK #2 -- the REDUNDANCY of the metric variable.
     W-43b forced it to ~1 bit by coupling the environment straight to sigma_z^mass. Here the
     environment couples only to Z on cut links, and [Z_cut (x) 1, 1 (x) sigma_z] = 0 with a product
     initial state, so at lambda = 0 the mass factorises out of the joint state exactly and
     I(M:F) = 0 IDENTICALLY for every fragment. Any nonzero I(M:F) is therefore transmitted by the
     potential's effect on the gauge dynamics and by nothing else.

  NOT AT RISK -- I(gauge:F). W-36 already measured 0.047 with this same machinery; lambda can only
     move it.""")

# ==============================================================================================
print()
print("-" * 100)
print("  OPERATORS EARN MEASUREMENT (requirement 5)")
print("-" * 100)
audit("R_gauge = rim loop (x) 1", RG32)
audit("M_metric = 1 (x) sigma_z", M32)
audit("R_gauge . M_metric", RG32 @ M32)
audit("single plaquette (x) 1", np.kron(BP[0] / 2, I2))
print(f"    [R_gauge , M_metric] norm = {np.linalg.norm(RG32@M32-M32@RG32):.2e}   (they are compatible records)")

# ==============================================================================================
print()
print("-" * 100)
print("  DOES THE METRIC DO WORK?  (requirement 1, checked directly on the evolution)")
print("-" * 100)
LAM = 1.6
H0 = Hsys(0.0, 0.0, 0.0); H1 = Hsys(LAM, 0.0, 0.0)
print(f"    || H(lambda={LAM}) - H(lambda=0) ||  = {np.linalg.norm(H1-H0):.6f}")
# evolve the SAME gauge state in the two mass sectors and compare
g = np.random.default_rng(11)
w = g.normal(size=DS) + 1j * g.normal(size=DS); w /= np.linalg.norm(w)
for T in (1.0, 3.0):
    ups = []
    for msign in (0, 1):                       # |0> = sigma_z +1 , |1> = sigma_z -1
        v = np.zeros(DS * 2, complex); v[msign::2] = w
        v = expm(-1j * H1 * T) @ v
        ups.append(v[msign::2])
    ov = abs(np.vdot(ups[0], ups[1]))
    v0 = expm(-1j * H0 * T) @ np.concatenate([[a, 0] for a in w])
    print(f"    T={T:4.1f}   |<psi_gauge(n=0) | psi_gauge(n=1)>| = {ov:.6f}   "
          f"(1.000000 would mean the mass changes nothing)")

# ==============================================================================================
# (a)  REDUNDANCY -- explicit local environment, same machinery for both records
# ==============================================================================================
print()
print("=" * 100)
print("  (a) REDUNDANCY.  Holevo information a fragment of a LOCAL environment holds about")
print("      the GAUGE record, and about the METRIC record. Same environment, same coupling,")
print("      same estimator. Ceiling 1 bit each. Route: exact unitary evolution, dim 32 x 32 = 1024.")
print("=" * 100)

NQ = 5
DE = 2 ** NQ
plus = np.ones(2, complex) / np.sqrt(2.0)
def kron_list(ops):
    o = np.array([[1]], complex)
    for x in ops: o = np.kron(o, x)
    return o

NEAR = [0, vx(0, 0)]                      # the two links touching the mass's own vertex (0,0)

def joint_state(lam, kappa, T, seed=3, links=None):
    """system (32) (x) NQ env qubits. env couples to Z on gauge links only, never to the mass.
       eps = 0 and g2 = 0 here so BOTH records are exactly conserved -- the QND setting of W-36,
       which is what makes projecting the final state equal to conditioning on the initial value."""
    links = CUT if links is None else links
    H = np.kron(Hsys(lam, 0.0, 0.0), np.eye(DE, dtype=complex))
    for q in range(NQ):
        Zk = np.kron(Zop([links[q % len(links)]]), I2)
        ops = [I2] * NQ; ops[q] = sz
        H = H + kappa * np.kron(Zk, kron_list(ops))
    U = expm(-1j * H * T)
    gg = np.random.default_rng(seed)
    w = gg.normal(size=DS) + 1j * gg.normal(size=DS)
    Pp = (IdS + RG) / 2; Pm = (IdS - RG) / 2
    a = Pp @ w; b = Pm @ w; a /= np.linalg.norm(a); b /= np.linalg.norm(b)
    psiG = (a + b) / np.sqrt(2.0); psiG /= np.linalg.norm(psiG)        # <R_gauge> = 0 exactly
    psiS = np.kron(psiG, plus)                                          # <M_metric> = 0 exactly
    psiE = kron_list([plus.reshape(2, 1)] * NQ).reshape(-1)
    psi = U @ np.kron(psiS, psiE)
    return psi.reshape((DS * 2,) + (2,) * NQ)

def holevo(psiT, frag, Rop):
    Pp = (Id32 + Rop) / 2; Pm = (Id32 - Rop) / 2
    br = []
    for Proj in (Pp, Pm):
        v = np.tensordot(Proj, psiT, axes=([1], [0]))
        p = float(np.vdot(v, v).real)
        if p < 1e-14: br.append((0.0, None)); continue
        v = v / np.sqrt(p)
        keep = [0] + [1 + i for i in frag]
        tr = [ax for ax in range(1 + NQ) if ax not in keep]
        Mx = np.transpose(v, keep + tr).reshape(DS * 2 * 2 ** len(frag), -1)
        rho = (Mx @ Mx.conj().T).reshape(DS * 2, 2 ** len(frag), DS * 2, 2 ** len(frag))
        br.append((p, np.einsum('ijik->jk', rho)))
    avg = sum(p * m for p, m in br if m is not None)
    return vn(avg) - sum(p * vn(m) for p, m in br if m is not None)

def profile(psiT, Rop):
    out = []
    for f in range(NQ + 1):
        combos = list(itertools.combinations(range(NQ), f))[:20]
        out.append(float(np.mean([holevo(psiT, c, Rop) for c in combos])))
    return out

print("\n  SATURATION FIRST (as in W-36): scan the environment coupling until the WHOLE environment")
print("  holds the gauge bit; a fragment profile means nothing before that. Then use the SAME")
print(f"  (kappa,T) for the metric readout. lambda = {LAM}.")
print(f"    {'kappa':>7s} {'T':>6s} {'I(gauge:all 5)':>16s} {'I(metric:all 5)':>17s}")
print("    " + "-" * 50)
best = None
for kap, T in [(0.6, 6.0), (1.2, 6.0), (2.0, 8.0), (3.0, 8.0), (3.0, 12.0), (4.0, 10.0),
               (5.0, 12.0), (6.0, 14.0), (8.0, 16.0), (12.0, 20.0)]:
    psiT = joint_state(LAM, kap, T)
    ig = holevo(psiT, tuple(range(NQ)), RG32)
    im = holevo(psiT, tuple(range(NQ)), M32)
    print(f"    {kap:7.2f} {T:6.1f} {ig:16.6f} {im:17.6f}")
    if best is None or ig > best[2]: best = (kap, T, ig)
KAP, TT = best[0], best[1]
print(f"    -> using kappa={KAP}, T={TT}")

rows = []
for tag, lam, lk, lknm in [("METRIC ON  lambda=%.2f" % LAM, LAM, CUT, "CUT"),
                           ("CONTROL    lambda=0", 0.0, CUT, "CUT"),
                           ("METRIC ON, env AT the mass", LAM, NEAR, "the 2 links at the mass")]:
    psiT = joint_state(lam, KAP, TT, links=lk)
    pg = profile(psiT, RG32); pm = profile(psiT, M32)
    print(f"\n  {tag}    (n={NQ} env qubits on {lknm}, kappa={KAP}, T={TT})")
    print("      |F|            : " + "  ".join(f"{i:8d}" for i in range(NQ + 1)))
    print("      I(gauge  : F)  : " + "  ".join(f"{v:8.5f}" for v in pg))
    print("      I(metric : F)  : " + "  ".join(f"{v:8.5f}" for v in pm))
    rows.append((tag, lam, pg, pm))

print("\n  READING. Redundancy = a SINGLE fragment already holds most of the bit.")
print(f"  {'case':>24s} {'record':>8s} {'I(|F|=1)':>10s} {'I(|F|=3)':>10s} {'I(all)':>9s} {'I(1)/I(all)':>12s}")
print("  " + "-" * 80)
RED = {}
for tag, lam, pg, pm in rows:
    for nm, v in (("gauge", pg), ("metric", pm)):
        tot = v[-1]
        r = v[1] / tot if abs(tot) > 1e-9 else float('nan')
        RED[(tag, nm)] = (v[1], v[3], tot, r)
        print(f"  {tag:>24s} {nm:>8s} {v[1]:10.5f} {v[3]:10.5f} {tot:9.5f} {r:12.4f}")

_, _, pg0, pm0 = rows[1]
print(f"\n  CONTROL CHECK -- lambda = 0 must give EXACTLY the no-metric answer:")
print(f"      max |I(metric:F)| over all fragment sizes at lambda=0 : {max(abs(x) for x in pm0):.3e}")
# and the gauge profile at lambda=0 must equal the profile computed with NO mass factor at all
def joint_state_nomass(kappa, T, seed=3):
    H = np.kron(-MAG, np.eye(DE, dtype=complex))
    for q in range(NQ):
        Zk = Zop([CUT[q % len(CUT)]])
        ops = [I2] * NQ; ops[q] = sz
        H = H + kappa * np.kron(Zk, kron_list(ops))
    U = expm(-1j * H * T)
    gg = np.random.default_rng(seed)
    w = gg.normal(size=DS) + 1j * gg.normal(size=DS)
    Pp = (IdS + RG) / 2; Pm = (IdS - RG) / 2
    a = Pp @ w; b = Pm @ w; a /= np.linalg.norm(a); b /= np.linalg.norm(b)
    psiG = (a + b) / np.sqrt(2.0); psiG /= np.linalg.norm(psiG)
    psi = U @ np.kron(psiG, kron_list([plus.reshape(2, 1)] * NQ).reshape(-1))
    return psi.reshape((DS,) + (2,) * NQ)
def holevo_nomass(psiT, frag):
    Pp = (IdS + RG) / 2; Pm = (IdS - RG) / 2
    br = []
    for Proj in (Pp, Pm):
        v = np.tensordot(Proj, psiT, axes=([1], [0])); p = float(np.vdot(v, v).real)
        if p < 1e-14: br.append((0.0, None)); continue
        v = v / np.sqrt(p)
        keep = [0] + [1 + i for i in frag]; tr = [ax for ax in range(1 + NQ) if ax not in keep]
        Mx = np.transpose(v, keep + tr).reshape(DS * 2 ** len(frag), -1)
        rho = (Mx @ Mx.conj().T).reshape(DS, 2 ** len(frag), DS, 2 ** len(frag))
        br.append((p, np.einsum('ijik->jk', rho)))
    avg = sum(p * m for p, m in br if m is not None)
    return vn(avg) - sum(p * vn(m) for p, m in br if m is not None)
psiNM = joint_state_nomass(KAP, TT)
pgNM = [float(np.mean([holevo_nomass(psiNM, c)
                       for c in list(itertools.combinations(range(NQ), f))[:20]])) for f in range(NQ + 1)]
print(f"      I(gauge:F) at lambda=0 vs the carrier with NO mass factor at all, max difference : "
      f"{max(abs(a-b) for a,b in zip(pg0,pgNM)):.3e}")

# ==============================================================================================
# (b)  SELECTION -- the predictability sieve with the metric term ON. Nominate nothing.
# ==============================================================================================
print()
print("=" * 100)
print("  (b) SELECTION.  Predictability sieve: eigendecomposition of the Lindblad generator,")
print("      slowest left-eigenmodes = the observables the dynamics protects. Nothing nominated.")
print("      Route: EIGENDECOMPOSITION of the 1024 x 1024 generator (system dim 32). No expm.")
print("=" * 100)

NAMES = {}; OPS = {}
for r in range(5):
    for S in itertools.combinations(range(4), r):
        G = Move(compose([P[i] for i in S])) if S else IdS.copy()
        gname = "identity" if not S else ("RIM LOOP (all 4)" if r == 4 else
                (f"plaquette {S}" if r == 1 else f"{r}-plaquette loop {S}"))
        for mn, mo in (("1", I2), ("sz", sz), ("sx", sx), ("sy", sy)):
            O = np.kron(G, mo)
            key = (S, mn)
            NAMES[key] = f"{gname} (x) {mn}"
            OPS[key] = O / np.linalg.norm(O)

TGRID = np.logspace(-3.0, 4.0, 3000)

def survival_rates(w, U, opdict):
    """ESTIMATOR-FREE protection score, and rotation-blind.
       O(t) = U exp(w t) U^-1 O, so the autocorrelation is C(t) = sum_i a_i exp(w_i t) with
       sum_i a_i = 1. Oscillation lives in Im(w_i) and says nothing about protection, so score the
       AMPLITUDE survival  S(t) = sum_i |a_i| exp(-Gamma_i t) / sum_i |a_i|, which is monotone.
       rate = 1 / t(S = 1/e); an operator that never decays reports rate 0 exactly."""
    keys = list(opdict.keys())
    Vs = np.stack([opdict[k].reshape(-1) for k in keys], axis=1)
    C = np.linalg.solve(U, Vs)                          # coefficients in the left-eigenbasis
    W = (U.conj().T @ Vs)                               # <mode_i , O>
    A = np.conj(W) * C                                  # wrong pairing would break sum a_i = 1
    A = A / np.sum(A, axis=0, keepdims=True)
    resid = np.abs(np.sum(A, axis=0) - 1.0).max()
    Gam = -np.real(w)
    Amp = np.abs(A); Amp = Amp / Amp.sum(axis=0, keepdims=True)
    S = np.exp(-np.outer(TGRID, Gam)) @ Amp             # (nt, nops)
    out = {}
    for j, k in enumerate(keys):
        s = S[:, j]
        hit = np.nonzero(s <= np.exp(-1.0))[0]
        out[k] = 0.0 if len(hit) == 0 else 1.0 / TGRID[hit[0]]
    return out, resid

def sieve(lam, eps, g2, gam=0.5, links=None, topk=8, label=""):
    links = CUT if links is None else links
    H = Hsys(lam, eps, g2)
    Ls = [np.kron(Zop([k]), I2) for k in links]
    Gm = -1j * (np.kron(H, Id32) - np.kron(Id32, H.T))
    for Lj in Ls: Gm += gam * (np.kron(Lj, Lj.conj()) - np.kron(Id32, Id32))
    w, U = np.linalg.eig(Gm.conj().T)
    rate = -np.real(w)
    order = np.argsort(rate)
    print(f"\n  {label}   lambda={lam}  eps={eps}  g2={g2}  gamma={gam}  bath on "
          f"{'CUT' if links is CUT else 'RIM'}")
    print(f"    {'rate':>13s} {'|Im|':>11s}  best operator match")
    print("    " + "-" * 66)
    shown = 0
    for i in order:
        if shown >= topk: break
        O = U[:, i].reshape(DS * 2, DS * 2)
        n = np.linalg.norm(O)
        if n < 1e-12: continue
        O = O / n
        best = max(OPS.items(), key=lambda kv: abs(np.vdot(kv[1].reshape(-1), O.reshape(-1))))
        ov = abs(np.vdot(best[1].reshape(-1), O.reshape(-1)))
        tag = NAMES[best[0]] if ov > 0.3 else "mixed / no clean match"
        print(f"    {rate[i]:13.6e} {abs(np.imag(w[i])):11.4e}  {tag}   (overlap {ov:.3f})")
        shown += 1
    # two rankings, both over the same 63 named operators (the identity is excluded: it is the
    # trace, conserved for free, and is not a candidate record)
    cand = {k: v for k, v in OPS.items() if k != ((), "1")}
    sr, resid = survival_rates(w, U, cand)
    Un = U / np.linalg.norm(U, axis=0)
    rk = []
    for key, O in cand.items():
        ov = np.abs(Un.conj().T @ O.reshape(-1)); ov = ov / max(ov.sum(), 1e-30)
        rk.append((sr[key], float((ov * rate).sum()), NAMES[key], key))
    rk.sort(key=lambda r: (r[0], r[1]))
    print(f"    (eigenbasis reconstruction residual max|sum a_i - 1| = {resid:.2e})")
    return rk

print("\n  REFERENCE. The no-metric sieve on the bare 16-dim carrier, so the lambda=0 run has")
print("  something exact to reproduce.")
def sieve_bare(g2=0.01, gam=0.5):
    H = -MAG - g2 * ELEC
    Gm = -1j * (np.kron(H, IdS) - np.kron(IdS, H.T))
    for k in CUT:
        Lj = Zop([k]); Gm += gam * (np.kron(Lj, Lj.conj()) - np.kron(IdS, IdS))
    w, U = np.linalg.eig(Gm.conj().T); rate = -np.real(w)
    Un = U / np.linalg.norm(U, axis=0)
    cand = {}
    for r in range(1, 5):
        for S in itertools.combinations(range(4), r):
            G = Move(compose([P[i] for i in S])); cand[S] = G / np.linalg.norm(G)
    sr, resid = survival_rates(w, U, cand)
    rk = []
    for S, G in cand.items():
        ov = np.abs(Un.conj().T @ G.reshape(-1)); ov = ov / max(ov.sum(), 1e-30)
        nm = ("RIM LOOP (all 4)" if len(S) == 4 else
              (f"plaquette {S}" if len(S) == 1 else f"{len(S)}-plaquette loop {S}"))
        rk.append((sr[S], float((ov * rate).sum()), nm))
    rk.sort(key=lambda r: (r[0], r[1]))
    print(f"      {'1/e rate':>13s} {'W-34 mean rate':>16s}   operator")
    for a, b, nm in rk[:5]: print(f"      {a:13.6e} {b:16.6e}   {nm}")
    print(f"      margin (2nd / 1st), 1/e rate = {rk[1][0]/rk[0][0]:.2f}x ; "
          f"W-34 mean rate = {rk[1][1]/rk[0][1]:.2f}x")
    return rk
bare = sieve_bare()

print("\n  EXACT CONTROL (requirement 6). Two things must hold at lambda = 0 and must FAIL at")
print("  lambda > 0, or the metric does no work:")
def gen32(lam, eps, g2, gam=0.5):
    H = Hsys(lam, eps, g2)
    Gm = -1j * (np.kron(H, Id32) - np.kron(Id32, H.T))
    for k in CUT:
        Lj = np.kron(Zop([k]), I2); Gm += gam * (np.kron(Lj, Lj.conj()) - np.kron(Id32, Id32))
    return Gm
def gen16(g2, gam=0.5):
    H = -MAG - g2 * ELEC
    Gm = -1j * (np.kron(H, IdS) - np.kron(IdS, H.T))
    for k in CUT:
        Lj = Zop([k]); Gm += gam * (np.kron(Lj, Lj.conj()) - np.kron(IdS, IdS))
    return Gm
def specsort(v, nd=7):
    v = np.round(v, nd)
    return v[np.lexsort((np.imag(v), np.real(v)))]
e16 = specsort(np.repeat(np.linalg.eigvals(gen16(0.01)), 4))
for lam in (0.0, LAM):
    e = specsort(np.linalg.eigvals(gen32(lam, 0.0, 0.01)))
    print(f"      lambda={lam:4.2f}   max |spectrum(32-dim) - spectrum(16-dim) repeated x4| = "
          f"{np.abs(e-e16).max():.3e}")
print("      (the lambda=0 line is the control and must be ~0; the lambda>0 line must not be)")
print("  and the DECLARED DOUBLING at eps=0: vec(M.X) = (M (x) 1) vec(X) row-major, so the metric")
print("  degeneracy is exactly the statement that left-multiplication by M commutes with Gdag.")
LM = np.kron(M32, Id32)
for lam, eps in ((LAM, 0.0), (LAM, 0.15)):
    Gd = gen32(lam, eps, 0.01).conj().T
    print(f"      lambda={lam:4.2f} eps={eps:4.2f}   || [ Gdag , left-mult by M ] || = "
          f"{np.linalg.norm(Gd@LM-LM@Gd):.3e}")

print("\n  FORCED-CASE DEMONSTRATION (eps = 0). The argument above says the whole spectrum must be")
print("  exactly doubled and M must sit at rate 0. If that is what comes out, the argument holds and")
print("  the eps = 0 sieve carries no information about the metric.")
def show(rk, n=8):
    print(f"    ranking of the 63 named operators (lowest = best protected):")
    print(f"      {'1/e rate':>13s} {'W-34 mean rate':>16s}   operator")
    for a, b, nm, _ in rk[:n]: print(f"      {a:13.6e} {b:16.6e}   {nm}")

rk_f = sieve(LAM, 0.0, 0.01, label="FORCED CASE eps=0", topk=6)
show(rk_f, 6)

RUNS = {}
for lab, lam, eps in [("CONTROL lambda=0 (metric decoupled)", 0.0, 0.15),
                      ("METRIC ON lambda=1.6", LAM, 0.15),
                      ("METRIC ON lambda=3.2", 3.2, 0.15)]:
    rk = sieve(lam, eps, 0.01, label=lab, topk=8)
    RUNS[lab] = rk
    show(rk)
    mg = rk[1][0] / rk[0][0] if rk[0][0] > 0 else float('inf')
    print(f"      margin (2nd / 1st), 1/e rate = {mg:.3f}x")

print("\n  COUPLING SCAN. The winner is allowed to change with lambda; print the scan rather than")
print("  one chosen point.  eps = 0.15, g2 = 0.01, bath on the cut, 1/e rates.")
print(f"  {'lambda':>8s} {'winner':>26s} {'rate(RIM(x)1)':>14s} {'rate(1(x)sz)':>14s} {'margin':>9s}")
print("  " + "-" * 78)
import io, contextlib
for lam in [0.0, 0.4, 0.8, 1.2, 1.6, 2.4, 3.2, 4.8, 6.4, 9.6, 14.4]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rk = sieve(lam, 0.15, 0.01, topk=0, label="")
    d = {k: a for a, b, nm, k in rk}
    mg = rk[1][0] / rk[0][0] if rk[0][0] > 0 else float('inf')
    print(f"  {lam:8.2f} {rk[0][2]:>26s} {d[((0,1,2,3),'1')]:14.6e} {d[((),'sz')]:14.6e} {mg:9.3f}")

print("\n  The two rates cross inside the scan. Bisect on rate(1(x)sz) - rate(RIM(x)1), 12 steps,")
print("  at two values of the mass tunnelling eps, so the crossing is not read off one parameter.")
print("  (the 1/e rate is read off a log time grid of 3000 points, so lambda* carries that resolution)")
def gap(lam, eps):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rk = sieve(lam, eps, 0.01, topk=0, label="")
    d = {k: a for a, b, nm, k in rk}
    return d[((), 'sz')] - d[((0, 1, 2, 3), '1')]
print(f"    {'eps':>6s} {'lambda*':>10s} {'gap at lambda*':>16s}")
print("    " + "-" * 36)
for eps in (0.15, 0.30):
    lo, hi = 2.4, 9.6
    if gap(lo, eps) > 0 and gap(hi, eps) < 0:
        for _ in range(12):
            mid = 0.5 * (lo + hi)
            if gap(mid, eps) > 0: lo = mid
            else: hi = mid
        lam_star = 0.5 * (lo + hi)
        print(f"    {eps:6.2f} {lam_star:10.4f} {gap(lam_star, eps):16.3e}")
    else:
        print(f"    {eps:6.2f}   no sign change in [2.4, 9.6]:  gap(2.4) = {gap(lo,eps):+.3e}, "
              f"gap(9.6) = {gap(hi,eps):+.3e}")

print("\n  SIDE BY SIDE -- where the two records sit, run by run.  (1/e rates)")
print(f"  {'run':>36s} {'winner':>26s} {'RIM(x)1':>13s} {'1(x)sz':>13s} {'margin':>9s}")
print("  " + "-" * 104)
for lab, rk in RUNS.items():
    d = {k: a for a, b, nm, k in rk}
    rrim = d[((0, 1, 2, 3), "1")]; rmet = d[((), "sz")]
    mg = rk[1][0] / rk[0][0] if rk[0][0] > 0 else float('inf')
    print(f"  {lab:>36s} {rk[0][2]:>26s} {rrim:13.6e} {rmet:13.6e} {mg:9.3f}")

# ==============================================================================================
# (c)  TRANSPORT -- the decisive locality control, in both directions
# ==============================================================================================
print()
print("=" * 100)
print("  (c) TRANSPORT.  A probe charge hops the 8 rim vertices (W-37's construction) while the")
print("      SAME mass sources a diagonal site potential V*w_v*n on the probe's position, n being")
print("      the mass OCCUPATION. The signed-mass coupling V*w_v*sigma_z is run beside it.")
print("      Route: exact unitary evolution, dim 128 x 2 = 256.")
print("=" * 100)

RIMV = [vid[(0, 0)], vid[(1, 0)], vid[(2, 0)], vid[(2, 1)],
        vid[(2, 2)], vid[(1, 2)], vid[(0, 2)], vid[(0, 1)]]
def link_between(u, v):
    for k, (a, b) in enumerate(E):
        if (a, b) == (u, v) or (b, a) == (u, v): return k
    raise ValueError((u, v))
RIML = [link_between(RIMV[k], RIMV[(k + 1) % 8]) for k in range(8)]
def div(s, v):
    return (sum(s[k] for k, (a, b) in enumerate(E) if a == v)
          - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % 2
STATES = []
for p in range(8):
    for s in itertools.product(range(2), repeat=L):
        q = lambda v: (1 if v == RIMV[p] else 0) ^ (1 if v == CENTER else 0)
        if all(div(s, v) == q(v) for v in range(9)): STATES.append((p, s))
TIDX = {x: i for i, x in enumerate(STATES)}; DT = len(STATES)
def shiftt(s, links):
    t = list(s)
    for k in links: t[k] ^= 1
    return tuple(t)
def RopT():
    Mx = np.zeros((DT, DT), complex)
    for j, (p, s) in enumerate(STATES): Mx[TIDX[(p, shiftt(s, PERIM))], j] = 1.0
    return Mx
def HhopT(tau, skip=None):
    Mx = np.zeros((DT, DT), complex)
    for k in range(8):
        if skip is not None and k == skip: continue
        kn = (k + 1) % 8
        for j, (p, s) in enumerate(STATES):
            if p != k: continue
            i = TIDX.get((kn, shiftt(s, [RIML[k]])))
            if i is not None: Mx[i, j] -= tau
    return Mx + Mx.conj().T
RT = RopT()
WV = np.array([newton(V2[RIMV[p]]) for p in range(8)]); WV = WV / WV.max()
POT = np.diag([WV[p] for p, _ in STATES]).astype(complex)
print(f"      probe ring dim {DT}; with the mass qubit, {DT*2}")
print(f"      site-potential profile over the 8 rim vertices w_v = {np.round(WV,4)}")
RT2 = np.kron(RT, I2); MT2 = np.kron(np.eye(DT, dtype=complex), sz)
print("    operators:")
audit("R_gauge (ring) (x) 1", RT2)
audit("M_metric (ring) = 1 (x) sigma_z", MT2)
Idt2 = np.eye(DT * 2, dtype=complex)

def HT(tau, V, skip=None, coup="n"):
    """coup='n'  : potential sourced by the mass OCCUPATION  (1-sigma_z)/2, eigenvalues {0,1}
       coup='sz' : potential sourced by the signed polarisation sigma_z, eigenvalues {+1,-1}"""
    C = NHAT if coup == "n" else sz
    return np.kron(HhopT(tau, skip), I2) + V * np.kron(POT, C)
for cp in ("n", "sz"):
    Hx = HT(1, 1, coup=cp)
    print(f"    coup={cp:>2s}   ||[R_gauge, H(tau=1,V=1)]|| = "
          f"{np.linalg.norm(RT2@Hx-Hx@RT2):.2e}   ||[M_metric, H]|| = "
          f"{np.linalg.norm(MT2@Hx-Hx@MT2):.2e}")

gg = np.random.default_rng(5)
mask = np.array([1.0 if p == 0 else 0.0 for p, _ in STATES])
wv = (gg.normal(size=DT) + 1j * gg.normal(size=DT)) * mask
Pp = (np.eye(DT) + RT) / 2; Pm = (np.eye(DT) - RT) / 2
a = Pp @ wv; b = Pm @ wv; a /= np.linalg.norm(a); b /= np.linalg.norm(b)
psiG0 = (a + b); psiG0 /= np.linalg.norm(psiG0)
psi0T = np.kron(psiG0, plus)

def probe_rdm(psi):
    """trace out gauge AND mass, keep the probe's 8 positions.
       NOTE, structural: two different probe positions carry DIFFERENT Gauss-law sectors, so no
       gauge configuration is shared between them and this matrix is exactly diagonal. The probe
       reads a record through its position PROBABILITIES and through nothing else."""
    Mx = np.zeros((8, 8), complex)
    v = psi.reshape(DT, 2)
    for i, (p, s) in enumerate(STATES):
        for j, (q, t) in enumerate(STATES):
            if s == t: Mx[p, q] += v[i] @ np.conj(v[j])
    return Mx
def info_probe(psi, Rop):
    Pplus = (Idt2 + Rop) / 2; Pminus = (Idt2 - Rop) / 2
    br = []
    for Pr in (Pplus, Pminus):
        v = Pr @ psi; pr = float(np.vdot(v, v).real)
        br.append((pr, probe_rdm(v / np.sqrt(pr)) if pr > 1e-14 else None))
    avg = sum(p * m for p, m in br if m is not None)
    return vn(avg) - sum(p * vn(m) for p, m in br if m is not None)
print(f"    probe position rdm is exactly diagonal:  max |offdiag| at T=8 = "
      f"{np.abs(probe_rdm(expm(-1j*HT(1.0,1.5)*8.0)@psi0T)-np.diag(np.diag(probe_rdm(expm(-1j*HT(1.0,1.5)*8.0)@psi0T)))).max():.2e}")

print(f"\n  I(record : probe position) in bits.  tau = 1.0, V = 1.5, closed ring.")
print(f"  {'T':>6s} {'I(gauge:probe)':>16s} {'I(metric:probe)':>17s} {'I(metric:probe)':>17s}")
print(f"  {'':>6s} {'':>16s} {'mass OCCUPATION n':>17s} {'signed mass sz':>17s}")
print("  " + "-" * 60)
Hn = HT(1.0, 1.5, coup="n"); Hs_ = HT(1.0, 1.5, coup="sz")
for T in [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 20.0]:
    pn = expm(-1j * Hn * T) @ psi0T
    ps = expm(-1j * Hs_ * T) @ psi0T
    print(f"  {T:6.2f} {info_probe(pn, RT2):16.6f} {info_probe(pn, MT2):17.6f} "
          f"{info_probe(ps, MT2):17.3e}")

print("\n  WHY THE SIGNED-MASS COLUMN IS EXACTLY ZERO. Sigma = (-1)^(probe position) is unitary and")
print("  gives Sigma H_hop Sigma = -H_hop while leaving the potential alone, and H_hop and the")
print("  potential are both REAL in this basis. So H(+V) and H(-V) are related by an antiunitary")
print("  that fixes every probe position, and a probe launched from ONE site cannot see the sign of")
print("  a potential -- only its magnitude. A potential is even under path reversal; the reversal")
print("  flips the sign of the mass. Checks:")
Sig = np.diag([(-1.0) ** p for p, _ in STATES]).astype(complex)
Hh1 = HhopT(1.0)
print(f"    || Sigma H_hop Sigma + H_hop ||     = {np.linalg.norm(Sig@Hh1@Sig+Hh1):.2e}")
print(f"    max |Im H_hop| , max |Im POT|       = {np.abs(np.imag(Hh1)).max():.2e} , "
      f"{np.abs(np.imag(POT)).max():.2e}")
print("    the occupation coupling has no such symmetry: n = (1-sz)/2 is not traceless, so the two")
print("    branches are 'potential present' and 'potential absent', not '+V' and '-V'.")

print("\n  CONTROL A -- V = 0. The metric reading must be exactly 0; the gauge reading must survive.")
H_V0 = HT(1.0, 0.0)
for T in [4.0, 12.0]:
    psi = expm(-1j * H_V0 * T) @ psi0T
    print(f"    T={T:5.1f}   I(gauge:probe) = {info_probe(psi,RT2):.6f}   "
          f"I(metric:probe) = {info_probe(psi,MT2):.3e}")

print("\n  CONTROL B -- tau = 0, no transport at all. Both must be 0.")
H_t0 = HT(0.0, 1.5)
for T in [4.0, 12.0]:
    psi = expm(-1j * H_t0 * T) @ psi0T
    print(f"    T={T:5.1f}   I(gauge:probe) = {info_probe(psi,RT2):.3e}   "
          f"I(metric:probe) = {info_probe(psi,MT2):.3e}")

print("\n  CONTROL C -- THE DECISIVE ONE. Cut the ring: remove ONE hop, keeping every local")
print("  interaction and the potential, removing only the closed path.")
print("  A gauge record needs the topology. A potential does not. This separates them or it does not.")
H_cut = HT(1.0, 1.5, skip=4)
for T in [4.0, 8.0, 12.0, 20.0]:
    psi = expm(-1j * H_cut * T) @ psi0T
    print(f"    T={T:5.1f}   I(gauge:probe) = {info_probe(psi,RT2):.3e}   "
          f"I(metric:probe) = {info_probe(psi,MT2):.6f}")

print()
print("=" * 100)
print("  END OF NUMBERS.")
print("=" * 100)
