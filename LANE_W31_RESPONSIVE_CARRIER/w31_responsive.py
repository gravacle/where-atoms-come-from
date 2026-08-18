"""W-31.  THE EVAPORATING CAP -- a responsive carrier, run against its own worst controls.

CONSTRUCTION (one line).  Make the carrier's FIRST BETTI NUMBER dynamical.  A Z_2 gauge field
lives on an annulus of m quads whose central disk-CAP (hub vertex + m spokes) can evaporate.
n=1: the complex is a DISK, b1=0.  n=0: it is an ANNULUS, b1=1.  The record is the winding
parity of the electric flux around the hole, R_E = (-1)^(s_{a0}+s_{c0}).  While the cap exists
the cap triangles fail to commute with R_E, so R_E is WRITABLE.  When the cap goes, the only
operator that failed to commute with R_E is physically removed, so R_E is EXACTLY conserved.

FRAMING CORRECTION (adopted from the losing proposals; it changes what is being claimed).
  * The W-30 "no-go" is not a theorem that anything here breaks.  R unitary and [L,R]=0 gives
    [L^dag,R]=0 and hence Tr(R D[rho]) = 0 for ANY L, unitary or not -- a tautology about
    superselection: a conserved charge is conserved.  STEP 4b verifies the NON-unitary case
    explicitly, so the lane stops treating unitarity as the load-bearing hypothesis.
  * W-30's R = prod_p W_p is, by discrete Stokes on a DISK, a function of the local plaquettes:
    the dual-Ising spin-flip GENERATOR, not an order parameter.  Generators are conserved,
    hence unwritable.  That is exactly and only what W-30 measured.
  * This construction does not "break" that.  It makes R_E non-conserved for a finite window
    and then deletes the term that made it so.  R_E is a HOLONOMY of a nontrivial cycle, which
    exists only once the hole is open.  That is what opening the hole buys.

WHAT IS AT RISK OF BEING CONSTRAINT-FORCED (declared before any dynamics; settled in STEP 2).
  p_+(infinity), the frozen weight of the R_E=+1 sector.  Three ways it could be forced:
   (i)   unequal sector dimensions -- counting alone would fix p_+;
   (ii)  an evaporation GATE preferring a sector -- the disqualified (1+R)/2 projector smuggled
         into the geometry;
   (iii) at large g^2 the state is pinned near the flux vacuum, where the winding is trivially 0.
  (i) and (ii) are settled by traces taken before any time step: tr(R_E P_punctured) and
  tr(R_E * sum_k L_k^dag L_k) for EVERY junction family.  (iii) is settled by no commutator and
  is left to the g^2 scan in STEP 9, where it is expected to show up as a collapse of the spread.

THE LOAD-BEARING CHOICE, RUN AS A CONTROL RATHER THAN ADMIRED (STEP 7).
  Topology change needs a junction condition.  The two natural clearings of the vanishing spokes
  differ by the kernel {0, sum_i tau_i}, and sum_i tau_i IS the inner cycle, whose winding is 1 --
  so the ambiguity in the junction condition is EXACTLY the record.  Four junction families are
  therefore implemented and all four are run:
    A  NO-HAIR    : one jump, sqrt(Gd) iota^dag.  The cap may vanish only where the field left no
                    flux on it.  Field-GATED.  sum L^dag L = Gd * Pi_nohair.
    B  DRAG/BLIND : the full Kraus family {sqrt(Gd) V_p}, sum_p V_p^dag V_p = I_capped.  The cap
                    carries its flux away, at the SAME rate from every state -- evaporation BLIND
                    to the field.  Branch convention T_0 = 0.
    Bp FLIPPED    : the other branch, T_0 = 1.  Every V_p carries an extra M_inner.
    C  SYMMETRISED: half rate on each branch.  If the correct treatment symmetrises, the design's
                    own prediction is that R_E is randomised at every topology change.
  Nothing in the construction derives the choice.  The numbers say what each choice gives.

ROUTE.  m=2 (dim 24, superoperator 576x576): EXACT Lindblad superoperator exponentiation, expm by
scaling-and-squaring Taylor (scipy is absent).  m=3 (dim 80) and m=4 (dim 288): RK4 on rho, with a
convergence check against the m=2 exact result.  numpy only.

GAUGE STRUCTURE, HONESTLY.  Within each fiber the full Z_2 Gauss law holds exactly (the fiber IS
the divergence-free flux sector of that graph).  ACROSS fibers there is no single gauge group --
the hub's Gauss law exists only when the hub exists -- so this is a bundle of constraint surfaces,
one per carrier configuration, with no canonical isomorphism between them.  A SECOND CONSTRAINT is
added for junction A and is written down here: the NO-HAIR JUNCTION CONDITION, E_sigma_i|psi> = 0
on the vanishing links, implemented by iota with iota iota^dag = Pi_nohair.  Junctions B/Bp/C drop
it, which is exactly why they are run.
"""

import itertools, numpy as np

# ----------------------------------------------------------------------------------
# numerics: expm by scaling and squaring Taylor (no scipy anywhere in this file)
# ----------------------------------------------------------------------------------
def expm(A):
    nrm = np.linalg.norm(A, np.inf)
    n = max(0, int(np.ceil(np.log2(nrm))) + 1) if nrm > 0 else 0
    B = A / (2.0 ** n)
    X = np.eye(A.shape[0], dtype=complex)
    T = X.copy()
    for k in range(1, 100):
        T = T @ B / k
        X = X + T
        if np.linalg.norm(T, np.inf) < 1e-18 * max(1.0, np.linalg.norm(X, np.inf)):
            break
    for _ in range(n):
        X = X @ X
    return X


# ----------------------------------------------------------------------------------
# the carrier
# ----------------------------------------------------------------------------------
def graph(m, capped):
    """vertices: inner 0..m-1, outer m..2m-1, hub 2m.
       links   : a_i = i (inner cycle), c_i = m+i (outer cycle), r_i = 2m+i (rungs),
                 sigma_i = 3m+i (spokes, only when capped)."""
    E = [(i, (i + 1) % m) for i in range(m)]
    E += [(m + i, m + ((i + 1) % m)) for i in range(m)]
    E += [(i, m + i) for i in range(m)]
    if capped:
        E += [(2 * m, i) for i in range(m)]
    return E, 2 * m + (1 if capped else 0)


def cycle_space(E, NV):
    """Divergence-free Z_2 flux configs: the Gauss-law physical sector, as in w27_patch.py."""
    out = []
    for s in itertools.product((0, 1), repeat=len(E)):
        d = [0] * NV
        for k, (u, v) in enumerate(E):
            if s[k]:
                d[u] ^= 1
                d[v] ^= 1
        if not any(d):
            out.append(s)
    return out


def move(S, idx, cyc):
    """M_p |s> = |s + p>.  W-27 / W-30 'Move' convention, flux basis."""
    n = len(S)
    M = np.zeros((n, n), complex)
    for j, s in enumerate(S):
        t = list(s)
        for k in cyc:
            t[k] ^= 1
        M[idx[tuple(t)], j] = 1.0
    return M


def zdiag(S, k):
    return np.array([(-1.0) ** s[k] for s in S], dtype=complex)


class Model:
    def __init__(self, m):
        self.m = m
        E1, NV1 = graph(m, True)
        E0, NV0 = graph(m, False)
        self.S1 = cycle_space(E1, NV1)
        self.S0 = cycle_space(E0, NV0)
        self.D1, self.D0 = len(self.S1), len(self.S0)
        self.D = self.D1 + self.D0
        self.i1 = {s: j for j, s in enumerate(self.S1)}
        self.i0 = {s: j for j, s in enumerate(self.S0)}
        a = lambda i: i
        c = lambda i: m + i
        r = lambda i: 2 * m + i
        sg = lambda i: 3 * m + i
        self.a, self.c, self.r, self.sg = a, c, r, sg
        self.quad = lambda i: [a(i), r(i), c(i), r((i + 1) % m)]
        self.tri = lambda i: [a(i), sg(i), sg((i + 1) % m)]

        D, D1, D0 = self.D, self.D1, self.D0
        self.Mq1 = [move(self.S1, self.i1, self.quad(i)) for i in range(m)]
        self.Mt1 = [move(self.S1, self.i1, self.tri(i)) for i in range(m)]
        self.Mq0 = [move(self.S0, self.i0, self.quad(i)) for i in range(m)]
        self.Minner1 = move(self.S1, self.i1, [a(i) for i in range(m)])

        # record: winding parity across a radial cut gamma_j = {a_j, c_j}
        self.RE1_all = [np.array([(-1.0) ** (s[a(j)] + s[c(j)]) for s in self.S1]) for j in range(m)]
        self.RE0_all = [np.array([(-1.0) ** (s[a(j)] + s[c(j)]) for s in self.S0]) for j in range(m)]
        self.rd = np.concatenate([self.RE1_all[0], self.RE0_all[0]]).astype(complex)
        self.RE = np.diag(self.rd)
        I = np.eye(D, dtype=complex)
        self.Pp = (I + self.RE) / 2.0
        self.Pm = (I - self.RE) / 2.0
        self.PCAP = np.zeros((D, D), complex)
        self.PCAP[:D1, :D1] = np.eye(D1)

        # electric operators on the total space (spokes act as identity on the punctured fiber)
        self.zt = []
        for k in range(3 * m):
            self.zt.append(np.concatenate([zdiag(self.S1, k), zdiag(self.S0, k)]))
        for k in range(3 * m, 4 * m):
            self.zt.append(np.concatenate([zdiag(self.S1, k), np.ones(D0, dtype=complex)]))

        # iota : punctured -> capped, identity on flux, zero on the spokes
        iota = np.zeros((D1, D0), complex)
        for j, s in enumerate(self.S0):
            iota[self.i1[tuple(list(s) + [0] * m)], j] = 1.0
        self.iota = iota
        self.Pnh = iota @ iota.conj().T
        self.IOTA = np.zeros((D, D), complex)
        self.IOTA[:D1, D1:] = iota

        # ---- junction families -------------------------------------------------
        # tau_i = a_i + sigma_i + sigma_{i+1}, so applying the set T flips sigma_j exactly
        # T_j + T_{j-1} times.  Clearing spoke pattern p therefore needs T_j + T_{j-1} = p_j.
        # Consistency around the loop needs sum_j p_j = 0, which the hub Gauss law guarantees.
        # Two solutions, differing by the all-ones set, i.e. by prod_i M_tau_i = M_inner, whose
        # winding is 1 -- the junction ambiguity IS the record.  T_0 = 0 fixes the branch.
        spoke_pat = [tuple(s[sg(i)] for i in range(m)) for s in self.S1]
        self.pats = sorted(set(spoke_pat))
        Vs = []
        for p in self.pats:
            T = [0] * m
            for j in range(1, m):
                T[j] = T[j - 1] ^ p[j]
            U = np.eye(D1, dtype=complex)
            for i in range(m):
                if T[i]:
                    U = U @ self.Mt1[i]
            Pi = np.diag([1.0 if spoke_pat[j] == p else 0.0 for j in range(D1)]).astype(complex)
            Vs.append((iota.conj().T @ U @ Pi, Pi))          # D0 x D1, isometric on the p sector
        self.piso = max(np.linalg.norm(V.conj().T @ V - Pi) for V, Pi in Vs)
        self.Vs = [V for V, _ in Vs]
        self.Vs_flip = [V @ self.Minner1 for V in self.Vs]
        self.kraus_defect = np.linalg.norm(sum(V.conj().T @ V for V in self.Vs) - np.eye(D1))

    def emb_down(self, V):
        M = np.zeros((self.D, self.D), complex)
        M[self.D1:, :self.D1] = V
        return M

    def junction(self, name, Gd):
        if Gd == 0.0:
            return []
        if name == "A":
            return [np.sqrt(Gd) * self.emb_down(self.iota.conj().T)]
        if name == "B":
            return [np.sqrt(Gd) * self.emb_down(V) for V in self.Vs]
        if name == "Bp":
            return [np.sqrt(Gd) * self.emb_down(V) for V in self.Vs_flip]
        if name == "C":
            return ([np.sqrt(Gd / 2) * self.emb_down(V) for V in self.Vs] +
                    [np.sqrt(Gd / 2) * self.emb_down(V) for V in self.Vs_flip])
        raise ValueError(name)

    def elec(self, gE):
        return [np.sqrt(gE) * np.diag(d) for d in self.zt]

    def H(self, g2, lam=0.0, kap=1.0, mu=0.0):
        m, D1 = self.m, self.D1
        H1 = -sum(M + M.conj().T for M in self.Mq1) - kap * sum(M + M.conj().T for M in self.Mt1)
        H1 = H1 - g2 * sum(2 * np.diag(zdiag(self.S1, k)) for k in range(4 * m)) + mu * np.eye(D1)
        H0 = -sum(M + M.conj().T for M in self.Mq0)
        H0 = H0 - g2 * sum(2 * np.diag(zdiag(self.S0, k)) for k in range(3 * m))
        H = np.zeros((self.D, self.D), complex)
        H[:D1, :D1] = H1
        H[D1:, D1:] = H0
        return H - lam * (self.IOTA + self.IOTA.conj().T)

    def Hfib(self, g2, kap=1.0, mu=0.0):
        H = self.H(g2, 0.0, kap, mu)
        return H[:self.D1, :self.D1], H[self.D1:, self.D1:]

    def starts(self):
        """Initial states.  UNBIASED family: (|s> + e^{i th}|s + inner>)/sqrt2.  Adding the inner
        cycle flips a_0, hence flips R_E, so every member has p_+(0) = 1/2 EXACTLY regardless of
        s and th.  Different s and th differ only in quantum numbers R_E cannot see -- which is
        the whole point of the spread test.  BIASED members C,E are labelled as such."""
        m, D = self.m, self.D

        def kk(s):
            v = np.zeros(D, complex)
            v[self.i1[s]] = 1.0
            return v

        def pair(base, th):
            u = list(base)
            for k in range(m):
                u[k] ^= 1
            return (kk(tuple(base)) + np.exp(1j * th) * kk(tuple(u))) / np.sqrt(2)

        z = [0] * (4 * m)
        q0 = list(z)
        for k in self.quad(0):
            q0[k] ^= 1
        t0 = list(z)
        for k in self.tri(0):
            t0[k] ^= 1
        q1 = list(z)
        for k in self.quad(1 % m):
            q1[k] ^= 1
        out = {"A": pair(z, 0.0), "B": pair(q0, 0.0), "D": pair(z, np.pi / 2),
               "F": pair(t0, 0.0), "G": pair(q1, np.pi / 2)}
        rg = np.random.default_rng(4242 + m)
        vC = np.zeros(D, complex)
        vC[:self.D1] = rg.normal(size=self.D1) + 1j * rg.normal(size=self.D1)
        out["C"] = vC / np.linalg.norm(vC)
        vE = np.zeros(D, complex)
        vE[:self.D1] = rg.normal(size=self.D1) + 1j * rg.normal(size=self.D1)
        vE = self.Pp @ vE + 0.4 * (self.Pm @ vE)
        out["E"] = vE / np.linalg.norm(vE)
        return out


UNB = ("A", "B", "D", "F", "G")     # p_+(0) = 1/2 exactly
BIA = ("C", "E")                    # deliberately biased starts


# ----------------------------------------------------------------------------------
# Lindblad machinery
# ----------------------------------------------------------------------------------
def superop(H, Ls, D):
    """Vectorisation is ROW-major (numpy C order), to match rho.reshape(-1) below:
       vec(A X B) = (A kron B^T) vec(X).  NOTE: w30b used the column-major kron with a
       row-major reshape.  That mismatch is invisible whenever rho_0 is real symmetric --
       which every state W-30b evolved happened to be -- and wrong otherwise.  The complex
       initial states used here (relative phase pi/2 between winding sectors) expose it, so
       the convention is fixed here and cross-checked against a naive matrix RK4 in STEP 12."""
    Id = np.eye(D, dtype=complex)
    M = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for L in Ls:
        LdL = L.conj().T @ L
        M += np.kron(L, L.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T))
    return M


_UC = {}


def propagator(H, Ls, dt, D, key):
    if key not in _UC:
        _UC[key] = expm(superop(H, Ls, D) * dt)
    return _UC[key]


def evolve_exact(H, Ls, rho0, ts, D, key):
    U = propagator(H, Ls, ts[1] - ts[0], D, key)
    v = rho0.reshape(-1).astype(complex)
    out = []
    for n in range(len(ts)):
        if n > 0:
            v = U @ v
        r = v.reshape(D, D)
        out.append(r / np.trace(r).real)
    return out


def obs(md, r):
    return (np.trace(md.Pp @ r).real,
            np.linalg.norm(md.Pp @ r @ md.Pm),
            np.trace(md.PCAP @ r).real,
            np.trace(md.RE @ r).real)


def rho_of(v):
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def late(o, lo=30, hi=41, i=0):
    return float(np.mean([o[k][i] for k in range(lo, hi)]))


# ==================================================================================
print("=" * 104)
print("W-31  THE EVAPORATING CAP -- responsive carrier.  numpy only, no scipy.")
print("=" * 104)

# ==================================================================================
# STEP 1  --  THE SPACE, AND THE STATIC COMPARISON CASE
# ==================================================================================
print("\n[STEP 1] THE SPACE\n")
MD = {}
for m in (2, 3, 4):
    MD[m] = Model(m)
    x = MD[m]
    print(f"  m={m}   capped fiber (DISK, b1=0) dim {x.D1:4d} (pred 2^{2*m} = {2**(2*m):4d})"
          f"   punctured fiber (ANNULUS, b1=1) dim {x.D0:4d} (pred 2^{m+1} = {2**(m+1):4d})"
          f"   TOTAL {x.D:4d}   superoperator {x.D**2}x{x.D**2}")
    print(f"         junction Kraus completeness ||sum_p V_p^dag V_p - I|| = {x.kraus_defect:.3e}"
          f"   max_p ||V_p^dag V_p - Pi_p|| = {x.piso:.3e}"
          f"   #spoke patterns = {len(x.pats)} (pred 2^{m-1} = {2**(m-1)})"
          f"   ||iota^dag iota - I|| = "
          f"{np.linalg.norm(x.iota.conj().T @ x.iota - np.eye(x.D0)):.3e}")
md = MD[2]
D = md.D

# static comparison case: W-30's own carrier, rebuilt from scratch here
V2 = [(i, j) for j in range(3) for i in range(3)]
vid = {v: k for k, v in enumerate(V2)}
Ep = []
for j in range(3):
    for i in range(2):
        Ep.append((vid[(i, j)], vid[(i + 1, j)]))
for j in range(2):
    for i in range(3):
        Ep.append((vid[(i, j)], vid[(i, j + 1)]))
Sp = cycle_space(Ep, 9)
ip = {s: j for j, s in enumerate(Sp)}
Dp = len(Sp)
hid = lambda i, j: j * 2 + i
vx = lambda i, j: 6 + j * 3 + i
PLQ = [[hid(i, j), vx(i + 1, j), hid(i, j + 1), vx(i, j)] for j in range(2) for i in range(2)]
print(f"\n  STATIC COMPARISON CASE (W-30's carrier, rebuilt here): 3x3 planar patch Z_2,"
      f" physical dim {Dp}, superoperator {Dp**2}x{Dp**2}")
print(f"  ROUTE: m=2 EXACT superoperator exponential; m=3,4 RK4 on rho, convergence-checked"
      f" against the m=2 exact result.")

# ==================================================================================
# STEP 2  --  FORCED-OR-NOT, RUN FIRST
# ==================================================================================
print("\n[STEP 2] FORCED-OR-NOT CHECK  (run before any dynamics)\n")
print("  QUANTITY AT RISK: p_+(inf), the frozen weight of the R_E = +1 winding sector.")
print("  Forced if (i) the two winding sectors have unequal dimension, or (ii) any evaporation")
print("  gate sum_k L_k^dag L_k has nonzero overlap with R_E.")
print("  COUNTING ARGUMENT, declared in advance: w(s) = s_{a0}+s_{c0} is a LINEAR functional on")
print("  the Z_2 cycle space and is surjective (w(inner cycle) = 1), so its kernel has index")
print("  exactly 2 and the sectors must be equal.\n")
forced = False
for m in (2, 3, 4):
    x = MD[m]
    n0p = int(np.sum(np.real(x.RE0_all[0]) > 0))
    n1p = int(np.sum(np.real(x.RE1_all[0]) > 0))
    tr0 = float(np.sum(x.RE0_all[0]))
    tr1 = float(np.sum(x.RE1_all[0]))
    tnh = float(np.real(np.trace(np.diag(x.RE1_all[0]).astype(complex) @ x.Pnh)))
    print(f"  m={m}  punctured sectors {n0p}/{x.D0-n0p}  tr(R_E P_punc) = {tr0:+.1f}   "
          f"capped sectors {n1p}/{x.D1-n1p}  tr(R_E P_cap) = {tr1:+.1f}   "
          f"tr(R_E Pi_nohair) = {tnh:+.1f}")
    RT = np.diag(np.concatenate([x.RE1_all[0], x.RE0_all[0]])).astype(complex)
    for nm in ("A", "B", "Bp", "C"):
        G = sum(L.conj().T @ L for L in x.junction(nm, 1.0))
        ov = float(np.real(np.trace(RT @ G)))
        print(f"        junction {nm:2s}: tr(R_E * sum_k L^dag L) = {ov:+.3e}"
              f"   tr(sum_k L^dag L) = {float(np.real(np.trace(G))):.1f}")
        if abs(ov) > 1e-9:
            forced = True
    if abs(tr0) > 1e-9 or abs(tnh) > 1e-9:
        forced = True
print(f"\n  Any nonzero overlap found (record constraint-forced / gate sector-biased)?  {forced}")
if forced:
    print("  THE RECORD OBSERVABLE IS CONSTRAINT-FORCED OR THE GATE IS BIASED.  STOPPING HERE.")
    raise SystemExit(0)
print("  Residual risk (iii), large-g^2 pinning, is settled by NO commutator.  See STEP 9.")

# ==================================================================================
# STEP 3  --  DOES R EARN THE RIGHT TO BE MEASURED?
# ==================================================================================
print("\n[STEP 3] THE RECORD OPERATOR MUST EARN ITS MEASUREMENT\n")
for m in (2, 3, 4):
    x = MD[m]
    R = x.RE
    ev = np.linalg.eigvals(R)
    print(f"  m={m}  ||R_E|| = {np.linalg.norm(R):8.4f}   ||R R^dag - I|| = "
          f"{np.linalg.norm(R @ R.conj().T - np.eye(x.D)):.3e}   distinct eigenvalues = "
          f"{len(np.unique(np.round(ev, 9)))}   cut-dependence max_j||R(gamma_j)-R(gamma_0)||:"
          f"  punctured {max(np.linalg.norm(x.RE0_all[j]-x.RE0_all[0]) for j in range(m)):.3e}"
          f"   capped {max(np.linalg.norm(x.RE1_all[j]-x.RE1_all[0]) for j in range(m)):.3e}")

# ==================================================================================
# STEP 4  --  REPRODUCE THE STATIC NO-GO INSIDE THIS SCRIPT
# ==================================================================================
print("\n[STEP 4] STATIC NO-GO CONTROLS  (if these do not reproduce, nothing below means anything)")

def mv_p(cyc):
    M = np.zeros((Dp, Dp), complex)
    for j, s in enumerate(Sp):
        t = list(s)
        for k in cyc:
            t[k] ^= 1
        M[ip[tuple(t)], j] = 1.0
    return M

Pl = [mv_p(q) for q in PLQ]
acc = {}
for q in PLQ:
    for k in q:
        acc[k] = acc.get(k, 0) + 1
Rrim = mv_p([k for k, v in acc.items() if v % 2])
rng = np.random.default_rng(20260817)

def rand_rho(n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    r = A @ A.conj().T
    return r / np.trace(r).real

def parts(R, H, Ls, n, ntr=20):
    """returns (max |Tr(R * dissipator)| , max |Tr(R * full drho/dt)|) over random rho"""
    wd = wf = 0.0
    for _ in range(ntr):
        r = rand_rho(n)
        dd = sum(L @ r @ L.conj().T - 0.5 * (L.conj().T @ L @ r + r @ L.conj().T @ L) for L in Ls)
        df = -1j * (H @ r - r @ H) + dd
        wd = max(wd, abs(np.trace(R @ dd)))
        wf = max(wf, abs(np.trace(R @ df)))
    return wd, wf

print("\n  4a  W-30 REPLICA on the 3x3 patch, R = product of the 4 plaquettes (discrete Stokes),")
print("      jumps = the 4 plaquettes (unitary).  W-30 measured: at g^2 = 0 the record is frozen.")
print(f"      ||R|| = {np.linalg.norm(Rrim):.4f}   unitarity defect = "
      f"{np.linalg.norm(Rrim @ Rrim.conj().T - np.eye(Dp)):.3e}   distinct eigenvalues = "
      f"{len(np.unique(np.round(np.linalg.eigvals(Rrim), 9)))}")
for g2 in (0.0, 1.0):
    Hq = -sum(L + L.conj().T for L in Pl) - g2 * sum(2 * np.diag(zdiag(Sp, k)) for k in range(len(Ep)))
    wd, wf = parts(Rrim, Hq, Pl, Dp)
    print(f"      g^2={g2:4.1f}  ||[H,R]|| = {np.linalg.norm(Hq @ Rrim - Rrim @ Hq):9.3e}   "
          f"max||[L,R]|| = {max(np.linalg.norm(L @ Rrim - Rrim @ L) for L in Pl):.3e}   "
          f"max|Tr(R D[rho])| = {wd:.3e}   max|d<R>/dt| = {wf:.3e}")

print("\n  4b  TAUTOLOGY GRAFT: the no-go never needed unitary jumps.  [L,R]=0 with R unitary")
print("      gives [L^dag,R]=0 and Tr(R D[rho]) = 0 term by term, for ANY L.  Measured:")
Bc = rng.normal(size=(Dp, Dp)) + 1j * rng.normal(size=(Dp, Dp))
Ppp = (np.eye(Dp) + Rrim) / 2.0
Pmm = (np.eye(Dp) - Rrim) / 2.0
Hq0 = -sum(L + L.conj().T for L in Pl)
for lab, L in (("unitary plaquette", Pl[0]),
               ("NON-unitary, [L,R]=0", Ppp @ Bc @ Ppp + Pmm @ Bc @ Pmm),
               ("NON-unitary, [L,R]!=0", Ppp @ Bc @ Pmm)):
    wd, _ = parts(Rrim, Hq0, [L], Dp)
    print(f"      {lab:24s}  ||L L^dag - I|| = "
          f"{np.linalg.norm(L @ L.conj().T - np.eye(Dp)):9.3e}   ||[L,R]|| = "
          f"{np.linalg.norm(L @ Rrim - Rrim @ L):9.3e}   max|Tr(R D[rho])| = {wd:.3e}")

print("\n  4c  ELITZUR / STOKES CONTRAST inside this construction: is the record a GENERATOR of")
print("      the local plaquette algebra (then conserved, hence unwritable) or a HOLONOMY?")
for m in (2, 3, 4):
    x = MD[m]
    pr = np.eye(x.D1, dtype=complex)
    for M in x.Mq1 + x.Mt1:
        pr = pr @ M
    pr0 = np.eye(x.D0, dtype=complex)
    for M in x.Mq0:
        pr0 = pr0 @ M
    out1 = move(x.S1, x.i1, [x.c(i) for i in range(m)])
    out0 = move(x.S0, x.i0, [x.c(i) for i in range(m)])
    print(f"      m={m}  ||prod(all plaquettes) - M_outer||    capped (DISK) = "
          f"{np.linalg.norm(pr - out1):.3e}    punctured (ANNULUS) = {np.linalg.norm(pr0 - out0):.3e}")

print("\n  4d  COMMUTATORS OF R_E WITH THE FULL GENERATOR (g^2 = 0.7, kappa = 1):")
for m in (2, 3, 4):
    x = MD[m]
    H1, H0 = x.Hfib(0.7)
    R1 = np.diag(x.RE1_all[0]).astype(complex)
    R0 = np.diag(x.RE0_all[0]).astype(complex)
    cA = max(np.linalg.norm(x.RE @ L - L @ x.RE)
             for L in x.elec(0.5) + x.junction("A", 0.8) + x.junction("B", 0.8))
    cC = max(np.linalg.norm(x.RE @ L - L @ x.RE) for L in x.junction("C", 0.8))
    nw = sum(1 for M in x.Mt1 if np.linalg.norm(R1 @ M - M @ R1) > 1e-9)
    nq = sum(1 for M in x.Mq1 if np.linalg.norm(R1 @ M - M @ R1) > 1e-9)
    print(f"      m={m}  ||[R_E,H_capped]|| = {np.linalg.norm(R1 @ H1 - H1 @ R1):9.4f}   "
          f"||[R_E,H_punctured]|| = {np.linalg.norm(R0 @ H0 - H0 @ R0):.3e}   "
          f"max||[R_E,L]|| junctions A,B + electric = {cA:.3e}   junction C = {cC:.3e}")
    print(f"            HOW MANY TERMS ACTUALLY WRITE: cap triangles with [R_E,M_tau] != 0 : "
          f"{nw} of {m}    quads with [R_E,M_q] != 0 : {nq} of {m}    "
          f"total plaquette terms in H_capped: {2*m}")

print("\n  4e  CARRIER RESPONSE OFF (Gamma_d = 0, lambda = 0).")
print("      START PUNCTURED -- the carrier is static and R_E is then a conserved charge, so")
print("      d<R_E>/dt must be 0 to machine precision:")
ts = np.linspace(0, 40, 41)
vpun = np.zeros(D, complex)
vpun[md.D1:] = rng.normal(size=md.D0) + 1j * rng.normal(size=md.D0)
r_pun = rho_of(vpun)
for g2 in (0.0, 0.7, 3.0):
    o = [obs(md, r) for r in evolve_exact(md.H(g2), md.elec(0.5), r_pun, ts, D, ("off", g2))]
    print(f"        g^2={g2:4.1f}  p+(0) = {o[0][0]:.12f}   p+(40) = {o[-1][0]:.12f}   |delta| = "
          f"{abs(o[-1][0]-o[0][0]):.3e}   <R_E>(40) = {o[-1][3]:+.12f}")
ST = md.starts()
H07 = md.H(0.7)
print("      START CAPPED, carrier still frozen -- the writer never stops, so nothing freezes:")
for nm in UNB + BIA:
    o = [obs(md, r) for r in evolve_exact(H07, md.elec(0.5), rho_of(ST[nm]), ts, D, ("offc",))]
    print(f"        init {nm}: p+(0) = {o[0][0]:.6f}   p+(40) = {o[-1][0]:.6f}   "
          f"late-avg p+[30,40] = {late(o):.6f}   drift[30,40] = {abs(o[-1][0]-o[-11][0]):.3e}   "
          f"p_cap(40) = {o[-1][2]:.3e}")

# ==================================================================================
# STEP 5  --  CARRIER RESPONSE ON
# ==================================================================================
print("\n[STEP 5] CARRIER RESPONSE ON  (junction A, no-hair gate; Gamma_d=0.8, gamma_E=0.5,"
      " g^2=0.7)\n")
print("  DECOY WARNING, adopted in advance from Proposal 3: for a pure state")
print("  ||P+ psi psi^dag P-||_1 = (1/2) sqrt(1 - <R>^2), so the coherence column CANNOT by")
print("  itself witness selection -- phase randomisation drives it to zero with no weight")
print("  transfer at all.  It is printed as a CONSISTENCY CHECK only.  The evidence is p_+ and")
print("  its dependence on the initial state.\n")
Ls_on = md.elec(0.5) + md.junction("A", 0.8)
print(f"  {'init':>5s} {'p+(0)':>9s} {'<R>(0)':>9s} {'p+(1)':>9s} {'p+(2)':>9s} {'p+(5)':>9s} "
      f"{'p+(10)':>9s} {'p+(40)':>9s} {'lateavg':>9s} {'<R>late':>9s} {'drift':>9s} "
      f"{'coh(40)':>9s} {'pcap(40)':>9s}")
print("  " + "-" * 133)
p_on = {}
for nm in UNB + BIA:
    o = [obs(md, r) for r in evolve_exact(H07, Ls_on, rho_of(ST[nm]), ts, D, ("on", 0.7, "A", 0.8))]
    p_on[nm] = late(o)
    print(f"  {nm:>5s} {o[0][0]:9.6f} {o[0][3]:+9.6f} {o[1][0]:9.6f} {o[2][0]:9.6f} "
          f"{o[5][0]:9.6f} {o[10][0]:9.6f} {o[40][0]:9.6f} {late(o):9.6f} "
          f"{late(o,i=3):+9.6f} {abs(o[-1][0]-o[-11][0]):9.3e} {o[-1][1]:9.2e} {o[-1][2]:9.2e}")
o = [obs(md, r) for r in evolve_exact(H07, Ls_on, rho_of(ST["A"]), ts, D, ("on", 0.7, "A", 0.8))]
print("\n  <R_E>(t), init A:  " + "   ".join(f"t={ts[k]:.0f}: {o[k][3]:+.6f}"
                                             for k in (0, 1, 2, 3, 5, 10, 20, 40)))
sp_unb = max(p_on[k] for k in UNB) - min(p_on[k] for k in UNB)
sp_all = max(p_on.values()) - min(p_on.values())
print(f"  SPREAD over the five UNBIASED starts (p+(0) = 1/2 exactly): max-min = {sp_unb:.6f}")
print(f"  max-min including the two BIASED starts C,E (p+(0) = "
      f"{obs(md, rho_of(ST['C']))[0]:.4f}, {obs(md, rho_of(ST['E']))[0]:.4f}): {sp_all:.6f}")

# ==================================================================================
# STEP 6  --  LINDBLADIAN SPECTRUM
# ==================================================================================
print("\n[STEP 6] LINDBLADIAN STRUCTURE  (adopted from Proposal 2: report structure, not a drift)")
print("  This replaces 'a plateau value and a drift' with an object that needs no protocol, no")
print("  initial state and no time evolution.  In the Heisenberg picture <R_E>(t) = Tr(Y(t)^dag")
print("  rho_0) with Y(t) = exp(M^dag t)[R_E].  Y_inf = lim Y(t) is the ASYMPTOTICALLY CONSERVED")
print("  PART of the record.  Its PUNCTURED block is the already-written record and is trivially")
print("  conserved for every junction (nothing acts on the punctured fiber).  Its CAPPED block is")
print("  the whole question: it is nonzero exactly when a state that starts on the DISK -- no")
print("  hole, no record yet -- ends with a winding that depends on which state it was.\n")
vR = md.RE.reshape(-1)
nR = np.linalg.norm(vR)
for jn in ("A", "B", "Bp", "C"):
    Ls = md.elec(0.5) + md.junction(jn, 0.8)
    M = superop(H07, Ls, D)
    Y = (expm(M.conj().T * 400.0) @ vR).reshape(D, D)
    w = np.linalg.eigvals(M)
    nk = int(np.sum(np.abs(w) < 1e-8))
    print(f"  junction {jn:2s}:  dim ker(M) = {nk:3d}   Hermiticity defect of Y_inf = "
          f"{np.linalg.norm(Y - Y.conj().T):.2e}   tr(Y_inf) = {np.trace(Y).real:+.3e}")
    print(f"      ||Y_inf||/||R_E|| = {np.linalg.norm(Y)/nR:.6f}"
          f"    PUNCTURED block {np.linalg.norm(Y[md.D1:, md.D1:])/nR:.6f}"
          f"    CAPPED block {np.linalg.norm(Y[:md.D1, :md.D1])/nR:.6f}")
    print(f"      predicted <R_E>(inf) from Y_inf, unbiased starts:  " +
          "  ".join(f"{k}:{float(np.real(np.trace(Y.conj().T @ rho_of(ST[k])))):+.6f}" for k in UNB))
    dec = sorted([z.real for z in w if z.real < -1e-9], reverse=True)[:3]
    print(f"      three slowest decay rates Re(lambda): " + ", ".join(f"{r:.4e}" for r in dec))

# ==================================================================================
# STEP 7  --  JUNCTION CONTROLS AND NULLS
# ==================================================================================
print("\n[STEP 7] JUNCTION CONTROLS AND NULL CONTROLS  (g^2=0.7, gamma_E=0.5, Gamma_d=0.8)\n")
print("  DECLARED PREDICTIONS, made before this table was run:")
print("    A  no-hair, field-GATED     -> record forms (spread > 0)")
print("    B  blind/drag, branch T0=0  -> if the spread SURVIVES, the field-gating is DECORATION")
print("    Bp blind/drag, branch T0=1  -> mirrored record (same magnitude, opposite sign)")
print("    C  symmetrised over the two branches -> R_E randomised, spread -> 0\n")
print(f"  {'junction':>9s} {'#jumps':>7s} {'init':>5s} {'p+(0)':>9s} {'late p+':>10s} "
      f"{'<R_E>late':>11s} {'drift':>10s} {'pcap(40)':>9s}")
print("  " + "-" * 78)
jres = {}
for jn in ("A", "B", "Bp", "C"):
    Ls = md.elec(0.5) + md.junction(jn, 0.8)
    for nm in UNB:
        o = [obs(md, r) for r in evolve_exact(H07, Ls, rho_of(ST[nm]), ts, D, ("j", jn))]
        jres[(jn, nm)] = late(o)
        print(f"  {jn:>9s} {len(Ls):7d} {nm:>5s} {o[0][0]:9.6f} {late(o):10.6f} "
              f"{late(o,i=3):+11.6f} {abs(o[-1][0]-o[-11][0]):10.3e} {o[-1][2]:9.2e}")
    sp = max(jres[(jn, k)] for k in UNB) - min(jres[(jn, k)] for k in UNB)
    print(f"  {'':>9s} {'':>7s} SPREAD over the five unbiased starts = {sp:.6f}")

print("\n  NULL CONTROL N1 -- start ALREADY punctured, full bath, long time.  Must be exactly flat.")
for T in (40.0, 200.0):
    tsl = np.linspace(0, T, 5)
    o = [obs(md, r) for r in evolve_exact(md.H(3.0), md.elec(0.5) + md.junction("A", 0.8),
                                          r_pun, tsl, D, ("N1", T))]
    print(f"      t=0  p+ = {o[0][0]:.12f}    t={T:.0f}  p+ = {o[-1][0]:.12f}    |delta| = "
          f"{abs(o[-1][0]-o[0][0]):.3e}")

print("\n  NULL CONTROL N2 -- remove the WRITER (kappa = 0, no cap triangles in H).  Then nothing")
print("      can ever be written and every start must stay at its own initial value.")
Hk0 = md.H(0.7, kap=0.0)
print(f"      ||[R_E,H_capped]|| at kappa=0 = "
      f"{np.linalg.norm(np.diag(md.RE1_all[0]).astype(complex) @ Hk0[:md.D1,:md.D1] - Hk0[:md.D1,:md.D1] @ np.diag(md.RE1_all[0]).astype(complex)):.3e}")
for nm in UNB + BIA:
    o = [obs(md, r) for r in evolve_exact(Hk0, Ls_on, rho_of(ST[nm]), ts, D, ("N2",))]
    print(f"      init {nm}: p+(0) = {o[0][0]:.9f}   p+(40) = {o[-1][0]:.9f}   |delta| = "
          f"{abs(o[-1][0]-o[0][0]):.3e}")

print("\n  NULL CONTROL N3 -- p_cap(t) must never be quoted as evidence.  Printed so it cannot be:")
for g2 in (0.0, 0.7, 3.0):
    a10 = obs(md, evolve_exact(md.H(g2), md.elec(0.5) + md.junction("A", 0.8),
                               rho_of(ST["A"]), ts, D, ("N3a", g2))[10])[2]
    b10 = obs(md, evolve_exact(md.H(g2), md.elec(0.5) + md.junction("B", 0.8),
                               rho_of(ST["A"]), ts, D, ("N3b", g2))[10])[2]
    print(f"      g^2={g2:4.1f}  p_cap(10) with the field-gated junction A = {a10:.6e}   "
          f"with the blind junction B = {b10:.6e}")

# ==================================================================================
# STEP 8  --  SCAN THE RESPONSE STRENGTH, INCLUDING EXACTLY ZERO
# ==================================================================================
print("\n[STEP 8] SCAN OF THE CARRIER RESPONSE STRENGTH Gamma_d  (g^2=0.7, gamma_E=0.5).")
print("  Gamma_d = 0 is the static carrier: the W-30 situation reproduced inside this design.\n")
print(f"  {'Gamma_d':>8s} {'junc':>5s} " + " ".join(f"{'p+('+k+')':>10s}" for k in UNB) +
      f" {'SPREAD':>9s} {'drift(A)':>10s} {'pcap(40)':>10s}")
print("  " + "-" * 104)
for jn in ("A", "B", "C"):
    for Gd in (0.0, 0.05, 0.2, 0.8, 3.2):
        Ls = md.elec(0.5) + md.junction(jn, Gd)
        vals, dr, pc = {}, 0.0, 0.0
        for nm in UNB:
            o = [obs(md, r) for r in evolve_exact(H07, Ls, rho_of(ST[nm]), ts, D, ("S8", jn, Gd))]
            vals[nm] = late(o)
            if nm == "A":
                dr, pc = abs(o[-1][0] - o[-11][0]), o[-1][2]
        print(f"  {Gd:8.2f} {jn:>5s} " + " ".join(f"{vals[k]:10.6f}" for k in UNB) +
              f" {max(vals.values())-min(vals.values()):9.6f} {dr:10.3e} {pc:10.2e}")

# ==================================================================================
# STEP 9  --  THE ALPHA AXIS
# ==================================================================================
print("\n[STEP 9] g^2 SCAN  (junction A, Gamma_d=0.8, gamma_E=0.5).  BOTH metrics are printed.")
print("  The honest metric is the spread over UNBIASED starts.  The max-min INCLUDING the two")
print("  biased starts C,E is also printed because it is the metric a hostile reader would")
print("  reach for, and a biased start that is never written simply retains its bias -- which")
print("  inflates that column at large g^2 exactly where nothing is being written.\n")
print(f"  {'g^2':>7s} {'||[R,Hcap]||':>12s} " + " ".join(f"{'p+('+k+')':>9s}" for k in UNB) +
      f" {'p+(C)':>9s} {'p+(E)':>9s} {'SPREAD(unb)':>12s} {'maxmin(all)':>12s}")
print("  " + "-" * 124)
R1d = np.diag(md.RE1_all[0]).astype(complex)
for g2 in (0.0, 0.1, 0.3, 0.7, 1.0, 1.5, 3.0, 8.0, 20.0):
    Hg = md.H(g2)
    Hc = Hg[:md.D1, :md.D1]
    Ls = md.elec(0.5) + md.junction("A", 0.8)
    vals = {}
    for nm in UNB + BIA:
        o = [obs(md, r) for r in evolve_exact(Hg, Ls, rho_of(ST[nm]), ts, D, ("S9", g2))]
        vals[nm] = late(o)
    su = max(vals[k] for k in UNB) - min(vals[k] for k in UNB)
    sa = max(vals.values()) - min(vals.values())
    print(f"  {g2:7.2f} {np.linalg.norm(R1d @ Hc - Hc @ R1d):12.4f} " +
          " ".join(f"{vals[k]:9.6f}" for k in UNB) +
          f" {vals['C']:9.6f} {vals['E']:9.6f} {su:12.6f} {sa:12.6f}")

# ==================================================================================
# STEP 10 --  UNRAVELLING
# ==================================================================================
print("\n[STEP 10] TRAJECTORY UNRAVELLING  (adopted from Proposal 3).  Ensemble coherence going")
print("  to zero is compatible with nothing having been selected.  The per-branch witness is")
print("  E[<R_E>^2] over trajectories and the fraction of trajectories with |<R_E>| > 0.99.")
print("  For the initial states used here <R_E>(0) = 0, so E[<R_E>^2](0) = 0 exactly.\n")

def unravel(H, Ls, psi0, T, dt, ntraj, seed, RE):
    r = np.random.default_rng(seed)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    U = expm(-1j * Heff * dt)
    nst = int(round(T / dt))
    Lstack = np.array(Ls)
    vals = []
    for _ in range(ntraj):
        v = psi0.copy()
        for _ in range(nst):
            w = U @ v
            nn = np.vdot(w, w).real
            if r.random() < 1.0 - nn:
                cand = np.einsum('kij,j->ki', Lstack, v)
                pk = np.einsum('ki,ki->k', cand.conj(), cand).real
                tot = pk.sum()
                if tot <= 0:
                    v = w / np.sqrt(nn)
                    continue
                k = int(r.choice(len(pk), p=pk / tot))
                v = cand[k] / np.linalg.norm(cand[k])
            else:
                v = w / np.sqrt(nn)
        vals.append(float(np.real(np.vdot(v, RE @ v))))
    return np.array(vals)

NTR = 160
for jn in ("A", "B", "C"):
    Ls = md.elec(0.5) + md.junction(jn, 0.8)
    vv = unravel(H07, Ls, ST["A"], 25.0, 0.0125, NTR, 11, md.RE)
    ex = obs(md, evolve_exact(H07, Ls, rho_of(ST["A"]), ts, D, ("S10", jn))[25])[3]
    print(f"  junction {jn:2s}, unitary (infinite-T) electric bath, {NTR} trajectories to t=25,"
          f" jump unravelling:")
    print(f"      E[<R_E>] = {vv.mean():+.6f} +- {vv.std(ddof=1)/np.sqrt(NTR):.6f}   "
          f"(master-equation <R_E>(25) = {ex:+.6f})")
    print(f"      E[<R_E>^2] = {np.mean(vv**2):.6f}   frac |<R_E>| > 0.99 = "
          f"{np.mean(np.abs(vv) > 0.99):.3f}   frac > 0.5 = {np.mean(np.abs(vv) > 0.5):.3f}")

# ==================================================================================
# STEP 11 --  FINITE-TEMPERATURE DAVIES BATH
# ==================================================================================
print("\n[STEP 11] FINITE-TEMPERATURE DAVIES BATH  (adopted from Proposal 2).  Every unitary jump")
print("  operator in this lane, our own L = sqrt(gamma_E) Z_l included, makes the identity")
print("  stationary: an INFINITE-temperature bath, where nothing anywhere keeps a record.")
print("  A(omega) = sum_{E_b - E_a = omega} |a><a| A |b><b|.  Rate gamma(omega) = g0/(1+exp(-w/T)),")
print("  which satisfies KMS, gamma(-w)/gamma(w) = exp(-w/T), and stays BOUNDED as T grows, so")
print("  the temperatures below are comparable at fixed overall coupling.  T = 1e6 is then the")
print("  infinite-temperature Davies bath, the honest analogue of this lane's unitary jumps.\n")

def davies(Hm, As, gamma0, Tmp, dim, tol=1e-6):
    w, U = np.linalg.eigh(Hm)
    Ls = []
    for A in As:
        Ab = U.conj().T @ A @ U
        freqs = {}
        for aa in range(dim):
            for bb in range(dim):
                om = w[bb] - w[aa]                       # A(omega) lowers energy by omega
                freqs.setdefault(round(om / tol) * tol, []).append((aa, bb))
        for om, pairs in freqs.items():
            M = np.zeros((dim, dim), complex)
            for aa, bb in pairs:
                M[aa, bb] = Ab[aa, bb]
            if np.linalg.norm(M) < 1e-12:
                continue
            g = gamma0 / (1.0 + np.exp(-om / Tmp))
            Ls.append(np.sqrt(g) * (U @ M @ U.conj().T))
    return Ls

H1c, H0c = md.Hfib(0.7)
As0 = [np.diag(zdiag(md.S0, k)).astype(complex) for k in range(3 * md.m)]
R0d = np.diag(md.RE0_all[0]).astype(complex)
print("  VALIDATION on the closed punctured fiber: the Gibbs state must be exactly stationary.")
for Tmp in (0.5, 2.0, 10.0):
    Ld = davies(H0c, As0, 0.5, Tmp, md.D0)
    w0, U0 = np.linalg.eigh(H0c)
    gb = U0 @ np.diag(np.exp(-w0 / Tmp)) @ U0.conj().T
    gb = gb / np.trace(gb).real
    dd = -1j * (H0c @ gb - gb @ H0c) + sum(
        L @ gb @ L.conj().T - 0.5 * (L.conj().T @ L @ gb + gb @ L.conj().T @ L) for L in Ld)
    print(f"      T = {Tmp:6.2f}   #jumps = {len(Ld):4d}   ||d rho_Gibbs / dt|| = "
          f"{np.linalg.norm(dd):.3e}   max||[R_E, A(omega)]|| = "
          f"{max(np.linalg.norm(R0d @ L - L @ R0d) for L in Ld):.3e}")

Ast = [np.diag(d) for d in md.zt]
print("\n  FULL RUN with a Davies electric bath (junction A, Gamma_d=0.8, g^2=0.7).  The claimed")
print("  protection is HOMOLOGICAL, not energetic, so the frozen p_+ should not care about T.")
print(f"  {'T':>9s} {'#jumps':>7s} " + " ".join(f"{'p+('+k+')':>10s}" for k in UNB) +
      f" {'SPREAD':>9s} {'drift':>10s}")
print("  " + "-" * 100)
for Tmp in (0.5, 2.0, 10.0, 1.0e6):
    Ld = davies(H07, Ast, 0.5, Tmp, D)
    Ls = Ld + md.junction("A", 0.8)
    vals, dr = {}, 0.0
    for nm in UNB:
        o = [obs(md, r) for r in evolve_exact(H07, Ls, rho_of(ST[nm]), ts, D, ("S11", Tmp))]
        vals[nm] = late(o)
        if nm == "A":
            dr = abs(o[-1][0] - o[-11][0])
    print(f"  {Tmp:9.1f} {len(Ls):7d} " + " ".join(f"{vals[k]:10.6f}" for k in UNB) +
          f" {max(vals.values())-min(vals.values()):9.6f} {dr:10.3e}")

# ==================================================================================
# STEP 12 --  SIZE SCALING
# ==================================================================================
print("\n[STEP 12] SIZE SCALING m = 2,3,4.  The protection is CLAIMED homological, so it must be")
print("  m-independent.  If the frozen p_+ or the spread drifts systematically with m, the")
print("  topological story is decoration on a finite-size accident.\n")

def rk4(x, H, jumps, rho0, T, nst, report=False):
    """RK4 on rho.  Three exploited structures, all exact, none approximations:
       (a) rho and drho/dt are Hermitian and H is Hermitian, so [H,rho] costs ONE matmul;
       (b) the electric jumps are diagonal in the flux basis -> elementwise;
       (c) every junction jump is block-lower (capped -> punctured) with V of shape D0 x D1,
           and in every junction family sum_k L^dag L is DIAGONAL (Gd*Pi_nohair for A, Gd*I for
           B/Bp/C).  That diagonality is verified numerically here, not assumed."""
    D1 = x.D1
    dt = T / nst
    sd = np.zeros(H.shape[0], complex)
    dj, vj = [], []
    Sfull = np.zeros_like(H)
    for kind, X in jumps:
        if kind == 'd':
            dj.append(X)
            sd += np.abs(X) ** 2
            Sfull += np.diag(np.abs(X) ** 2)
        else:                                   # X is V, shape D0 x D1
            vj.append(X)
            blk = X.conj().T @ X
            Sfull[:D1, :D1] += blk
            sd[:D1] += np.diag(blk)
    offdiag = np.linalg.norm(Sfull - np.diag(np.diag(Sfull)))
    if report:
        print(f"      structural check: ||sum_k L^dag L - diag(sum_k L^dag L)|| = {offdiag:.3e}")
    if offdiag > 1e-10:
        raise RuntimeError("sum_k L^dag L is not diagonal; the fast RK4 path is invalid")
    anti = -0.5 * (sd[:, None] + sd[None, :])

    def rhs(r):
        X1 = H @ r
        out = -1j * (X1 - X1.conj().T)
        for d in dj:
            out += np.outer(d, d.conj()) * r
        for V in vj:
            out[D1:, D1:] += V @ r[:D1, :D1] @ V.conj().T
        return out + anti * r

    r = rho0.copy()
    for _ in range(nst):
        k1 = rhs(r)
        k2 = rhs(r + 0.5 * dt * k1)
        k3 = rhs(r + 0.5 * dt * k2)
        k4 = rhs(r + dt * k3)
        r = r + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return r / np.trace(r).real

def jumps_of(x, gE=0.5, Gd=0.8):
    return ([('d', np.sqrt(gE) * d) for d in x.zt] +
            [('v', np.sqrt(Gd) * x.iota.conj().T)])

print("  CONVERGENCE CHECK of the RK4 route against the m=2 exact superoperator exponential,")
print("  and a third, assumption-free NAIVE matrix RK4 that uses none of the fast-path structure.")
print("  This third route exists because the first version of this script had the two disagree:")
print("  the superoperator was built with the column-major kron and fed a row-major reshape, an")
print("  error that is invisible for real symmetric rho_0 and wrong for complex ones.")

def naive_rhs(H, Ls, r):
    o = -1j * (H @ r - r @ H)
    for L in Ls:
        o += L @ r @ L.conj().T - 0.5 * (L.conj().T @ L @ r + r @ L.conj().T @ L)
    return o

def naive_rk4(H, Ls, r, T, n):
    dt = T / n
    for _ in range(n):
        k1 = naive_rhs(H, Ls, r)
        k2 = naive_rhs(H, Ls, r + 0.5 * dt * k1)
        k3 = naive_rhs(H, Ls, r + 0.5 * dt * k2)
        k4 = naive_rhs(H, Ls, r + dt * k3)
        r = r + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return r / np.trace(r).real

for nm in ("A", "D"):
    exA = obs(md, evolve_exact(H07, Ls_on, rho_of(ST[nm]), ts, D, ("on", 0.7, "A", 0.8)
                               if nm == "A" else ("cv", nm))[-1])[0]
    nv = obs(md, naive_rk4(H07, Ls_on, rho_of(ST[nm]).copy(), 40.0, 8000))[0]
    print(f"      init {nm} (rho_0 symmetry defect ||rho-rho^T|| = "
          f"{np.linalg.norm(rho_of(ST[nm]) - rho_of(ST[nm]).T):.3e}):"
          f"  exact expm = {exA:.10f}   naive RK4 = {nv:.10f}   diff = {abs(exA-nv):.3e}")
    for j, nst in enumerate((400, 1600, 3200)):
        r = rk4(md, H07, jumps_of(md), rho_of(ST[nm]), 40.0, nst, report=(j == 0 and nm == "A"))
        print(f"          fast RK4 steps = {nst:5d}   p+(40) = {obs(md, r)[0]:.10f}"
              f"   diff from exact = {abs(obs(md, r)[0]-exA):.3e}")

print("\n  m-scaling (junction A, g^2=0.7, gamma_E=0.5, Gamma_d=0.8).  Each m is run to T=40 AND to")
print("  T=160, because at larger m the cap has not finished evaporating by T=40 and a spread")
print("  that shrinks with m would otherwise be confounded with an unfinished run.")
print(f"  {'m':>3s} {'dim':>5s} {'T':>4s} " + " ".join(f"{'p+('+k+')':>10s}" for k in UNB) +
      f" {'SPREAD':>9s} {'<R_E>(A)':>10s} {'pcap':>9s}")
print("  " + "-" * 103)
for m in (2, 3, 4):
    x = MD[m]
    st = x.starts()
    Hx = x.H(0.7)
    jm = jumps_of(x)
    RT = np.diag(np.concatenate([x.RE1_all[0], x.RE0_all[0]])).astype(complex)
    Px = (np.eye(x.D) + RT) / 2.0
    Cx = np.zeros((x.D, x.D), complex)
    Cx[:x.D1, :x.D1] = np.eye(x.D1)
    for T, nst in ((40.0, 1600), (160.0, 6400)):
        vals, rA, pcA = {}, 0.0, 0.0
        for nm in UNB:
            r0 = rho_of(st[nm])
            p0 = np.trace(Px @ r0).real
            if abs(p0 - 0.5) > 1e-12:
                print(f"      WARNING m={m} start {nm} is not unbiased: p+(0) = {p0:.12f}")
            r = rk4(x, Hx, jm, r0, T, nst)
            vals[nm] = np.trace(Px @ r).real
            if nm == "A":
                rA = np.trace(RT @ r).real
                pcA = np.trace(Cx @ r).real
        print(f"  {m:3d} {x.D:5d} {T:4.0f} " + " ".join(f"{vals[k]:10.6f}" for k in UNB) +
              f" {max(vals.values())-min(vals.values()):9.6f} {rA:+10.6f} {pcA:9.2e}")

print("\n  m-scaling with the SYMMETRISED junction C (the null), same parameters:")
print(f"  {'m':>3s} {'dim':>5s} " + " ".join(f"{'p+('+k+')':>10s}" for k in UNB) + f" {'SPREAD':>9s}")
print("  " + "-" * 74)
for m in (2, 3, 4):
    x = MD[m]
    st = x.starts()
    Hx = x.H(0.7)
    jm = ([('d', np.sqrt(0.5) * d) for d in x.zt] +
          [('v', np.sqrt(0.4) * V) for V in x.Vs] +
          [('v', np.sqrt(0.4) * V) for V in x.Vs_flip])
    RT = np.diag(np.concatenate([x.RE1_all[0], x.RE0_all[0]])).astype(complex)
    Px = (np.eye(x.D) + RT) / 2.0
    vals = {}
    for nm in UNB:
        r = rk4(x, Hx, jm, rho_of(st[nm]), 40.0, 1600)
        vals[nm] = np.trace(Px @ r).real
    print(f"  {m:3d} {x.D:5d} " + " ".join(f"{vals[k]:10.6f}" for k in UNB) +
          f" {max(vals.values())-min(vals.values()):9.6f}")

print("\n" + "=" * 104)
print("END OF RUN.  Every prediction above is labelled as a prediction and was fixed before the")
print("number beside it was computed.  No verdict is hard-coded.")
print("=" * 104)
