"""W-31.  THE EVAPORATING CAP -- a responsive carrier, run against its own worst controls.

CONSTRUCTION (one line).  Make the carrier's FIRST BETTI NUMBER dynamical.  A Z_2 gauge field
lives on an annulus of m quads whose central disk-CAP (hub vertex + m spokes) can evaporate.
n=1: the complex is a DISK, b1=0.  n=0: it is an ANNULUS, b1=1.  The record is the winding
parity of the electric flux around the hole, R_E = (-1)^(s_{a0}+s_{c0}).  While the cap exists
the cap triangles fail to commute with R_E, so R_E is WRITABLE.  When the cap goes, the only
operator that failed to commute with R_E is physically removed, so R_E is EXACTLY conserved.

FRAMING CORRECTION (adopted from the losing proposals, and it matters).
  * The W-30 "no-go" is not a theorem that anything breaks.  R unitary and [L,R]=0 gives
    [L^dag,R]=0 and hence Tr(R D[rho]) = 0 for ANY L, unitary or not -- it is a tautology about
    superselection: a conserved charge is conserved.  This script verifies the non-unitary case
    explicitly (STEP 4b) so the lane stops treating unitarity as the load-bearing hypothesis.
  * W-30's R = prod_p W_p is, by discrete Stokes on a DISK, a function of the local plaquettes:
    it is the dual-Ising spin-flip GENERATOR, not an order parameter.  Generators are conserved,
    hence unwritable.  That is exactly and only what W-30 measured.
  * This construction does NOT break that.  It makes R_E non-conserved for a finite window and
    then deletes the term that made it so.  R_E is a HOLONOMY of a nontrivial cycle, which only
    exists once the hole is open; that is what opening the hole buys.

WHAT IS AT RISK OF BEING CONSTRAINT-FORCED (declared before any dynamics, STEP 2).
  p_+(infinity), the frozen weight of the R_E=+1 sector.  Three ways it could be forced:
   (i)  unequal sector dimensions -- counting alone would fix p_+;
   (ii) the evaporation GATE preferring a sector -- this is the disqualified (1+R)/2 projector
        smuggled into the geometry;
   (iii) at large g^2 the state is pinned near the flux vacuum where the winding is trivially 0.
  (i) and (ii) are settled by two traces before any time step: tr(R_E * P_punctured) and
  tr(R_E * sum_k L_k^dag L_k) for EVERY junction family.  (iii) is not settled by any commutator
  and is left to the g^2 scan, where it is expected to appear as a collapse of the spread.

THE LOAD-BEARING CHOICE, AND IT IS RUN HERE AS A CONTROL, NOT ADMIRED (STEP 7).
  Topology change needs a junction condition.  The two natural clearings of the vanishing spokes
  differ by the kernel {0, sum_i tau_i}, and sum_i tau_i IS the inner cycle, whose winding is 1 --
  so the ambiguity in the junction condition is EXACTLY the record.  Three junction families are
  therefore implemented and all three are run:
    A  NO-HAIR   : one jump, sqrt(Gd) iota^dag.  The cap may only vanish where the field put no
                   flux on it.  Field-gated.
    B  DRAG/BLIND: a full Kraus family {sqrt(Gd) V_p}, sum_p V_p^dag V_p = I.  The cap carries its
                   flux away, at the SAME rate from every state -- evaporation blind to the field.
                   Branch convention T_0 = 0.
    B' FLIPPED   : the other branch, T_0 = 1.  Every V_p carries an extra M_inner.
    C  SYMMETRISED: half rate on each branch.  If the correct treatment symmetrises, the design
                   predicts R_E is randomised at every topology change and the record dies.
  Nothing derives the choice.  The numbers below say what each choice gives.

ROUTE.  m=2 (dim 24, superoperator 576x576): EXACT Lindblad superoperator exponentiation, expm
implemented here by scaling-and-squaring Taylor (scipy absent).  m=3 (dim 80) and m=4 (dim 288):
RK4 on rho, validated against the m=2 exact result.  numpy only.

GAUGE STRUCTURE, HONESTLY.  Within each fiber the full Z_2 Gauss law holds exactly (the fiber IS
the divergence-free flux sector of that graph).  ACROSS fibers there is no single gauge group --
the hub's Gauss law exists only when the hub exists -- so this is a bundle of constraint surfaces,
one per carrier configuration.  A SECOND CONSTRAINT is added for junction A and must be written
down: the NO-HAIR JUNCTION CONDITION, E_sigma_i |psi> = 0 on the vanishing links, implemented as
iota with iota iota^dag = Pi_nohair.  Junctions B/B'/C drop it and are run for that reason.
"""

import itertools, numpy as np

np.set_printoptions(linewidth=200)

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
    """Divergence-free Z_2 flux configs: the Gauss-law physical sector, exactly as w27_patch.py."""
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
    """M_p |s> = |s + p>.  W-27/W-30 'Move' convention, flux basis."""
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

        # record operator: winding parity across a radial cut gamma_j = {a_j, c_j}
        self.RE1_all = [np.array([(-1.0) ** (s[a(j)] + s[c(j)]) for s in self.S1]) for j in range(m)]
        self.RE0_all = [np.array([(-1.0) ** (s[a(j)] + s[c(j)]) for s in self.S0]) for j in range(m)]
        self.rd = np.concatenate([self.RE1_all[0], self.RE0_all[0]]).astype(complex)  # diagonal of R_E
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
        self.Pnh = iota @ iota.conj().T                      # Pi_nohair on the capped fiber
        self.IOTA = np.zeros((D, D), complex)
        self.IOTA[:D1, D1:] = iota                           # punctured -> capped

        # ---- junction families ------------------------------------------------
        spoke_pat = [tuple(s[sg(i)] for i in range(m)) for s in self.S1]
        pats = sorted(set(spoke_pat))
        self.pats = pats
        Vs = []
        for p in pats:
            T = [0] * m
            for j in range(1, m):
                T[j] = T[j - 1] ^ p[j - 1]
            U = np.eye(D1, dtype=complex)
            for i in range(m):
                if T[i]:
                    U = U @ self.Mt1[i]
            Pi = np.diag([1.0 if spoke_pat[j] == p else 0.0 for j in range(D1)]).astype(complex)
            Vs.append(iota.conj().T @ U @ Pi)                # D0 x D1
        self.Vs = Vs
        self.Vs_flip = [V @ self.Minner1 for V in Vs]

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
        m, D1, D0 = self.m, self.D1, self.D0
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


# ----------------------------------------------------------------------------------
# Lindblad machinery
# ----------------------------------------------------------------------------------
def superop(H, Ls, D):
    Id = np.eye(D, dtype=complex)
    M = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for L in Ls:
        LdL = L.conj().T @ L
        M += np.kron(L.conj(), L) - 0.5 * (np.kron(Id, LdL) + np.kron(LdL.T, Id))
    return M


def evolve_exact(H, Ls, rho0, ts, D):
    """uniform grid ts; one expm, then repeated application."""
    dt = ts[1] - ts[0]
    U = expm(superop(H, Ls, D) * dt)
    v = rho0.reshape(-1).astype(complex)
    out = []
    for n in range(len(ts)):
        if n > 0:
            v = U @ v
        r = v.reshape(D, D)
        r = r / np.trace(r).real
        out.append(r)
    return out


def obs(md, r):
    return (np.trace(md.Pp @ r).real,
            np.linalg.norm(md.Pp @ r @ md.Pm),
            np.trace(md.PCAP @ r).real,
            np.trace(md.RE @ r).real)


def rho_of(v):
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


# ----------------------------------------------------------------------------------
print("=" * 100)
print("W-31  THE EVAPORATING CAP -- responsive carrier.  numpy only, no scipy.")
print("=" * 100)

# ==================================================================================
# STEP 1  --  THE SPACE, AND THE STATIC COMPARISON CASE
# ==================================================================================
print("\n[STEP 1] THE SPACE\n")
MD = {}
for m in (2, 3, 4):
    MD[m] = Model(m)
    md = MD[m]
    print(f"  m={m}   capped fiber (DISK, b1=0) dim {md.D1:4d} (pred 2^{2*m}={2**(2*m):4d})"
          f"   punctured fiber (ANNULUS, b1=1) dim {md.D0:4d} (pred 2^{m+1}={2**(m+1):4d})"
          f"   TOTAL {md.D:4d}   superoperator {md.D**2}x{md.D**2}")
md = MD[2]

# static comparison case: W-30's own carrier, built here from scratch
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
      f" physical dim {Dp}   superoperator {Dp**2}x{Dp**2}")
print(f"  ROUTE: m=2 exact superoperator exponential (576x576); m=3,4 RK4 on rho, cross-checked"
      f" against m=2 exact.")

# ==================================================================================
# STEP 2  --  FORCED-OR-NOT, RUN FIRST
# ==================================================================================
print("\n[STEP 2] FORCED-OR-NOT CHECK  (run before any dynamics)\n")
print("  QUANTITY AT RISK: p_+(inf), the frozen weight of the R_E=+1 winding sector.")
print("  Forced if (i) the sectors have unequal dimension, or (ii) any evaporation gate")
print("  sum_k L_k^dag L_k has nonzero overlap with R_E.  Counting argument declared in advance:")
print("  w(s)=s_{a0}+s_{c0} is a LINEAR functional on the Z_2 cycle space and is surjective")
print("  (w(inner cycle)=1), so its kernel has index exactly 2 and the sectors are equal.\n")
forced = False
for m in (2, 3, 4):
    x = MD[m]
    n1p = int(np.sum(np.real(x.RE1_all[0]) > 0))
    n1m = x.D1 - n1p
    n0p = int(np.sum(np.real(x.RE0_all[0]) > 0))
    n0m = x.D0 - n0p
    tr0 = float(np.sum(x.RE0_all[0]))
    tr1 = float(np.sum(x.RE1_all[0]))
    tnh = float(np.real(np.trace(np.diag(x.RE1_all[0]).astype(complex) @ x.Pnh)))
    print(f"  m={m}  punctured sectors {n0p}/{n0m}   tr(R_E*P_punc) = {tr0:+.1f}"
          f"    capped sectors {n1p}/{n1m}   tr(R_E*P_cap) = {tr1:+.1f}"
          f"    tr(R_E*Pi_nohair) = {tnh:+.1f}")
    for nm in ("A", "B", "Bp", "C"):
        Ls = x.junction(nm, 1.0)
        G = sum(L.conj().T @ L for L in Ls)
        ov = float(np.real(np.trace(np.diag(np.concatenate([x.RE1_all[0], x.RE0_all[0]])).astype(complex) @ G)))
        print(f"        junction {nm:2s}: sum_k L^dag L  ->  tr(R_E * gate) = {ov:+.3e}"
              f"   ||sum L^dag L - Gd*Pi|| info: trace(gate) = {float(np.real(np.trace(G))):.1f}")
        if abs(ov) > 1e-9:
            forced = True
    if abs(tr0) > 1e-9 or abs(tnh) > 1e-9:
        forced = True
print(f"\n  ANY nonzero overlap found (record constraint-forced / gate biased)?  {forced}")
if forced:
    print("  THE RECORD OBSERVABLE IS CONSTRAINT-FORCED OR THE GATE IS BIASED.  STOPPING.")
    raise SystemExit(0)
print("  (residual risk (iii), large-g^2 pinning, is NOT settled by any commutator -- see STEP 9)")

# ==================================================================================
# STEP 3  --  DOES R EARN THE RIGHT TO BE MEASURED?
# ==================================================================================
print("\n[STEP 3] THE RECORD OPERATOR MUST EARN ITS MEASUREMENT\n")
for m in (2, 3, 4):
    x = MD[m]
    R = x.RE
    nR = np.linalg.norm(R)
    uni = np.linalg.norm(R @ R.conj().T - np.eye(x.D))
    ev = np.linalg.eigvals(R)
    nd = len(np.unique(np.round(ev, 9)))
    gi0 = max(np.linalg.norm(x.RE0_all[j] - x.RE0_all[0]) for j in range(m))
    gi1 = max(np.linalg.norm(x.RE1_all[j] - x.RE1_all[0]) for j in range(m))
    print(f"  m={m}  ||R_E|| = {nR:8.4f}   ||R R^dag - I|| = {uni:.3e}   distinct eigenvalues = {nd}"
          f"   cut-dependence: punctured {gi0:.3e}   capped {gi1:.3e}")

# ==================================================================================
# STEP 4  --  REPRODUCE THE STATIC NO-GO INSIDE THIS SCRIPT
# ==================================================================================
print("\n[STEP 4] STATIC NO-GO CONTROLS  (if these do not reproduce, nothing below means anything)\n")

# --- 4a: literal W-30 replica on the 3x3 patch -------------------------------------
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
Hp = -sum(L + L.conj().T for L in Pl) - 1.0 * sum(2 * np.diag(zdiag(Sp, k)) for k in range(len(Ep)))
worst_u = 0.0
for _ in range(20):
    A = rng.normal(size=(Dp, Dp)) + 1j * rng.normal(size=(Dp, Dp))
    rr = A @ A.conj().T
    rr = rr / np.trace(rr).real
    d = -1j * (Hp @ rr - rr @ Hp) + sum(L @ rr @ L.conj().T - rr for L in Pl)
    worst_u = max(worst_u, abs(np.trace(Rrim @ d)))
print(f"  4a  W-30 REPLICA, 3x3 patch dim {Dp}, R = product of the 4 plaquettes (discrete Stokes):")
print(f"      ||R|| = {np.linalg.norm(Rrim):.4f}   unitarity defect = "
      f"{np.linalg.norm(Rrim@Rrim.conj().T-np.eye(Dp)):.3e}   distinct eigenvalues = "
      f"{len(np.unique(np.round(np.linalg.eigvals(Rrim),9)))}")
print(f"      ||[H,R]|| = {np.linalg.norm(Hp@Rrim-Rrim@Hp):.3e}   "
      f"max||[L,R]|| = {max(np.linalg.norm(L@Rrim-Rrim@L) for L in Pl):.3e}   "
      f"max |d<R>/dt| over 20 random rho = {worst_u:.3e}")

# --- 4b: the tautology graft -- non-unitary L too -----------------------------------
Bc = rng.normal(size=(Dp, Dp)) + 1j * rng.normal(size=(Dp, Dp))
Pplus_p = (np.eye(Dp) + Rrim) / 2.0
Pminus_p = (np.eye(Dp) - Rrim) / 2.0
Lcom = Pplus_p @ Bc @ Pplus_p + Pminus_p @ Bc @ Pminus_p       # commutes with Rrim, NOT unitary
Lnon = Pplus_p @ Bc @ Pminus_p                                  # does not commute
def dR(L, R, H):
    w = 0.0
    for _ in range(10):
        A = rng.normal(size=(Dp, Dp)) + 1j * rng.normal(size=(Dp, Dp))
        rr = A @ A.conj().T
        rr = rr / np.trace(rr).real
        d = -1j * (H @ rr - rr @ H) + (L @ rr @ L.conj().T
                                       - 0.5 * (L.conj().T @ L @ rr + rr @ L.conj().T @ L))
        w = max(w, abs(np.trace(R @ d)))
    return w
print(f"  4b  TAUTOLOGY GRAFT (the no-go does not need unitary jumps):")
print(f"      L non-unitary, [L,R]=0 : ||L L^dag - I|| = "
      f"{np.linalg.norm(Lcom@Lcom.conj().T-np.eye(Dp)):.3e}   ||[L,R]|| = "
      f"{np.linalg.norm(Lcom@Rrim-Rrim@Lcom):.3e}   max|d<R>/dt| = {dR(Lcom,Rrim,Hp):.3e}")
print(f"      L non-unitary, [L,R]!=0: ||[L,R]|| = "
      f"{np.linalg.norm(Lnon@Rrim-Rrim@Lnon):.3e}   max|d<R>/dt| = {dR(Lnon,Rrim,Hp):.3e}")

# --- 4c/4d/4e: the same disease inside this construction ----------------------------
print(f"  4c  ELITZUR / STOKES CONTRAST inside this construction "
      f"(is the record a GENERATOR or a HOLONOMY?):")
for m in (2, 3, 4):
    x = MD[m]
    pr = np.eye(x.D1, dtype=complex)
    for M in x.Mq1 + x.Mt1:
        pr = pr @ M
    out1 = move(x.S1, x.i1, [x.c(i) for i in range(m)])
    pr0 = np.eye(x.D0, dtype=complex)
    for M in x.Mq0:
        pr0 = pr0 @ M
    out0 = move(x.S0, x.i0, [x.c(i) for i in range(m)])
    print(f"      m={m}  ||prod(all plaquettes) - M_outer||   capped(DISK) = "
          f"{np.linalg.norm(pr-out1):.3e}    punctured(ANNULUS) = {np.linalg.norm(pr0-out0):.3e}")

print(f"  4d  COMMUTATORS OF R_E WITH THE FULL GENERATOR (g^2=0.7, kappa=1):")
for m in (2, 3, 4):
    x = MD[m]
    H1, H0 = x.Hfib(0.7)
    R1 = np.diag(x.RE1_all[0]).astype(complex)
    R0 = np.diag(x.RE0_all[0]).astype(complex)
    cL = 0.0
    for L in x.elec(0.5) + x.junction("A", 0.8):
        cL = max(cL, np.linalg.norm(x.RE @ L - L @ x.RE))
    print(f"      m={m}  ||[R_E,H_capped]|| = {np.linalg.norm(R1@H1-H1@R1):9.4f}   "
          f"||[R_E,H_punctured]|| = {np.linalg.norm(R0@H0-H0@R0):.3e}   "
          f"max||[R_E,L_k]|| over ALL jumps (junction A) = {cL:.3e}")

print(f"  4e  CARRIER RESPONSE OFF (Gamma_d=0, lambda=0).  START PUNCTURED: the carrier is now")
print(f"      static and R_E is a conserved charge -- d<R_E>/dt must be 0 to machine precision.")
D = md.D
ts = np.linspace(0, 40, 41)
vpun = np.zeros(D, complex)
vpun[md.D1:] = rng.normal(size=md.D0) + 1j * rng.normal(size=md.D0)
r_pun = rho_of(vpun)
for g2 in (0.0, 0.7, 3.0):
    H = md.H(g2)
    tr = evolve_exact(H, md.elec(0.5), r_pun, ts, D)
    o0, o1 = obs(md, tr[0]), obs(md, tr[-1])
    print(f"        g^2={g2:4.1f}  p+(0) = {o0[0]:.12f}   p+(40) = {o1[0]:.12f}   "
          f"|delta| = {abs(o1[0]-o0[0]):.3e}   <R_E>(40) = {o1[3]:+.12f}")
print(f"      START CAPPED (carrier still frozen): the writer never stops, so nothing freezes.")
zero = tuple([0] * (4 * md.m))
inner = tuple([1 if k < md.m else 0 for k in range(4 * md.m)])
q0 = [0] * (4 * md.m)
for k in md.quad(0):
    q0[k] ^= 1
q0 = tuple(q0)
q0a = list(q0)
for k in range(md.m):
    q0a[k] ^= 1
q0a = tuple(q0a)
q1 = [0] * (4 * md.m)
for k in md.quad(1 % md.m):
    q1[k] ^= 1
q1 = tuple(q1)
q1a = list(q1)
for k in range(md.m):
    q1a[k] ^= 1
q1a = tuple(q1a)

def ketc(s):
    v = np.zeros(D, complex)
    v[md.i1[s]] = 1.0
    return v

psi = {
    "A": (ketc(zero) + ketc(inner)) / np.sqrt(2),
    "B": (ketc(q0) + ketc(q0a)) / np.sqrt(2),
    "D": (ketc(q1) + ketc(q1a)) / np.sqrt(2),
}
vC = np.zeros(D, complex)
vC[:md.D1] = rng.normal(size=md.D1) + 1j * rng.normal(size=md.D1)
psi["C"] = vC / np.linalg.norm(vC)
# an initial state deliberately biased, to see whether bias survives with the carrier frozen
vE = np.zeros(D, complex)
vE[:md.D1] = rng.normal(size=md.D1) + 1j * rng.normal(size=md.D1)
vE = md.Pp @ vE + 0.4 * (md.Pm @ vE)
psi["E"] = vE / np.linalg.norm(vE)

H07 = md.H(0.7)
for nm in ("A", "B", "C", "E"):
    tr = evolve_exact(H07, md.elec(0.5), rho_of(psi[nm]), ts, D)
    o = [obs(md, r) for r in tr]
    late = np.mean([o[k][0] for k in range(30, 41)])
    print(f"        init {nm}: p+(0) = {o[0][0]:.6f}   p+(40) = {o[-1][0]:.6f}   "
          f"late-avg p+[30,40] = {late:.6f}   drift[30,40] = {abs(o[-1][0]-o[-11][0]):.3e}   "
          f"p_cap(40) = {o[-1][2]:.3e}")

# ==================================================================================
# STEP 5  --  CARRIER RESPONSE ON
# ==================================================================================
print("\n[STEP 5] CARRIER RESPONSE ON  (junction A, no-hair gate, Gamma_d=0.8, gamma_E=0.5,"
      " g^2=0.7)\n")
print("  DECOY WARNING, adopted in advance from Proposal 3: for a pure state")
print("  ||P+ psi psi^dag P-||_1 = (1/2) sqrt(1 - <R>^2), so the coherence column CANNOT by itself")
print("  witness selection -- phase randomisation drives it to zero with no weight transfer.")
print("  The coherence column below is a consistency check only.  The evidence is p_+ and its")
print("  dependence on the initial state.\n")
Ls_on = md.elec(0.5) + md.junction("A", 0.8)
print(f"  {'init':>5s} {'p+(0)':>10s} {'<R_E>(0)':>10s} {'p+(2)':>10s} {'p+(5)':>10s} "
      f"{'p+(10)':>10s} {'p+(20)':>10s} {'p+(40)':>10s} {'late-avg':>10s} {'drift':>10s} "
      f"{'coh(40)':>9s} {'pcap(40)':>9s}")
print("  " + "-" * 128)
p_on = {}
for nm in ("A", "B", "D", "C", "E"):
    tr = evolve_exact(H07, Ls_on, rho_of(psi[nm]), ts, D)
    o = [obs(md, r) for r in tr]
    late = np.mean([o[k][0] for k in range(30, 41)])
    p_on[nm] = late
    print(f"  {nm:>5s} {o[0][0]:10.6f} {o[0][3]:+10.6f} {o[2][0]:10.6f} {o[5][0]:10.6f} "
          f"{o[10][0]:10.6f} {o[20][0]:10.6f} {o[40][0]:10.6f} {late:10.6f} "
          f"{abs(o[-1][0]-o[-11][0]):10.3e} {o[-1][1]:9.2e} {o[-1][2]:9.2e}")
print(f"\n  <R_E>(t) for init A: " +
      "  ".join(f"t={ts[k]:.0f}:{obs(md, evolve_exact(H07, Ls_on, rho_of(psi['A']), ts, D)[k])[3]:+.6f}"
               for k in (0, 1, 2, 5, 10, 20, 40)))
print(f"  SPREAD over the three UNBIASED starts (p+(0)=1/2 exactly): "
      f"max-min over A,B,D = {max(p_on[k] for k in 'ABD')-min(p_on[k] for k in 'ABD'):.6f}")
print(f"  SPREAD including the two BIASED starts C,E (p+(0) = "
      f"{obs(md, rho_of(psi['C']))[0]:.4f}, {obs(md, rho_of(psi['E']))[0]:.4f}): "
      f"max-min over A,B,D,C,E = "
      f"{max(p_on[k] for k in 'ABDCE')-min(p_on[k] for k in 'ABDCE'):.6f}")

# ==================================================================================
# STEP 6  --  LINDBLADIAN SPECTRUM: LIFETIME, NOT A PLATEAU
# ==================================================================================
print("\n[STEP 6] LINDBLADIAN SPECTRUM  (adopted from Proposal 2: report a lifetime, not a drift)\n")
for nm in ("A", "B", "C"):
    Ls = md.elec(0.5) + md.junction(nm, 0.8)
    Lsup = superop(H07, Ls, D)
    w, Vv = np.linalg.eig(Lsup)
    nk = int(np.sum(np.abs(w) < 1e-8))
    ov = []
    for k in range(len(w)):
        Vk = Vv[:, k].reshape(D, D)
        nv = np.linalg.norm(Vk)
        if nv < 1e-12:
            continue
        o = abs(np.trace(md.RE @ Vk)) / nv
        if o > 1e-6 and w[k].real < -1e-9:
            ov.append((w[k].real, o))
    ov.sort(key=lambda t: -t[0])
    slow = ov[0][0] if ov else None
    print(f"  junction {nm:2s}: dim ker(Lindbladian) = {nk:3d}  (>1 means a multi-dimensional"
          f" steady-state manifold: information is retained)")
    if slow is None:
        print(f"              NO decaying mode has nonzero overlap with R_E  -> tau_erase = infinite"
              f" within numerical resolution")
    else:
        print(f"              slowest decaying mode with nonzero R_E overlap: Re(lambda) = {slow:.3e}"
              f"  -> tau_erase = {-1.0/slow:.3e}")
        print(f"              next three: " + ", ".join(f"{r:.3e}" for r, _ in ov[1:4]))

# ==================================================================================
# STEP 7  --  THE JUNCTION CONTROLS (the load-bearing choice, run not admired)
#             and the NULL CONTROLS
# ==================================================================================
print("\n[STEP 7] JUNCTION CONTROLS AND NULLS  (g^2=0.7, gamma_E=0.5, Gamma_d=0.8)\n")
print(f"  {'junction':>9s} {'#jumps':>7s} {'init':>5s} {'p+(0)':>9s} {'late-avg p+':>12s} "
      f"{'drift':>10s} {'<R_E>late':>11s} {'pcap(40)':>9s}")
print("  " + "-" * 90)
junc_res = {}
for jn in ("A", "B", "Bp", "C"):
    Ls = md.elec(0.5) + md.junction(jn, 0.8)
    for nm in ("A", "B", "D"):
        tr = evolve_exact(H07, Ls, rho_of(psi[nm]), ts, D)
        o = [obs(md, r) for r in tr]
        late = np.mean([o[k][0] for k in range(30, 41)])
        lateR = np.mean([o[k][3] for k in range(30, 41)])
        junc_res[(jn, nm)] = late
        print(f"  {jn:>9s} {len(Ls):7d} {nm:>5s} {o[0][0]:9.6f} {late:12.6f} "
              f"{abs(o[-1][0]-o[-11][0]):10.3e} {lateR:+11.6f} {o[-1][2]:9.2e}")
    sp = max(junc_res[(jn, k)] for k in "ABD") - min(junc_res[(jn, k)] for k in "ABD")
    print(f"  {'':>9s} {'':>7s} SPREAD over unbiased starts A,B,D = {sp:.6f}")
print("\n  DECLARED PREDICTIONS FOR THIS TABLE, made before it was run:")
print("    A  no-hair gate           -> record forms (spread > 0)")
print("    B  blind/drag, branch T0=0 -> if the spread survives, the field-gating is DECORATION")
print("    Bp blind/drag, branch T0=1 -> mirrored record")
print("    C  symmetrised             -> R_E randomised at every topology change, spread -> 0")

print("\n  NULL CONTROL N1 -- start ALREADY punctured, full bath, long time.  Must be exactly flat.")
for T in (40.0, 200.0):
    tsl = np.linspace(0, T, 5)
    tr = evolve_exact(md.H(3.0), md.elec(0.5) + md.junction("A", 0.8), r_pun, tsl, D)
    o0, o1 = obs(md, tr[0]), obs(md, tr[-1])
    print(f"      t=0 p+ = {o0[0]:.12f}   t={T:.0f} p+ = {o1[0]:.12f}   |delta| = "
          f"{abs(o1[0]-o0[0]):.3e}")

print("\n  NULL CONTROL N2 -- kill the record operator's writer (kappa=0): nothing should ever")
print("      be written, so every initial state must freeze at its own starting value.")
Hk0 = md.H(0.7, kap=0.0)
for nm in ("A", "B", "D"):
    tr = evolve_exact(Hk0, Ls_on, rho_of(psi[nm]), ts, D)
    o = [obs(md, r) for r in tr]
    print(f"      init {nm}: p+(0) = {o[0][0]:.9f}   p+(40) = {o[-1][0]:.9f}   |delta| = "
          f"{abs(o[-1][0]-o[0][0]):.3e}")

print("\n  NULL CONTROL N3 -- p_cap(t) is NOT evidence.  Reported here so it is never quoted:")
for g2 in (0.0, 0.7, 3.0):
    tr = evolve_exact(md.H(g2), md.elec(0.5) + md.junction("A", 0.8), rho_of(psi["A"]), ts, D)
    trB = evolve_exact(md.H(g2), md.elec(0.5) + md.junction("B", 0.8), rho_of(psi["A"]), ts, D)
    print(f"      g^2={g2:4.1f}  p_cap(10) gate A = {obs(md,tr[10])[2]:.6e}   gate B (blind) = "
          f"{obs(md,trB[10])[2]:.6e}")

# ==================================================================================
# STEP 8  --  SCAN THE RESPONSE STRENGTH, INCLUDING EXACTLY ZERO
# ==================================================================================
print("\n[STEP 8] SCAN OF THE CARRIER RESPONSE STRENGTH Gamma_d  (g^2=0.7, gamma_E=0.5)\n")
print(f"  {'Gamma_d':>8s} {'junc':>5s} {'late p+ (A)':>12s} {'late p+ (B)':>12s} "
      f"{'late p+ (D)':>12s} {'SPREAD':>9s} {'drift(A)':>10s} {'pcap(40)':>9s}")
print("  " + "-" * 86)
for jn in ("A", "B", "C"):
    for Gd in (0.0, 0.05, 0.2, 0.8, 3.2):
        Ls = md.elec(0.5) + md.junction(jn, Gd)
        vals = {}
        dr = 0.0
        pc = 0.0
        for nm in ("A", "B", "D"):
            tr = evolve_exact(H07, Ls, rho_of(psi[nm]), ts, D)
            o = [obs(md, r) for r in tr]
            vals[nm] = np.mean([o[k][0] for k in range(30, 41)])
            if nm == "A":
                dr = abs(o[-1][0] - o[-11][0])
                pc = o[-1][2]
        sp = max(vals.values()) - min(vals.values())
        print(f"  {Gd:8.2f} {jn:>5s} {vals['A']:12.6f} {vals['B']:12.6f} {vals['D']:12.6f} "
              f"{sp:9.6f} {dr:10.3e} {pc:9.2e}")

# ==================================================================================
# STEP 9  --  THE ALPHA AXIS: g^2 SCAN, BOTH SPREAD METRICS
# ==================================================================================
print("\n[STEP 9] g^2 SCAN  (junction A, Gamma_d=0.8, gamma_E=0.5).  BOTH metrics are printed:")
print("  the spread over UNBIASED starts (the honest one) and the max-min including the two")
print("  BIASED starts C,E -- because a biased start that is never written simply retains its")
print("  bias, which inflates max-min at large g^2 exactly where nothing is being written.\n")
print(f"  {'g^2':>7s} {'||[R,Hcap]||':>13s} {'p+(A)':>10s} {'p+(B)':>10s} {'p+(D)':>10s} "
      f"{'p+(C)':>10s} {'p+(E)':>10s} {'SPREAD(ABD)':>12s} {'maxmin(ABDCE)':>14s}")
print("  " + "-" * 106)
H1r = np.diag(md.RE1_all[0]).astype(complex)
for g2 in (0.0, 0.1, 0.3, 0.7, 1.0, 1.5, 3.0, 8.0, 20.0):
    Hg = md.H(g2)
    Hc = Hg[:md.D1, :md.D1]
    cw = np.linalg.norm(H1r @ Hc - Hc @ H1r)
    Ls = md.elec(0.5) + md.junction("A", 0.8)
    vals = {}
    for nm in ("A", "B", "D", "C", "E"):
        tr = evolve_exact(Hg, Ls, rho_of(psi[nm]), ts, D)
        o = [obs(md, r) for r in tr]
        vals[nm] = np.mean([o[k][0] for k in range(30, 41)])
    s3 = max(vals[k] for k in "ABD") - min(vals[k] for k in "ABD")
    s5 = max(vals.values()) - min(vals.values())
    print(f"  {g2:7.2f} {cw:13.4f} {vals['A']:10.6f} {vals['B']:10.6f} {vals['D']:10.6f} "
          f"{vals['C']:10.6f} {vals['E']:10.6f} {s3:12.6f} {s5:14.6f}")

# ==================================================================================
# STEP 10 --  UNRAVELLING: IS THE FROZEN VALUE PER-BRANCH DEFINITE?
# ==================================================================================
print("\n[STEP 10] TRAJECTORY UNRAVELLING  (adopted from Proposal 3).  The ensemble coherence")
print("  going to zero is compatible with nothing having been selected.  The per-branch witness")
print("  is E[<R_E>^2] over trajectories and the fraction of trajectories with |<R_E>| > 0.99.\n")

def unravel(H, Ls, psi0, T, dt, ntraj, seed):
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
        vals.append(float(np.real(np.vdot(v, md.RE @ v))))
    return np.array(vals)

for jn in ("A", "C"):
    Ls = md.elec(0.5) + md.junction(jn, 0.8)
    vv = unravel(H07, Ls, psi["A"], 25.0, 0.0125, 120, 11)
    print(f"  junction {jn:2s}, unitary(infinite-T) electric bath, 120 trajectories to t=25:")
    print(f"      E[<R_E>] = {vv.mean():+.6f}   E[<R_E>^2] = {np.mean(vv**2):.6f}   "
          f"frac |<R_E>|>0.99 = {np.mean(np.abs(vv)>0.99):.3f}   "
          f"frac |<R_E>|>0.5 = {np.mean(np.abs(vv)>0.5):.3f}")

# ==================================================================================
# STEP 11 --  FINITE-TEMPERATURE DAVIES BATH (the infinite-T defect of the whole lane)
# ==================================================================================
print("\n[STEP 11] FINITE-TEMPERATURE DAVIES BATH  (adopted from Proposal 2).  Every unitary jump")
print("  operator in this lane, including our own L = sqrt(gamma_E) Z_l, makes the identity")
print("  stationary: an INFINITE-temperature bath.  Replace it with a Davies bath and re-run.\n")

def davies(Hm, As, gamma0, Tmp, dim, tol=1e-7):
    w, U = np.linalg.eigh(Hm)
    Ls = []
    for A in As:
        Ab = U.conj().T @ A @ U
        freqs = {}
        for aa in range(dim):
            for bb in range(dim):
                om = w[aa] - w[bb]
                key = round(om / tol) * tol
                freqs.setdefault(key, []).append((aa, bb))
        for om, pairs in freqs.items():
            M = np.zeros((dim, dim), complex)
            for aa, bb in pairs:
                M[aa, bb] = Ab[aa, bb]
            if np.linalg.norm(M) < 1e-12:
                continue
            if abs(om) < 1e-12:
                g = gamma0 * Tmp
            else:
                g = gamma0 * om / (1.0 - np.exp(-om / Tmp))
            Ls.append(np.sqrt(abs(g)) * (U @ M @ U.conj().T))
    return Ls

# validate on the closed punctured fiber: the Gibbs state must be stationary
H1c, H0c = md.Hfib(0.7)
As0 = [np.diag(zdiag(md.S0, k)).astype(complex) for k in range(3 * md.m)]
for Tmp in (0.5, 2.0, 10.0):
    Ld = davies(H0c, As0, 0.5, Tmp, md.D0)
    w0, U0 = np.linalg.eigh(H0c)
    gb = U0 @ np.diag(np.exp(-w0 / Tmp)) @ U0.conj().T
    gb = gb / np.trace(gb).real
    dd = -1j * (H0c @ gb - gb @ H0c) + sum(L @ gb @ L.conj().T
                                           - 0.5 * (L.conj().T @ L @ gb + gb @ L.conj().T @ L)
                                           for L in Ld)
    R0d = np.diag(md.RE0_all[0]).astype(complex)
    cR = max(np.linalg.norm(R0d @ L - L @ R0d) for L in Ld)
    print(f"  T={Tmp:5.2f}  #Davies jumps on the punctured fiber = {len(Ld):4d}   "
          f"||D[rho_Gibbs]|| = {np.linalg.norm(dd):.3e}   max||[R_E, A(omega)]|| = {cR:.3e}")

# full-space Davies + evaporation
Ast = [np.diag(d) for d in md.zt]
print(f"\n  FULL RUN with a Davies electric bath (junction A, Gamma_d=0.8, g^2=0.7):")
print(f"  {'T':>7s} {'#jumps':>7s} {'p+(0)':>10s} {'late p+ (A)':>12s} {'late p+ (B)':>12s} "
      f"{'late p+ (D)':>12s} {'SPREAD':>9s} {'drift':>10s}")
print("  " + "-" * 84)
for Tmp in (0.5, 2.0, 10.0, 1e6):
    Ld = davies(H07, Ast, 0.5, Tmp, D)
    Ls = Ld + md.junction("A", 0.8)
    vals = {}
    dr = 0.0
    for nm in ("A", "B", "D"):
        tr = evolve_exact(H07, Ls, rho_of(psi[nm]), ts, D)
        o = [obs(md, r) for r in tr]
        vals[nm] = np.mean([o[k][0] for k in range(30, 41)])
        if nm == "A":
            dr = abs(o[-1][0] - o[-11][0])
            p0 = o[0][0]
    sp = max(vals.values()) - min(vals.values())
    print(f"  {Tmp:7.1f} {len(Ls):7d} {p0:10.6f} {vals['A']:12.6f} {vals['B']:12.6f} "
          f"{vals['D']:12.6f} {sp:9.6f} {dr:10.3e}")

# ==================================================================================
# STEP 12 --  SIZE SCALING.  The protection is CLAIMED homological, so it must be m-independent.
# ==================================================================================
print("\n[STEP 12] SIZE SCALING  m = 2,3,4.  m=2 exact; m=3,4 RK4 on rho with a convergence check.")
print("  If the frozen p_+ or the spread drifts systematically with m, the topological story is")
print("  decoration on a finite-size accident.\n")

def rk4(H, jumps, rho0, T, nst):
    """jumps: list of ('d', diag_vector) or ('m', matrix).  RK4 on rho."""
    dt = T / nst
    Sgg = np.zeros(H.shape[0], complex)
    Smat = np.zeros_like(H)
    dj, mj = [], []
    for kind, X in jumps:
        if kind == 'd':
            dj.append(X)
            Sgg += np.abs(X) ** 2
        else:
            mj.append(X)
            Smat += X.conj().T @ X
    Sd = np.diag(Sgg) + Smat

    def rhs(r):
        out = -1j * (H @ r - r @ H)
        for d in dj:
            out += np.outer(d, d.conj()) * r
        for L in mj:
            out += L @ r @ L.conj().T
        out -= 0.5 * (Sd @ r + r @ Sd)
        return out

    r = rho0.copy()
    for _ in range(nst):
        k1 = rhs(r)
        k2 = rhs(r + 0.5 * dt * k1)
        k3 = rhs(r + 0.5 * dt * k2)
        k4 = rhs(r + dt * k3)
        r = r + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return r / np.trace(r).real


def unbiased_starts(x):
    m = x.m
    z = tuple([0] * (4 * m))
    inn = tuple([1 if k < m else 0 for k in range(4 * m)])
    out = {}
    def kk(s):
        v = np.zeros(x.D, complex)
        v[x.i1[s]] = 1.0
        return v
    out["A"] = (kk(z) + kk(inn)) / np.sqrt(2)
    for lab, qi in (("B", 0), ("D", 1 % m)):
        t = [0] * (4 * m)
        for k in x.quad(qi):
            t[k] ^= 1
        t = tuple(t)
        u = list(t)
        for k in range(m):
            u[k] ^= 1
        out[lab] = (kk(t) + kk(tuple(u))) / np.sqrt(2)
    return out


print("  CONVERGENCE CHECK of the RK4 route against the m=2 exact superoperator exponential:")
for nst in (400, 800, 1600):
    jm = [('d', np.sqrt(0.5) * d) for d in md.zt] + [('m', np.sqrt(0.8) * md.emb_down(md.iota.conj().T))]
    r = rk4(H07, jm, rho_of(psi["A"]), 40.0, nst)
    ex = evolve_exact(H07, Ls_on, rho_of(psi["A"]), ts, D)[-1]
    print(f"      steps={nst:5d}  RK4 p+(40) = {obs(md,r)[0]:.9f}   exact = {obs(md,ex)[0]:.9f}"
          f"   diff = {abs(obs(md,r)[0]-obs(md,ex)[0]):.3e}")

print("\n  m-scaling (junction A, g^2=0.7, gamma_E=0.5, Gamma_d=0.8, T=40, RK4 1600 steps):")
print(f"  {'m':>3s} {'dim':>5s} {'init':>5s} {'p+(0)':>10s} {'p+(40)':>10s} {'<R_E>(40)':>11s} "
      f"{'pcap(40)':>10s} {'SPREAD(ABD)':>12s}")
print("  " + "-" * 76)
for m in (2, 3, 4):
    x = MD[m]
    st = unbiased_starts(x)
    Hx = x.H(0.7)
    jm = [('d', np.sqrt(0.5) * d) for d in x.zt] + [('m', np.sqrt(0.8) * x.emb_down(x.iota.conj().T))]
    Px = (np.eye(x.D) + np.diag(np.concatenate([x.RE1_all[0], x.RE0_all[0]]))) / 2.0
    Rx = np.diag(np.concatenate([x.RE1_all[0], x.RE0_all[0]])).astype(complex)
    Cx = np.zeros((x.D, x.D), complex)
    Cx[:x.D1, :x.D1] = np.eye(x.D1)
    vals = {}
    for nm in ("A", "B", "D"):
        r0 = rho_of(st[nm])
        r = rk4(Hx, jm, r0, 40.0, 1600)
        vals[nm] = np.trace(Px @ r).real
        print(f"  {m:3d} {x.D:5d} {nm:>5s} {np.trace(Px@r0).real:10.6f} {vals[nm]:10.6f} "
              f"{np.trace(Rx@r).real:+11.6f} {np.trace(Cx@r).real:10.3e} "
              f"{'':>12s}")
    print(f"  {'':>3s} {'':>5s} {'':>5s} {'':>10s} {'':>10s} {'':>11s} {'':>10s} "
          f"{max(vals.values())-min(vals.values()):12.6f}")

print("\n" + "=" * 100)
print("END OF RUN.  No verdicts printed above are hard-coded; every declared prediction is labelled")
print("as a prediction and every number is measured.")
print("=" * 100)
