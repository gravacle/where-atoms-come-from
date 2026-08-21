"""THE ARROW LAYER -- the family PROOF_V002 section 3 says is owed a layer, folded into
the model (T-54/T-55, the T-46 pattern).  The principal's directive, 2026-08-21: the URM is
the world model -- the framework new observations are added INTO.  This module is therefore
two things at once: the ported sealed machinery of the arrow results, and the ENTRY GATE
through which a NEW bath/fragment observation of this family's kind is scored
(score_bath_observation, D-25 provenance enforced).

Everything here is PORTED from the sealed lanes, with fidelity to sealed behavior over
elegance, or DELEGATED to record_model machinery the lanes' numbers were already re-derived
through (validate_formation.py).  Floats throughout -- these are eigh lanes, not F_2 lanes;
the one exact-integer stretch (the homology construction of the logicals) is ported intact.

WHAT EACH SECTION IMPLEMENTS, WITH ITS CLAIM ROW AND SEALED SOURCE:
  1  the toric 2x2 carrier, logicals computed     (shared venue of every arrow lane;
     from homology, NEVER nominated                LANE_F7_OCCUPANCY/f7_davies.py header,
     (W-62 convention)                             exec'd by both arrow lanes; integer-exact)
  2  HOLEVO INSTRUMENT chi(O:B) + I(S:B)          (shared; LANE_F1_ARROW/f1_arrow.py,
     for an arbitrary finite bath                  f1b_invariance.py.
                                                   record_model.Environment.holevo is the
                                                   qubit-bath twin of the same instrument)
  3  THE STATIC ARROW: mean-force state,   F-17   (LANE_F1_ARROW/f1_arrow.py part 4: chi
     the weight-threshold sweep                    0.00000000 over ALL 24 weight-1
                                                   observables, 0.11448276 at weight 2 = d;
                                                   closed form from Z_B(+-1) exact)
  4  ENTANGLEMENT WITHOUT RECORD BITS       F-18   (LANE_F1_ARROW/f1b_invariance.py (c):
     -- the ledger separating decoherence          weight-1 coupling I(S:B) = 0.04549256
     from the record's arrow                       with chi(Zbar:B) = 0.00000000)
  5  IRREVERSIBILITY FROM INSIDE            F-19   (LANE_F1_ARROW/f1b_invariance.py (a,b):
                                                   I(S:B) invariant to 3.686e-14 under 12
                                                   system-only unitaries, covariance
                                                   9.992e-16; chi about the FIXED label
                                                   moves by 1.145e-01 -- the copy is not
                                                   erasable from inside, only relocatable)
  6  THE HISTORY -- the arrow as formation  F-20-  (LANE_PF2_DYNAMICAL/pf2_history.py:
     from a product state, one eigen-      adjacent chi grows 0 -> 0.97527192 (t=1) while
     decomposition serving every time              <Zbar> is exactly constant; full
                                                   reversal control.  RecordModel.formation
                                                   is the sealed-equal single-time route,
                                                   validate_formation.py)
  7  FRAGMENT REDUNDANCY                    F-21   (LANE_PF2_DYNAMICAL parts 3-4: fragments
     -- WIRES RecordModel.redundancy,              hold 0.789366/0.048377/0.678602 bits
     which no validator called before T-55         under the weight-d coupling, exactly
                                                   zero under weight-1)
  8  OBSERVATION ENTRY: score a NEW               (T-54/T-55; the D-25 gate mirrored from
     bath/fragment observation through             project_model.URM.surface -- world-tier
     the family's own functions                    bath observations require provenance,
                                                   corner baths must self-declare DEF-A)

SCOPE, STATED PLAINLY (the sealed T-9 audit): F-17/F-18/F-19 are SINGLE-CARRIER
(toric-2x2); F-21's fragment bits are toric-only -- the T-9 battery replicated only the
whole-bath weight-1 null on [[8,1,2]] and [[4,2,2]], not fragment redundancy.  F-20 is
TWO-CARRIER (toric-2x2; bouquet).  Nothing here upgrades those scopes; the module
reproduces the sealed venue and gives new observations a scored way in.

BORROWED INSTRUMENTS, OWNERS NAMED: the Holevo quantity chi is Holevo (1973); the
fragment-redundancy apparatus is quantum-Darwinism-style, Zurek / Blume-Kohout-Zurek
(F-21 owner: ASSEMBLED -- the threshold-in-fragments finding is ours); the invariance of
I(S:B) under system-local unitaries is textbook quantum information (F-19 owner: BORROWED
-- the relative-arrow reading for records is ours, T-III.6); the toric carrier is Kitaev
(quant-ph/9707021); the mean-force Gibbs state is standard open-system equilibrium.
OURS (per the register rows): the threshold-at-d law F-17 (owner: ORIGINAL), the measured
decoherence/record separation F-18 (owner: ORIGINAL), the read-not-written copy mechanism
checks F-20 (owner: ASSEMBLED).

PACKAGING NOTES (the two adaptations beyond re-organisation, both identical maps, never
approximations): (1) the lanes apply the conditioning projector as an explicit
kron(P, I_bath) matrix product at joint dimension; here the SAME projector is applied
blockwise on the reshaped state, so the 252-observable weight-2 sweep stays inside the
check budget.  (2) arrow_history diagonalises the venue's joint Hamiltonian on the real
path when its imaginary part is EXACTLY zero (it is, on this venue) -- ~7x cheaper at
dim 2048; the check block gates that route against the complex path inside
RecordModel.redundancy.  Every gated number still equals its sealed print.  The one
number these economies visibly move is a machine-epsilon-scale instrument bound: the
covariance worst case reads 3.9e-15 here vs the lane's 9.992e-16 -- both are the same
statement (zero at working precision) and the gate is the lane's own 1e-8 bound.

The sealed lanes remain the source of truth; checks_arrow.py gates every number this
module reproduces against its SEALED value.  Returns are DATA."""
import os as _os
import sys as _sys
import itertools

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

import numpy as np


# =====================================================================================
# 1  THE CARRIER (LANE_F7_OCCUPANCY/f7_davies.py header, verbatim-in-substance; the
#    construction both arrow lanes exec'd.  Logicals computed from homology -- the W-62
#    convention: NEVER nominated.  Integer-exact F_2 bitmask arithmetic, as the lane.)
# =====================================================================================
I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def op(p, L):
    """Tensor product over L qubits with single-site factors given in dict p."""
    M = np.array([[1]], dtype=complex)
    for l in range(L):
        M = np.kron(M, p.get(l, I2))
    return M


def _rref_bits(rows):
    piv = {}
    for r in rows:
        c = r
        for j in sorted(piv, reverse=True):
            if (c >> j) & 1:
                c ^= piv[j]
        if c:
            piv[c.bit_length() - 1] = c
    return piv


def _inspan(v, piv):
    for j in sorted(piv, reverse=True):
        if (v >> j) & 1:
            v ^= piv[j]
    return v == 0


def _nullspace(M):
    M = M.copy() % 2
    rows, cols = M.shape
    pc = []
    r = 0
    for c in range(cols):
        pv = next((i for i in range(r, rows) if M[i, c]), None)
        if pv is None:
            continue
        M[[r, pv]] = M[[pv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        pc.append(c)
        r += 1
    out = []
    for fc in [c for c in range(cols) if c not in pc]:
        v = np.zeros(cols, dtype=np.int8)
        v[fc] = 1
        for i, pcc in enumerate(pc):
            v[pcc] = M[i, fc]
        out.append(v)
    return out


def _toint(v):
    return int(''.join(map(str, v)), 2)


_CARRIER = None


def carrier():
    """The sealed venue: toric code on the 2x2 torus, L=8 qubits, dim 256, ground space 4.
       Logicals Zbar, Zbar2 (weight 2 = d), Xbar computed from the chain complex; Ze is the
       sealed single-site (weight-1) coupling op({ind[('h',0,0)]: Z}).  Built once, cached.
       Sealed source: LANE_F7_OCCUPANCY/f7_davies.py header (exec'd by LANE_F1_ARROW and
       LANE_PF2_DYNAMICAL); owners: carrier Kitaev, construction the lane's."""
    global _CARRIER
    if _CARRIER is not None:
        return _CARRIER
    nx = ny = 2
    L = 2 * nx * ny
    ind = {}
    k = 0
    for j in range(ny):
        for i in range(nx):
            ind[('h', i, j)] = k; k += 1
            ind[('v', i, j)] = k; k += 1
    STAR = [[ind[('h', i, j)], ind[('h', (i - 1) % nx, j)], ind[('v', i, j)],
             ind[('v', i, (j - 1) % ny)]] for j in range(ny) for i in range(nx)]
    PLAQ = [[ind[('h', i, j)], ind[('v', (i + 1) % nx, j)], ind[('h', i, (j + 1) % ny)],
             ind[('v', i, j)]] for j in range(ny) for i in range(nx)]
    H0 = -sum(op({l: X for l in st}, L) for st in STAR) \
         - sum(op({l: Z for l in p}, L) for p in PLAQ)
    # ---- logicals computed from structure, not nominated (W-62 convention) ----
    EDGES = [None] * L
    for j in range(ny):
        for i in range(nx):
            EDGES[ind[('h', i, j)]] = (j * nx + i, j * nx + (i + 1) % nx)
            EDGES[ind[('v', i, j)]] = (j * nx + i, ((j + 1) % ny) * nx + i)
    NV = nx * ny
    d1 = np.zeros((NV, L), dtype=np.int8)
    for k_, (a, b) in enumerate(EDGES):
        d1[a, k_] ^= 1
        d1[b, k_] ^= 1
    d2 = np.zeros((L, len(PLAQ)), dtype=np.int8)
    for k_, pl in enumerate(PLAQ):
        for e in pl:
            d2[e, k_] ^= 1
    assert not ((d1 @ d2) % 2).any(), "d1 d2 != 0"
    Z1 = [int(_toint(v)) for v in _nullspace(d1)]
    B1 = _rref_bits([int(_toint(d2[:, c])) for c in range(d2.shape[1])])
    Zcand = [c for c in
             [x for m in range(1, 1 << len(Z1))
              for x in [int(np.bitwise_xor.reduce([Z1[i] for i in range(len(Z1))
                                                   if (m >> i) & 1]))]]
             if not _inspan(c, B1)]
    Zc1 = min(Zcand, key=lambda c: bin(c).count('1'))
    Zsp = _rref_bits([Zc1] + [int(_toint(d2[:, c])) for c in range(d2.shape[1])])
    Zc2 = min([c for c in Zcand if not _inspan(c, Zsp)], key=lambda c: bin(c).count('1'))
    Zdual = [int(_toint(v)) for v in _nullspace(d2.T)]
    Bd = _rref_bits([int(_toint(d1[r, :])) for r in range(d1.shape[0])])
    Xcand = [c for c in
             [x for m in range(1, 1 << len(Zdual))
              for x in [int(np.bitwise_xor.reduce([Zdual[i] for i in range(len(Zdual))
                                                   if (m >> i) & 1]))]]
             if not _inspan(c, Bd)]

    def ov(a, b):
        return bin(a & b).count('1') % 2

    Xc1 = min([c for c in Xcand if ov(c, Zc1) == 1 and ov(c, Zc2) == 0],
              key=lambda c: bin(c).count('1'))

    def zopc(c):
        return op({l: Z for l in range(L) if (c >> (L - 1 - l)) & 1}, L)

    def xopc(c):
        return op({l: X for l in range(L) if (c >> (L - 1 - l)) & 1}, L)

    Zbar, Zbar2, Xbar = zopc(Zc1), zopc(Zc2), xopc(Xc1)
    E0 = np.linalg.eigvalsh(H0)
    gs = int(np.sum(np.abs(E0 - E0[0]) < 1e-9))
    _CARRIER = dict(
        L=L, nS=2 ** L, ind=ind, STAR=STAR, PLAQ=PLAQ, H0=H0,
        Zbar=Zbar, Zbar2=Zbar2, Xbar=Xbar,
        Ze=op({ind[('h', 0, 0)]: Z}, L),
        weights=(bin(Zc1).count('1'), bin(Zc2).count('1'), bin(Xc1).count('1')),
        dimH1=len(Z1) - len(B1), gs=gs)
    return _CARRIER


# =====================================================================================
# 2  THE HOLEVO INSTRUMENT (LANE_F1_ARROW/f1_arrow.py chi(), f1b_invariance.py mutual(),
#    verbatim-in-substance; blockwise projector application per the header's packaging
#    note).  record_model.Environment.holevo is the same instrument specialised to the
#    qubit bath -- checks_arrow.py gates the two routes against each other.
# =====================================================================================
def vN(r):
    """von Neumann entropy in bits (the lanes' own clamp at 1e-13)."""
    e = np.linalg.eigvalsh(r)
    e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())


def _conditional_bath(r, O, nS, nB, s):
    """(p_s, bath state conditioned on O = s): the lane's kron(P,I) r kron(P,I) followed
       by the partial trace over the system, applied blockwise (identical linear map)."""
    P = (np.eye(nS) + s * O) / 2
    rr = r.reshape(nS, nB, nS, nB)
    t1 = np.tensordot(P, rr, axes=([1], [0]))          # [a, i, d, j]
    blk = np.moveaxis(np.tensordot(t1, P, axes=([2], [0])), 3, 2)   # [a, i, c, j]
    p = float(np.real(np.einsum('aiai->', blk)))
    rB = np.einsum('aiaj->ij', blk)
    return p, rB


def chi(r, O, nS, nB):
    """Holevo information chi(O:B) the bath holds about the +-1 observable O (O^2 = I).
       Claim rows F-17/F-18; sealed source LANE_F1_ARROW/f1_arrow.py chi() (its convention:
       an empty branch contributes (0, 0-matrix)).  Owner: Holevo (1973); the use is the
       lane's."""
    out = []
    for s in (+1, -1):
        p, rB = _conditional_bath(r, O, nS, nB, s)
        if p < 1e-12:
            out.append((0.0, np.zeros((nB, nB), dtype=complex)))
            continue
        out.append((p, rB / p))
    av = sum(p * rB for p, rB in out)
    return max(vN(av) - sum(p * vN(rB) for p, rB in out), 0.0)


def mutual(r, nS, nB):
    """I(S:B) of a joint state -- the F-19 invariant.  Sealed source
       LANE_F1_ARROW/f1b_invariance.py mutual(); owner: textbook quantum information."""
    rS = r.reshape(nS, nB, nS, nB).trace(axis1=1, axis2=3)
    rB = r.reshape(nS, nB, nS, nB).trace(axis1=0, axis2=2)
    return vN(rS) + vN(rB) - vN(r)


# =====================================================================================
# 3-5  THE STATIC ARROW (LANE_F1_ARROW): the 4-level bath, the mean-force Gibbs state,
#      the threshold sweep (F-17), the ledger (F-18), the invariance theorem (F-19).
# =====================================================================================
EB = np.array([0.0, 0.7, 1.3, 2.1])      # the sealed 4-level bath energies
BB = np.array([1.0, 0.3, -0.2, -0.9])    # the sealed bath-coupling eigenvalues
NB4 = 4
BETA = 2.0
HB4 = np.diag(EB)
BOP4 = np.diag(BB)


def mean_force_state(A, lam=0.8, beta=BETA):
    """Joint Gibbs state of H_tot = H0 x I + I x HB + lam A x Bop at the sealed beta = 2.
       Sealed source: LANE_F1_ARROW rho_SB(), verbatim.  Owner: mean-force Gibbs is
       standard open-system equilibrium."""
    C = carrier()
    nS = C['nS']
    Ht = np.kron(C['H0'], np.eye(NB4)) + np.kron(np.eye(nS), HB4) + lam * np.kron(A, BOP4)
    w, U = np.linalg.eigh(Ht)
    w = w - w.min()
    M = (U * np.exp(-beta * w)) @ U.conj().T
    return M / np.trace(M)


def closed_form_chi(lam, beta=BETA):
    """The INDEPENDENT instrument for chi(Zbar:B) under the Zbar coupling: conditional
       bath partition functions Z_B(+-1) -- no joint diagonalisation, no chi() call.
       Sealed source: LANE_F1_ARROW/f1_arrow.py 'closed form' column (exact match row
       F-17).  A computed comparison in the D-8 sense: the check block gates chi against
       THIS, not against a transcribed number alone."""
    pB = np.exp(-beta * (EB + lam * BB))
    pM = np.exp(-beta * (EB - lam * BB))
    Zp, Zm = pB.sum(), pM.sum()
    pp, pm = Zp / (Zp + Zm), Zm / (Zp + Zm)
    av = pp * (pB / Zp) + pm * (pM / Zm)

    def H(v):
        v = v[v > 1e-13]
        return float(-(v * np.log2(v)).sum())

    return H(av) - pp * H(pB / Zp) - pm * H(pM / Zm)


def weight_sweep(r, weight):
    """F-17's deciding instrument: sweep ALL observables of the given weight (every site
       combination x every Pauli letter) in the SAME state, return the count, the maximum
       chi, and where it was attained.  Sealed source: LANE_F1_ARROW/f1_arrow.py part 4
       (24 weight-1 observables all at 0.00000000; 252 weight-2, max 0.11448276)."""
    C = carrier()
    L, nS = C['L'], C['nS']
    P1 = {'X': X, 'Y': 1j * (X @ Z), 'Z': Z}
    mx = 0.0
    arg = None
    n = 0
    for sites in itertools.combinations(range(L), weight):
        for lets in itertools.product('XYZ', repeat=weight):
            O = op({s: P1[c] for s, c in zip(sites, lets)}, L)
            n += 1
            v = chi(r, O, nS, NB4)
            if v > mx:
                mx = v
                arg = ''.join(lets) + str(list(sites))
    return dict(weight=weight, n_swept=n, max_chi=mx, argmax=arg)


def arrow_threshold(lam=0.8, weights=(1, 2), coupling=None):
    """F-17: THE ARROW CARRIES THE RECORD'S OWN THRESHOLD.  In the mean-force state under
       the weight-d coupling (default Zbar), chi(O:B) over every observable of each listed
       weight, plus the logical row.  Sealed: 0.00000000 at weight 1 (all 24), 0.11448276
       at weight 2 = d, closed form exact.  Owner: ORIGINAL (T-36).  Scope: SINGLE-CARRIER
       (toric-2x2, sealed T-9 audit)."""
    C = carrier()
    A = C['Zbar'] if coupling is None else coupling
    r = mean_force_state(A, lam=lam)
    out = {w: weight_sweep(r, w) for w in weights}
    out['logical'] = dict(chi_Zbar=chi(r, C['Zbar'], C['nS'], NB4),
                          chi_Zbar2=chi(r, C['Zbar2'], C['nS'], NB4),
                          closed_form=closed_form_chi(lam))
    return out


def arrow_ledger(lam=0.8):
    """F-18: THE RECORD'S ARROW IS NOT AMBIENT DECOHERENCE.  The sealed four-row ledger --
       coupling, weight, I(S:B), chi(Zbar:B).  The weight-1 row entangles (I = 0.04549256)
       and transfers ZERO record bits.  Sealed source: LANE_F1_ARROW/f1b_invariance.py (c).
       Owner: ORIGINAL (T-36).  Scope: SINGLE-CARRIER."""
    C = carrier()
    nS = C['nS']
    rows = []
    for nm, A, wt in (("Zbar", C['Zbar'], 2), ("Zbar2", C['Zbar2'], 2),
                      ("Ze", C['Ze'], 1), ("identity", np.eye(nS), 0)):
        r = mean_force_state(A, lam=lam)
        rows.append(dict(coupling=nm, weight=wt,
                         I_SB=mutual(r, nS, NB4),
                         chi_record=chi(r, C['Zbar'], nS, NB4)))
    return rows


def arrow_invariance(n_unitaries=12, seed=5, lam=0.8):
    """F-19: IRREVERSIBILITY FROM INSIDE.  Random system-only unitaries on the mean-force
       state: covariance of the instrument (must vanish), invariance of I(S:B) (the
       theorem; sealed 3.686e-14 over 12 unitaries at seed 5), movement of chi about the
       FIXED label (may move -- and does, 1.145e-01: the copy is relocatable, never
       erasable, from inside).  Sealed source: LANE_F1_ARROW/f1b_invariance.py (a,b).
       Owner: BORROWED (textbook invariance); the relative-arrow reading ours (T-III.6)."""
    C = carrier()
    nS = C['nS']
    r = mean_force_state(C['Zbar'], lam=lam)
    c0 = chi(r, C['Zbar'], nS, NB4)
    I0 = mutual(r, nS, NB4)
    rng = np.random.default_rng(seed)
    wc = wi = wfixed = 0.0
    for _ in range(n_unitaries):
        M = rng.normal(size=(nS, nS)) + 1j * rng.normal(size=(nS, nS))
        Q, _ = np.linalg.qr(M)
        Us = np.kron(Q, np.eye(NB4))
        rp = Us @ r @ Us.conj().T
        wc = max(wc, abs(chi(rp, Q @ C['Zbar'] @ Q.conj().T, nS, NB4) - c0))
        wi = max(wi, abs(mutual(rp, nS, NB4) - I0))
        wfixed = max(wfixed, abs(chi(rp, C['Zbar'], nS, NB4) - c0))
    return dict(chi0=c0, I0=I0, covariance_worst=wc, mutual_worst=wi,
                fixed_label_worst=wfixed, n_unitaries=n_unitaries, seed=seed)


# =====================================================================================
# 6-7  THE DYNAMICAL ARROW (LANE_PF2_DYNAMICAL): the history from a product state, and
#      fragment redundancy -- DELEGATING to record_model machinery wherever it exists
#      (Environment, holevo, RecordModel.redundancy: the wiring T-54 called for).
# =====================================================================================
_MODEL = None
_ENV = None


def pf2_env():
    """The sealed PF-2 bath IS record_model.Environment's default: 3 qubits, energies
       (1.0, 1.4, 0.7), beta 2.0, probe = sum_j X_j.  Built once, cached."""
    global _ENV
    if _ENV is None:
        from record_model import Environment
        _ENV = Environment()
    return _ENV


def record_model_instance():
    """RecordModel on the carrier with no noise -- the object whose .formation and
       .redundancy the T-7 validator proved sealed-equal.  Built once, cached."""
    global _MODEL
    if _MODEL is None:
        from record_model import RecordModel
        _MODEL = RecordModel(carrier()['H0'], [])
    return _MODEL


def arrow_history(times, coupling=None, lam=0.8, env=None, keep_states=()):
    """F-20-adjacent: THE ARROW AS A HISTORY.  From the product state (maximally mixed
       code space) x (thermal bath), one eigendecomposition serves every requested time --
       the lane's own economy (LANE_PF2_DYNAMICAL/pf2_history.py run()).  Per time:
       <Zbar> (exactly conserved under the Zbar coupling -- the record is READ, not
       written) and chi(Zbar:B) through record_model.Environment.holevo.  Negative times
       are legal and are the sealed reversal control.  keep_states: times whose joint
       state to return (for fragment or I(S:B) readouts without another eigh).
       Sealed anchors: chi = 0 at t=0; 0.40660635 / 0.81447230 / 0.97527192 / 0.78665760 /
       0.90811968 at t = 0.25/0.5/1/2/4.  RecordModel.formation is the single-time twin
       (validate_formation.py gates it).  Owner: ASSEMBLED; scope TWO-CARRIER for the
       mechanism (F-20: toric-2x2; bouquet), this venue toric."""
    C = carrier()
    env = pf2_env() if env is None else env
    m = record_model_instance()
    A = C['Zbar'] if coupling is None else coupling
    nS, nB = C['nS'], env.dim
    Pg, k = m.ground_space()
    r0 = np.kron(Pg / k, env.thermal())
    Ht = np.kron(C['H0'], np.eye(nB)) + np.kron(np.eye(nS), env.HB) \
        + lam * np.kron(A, env.probe)
    # Packaging economy, second header note: on this venue Ht is EXACTLY real symmetric
    # (X/Z-type operators, real bath), and the real-path eigh is ~7x cheaper at dim 2048.
    # Same operator, same spectrum; the check block gates this route against the complex
    # path inside RecordModel.redundancy at the sealed print precision.
    if np.abs(Ht.imag).max() == 0.0:
        Ht = Ht.real
    w, U = np.linalg.eigh(Ht)
    Uc = U.conj().T @ r0 @ U
    out = []
    states = {}
    for t in times:
        ph = np.exp(-1j * w * t)
        r = U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
        rS = r.reshape(nS, nB, nS, nB).trace(axis1=1, axis2=3)
        out.append(dict(t=t,
                        value=float(np.real(np.trace(rS @ C['Zbar']))),
                        chi=env.holevo(r, C['Zbar'], nS)))
        if t in keep_states:
            states[t] = r
    return out, states


def arrow_redundancy(coupling=None, lam=0.8, t=4.0, env=None):
    """F-21: REDUNDANCY CARRIES THE RECORD'S THRESHOLD.  Delegates to
       RecordModel.redundancy -- whole-bath chi and every single-qubit fragment from ONE
       evolution -- which no validator called before T-55; checks_arrow.py now gates it.
       Sealed anchors (weight-d coupling, t=4): whole 0.90811968, fragments
       0.789366 / 0.048377 / 0.678602; weight-1 coupling: all EXACTLY ZERO.
       Owner: ASSEMBLED (quantum-Darwinism-style apparatus; threshold-in-fragments ours).
       Scope: fragment bits SINGLE-CARRIER toric-only; the T-9 battery replicated only the
       whole-bath weight-1 null on [[8,1,2]]/[[4,2,2]] (sealed audit row F-21)."""
    C = carrier()
    env = pf2_env() if env is None else env
    m = record_model_instance()
    A = C['Zbar'] if coupling is None else coupling
    whole, parts = m.redundancy(C['Zbar'], A, env, lam=lam, t=t)
    return dict(whole=whole, fragments=[float(p) for p in parts])


# =====================================================================================
# 8  OBSERVATION ENTRY -- how a NEW bath/fragment observation enters the URM through
#    this layer (T-54/T-55; the principal 2026-08-21: the framework observations are
#    added INTO).  D-25 mirrored from project_model.URM.surface.
# =====================================================================================
def score_bath_observation(env, coupling, record=None, model=None, lam=0.8, t=4.0,
                           tier="world", provenance=None, tol=1e-9):
    """Score a NEW bath/fragment observation through the arrow family's own instruments.

       env: a record_model.Environment (the bath spec: qubit count, energies, beta) --
       the observed environment.  coupling: the system operator it couples through (or
       the distributed/full forms RecordModel.formation accepts).  record: the +-1 record
       observable (default: the carrier's Zbar only when the default toric model is
       used).  model: a RecordModel (default: the sealed toric venue).  A caller that
       supplies a custom model must also supply that model's record explicitly; the
       toric record is never silently paired with a different Hilbert space.

       D-25 AT THE GATE: a world-tier observation REFUSES to enter without provenance --
       the real bath it models and its constants' pinned sources; a corner bath must
       self-declare provenance='DEF-A'.  The exact idealisation may never silently pose
       as the world.

       Returns the scored entry: I(S:B) (does it entangle at all), chi_whole (does it
       hold record bits -- the F-17/F-18 discriminator), per-fragment chi and the
       redundancy count (F-21's instrument), and the verdicts.  every outcome registers:
       entangled_without_record = True is a RESULT (the F-18 class), not a failure."""
    if tier == "corner":
        if provenance != "DEF-A":
            raise ValueError(
                "ARROW GATE REFUSES: a corner bath must self-declare provenance='DEF-A' "
                "-- the exact idealisation may never silently pose as the world (D-25).")
    else:
        if not provenance or not str(provenance).strip():
            raise ValueError(
                "ARROW GATE REFUSES: a world-tier bath observation requires PROVENANCE -- "
                "the real environment it models and its constants' pinned sources (D-25, "
                "the principal 2026-08-20: the model is grounded in real record data, "
                "never the toy category).")
    if model is None:
        C = carrier()
        m = record_model_instance()
        R = C['Zbar'] if record is None else np.asarray(record)
    else:
        m = model
        if record is None:
            raise ValueError(
                "ARROW GATE REFUSES: a custom RecordModel requires its record observable "
                "explicitly; the toric default record belongs only to the default model.")
        R = np.asarray(record)
    if R.shape != (m.n, m.n):
        raise ValueError(
            f"ARROW GATE REFUSES: record shape {R.shape} does not match model dimension "
            f"{m.n}x{m.n}.")
    r = m.evolve(coupling, env, lam=lam, t=t)
    nS = m.n
    chi_whole = env.holevo(r, R, nS)
    frags = [env.holevo(r, R, nS, fragment=[j]) for j in range(env.nq)]
    I_SB = mutual(r, nS, env.dim)
    return dict(tier=tier, provenance=provenance, t=t, lam=lam,
                I_SB=I_SB, chi_whole=chi_whole, fragments=frags,
                holds_record_bits=bool(chi_whole > tol),
                entangled_without_record=bool(I_SB > tol and chi_whole <= tol),
                redundant_fragments=int(sum(f > tol for f in frags)))


if __name__ == "__main__":
    # Standalone smoke run: the carrier and the two headline instruments, briefly.
    C = carrier()
    print(f"carrier: dim {C['nS']}, ground degeneracy {C['gs']}, dim H_1 = {C['dimH1']}, "
          f"logical weights {C['weights']}")
    r = mean_force_state(C['Zbar'])
    print(f"chi(Zbar:B) mean-force = {chi(r, C['Zbar'], C['nS'], NB4):.8f}   "
          f"closed form = {closed_form_chi(0.8):.8f}")
    print("full gates: python3 checks_arrow.py")
