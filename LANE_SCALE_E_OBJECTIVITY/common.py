"""LANE_SCALE_E_OBJECTIVITY -- shared machinery.

CARRIER FAMILY: [[n, n-2, 2]], n even, stabilisers X^(x)n and Z^(x)n, k = n-2 logical qubits.
    H = -(X^(x)n + Z^(x)n)   ->  eigenvalues -2 (mult 2^k), 0 (mult 2^(k+1)), +2 (mult 2^k)
    ground space = the (+1,+1) stabiliser sector, dim 2^k.

RECORDS: the Z_i of the conjugate pairs returned by symplectic_logicals.  NEVER NOMINATED.

THE EXACT FACTORISATION USED BY THE FAST PATH (and verified against record_model in s2):
  couple through the records themselves,
        H_tot = H (x) I + I (x) H_B + lam * sum_{i,j} W[i,j] R_i (x) X_j .
  Every system operator appearing (H, R_1..R_k) commutes with every other, so H_tot is block
  diagonal in the joint eigenbasis of (H, R_1..R_k).  The initial state P_g/2^k is exactly the
  uniform mixture over the 2^k joint sign sectors of the ground space (each 1-dimensional,
  self-checked), so

        rho(t) = sum_r 2^-k |r><r| (x) rho_B^r(t),       r in {+1,-1}^k

  a classical-quantum broadcast state.  Inside a sector the bath Hamiltonian is
        H_B + lam * sum_j ( sum_i W[i,j] r_i ) X_j = sum_j [ e_j Z_j + h_j(r) X_j ]
  which is a SUM OF SINGLE-SITE terms, so rho_B^r(t) = (x)_j rho_j^{h_j(r)}(t) exactly.

  This is not an approximation and not a shortcut around the model: s2 reproduces
  record_model.evolve + Environment.holevo to ~1e-12 at three sizes.  It is what makes
  k = 10 records with a 6-qubit bath reachable at all.
"""
import numpy as np
import sys, itertools

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import (RecordModel, Environment, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)

I2 = np.eye(2, dtype=complex)
XM = np.array([[0, 1], [1, 0]], dtype=complex)
ZM = np.array([[1, 0], [0, -1]], dtype=complex)

ENERGIES = (1.0, 1.4, 0.7, 1.1, 0.9, 1.25, 0.85, 1.55, 1.05, 0.75)
BETA = 2.0
TIMES = np.linspace(1.0, 13.0, 25)          # 25 times in [1,13]  -- time-average chi


# ------------------------------------------------------------------ the carrier
def code_stabilisers(n):
    """[[n, n-2, 2]] as (x|z) vectors over F_2^{2n}."""
    return [[1] * n + [0] * n, [0] * n + [1] * n]


def sp_form(a, b, n):
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2


def carrier(n):
    """Return dict with stabilisers, symplectic pairs, the k record (x|z) vectors, and checks."""
    stab = code_stabilisers(n)
    pairs = symplectic_logicals(stab, n)
    k = n - 2
    checks = {}
    checks['n_pairs'] = len(pairs)
    checks['n_pairs_ok'] = (len(pairs) == k)
    # SELF-CHECK 1: the symplectic pairing matrix must be the identity (non-degenerate).
    G = np.array([[sp_form(pairs[i][0], pairs[j][1], n) for j in range(len(pairs))]
                  for i in range(len(pairs))], dtype=int)
    checks['pairing_matrix'] = G
    checks['pairing_is_identity'] = bool(np.array_equal(G, np.eye(len(pairs), dtype=int)))
    # SELF-CHECK 2: the chosen records (the Z_i) must mutually commute.
    recs = [p[1] for p in pairs]
    ZZ = np.array([[sp_form(recs[i], recs[j], n) for j in range(len(recs))]
                   for i in range(len(recs))], dtype=int)
    checks['records_commute'] = bool(not ZZ.any())
    # SELF-CHECK 3: each record commutes with both stabilisers, and is not in <S>.
    checks['commute_with_stab'] = all(sp_form(r, s, n) == 0 for r in recs for s in stab)
    span = set()
    for a in (0, 1):
        for b in (0, 1):
            v = tuple((a * stab[0][t] + b * stab[1][t]) % 2 for t in range(2 * n))
            span.add(v)
    checks['not_a_stabiliser'] = all(tuple(r) not in span for r in recs)
    return dict(n=n, k=k, stab=stab, pairs=pairs, recs_xz=recs, checks=checks)


def record_matrices(car):
    n = car['n']
    return [xz_to_matrix(v, n) for v in car['recs_xz']]


def code_hamiltonian(n):
    Sx = xz_to_matrix(code_stabilisers(n)[0], n)
    Sz = xz_to_matrix(code_stabilisers(n)[1], n)
    return -(Sx + Sz)


# ------------------------------------------------------------------ coupling geometries
def weights(kind, k, nq, seed=7):
    """W[i,j] = weight with which record i couples to bath site j.

    EVERY GEOMETRY IS NORMALISED TO THE SAME TOTAL COUPLING PER SITE, sum_i W[i,j]^2 = 1, so
    that 'crowded' and 'separate' spend the same coupling budget at every site.  (D-16: the
    comparison is against the SPREAD control, never against an 'alone' value.)

      'crowded'   every record couples to every site, generic weights   -> a site can hold relations
      'sym'       every record couples to every site, EQUAL weights     -> degenerate venue (D-17)
      'separate'  site j couples to record (j mod k) ALONE              -> CONTROL: k independent
                  carriers with k disjoint baths; no relation can exist by construction
    """
    W = np.zeros((k, nq))
    if kind == 'crowded':
        rng = np.random.default_rng(seed)
        W = rng.normal(size=(k, nq))
    elif kind == 'sym':
        W = np.ones((k, nq))
    elif kind == 'separate':
        for j in range(nq):
            W[j % k, j] = 1.0
    else:
        raise ValueError(kind)
    W = W / np.sqrt((W ** 2).sum(axis=0, keepdims=True))
    return W


# ------------------------------------------------------------------ THE FAST PATH
def sign_patterns(k):
    """(2^k, k) array of +-1 sign patterns; bit i of the row index is record i."""
    idx = np.arange(2 ** k)
    bits = ((idx[:, None] >> np.arange(k)[None, :]) & 1)
    return 1 - 2 * bits                      # bit 0 -> +1, bit 1 -> -1


def site_states(nq, W, lam, times, energies=ENERGIES, beta=BETA):
    """rho_j^{h_j(r)}(t) for every sign pattern r, site j, time t.  Shape (nP, nq, nT, 2, 2)."""
    k = W.shape[0]
    S = sign_patterns(k).astype(float)                 # (nP, k)
    fields = lam * (S @ W)                             # (nP, nq)
    e = np.array(energies[:nq], dtype=float)
    nP = fields.shape[0]
    Hs = e[None, :, None, None] * ZM[None, None, :, :] + fields[:, :, None, None] * XM[None, None, :, :]
    w, U = np.linalg.eigh(Hs)                          # (nP,nq,2), (nP,nq,2,2)
    # thermal single-site state of e_j Z_j at inverse temperature beta
    tau = np.zeros((nq, 2, 2), dtype=complex)
    for j in range(nq):
        p = np.exp(-beta * np.array([e[j], -e[j]]))
        p = p / p.sum()
        tau[j] = np.diag(p).astype(complex)
    T = np.asarray(times, dtype=float)
    ph = np.exp(-1j * w[:, :, None, :] * T[None, None, :, None])     # (nP,nq,nT,2)
    Ue = U[:, :, None, :, :] * ph[:, :, :, None, :]                  # U diag(ph)
    Uc = np.conj(np.swapaxes(U, -1, -2))[:, :, None, :, :]
    Ut = Ue @ Uc                                                     # (nP,nq,nT,2,2)
    tt = tau[None, :, None, :, :]
    rho = Ut @ tt @ np.conj(np.swapaxes(Ut, -1, -2))
    return rho                                                        # (nP,nq,nT,2,2)


def kron_sites(rho, sites, ti):
    """Batched kron over the listed bath sites at time index ti.  rho: (nP,nq,nT,2,2)."""
    out = rho[:, sites[0], ti]                                       # (nP,2,2)
    for j in sites[1:]:
        a = out
        b = rho[:, j, ti]
        d = a.shape[-1]
        out = np.einsum('pij,pkl->pikjl', a, b).reshape(a.shape[0], 2 * d, 2 * d)
    return out


def entropies(mats):
    """von Neumann entropy in bits of a stack of density matrices."""
    ev = np.linalg.eigvalsh(mats)
    ev = np.where(ev > 1e-13, ev, 1.0)
    return -(ev * np.log2(ev)).sum(axis=-1)


def averages(rhoF, masks):
    """masks: (nMask, nP) row-stochastic weights.  Returns (nMask, d, d)."""
    nP, d, _ = rhoF.shape
    return (masks.astype(complex) @ rhoF.reshape(nP, d * d)).reshape(-1, d, d)


class Broadcast:
    """The classical-quantum broadcast state for one (carrier k, bath nq, W, lam), all times."""

    def __init__(self, k, nq, W, lam, times=TIMES, energies=ENERGIES, beta=BETA):
        self.k, self.nq, self.W, self.lam = k, nq, W, lam
        self.times = np.asarray(times, dtype=float)
        self.S = sign_patterns(k)                    # (nP,k)  +-1
        self.nP = self.S.shape[0]
        self.rho = site_states(nq, W, lam, self.times, energies, beta)

    # --- conditioning masks ------------------------------------------------
    def _mask(self, cond):
        """cond: dict {record_index: +-1} or callable(S)->bool array."""
        if callable(cond):
            sel = cond(self.S)
        else:
            sel = np.ones(self.nP, dtype=bool)
            for i, s in cond.items():
                sel &= (self.S[:, i] == s)
        m = sel.astype(float)
        return m / m.sum()

    def chi(self, sites, conds, ti):
        """Holevo of the observable whose outcomes are the given list of conditions.
           conds: list of (prior, cond).  Returns S(avg) - sum prior * S(cond state)."""
        rhoF = kron_sites(self.rho, sites, ti)
        masks = [self._mask(c) for _, c in conds]
        pri = np.array([p for p, _ in conds], dtype=float)
        allm = np.vstack([pri @ np.vstack(masks)] + masks)
        st = averages(rhoF, allm)
        S = entropies(st)
        return float(max(S[0] - float((pri * S[1:]).sum()), 0.0))

    # --- the standard readouts --------------------------------------------
    def chi_single(self, sites, i, ti):
        return self.chi(sites, [(0.5, {i: +1}), (0.5, {i: -1})], ti)

    def chi_pair(self, sites, i, j, ti):
        return self.chi(sites, [(0.25, {i: a, j: b}) for a in (+1, -1) for b in (+1, -1)], ti)

    def chi_parity(self, sites, i, j, ti):
        f = lambda p: (lambda S: S[:, i] * S[:, j] == p)
        return self.chi(sites, [(0.5, f(+1)), (0.5, f(-1))], ti)

    def chi_triple(self, sites, i, j, l, ti):
        return self.chi(sites, [(0.125, {i: a, j: b, l: c})
                                for a in (+1, -1) for b in (+1, -1) for c in (+1, -1)], ti)

    def chi_all(self, sites, ti):
        return self.chi(sites, [(2.0 ** -self.k, {i: int(s) for i, s in enumerate(row)})
                                for row in self.S], ti)


def time_average(fn, nT):
    return float(np.mean([fn(ti) for ti in range(nT)]))


def fmt(x, w=7, p=4):
    return f"{x:{w}.{p}f}"


# ------------------------------------------------------------------ batched readouts
def chi_batch(B, sites, ti, specs):
    """MANY Holevo quantities on ONE fragment at ONE time, with a single kron, a single masked
       average and a single batched eigvalsh.

       specs: list of (name, [(prior, cond), ...]).  For every observable used in this lane the
       priors average the conditional masks back to the UNIFORM mask, so S(avg) is computed once.
       That identity is asserted, not assumed."""
    rhoF = kron_sites(B.rho, sites, ti)
    uniform = np.full(B.nP, 1.0 / B.nP)
    rows = [uniform]
    layout = []
    for name, conds in specs:
        idx = []
        acc = np.zeros(B.nP)
        for pri, c in conds:
            m = B._mask(c)
            acc += pri * m
            idx.append(len(rows)); rows.append(m)
        assert np.max(np.abs(acc - uniform)) < 1e-12, f"{name}: priors do not average to uniform"
        layout.append((name, [p for p, _ in conds], idx))
    S = entropies(averages(rhoF, np.vstack(rows)))
    out = {'_Savg': float(S[0])}
    for name, pri, idx in layout:
        out[name] = float(max(S[0] - float(sum(p * S[q] for p, q in zip(pri, idx))), 0.0))
    return out


def thermal_entropies(nq, energies=ENERGIES, beta=BETA):
    """S(tau_j) per bath site, in bits.

    EVERY conditional bath state rho_B^r(t) is a UNITARY image of the thermal state -- the
    conditional dynamics inside a sector is unitary -- so its entropy is exactly S(tau) for
    every r and every t.  Hence the Holevo of the FULL k-record register is
        chi_ALL(F) = S(rho-bar_F) - sum_{j in F} S(tau_j),
    with no 2^k-fold averaging.  s3 checks this identity against the brute-force sum.
    """
    out = []
    for j in range(nq):
        e = energies[j]
        p = np.exp(-beta * np.array([e, -e])); p = p / p.sum()
        out.append(float(-(p * np.log2(p)).sum()))
    return np.array(out)


def subset_cond(B, subset, value):
    """condition on the PRODUCT of records in `subset` equalling `value` (+-1).
       subset = () is the identity element and is never used as an observable."""
    sub = list(subset)
    return lambda S: np.prod(S[:, sub], axis=1) == value


def spec_group_element(B, subset):
    return [(0.5, subset_cond(B, subset, +1)), (0.5, subset_cond(B, subset, -1))]
