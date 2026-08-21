"""THE COUNT-LAW LAYER -- the C-86 census machinery, folded into the model (T-54, family countlaw).

THE CLAIM ROW THIS MODULE HOMES: C-86 (grade FORMAL) -- THE SURVIVING-RECORD COUNT LAW.
One criterion -- clause (ii') |lambda_record| <= 1/t_m on the record's OWN Liouvillian mode --
yields BOTH widths with zero adjustable content:

    delta_pop(t_m) = kT ln( expm1( B/kT - ln(f0 t_m) ) )      (population face, diagonal record)
                   = [B - kT ln(f0 t_m)] - kT ln(1 + e^{-dE*/kT})   (naive width, derived correction)
    delta_coh(t_m) = hbar / t_m                               (coherence face, off-diagonal record)

and the dated staircase

    k(t_m) = #{ i : dE_i <= delta_pop,i(t_m) } = #{ i : g_u(i) + g_l(i) <= 1/t_m },
    t*_i   = f0^{-1} exp((B_i - dE_i)/kT) / (1 + e^{-dE_i/kT})   (parameter-free drop times;
             each record dies at its SHALLOWER value's escape).

Everything here is PORTED from the sealed lanes, with fidelity to sealed behavior over elegance:
  1  THE TWO WIDTHS + the one-modulus unification   C-86   (LANE_T47_A_WIDTH/t47_a_width.py)
  2  THE RECORD'S OWN MODE on a RecordSurface       C-86   (LANE_T47_A_WIDTH rate_num;
     (instrument path; closed form is the check)            grounded.clause_ii, the sealed instrument)
  3  THE CENSUS census(surfaces, t_m)               C-86   (LANE_T47_A_WIDTH section 6,
     -- the ProjectModel-ready layer method                 LANE_T47_B_STAIRCASE sections C/D/G)
  4  THE T-31 MULTI-SHELL CONTROL CARRIER           C-86 / (LANE_T47_B_STAIRCASE section F --
     [[4,2,2]]-class, Davies bath, basin-lumped     C-76    the kill becomes a staircase; register
     record modes                                           note: clause_ii misreads multi-shell
                                                            carriers, the basin instrument is the
                                                            record's own mode there)

THE CORNER IS REFERENCED, NOT DUPLICATED: the exact corner law k = min_E v2(m_E) (C-14, PARTIAL --
recovered as the zero-splitting corner of C-86) is homed at model/count_law.py. Nothing in this
module offers a spectral-multiplicity count; the private _v2/_exact_multiplicities helpers exist
ONLY so the T-31 control can reproduce the corner instrument's reading beside the record-mode
census (C-76's kill exhibited, never adopted).

THE C-76 GATE, STRUCTURAL: C-76 killed the first k(t_m) because its clustering width was CHOSEN
(an unowned 1e-2 while the actual derived width was ~1e-39x smaller). The successor width is the
carrier's own, and this module's API carries the lesson in its shape: census(surfaces, t_m) takes
NOTHING else -- no width, no tolerance, no margin, no clustering knob -- and no public function
here accepts one. checks_countlaw.py gates that unreachability.

CONVENTION (declared, the T-47 lanes' own): B is the barrier measured from the LOWER well;
project_model.RecordSurface's E_b is the standard Arrhenius activation from the METASTABLE (upper)
well, so B = E_b + dE and B - dE IS the model's E_b. The map is exact and is used, not re-derived.
H-UNITS (register note, T-47): grounded.liouvillian takes H in generator (angular-frequency)
units; RecordSurface.open_system returns H in joules. Harmless on the population face ([H, R] = 0,
the Hamiltonian part drops out of the record's mode -- the sealed numbers came down this exact
path); coherence-face callers must pass H/hbar, as the sealed lanes did.

DECLINE, NOT A NUMBER: a surface that is not thermally activated, or whose Arrhenius rates
underflow floats, DECLINES through open_system() -- the census reports it in `declined`, never
silently counts or drops it. (LANE_T47_A's sweep shim returned rate 0.0 there; the model's
registered convention is to decline, and the census follows the model.)

ASSUMPTIONS (declared, the lanes'): shared Arrhenius prefactor f0 for both wells (equal Kramers
attempt frequencies; detailed balance is then a consequence); Markovian GKSL thermal activation;
f0 t_m >> 1. Scope: on multi-pathway/multi-shell carriers the one-f0 closed form leaves
ln-residuals and the staircase is carried by the record-mode instrument (C-86 scope clause; the
T-31 section here is that instrument).

OWNERS, NAMED WHERE THE SEALED LANES NAMED THEM (ownership PARTIAL per C-86): Neel 1949 /
Street-Woolley 1949 / Sharrock 1994 own the activation window kT ln(f0 t) as a remanence-DECAY
device; Charap-Lu-He 1997 and Weller-Moser 1999 own the dE = 0 corner criterion; field-tilted
Neel-Brown and Korman-Mayergoyz Preisach-Arrhenius own (B - dE)-type barriers at the RATE level;
textbook two-state kinetics owns the (1 + e^{-dE/kT}) both-rates factor; Alicki-Fannes-Horodecki /
Bravyi-Terhal own coherence-face storage-time bounds; Davies owns the weak-coupling generator used
by the T-31 bath. WHOLLY PROGRAM-OWNED: the width DERIVED from the record definition, the
diagonal/coherence unification under one clause, the margin-free integer CENSUS as a law, and the
departure term sum_over_dead tanh(dE_i/2kT).

The sealed lanes remain the source of truth; checks_countlaw.py gates every number this module
reproduces against its SEALED value, plus API-fidelity probes beyond the gated range.
Returns are DATA."""
import sys as _sys
import os as _os

_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in _sys.path:
    _sys.path.insert(0, _HERE)

import numpy as _np
import grounded as G

EV = 1.602176634e-19       # J, exact SI (the lanes' own constant)
YEAR = 3.156e7             # s, the house convention (T-29, T-31, T-47)


# =====================================================================================
# 1  THE TWO WIDTHS  (C-86; LANE_T47_A_WIDTH/t47_a_width.py delta_exact + section 4/5)
# =====================================================================================
def delta_pop(B, T, f0, t_m):
    """C-86 population-face width, the EXACT derived form (LANE_T47_A_WIDTH, sealed):
       delta_pop(t_m) = kT ln(expm1(B/kT - ln(f0 t_m))), overflow-safe. B is the barrier
       above the LOWER well, in J; T in K; f0 in Hz; t_m in s. Returns None when no
       dE >= 0 crossing exists -- expm1(y) <= 1, i.e. t_m beyond the symmetric bound
       exp(B/kT)/(2 f0); the no-crossing condition IS the symmetric bound (C-86 corner
       limit iii). Zero adjustable content: every symbol is the carrier's own or SI."""
    kT = G.KB * T
    y = B / kT - _np.log(f0 * t_m)
    if y <= _np.log(2.0):
        return None                      # dE* <= 0: no record with this B survives t_m
    with _np.errstate(over='ignore'):    # huge-barrier surfaces: expm1 -> inf, delta -> inf
        return float(kT * _np.log(_np.expm1(y)))


def delta_coh(t_m):
    """C-86 coherence-face width: delta_coh(t_m) = hbar/t_m, the SAME clause-(ii')
       criterion on the off-diagonal record mode (LANE_T47_A_WIDTH section 4, sealed)."""
    return G.HBAR / t_m


def coh_modulus(g_total, dE):
    """C-86's one-expression unification (LANE_T47_A_WIDTH section 5, sealed): with
       dissipation on, the coherence eigenvalue modulus is
           |lambda| = sqrt( ((g_u+g_l)/2)^2 + (dE/hbar)^2 )
       whose two corners are the two widths -- dissipation-dominated -> the Arrhenius
       face, rotation-dominated -> hbar/t_m. g_total = g_u + g_l in 1/s, dE in J."""
    return float(_np.sqrt((g_total / 2.0) ** 2 + (dE / G.HBAR) ** 2))


# =====================================================================================
# 2  THE RECORD'S OWN MODE ON A RecordSurface  (the instrument path; LANE_T47_A rate_num)
# =====================================================================================
def record_rate(s):
    """The record's own Liouvillian-mode rate FROM THE SEALED INSTRUMENT
       (grounded.clause_ii's Rayleigh quotient -- provably exact on the two-well carrier,
       register note T-47), never from the closed form. Returns 1/s, or None where the
       surface DECLINES (not thermally activated, or rates below float range).
       s: a project_model.RecordSurface (E_b in J from the metastable well, per the
       solidity-review convention)."""
    _require_urm_surface(s)
    opened = s.open_system()
    if opened is None:
        return None
    H, Ls, R = opened
    return G.clause_ii(H, Ls, R, _np.inf)['rate']


def drop_time(s):
    """t*_i by the INSTRUMENT: 1 / (the record's own mode rate). None where the model
       declines. This is the measurement path; drop_time_formula is the check."""
    r = record_rate(s)
    return None if (r is None or r <= 0.0) else 1.0 / r


def drop_time_formula(s):
    """C-86's parameter-free drop time, closed form (LANE_T47_A_WIDTH section 6, sealed):
           t*_i = f0^{-1} exp((B_i - dE_i)/kT) / (1 + e^{-dE_i/kT})
       with B - dE = E_b (the model's own convention map, exact). THE CHECK, NEVER THE
       SOURCE: the census measures through record_rate; this form gates it."""
    _require_urm_surface(s)
    kT = G.KB * s.T
    with _np.errstate(over='ignore'):
        return float(_np.exp(s.E_b / kT) / (s.f0 * (1.0 + _np.exp(-s.dE / kT))))


# =====================================================================================
# 3  THE CENSUS  (C-86; the ProjectModel-ready layer method)
# =====================================================================================
def _require_urm_surface(s):
    """Enforce D-25 at every public surface-consuming path, even after a bypass.

    Construction through URM.surface is the normal entry path, but provenance must be
    defended where the observation is consumed as well: RecordSurface is intentionally
    importable for sealed legacy lanes, so constructor-only enforcement is bypassable.
    """
    tier = getattr(s, 'tier', None)
    provenance = getattr(s, 'provenance', None)
    if tier == "corner":
        if provenance != "DEF-A":
            raise ValueError(
                "COUNTLAW REFUSES: a corner census surface must carry provenance='DEF-A' "
                "from URM.surface; an exact idealisation may not silently pose as world data "
                "(D-25).")
        return
    if tier == "world":
        if not provenance or not str(provenance).strip():
            raise ValueError(
                "COUNTLAW REFUSES: a world-tier census surface requires pinned provenance "
                "from URM.surface (D-25).")
        return
    raise ValueError(
        "COUNTLAW REFUSES: every census surface must enter through URM.surface with an "
        "explicit 'world' or 'corner' tier and the corresponding provenance (D-25).")


def census(surfaces, t_m):
    """THE SURVIVING-RECORD COUNT LAW k(t_m) -- C-86, the URM's wholly-owned falsifiable
       count law (LANE_T47_A_WIDTH section 6 + LANE_T47_B_STAIRCASE sections C/D/G, all
       sealed; registered via LANE_T47_D_REGISTER).

       surfaces: a LIST of RecordSurface objects -- a real device census enters here,
       each record surface built through URM.surface() so D-25 provenance is already on
       it (this is the machinery the Saira/Woodside grounding lanes call).  The census
       independently rechecks tier/provenance and REFUSES bypassed raw or mutated
       surfaces; constructor-only enforcement is not treated as a security boundary.
       t_m: the retention spec, in s.

       THE SIGNATURE IS THE C-76 GATE: (surfaces, t_m) and NOTHING else. No width,
       tolerance, margin, or clustering parameter exists on this path -- the width is
       the carrier's own, computed inside. checks_countlaw.py gates the unreachability.

       Returns DATA, both routes, so agreement is CHECKED, never assumed (the corner()
       pattern):
         k            the count by the INSTRUMENT (record-mode rate <= 1/t_m)
         k_formula    the count by the derived width (dE_i <= delta_pop,i(t_m))
         schedule     one row per readable surface: constants, t_star (instrument),
                      t_star_formula (the check), delta_pop at this t_m, alive both
                      routes, m_eq = tanh(dE/2kT), provenance/tier carried through
         drop_order   surface indices, first-to-die first -- the dated staircase
         departure    sum over dead rows of m_eq: the C-86 departure term (remanence
                      persists while records die; 0 on symmetric carriers)
         delta_coh    hbar/t_m, the coherence-face width of the same spec
         declined     surfaces the model cannot read, declared, never silently counted"""
    inv_tm = 1.0 / float(t_m)
    rows, declined = [], []
    for i, s in enumerate(surfaces):
        _require_urm_surface(s)
        rate = record_rate(s)
        if rate is None:
            declined.append(dict(index=i, name=s.name,
                                 why="not a thermally activated two-state record (or rates "
                                     "below float range); the model declines"))
            continue
        kT = G.KB * s.T
        B = s.E_b + s.dE                      # exact convention map: B from the lower well
        d = delta_pop(B, s.T, s.f0, t_m)
        rows.append(dict(index=i, name=s.name, mechanism=s.mechanism,
                         tier=getattr(s, 'tier', None),
                         provenance=getattr(s, 'provenance', None),
                         B=B, dE=s.dE, E_b=s.E_b, T=s.T, f0=s.f0,
                         t_star=1.0 / rate,
                         t_star_formula=drop_time_formula(s),
                         delta_pop=d,
                         alive=bool(rate <= inv_tm),
                         alive_formula=bool(d is not None and s.dE <= d),
                         m_eq=float(_np.tanh(s.dE / (2.0 * kT)))))
    order = sorted(range(len(rows)), key=lambda j: rows[j]['t_star'])
    return dict(t_m=float(t_m),
                k=sum(1 for r in rows if r['alive']),
                k_formula=sum(1 for r in rows if r['alive_formula']),
                delta_coh=delta_coh(t_m),
                schedule=rows,
                drop_order=[rows[j]['index'] for j in order],
                departure=float(sum(r['m_eq'] for r in rows if not r['alive'])),
                declined=declined)


# =====================================================================================
# 4  THE T-31 MULTI-SHELL CONTROL CARRIER  (LANE_T47_B_STAIRCASE section F, sealed;
#    ported verbatim-in-substance).  The [[4,2,2]]-class 16-dim carrier on which T-31
#    recorded the kill (generic asymmetry -> exact-multiplicity count 0 forever) and
#    C-86 measured the staircase that replaced it (2 -> 1 -> 0).  Carrier units, all
#    DECLARED as the lane declared them: hbar = 1, f0 = 1 (Metropolis attempt rate),
#    kT = 0.20 x the carrier's own unit; bath = Davies-Metropolis single-qubit X_j, Z_j
#    couplings, detailed balance.  The two numerical floors below are the sealed lane's
#    own declared INSTRUMENT floors (not physical widths, and not caller-adjustable):
T31_KT = 0.20
T31_F0 = 1.0
_T31_BOHR_FLOOR = 1e-9      # Davies Bohr-frequency grouping / exact-degeneracy floor (sealed)
_T31_ELEM_FLOOR = 1e-12     # matrix-element zero floor inside the Davies build (sealed)

_I2 = _np.eye(2, dtype=complex)
_PX = _np.array([[0, 1], [1, 0]], dtype=complex)
_PZ = _np.array([[1, 0], [0, -1]], dtype=complex)


def _t31_word(spec):
    M = _np.array([[1.0]], dtype=complex)
    for c in spec:
        M = _np.kron(M, {'I': _I2, 'X': _PX, 'Z': _PZ}[c])
    return M


_T31_HSYM = -(_t31_word('XXXX') + _t31_word('ZZZZ'))
_T31_COUPLE = [_t31_word(''.join(p if i == j else 'I' for i in range(4)))
               for j in range(4) for p in 'XZ']
_T31_ZBAR = [_t31_word('ZZII'), _t31_word('ZIZI')]
_T31_BASINS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def t31_carrier(eps):
    """The T-31 carrier: Hsym = -(XXXX+ZZZZ) plus GENERIC single-site Z fields
       eps*(1 + 0.618 j) -- the lane's own perturbation, golden-ratio coefficients so no
       accidental degeneracy survives (LANE_T47_B_STAIRCASE section F, sealed)."""
    return _T31_HSYM + eps * sum(
        (1.0 + 0.6180339887 * j) *
        _t31_word(''.join('Z' if i == j else 'I' for i in range(4)))
        for j in range(4))


def _t31_davies(H):
    """Davies-Metropolis bath on the carrier's own Bohr frequencies (sealed port)."""
    E, V = _np.linalg.eigh(H)
    Ls = []
    for A in _T31_COUPLE:
        At = V.conj().T @ A @ V
        groups = {}
        for m in range(16):
            for n in range(16):
                if abs(At[m, n]) < _T31_ELEM_FLOOR:
                    continue
                key = round(float(E[m] - E[n]), 9)
                groups.setdefault(key, _np.zeros((16, 16), complex))[m, n] = At[m, n]
        for wk, Mat in groups.items():
            gam = T31_F0 * _np.exp(-max(wk, 0.0) / T31_KT)
            Ls.append(_np.sqrt(gam) * (V @ Mat @ V.conj().T))
    return E, V, Ls


def _t31_basin_instrument(H, Ls):
    """The carrier's own slow sector, no formula and no rate theory (sealed port):
       basins = the joint (Zbar_1, Zbar_2) sectors; the 4-state generator Q is read off
       the FULL 256-dim Liouvillian's own 4 slowest eigenmodes through the basin readout
       -- exact at slow times by construction. On multi-shell carriers THIS is the
       record's own mode (register note: clause_ii's uniform trace weighting misreads
       them). Q must come out a CLASSICAL MARKOV GENERATOR -- gen_dev measures the
       deviation, gated by the checks, so the classical 4-well structure is MEASURED,
       not assumed. Per-record lifetime = the record's own mode of its lumped two-state
       chain, lumping weighted by the carrier's own conditional Gibbs."""
    n = H.shape[0]
    Ev, Vv = _np.linalg.eigh(H)
    gib = Vv @ _np.diag(_np.exp(-(Ev - Ev.min()) / T31_KT)) @ Vv.conj().T
    Pi = {(s1, s2): 0.25 * ((_np.eye(n) + s1 * _T31_ZBAR[0]) @
                            (_np.eye(n) + s2 * _T31_ZBAR[1]))
          for s1 in (1, -1) for s2 in (1, -1)}
    L = G.liouvillian(H, Ls)
    w, U = _np.linalg.eig(L)
    idx = _np.argsort(_np.abs(_np.real(w)))[:4]
    lam = _np.real(w[idx])
    imag_dev = float(_np.max(_np.abs(_np.imag(w[idx]))))
    Y = _np.zeros((4, 4))
    for col, j in enumerate(idx):
        Mj = U[:, j].reshape(n, n, order='F')
        for i, b in enumerate(_T31_BASINS):
            Y[i, col] = float(_np.real(_np.trace(Pi[b] @ Mj)))
    Q = Y @ _np.diag(lam) @ _np.linalg.inv(Y)
    gen_dev = max(float(_np.max(_np.abs(Q.sum(axis=0)))),
                  float(-min(0.0, float(_np.min(Q - _np.diag(_np.diag(Q)))))))
    pb = _np.array([float(_np.real(_np.trace(Pi[b] @ gib))) for b in _T31_BASINS])
    pb = pb / pb.sum()
    taus = []
    for rec in (0, 1):
        rate = 0.0
        for sv in (1, -1):          # escape from value sv, spectator at conditional equil.
            jdx = [j for j, b in enumerate(_T31_BASINS) if b[rec] == sv]
            wts = pb[jdx] / pb[jdx].sum()
            out_rate = [sum(Q[i, j] for i, b in enumerate(_T31_BASINS) if b[rec] != sv)
                        for j in jdx]
            rate += float(_np.dot(wts, out_rate))
        taus.append(1.0 / rate)
    return Q, taus, gen_dev, _np.sort(_np.abs(lam)), imag_dev


def _t31_shell_data(E, V):
    """(dE_i, B_i - dE_i) read off the carrier's OWN spectrum (sealed port): bottom
       shell labeled by the exact conserved Zbar_i; saddle = lowest excited-shell level."""
    order = _np.argsort(E)
    bot, sad = order[:4], E[order[4]]
    out = []
    for Rm in _T31_ZBAR:
        Elo = {+1: _np.inf, -1: _np.inf}
        for nn in bot:
            sgn = int(round(float(_np.real(V[:, nn].conj() @ (Rm @ V[:, nn])))))
            Elo[sgn] = min(Elo[sgn], E[nn])
        dE = abs(Elo[+1] - Elo[-1])
        Eup = max(Elo[+1], Elo[-1])
        out.append(dict(dE=float(dE), B_minus_dE=float(sad - Eup)))
    return out


def _v2(m):
    """2-adic valuation. THE CORNER LAW k = min_E v2(m_E) IS HOMED AT model/count_law.py
       (C-14, recovered as C-86's zero-splitting corner) -- referenced, not duplicated:
       this private helper exists only so the T-31 control below can reproduce the
       corner instrument's reading BESIDE the record-mode census. It is not exported
       through census() or any public path."""
    k = 0
    while m % 2 == 0 and m > 0:
        m //= 2
        k += 1
    return k


def _exact_multiplicities(levels):
    """Exact spectral multiplicities at the sealed lane's declared NUMERICAL degeneracy
       floor (_T31_BOHR_FLOOR = 1e-9 -- an instrument floor against fp noise, NOT a
       physical clustering width; C-76 is about exactly the difference). Private: the
       chosen-width count is not a public capability of this module."""
    lv = _np.sort(_np.asarray(levels))
    out = [[lv[0], 1]]
    for x in lv[1:]:
        if x - out[-1][0] <= _T31_BOHR_FLOOR:
            out[-1][1] += 1
            out[-1][0] = x
        else:
            out.append([x, 1])
    return [m for _, m in out]


_T31_MEMO = {}


def t31_basin(eps):
    """The T-31 control, one call (LANE_T47_B_STAIRCASE section F, sealed): builds the
       carrier at asymmetry eps, reads the basin-lumped record modes, and returns beside
       them the exact-multiplicity corner reading (the C-76 kill: 0 at every eps > 0)
       so the checks can gate the kill AND what survived it. Returns DATA:
       taus (the two records' lifetimes, carrier units), Q, gen_dev, imag_dev,
       slow_rates, dE / B_minus_dE (read off the carrier's own spectrum),
       multiplicities and v2_exact (the corner instrument's reading, see _v2)."""
    if eps in _T31_MEMO:
        return _T31_MEMO[eps]
    H = t31_carrier(eps)
    E, V, Ls = _t31_davies(H)
    Q, taus, gen_dev, slow, imag_dev = _t31_basin_instrument(H, Ls)
    if eps == 0.0:
        Es = _np.sort(E)
        sd = [dict(dE=0.0, B_minus_dE=float(Es[4] - Es[0]))] * 2
    else:
        sd = _t31_shell_data(E, V)
    mults = _exact_multiplicities(_np.linalg.eigvalsh(H))
    out = dict(eps=float(eps), taus=[float(t) for t in taus], Q=Q,
               gen_dev=float(gen_dev), imag_dev=float(imag_dev),
               slow_rates=[float(x) for x in slow],
               dE=[s['dE'] for s in sd], B_minus_dE=[s['B_minus_dE'] for s in sd],
               multiplicities=mults, v2_exact=min(_v2(m) for m in mults))
    _T31_MEMO[eps] = out
    return out


def t31_staircase(eps, t_grid):
    """k(t_m) on the T-31 carrier over t_grid (carrier units): the record-mode count
       from the basin instrument -- the staircase that replaced T-31's binary kill."""
    taus = t31_basin(eps)['taus']
    return [int(sum(1 for tau in taus if tau >= t)) for t in t_grid]


# =====================================================================================
if __name__ == '__main__':
    # standalone demonstration: the sealed six-record ensemble (LANE_T47_A_WIDTH sec. 6)
    from project_model import URM
    ens = [URM.surface("r%d" % i, "thermal two-well (declared corner grid)",
                       dE * EV, (1.2 - dE) * EV, 300.0, 1e9,
                       provenance="DEF-A", tier="corner")
           for i, dE in enumerate((0.0, 0.05, 0.10, 0.15, 0.20, 0.25))]
    for t_m, label in ((3.156e4, "~9 h"), (3.156e7, "1 y"), (3.156e8, "10 y"),
                       (3.156e9, "100 y"), (3.156e11, "10 ky")):
        c = census(ens, t_m)
        print(f"  t_m = {label:>6}: k = {c['k']} (formula route {c['k_formula']}), "
              f"drop order {c['drop_order']}, departure {c['departure']:.4f}")
    print("  drop schedule (s):",
          ["%.4e" % r['t_star'] for r in census(ens, 3.156e8)['schedule']])
