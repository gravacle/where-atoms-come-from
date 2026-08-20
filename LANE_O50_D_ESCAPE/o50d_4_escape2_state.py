"""O-50 D  PART 4 -- ESCAPE (2): A FUNCTIONAL OF THE STATE RATHER THAN OF THE CONFIGURATION.

Off-diagonal elements, coherences between configurations, entanglement across a cut.  These
are NOT functions of the +-1 values, so the configuration theorem does not reach them.
Everything below is on the torus (D-23).  Dense at L=2 (dim 256); exact F_2 stabiliser
entropies to L=6.
"""
import sys, itertools, math
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE")
from o50d_common import *

say("=" * 104)
say("O-50 D  PART 4   ESCAPE (2)  STATE FUNCTIONALS")
say("=" * 104)

# ---------------------------------------------------------------- setup at L=2
T = Torus(2); n = T.nq; N = 2 ** n
pairs = symplectic_logicals(T.stab, n); basisv = [x for p in pairs for x in p]
def comb(coef):
    v = [0] * (2 * n)
    for c, b in zip(coef, basisv):
        if c: v = [(x + y) % 2 for x, y in zip(v, b)]
    return v
R1v = comb((0, 0, 0, 1)); R2v = comb((0, 1, 0, 0))
R1 = dense(R1v, n); R2 = dense(R2v, n)
H = -sum(dense(s, n) for s in T.stab)
Pg = np.eye(N, dtype=complex)
for s in T.stab: Pg = Pg @ (np.eye(N) + dense(s, n)) / 2
w, V = np.linalg.eigh(Pg); V = V[:, w > 0.5]            # 256 x 4 code basis
# configuration basis: joint eigenvectors of R1,R2 inside the code space
r1 = V.conj().T @ R1 @ V; r2 = V.conj().T @ R2 @ V
ww, U1 = np.linalg.eigh(r1 + 3.0 * r2)
CB = V @ U1                                              # columns = configuration states
cfg = [(int(round(np.real(CB[:, a].conj() @ R1 @ CB[:, a]))),
        int(round(np.real(CB[:, a].conj() @ R2 @ CB[:, a])))) for a in range(4)]
say("")
say(f"1. SETUP.  L=2 toric code, dim {N}; code space dim {CB.shape[1]}; configuration labels "
    f"(s1,s2) = {cfg}")
def logmat(O): return CB.conj().T @ O @ CB           # 4x4 logical matrix
LOGICALS = {}
for coef in itertools.product((0, 1), repeat=4):
    LOGICALS[coef] = dense(comb(coef), n)
W1 = [(c, M) for c, M in LOGICALS.items() if T.sp(comb(c), R1v) == 1 and T.sp(comb(c), R2v) == 0]
W2 = [(c, M) for c, M in LOGICALS.items() if T.sp(comb(c), R2v) == 1 and T.sp(comb(c), R1v) == 0]
say(f"   writers of R1 SEARCHED: {len(W1)} classes {[c for c,_ in W1]}")
for c, M in W1:
    assert np.linalg.norm(M @ H - H @ M) < 1e-9 and np.linalg.norm(M @ R1 + R1 @ M) < 1e-9
say("   every one verified: [W,H] = 0 and {W,R1} = 0.  ALL FOUR ARE EQUALLY ADMISSIBLE WRITERS.")

# ---------------------------------------------------------------- states
def code_state(amps):
    v = CB @ np.array(amps, dtype=complex); return v / np.linalg.norm(v)
def rho(v): return np.outer(v, v.conj())
def dephase_logical(rL, p):
    """INSERTED logical dephasing on the 4x4 logical density matrix"""
    out = rL.copy()
    for a in range(4):
        for b in range(4):
            if a != b: out[a, b] *= (1 - p)
    return out

# ---------------------------------------------------------------- functionals
def partial_trace(r, keep, nq):
    t = r.reshape([2] * nq + [2] * nq)
    drop = [q for q in range(nq) if q not in keep]
    for q in reversed(drop):
        t = np.trace(t, axis1=q, axis2=q + t.ndim // 2)
    d = 2 ** len(keep); return t.reshape(d, d)
def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-12]
    return float(-(e * np.log2(e)).sum())

supp_R1 = T.support(R1v); supp_R2 = T.support(R2v)
regA = supp_R1                                            # NON-contractible: it carries R1
regB = [q for q in range(n) if q not in supp_R1 and q not in supp_R2][:2]
say(f"   region A = supp(R1) = {regA} (wraps the torus, carries a logical); region B = {regB}.")
say("   NOTE, honestly: at L=2 the code distance is 2 and NO region of size >= 2 is contractible.")
say("   The contractible case is done exactly in section 5 at L=3..6.")

def functionals(r_full):
    """r_full is a 256x256 density matrix.  Returns the battery."""
    rL = CB.conj().T @ r_full @ CB
    ex = lambda O: float(np.real(np.trace(r_full @ O)))
    out = {}
    out['<R1>'] = ex(R1)
    out['<R2>'] = ex(R2)
    out['Re<Xbar1>'] = ex(LOGICALS[(0, 0, 1, 0)])
    out['|rho_L 01|'] = float(abs(rL[0, 1]))
    out['purity(rho_L)'] = float(np.real(np.trace(rL @ rL)))
    out['S(rho_L)'] = vN(rL / max(np.real(np.trace(rL)), 1e-300))
    out['S_A noncontr'] = vN(partial_trace(r_full, regA, n))
    out['S_B contr'] = vN(partial_trace(r_full, regB, n))
    bl = 0.0
    for c in [(0, 0, 1, 0), (0, 0, 0, 1), (0, 0, 1, 1)]:
        bl += ex(LOGICALS[c]) ** 2
    out['Bloch radius r_1'] = float(bl)
    integ = abs(ex(R1)) + abs(ex(R2))
    out['integrity S|<R_i>|'] = float(integ)
    return out

KEYS = ['<R1>', '<R2>', 'Re<Xbar1>', '|rho_L 01|', 'purity(rho_L)', 'S(rho_L)',
        'S_A noncontr', 'S_B contr', 'Bloch radius r_1', 'integrity S|<R_i>|']

STATES = {
    'basis |s1=+1,s2=+1>': rho(CB[:, cfg.index((1, 1))]),
    'cat (|++> + |-->)/r2': rho(code_state([1, 0, 0, 1] if cfg[0] == (1, 1) else [1, 0, 0, 1])),
    'generic pure code st': rho(code_state([0.6, 0.3 + 0.2j, -0.4, 0.5j])),
    'dephased p=0.5 (INS)': None,
}
gp = code_state([0.6, 0.3 + 0.2j, -0.4, 0.5j])
rLg = np.outer(CB.conj().T @ gp, (CB.conj().T @ gp).conj())
STATES['dephased p=0.5 (INS)'] = CB @ dephase_logical(rLg, 0.5) @ CB.conj().T

say("")
say("2. THE BATTERY, AND WHETHER A WRITE MOVES IT.  For each functional: its value on a generic")
say("   code state, and the SPREAD of its value over the four equally admissible writers of R1.")
say("")
say(f"   {'functional':<22}{'value':>12}{'after Xbar1':>13}{'after Ybar1':>13}"
    f"{'max|delta| over':>17}{'orbit mean':>13}{'a fn of s?':>12}")
say(f"   {'':<22}{'':>12}{'':>13}{'':>13}{'the 4 writers':>17}{'over G_W':>13}{'':>12}")
r0 = STATES['generic pure code st']
Xb = LOGICALS[(0, 0, 1, 0)]
Yb = None
for c, M in W1:
    if c != (0, 0, 1, 0):
        Yb = M; Yc = c; break
vals0 = functionals(r0)
valsX = functionals(Xb @ r0 @ Xb.conj().T)
valsY = functionals(Yb @ r0 @ Yb.conj().T)
allw = [functionals(M @ r0 @ M.conj().T) for _, M in W1]
grp = [functionals(M @ r0 @ M.conj().T) for M in LOGICALS.values()]
ISFN = {'<R1>': 'YES', '<R2>': 'YES', 'Re<Xbar1>': 'no', '|rho_L 01|': 'no',
        'purity(rho_L)': 'no', 'S(rho_L)': 'no', 'S_A noncontr': 'no', 'S_B contr': 'no',
        'Bloch radius r_1': 'no', 'integrity S|<R_i>|': 'no'}
for kk in KEYS:
    spread = max(abs(v[kk] - vals0[kk]) for v in allw)
    om = float(np.mean([v[kk] for v in grp]))
    say(f"   {kk:<22}{vals0[kk]:>12.6f}{valsX[kk]:>13.6f}{valsY[kk]:>13.6f}"
        f"{spread:>17.6f}{om:>13.6f}{ISFN[kk]:>12}")
say(f"   (Ybar1 here is the writer class {Yc}; it is as admissible as Xbar1 -- see part 1.)")

say("")
say("   THE SAME BATTERY ACROSS A STATE ENSEMBLE -- max |change| over ALL 16 admissible writers:")
say(f"   {'state':<24}" + "".join(f"{k:>21}" for k in
    ['<R1>', 'Re<Xbar1>', '|rho_L 01|', 'S_A noncontr', 'integrity S|<R_i>|']))
for nm, rr in STATES.items():
    row = f"   {nm:<24}"
    for kk in ['<R1>', 'Re<Xbar1>', '|rho_L 01|', 'S_A noncontr', 'integrity S|<R_i>|']:
        base = functionals(rr)[kk]
        mx = max(abs(functionals(M @ rr @ M.conj().T)[kk] - base) for M in LOGICALS.values())
        row += f"{base:>12.4f}/{mx:>8.4f}"
    say(row)
say("   format: value / max change under a write.  <R1> moves (it is the record).  The coherence")
say("   Re<Xbar1> moves by twice its value and averages to zero.  |rho_L 01| moves only because a")
say("   write PERMUTES the configuration labels, so a NAMED matrix element is carried to another")
say("   one -- the multiset of |rho_ab| is invariant (section 4).  S_A and the integrity do not")
say("   move at all, for any state in the ensemble.")

say("")
say("3. THE TWO THINGS THAT TABLE SHOWS")
say("   (a) EVERY functional that MOVES under a write has orbit mean 0 (or, for the ones built")
say("       from squares and absolute values, moves not at all).  Part 1 proved why: the only")
say("       logical observable fixed by all of G_W is the IDENTITY, so every linear state")
say("       functional -- diagonal or OFF-diagonal -- averages to exactly zero over the writers.")
say("   (b) THE WRITER AMBIGUITY.  Xbar1 and Ybar1 are BOTH admissible writers of record 1 and")
say("       both leave record 2 alone.  They send the same coherence to OPPOSITE values.  So a")
say("       coherence is not even a well-defined function of 'record 1 was flipped'.")

say("")
say("4. THE COMPLETE INVARIANT: which state functionals survive every write?")
rng = np.random.default_rng(4)
tests = []
for _ in range(200):
    a = rng.normal(size=4) + 1j * rng.normal(size=4)
    r = rho(code_state(a))
    r = 0.8 * r + 0.2 * (CB @ (np.eye(4) / 4) @ CB.conj().T)
    tests.append(r)
maxdev = {}
for c, M in LOGICALS.items():
    for r in tests[:20]:
        rr = M @ r @ M.conj().T
        for cc, MM in LOGICALS.items():
            k = 'abs<P>'
            d = abs(abs(np.real(np.trace(rr @ MM))) - abs(np.real(np.trace(r @ MM))))
            maxdev[k] = max(maxdev.get(k, 0.0), d)
        maxdev['purity'] = max(maxdev.get('purity', 0.0),
                               abs(np.real(np.trace(rr @ rr)) - np.real(np.trace(r @ r))))
say(f"   over 20 random mixed code states x all 16 writers: max change in |<P>| for EVERY logical")
say(f"   Pauli P = {maxdev['abs<P>']:.2e}; max change in purity = {maxdev['purity']:.2e}")
# separation: do the |<P>| determine the G_W orbit?
def sig(r): return tuple(round(abs(np.real(np.trace(r @ M))), 9) for M in LOGICALS.values())
same = sum(1 for r in tests[:40] for _, M in list(LOGICALS.items())[1:]
           if sig(M @ r @ M.conj().T) == sig(r))
say(f"   and the signature |<P>| is constant on every writer orbit in {same}/{40*15} checks.")
say("   EXACT CHARACTERISATION: conjugation by G_W = the full logical Pauli group multiplies each")
say("   logical Pauli coefficient by +-1 and fixes the identity coefficient.  Hence THE COMPLETE")
say("   ALGEBRA OF WRITE-INVARIANT STATE FUNCTIONALS IS THE ALGEBRA OF FUNCTIONS OF THE ABSOLUTE")
say("   VALUES |<Pbar>|.  Every one of them is blind to WHICH configuration the state is in and")
say("   sensitive only to HOW SHARP the record is.  Escape (2) yields exactly one kind of")
say("   survivor: an INTEGRITY functional, never a VALUE functional.")

# ---------------------------------------------------------------- entanglement, exact, all L
say("")
say("5. ENTANGLEMENT ACROSS A CUT -- EXACT, AT EVERY L, FOR EVERY CONFIGURATION.")
say("   For a stabiliser codeword, S_A = |A| - dim{ g in Stab : supp(g) subset A }, and that")
say("   dimension is computed from the (x|z) LABELS ONLY.  The configuration enters the")
say("   stabiliser group only through SIGNS.  Therefore S_A IS IDENTICAL FOR ALL 2^k")
say("   CONFIGURATIONS, at every L, for every region.  EXACT -- no computation can disagree.")
def S_A(TT, gens, region):
    """S_A = |A| - dim{ g in Stab : supp(g) subset A }.  THE GENERATORS MUST BE INDEPENDENT:
       the 2L^2 star/plaquette operators carry TWO relations, and feeding them in raw inflates
       the dimension by 2 and returns negative entropies.  This lane hit exactly that."""
    nq = TT.nq; inside = set(region)
    G, _ = rref(gens, 2 * nq)                      # independent generating set
    assert len(G) == nq, f"expected {nq} independent generators, got {len(G)}"
    cons = []
    for e in range(nq):
        if e in inside: continue
        cons.append([g[e] for g in G]); cons.append([g[nq + e] for g in G])
    ns = nullspace2(cons, len(G)) if cons else [[1 if i == j else 0 for i in range(len(G))]
                                               for j in range(len(G))]
    return len(region) - rank2(ns, len(G))
say("")
say(f"   {'L':>3}{'region':>26}{'|A|':>5}{'S_A (bits)':>12}{'same for all 4 configs':>25}")
for L in (2, 3, 4, 5, 6):
    TT = Torus(L)
    pr = symplectic_logicals(TT.stab, TT.nq)
    bs = [x for p in pr for x in p]
    def cb(coef, TT=TT, bs=bs):
        v = [0] * (2 * TT.nq)
        for c, b in zip(coef, bs):
            if c: v = [(x + y) % 2 for x, y in zip(v, b)]
        return v
    RA = cb((0, 0, 0, 1)); RB = cb((0, 1, 0, 0))
    gens = list(TT.stab) + [RA, RB]
    for name, reg in (("contractible 2x2 patch",
                       sorted({TT.h(i, j) for i in range(min(2, L)) for j in range(min(2, L))} |
                              {TT.v(i, j) for i in range(min(2, L)) for j in range(min(2, L))})),
                      ("non-contractible: supp(R1)", TT.support(RA))):
        s = S_A(TT, gens, reg)
        say(f"   {L:>3}{name:>26}{len(reg):>5}{s:>12}{'YES (signs do not enter)':>25}")

say("")
say("   DENSE CONFIRMATION AT L=2, all four configuration states, both regions:")
for a in range(4):
    r = rho(CB[:, a])
    say(f"     config {cfg[a]}:  S_A(noncontr) = {vN(partial_trace(r, regA, n)):.6f}   "
        f"S_B(contr) = {vN(partial_trace(r, regB, n)):.6f}")
say("   Identical across configurations, as the exact argument requires.")

say("")
say("6. CLAUSE (v) IS ITSELF A RECORD-BLINDNESS THEOREM FOR LOCAL STATE FUNCTIONALS.")
say("   rho_A is fixed by the expectations of the Paulis supported in A.  For a stabiliser")
say("   codeword those vanish except on the stabiliser group of the state, Stab = <S, R1, R2>.")
say("   CLAUSE (v) SAYS EXACTLY: for a contractible A, every element of Stab supported in A")
say("   already lies in S -- so its expectation is +1 in EVERY codeword and the signs, which are")
say("   the record values, never appear.  THEREFORE rho_A IS IDENTICAL FOR EVERY CONFIGURATION,")
say("   AND FOR EVERY SUPERPOSITION OF THEM.  Every functional of a contractible region is")
say("   RECORD-BLIND, by protection itself.  Measured, exactly:")
say(f"   {'L':>3}{'region':>22}{'contractible':>14}{'dim Stab_A in S':>17}"
    f"{'dim Stab_A in <S,R>':>21}{'rho_A sees s?':>15}")
for L in (3, 4, 5, 6):
    TT = Torus(L)
    pr = symplectic_logicals(TT.stab, TT.nq); bs = [x for p in pr for x in p]
    def cb2(coef, TT=TT, bs=bs):
        v = [0] * (2 * TT.nq)
        for c, b in zip(coef, bs):
            if c: v = [(x + y) % 2 for x, y in zip(v, b)]
        return v
    RA = cb2((0, 0, 0, 1)); RB = cb2((0, 1, 0, 0))
    def dim_in(gens, reg, TT=TT):
        nq = TT.nq; G, _ = rref(gens, 2 * nq); inside = set(reg); cons = []
        for e in range(nq):
            if e in inside: continue
            cons.append([g[e] for g in G]); cons.append([g[nq + e] for g in G])
        ns = nullspace2(cons, len(G)) if cons else [[1 if i == j else 0 for i in range(len(G))]
                                                   for j in range(len(G))]
        return rank2(ns, len(G))
    regs = [("contractible 2x2", sorted({TT.h(i, j) for i in range(2) for j in range(2)} |
                                        {TT.v(i, j) for i in range(2) for j in range(2)}), True),
            ("contractible (L-1)^2", sorted({TT.h(i, j) for i in range(L - 1) for j in range(L - 1)} |
                                            {TT.v(i, j) for i in range(L - 1) for j in range(L - 1)}), True),
            ("wrapping: supp(R1)", TT.support(RA), False)]
    for nm, reg, contr in regs:
        d0 = dim_in(list(TT.stab), reg); d1 = dim_in(list(TT.stab) + [RA, RB], reg)
        say(f"   {L:>3}{nm:>22}{str(contr):>14}{d0:>17}{d1:>21}"
            f"{('NO' if d0 == d1 else 'YES'):>15}")
say("   For every CONTRACTIBLE region the two dimensions agree: no element of Stab_A carries a")
say("   record sign, so rho_A is configuration-independent.  For the WRAPPING region they differ")
say("   by exactly 1 -- the record itself is inside -- and rho_A then DOES see s.  THE INSTRUMENT")
say("   REGISTERS CONFIGURATION-DEPENDENCE WHERE THERE IS SOME (D-15).")
say("")
say("   DENSE POSITIVE CONTROL AT L=2 on the wrapping region A = supp(R1):")
for a in range(4):
    rA = partial_trace(rho(CB[:, a]), regA, n)
    say(f"     config {str(cfg[a]):>9}:  Tr(rho_A . R1|_A) = "
        f"{np.real(np.trace(rA @ dense_from_sets(set(), set(range(len(regA))), len(regA)))):+.6f}   "
        f"S_A = {vN(rA):.6f}")
say("     the reduced state MOVES with the configuration while its ENTROPY does not.  The")
say("     record-dependence of a non-contractible region is carried entirely by <R_i> -- i.e. by")
say("     the configuration -- so it is back inside the cancellation law of part 2.")
say("")
say("   AND THE C-61 CONTROL, for completeness: the zero-record carrier of C-61 has a PRODUCT")
say("   ground state, S_A = 0 for every region.  It differs from the toric code in far more than")
say("   record content, so it cannot by itself certify record-sensitivity.  The sharp control is")
say("   the contractible/wrapping comparison above, which differs ONLY in whether the record is")
say("   inside the region.")
say("")
say("   VERDICT ON ESCAPE (2): state functionals split exactly as the configuration functionals")
say("   did.  Responsive ones (coherences, off-diagonals, <Pbar>) have orbit mean EXACTLY zero")
say("   AND are not even well defined given only 'record i was flipped' -- two equally admissible")
say("   writers give opposite signs.  Invariant ones are exactly the functions of |<Pbar>|:")
say("   record INTEGRITY, not record VALUE, and they do not respond to writing at all.")
say("   NOT AN ESCAPE from responsiveness.  It IS the cleanest source of a non-record-blind,")
say("   sign-definite, non-cancelling quantity -- which part 7 scores against all five standards.")
say("=" * 104)
