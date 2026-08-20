"""O-50 D  PART 3 -- ESCAPE (1): A FUNCTIONAL OF THE CONFIGURATION *AND* THE CARRIER.

On the torus a record is a HOMOLOGY CLASS; each class has 2^rank(S) REPRESENTATIVES, which
differ in weight and in where they lie.  Does anything gauge-invariant depend on the
representative?  And can carrier geometry -- WHERE the records sit -- rescue a functional
from the cancellation law?
"""
import sys, itertools, math
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE")
from o50d_common import *

say("=" * 104)
say("O-50 D  PART 3   ESCAPE (1)  CONFIGURATION x CARRIER")
say("=" * 104)

T = Torus(2); n = T.nq
pairs = symplectic_logicals(T.stab, n)
basis = [x for p in pairs for x in p]
def comb(coef):
    v = [0] * (2 * n)
    for c, b in zip(coef, basis):
        if c: v = [(x + y) % 2 for x, y in zip(v, b)]
    return v
R1v = comb((0, 0, 0, 1)); R2v = comb((0, 1, 0, 0))          # as SEARCHED in part 1
sb, _ = rref(T.stab, 2 * n)

# ---------------------------------------------------------------- 1. the gauge orbit
say("")
say("1. THE REPRESENTATIVE IS A GAUGE ORBIT, AND THE GAUGE GROUP IS INSIDE THE ADMISSIBLE GROUP.")
# CAUTION, and this lane hit it: the F_2 XOR of two Pauli vectors loses the PHASE.  A
# representative is the OPERATOR PRODUCT R.s, not the Pauli labelled by the XOR -- those can
# differ by a sign, and a sign here is the whole question.  Build the products as MATRICES.
Sgrp = []
for coef in itertools.product((0, 1), repeat=len(sb)):
    Mx = np.eye(2 ** n, dtype=complex); w = [0] * (2 * n)
    for c, s in zip(coef, sb):
        if c:
            Mx = Mx @ dense(s, n); w = [(x + y) % 2 for x, y in zip(w, s)]
    Sgrp.append((Mx, w))
R1m = dense(R1v, n)
reps = [(R1m @ Mx, [(x + y) % 2 for x, y in zip(R1v, w)]) for Mx, w in Sgrp]
wts = sorted({T.weight(w) for _, w in reps})
say(f"   L=2, record R1: {len(reps)} representatives of the SAME class; weights present {wts}; "
    f"min {min(T.weight(w) for _, w in reps)} = d = L")
Pg = np.eye(2 ** n, dtype=complex)
for s in T.stab: Pg = Pg @ (np.eye(2 ** n) + dense(s, n)) / 2
gs = int(round(np.real(np.trace(Pg))))
say(f"   code space dimension {gs}")
say(f"   sanity: (1/|S|) sum_s s  vs  P_code : "
    f"{np.linalg.norm(sum(Mx for Mx, _ in Sgrp) / len(Sgrp) - Pg):.2e}")
M0 = Pg @ R1m @ Pg
diffs = [np.linalg.norm(Pg @ Rm @ Pg - M0) for Rm, _ in reps]
diffs_xor = [np.linalg.norm(Pg @ dense(w, n) @ Pg - M0) for _, w in reps]
say(f"   max over ALL {len(reps)} representatives of || P_code R' P_code - P_code R P_code || = "
    f"{max(diffs):.3e}    (operator products R.s)")
say(f"   the same quantity computed from the F_2 XOR labels instead: {max(diffs_xor):.3e} -- the")
say(f"   XOR drops the Pauli phase and manufactures a spurious sign.  D-19 in miniature.")
say("   EVERY representative acts IDENTICALLY on the code space.  EXACT: R' = R s with s in S,")
say("   and s P_code = P_code.  No code-space quantity can see the representative.")

say("")
say("   OFF the code space they DIFFER -- so ask what the difference is worth:")
rng = np.random.default_rng(1)
psi = rng.normal(size=2 ** n) + 1j * rng.normal(size=2 ** n); psi /= np.linalg.norm(psi)
vals = np.array([np.real(psi.conj() @ Rm @ psi) for Rm, _ in reps])
gauge_avg = vals.mean()
code_val = np.real(psi.conj() @ (R1m @ Pg) @ psi)
say(f"   random state in the FULL 256-dim space: <R'> over representatives ranges "
    f"[{vals.min():+.4f}, {vals.max():+.4f}]")
say(f"   average over the gauge orbit = {gauge_avg:+.8f}   <R P_code> = {code_val:+.8f}   "
    f"difference {abs(gauge_avg - code_val):.2e}")
say("   EXACT IDENTITY: (1/|S|) sum_{s in S} s = P_code, so the gauge average of ANY")
say("   representative-dependent quantity IS its code-space value.  Representative-dependence")
say("   carries NO information beyond the class.  ESCAPE (1a) CLOSED, exactly.")
say("   And the orbit-averaging lemma of Part 2 applies a SECOND time here, with G = S: a")
say("   representative-dependent quantity is either gauge-invariant (blind) or averages away.")

# ---------------------------------------------------------------- 2. class invariants
say("")
say("2. WHAT THE CARRIER *DOES* CONTRIBUTE: class invariants.  Are any of them record-DEPENDENT?")
say(f"   {'L':>3}{'min weight of R1 class':>24}{'min weight of R2 class':>24}"
    f"{'depends on the VALUE s?':>26}")
def min_class_weight(TT, rep, gens):
    best = 10 ** 9
    for coef in itertools.product((0, 1), repeat=len(gens)):
        sset = set(rep)
        for c, g in zip(coef, gens):
            if c: sset ^= set(g)
        if sset: best = min(best, len(sset))
    return best
for L in (2, 3, 4):
    TT = Torus(L)
    mz = min_class_weight(TT, [TT.h(0, j) for j in range(L)], TT.plaq)
    mx = min_class_weight(TT, [TT.v(0, j) for j in range(L)], TT.star)
    say(f"   {L:>3}{mz:>24}{mx:>24}{'NO -- s does not enter':>26}")
say("   The code distance, the minimum weight, the homology class, the writer's minimum weight:")
say("   all are properties of the CLASS and of the CARRIER.  None of them contains s.  They are")
say("   invariant under every write -- they are in the Pf branch of the lemma, non-responsive.")

# ---------------------------------------------------------------- 3. geometry-weighted functionals
say("")
say("3. CAN CARRIER GEOMETRY RESCUE A RESPONSIVE FUNCTIONAL?  m disjoint tori, records at")
say("   separations r_ij, functional weighted by a POWER LAW -- criterion (e)'s own form.")
say("   (D-22: the torus has a metric and no permutation symmetry, so separation is readable.)")
def geo_table(m, L=3, delta=4.0, p=1.0):
    k = 2 * m
    pos = [delta * (i // 2) + 0.5 * (i % 2) for i in range(k)]
    W = np.zeros((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            r = max(1.0, abs(pos[i] - pos[j]))
            W[i, j] = W[j, i] = r ** (-p)
    cfgs = np.array(list(itertools.product((1, -1), repeat=k))) if k <= 16 else None
    tot = sum(W[i, j] for i in range(k) for j in range(i + 1, k))
    if cfgs is not None:
        F = np.array([sum(W[i, j] * s[i] * s[j] for i in range(k) for j in range(i + 1, k))
                      for s in cfgs])
    else:
        rng = np.random.default_rng(0)
        S = rng.integers(0, 2, size=(200000, k)) * 2 - 1
        F = np.einsum('ai,ij,aj->a', S, np.triu(W, 1), S)
    return k, tot, F

say(f"   {'m tori':>7}{'k records':>10}{'mean of F':>18}{'min F':>10}{'max F':>10}"
    f"{'sign-definite':>15}{'coherence mean|F|/sum|w|':>26}{'method':>12}")
for m in (1, 2, 3, 4, 5, 6, 8, 12, 20):
    k, tot, F = geo_table(m)
    sd = (F.min() * F.max() > 0)
    meth = 'EXACT' if k <= 16 else 'sampled 2e5'
    say(f"   {m:>7}{k:>10}{F.mean():>18.10f}{F.min():>10.4f}{F.max():>10.4f}{str(sd):>15}"
        f"{np.abs(F).mean() / tot:>26.6f}{meth:>12}")
say("   MEAN EXACTLY ZERO at every m, both signs always present, and the coherence FALLS as")
say("   records are added.  Putting a power law in by hand (INSERTED, not induced) does not help:")
say("   the geometry sits in the COEFFICIENTS w_ij, and the cancellation is a fact about the")
say("   CHARACTERS chi_ij(s) = s_i s_j, which no choice of coefficients can make sign-definite.")

# ---------------------------------------------------------------- 4. the sign-definite dodge
say("")
say("4. THE ONLY WAY TO MAKE IT SIGN-DEFINITE IS TO ADD A CONSTANT -- AND THE CONSTANT IS")
say("   EXACTLY THE RECORD-BLIND PART.  This is the whole trichotomy in one table.")
k, tot, F = geo_table(4)
G = 0.5 * (tot + F)                                    # sum_ij (1 + s_i s_j)/2 / r_ij  >= 0
say(f"   G(s) = sum_{{i<j}} (1 + s_i s_j) / (2 r_ij)   on k={k} records")
say(f"     G >= 0 everywhere: {bool(G.min() >= -1e-12)}   min {G.min():.6f}  max {G.max():.6f}  "
    f"mean {G.mean():.6f}")
say(f"     invariant part (orbit average) = {G.mean():.6f} = tot/2 = {tot/2:.6f}  -- contains NO s")
say(f"     responsive part G - Pf: mean {np.mean(G - G.mean()):+.3e}  min {np.min(G-G.mean()):+.4f}"
    f"  max {np.max(G-G.mean()):+.4f}  sign-definite: "
    f"{bool(np.min(G-G.mean())*np.max(G-G.mean())>0)}")
say("")
say("   CONTROL (C-61).  The same G evaluated on the +-1 LABELS of the zero-record control carrier")
say("   -- non-degenerate spectrum, ZERO records by P-1 -- with the same geometric weights:")
Hc, J, h, diag = control_carrier(8, seed=5)
say(f"     control carrier n=8: {len(set(np.round(diag,9)))}/{len(diag)} distinct energies -> "
    f"NON-DEGENERATE -> zero records (P-1)")
kk = 8
posc = [4.0 * (i // 2) + 0.5 * (i % 2) for i in range(kk)]
Wc = np.zeros((kk, kk))
for i in range(kk):
    for j in range(i + 1, kk):
        Wc[i, j] = Wc[j, i] = max(1.0, abs(posc[i] - posc[j])) ** -1.0
labels = np.array([[1 - 2 * ((b >> (kk - 1 - q)) & 1) for q in range(kk)] for b in range(2 ** kk)])
Gc = np.array([0.5 * (Wc.sum() / 2 + sum(Wc[i, j] * s[i] * s[j]
                                         for i in range(kk) for j in range(i + 1, kk)))
               for s in labels])
say(f"     G on the CONTROL's labels: min {Gc.min():.6f}  max {Gc.max():.6f}  mean {Gc.mean():.6f}")
say(f"     G on the RECORD carrier   : min {G.min():.6f}  max {G.max():.6f}  mean {G.mean():.6f}")
say("     IDENTICAL IN FORM.  The sign-definite, non-cancelling part of G is present on a carrier")
say("     that holds ZERO records.  IT IS RECORD-BLIND.  ESCAPE (1b) CLOSED.")
say("")
say("   VERDICT ON ESCAPE (1): a functional of the configuration AND the carrier splits into a")
say("   carrier-only part (record-blind: it survives the C-61 control) and a configuration-")
say("   responsive part (mean exactly zero, never sign-definite).  Carrier data changes the")
say("   COEFFICIENTS and never the CHARACTERS.  NOT AN ESCAPE.")
say("=" * 104)
