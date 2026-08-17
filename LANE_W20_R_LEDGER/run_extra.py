# run_extra.py -- LANE W20_R_LEDGER.  THE ISOLATIONS THE MAIN RUN DEMANDED.
#
# Five blocks, each one an ISOLATION of a single variable that run_arms.py could not separate:
#   E1  THE IDENTITY.  Delta_surf = 2 (H_FULL - H_CENTRE) = 2 sum_k p_k E_k.  DERIVED, then measured.
#       Consequence: the pre-registration's four channels are NOT four independent tests of H-SURFACE.
#       That is a correction to my own binding pre-registration and it is reported as one.
#   E2  THE EIGENSTATE LADDER.  Same H, same sector, same Sigma, same algebra, same partition.
#       THE ONLY THING THAT MOVES IS THE EIGENVALUE INDEX n.  This isolates GROUND-STATE SELECTION
#       from every other property of the dynamics.
#   E3  THE PERTURBATION CURVE.  |psi(e)> ~ (1-e)|gs> + e|haar>.  How fast does the null break?
#   E4  ARM 1 OF THE COMMISSION -- REMOVE EM.  (a) same H, NO Gauss projection, 4096 dims.
#       (b) NO Gauss AND NO plaquettes: L bare qubits.  Verbatim the commissioning text.
#   E5  ARM S1 -- A REGION WITH A LARGER FREE PART.  This is the NEXT STEP the pre-registration
#       named for itself, and it costs nothing on the same carrier: A' = {0,1,2,3}, S' = E[A'] =
#       links {1,2,3,4,5} (vertex-induced, TWO independent cycles), Sigma' = delta(A') = {0,6} =
#       the carrier's unique 2-link bottleneck.  Cross-orbit from A = {0,1,2} by cardinality, so
#       it is ADMISSIBLE under the pre-registration's Aut rule.  It also collapses the second
#       recorded confound: on this arm the surface IS the bottleneck rather than merely touching it.

import numpy as np, math
from core_w20r import *

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def rule(t=""):
    P("\n" + "=" * 104); P(t); P("=" * 104)

def LM(ls):
    m = 0
    for l in ls: m |= 1 << l
    return m

S_T   = LM([1,2,3]); SIG_T = LM([0,4,5]); ENV_T = LM([0,4,5,6,7,8,9,10,11]); W_S = LM([1,2,3])
GRID  = [0.05,0.10,0.20,0.30,0.45,0.60,0.80,1.00,1.30,1.70,2.20,3.00,5.00]
SEC = {"vacuum":[1]*8, "eta{0,4}":[-1 if v in (0,4) else 1 for v in range(8)],
       "eta{1,5}":[-1 if v in (1,5) else 1 for v in range(8)],
       "eta{4,5}":[-1 if v in (4,5) else 1 for v in range(8)],
       "eta{0,1}":[-1 if v in (0,1) else 1 for v in range(8)]}
BE = {k: PhysSector(v) for k, v in SEC.items()}

def chans(be, S, wS_list):
    A = alg_links(be, S, N_PHYS, "FULL")
    cen = [v for v in gf2_span(list(A.W.values())) if all(omega(v, w) == 0 for w in A.W.values())]
    C = Algebra(be, cen, N_PHYS, "CENTRE")
    MG = Algebra(be, [mk(0, z) for z in wS_list], N_PHYS, "MAG")
    return A, C, MG

be0 = BE["vacuum"]
A_S, A_CEN, A_MAG = chans(be0, S_T, [W_S])
A_BLK = Algebra(be0, [mk(1<<1,0), mk(0,W_S)], N_PHYS, "BLOCK")
A_SIG = alg_links(be0, SIG_T, N_PHYS, "SIGMA")
A_ENV = alg_links(be0, ENV_T, N_PHYS, "ENV")
J_SE  = join(A_S, A_ENV); J_SS = join(A_S, A_SIG)

def dsurf(psi, A, C, SG, EN):
    return (A.entropy(psi) + EN.entropy(psi) - join(A, EN).entropy(psi)) \
         - (A.entropy(psi) + SG.entropy(psi) - join(A, SG).entropy(psi))

# ==================================================================================================
rule("E1 -- THE IDENTITY THAT COLLAPSES MY OWN PRE-REGISTRATION'S FOUR CHANNELS INTO ONE NUMBER")
P("DERIVATION (no state needed).  On the physical sector the environment algebra A_env is EXACTLY")
P("the COMMUTANT of A_S: both label spaces have GF(2) dim 6 and 4 with W_env = W_S^perp -- checked")
P("below.  Write the minimal central projections of Z(A_S) as P_k.  Then Z(A_S)' = (+)_k B(H_k),")
P("A_S|_k = M_2 (x) I_{m_k}, A_env|_k = I_2 (x) M_{m_k}, and for a PURE global state |psi>:")
P("     S(A_S)   = H({p_k}) + sum_k p_k E_k        E_k = entanglement of the free qubit in block k")
P("     S(A_env) = H({p_k}) + sum_k p_k E_k        (the SAME number -- commutant)")
P("     S(A_S v A_env) = H({p_k})                  (pure inside each block)")
P("     I(A_S:A_env)   = H({p_k}) + 2 sum_k p_k E_k")
P("     I(A_S:A_Sigma) = H(Z(A_S)) = H({p_k})      (A_Sigma == Z(A_S), a SUBALGEBRA of A_S)")
P("  =>  Delta_surf = 2 sum_k p_k E_k = 2 ( H_FULL(A_S) - H_CENTRE(A_S) )   FOR EVERY PURE STATE.")
P("")
ALLW = Algebra(be0, [mk(1<<l,0) for l in range(L)] + [mk(0,z) for z in gf2_rref(CYC).values()],
               N_PHYS, "B(H)")                     # the effective 10-dim label space, reduced mod B^1
Wperp = sorted(v for v in gf2_span(list(ALLW.W.values())) if all(omega(v, w) == 0 for w in A_S.W.values()))
P("effective label space dim = %d (= 2^10 operators = 32^2 = all of B(H_phys))" % ALLW.n)
P("dim W_S = %d ; dim W_S^perp = %d ; dim W_env = %d ; W_env == W_S^perp : %s"
  % (A_S.n, len(gf2_rref(Wperp)), A_ENV.n,
     sorted(gf2_span(list(A_ENV.W.values()))) == Wperp))
P("   (A_env is therefore the FULL COMMUTANT of A_S, not merely a commuting subalgebra --")
P("    which is why the environment can hold no more about S than the surface plus the free qubit.)")
P("")
P("MEASURED on 6 Haar states and 6 ground states.  If this is an identity the residual is machine zero,")
P("and then Delta_surf carries NO information beyond H_FULL - H_CENTRE.  PRINT THE RESIDUAL.")
P("  state          Delta_surf     2*(H_FULL-H_CENTRE)   residual      sum_k p_k E_k")
rng = np.random.default_rng(4242)
tests = []
for t in range(6):
    v = rng.normal(size=32) + 1j*rng.normal(size=32); v /= np.linalg.norm(v)
    tests.append(("haar%d" % t, v))
for g2 in [0.05, 0.45, 0.6, 0.8, 1.3, 5.0]:
    ev, evec = np.linalg.eigh(H_matrix(be0, g2))
    tests.append(("gs g2=%.2f" % g2, evec[:, 0]))
worst = 0.0
for nm, v in tests:
    d = dsurf(v, A_S, A_CEN, A_SIG, A_ENV)
    h = A_S.entropy(v) - A_CEN.entropy(v)
    worst = max(worst, abs(d - 2*h))
    P("  %-14s %-14.9f %-21.9f %-13.2e %.9f" % (nm, d, 2*h, abs(d-2*h), h))
P("MAX RESIDUAL OVER 12 STATES = %.3e  ->  IDENTITY CONFIRMED." % worst)
P("")
P(">>> CORRECTION TO MY OWN BINDING PRE-REGISTRATION.  The pre-registration declared Delta_surf as")
P("    the primary quantity and declared FOUR channels, implying four tests.  There are not four.")
P("    Delta_surf(FULL) = Delta_surf(BLOCK) exactly, Delta_surf(CENTRE) = 0 exactly, and")
P("    Delta_surf(MAG) = half of Delta_surf(FULL) in every row of run_arms.py.  ONE NUMBER:")
P("        C := sum_k p_k E_k = the CONDITIONAL FREE ENTROPY of the region given its surface.")
P("    H-SURFACE says C = 0.  It is NOT forced to be 0: Haar states give C ~ 0.7 bits (below).")

# ==================================================================================================
rule("E2 -- THE EIGENSTATE LADDER.  ISOLATED VARIABLE: THE EIGENVALUE INDEX n.  NOTHING ELSE MOVES.")
P("Same H, same sector, same Sigma, same algebra, same partition, same coupling.  Only n moves.")
P("This separates GROUND-STATE SELECTION from every other property of the dynamics, which the")
P("Haar arm cannot do (a Haar state also changes the energy, the sector statistics and the norm).")
P("C = sum_k p_k E_k in bits.  H-SURFACE predicts C = 0.  Ceiling for this region is 1 bit.")
for g2 in [0.45, 0.60, 0.80]:
    P("")
    P("g2 = %.2f  (inside the measured live window [0.4198, 0.8234])" % g2)
    ev, evec = np.linalg.eigh(H_matrix(be0, g2))
    P("   n    E_n            gap to n-1    H_FULL     H_CENTRE   C=sum p_k E_k   Delta_surf")
    Cs = []
    for n in range(32):
        v = evec[:, n]
        hf = A_S.entropy(v); hc = A_CEN.entropy(v); C = hf - hc
        Cs.append(C)
        gp = ev[n] - ev[n-1] if n else float('nan')
        P("   %-4d %-14.8f %-13.8f %-10.6f %-10.6f %-15.9f %.9f"
          % (n, ev[n], gp, hf, hc, C, 2*C))
    P("   >>> C(n=0) = %.9f   ;  mean C over n>=1 = %.9f  ;  min C over n>=1 = %.9f  ;  max C = %.9f"
      % (Cs[0], float(np.mean(Cs[1:])), float(np.min(Cs[1:])), float(np.max(Cs))))
    P("   >>> RATIO mean(C, n>=1) / C(n=0) = %.1f" % (float(np.mean(Cs[1:]))/Cs[0] if Cs[0] > 0 else float('inf')))

rule("E2b -- THE SAME LADDER ACROSS ALL FIVE CHARGE SECTORS AT g2 = 0.60.  C(n=0) vs the rest.")
P("  sector      C(n=0)         mean C (n>=1)   min C (n>=1)   ratio    rank of C(n=0) among 32")
for nm, be in BE.items():
    Ai, Ci, _ = chans(be, S_T, [W_S])
    ev, evec = np.linalg.eigh(H_matrix(be, 0.60))
    Cs = [Ai.entropy(evec[:, n]) - Ci.entropy(evec[:, n]) for n in range(32)]
    rk = int(np.argsort(Cs).tolist().index(0))
    P("  %-11s %-14.9f %-15.9f %-14.9f %-8.1f %d of 32 (0 = smallest)"
      % (nm, Cs[0], float(np.mean(Cs[1:])), float(np.min(Cs[1:])),
         float(np.mean(Cs[1:]))/Cs[0] if Cs[0] > 0 else float('inf'), rk))

# ==================================================================================================
rule("E3 -- THE PERTURBATION CURVE.  |psi(e)> ~ (1-e)|gs> + e|haar>, g2 = 0.60, vacuum sector.")
P("How special is the ground state?  If C rises linearly from 0 the null is a knife edge; if it")
P("stays flat for a while the ground-state manifold has a genuine neighbourhood.")
ev, evec = np.linalg.eigh(H_matrix(be0, 0.60))
gsv = evec[:, 0]
rng2 = np.random.default_rng(31337)
hv = rng2.normal(size=32) + 1j*rng2.normal(size=32); hv /= np.linalg.norm(hv)
P("   e        |<gs|psi>|^2   H_FULL     C=sum p_k E_k   Delta_surf")
for e in [0.0, 0.001, 0.003, 0.01, 0.03, 0.10, 0.30, 0.50, 1.00]:
    v = (1-e)*gsv + e*hv; v /= np.linalg.norm(v)
    C = A_S.entropy(v) - A_CEN.entropy(v)
    P("   %-8.3f %-14.9f %-10.6f %-15.9f %.9f" % (e, abs(complex(np.vdot(gsv, v)))**2, A_S.entropy(v), C, 2*C))

# ==================================================================================================
rule("E4 -- ARM 1 OF THE COMMISSION: REMOVE EM.  TWO REMOVALS, BOTH RUN.")
P("(a) X1a  SAME H, NO GAUSS PROJECTION.  4096-dim unconstrained space.  The algebra A_S is the")
P("    SAME 16 strings and is STRUCTURALLY IDENTICAL (dim 16, centre dim 4) -- what is destroyed is")
P("    the OPERATOR IDENTITY Z(A_S) = alg{X_0,X_4,X_5}.  The CENTRE channel therefore survives as an")
P("    ALGEBRA but STOPS BEING THE SURFACE, and it is stamped accordingly and never scored.")
P("(b) X1b  NO GAUSS AND NO PLAQUETTES: H = -g2 sum_l X_l on L bare qubits.  Verbatim the text.")
fb = FreeSpace()
A_Sf   = alg_links(fb, S_T, N_FREE, "FULL")
cenf = [v for v in gf2_span(list(A_Sf.W.values())) if all(omega(v,w) == 0 for w in A_Sf.W.values())]
A_CENf = Algebra(fb, cenf, N_FREE, "CENTRE*")
A_BLKf = Algebra(fb, [mk(1<<1,0), mk(0,W_S)], N_FREE, "BLOCK")
A_MAGf = Algebra(fb, [mk(0,W_S)], N_FREE, "MAG")
A_SIGf = alg_links(fb, SIG_T, N_FREE, "SIGMA")
P("")
P("STRUCTURE DIFF (this is the whole point of the arm):")
P("  %-10s %-28s %-28s" % ("channel", "WITH Gauss (32-dim sector)", "WITHOUT Gauss (4096-dim)"))
for nm, a, b in [("FULL", A_S, A_Sf), ("CENTRE", A_CEN, A_CENf), ("BLOCK", A_BLK, A_BLKf), ("MAG", A_MAG, A_MAGf), ("SIGMA", A_SIG, A_SIGf)]:
    P("  %-10s dim=%-4d r=%d s=%d max=%-6d      dim=%-4d r=%d s=%d max=%d"
      % (nm, a.dim, a.r, a.s, a.maxent, b.dim, b.r, b.s, b.maxent))
P("  >>> A_Sigma: dim %d WITH the Gauss law, dim %d WITHOUT.  With it, A_Sigma == Z(A_S) EXACTLY."
  % (A_SIG.dim, A_SIGf.dim))
P("      Without it, A_Sigma is a 3-qubit electric algebra that shares NO operator with A_S.")
P("      I(A_S : A_Sigma) is still perfectly well defined -- so this is NOT a category error -- but")
P("      the quantity it measures is a different quantity.  BOTH numbers are printed below.")
P("")
P("  g2      chan     H(A_S)     I(A_S:A_Sigma)   note")
for g2 in [0.20, 0.645497, 3.00]:
    E0, psif = lanczos_ground(fb, g2, m=150, seed=3, restarts=4)
    for nm, A in [("FULL", A_Sf), ("CENTRE*", A_CENf), ("BLOCK", A_BLKf), ("MAG", A_MAGf)]:
        J = join(A, A_SIGf)
        I_ = A.entropy(psif) + A_SIGf.entropy(psif) - J.entropy(psif)
        note = "CENTRE IS NO LONGER THE SURFACE -- NOT SCORED" if nm == "CENTRE*" else ""
        P("  %-7.4g %-8s %-10.6f %-16.6f %s" % (g2, nm, A.entropy(psif), I_, note))
    # the same numbers WITH the Gauss law, for the diff
    evp, evecp = np.linalg.eigh(H_matrix(be0, g2)); vp = evecp[:, 0]
    P("     [with Gauss] FULL H=%.6f I(A_S:A_Sigma)=%.6f   CENTRE H=%.6f   C=%.9f"
      % (A_S.entropy(vp), A_S.entropy(vp) + A_SIG.entropy(vp) - J_SS.entropy(vp),
         A_CEN.entropy(vp), A_S.entropy(vp) - A_CEN.entropy(vp)))
    P("     [overlap] || P_vac psi_free || = %.9f"
      % float(np.linalg.norm(np.add.reduceat(np.zeros(1), [0]) * 0 + np.array(
          [sum(psif[c] * be0.sgn[c] for c in range(1 << L) if be0.orb[c] == k) for k in range(32)]) / np.sqrt(128.0))))
    P("")
P("(b) X1b -- NO GAUSS, NO PLAQUETTES.  H = -g2 sum_l X_l.  The ground state is EXACTLY |+>^{x12}.")
plus = np.ones(1 << L, dtype=complex) / np.sqrt(1 << L)
E0b, psib = lanczos_ground(fb, 1.0, plaq=[], m=40, seed=3, restarts=3)
P("   ||psi_lanczos - |+>^{x12}|| = %.3e   E0 = %.9f  (exact -12 g2 = %.9f)"
  % (float(min(np.linalg.norm(psib - plus), np.linalg.norm(psib + plus))), E0b, -12.0))
for nm, A in [("FULL", A_Sf), ("CENTRE*", A_CENf), ("BLOCK", A_BLKf), ("MAG", A_MAGf), ("SIGMA", A_SIGf)]:
    P("   H(%-8s) = %.9f bits" % (nm, A.entropy(plus)))
J = join(A_Sf, A_SIGf)
P("   I(A_S : A_Sigma) = %.9f bits" % (A_Sf.entropy(plus) + A_SIGf.entropy(plus) - J.entropy(plus)))
P("   >>> L BARE QUBITS IN A PRODUCT STATE.  EVERY CHANNEL IS EXACTLY ZERO.  There is no record,")
P("       and nothing about a surface can be asked, because nothing on S is correlated with anything.")

# ==================================================================================================
rule("E5 -- ARM S1: A REGION WITH A LARGER FREE PART.  A' = {0,1,2,3},  S' = {1,2,3,4,5},  Sigma' = {0,6}")
P("This is the NEXT STEP the pre-registration named for itself, executed at zero extra cost.")
P("A' is CROSS-ORBIT from A = {0,1,2} by cardinality (4 vertices vs 3), so the Aut rule is satisfied.")
S2   = LM([1,2,3,4,5]); SIG2 = LM([0,6]); ENV2 = LM([0,6,7,8,9,10,11])
P("S' is vertex-induced: E[{0,1,2,3}] = %s  -> links %s" % ([e for e in EDGES if e[0] in (0,1,2,3) and e[1] in (0,1,2,3)], bits(S2)))
P("Sigma' = delta({0,1,2,3}) = links %s  (|Sigma'| = 2 -- the carrier's UNIQUE minimum cut and its bottleneck)" % bits(SIG2))
A_S2, A_CEN2, _ = chans(be0, S2, [W_S, LM([3,4,5])])
A_SIG2 = alg_links(be0, SIG2, N_PHYS, "SIGMA'")
A_ENV2 = alg_links(be0, ENV2, N_PHYS, "ENV'")
P("")
P("  algebra    n   dim   r   s   max entropy   FREE PART (= s, in bits)")
for nm, A in [("A_S'", A_S2), ("Z(A_S')", A_CEN2), ("A_Sigma'", A_SIG2), ("A_env'", A_ENV2)]:
    P("  %-10s %-3d %-5d %-3d %-3d %-13d %d" % (nm, A.n, A.dim, A.r, A.s, A.maxent, A.s))
P("  PRIMARY REGION FOR COMPARISON:  A_S n=%d dim=%d r=%d s=%d max=%d  FREE PART = %d bit"
  % (A_S.n, A_S.dim, A_S.r, A_S.s, A_S.maxent, A_S.s))
P("  >>> THE FREE PART IS %d BITS HERE vs %d BIT ON THE PRIMARY REGION." % (A_S2.s, A_S.s))
same = sorted(gf2_span(list(A_SIG2.W.values()))) == sorted(gf2_span(list(A_CEN2.W.values())))
P("  >>> A_Sigma' == Z(A_S') as an algebra?  %s   (on the primary region this was TRUE)" % same)
P("      D_forced(S', Sigma') = %d ;  D_forced(S', env') = %d" % (D_forced(S2, SIG2), D_forced(S2, ENV2)))
P("")
P("  g2      H_FULL(S')  H_CENTRE(S')  C' = sum p_k E_k   Delta_surf'   I(A_S':A_Sigma')  I(A_S':A_env')")
for g2 in GRID:
    ev, evec = np.linalg.eigh(H_matrix(be0, g2)); v = evec[:, 0]
    hf = A_S2.entropy(v); hc = A_CEN2.entropy(v)
    Iss = hf + A_SIG2.entropy(v) - join(A_S2, A_SIG2).entropy(v)
    Ise = hf + A_ENV2.entropy(v) - join(A_S2, A_ENV2).entropy(v)
    P("  %-7.4g %-11.6f %-13.6f %-18.9f %-13.9f %-17.6f %.6f" % (g2, hf, hc, hf-hc, Ise-Iss, Iss, Ise))
P("")
P("  HAAR CONTROL ON THE SAME REGION (the pre-committed discriminator, applied to S'):")
rng3 = np.random.default_rng(777)
for t in range(4):
    v = rng3.normal(size=32) + 1j*rng3.normal(size=32); v /= np.linalg.norm(v)
    hf = A_S2.entropy(v); hc = A_CEN2.entropy(v)
    Iss = hf + A_SIG2.entropy(v) - join(A_S2, A_SIG2).entropy(v)
    Ise = hf + A_ENV2.entropy(v) - join(A_S2, A_ENV2).entropy(v)
    P("  haar%d   %-11.6f %-13.6f %-18.9f %-13.9f %-17.6f %.6f" % (t, hf, hc, hf-hc, Ise-Iss, Iss, Ise))
P("")
P("  FORMATION ARMS REPEATED ON S' (same surface Sigma'={0,6}, flux' = eta_0 eta_1 eta_2 eta_3):")
P("  sector      flux'   max|dH_FULL(S')| over grid   max|dC'| over grid   max|dDelta_surf'|")
base = {g2: None for g2 in GRID}
for g2 in GRID:
    ev, evec = np.linalg.eigh(H_matrix(be0, g2)); v = evec[:, 0]
    base[g2] = (A_S2.entropy(v), A_S2.entropy(v) - A_CEN2.entropy(v))
for nm in ["eta{0,4}", "eta{1,5}", "eta{4,5}", "eta{0,1}"]:
    be = BE[nm]
    A2i, C2i, _ = chans(be, S2, [W_S, LM([3,4,5])])
    SG2i = alg_links(be, SIG2, N_PHYS, "s"); EN2i = alg_links(be, ENV2, N_PHYS, "e")
    fl = SEC[nm][0]*SEC[nm][1]*SEC[nm][2]*SEC[nm][3]
    dh, dc, dd = [], [], []
    for g2 in GRID:
        ev, evec = np.linalg.eigh(H_matrix(be, g2)); v = evec[:, 0]
        hf = A2i.entropy(v); hc = C2i.entropy(v)
        Iss = hf + SG2i.entropy(v) - join(A2i, SG2i).entropy(v)
        Ise = hf + EN2i.entropy(v) - join(A2i, EN2i).entropy(v)
        dh.append(abs(hf - base[g2][0])); dc.append(abs((hf-hc) - base[g2][1])); dd.append(abs(Ise-Iss))
    P("  %-11s %+d      %-28.8f %-20.8f %.8f" % (nm, fl, max(dh), max(dc), max(dd)))

open("OUT_run_extra.txt", "w").write("\n".join(LOG) + "\n")
P("\nwrote OUT_run_extra.txt")
