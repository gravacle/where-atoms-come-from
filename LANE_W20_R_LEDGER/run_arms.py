# run_arms.py -- LANE W20_R_LEDGER.  ALL ARMS.  ONE VARIABLE PER ARM, EVERY DIFF PRINTED.
#
# The pre-registration is BINDING and is reproduced in verify_core.py.  Nothing here adjusts a
# criterion, a partition, a fragment, a channel, a coupling grid or a threshold.
#
# ONE ARM IS ADDED BEYOND THE PRE-REGISTRATION AND IT IS FLAGGED EVERYWHERE IT APPEARS:
#   R4  ADIABATIC vs QUENCH.  It is added because the COMMISSIONING TEXT of ARM 3 names
#       "adiabatic versus quench" explicitly as a formation route, and the pre-registration's
#       R1/R2/R3 are charge-PLACEMENT routes only.  It is declared here BEFORE any number is seen:
#           quench: start in the ground state of H(g2 = 5.00), evolve under H(g2_final).
#           report record content at T = 10.0 and the time-average over 200 points in T in [0,50].
#       No other quench parameters will be tried.

import numpy as np, math, sys
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

# ------------------------------------------------------------------ PINNED OBJECTS
S_T    = LM([1, 2, 3])
SIG_T  = LM([0, 4, 5])
ENV_T  = LM([0, 4, 5, 6, 7, 8, 9, 10, 11])
WIDE_T = LM([0, 1, 2, 3, 4, 5])
W_S    = LM([1, 2, 3])
PRIM = [("F1", [0,4,5]), ("F2", [7,8,9]), ("F3", [11]), ("F4", [6,10])]
SECO = [("G1", [0,4,5]), ("G2", [6,9]),   ("G3", [7,11]), ("G4", [8,10])]
GRID = [0.05, 0.10, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00, 1.30, 1.70, 2.20, 3.00, 5.00]
G2_SD = math.sqrt(5.0/12.0)
DELTA = 0.10
GATE  = 0.10

SECTORS = {
    "vacuum":  [1]*8,
    "eta{0,4}":[-1 if v in (0,4) else 1 for v in range(8)],
    "eta{1,5}":[-1 if v in (1,5) else 1 for v in range(8)],
    "eta{4,5}":[-1 if v in (4,5) else 1 for v in range(8)],
    "eta{0,1}":[-1 if v in (0,1) else 1 for v in range(8)],
}
def flux(eta):  # X^Sigma = eta_0 eta_1 eta_2
    return eta[0]*eta[1]*eta[2]

_BE = {}
def backend(name):
    if name not in _BE: _BE[name] = PhysSector(SECTORS[name])
    return _BE[name]

_ALG = {}
def algs(name, wide=False, partition=PRIM):
    key = (name, wide, tuple(k for k,_ in partition))
    if key in _ALG: return _ALG[key]
    be = backend(name)
    if wide:
        chan = {"FULL": alg_links(be, WIDE_T, N_PHYS, "WIDE FULL")}
        A = chan["FULL"]
        cen = [v for v in gf2_span(list(A.W.values())) if all(omega(v,w)==0 for w in A.W.values())]
        chan["CENTRE"] = Algebra(be, cen, N_PHYS, "WIDE CENTRE")
        chan["BLOCK"]  = Algebra(be, [mk(1<<1,0), mk(0,W_S)], N_PHYS, "WIDE-ctx BLOCK")
        chan["MAG"]    = Algebra(be, [mk(0,W_S), mk(0,LM([3,4,5]))], N_PHYS, "WIDE MAG")
    else:
        chan = {
            "FULL":   alg_links(be, S_T, N_PHYS, "FULL"),
            "CENTRE": Algebra(be, [mk(LM([1,2]),0), mk(LM([1,3]),0)], N_PHYS, "CENTRE"),
            "BLOCK":  Algebra(be, [mk(1<<1,0), mk(0,W_S)], N_PHYS, "BLOCK"),
            "MAG":    Algebra(be, [mk(0,W_S)], N_PHYS, "MAG"),
        }
    frags = {k: alg_links(be, LM(ls), N_PHYS, k) for k, ls in partition}
    frags["SIGMA"] = alg_links(be, SIG_T, N_PHYS, "SIGMA")
    frags["ENV"]   = alg_links(be, ENV_T, N_PHYS, "ENV")
    joins = {(c,f): join(chan[c], frags[f]) for c in chan for f in frags}
    _ALG[key] = (chan, frags, joins)
    return _ALG[key]

CH = ["FULL", "CENTRE", "BLOCK", "MAG"]

def report(name, psi, wide=False, partition=PRIM):
    chan, frags, joins = algs(name, wide, partition)
    out = {"H": {}, "I": {}, "R": {}, "DS": {}, "HF": {}}
    for f in frags: out["HF"][f] = frags[f].entropy(psi)
    for c in CH:
        hs = chan[c].entropy(psi)
        out["H"][c] = hs
        for f in frags:
            out["I"][(c,f)] = hs + out["HF"][f] - joins[(c,f)].entropy(psi)
        out["R"][c] = sum(1 for k,_ in partition if out["I"][(c,k)] >= (1-DELTA)*hs) if hs > 1e-12 else 0
        out["DS"][c] = out["I"][(c,"ENV")] - out["I"][(c,"SIGMA")]
    return out

def gs(name, g2, plaq=None):
    H = H_matrix(backend(name), g2, plaq)
    ev, evec = np.linalg.eigh(H)
    return evec[:, 0], ev[1]-ev[0]

def hdr(partition=PRIM):
    return ("  g2      chan   H(A_S)   " + "".join("I:%-7s" % k for k,_ in partition)
            + "I:SIGMA  I:ENV    Delta_surf  Rd  gate")

def line(g2, o, c, partition=PRIM):
    g = "VACUOUS" if o["H"][c] < GATE else "live"
    g2s = g2 if isinstance(g2, str) else "%.4g" % g2
    return ("  %-7s %-6s %-8.5f " % (g2s, c, o["H"][c])
            + "".join("%-9.5f" % o["I"][(c,k)] for k,_ in partition)
            + "%-8.5f %-8.5f %-11.5f %-3d %s" % (o["I"][(c,"SIGMA")], o["I"][(c,"ENV")], o["DS"][c], o["R"][c], g))

# ==================================================================================================
rule("BLOCK 0 -- THE MEASURED CROSSOVER g2_*, AND THE LIVE WINDOW.  FIRST NUMBERS OF THE RUN.")
P("g2_* is MEASURED as argmax over a 200-point log grid of Var_psi(sum_p W_p) in the vacuum sector.")
P("It is never coined.  The a-priori term-balance marker g2_sd = sqrt(5/12) = %.6f is printed beside it." % G2_SD)
be0 = backend("vacuum")
PLAQOP = None
def varW(psi):
    v = np.zeros_like(psi)
    for p in PLAQ: v = v + be0.apply(psi, 0, p)
    m1 = complex(np.vdot(psi, v)).real
    m2 = complex(np.vdot(v, v)).real
    return m2 - m1*m1
fine = np.exp(np.linspace(math.log(0.03), math.log(8.0), 200))
vs, hb, hm = [], [], []
chanV, fragV, joinV = algs("vacuum")
for g2 in fine:
    psi, _ = gs("vacuum", g2)
    vs.append(varW(psi))
    hb.append(chanV["BLOCK"].entropy(psi)); hm.append(chanV["MAG"].entropy(psi))
vs = np.array(vs); hb = np.array(hb); hm = np.array(hm)
istar = int(np.argmax(vs))
G2STAR = float(fine[istar])
P("MEASURED g2_* = %.6f   (Var = %.6f)      g2_sd = %.6f   ratio g2_*/g2_sd = %.4f"
  % (G2STAR, vs[istar], G2_SD, G2STAR/G2_SD))
live = (hb >= GATE) & (hm >= GATE)
if live.any():
    lo, hi = float(fine[live][0]), float(fine[live][-1])
    P("LIVE WINDOW (H_BLOCK >= %.2f AND H_MAG >= %.2f bits): g2 in [%.6f, %.6f]" % (GATE, GATE, lo, hi))
    P("   width in log10(g2) = %.6f decades;  %d of 200 grid points live (%.1f%%)"
      % (math.log10(hi/lo), int(live.sum()), 100.0*live.sum()/200))
    P("   window contains g2_* = %s ;  contains g2_sd = %s" % (lo <= G2STAR <= hi, lo <= G2_SD <= hi))
    P("   max H_BLOCK over the sweep = %.6f bits at g2 = %.6f  (ceiling is 1 bit)"
      % (hb.max(), fine[int(np.argmax(hb))]))
    P("   max H_MAG   over the sweep = %.6f bits at g2 = %.6f  (ceiling is 1 bit)"
      % (hm.max(), fine[int(np.argmax(hm))]))
else:
    lo = hi = None
    P("LIVE WINDOW IS EMPTY.  Every point of the 200-grid fails the AXIS 6 gate in BLOCK or MAG.")
    P("   max H_BLOCK = %.6f at g2=%.6f ; max H_MAG = %.6f at g2=%.6f"
      % (hb.max(), fine[int(np.argmax(hb))], hm.max(), fine[int(np.argmax(hm))]))
P("")
P("  %-10s %-10s %-10s %-10s %-10s" % ("g2", "Var(sumW)", "H_BLOCK", "H_MAG", "H_FULL"))
for g2 in GRID:
    psi, _ = gs("vacuum", g2)
    P("  %-10.4g %-10.5f %-10.6f %-10.6f %-10.6f"
      % (g2, varW(psi), chanV["BLOCK"].entropy(psi), chanV["MAG"].entropy(psi), chanV["FULL"].entropy(psi)))

# ==================================================================================================
rule("ARM A0 -- REFERENCE.  vacuum sector, primary partition, primary PLAQ, 13-point grid.")
P("D(F) for the primary partition (pre-registration's vacuity numbers, recomputed): F1=%d F2=%d F3=%d F4=%d ; SIGMA=%d ENV=%d"
  % tuple([D_forced(S_T, LM(ls)) for _,ls in PRIM] + [D_forced(S_T,SIG_T), D_forced(S_T,ENV_T)]))
P("STRUCTURAL FACT MEASURED IN verify_core V8/V11 AND USED THROUGHOUT:")
P("   A_Sigma == Z(A_S) as an ALGEBRA (identical label spaces), hence A_Sigma is a SUBALGEBRA OF A_S,")
P("   hence join(A_S, A_Sigma) = A_S and  I(A_S : A_Sigma) = H(Z(A_S))  IDENTICALLY.  The 'surface")
P("   fragment' column below is therefore a THEOREM in the FULL and CENTRE channels, not a measurement.")
P(hdr())
A0 = {}
for g2 in GRID:
    psi, gap = gs("vacuum", g2)
    A0[g2] = report("vacuum", psi)
    for c in CH: P(line(g2, A0[g2], c))
    P("")

# ==================================================================================================
rule("ARM A1 -- CHARGE SECTOR.  eta = -1 at {0,4}.  ISOLATED VARIABLE: the sign pattern, nothing else.")
P("flux(Sigma) = %+d in A0(vacuum)   ,   %+d in A1(eta{0,4})" % (flux(SECTORS["vacuum"]), flux(SECTORS["eta{0,4}"])))
P(hdr())
A1 = {}
for g2 in GRID:
    psi, _ = gs("eta{0,4}", g2)
    A1[g2] = report("eta{0,4}", psi)
    for c in CH: P(line(g2, A1[g2], c))
    P("")
rule("ARM A1 DIFF -- A1 minus A0.  PRINTED, NOT SUMMARISED.")
P("  g2      chan   dH(A_S)    dI:F1      dI:F2      dI:F3      dI:F4      dI:SIGMA   dI:ENV     dDelta_surf")
for g2 in GRID:
    for c in CH:
        d = lambda k: A1[g2]["I"][(c,k)] - A0[g2]["I"][(c,k)]
        P("  %-7.4g %-6s %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f"
          % (g2, c, A1[g2]["H"][c]-A0[g2]["H"][c], d("F1"), d("F2"), d("F3"), d("F4"),
             d("SIGMA"), d("ENV"), A1[g2]["DS"][c]-A0[g2]["DS"][c]))
    P("")

# ==================================================================================================
rule("ARM A2 -- BOUNDARY ALGEBRA.  WIDE (S u Sigma) vs NARROW (=CENTER on the physical sector).")
P("ISOLATED VARIABLE: the algebra assigned to S.  State, partition, coupling, sector all fixed.")
chW, frW, joW = algs("vacuum", wide=True)
P("WIDE FULL   : n=%d dim=%d r=%d s=%d maxent=%d bits" % (chW["FULL"].n, chW["FULL"].dim, chW["FULL"].r, chW["FULL"].s, chW["FULL"].maxent))
P("NARROW FULL : n=%d dim=%d r=%d s=%d maxent=%d bits" % (chanV["FULL"].n, chanV["FULL"].dim, chanV["FULL"].r, chanV["FULL"].s, chanV["FULL"].maxent))
P("  g2      H_FULL(WIDE)  H_FULL(NARROW)  DISAGREEMENT   H_MAG(WIDE)  H_MAG(NARROW)  dMAG   DS(W)-DS(N) FULL")
A2 = {}
for g2 in GRID:
    psi, _ = gs("vacuum", g2)
    w = report("vacuum", psi, wide=True); n = A0[g2]
    A2[g2] = w
    P("  %-7.4g %-13.6f %-15.6f %-14.6f %-12.6f %-14.6f %-6.6f %-.6f"
      % (g2, w["H"]["FULL"], n["H"]["FULL"], w["H"]["FULL"]-n["H"]["FULL"],
         w["H"]["MAG"], n["H"]["MAG"], w["H"]["MAG"]-n["H"]["MAG"], w["DS"]["FULL"]-n["DS"]["FULL"]))

# ==================================================================================================
rule("ARM 2 OF THE COMMISSION -- REMOVE THE COUPLING.  Fix g2 at ONE value, delete the sweep.")
P("Does the boundary-algebra disagreement stop varying, or was the phase-indexing an artefact?")
P("The disagreement above is a FUNCTION of g2.  Deleting the sweep does not change any number; it")
P("deletes the ability to SEE the number vary.  Reported as the range of the disagreement:")
dis = [A2[g]["H"]["FULL"] - A0[g]["H"]["FULL"] for g in GRID]
dmag = [A2[g]["H"]["MAG"] - A0[g]["H"]["MAG"] for g in GRID]
P("   WIDE-NARROW FULL disagreement over the grid: min = %.6f at g2=%.4g ; max = %.6f at g2=%.4g ; range = %.6f bits"
  % (min(dis), GRID[int(np.argmin(dis))], max(dis), GRID[int(np.argmax(dis))], max(dis)-min(dis)))
P("   WIDE-NARROW MAG  disagreement over the grid: min = %.6f ; max = %.6f ; range = %.6f bits"
  % (min(dmag), max(dmag), max(dmag)-min(dmag)))
P("   at the single fixed point g2 = g2_sd = %.6f the disagreement is a CONSTANT:" % G2_SD)
psi_sd, _ = gs("vacuum", G2_SD)
wsd = report("vacuum", psi_sd, wide=True); nsd = report("vacuum", psi_sd)
P("       H_FULL(WIDE) = %.6f   H_FULL(NARROW) = %.6f   disagreement = %.6f bits"
  % (wsd["H"]["FULL"], nsd["H"]["FULL"], wsd["H"]["FULL"]-nsd["H"]["FULL"]))

# ==================================================================================================
rule("ARM 3 OF THE COMMISSION -- VARY THE BOUNDARY FORMATION.  R1 / R2 / R3, THE PRE-REGISTERED PAIRS.")
P("The boundary's EXISTENCE is held byte-identical in every pair: Sigma = {0,4,5}, same links, same")
P("algebra, same partition, same Hamiltonian, same coupling grid.  Only the ROUTE moves.")
P("")
FORM = [("R1 ALPHA", "eta{0,4}", "eta{1,5}"), ("R2 BETA", "vacuum", "eta{4,5}"), ("R3 GAMMA", "vacuum", "eta{0,1}")]
STATES = {}
for nm in SECTORS:
    STATES[nm] = {}
    for g2 in GRID:
        psi, _ = gs(nm, g2)
        STATES[nm][g2] = report(nm, psi)
for tag, a, b in FORM:
    P("-"*104)
    P("%s :  %s  vs  %s      flux(Sigma) = %+d  vs  %+d   -> IDENTICAL: %s"
      % (tag, a, b, flux(SECTORS[a]), flux(SECTORS[b]), flux(SECTORS[a])==flux(SECTORS[b])))
    P("  g2      chan   H(A_S)_a   H(A_S)_b   dH        I:ENV_a    I:ENV_b    dI:ENV     DS_a      DS_b      dDS")
    for g2 in GRID:
        for c in CH:
            oa, ob = STATES[a][g2], STATES[b][g2]
            P("  %-7.4g %-6s %-10.6f %-10.6f %-9.6f %-10.6f %-10.6f %-10.6f %-9.6f %-9.6f %-9.6f"
              % (g2, c, oa["H"][c], ob["H"][c], ob["H"][c]-oa["H"][c],
                 oa["I"][(c,"ENV")], ob["I"][(c,"ENV")], ob["I"][(c,"ENV")]-oa["I"][(c,"ENV")],
                 oa["DS"][c], ob["DS"][c], ob["DS"][c]-oa["DS"][c]))
        P("")

rule("ARM 3 -- MAXIMUM FORMATION DIFFERENCE PER CHANNEL, OVER THE WHOLE GRID.  THE HEADLINE NUMBERS.")
P("  pair       chan   max|dH(A_S)|  max|dI:ENV|  max|dI:F3|   max|dDelta_surf|   both-live-at-g2")
for tag, a, b in FORM:
    for c in CH:
        dh = [abs(STATES[b][g]["H"][c]-STATES[a][g]["H"][c]) for g in GRID]
        de = [abs(STATES[b][g]["I"][(c,"ENV")]-STATES[a][g]["I"][(c,"ENV")]) for g in GRID]
        d3 = [abs(STATES[b][g]["I"][(c,"F3")]-STATES[a][g]["I"][(c,"F3")]) for g in GRID]
        dd = [abs(STATES[b][g]["DS"][c]-STATES[a][g]["DS"][c]) for g in GRID]
        lg = [g for g in GRID if STATES[a][g]["H"][c] >= GATE and STATES[b][g]["H"][c] >= GATE]
        P("  %-10s %-6s %-13.8f %-12.8f %-12.8f %-18.8f %d/13" % (tag, c, max(dh), max(de), max(d3), max(dd), len(lg)))

# ==================================================================================================
rule("ARM R4 -- ADIABATIC vs QUENCH.  *** NOT IN THE PRE-REGISTRATION *** -- added because ARM 3 of")
P("the commissioning text names it.  Parameters declared before the numbers: start = ground state of")
P("H(g2=5.00) [deep electric], evolve under H(g2_final).  Report T = 10.0 and the average over 200")
P("points T in [0,50].  Same sector (vacuum), same Sigma, same flux(+1), same partition, same algebra.")
P("  g2_f    chan   H_gs       H_T=10     H_Tavg     dH(T=10)   DS_gs     DS_T=10   dDS(T=10)  DS_Tavg")
psi0, _ = gs("vacuum", 5.00)
Tgrid = np.linspace(0.0, 50.0, 200)
R4 = {}
for g2 in GRID:
    H = H_matrix(be0, g2)
    ev, evec = np.linalg.eigh(H)
    c0 = evec.conj().T @ psi0
    psiT = evec @ (np.exp(-1j*ev*10.0) * c0)
    oT = report("vacuum", psiT)
    accH = {c: 0.0 for c in CH}; accD = {c: 0.0 for c in CH}
    for T in Tgrid:
        pv = evec @ (np.exp(-1j*ev*T) * c0)
        o = report("vacuum", pv)
        for c in CH: accH[c] += o["H"][c]/len(Tgrid); accD[c] += o["DS"][c]/len(Tgrid)
    R4[g2] = (oT, accH, accD)
    for c in CH:
        P("  %-7.4g %-6s %-10.6f %-10.6f %-10.6f %-10.6f %-9.6f %-9.6f %-10.6f %-9.6f"
          % (g2, c, A0[g2]["H"][c], oT["H"][c], accH[c], oT["H"][c]-A0[g2]["H"][c],
             A0[g2]["DS"][c], oT["DS"][c], oT["DS"][c]-A0[g2]["DS"][c], accD[c]))
    P("")

# ==================================================================================================
rule("ARM D1 / ARM 4 OF THE COMMISSION -- REMOVE THE DYNAMICS.")
P("ground state  vs  Haar-random physical state  vs  product states.  Same sector, same partition,")
P("same algebra.  ISOLATED VARIABLE: the existence of a state-selecting Hamiltonian.")
P("PRODUCT STATES, declared: PROD-Z = |0>^{x12} projected into the vacuum sector (a single coset")
P("basis vector); PROD-X = |+>^{x12} projected (the uniform superposition, = the g2 -> inf ground state).")
rng = np.random.default_rng(20260817)
HAAR = []
for t in range(4):
    v = rng.normal(size=32) + 1j*rng.normal(size=32); v /= np.linalg.norm(v)
    HAAR.append(v)
pz = np.zeros(32, dtype=complex); pz[be0.orb[0]] = 1.0
px = np.ones(32, dtype=complex)/np.sqrt(32.0)
P(hdr())
NOSTATE = {}
for nm, v in [("HAAR0", HAAR[0]), ("HAAR1", HAAR[1]), ("HAAR2", HAAR[2]), ("HAAR3", HAAR[3]),
              ("PROD-Z", pz), ("PROD-X", px)]:
    o = report("vacuum", v); NOSTATE[nm] = o
    for c in CH: P(line(nm, o, c))
    P("")
rule("ARM D1 DIFF -- HAAR MEAN minus GROUND STATE AT EACH COUPLING.  THE STANDING RULE'S NUMBER.")
P("STANDING RULE: a plateau a Haar state also produces is NOT reported as a result.")
P("  g2      chan   H_gs      H_haar    dH        DS_gs     DS_haar   dDS       I:F3_gs   I:F3_haar dI:F3")
hm_ = {c: float(np.mean([NOSTATE["HAAR%d"%t]["H"][c] for t in range(4)])) for c in CH}
hd_ = {c: float(np.mean([NOSTATE["HAAR%d"%t]["DS"][c] for t in range(4)])) for c in CH}
h3_ = {c: float(np.mean([NOSTATE["HAAR%d"%t]["I"][(c,"F3")] for t in range(4)])) for c in CH}
for g2 in GRID:
    for c in CH:
        P("  %-7.4g %-6s %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f %-9.6f"
          % (g2, c, A0[g2]["H"][c], hm_[c], hm_[c]-A0[g2]["H"][c],
             A0[g2]["DS"][c], hd_[c], hd_[c]-A0[g2]["DS"][c],
             A0[g2]["I"][(c,"F3")], h3_[c], h3_[c]-A0[g2]["I"][(c,"F3")]))
    P("")

# ==================================================================================================
rule("ARM P1 -- PARTITION.  secondary partition G1..G4 vs primary F1..F4.  ISOLATED: the partition.")
P("D(G) recomputed: G1=%d G2=%d G3=%d G4=%d" % tuple(D_forced(S_T, LM(ls)) for _,ls in SECO))
P(hdr(SECO))
P1 = {}
for g2 in GRID:
    psi, _ = gs("vacuum", g2)
    P1[g2] = report("vacuum", psi, partition=SECO)
    for c in CH: P(line(g2, P1[g2], c, SECO))
    P("")
rule("ARM P1 DIFF -- R_delta AND Delta_surf UNDER THE TWO PARTITIONS.")
P("  g2      chan   R_prim  R_seco  dR   DS_prim   DS_seco   dDS   (Delta_surf must be partition-INDEPENDENT)")
for g2 in GRID:
    for c in CH:
        P("  %-7.4g %-6s %-7d %-7d %-4d %-9.6f %-9.6f %-.2e"
          % (g2, c, A0[g2]["R"][c], P1[g2]["R"][c], P1[g2]["R"][c]-A0[g2]["R"][c],
             A0[g2]["DS"][c], P1[g2]["DS"][c], abs(P1[g2]["DS"][c]-A0[g2]["DS"][c])))
    P("")

open("OUT_run_arms.txt", "w").write("\n".join(LOG) + "\n")
P("\nwrote OUT_run_arms.txt")
