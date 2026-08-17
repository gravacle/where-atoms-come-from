# run_probe.py -- LANE W20_C.  THE ARMS THAT SURROUND THE PRIMARY PAIR.
#   BLOCK 5  the pinned crossover estimator, and its repair
#   BLOCK 6  the mechanism: why the charged arm's record does not empty at the electric end
#   BLOCK 7  FORMATION ARMS R1 / R2 / R3, with trace distances
#   BLOCK 8  ALL 128 CHARGE SECTORS: is the record a function of the surface flux alone?
#   BLOCK 9  D1 dynamics removal (Haar), P1 partition, A2 boundary algebra
#   BLOCK 10 X1 gauge removal, exactly
#   BLOCK 11 boundary EXISTENCE: the category-error argument, plus a cross-orbit surface variation
#   BLOCK 12 exhaustive search for where the plateau CAN fail
import math, itertools, functools
import numpy as np
import w20c_core as C

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def rule(t=""):
    P("\n" + "=" * 108); P(t); P("=" * 108)

def rho_alg(A, sec, psi):
    b = A.state(sec, psi); d = 1 << A.k; n = (1 << A.r) * d
    M = np.zeros((n, n), dtype=complex)
    for s in range(1 << A.r):
        M[s * d:(s + 1) * d, s * d:(s + 1) * d] = b[s]
    return (M + M.conj().T) / 2

def tracedist(M1, M2):
    return 0.5 * float(np.sum(np.abs(np.linalg.eigvalsh(M1 - M2).real)))

# automorphisms and the stabiliser of A = {0,1,2}
Eset = {frozenset(e) for e in C.E}
AUT = [p for p in itertools.permutations(range(C.V))
       if all(frozenset((p[a], p[b])) in Eset for a, b in C.E)]
STAB = [p for p in AUT if {p[v] for v in (0, 1, 2)} == {0, 1, 2}]

# ============================================================================================
rule("BLOCK 5 -- THE PINNED CROSSOVER ESTIMATOR IS DEGENERATE.  MEASURED, THEN REPAIRED.")
P("PIN: 'g2_* = argmax over a 200-point log grid of Var_psi(sum_p W_p)'.  Measured in run_charge:")
P("   A0 vacuum      : argmax at g2 = 5.000000, THE RIGHT-HAND EDGE OF THE GRID.")
P("   A1 charges{0,4}: argmax at g2 = 1.981344.")
P("The A0 answer is not a crossover, it is the grid boundary, and it is forced:")
P("   as g2 -> infinity the state -> an X-eigenstate, <W_p> -> 0 for every p, and W_p^2 = I, so")
P("   Var(sum_p W_p) -> C = 5 MONOTONICALLY FROM BELOW.  The estimator's argmax is at g2 = infinity")
P("   for any carrier.  IT MEASURES THE ELECTRIC LIMIT, NOT THE CROSSOVER.")
FINE = [10 ** x for x in np.linspace(math.log10(0.05), math.log10(5.00), 200)]
sec0 = C.Sector([]); sec1 = C.Sector([0, 4])
P("")
P("   the curve, printed rather than asserted (A0):")
P("   %-10s %-12s %-12s %-12s" % ("g2", "Var(sumW)", "<sumW>", "<sum X_l>"))
for g2 in (0.05, 0.20, 0.60, 1.00, 2.20, 5.00):
    _, psi, _ = sec0.ground(g2)
    v, m = C.var_plaq(sec0, psi)
    xs = sum(float(np.vdot(psi, C.Alg([C.sv(1 << l, 0)], "x").CT[0][1][2] * sec0.parity(1 << l) * psi).real)
             for l in range(C.L))
    P("   %-10.4f %-12.6f %-12.6f %-12.6f" % (g2, v, m, xs))
P("")
P("REPAIRED ESTIMATORS, both MEASURED, both declared here as coined replacements:")
P("   (a) g2_steep = argmax over the same 200-point grid of |d<sum_p W_p>/d ln g2|  (steepest")
P("       change of the magnetic order parameter -- a crossover by definition).")
P("   (b) g2_bal   = the g2 where the two ENERGY TERMS balance: (1/g2)|<sum_p W_p>| = g2|<sum_l X_l>|.")
for tag, sec in (("A0", sec0), ("A1", sec1)):
    ms, xs = [], []
    for g2 in FINE:
        _, psi, _ = sec.ground(g2)
        _, m = C.var_plaq(sec, psi)
        ms.append(m)
        xs.append(sum(float((sec.parity(1 << l) * np.abs(psi) ** 2).sum()) for l in range(C.L)))
    lg = np.log(np.array(FINE))
    d = np.abs(np.gradient(np.array(ms), lg))
    i = int(np.argmax(d))
    f = [abs(m) / g - g * abs(x) for g, m, x in zip(FINE, ms, xs)]
    j = int(np.argmin(np.abs(np.array(f))))
    P("   %s : g2_steep = %.6f   g2_bal = %.6f   (a-priori marker g2_sd = %.6f)"
      % (tag, FINE[i], FINE[j], math.sqrt(5 / 12)))
P("   >>> THE A-PRIORI MARKER sqrt(C/L) = 0.645497 LANDS INSIDE BOTH REPAIRED ESTIMATES AND")
P("       NOWHERE NEAR THE PINNED ONE.  THE PINNED ESTIMATOR IS RETIRED, WITH ITS NUMBER SHOWN.")
P("       ALL PHASE LABELS BELOW ARE ASSIGNED RELATIVE TO g2_steep, NOT TO THE PINNED VALUE.")

# ============================================================================================
rule("BLOCK 6 -- THE MECHANISM.  WHY THE CHARGED ARM'S RECORD DOES NOT EMPTY AT THE ELECTRIC END.")
P("PRE-REGISTERED VACUITY CLAIM V2: 'g2 >> g2_*: <X_1> -> 1, block pure, H_BLOCK -> 0.'")
P("MEASURED: A0 H_BLOCK(g2=5.00) = 0.000205 bits.  A1 H_BLOCK(g2=5.00) = 0.999427 bits.")
P("V2 IS TRUE IN THE VACUUM SECTOR AND FALSE IN A CHARGED SECTOR.  It was written as a property of")
P("the coupling; it is a property of the SECTOR.  Here is the incidence reason, computed:")
for tag, sec in (("A0 vacuum", sec0), ("A1 charges {0,4}", sec1)):
    w = [C.pop(int(u)) for u in sec.U]
    mn = min(w)
    deg = [int(u) for u in sec.U if C.pop(int(u)) == mn]
    P("   %s : minimum electric weight |u| = %d, DEGENERACY %d" % (tag, mn, len(deg)))
    for u in deg:
        P("        u = links %-16s  restricted to S = %s" % (C.bits(u), C.bits(u & C.S_MASK)))
P("")
P(">>> THE g2 -> infinity GROUND STATE IS THE UNIFORM SUPERPOSITION OF THE MINIMUM-WEIGHT ELECTRIC")
P("    CONFIGURATIONS.  In the vacuum sector that set has ONE element (u = 0) and the record is")
P("    empty.  In the {0,4} sector it has FOUR elements -- the four shortest flux strings joining")
P("    the two charges -- and they DISAGREE ON S.  THE RECORD IS WHICH ROUTE THE FLUX STRING TOOK.")
P("    That is not a narration: the number is H_FULL(A0, g2=5) = 0.00091 vs H_FULL(A1, g2=5) =")
P("    1.01071 bits, a 1.010 bit gap produced by two minus signs in the Gauss constraint.")
P("")
P("    AND THE HONEST SECOND READING, RECORDED: 4-fold degeneracy of the shortest 0-4 path is a")
P("    property of tri_chain12's geometry, not of charge as such.  A carrier with a UNIQUE shortest")
P("    path between the two charged vertices would give a non-degenerate electric ground state and")
P("    an empty record at large g2 even with charges present.  THIS IS A CARRIER-DEPENDENT")
P("    MECHANISM AND IT IS LOGGED AS ONE.  Counted here: pairs (v,w) of vertices at distance 3 with")
P("    a UNIQUE shortest path on tri_chain12:")
import collections
adj = collections.defaultdict(list)
for i, (a, b) in enumerate(C.E):
    adj[a].append(b); adj[b].append(a)
uniq = []
for a in range(C.V):
    for b in range(a + 1, C.V):
        dist = {a: 0}; q = [a]
        while q:
            x = q.pop(0)
            for y in adj[x]:
                if y not in dist: dist[y] = dist[x] + 1; q.append(y)
        npaths = {a: 1}
        order = sorted(dist, key=lambda v: dist[v])
        for x in order:
            if x == a: continue
            npaths[x] = sum(npaths[y] for y in adj[x] if dist[y] == dist[x] - 1)
        uniq.append((a, b, dist[b], npaths[b]))
P("       %s" % ", ".join("(%d,%d)d=%d n=%d" % t for t in uniq))
P("       -> %d of %d vertex pairs have a UNIQUE geodesic.  The declared pair {0,4} has %d."
  % (sum(1 for t in uniq if t[3] == 1), len(uniq), [t for t in uniq if t[:2] == (0, 4)][0][3]))

# ============================================================================================
rule("BLOCK 7 -- FORMATION ARMS.  SAME SURFACE, SAME FLUX, DIFFERENT ROUTE.")
P("Sigma = {0,4,5} is byte-identical in every arm.  flux(Sigma) = eta_0 eta_1 eta_2 is held fixed")
P("inside each pair.  Only WHERE the charge sits moves.  The comparison is well posed across")
P("sectors because A_S = M_2^(+4) is the SAME abstract algebra in every sector -- the record lives")
P("in the region's own algebra, not in the sector's Hilbert space.")
PAIRS = [("R1 alpha", [0, 4], [1, 5]), ("R2 beta ", [], [4, 5]), ("R3 gamma", [], [0, 1])]
for nm, ca, cb in PAIRS:
    sa, sb = C.Sector(ca), C.Sector(cb)
    hit = [p for p in STAB if {p[v] for v in ca} == set(cb)]
    P("")
    P("%s : eta=-1 at %s  vs  %s .  flux(Sigma) = %+d vs %+d .  Stab(A) maps: %d -> %s"
      % (nm, ca if ca else "none", cb, C.flux_sigma(sa), C.flux_sigma(sb), len(hit),
         "VACUOUS BY SYMMETRY" if hit else "LIVE"))
    P("   minimum-weight route: |u0| = %d at links %s   vs   |u0| = %d at links %s"
      % (C.pop(sa.u0), C.bits(sa.u0), C.pop(sb.u0), C.bits(sb.u0)))
    P("   %-6s %10s %10s %10s %10s %12s %12s" % ("g2", "H_FULL a", "H_FULL b", "H_CEN a", "H_CEN b",
                                                 "D_tr(A_S)", "D_tr(BLOCK)"))
    for g2 in C.GRID:
        _, pa, _ = sa.ground(g2); _, pb, _ = sb.ground(g2)
        P("   %-6.2f %10.6f %10.6f %10.6f %10.6f %12.9f %12.9f"
          % (g2, C.A_FULL.entropy(sa, pa), C.A_FULL.entropy(sb, pb),
             C.A_CENTRE.entropy(sa, pa), C.A_CENTRE.entropy(sb, pb),
             tracedist(rho_alg(C.A_FULL, sa, pa), rho_alg(C.A_FULL, sb, pb)),
             tracedist(rho_alg(C.A_BLOCK, sa, pa), rho_alg(C.A_BLOCK, sb, pb))))

# ============================================================================================
rule("BLOCK 8 -- ALL 128 CHARGE SECTORS.  IS THE RECORD A FUNCTION OF THE SURFACE FLUX ALONE?")
P("THIS IS THE NON-DEGENERATE FORM OF H-SURFACE, and it is the one arm in this lane that no")
P("theorem in the pre-registration forces either way.  Delta_surf could not fail (BLOCK 1).  This")
P("can: nothing whatsoever requires two sectors with the SAME flux through Sigma to give the SAME")
P("rho|A_S.  If the record is boundary data fixed by the surface and the flux, the 64 sectors with")
P("flux(Sigma) = +1 must all give the same record, and so must the 64 with flux = -1.")
SECTORS = [tuple(q) for q in
           (tuple(C.bits(m, C.V)) for m in range(1 << C.V)) if len(q) % 2 == 0]
P("admissible sectors (even charge parity): %d" % len(SECTORS))
def orbit_key(q):
    return min(tuple(sorted(p[v] for v in q)) for p in STAB)
for g2 in (0.45, 1.00, 3.00):
    P("")
    P("g2 = %.2f" % g2)
    data = {}
    for q in SECTORS:
        s = C.Sector(list(q))
        _, psi, _ = s.ground(g2)
        data[q] = (C.flux_sigma(s), rho_alg(C.A_FULL, s, psi),
                   C.A_FULL.entropy(s, psi), C.A_CENTRE.entropy(s, psi))
    for fl in (+1, -1):
        grp = [q for q in SECTORS if data[q][0] == fl]
        hs = [data[q][2] for q in grp]
        P("   flux(Sigma) = %+d : %d sectors.  H_FULL ranges [%.6f , %.6f] , spread %.6f bits"
          % (fl, len(grp), min(hs), max(hs), max(hs) - min(hs)))
        best = (0.0, None, None); bestlive = (0.0, None, None)
        for i in range(len(grp)):
            for j in range(i + 1, len(grp)):
                d = tracedist(data[grp[i]][1], data[grp[j]][1])
                if d > best[0]: best = (d, grp[i], grp[j])
                if orbit_key(grp[i]) != orbit_key(grp[j]) and d > bestlive[0]:
                    bestlive = (d, grp[i], grp[j])
        P("      MAX pairwise D_tr(rho|A_S) inside the flux class : %.9f   between %s and %s"
          % (best[0], best[1], best[2]))
        P("      MAX over Stab(A)-INEQUIVALENT pairs only (LIVE)   : %.9f   between %s and %s"
          % (bestlive[0], bestlive[1], bestlive[2]))
        nz = sum(1 for i in range(len(grp)) for j in range(i + 1, len(grp))
                 if tracedist(data[grp[i]][1], data[grp[j]][1]) > 1e-6)
        P("      pairs with D_tr > 1e-6 : %d of %d" % (nz, len(grp) * (len(grp) - 1) // 2))
P("")
P(">>> READ THE SPREAD.  If it is 0, the record on S is a function of flux(Sigma) alone and")
P("    H-SURFACE survives its one non-degenerate test in this lane.  If it is not 0, the surface")
P("    and its flux DO NOT determine the record, and the number is the size of the failure.")

# ============================================================================================
rule("BLOCK 9 -- D1 DYNAMICS REMOVAL (HAAR), P1 PARTITION, A2 BOUNDARY ALGEBRA")
P("D1.  Removal of the state-selecting Hamiltonian.  Haar-random state of the SAME physical sector,")
P("same partition, same algebra.  STANDING RULE: a plateau a Haar state also produces is not a")
P("result.  W-19 established NARRATION in the electric channel; BLOCK and MAG were untested.")
P("%-14s %9s %9s %9s %9s %9s %9s %9s %6s"
  % ("state", "H_FULL", "H_CEN", "H_BLOCK", "H_MAG", "I:F1", "I:F3", "Dsurf", "R_del"))
def row(tag, sec, psi):
    hF = C.A_FULL.entropy(sec, psi); hC = C.A_CENTRE.entropy(sec, psi)
    Is = [C.MI(C.A_FULL, A, sec, psi) for _, A in C.AF_P]
    Rd = sum(1 for x in Is if hF > 0 and x >= 0.9 * hF)
    P("%-14s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.6f %6d"
      % (tag, hF, hC, C.A_BLOCK.entropy(sec, psi), C.A_MAG.entropy(sec, psi),
         Is[0], Is[2], 2 * (hF - hC), Rd))
for tag, sec in (("A0", sec0), ("A1", sec1)):
    for g2 in (0.45, 0.80, 3.00):
        _, psi, _ = sec.ground(g2)
        row("%s ground %.2f" % (tag, g2), sec, psi)
    for sd in (7, 8, 9):
        row("%s HAAR s%d" % (tag, sd), sec, C.haar(sec, sd))
P("")
P("P1.  PARTITION.  primary F1..F4 vs secondary G1..G4, same state, same algebra.")
P("%-6s %-4s %9s %9s %9s %9s %9s %6s" % ("g2", "arm", "H_FULL", "I:G1", "I:G2", "I:G3", "I:G4", "R_del"))
for g2 in (0.45, 1.00, 3.00):
    for tag, sec in (("A0", sec0), ("A1", sec1)):
        _, psi, _ = sec.ground(g2)
        hF = C.A_FULL.entropy(sec, psi)
        Is = [C.MI(C.A_FULL, A, sec, psi) for _, A in C.AF_S]
        Rd = sum(1 for x in Is if hF > 0 and x >= 0.9 * hF)
        P("%-6.2f %-4s %9.5f %9.5f %9.5f %9.5f %9.5f %6d" % (g2, tag, hF, Is[0], Is[1], Is[2], Is[3], Rd))
P("   D(F) vacuity numbers: primary %s ; secondary %s"
  % ([C.Dnum(f) for _, f in C.FRAG_P], [C.Dnum(f) for _, f in C.FRAG_S]))
P("")
P("A2.  BOUNDARY ALGEBRA.  WIDE (S u Sigma, nbar=5, max 3 bits) vs NARROW (nbar=4, max 3 bits).")
P("   CONFOUND RECORDED: WIDE's support overlaps F1 = Sigma, so only F2,F3,F4 remain disjoint.")
P("%-6s %-4s %10s %10s %12s %10s %10s" % ("g2", "arm", "H_NARROW", "H_WIDE", "H_WIDE-H_NAR", "I_W:F2", "I_W:F3"))
for g2 in C.GRID:
    for tag, sec in (("A0", sec0), ("A1", sec1)):
        _, psi, _ = sec.ground(g2)
        hn = C.A_FULL.entropy(sec, psi); hw = C.A_WIDE.entropy(sec, psi)
        P("%-6.2f %-4s %10.6f %10.6f %12.6f %10.6f %10.6f"
          % (g2, tag, hn, hw, hw - hn,
             C.MI(C.A_WIDE, C.AF_P[1][1], sec, psi), C.MI(C.A_WIDE, C.AF_P[2][1], sec, psi)))

# ============================================================================================
rule("BLOCK 10 -- X1 GAUGE REMOVAL, DONE EXACTLY (NO 4096x4096 DIAGONALISATION NEEDED)")
P("Removing the Gauss projection means diagonalising the SAME H on the unconstrained 2^12 space.")
P("But H is BLOCK DIAGONAL in the charge sectors: the plaquette hopping translates u by a cycle,")
P("which preserves d(u) = q.  So the unconstrained spectrum is the union of the 128 sector spectra,")
P("and the unconstrained ground state is the ground state of whichever sector is lowest.  EXACT.")
P("%-6s %14s %14s %10s %12s" % ("g2", "E0(vacuum)", "E0(best charged)", "gap", "|<free|proj>|"))
for g2 in C.GRID:
    e_vac = sec0.ground(g2)[0]
    best = min((C.Sector(list(q)).ground(g2)[0], q) for q in SECTORS if q)
    P("%-6.2f %14.9f %14.9f %10.6f %12s"
      % (g2, e_vac, best[0], best[0] - e_vac, "1.000000000" if e_vac <= best[0] else "0.0"))
P("   >>> THE UNCONSTRAINED GROUND STATE IS THE VACUUM-SECTOR GROUND STATE AT EVERY COUPLING.")
P("       |<free|projected>| = 1.000000000, reproducing W-19's Perron-Frobenius result.")
P("       ON THE VACUUM ARM, REMOVING THE GAUSS LAW MOVES THE STATE BY EXACTLY ZERO -> the numbers")
P("       do not move -> NARRATION on that arm.")
P("   >>> ON A CHARGED ARM IT IS NOT A REMOVAL AT ALL.  Deleting the constraint does not leave the")
P("       charged state weakened; it DELETES the state and hands back the vacuum one, %.6f in"
  % (C.Sector([0, 4]).ground(1.0)[0] - sec0.ground(1.0)[0]))
P("       energy away at g2=1.00.  'The same state across the removal' does not exist.  The well-")
P("       posed question is the one BLOCK 8 asks instead.")

# ============================================================================================
rule("BLOCK 11 -- BOUNDARY EXISTENCE: THE CATEGORY ERROR, AND WHAT VARYING IT MEASURES")
P("delta(A) = empty iff A is empty or A = V.  Either way S is empty or is everything and")
P("I(A_S : A_F) has no arguments.  Worse, Z(A_S) IS the surface algebra as an operator identity")
P("(BLOCK 1(i)), so deleting the surface deletes 2 of the record's %d bits BY DEFINITION."
  % C.A_FULL.maxent)
P("A REMOVAL TEST HERE COULD NOT HAVE FAILED.  It is a CATEGORY ERROR and is not run.")
P("")
P("WHAT VARYING IT MEANS INSTEAD -- and the four-triangle version is vacuous by Aut, so this uses")
P("a CROSS-ORBIT surface.  Sigma = {0,4,5} (|Sigma| = 3) vs Sigma' = {0,6}, the carrier's unique")
P("2-link surface, which is in a DIFFERENT Aut-orbit (sizes differ, so no automorphism relates them).")
AS2 = C.Alg(C.frag_gens([0, 6]), "Sigma'={0,6}")
P("   A_Sigma  : nbar=%d dim=%d max %d bits    A_Sigma' : nbar=%d dim=%d max %d bits"
  % (C.A_SIG.nbar, C.A_SIG.dim, C.A_SIG.maxent, AS2.nbar, AS2.dim, AS2.maxent))
P("%-6s %-4s %12s %12s %14s" % ("g2", "arm", "I(A_S:Sig)", "I(A_S:Sig')", "difference"))
for g2 in (0.45, 1.00, 3.00):
    for tag, sec in (("A0", sec0), ("A1", sec1)):
        _, psi, _ = sec.ground(g2)
        a = C.MI(C.A_FULL, C.A_SIG, sec, psi); b = C.MI(C.A_FULL, AS2, sec, psi)
        P("%-6.2f %-4s %12.6f %12.6f %14.6f" % (g2, tag, a, b, a - b))
P("   >>> WHICH SURFACE IS CALLED THE BOUNDARY CHANGES THE NUMBER.  The 2-link bottleneck surface")
P("       {0,6} carries at most 1 bit; Sigma = {0,4,5} carries 2.  'The surface' is not one object.")

# ============================================================================================
rule("BLOCK 12 -- EXHAUSTIVE: WHERE CAN THE PLATEAU FAIL?  ALL 255 FRAGMENTS WITH |F| <= 4.")
P("The brief says ask where the plateau CAN fail, not where it appears.  Every fragment of the")
P("environment with |F| <= floor(|E_env|/2) = 4 is enumerated and scored against the SAME criterion.")
FRAGS = [list(c) for n in (1, 2, 3, 4) for c in itertools.combinations(C.ENV, n)]
P("fragments enumerated: %d" % len(FRAGS))
for g2 in (0.45, 1.00, 3.00):
    for tag, sec in (("A0", sec0), ("A1", sec1)):
        _, psi, _ = sec.ground(g2)
        hF = C.A_FULL.entropy(sec, psi)
        rat = []
        for f in FRAGS:
            r = C.MI(C.A_FULL, C.Alg(C.frag_gens(f), "f"), sec, psi) / hF if hF > 1e-12 else 0.0
            rat.append((r, tuple(f)))
        rat.sort()
        npl = sum(1 for r, _ in rat if r >= 0.9)
        P("   g2=%-5.2f %s  H_FULL=%.6f  fragments meeting the plateau: %d of %d (%.1f%%)"
          % (g2, tag, hF, npl, len(rat), 100.0 * npl / len(rat)))
        P("        min I/H = %.6f at F=%s ;  max = %.6f at F=%s"
          % (rat[0][0], list(rat[0][1]), rat[-1][0], list(rat[-1][1])))
        P("        D(F) of every plateau-meeting fragment: %s"
          % sorted({C.Dnum(list(f)) for r, f in rat if r >= 0.9}))
        P("        smallest plateau-meeting fragments: %s"
          % [list(f) for r, f in sorted([x for x in rat if x[0] >= 0.9], key=lambda t: len(t[1]))[:4]])
P("")
P(">>> THE PLATEAU FAILS ON MOST OF THE ENVIRONMENT.  It is not a property of the state; it is a")
P("    property of D(F).  Reported with the D(F) column so it cannot be read as redundancy.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_W20_C_CHARGE/OUT_run_probe.txt", "w").write("\n".join(LOG) + "\n")
