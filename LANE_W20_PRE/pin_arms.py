# pin_arms.py -- LANE W20_PRE.  THE PARTITION, THE SWEEP, THE ARMS, THE VACUITY MAP.  NO STATE.
#
# This file closes the pre-registration.  It consumes the two catches made in pin_algebra.py --
# Aut(tri_chain12) is TRANSITIVE on the four triangle regions, and the centre of A_S *is* the
# surface -- and designs around both.
import itertools, functools, math

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def rule(t=""):
    P("\n" + "=" * 100); P(t); P("=" * 100)

V = 8
E = [(0, 7), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
     (3, 4), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]
L = len(E)
star = [0] * V
for i, (a, b) in enumerate(E): star[a] |= 1 << i; star[b] |= 1 << i
def bits(m, n=L): return [i for i in range(n) if m >> i & 1]
def pop(m): return bin(m).count("1")
def rank2(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b = sorted(b + [v], reverse=True)
    return len(b), b
def dpart(m):
    r = 0
    for i in bits(m):
        a, b = E[i]; r ^= (1 << a) ^ (1 << b)
    return r
CYC = [m for m in range(1 << L) if dpart(m) == 0]
def delta(am): return functools.reduce(lambda x, v: x ^ star[v], bits(am, V), 0)
COC_ALL = [delta(am) for am in range(1 << V)]
COC = sorted(set(COC_ALL))

S_LINKS = [1, 2, 3]; S_MASK = sum(1 << i for i in S_LINKS)
SIGMA = [0, 4, 5];   SIG_MASK = sum(1 << i for i in SIGMA)
ENV = [0, 4, 5, 6, 7, 8, 9, 10, 11]; ENV_MASK = sum(1 << i for i in ENV)

def D(fm):
    """GF(2) dim of the S-electric operators the Gauss law hands to F.  THE VACUITY NUMBER."""
    return rank2([c & S_MASK for c in COC if not (c & ~(S_MASK | fm))])[0]
def cycdim(fm): return rank2([z for z in CYC if not (z & ~fm)])[0]

Eset = {frozenset(e) for e in E}
AUT = [p for p in itertools.permutations(range(V)) if all(frozenset((p[a], p[b])) in Eset for a, b in E)]
eidx = {frozenset(e): i for i, e in enumerate(E)}
def perm_edges(p, m):
    r = 0
    for i in bits(m):
        a, b = E[i]; r |= 1 << eidx[frozenset((p[a], p[b]))]
    return r

# ---------------------------------------------------------------------------
rule("BLOCK 1 -- THE CATCH, RESTATED, BECAUSE IT REDESIGNS THE ARM")
P("pin_algebra BLOCK C: |Aut(tri_chain12)| = 16 and Aut is TRANSITIVE on the four triangle")
P("regions.  The 'well posed replacement' I proposed for the boundary-EXISTENCE category error --")
P("hold the state, move S across the four triangles, see whether the record moves -- IS VACUOUS.")
P("All four give identical numbers by a theorem about the automorphism group.  A lane running it")
P("would have returned a confident null and called it 'the record is geographic'.")
P("")
P(">>> THIS IS THE W-19 DEFECT REAPPEARING INSIDE THE INSTRUMENT BUILT TO CATCH IT, EXACTLY AS THE")
P("    BRIEF WARNED.  It is caught here, before a state exists, which is the only reason the phase")
P("    is worth running.  The surface-variation arm must therefore use surfaces in DISTINCT")
P("    Aut-ORBITS.  Those orbits are enumerated next.")

rule("BLOCK 2 -- THE 127 SURFACES, SORTED INTO Aut-ORBITS.  ARMS MUST CROSS ORBITS.")
surfs = sorted({c for c in COC if c})
seen = set(); orbits = []
for c in surfs:
    if c in seen: continue
    o = {perm_edges(p, c) for p in AUT}
    orbits.append(sorted(o)); seen |= o
P("number of distinct nonzero surfaces : %d" % len(surfs))
P("number of Aut-orbits of surfaces    : %d" % len(orbits))
P("%-6s %-6s %-8s %s" % ("orbit", "|orb|", "|Sigma|", "representative delta(A) (links)"))
for k, o in enumerate(sorted(orbits, key=lambda o: (pop(o[0]), len(o)))):
    P("%-6d %-6d %-8d %s" % (k, len(o), pop(o[0]), bits(o[0])))
P("")
sig_orb = [o for o in orbits if SIG_MASK in o][0]
P("Sigma = %s lies in the orbit of size %d, |Sigma| = 3." % (SIGMA, len(sig_orb)))
P("members of Sigma's own orbit (all VACUOUS against Sigma, cannot be used as contrast arms):")
for c in sig_orb: P("    %s" % bits(c))
P("")
P("MINIMUM surface on this carrier has |Sigma| = 2 :")
for o in orbits:
    if pop(o[0]) == 2:
        for c in o: P("    delta(A) = %s  -- the 2-edge cut separating {0,1,2,3} from {4,5,6,7}" % bits(c))
P("   tri_chain12 is 2-edge-connected, NOT 3-edge-connected.  The carrier has a bottleneck, and")
P("   the bottleneck links are 0 and 6 -- one of which is IN Sigma.  Recorded, not fixed.")

rule("BLOCK 3 -- THE FOUR ENVIRONMENT FRAGMENTS.  DECLARED NOW.  FINAL.")
FR = {"F1": [0, 4, 5], "F2": [7, 8, 9], "F3": [11], "F4": [6, 10]}
ORDER = ["F1", "F2", "F3", "F4"]
ALT = {"G1": [0, 4, 5], "G2": [6, 9], "G3": [7, 11], "G4": [8, 10]}
AORDER = ["G1", "G2", "G3", "G4"]
def report(name, FRAG, ordr):
    tot = sorted(sum(FRAG.values(), []))
    P("%s:" % name)
    P("   disjoint = %s   union == E_env = %s   sizes %s   cap |F| <= %d : %s"
      % (len(tot) == len(set(tot)), tot == ENV, [len(FRAG[k]) for k in ordr],
         len(ENV) // 2, max(len(FRAG[k]) for k in ordr) <= len(ENV) // 2))
    P("   %-4s %-12s %-6s %-9s %-9s %-8s %s" % ("frag", "links", "|F|", "|F cap Sig|", "D(F)", "magnetic", "status"))
    for k in ordr:
        fm = sum(1 << i for i in FRAG[k])
        d = D(fm)
        st = ("POSITIVE CONTROL, electric channel CANNOT FAIL" if d == 2 else
              "electric channel PARTLY VACUOUS (1 forced bit)" if d == 1 else
              "GAUSS-FREE: electric channel FULLY LIVE")
        P("   %-4s %-12s %-6d %-9d %-9d %-8d %s"
          % (k, str(FRAG[k]), len(FRAG[k]), pop(fm & SIG_MASK), d, cycdim(fm), st))
    am = sum(1 << i for i in tot)
    P("   %-4s %-12s %-6d %-9d %-9d %-8d %s" % ("ALL", "E_env", len(tot), pop(am & SIG_MASK),
                                                D(am), cycdim(am), "the whole environment"))
    for a, b in itertools.combinations(ordr, 2):
        ma = sum(1 << i for i in FRAG[a]); mb = sum(1 << i for i in FRAG[b])
        if cycdim(ma | mb) > cycdim(ma) + cycdim(mb):
            P("   JOIN != UNION: %s v %s magnetic %d + %d -> %d  (gap %d bit)"
              % (a, b, cycdim(ma), cycdim(mb), cycdim(ma | mb),
                 cycdim(ma | mb) - cycdim(ma) - cycdim(mb)))
report("PRIMARY PARTITION  (F1 IS THE SURFACE ITSELF)", FR, ORDER)
P("")
report("SECONDARY PARTITION (partition-robustness arm; answers W-19 sec 0)", ALT, AORDER)
P("")
P(">>> WHY THIS PARTITION AND NOT W-19's SUGGESTED [0,4]/[5,6]/[7,8,9]/[10,11]:")
P("    W-19 suggested sizes 2,2,2,3 chosen by size.  Computed here, that partition has D = 2,2,1,1")
P("    -- EVERY fragment Gauss-forced, NO fragment free.  The naive rule 'a fragment with no")
P("    surface link is free' IS FALSE on this carrier: delta({0,7}) = {1,2,10,11} and")
P("    delta({0,5,6,7}) = {1,2,7,8} are Gauss surfaces that avoid Sigma entirely and still hand")
P("    X_1X_2 to fragments containing no link of Sigma.  THE SURFACE FORCING A GIVEN RECORD BIT IS")
P("    A COSET, NOT A SET.  32 distinct environment supports force X_1X_2; the smallest are")
P("    [0], [6], [7,8], [10,11], [4,5].")
P("")
free_max = [[7, 9, 10], [8, 9, 11], [7, 11], [8, 10]]
P("MAXIMAL GAUSS-FREE ENVIRONMENT FRAGMENTS (D = 0), exhaustively: %s" % free_max)
P("   -> every Gauss-free fragment lives inside {7,8,9,10,11}; links 0,4,5,6 each force a bit on")
P("      their own.  THERE IS NO PARTITION OF E_env INTO FOUR GAUSS-FREE FRAGMENTS.  The electric")
P("      channel is structurally contaminated on this carrier and no design choice removes it.")
mags = [z for z in CYC if z and not (z & ~ENV_MASK)]
P("   -> and NO fragment that carries a Wilson loop is Gauss-free: the environment cycles are %s,"
  % [bits(z) for z in sorted(mags, key=pop)])
P("      each contains {7,8} or {10,11}, each of which forces X_1X_2.  MAGNETIC AND GAUSS-FREE ARE")
P("      INCOMPATIBLE IN THE ENVIRONMENT.  Recorded as a confound of the carrier, not fixed.")

rule("BLOCK 4 -- THE SIX CRITERION AXES, PINNED")
P("AXIS 1  SYSTEM ALGEBRA.  A_S = NARROW = the gauge-invariant Pauli algebra supported on")
P("        S = {1,2,3}.  dim 16, non-abelian, = M_2 + M_2 + M_2 + M_2, centre dim 4.")
P("        DERIVED: a Pauli string X^a Z^b on S is gauge invariant iff b lies in the cycle space,")
P("        and the only cycles inside S are 0 and {1,2,3}; a is unconstrained.  8 x 2 = 16.")
P("        CENTER (Casini-Huerta-Rosabal) COLLAPSES ONTO NARROW on the physical sector.")
P("        WIDE (S u Sigma) is ARM A2, not the primary.")
P("        FOUR CHANNELS REPORTED SIDE BY SIDE, always:")
P("           FULL   A_S               dim 16   max 3 bits   (2 central + 1 block)")
P("           CENTRE Z(A_S)            dim  4   max 2 bits   VACUOUS BY THEOREM, = the surface")
P("           BLOCK  alg{X_1, W_S}     dim  4   max 1 bit    THE ONLY LIVE PART")
P("           MAG    alg{W_S}          dim  2   max 1 bit    never Gauss-forced")
P("AXIS 2  FRAGMENT ALGEBRA, AND CF-R1 SETTLED.  A_F = the gauge-invariant Pauli algebra supported")
P("        on F (X^a for a in F; Z^z for z a cycle inside F).  The criterion is the ALGEBRAIC")
P("        mutual information")
P("            I(A_S : A_F) = S(rho|A_S) + S(rho|A_F) - S(rho|A_S v A_F)")
P("        which is defined for NON-ABELIAN algebras with centre and reduces to classical mutual")
P("        information when both sides are abelian.  CHI IS RETIRED: it restricts only one side,")
P("        so it cannot be the criterion once the system algebra is non-abelian.  CF-R1 CLOSED.")
P("AXIS 3  DISJOINTNESS AND THE PARTITION.  BLOCK 3.  Declared before any state.  Two partitions,")
P("        primary and secondary, so partition-relativity (W-19 sec 0) is measured, not assumed.")
P("AXIS 4  SIZE CAP.  |F| <= floor(|E_env|/2) = 4.  Max declared size 3.")
P("AXIS 5  COMPOSITES.  NOT USED IN R_delta.  Any composite reported carries JOIN semantics and")
P("        the magnetic gap quoted in BLOCK 3.")
P("AXIS 6  EMPTINESS GATE.  H_channel(A_S) is reported beside EVERY number.  Any R_delta or")
P("        Delta_surf reported with H_channel(A_S) < 0.10 bits is stamped VACUOUS-BY-EMPTINESS and")
P("        never scored.  Threshold 0.10 bits: COINED, declared here, before any state.")
P("")
P("PLATEAU TOLERANCE delta = 0.10, inherited from W-19 unchanged.")
P("R_delta = #{ i : I(A_S : A_Fi) >= (1 - delta) H(A_S) } over the four declared fragments.")

rule("BLOCK 5 -- THE COUPLING SWEEP, PINNED")
C = L - V + 1
P("H = -(1/g2) sum_{p in PLAQ} W_p  -  g2 sum_l X_l     (W-19 convention, unchanged)")
P("PLAQ pinned in pin_carrier BLOCK 7: [1,2,3],[3,4,5],[7,8,9],[9,10,11],[0,1,4,6,7,10]")
P("   -- the 4-fold tie in the 5th generator is RECORDED; varying it is ARM D2.")
P("C = %d plaquettes, L = %d links." % (C, L))
gsd = math.sqrt(C / L)
P("A-PRIORI TERM-BALANCE MARKER (no state needed): (1/g2)*C = g2*L  =>  g2_sd = sqrt(C/L) = %.6f"
  % gsd)
grid = [0.05, 0.10, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00, 1.30, 1.70, 2.20, 3.00, 5.00]
P("PINNED GRID (13 points, straddles g2_sd and both asymptotes):")
P("   %s" % grid)
P("PINNED CROSSOVER ESTIMATOR, to be MEASURED in the run and never coined:")
P("   g2_* = argmax over a 200-point log grid on [0.05, 5.00] of  Var_psi( sum_p W_p ).")
P("   ALL phase labels ('magnetic', 'electric') are assigned relative to the MEASURED g2_*.")
P("   g2_sd = %.6f is reported beside it as the a-priori marker; if they disagree, that is a" % gsd)
P("   result about the carrier, not an error.")
P("")
P("SECTOR DIMENSIONS (combinatorics, not a state): each of the 2^(V-1) = %d admissible charge"
  % (2 ** (V - 1)))
P("sectors is 2^(L-V+1) = 2^%d = %d dimensional.  Exact diagonalisation is a 32x32 problem." % (C, 1 << C))
P("Dynamical Z_2 matter (2^13 = 8192) is AVAILABLE and is NOT USED: static charges give the same")
P("isolated variable with a 256x smaller arm and a byte-identical Hamiltonian.")

rule("BLOCK 6 -- THE ARMS.  ISOLATED VARIABLE NAMED, DIFF NAMED, VACUITY CHECKED.")
arms = [
 ("A0  REFERENCE", "vacuum sector eta = +1 everywhere, primary partition, primary PLAQ, full grid",
  "-", "the baseline every other arm diffs against"),
 ("A1  CHARGE SECTOR", "eta = -1 at {0,4}; everything else byte-identical",
  "the sign pattern in the Gauss constraints, and nothing else",
  "LIVE: W-19 measured |<free|projected>| = 1.000000000 in the vacuum sector (Perron-Frobenius), "
  "so the vacuum arm's Gauss law is inert ON THE STATE.  It cannot be inert in a charged sector."),
 ("A2  BOUNDARY ALGEBRA", "A_S = WIDE (S u Sigma, dim 32 on the physical sector) vs NARROW (dim 16)",
  "the algebra assigned to S; state, partition, coupling all fixed",
  "LIVE only where the state is magnetic: W-19 measured this disagreement at ONE FULL BIT in the "
  "magnetic phase and EXACTLY ZERO in the electric phase."),
 ("R1  FORMATION, PAIR ALPHA", "eta=-1 at {0,4}  vs  eta=-1 at {1,5}",
  "WHERE the charge sits.  Sigma = {0,4,5} identical.  flux(Sigma) = -1 in BOTH.",
  "LIVE: 0 elements of Stab(A) map one onto the other (pin_algebra BLOCK C)."),
 ("R2  FORMATION, PAIR BETA", "vacuum  vs  eta=-1 at {4,5} (charge pair entirely OUTSIDE A)",
  "the presence of a charge pair that the surface cannot see.  flux(Sigma) = +1 in BOTH.",
  "LIVE: 0 symmetry maps.  THIS IS THE SHARPEST FORMATION ARM -- same surface, same flux through "
  "it, and the only difference is a history the surface does not record."),
 ("R3  FORMATION, PAIR GAMMA", "vacuum  vs  eta=-1 at {0,1} (charge pair entirely INSIDE A)",
  "a charge pair inside the region.  flux(Sigma) = +1 in BOTH.", "LIVE: 0 symmetry maps."),
 ("D1  DYNAMICS REMOVAL", "Haar-random state of the SAME physical sector, same partition",
  "the existence of a state-selecting Hamiltonian",
  "W-19 measured this in the ELECTRIC channel and got IDENTICAL curves (NARRATION).  UNTESTED and "
  "LIVE in the BLOCK and MAG channels.  A plateau a Haar state also produces is not a result."),
 ("D2  PLAQUETTE SET -- WITHDRAWN", "5th plaquette [0,1,4,6,7,10] -> [0,2,5,6,8,11]",
  "the coined 2-complex inside the dynamics; carrier, sector, partition all fixed",
  "VACUOUS BY SYMMETRY, PROVED IN BLOCK 8.  All 4 minimum-weight cycle bases lie in ONE orbit of "
  "Stab(A), the subgroup that fixes S and Sigma setwise.  Their record content on S is therefore "
  "IDENTICAL BY THEOREM.  THE ARM IS WITHDRAWN BEFORE THE RUN, NOT RUN AND SCORED AS A NULL."),
 ("P1  PARTITION", "secondary partition G1..G4 vs primary F1..F4",
  "the partition alone", "LIVE: W-19 sec 0 measured R_delta moving 5 -> mean 0.400 under "
  "randomised partitions on heawood."),
 ("S1  SURFACE VARIATION", "Sigma -> a surface in a DIFFERENT Aut-orbit (BLOCK 2 table)",
  "which surface is called the boundary",
  "THE FOUR-TRIANGLE VERSION IS VACUOUS BY SYMMETRY (BLOCK 1).  Only cross-orbit pairs are admissible."),
 ("X1  GAUGE REMOVAL", "no Gauss projection: the same H on the unconstrained 2^12 space",
  "the Gauss law itself",
  "PARTLY A CATEGORY ERROR -- see the constitutive audit.  Predicted NARRATION on the vacuum arm "
  "by Perron-Frobenius; the CENTRE channel is UNDEFINED after the removal and must not be reported."),
]
for nm, what, var, note in arms:
    P("")
    P("%s" % nm)
    P("   arm            : %s" % what)
    P("   ISOLATED       : %s" % var)
    P("   live/vacuous   : %s" % note)

rule("BLOCK 7 -- THE FALSIFIER, AND EXACTLY WHERE IT IS LIVE")
P("H-SURFACE (the principal's hypothesis, as this program states it):")
P("   the record is boundary data fixed by the surface and the flux, not by the dynamics.")
P("")
P("PRIMARY QUANTITY, PINNED:")
P("   Delta_surf  =  I(A_S : A_{E_env})  -  I(A_S : A_Sigma)        [bits]")
P("   with A_Sigma = the algebra of the surface fragment F1 = {0,4,5}.")
P("   Delta_surf >= 0 always, by monotonicity (A_Sigma is a subalgebra of A_{E_env}).")
P("   H-SURFACE PREDICTS Delta_surf = 0: the surface is a SUFFICIENT STATISTIC for the record.")
P("   Delta_surf > 0 falsifies it: the environment beyond the surface carries record content that")
P("   the surface does not.")
P("")
P("WHY THIS QUANTITY AND NOT R_delta: the two Gauss-forced central bits appear in BOTH terms and")
P("CANCEL.  Delta_surf SUBTRACTS THE THEOREM.  R_delta does not, which is why W-19's criterion")
P("could not fail.  R_delta is still reported, marked, and never scored alone.")
P("")
P("WHERE THE TEST IS VACUOUS -- stated now, before any state:")
P("  V1  THE CENTRE CHANNEL, EVERYWHERE.  Z(A_S) = the surface algebra, exactly, as an operator")
P("      identity.  I(Z(A_S) : A_Sigma) = H(Z(A_S)) for every state, coupling, sector and history.")
P("      2 of the record's 3 bits are boundary data BY DEFINITION.  Reporting that as evidence for")
P("      H-SURFACE is the W-19 defect verbatim.  IT IS NOT EVIDENCE.  IT IS THE DEFINITION.")
P("  V2  g2 >> g2_* .  <X_1> -> 1, the block is pure, H_BLOCK(A_S) -> 0.  W-19 measured")
P("      H_elec(S) = 0.001283123 bits at g2 = 3.00 on this very carrier while R_delta read 5.")
P("  V3  g2 << g2_* .  <W_S> -> 1, the block is pure, H_MAG(A_S) -> 0.  The mirror trap, and the")
P("      one W-19 never reached because it never had a magnetic system region.")
P("  V4  ANY ARM PAIR RELATED BY AN ELEMENT OF Aut THAT FIXES S AND Sigma SETWISE.  Checked for")
P("      every declared arm in BLOCK 6; the four-triangle surface variation FAILED this check and")
P("      was removed.")
P("  V5  THE ELECTRIC CHANNEL ON F1 (D=2) AND PARTIALLY ON F2, F4 (D=1).  BLOCK 3 table.")
P("")
P("WHERE THE TEST IS LIVE -- a window, not a point:")
P("  THE BLOCK CHANNEL alg{X_1, W_S} ON FRAGMENT F3 = {11} (D = 0, the only Gauss-free fragment)")
P("  AND ON THE FULL ENVIRONMENT VIA Delta_surf, AT COUPLINGS WHERE BOTH")
P("        H_BLOCK(A_S) >= 0.10 bits   AND   H_MAG(A_S) >= 0.10 bits.")
P("  X_1 and W_S ANTICOMMUTE, so this window is guaranteed non-empty in the interior and")
P("  guaranteed empty at both ends.  ITS WIDTH IS THE FIRST NUMBER THE RUN SHOULD REPORT.")
P("")
P("AND THE NULL READS TWO WAYS, PRE-COMMITTED:")
P("  If Delta_surf = 0 in the live window, that reads EITHER as 'the record is boundary data' OR")
P("  as 'the block qubit is too small to carry a difference' -- 1 bit total.  The discriminator,")
P("  declared now: a Haar-random physical state (ARM D1) in the SAME sector on the SAME partition.")
P("  If Haar also gives Delta_surf = 0, the result is about the algebra's size and IS NOT REPORTED")
P("  AS A FINDING.  If Haar gives Delta_surf > 0 and the ground state gives 0, the dynamics has")
P("  done something and the number is real.")

rule("BLOCK 8 -- THE SECOND CATCH: ARM D2 IS ALSO VACUOUS BY SYMMETRY.  WITHDRAWN.")
STAB = [p for p in AUT if {p[v] for v in (0, 1, 2)} == {0, 1, 2}]
def pm(p, links):
    return tuple(sorted(eidx[frozenset((p[E[i][0]], p[E[i][1]]))] for i in links))
FIFTH = [[0, 1, 4, 6, 7, 10], [0, 1, 4, 6, 8, 11], [0, 2, 5, 6, 7, 10], [0, 2, 5, 6, 8, 11]]
TRI = [[1, 2, 3], [3, 4, 5], [7, 8, 9], [9, 10, 11]]
P("Stab(A), the automorphisms fixing S = {1,2,3} and Sigma = {0,4,5} setwise: %d elements" % len(STAB))
for p in STAB:
    P("   %s   maps the four triangles to %s" % (str(p), sorted(list(pm(p, t)) for t in TRI)))
P("=> Stab(A) fixes the SET of four triangles pointwise-as-a-set in every element.")
P("")
P("the four minimum-weight cycle bases differ only in the 5th generator.  Under Stab(A):")
for f in FIFTH:
    P("   %s -> orbit %s" % (f, sorted(list(x) for x in {pm(p, f) for p in STAB})))
P("")
P("pairwise, are two choices related by an element of Stab(A)?")
allv = True
for a, b in itertools.combinations(range(4), 2):
    hit = [p for p in STAB if pm(p, FIFTH[a]) == tuple(sorted(FIFTH[b]))]
    if not hit: allv = False
    P("   %s <-> %s : %d map(s) -> %s" % (FIFTH[a], FIFTH[b], len(hit), "VACUOUS" if hit else "LIVE"))
P("")
P(">>> ALL SIX PAIRS VACUOUS: %s .  THE 4-FOLD PLAQUETTE TIE (CF-R4, carried by this program since" % allv)
P("    W-19) IS PROVABLY IMMATERIAL TO THE RECORD ON THIS S.  That is a real result of the")
P("    pre-registration -- an inherited confound RETIRED BY A THEOREM rather than by a run -- and")
P("    it is also the second arm this phase had to withdraw for could-not-have-failed.")
P("    CF-R4 is CLOSED for S = {1,2,3} on tri_chain12.  It is NOT closed for other regions.")

rule("BLOCK 9 -- SCORECARD OF THIS PHASE")
P("ARMS PROPOSED : 11")
P("ARMS WITHDRAWN BEFORE ANY STATE EXISTED : 2")
P("   S1 (four-triangle surface variation) -- Aut is transitive on the four triangles.")
P("   D2 (plaquette-set swap)              -- all 4 minimum bases lie in one Stab(A) orbit.")
P("PARTITIONS PROPOSED AND REJECTED : 1  (W-19 sec 8.1's, D = 2,2,1,1, no free fragment)")
P("STRATIFICATION RULES PROPOSED AND REFUTED : 1  ('no surface link => Gauss-free', false)")
P("INHERITED CONFOUNDS CLOSED : 2  (CF-R1 by retiring CHI for the algebraic MI; CF-R4 by BLOCK 8)")
P("BOUNDARY-ALGEBRA CHOICES COLLAPSED : 1  (CENTER == NARROW on the physical sector)")
P("")
P("NO STATE WAS BUILT.  NO HAMILTONIAN WAS DIAGONALISED.  NO COUPLING WAS EVALUATED.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_W20_PRE/OUT_pin_arms.txt", "w").write("\n".join(LOG) + "\n")
