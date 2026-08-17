# W10-D REFUTE-2  LENS 2 = COMPLETENESS.  LEG 2.
#
# FOUR REGISTERED LOAD-BEARING CLAIMS THAT APPEAR IN NO ROW OF LANE D's SCOPE TABLE, EACH
# CHECKED HERE ON THE SAME FOUR-CLASS CARRIERS LANE D USED.
#
#  2A  W-02 (REGISTER:118): "THE TRAP WAS REAL AND IS DISARMED BY COMPUTATION.  Repeated
#      circuits of one loop span 3 DIMENSIONS at N=1 and 3 at N=100."  S3:264-271 states it as
#      "circuits generate EXACTLY the three-way split algebra C^3 of sec2.3 and nothing more."
#      Lane D's table has no row for it.  It is item 2 of S3's own headline and it is what
#      licenses "the directed system is forced to be multiplicative, not cellular".
#
#  2B  W-02 (REGISTER:139): "P-9 is resolved: ||omega_F^N - omega_C^N|| -> 2.000000000000."
#      Lane D row 2.9 marks ALL NINE "UNDETERMINED ... to determine their scope one must re-run
#      S3 sec5 with V = 9 and V = 6, which requires S3's code -- S3 published no lane directory."
#      S3 sec5.9 gives P-9 IN CLOSED FORM: ||omega_F^N - omega_C^N|| = 2 sqrt(1 - |Omega_N|^2).
#      |Omega_N| is the object lane D itself computed on BOTH four-class carriers (its legs 2B,
#      3B, 3E).  No S3 code is required for P-9, and none was required for P-1 or P-3 either.
#
#  2C  S3 sec5.9's exception clause: "at finite N on W-01's FIRING LOCUS, P-9 HOLDS exactly."
#      Lane D's own leg 1D proves neither four-class carrier has a torus zero at its own
#      published weights.  The clause is therefore EMPTY on both, and lane D never says so.
#
#  2E  S4 sec4.1: "B5, the degenerate end, has b1 = 0: no free cycle exists, so gamma_C cannot
#      be designated and THE FORMATION DATUM DOES NOT EXIST ON IT.  This is the one place in S4
#      where topology, and nothing else, decides an outcome."  Lane D row 3.4 rules that the
#      topology-inertness row is "UNTESTABLE BY CONSTRUCTION ... no experiment on any carrier
#      can bear on this row."  The corpus owns the carrier that bears on it.
#
# THE ONE VARIABLE PER LEG IS NAMED AT THE HEAD OF THAT LEG, AND EVERY LEG PRINTS ITS ARMS.
# PRECISION: numpy float64; the rank statistics are checked at two tolerances and the
# |Omega_N| products are additionally accumulated in logs to avoid underflow.

import numpy as np
from collections import Counter

np.set_printoptions(precision=12)

# ---------------------------------------------------------------- carriers, from S4:575-590
# class label = (is v in gamma_F ?)(is v in gamma_C ?);  chi(00)=1, chi(10)=u, chi(01)=v, chi(11)=uv
CARRIERS = {
    'B1  K1 as handed  ': {'01': 2, '10': 2, '11': 1},
    'B1q + spectator   ': {'00': 1, '01': 3, '10': 3},
    'B1p bridged       ': {'01': 3, '10': 3},
    'B1s subdivided    ': {'01': 5, '10': 5, '11': 1},
    'B0a torus disjoint': {'00': 2, '01': 3, '10': 4},
    'B0b torus meeting ': {'00': 4, '01': 1, '10': 2, '11': 2},
    'B4  spindle       ': {'00': 1, '01': 1, '10': 1, '11': 3},
}
ORDER = ['00', '10', '01', '11']


def vertex_classes(mult):
    out = []
    for c in ORDER:
        out += [c] * mult.get(c, 0)
    return out


def weights(mult):
    V = sum(mult.values())
    return np.array([mult.get(c, 0) / V for c in ORDER])


print("=" * 104)
print("ARM DIFF FIRST.  The carriers, their vertex-class strings, and their SENSE-U weight")
print("vectors.  Distinctness is checked on the class STRING, which is what every leg below")
print("actually consumes -- not on the label.")
print("=" * 104)
strs = {}
for nm, mult in CARRIERS.items():
    vc = vertex_classes(mult)
    strs[nm] = "".join(vc)
    print(f"  {nm}  V={sum(mult.values()):2d}  classes = {''.join(vc)}   "
          f"p = {np.array2string(weights(mult), precision=6)}")
dups = [(a, b) for i, a in enumerate(strs) for b in list(strs)[i + 1:] if strs[a] == strs[b]]
print(f"  DISTINCT ARMS: {len(set(strs.values()))} of {len(strs)}" +
      ("" if not dups else f"   !! IDENTICAL: {dups}"))
print("  NOTE, DECLARED: B0b and B4 have DIFFERENT class strings, so for the legs below they")
print("  are two arms; for W-01's hull criterion they are one arm (lane D's own finding, kept).")

# =============================================================== 2A  the circuit-algebra trap
print("\n" + "=" * 104)
print("== 2A  W-02's TRAP FIGURE: 'repeated circuits span 3 DIMENSIONS at N=1 and 3 at N=100' ==")
print("=" * 104)
print("  ONE VARIABLE: the carrier's vertex-class string.  Connection, rank routine, tolerance")
print("  and N-list identical in every row.  W_F = e^{i}, W_C = e^{i sqrt2} (generic, fixed).")
WF, WC = np.exp(1j * 1.0), np.exp(1j * np.sqrt(2))


def circuit_span(vc, N):
    """rank of span{ M_dF^n, M_c^n : n = 0..N } inside M_V(C), as S3:264 measures it."""
    a = np.array([1.0 if c[0] == '1' else 0.0 for c in vc])
    b = np.array([1.0 if c[1] == '1' else 0.0 for c in vc])
    rows = []
    for n in range(0, N + 1):
        rows.append(WF ** (n * a))
        rows.append(WC ** (n * b))
    M = np.array(rows)
    return np.linalg.matrix_rank(M, tol=1e-9), np.linalg.matrix_rank(M, tol=1e-6)


def generated_algebra(vc, N):
    """rank of span{ M_dF^i M_c^j } -- the ALGEBRA the two transports generate."""
    a = np.array([1.0 if c[0] == '1' else 0.0 for c in vc])
    b = np.array([1.0 if c[1] == '1' else 0.0 for c in vc])
    rows = [WF ** (i * a) * WC ** (j * b) for i in range(N + 1) for j in range(N + 1)]
    M = np.array(rows)
    return np.linalg.matrix_rank(M, tol=1e-9)


print(f"\n  {'carrier':20s} {'V':>3s} {'#classes':>9s}" +
      "".join(f"{'span N='+str(N):>12s}" for N in (1, 2, 5, 10, 25, 100)) +
      f"{'ALGEBRA':>9s}")
for nm, mult in CARRIERS.items():
    vc = vertex_classes(mult)
    row = f"  {nm:20s} {len(vc):3d} {len(mult):9d}"
    for N in (1, 2, 5, 10, 25, 100):
        r9, r6 = circuit_span(vc, N)
        row += f"{r9:12d}" if r9 == r6 else f"{str(r9)+'/'+str(r6):>12s}"
    row += f"{generated_algebra(vc, 6):9d}"
    print(row)
print("\n  S3:264-271 REPRODUCED ON K1: dim = 3, constant in N, and equal to the class algebra")
print("  C^3 -- 'circuits generate EXACTLY the three-way split algebra of sec2.3'.")
print("  OFF K1 THE SENTENCE SPLITS IN TWO AND ONLY ONE HALF SURVIVES:")
print("    - 'constant in N / circuits grow no algebra'  -- CARRIER-INDEPENDENT (every row flat).")
print("    - 'the dimension is 3'                        -- THREE-CLASS-SCOPED (it is the number")
print("      of occupied classes on K1 and B1q; it is 2 on B1p).")
print("    - 'EXACTLY the class algebra, and nothing more' -- FALSE AT FOUR CLASSES: the span of")
print("      PURE POWERS is 3 on B0b and B4 while the algebra the same two transports GENERATE")
print("      is 4.  On a four-class carrier the circuits span a PROPER subspace of the algebra")
print("      they generate, and the coincidence S3 reports is an artefact of p00 = 0.")
print("  (Reason, one line: the four class-indicator combinations reachable by pure powers are")
print("   e10+e11, e00+e01, e01+e11, e00+e10, whose sum-relation (e10+e11)+(e00+e01) =")
print("   (e01+e11)+(e00+e10) drops the rank to 3; with class 00 empty the relation is vacuous.)")

# =============================================================== 2B  P-9 in closed form
print("\n" + "=" * 104)
print("== 2B  P-9 IS A CLOSED FORM IN |Omega_N| -- NO S3 CODE IS REQUIRED TO SCOPE IT ==")
print("=" * 104)
print("  S3 sec5.9 verbatim: 'The state-norm distance is ||omega_F^N - omega_C^N|| =")
print("  2 sqrt(1 - |Omega_N|^2), maximal at 2.'  Step one: reproduce S3's OWN published table")
print("  on K1 at ITS OWN test connection f=2.0, c=1.1, p=(0.4,.15,.15,.15,.15).")


def omega_log(p, f, c, N):
    k = np.arange(1, N + 1)
    Z = p[0] + p[1] * np.exp(-1j * f * k) + p[2] * np.exp(1j * c * k) + p[3] * np.exp(1j * (c - f) * k)
    return np.cumsum(np.log(np.abs(Z))), np.abs(Z)


S3_TABLE = {1: 4.112706e-01, 2: 1.934559e-01, 5: 1.772339e-02, 10: 6.540411e-04,
            20: 2.514545e-07, 42: 2.505486e-15, 100: 4.036647e-34}
pS3 = np.array([0.0, 0.3, 0.3, 0.4])
L, _ = omega_log(pS3, 2.0, 1.1, 100)
print(f"    {'N':>4s} {'|Omega_N| mine':>16s} {'S3 published':>16s} {'rel dev':>10s} "
      f"{'||om_F-om_C|| mine':>20s}")
for N, v in S3_TABLE.items():
    mine = np.exp(L[N - 1])
    nrm = 2 * np.sqrt(max(0.0, 1 - mine ** 2))
    print(f"    {N:4d} {mine:16.6e} {v:16.6e} {abs(mine-v)/v:10.2e} {nrm:20.12f}")
print("  S3's table REPRODUCED from the class weight vector alone.  P-9's whole content at every")
print("  finite N is |Omega_N|, and its limit statement is |Omega_N| -> 0, i.e. G != {1}.")
print("\n  STEP TWO -- THE SCOPE QUESTION LANE D DECLINED, ANSWERED ON BOTH FOUR-CLASS CARRIERS.")
print("  ONE VARIABLE: the class weight vector.  Same (f,c) = (2.0, 1.1), same code path.")
print(f"    {'carrier':20s}" + "".join(f"{'N='+str(N):>16s}" for N in (1, 5, 20, 100)) +
      f"{'||.|| at N=100':>16s}")
for nm, mult in CARRIERS.items():
    p = weights(mult)
    L, _ = omega_log(p, 2.0, 1.1, 100)
    row = f"    {nm:20s}"
    for N in (1, 5, 20, 100):
        row += f"{np.exp(L[N-1]):16.6e}"
    om = np.exp(L[99])
    row += f"{2*np.sqrt(max(0.0,1-om**2)):16.12f}"
    print(row)
print("  P-9 HOLDS IN THE LIMIT ON EVERY CARRIER WITH G != {1}, FOUR CLASSES INCLUDED, AND ITS")
print("  SCOPE IS EXACTLY W-08's G != {1} CRITERION -- which lane D itself scoped")
print("  CARRIER_INDEPENDENT at its rows 2.1 and 8.4.  The verdict was already in the lane's")
print("  own hands; the row says UNDETERMINED.")

# =============================================================== 2C  the firing-locus clause
print("\n" + "=" * 104)
print("== 2C  P-9's FINITE-STAGE EXCEPTION CLAUSE IS EMPTY ON BOTH FOUR-CLASS CARRIERS ==")
print("=" * 104)
print("  S3 sec5.9: 'At finite N on W-01's firing locus: if Z_{k_n} = 0 for some n -- exactly")
print("  W-01's convex-hull criterion -- then Omega_N = 0 EXACTLY for all N >= n.  On K1's own")
print("  published connection this occurs at the FIRST CELL.  HOLDS exactly, at finite stage.'")
print("  ONE VARIABLE: the class weight vector.  Torus-zero test is exact (a linear equation in")
print("  cos t); the grid minimum is a 2000x2000 cross-check on the same weights.")
print(f"\n  {'carrier':20s} {'A':>11s} {'B':>11s} {'torus zero?':>12s} {'min|Z_1| on 2000^2 grid':>25s}")
g = np.linspace(0, 2 * np.pi, 2000, endpoint=False)
FF, CC = np.meshgrid(g, g, indexing='ij')
for nm, mult in CARRIERS.items():
    p = weights(mult)
    p00, p10, p01, p11 = p
    A = p00 ** 2 + p01 ** 2 - p10 ** 2 - p11 ** 2
    B = 2 * (p00 * p01 - p10 * p11)
    if abs(B) > 1e-300:
        tz = -1.0 <= -A / B <= 1.0
    else:
        tz = abs(A) < 1e-15
    Z1 = p00 + p10 * np.exp(-1j * FF) + p01 * np.exp(1j * CC) + p11 * np.exp(1j * (CC - FF))
    print(f"  {nm:20s} {A:11.6f} {B:11.6f} {str(tz):>12s} {np.abs(Z1).min():25.9f}")
print("\n  BOTH FOUR-CLASS CARRIERS: NO TORUS ZERO AT THEIR OWN PUBLISHED WEIGHTS, so there is")
print("  NO connection at which W-01's criterion fires, so P-9's finite-stage exception clause")
print("  is VACUOUS on them.  On K1 the clause is not vacuous and S1's published connection is")
print("  in it.  Lane D's leg 1D computes the same torus-zero fact and never connects it to")
print("  P-9, because P-9 sits inside a block it marked UNDETERMINED.")

# =============================================================== 2D  P-1 and P-3
print("\n" + "=" * 104)
print("== 2D  P-1 AND P-3 ARE ALSO ALREADY DECIDED BY ROWS LANE D SCOPED ==")
print("=" * 104)
print("  P-3 'thresholded non-return (EXACT, by monotonicity, not asymptotic)' is W-08's")
print("  |Z_k| <= 1, which lane D scoped CARRIER_INDEPENDENT at its row 8.1.  Re-exhibited:")
for nm, mult in CARRIERS.items():
    p = weights(mult)
    _, Zabs = omega_log(p, 2.0, 1.1, 200000)
    print(f"    {nm:20s} max_k(|Z_k| - 1) over k<=2e5 = {(Zabs-1).max():+.3e}"
          f"   monotone non-increasing: {bool((Zabs<=1+1e-12).all())}")
print("  P-1 'holds in the limit, FAILS at every finite N NECESSARILY' -- S3 sec5.3's own")
print("  statement is that this 'is a theorem about ALL finite constructions', i.e. it is")
print("  carrier-independent by the wording of the claim itself, with no computation needed.")
print("  THREE OF THE NINE ARE THEREFORE DETERMINABLE WITHOUT S3's CODE.  A fourth, P-7")
print("  sector-hood, is Hepp disjointness of the two limit states, which is again |Omega_N|")
print("  -> 0.  Lane D's blanket UNDETERMINED over all nine is over-broad by at least three.")

# =============================================================== 2E  the topological premise
print("\n" + "=" * 104)
print("== 2E  THE TOPOLOGICAL PRECONDITION NO ROW OF THE SCOPE TABLE STATES ==")
print("=" * 104)
print("  S4:519's own table, quoted.  b1 is the number of independent flat cycles; rank(d2) is")
print("  the number of independent curvatures.  gamma_F needs rank(d2) >= 1 (it must BOUND) and")
print("  gamma_C needs b1 >= 1 (it must NOT).")
S4ROWS = [  # name, V,E,F, chi, b0,b1,b2, has gamma_C
    ('B0a ring torus disjoint', 9, 18, 9, 0, 1, 2, 1, True),
    ('B0b ring torus meeting ', 9, 18, 9, 0, 1, 2, 1, True),
    ('B3  horn torus         ', 5, 12, 8, 1, 1, 1, 1, True),
    ('B1  K1 as handed       ', 5, 6, 1, 0, 1, 1, 0, True),
    ('B4  spindle            ', 6, 8, 4, 2, 1, 1, 2, True),
    ('B5  double sphere      ', 4, 4, 2, 2, 1, 0, 1, False),
    ('B2  K1 both filled     ', 5, 6, 2, 1, 1, 0, 0, True),
    ('B1p K1-bridged         ', 6, 7, 1, 0, 1, 1, 0, True),
    ('B1q + spectator        ', 7, 8, 1, 0, 1, 1, 0, True),
    ('B1s K1 subdivided      ', 11, 12, 1, 0, 1, 1, 0, True),
]
print("  MY OWN CONFOUND, RECORDED NOT PATCHED.  The first version of this leg tested the single")
print("  condition b1 >= 1 and flagged BOTH B5 and B2 as 'formation datum does not exist'.  THAT")
print("  IS WRONG ABOUT B2: S4's designated-loop table gives B2 'gC bounds = TRUE, independent =")
print("  True', i.e. B2 HAS two designable loops and S4 calls it 'the deliberate exception' --")
print("  both bound, so its second holonomy is a second CURVATURE, not a flat one.  There are")
print("  TWO preconditions, not one, and the corpus's own ten-carrier table breaks each of them")
print("  on a different carrier.  Both statistics are printed.")
print(f"  {'carrier':25s} {'V':>3s} {'cyc rk':>7s} {'b1':>3s} {'rk d2':>6s} "
      f"{'2 loops?':>9s} {'gC flat?':>9s} {'(1st run: b1>=1)':>17s}")
for nm, V, E, F, chi, b0, b1, b2, hasC in S4ROWS:
    cyc = E - V + b0
    rkd2 = F - b2
    note = ""
    if cyc < 2:
        note = "  <-- S4: 'the formation datum DOES NOT EXIST on it'"
    elif b1 == 0:
        note = "  <-- S4: 'the deliberate exception -- no flat holonomy at all'"
    print(f"  {nm:25s} {V:3d} {cyc:7d} {b1:3d} {rkd2:6d} {str(cyc >= 2):>9s} "
          f"{str(b1 >= 1):>9s} {str(b1 >= 1):>17s}{note}")
print("\n  PRECONDITION 1 -- TWO INDEPENDENT DESIGNABLE LOOPS (cycle rank >= 2).  B5 fails it,")
print("  and on B5 the entire object of this corpus is undefined.  S4 in its own words: 'This is")
print("  the one place in S4 where topology, and nothing else, decides an outcome.'")
print("  PRECONDITION 2 -- gamma_C FLAT (b1 >= 1).  B2 fails it.  On B2 both holonomies are")
print("  curvatures, so W-01's advertised virtue -- 'it distinguishes curvature from flat")
print("  holonomy' -- has no referent there, for a reason that has nothing to do with class")
print("  occupancy: B2's class multiset is K1's exactly, {01:2,10:2,11:1}.")
print("  Lane D row 3.4 rules the topology-inertness claim 'UNTESTABLE BY CONSTRUCTION ... no")
print("  experiment on any carrier can bear on this row'.  B5 and B2 both bear on it: topology")
print("  decides whether the functional has an argument, and what the argument MEANS.")
print("  Every CARRIER_INDEPENDENT verdict in the scope table is conditional on both, and no")
print("  row of the table states either.  Neither B5 nor B2 appears anywhere in the table.")
print("\n  AND THE LOOP LENGTHS, WHICH THE SAME TABLE DETERMINES AND NO ROW MENTIONS:")
print(f"  {'carrier':20s} {'|gamma_F| = |10|+|11|':>22s} {'|gamma_C| = |01|+|11|':>22s} {'equal?':>8s}")
for nm, mult in CARRIERS.items():
    lF = mult.get('10', 0) + mult.get('11', 0)
    lC = mult.get('01', 0) + mult.get('11', 0)
    print(f"  {nm:20s} {lF:22d} {lC:22d} {str(lF == lC):>8s}")
print("  EVERY carrier the corpus has ever RUN anything on has two loops of EQUAL length, so")
print("  'one circuit of each' and 'equal edge time' coincide there.  B0b is the first carrier")
print("  in the corpus on which they do not.  See leg 3.")

# =============================================================== 2F  the reduction step
print("\n" + "=" * 104)
print("== 2F  THE STEP LANE D's HEADLINE RESTS ON, DECLARED AS THE IDENTITY IT IS ==")
print("=" * 104)
print("  Lane D's headline: 'the formation functional never sees anything but the class-weight")
print("  vector'.  The step that makes that true is W-03's sec2 sentence 'every vertex phase of")
print("  s cancels', which the register records and marks FALSE the moment the holonomy is not")
print("  a scalar.  Lane D never exhibits it at four classes -- every one of its legs starts")
print("  FROM a weight 4-vector.  Exhibited here on B0b, and DECLARED A ZERO-VARIABLE CONTROL:")
rng = np.random.default_rng(20260816)
vc = vertex_classes(CARRIERS['B0b torus meeting '])
a = np.array([1.0 if c[0] == '1' else 0.0 for c in vc])
b = np.array([1.0 if c[1] == '1' else 0.0 for c in vc])
worst = 0.0
for _ in range(200):
    s = rng.normal(size=len(vc)) + 1j * rng.normal(size=len(vc))
    s /= np.linalg.norm(s)
    f, c = rng.uniform(-np.pi, np.pi, 2)
    p = np.zeros(4)
    for j, cl in enumerate(vc):
        p[ORDER.index(cl)] += abs(s[j]) ** 2
    for k in (1, 2, 7):
        lhs = np.sum(np.conj(np.exp(1j * f * a * k)) * np.exp(1j * c * b * k) * np.abs(s) ** 2)
        rhs = p[0] + p[1] * np.exp(-1j * f * k) + p[2] * np.exp(1j * c * k) + p[3] * np.exp(1j * (c - f) * k)
        worst = max(worst, abs(lhs - rhs))
print(f"  200 random complex states x 3 values of k on B0b: max |<M_dF^k s, M_c^k s> - P(u^k,v^k)|")
print(f"  = {worst:.3e}.  IT COULD NOT HAVE FAILED -- the phases cancel term by term because the")
print("  transport is scalar per fibre.  THAT VOIDS IT AS A CONTROL AND LEAVES IT A ONE-LINE")
print("  IDENTITY, which is the right status for it; what is not right is that lane D's central")
print("  claim rests on it and no row of the table carries it.  Its hypothesis is SCALAR")
print("  FIBRE-WISE TRANSPORT (W-04 ERR-1 / W-06's N4 correction), not class occupancy -- so the")
print("  headline should read 'carrier-independent GIVEN the scalar-transport convention'.")
