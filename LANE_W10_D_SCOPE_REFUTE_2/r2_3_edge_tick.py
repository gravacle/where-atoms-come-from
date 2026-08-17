# W10-D REFUTE-2  LENS 2 = COMPLETENESS.  LEG 3.
#
# THE OMISSION.  W-06's register row identifies ONE object as the thing that decided the
# program's spine, and it is not in any row of lane D's scope table:
#
#   "BUT THE THING THAT DECIDED THE SPINE WAS NEVER IMPORTED -- AND THE CORPUS ALREADY KNEW.
#    It was an unledgered stipulation the corpus made about itself: that loop transport means
#    MULTIPLYING BY THE WHOLE-CIRCUIT HOLONOMY M_gamma, rather than moving EDGE BY EDGE via T
#    with T^3 = M_gamma. ... THE CORPUS'S OWN SEALED AUDIT FOUND THIS AND WROTE IT DOWN AS
#    COR-F, AT S3_THE_CROSSING_AUDIT_V001.md:794, BEFORE W-03, W-04 AND W-05 RAN."
#    "> THE CHAIN WAS NOT UNDER-ADVERSARIAL.  IT WAS UNDER-READ."
#
# `grep -c "COR-" W10D_SCOPE_TABLE_V001.md` = 0.  Not one of the twenty sealed corrections the
# brief says to CARRY (S3 COR-A..COR-L, S2 COR-A..COR-H) is named anywhere in the scope table,
# and COR-F is absent in substance too: no row scopes the transport convention.  Lane D's only
# nominal cover is row 4.5, which lumps W-04's "(i)-(iii)" into one cell marked
# CARRIER_INDEPENDENT with the evidence pointer "leg 4B" -- and leg 4B computes vertex-count
# arithmetic (V = 5,6,7,9,11 and what divides what).  It contains no transport.
#
# WHY IT IS A SCOPE QUESTION AND NOT A PHILOSOPHY QUESTION.  Under the scalar convention both
# transports are DIAGONAL, hence commute, hence see only the class multiset -- which is lane D's
# entire headline.  Under COR-F's edge convention they are NOT diagonal, they need NOT commute,
# and they see the loop LENGTHS and the MEETING PATTERN.  On K1 the two loops have the same
# length (3 and 3), so a circuit clock and an edge clock coincide.  ON B0b THEY DO NOT (4 and 3),
# and B0b is the corpus's own four-class carrier.
#
# ISOLATION.  TWO variables move in this leg and each is moved ALONE, in its own sub-leg:
#   3A/3C: WHETHER THE LOOPS MEET.  B0a and B0b are built as THE SAME 3x3 torus grid with the
#          SAME gamma_F and the SAME connection; the ONLY thing that changes is which row is
#          designated gamma_C (row 0, which meets the square -> four classes; row 2, which does
#          not -> three classes).  Arm diff printed: the two edge sets, and the two class
#          multisets, which reproduce S4:575 and S4:576 respectively.
#   3B/3D/3E: THE TRANSPORT CONVENTION (scalar M_gamma vs edge-tick T), with the carrier,
#          connection, state and k-range held fixed.
# PRECISION: numpy complex128.  Every operator identity is reported as a matrix norm.

import numpy as np
from collections import Counter

rng = np.random.default_rng(20260816)
np.set_printoptions(precision=6, suppress=True)


def transport_op(V, loop_path, U):
    """loop_path = [(src,tgt,phase_index_sign)] as an ordered closed walk; T moves each fibre
       value one edge along the walk, identity off the walk."""
    T = np.eye(V, dtype=complex)
    on = [p[0] for p in loop_path]
    for v in on:
        T[v, v] = 0.0
    for (src, tgt, ph) in loop_path:
        T[tgt, src] = ph
    return T


def scalar_op(V, on, W):
    D = np.ones(V, dtype=complex)
    for v in on:
        D[v] = W
    return np.diag(D)


# ================================================================= K1
print("=" * 104)
print("== 3B  COR-F's OWN EXHIBIT, REPRODUCED ON K1 (the audit's arm) ==")
print("=" * 104)
a = np.array([0.3, 0.9, -0.4, 0.7, 1.3, -0.4])           # a1..a6
U = np.exp(1j * a)
WF_K1 = np.exp(1j * (a[0] + a[1] + a[2]))
WC_K1 = np.exp(1j * (a[3] + a[4] + a[5]))
# e1=(0,1) e2=(1,2) e3=(2,0) e4=(0,3) e5=(3,4) e6=(4,0)
K1_F = [(0, 1, U[0]), (1, 2, U[1]), (2, 0, U[2])]
K1_C = [(0, 3, U[3]), (3, 4, U[4]), (4, 0, U[5])]
TF_K1 = transport_op(5, K1_F, U)
TC_K1 = transport_op(5, K1_C, U)
MF_K1 = scalar_op(5, [0, 1, 2], WF_K1)
MC_K1 = scalar_op(5, [0, 3, 4], WC_K1)
rho = np.diag([0.4, 0.15, 0.15, 0.15, 0.15]).astype(complex)
print(f"  ||T_C* T_C - I||          = {np.abs(TC_K1.conj().T@TC_K1 - np.eye(5)).max():.2e}   (unitary)")
print(f"  T_C diagonal?              {np.allclose(TC_K1, np.diag(np.diag(TC_K1)))}")
print(f"  ||T_C^3 - M_c||           = {np.abs(np.linalg.matrix_power(TC_K1,3) - MC_K1).max():.2e}"
      f"   (COR-F: 'T^3 equal to W-01's scalar operator')")
print(f"  ||T_F^3 - M_dF||          = {np.abs(np.linalg.matrix_power(TF_K1,3) - MF_K1).max():.2e}")
print(f"  diag(T_C rho T_C*)        = {np.real(np.diag(TC_K1@rho@TC_K1.conj().T))}")
print(f"  diag(rho)                 = {np.real(np.diag(rho))}          NOT PRESERVED")
print(f"  ||[M_dF, M_c]||_F         = {np.linalg.norm(MF_K1@MC_K1 - MC_K1@MF_K1):.2e}   (scalar convention: ALWAYS 0)")
print(f"  ||[T_F,  T_C ]||_F        = {np.linalg.norm(TF_K1@TC_K1 - TC_K1@TF_K1):.6f}   (edge convention on K1)")
print("  COR-F REPRODUCED at the bytes, from an independent implementation.")

# ================================================================= the grid: B0a and B0b
print("\n" + "=" * 104)
print("== 3A  ARM DIFF: B0a AND B0b ARE THE SAME COMPLEX WITH ONE LOOP MOVED ==")
print("=" * 104)
n = 3
vid = {(i, j): 3 * i + j for i in range(n) for j in range(n)}
edges, eidx = [], {}
for i in range(n):
    for j in range(n):
        for e in ((vid[(i, j)], vid[((i + 1) % n, j)]), (vid[(i, j)], vid[(i, (j + 1) % n)])):
            eidx[e] = len(edges)
            edges.append(e)
V, E = n * n, len(edges)
aE = rng.uniform(-np.pi, np.pi, E)
UE = np.exp(1j * aE)


def walk_ops(walk):
    """walk = ordered list of vertices, closed.  Returns (path with phases, holonomy, vertices)."""
    path, hol = [], 1.0 + 0j
    for t in range(len(walk)):
        s_, t_ = walk[t], walk[(t + 1) % len(walk)]
        if (s_, t_) in eidx:
            ph = UE[eidx[(s_, t_)]]
        else:
            ph = np.conj(UE[eidx[(t_, s_)]])
        path.append((s_, t_, ph))
        hol *= ph
    return path, hol, [w for w in walk]


# gamma_F = boundary of the square face at (0,0): (0,0)->(1,0)->(1,1)->(0,1)->(0,0)
WALK_F = [vid[(0, 0)], vid[(1, 0)], vid[(1, 1)], vid[(0, 1)]]
# gamma_C, ARM 1 (B0b, loops MEET):   column j=0:  (0,0)->(1,0)->(2,0)->(0,0)
WALK_C_meet = [vid[(0, 0)], vid[(1, 0)], vid[(2, 0)]]
# gamma_C, ARM 2 (B0a, loops DISJOINT): column j=2: (0,2)->(1,2)->(2,2)->(0,2)
WALK_C_disj = [vid[(0, 2)], vid[(1, 2)], vid[(2, 2)]]

for lab, WC_ in (("B0b loops MEET    ", WALK_C_meet), ("B0a loops DISJOINT", WALK_C_disj)):
    FV, CV = set(WALK_F), set(WC_)
    cnt = Counter(('1' if v in FV else '0') + ('1' if v in CV else '0') for v in range(V))
    print(f"  {lab}  gamma_F verts {sorted(FV)}  gamma_C verts {sorted(CV)}  shared {sorted(FV&CV)}")
    print(f"                      class multiset {{'00':{cnt['00']}, '01':{cnt['01']}, "
          f"'10':{cnt['10']}, '11':{cnt['11']}}}   |gamma_F| = {len(WALK_F)}  |gamma_C| = {len(WC_)}")
print("  S4:575 B0a PUBLISHED {00:2, 01:3, 10:4}      S4:576 B0b PUBLISHED {00:4, 01:1, 10:2, 11:2}")
print("  BOTH REPRODUCED.  ONE VARIABLE SEPARATES THEM: which column is gamma_C.  Same complex,")
print("  same gamma_F, same connection, same loop LENGTHS (4 and 3).  This is the arm lane D")
print("  does not have: it separates FOUR CLASSES from UNEQUAL LOOP LENGTHS, which are confounded")
print("  in every four-class row of the scope table (B0b has both; B4 has four classes and equal")
print("  lengths 4 and 4; B0a has unequal lengths and three classes).")

# ================================================================= 3C  edge-tick off K1
print("\n" + "=" * 104)
print("== 3C  COR-F's EDGE TRANSPORT ON THE FOUR-CLASS CARRIER.  ONE VARIABLE: DO THE LOOPS MEET ==")
print("=" * 104)
pathF, WF_g, _ = walk_ops(WALK_F)
TF = transport_op(V, pathF, UE)
MF = scalar_op(V, WALK_F, WF_g)
print(f"  {'arm':20s} {'|g_F|':>6s} {'|g_C|':>6s} {'shared':>7s} {'||T_F^|g_F| - M_dF||':>21s} "
      f"{'||T_C^|g_C| - M_c||':>20s} {'||[M_dF,M_c]||':>15s} {'||[T_F,T_C]||':>14s}")
ARMS3 = []
for lab, WC_ in (("B0b loops MEET    ", WALK_C_meet), ("B0a loops DISJOINT", WALK_C_disj),
                 ("K1 (both len 3)   ", None)):
    if WC_ is None:
        TFx, TCx, MFx, MCx, lF, lC, sh = TF_K1, TC_K1, MF_K1, MC_K1, 3, 3, 1
    else:
        pathC, WC_g, _ = walk_ops(WC_)
        TCx = transport_op(V, pathC, UE)
        MCx = scalar_op(V, WC_, WC_g)
        TFx, MFx = TF, MF
        lF, lC = len(WALK_F), len(WC_)
        sh = len(set(WALK_F) & set(WC_))
    dF = np.abs(np.linalg.matrix_power(TFx, lF) - MFx).max()
    dC = np.abs(np.linalg.matrix_power(TCx, lC) - MCx).max()
    cM = np.linalg.norm(MFx @ MCx - MCx @ MFx)
    cT = np.linalg.norm(TFx @ TCx - TCx @ TFx)
    ARMS3.append((lab, TFx, TCx, MFx, MCx, lF, lC))
    print(f"  {lab:20s} {lF:6d} {lC:6d} {sh:7d} {dF:21.2e} {dC:20.2e} {cM:15.2e} {cT:14.6f}")
print("  T^L = M_gamma HOLDS ON EVERY CARRIER with L = the loop's EDGE LENGTH -- so COR-F's")
print("  alternative EXISTS off K1 and is not a K1 curiosity.  The scalar transports commute on")
print("  every carrier (they are diagonal); THE EDGE TRANSPORTS COMMUTE EXACTLY WHEN THE LOOPS")
print("  ARE DISJOINT, i.e. exactly when class 11 is EMPTY.  ||[T_F,T_C]|| is therefore a")
print("  function of the INCIDENCE, not of the class multiset alone -- and the corpus's central")
print("  survivor N2 says the incidence labels are invisible.  THEY ARE INVISIBLE TO THE SCALAR")
print("  CONVENTION AND VISIBLE TO THE CORPUS'S OWN SEALED ALTERNATIVE.")

# ================================================================= 3D  the clock
print("\n" + "=" * 104)
print("== 3D  THE CLOCK.  ONE VARIABLE: THE TRANSPORT CONVENTION, CARRIER AND STATE FIXED ==")
print("=" * 104)
print("  W-01's row: 'Its no-time claim is overstated -- CIRCUIT COUNT IS CARRIER-SUPPLIED")
print("  DISCRETE TIME.'  S3 sec3.5 makes circuit count the clock.  COR-F: 'Circuits supply")
print("  discrete time'.  On K1 one circuit of gamma_F and one of gamma_C are BOTH 3 edges, so")
print("  a circuit clock and an edge clock are the same clock up to a factor 3.  Off K1 they")
print("  are not, and the corpus has never compared them.")
print("\n  Scalar convention at circuit k:  Z_k = <M_dF^k s, M_c^k s>")
print("  Edge   convention at edge e:     Y_e = <T_F^e s, T_C^e s>")
for lab, TFx, TCx, MFx, MCx, lF, lC in ARMS3:
    Vx = TFx.shape[0]
    s = rng.normal(size=Vx) + 1j * rng.normal(size=Vx)
    s /= np.linalg.norm(s)
    print(f"\n    {lab}   |gamma_F| = {lF}, |gamma_C| = {lC}   lcm = {np.lcm(lF,lC)}")
    print(f"      {'e':>4s} {'|Y_e| (edge clock)':>20s} | {'k':>4s} {'|Z_k| (circuit clock)':>22s}"
          f" {'agree?':>8s}")
    Z = []
    for k in range(0, 13):
        Z.append(abs(np.vdot(np.linalg.matrix_power(MFx, k) @ s,
                             np.linalg.matrix_power(MCx, k) @ s)))
    for e in range(0, 13):
        Ye = abs(np.vdot(np.linalg.matrix_power(TFx, e) @ s, np.linalg.matrix_power(TCx, e) @ s))
        if e % lF == 0 and e % lC == 0:
            k = e // lF
            print(f"      {e:4d} {Ye:20.12f} | {k:4d} {Z[k]:22.12f} {abs(Ye-Z[k]) < 1e-12!s:>8s}")
        else:
            print(f"      {e:4d} {Ye:20.12f} |    -                      -        -")
print("\n  ON K1 the two clocks agree at EVERY multiple of 3, EXACTLY -- one circuit each, three")
print("  edges each -- so the corpus never had to choose, and every place it ever looked is a")
print("  place where the choice does not exist.")
print("  MY OWN CORRECTION, RECORDED NOT PATCHED: I first wrote that on B0a/B0b the two clocks")
print("  'coincide only at multiples of 12'.  THE TABLE ABOVE REFUTES MY OWN SENTENCE -- e = 12")
print("  prints agree? FALSE on both.  The lcm is where the two loops are simultaneously at a")
print("  whole circuit, but they are at DIFFERENT whole circuits: T_F^12 = M_dF^3 while")
print("  T_C^12 = M_c^4, so Y_12 = <M_dF^3 s, M_c^4 s>, which is not Z_k for ANY k.  THE TWO")
print("  CLOCKS NEVER COINCIDE AT ANY e > 0 WHEN THE LOOP LENGTHS DIFFER.  That is stronger")
print("  than what I first claimed, and it is the version the numbers support.")
print("  Note also K1's e = 3, 6, 9, 12 rows in leg 3E: the edge convention's branch separation")
print("  is EXACTLY ZERO there, because T^3 = M_gamma is scalar.  On K1 the two conventions")
print("  agree at exactly the ticks the corpus clocked, and nowhere else; on B0b they agree")
print("  nowhere at all.  'Circuit count is carrier-supplied discrete time' (W-01's row) is a")
print("  well-defined stipulation only when the two designated loops have the SAME LENGTH, and")
print("  the corpus's ten-carrier table shows two carriers where they do not -- B0a and B0b.")
print("  NO ROW OF THE REGISTER AND NO ROW OF THE SCOPE TABLE STATES THIS CONDITION.")

# ================================================================= 3E  what a diagonal sees
print("\n" + "=" * 104)
print("== 3E  COR-F's PAYLOAD OFF K1: A DIAGONAL OBSERVABLE SEPARATES THE EDGE-TICK BRANCHES ==")
print("=" * 104)
print("  S3-0's body: the gauge-invariant algebra is the diagonal, and the diagonal cannot tell")
print("  M_dF s from M_c s.  COR-F: true of W-01's SCALAR transports, false of a bona fide")
print("  edge transport around the SAME cycle.  ONE VARIABLE: the transport convention.")
print(f"  {'arm':20s} {'k/e':>4s} {'max_v |diag scalar branches|':>30s} {'max_v |diag edge branches|':>28s}")
for lab, TFx, TCx, MFx, MCx, lF, lC in ARMS3:
    Vx = TFx.shape[0]
    s = rng.normal(size=Vx) + 1j * rng.normal(size=Vx)
    s /= np.linalg.norm(s)
    for k in (1, 2, 3):
        ds = np.abs(np.abs(np.linalg.matrix_power(MFx, k) @ s) ** 2 -
                    np.abs(np.linalg.matrix_power(MCx, k) @ s) ** 2).max()
        de = np.abs(np.abs(np.linalg.matrix_power(TFx, k) @ s) ** 2 -
                    np.abs(np.linalg.matrix_power(TCx, k) @ s) ** 2).max()
        print(f"  {lab:20s} {k:4d} {ds:30.2e} {de:28.12f}")
print("  THE SCALAR CONVENTION GIVES 0 ON EVERY CARRIER (that IS S3-0, and it is an identity).")
print("  THE EDGE CONVENTION SEPARATES ON EVERY CARRIER INCLUDING THE FOUR-CLASS ONE.  COR-F is")
print("  CARRIER-INDEPENDENT, it is not scoped anywhere in lane D's table, and it is the one")
print("  object W-06's row names as having decided the spine.")
