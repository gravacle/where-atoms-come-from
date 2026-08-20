"""T48-C driver -- PROBE 3: THE WORLD VENUE COMPUTATION (O-58 N2, world tier).

QUESTION (the judge's N2, verbatim hypothesis under test, never assumed):
  mu_c = 1/deg on every computed venue, and mu = 1/deg is exactly the stochastic
  measure-conserving normalization of one writer step; if Gamma's writer kernel
  conserves measure, mu = mu_c IDENTICALLY -- masslessness is measure conservation,
  not tuning.
TASK: on the census access geometry's deg-6 grain lattice, compute the INDUCED
per-link amplitude of one writer step from the world ensemble the model actually
has -- the (iv') energy-conserving dilation writers on written media (C-71/C-72
occupancy and orientation encodings), with the activation convention of
model/project_model.py (corrected form):
    gu = f0 exp(-E_b/kT)         escape FROM the metastable well
    gd = f0 exp(-(E_b+dE)/kT)    reverse, INTO the metastable well
Does detailed balance at temperature T induce the measure-conserving normalization
1/deg?  Only at dE = 0?  Never?  dE is the physical bias dial; if measure
conservation breaks, compute WHERE the induced amplitude lands relative to
mu_c = 1/6 (mass gap ln(mu_c/mu)) as a function of dE/kT.

AMBIGUITY GUARD (binding): more than one honest construction of "the writer
ensemble" exists.  ALL are declared here, BEFORE measurement, and ALL are computed;
none is silently preferred.  Per-attempt probabilities (the model's own attempt
clock f0; one attempt per 1/f0):
    u := exp(-E_b/kT)      per-attempt flip probability, metastable -> stable
    v := u * b             per-attempt flip probability, stable -> metastable,
    b := exp(-dE/kT)       the Boltzmann factor of the well asymmetry (detailed
                           balance: v/u = b; gated in S2 from the model's own
                           two-state kernel, not asserted)
  E1  TRANSPORT (iv' read literally): one writer step = ONE energy-conserving
      dilation event that transports the written unit across one grain boundary
      (erase at x releases dE, write at y costs dE; net configuration-energy
      change COMPUTED in S2).  Per-link amplitude = the single event's activation
      factor a; a is swept over BOTH readings a = u and a = u*b, because the
      composite saddle is not pinned by the model -- the result must be (and is
      gated) independent of that choice.  Stay = 1 - 6a (nothing happens this
      attempt).  This is the mediation that leaves no written trail (the
      R1-transparent reading).
  E2  TRAIL WITH RETREAT (the raw two rates, strings counted WITH backtracking,
      exactly H1's ensemble): tip state (grain, last direction); per attempt the
      tip writes a fresh neighbor (v each, 5 directions), or the tip grain erases
      and the string retreats (u, 1 direction -- the backtracking link), or
      nothing.  Measure kept or lost is COMPUTED, not declared.
  E3a TRAIL WITH DECAY, H1-matched counting (extension-only strings, the walk
      ensemble of the sealed T-44 coupling with the model's own erase channel as
      absorption): per attempt, extend to any of the 6 neighbors (v each) or the
      tip grain decays and the string leaves the ensemble (probability u = the
      model's own erase probability), or nothing.
  E3b TRAIL WITH DECAY, non-backtracking counting: extend fresh (v each, 5
      directions) or decay out (u), or nothing; criticality reference for this
      counting = the venue's own directed-edge (non-backtracking) operator,
      earned by the same row-sum instrument.

DECLARED CONSTANTS AND GATES (fixed before measurement; D-8: every verdict a
computed boolean; exact-equality gates have tolerance ZERO on Fractions):
  N = 8 venue side; K_WRAP = 12 (wrap identity, every cell, every k <= 12);
  K_AXIS = 10 (axis-split vs DP, all endpoints); K_BOX = 6, BOX_B = 7
  (block-interior vs free); K_MEAS = 12 (exact measure evolution);
  K_SUM = 12 (propagator partial sums).
  U_SAMPLES = [1/20, 1/100]; B_SAMPLES = [1, 9/10, 3/4, 1/2, 1/4, 1/10]
  (declared rational sample points of the Boltzmann factors; dE/kT = -ln b exists
  only in the labeled COMPARISON section).  b = 1 is the orientation encoding's
  field-free degenerate wells (C-71 CoCrPt, dE = 0); b < 1 is the occupancy
  encoding (C-71 NAND, dE > 0).
  MU_C_REF := 1/deg_venue with deg_venue COMPUTED in G01 (comparison-to-sealed:
  PUBLISHED_CONVENTIONS' deg-6, T44-B's mu_c = 1/6).
  SEALED_MU_SUB_ROW := 1/8 (T44-B's declared subcritical sweep row; comparison).
  COMPARISON float tolerance: 1e-12 (labeled; after results).
VERDICT SEMANTICS (declared): "measure conserved" = every row sum exactly 1;
"induced amplitude" = the m in I - W = c(I - sum_e m_e S_e), extracted entrywise
and gated for uniformity across cells and directions (uniformity EARNED or shown
to fail, never inserted); "AT criticality" = m uniform and equal to the venue's
computed 1/deg (E1/E2/E3a) or 1/deg_NB (E3b); "mass gap" = ln(mu_c/mu), reported
as the exact rational ratio mu_c/mu on the measurement path.
"""
import json
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t48c_lib import (DIRS, OPP, torus3, bfs, l1_wrap, torus_counts, dp3, N3_free,
                      kernel_pos, kernel_edge, nb_edge_adjacency, row_sums, col_sums,
                      apply_forward, extract_pos, extract_edge, op_identity_pos,
                      drift_per_state)

N = 8
K_WRAP = 12
K_AXIS = 10
K_BOX = 6
BOX_B = 7
K_MEAS = 12
K_SUM = 12
U_SAMPLES = [F(1, 20), F(1, 100)]
B_SAMPLES = [F(1), F(9, 10), F(3, 4), F(1, 2), F(1, 4), F(1, 10)]
SEALED_DEG = 6              # comparison-to-sealed only
SEALED_MU_SUB_ROW = F(1, 8)  # comparison-to-sealed only
CMP_TOL = 1e-12             # labeled COMPARISON tolerance (floats, after results)

GATES = []


def gate(gid, ok, desc, val=""):
    GATES.append(bool(ok))
    line = "[%s] %s  %s" % (gid, "PASS" if ok else "FAIL", desc)
    if val != "":
        line += "  :: %s" % (val,)
    print(line)
    return bool(ok)


print(__doc__)
print("=" * 78)
print("SECTION 1 -- VENUE REBUILD GATES (identities gated against brute force"
      " BEFORE any kernel is built)")
print("=" * 78)

cells, idx, nbr = torus3(N)
NC = len(cells)

# G01 row sums / degree
degs = {len(set(row)) for row in nbr} | {len(row) for row in nbr}
g01 = (len(degs) == 1)
DEG = nbr and len(nbr[0])
gate("G01", g01 and all(len(set(row)) == DEG for row in nbr),
     "venue row sums: every grain has the same face-degree, no multi-edges",
     "deg_venue = %d (sealed conventions say %d: %s)" % (DEG, SEALED_DEG, DEG == SEALED_DEG))
MU_C_REF = F(1, DEG)

# G02 BFS == L1 wrap distance
d_bfs = bfs(nbr, 0)
g02 = all(d_bfs[i] == l1_wrap(cells[i], N) for i in range(NC))
gate("G02", g02, "BFS distance == earned L1 (wrap) separation on every grain")

# G03 block-interior identity
free = dp3(K_WRAP)
box = dp3(K_BOX, inside=lambda p: max(abs(c) for c in p) <= BOX_B)
g03 = all(box[key] == free[key] for key in box if key[0] <= K_BOX)
gate("G03", g03, "block-interior DP == free-sector DP, k <= %d, box |x|inf <= %d"
     % (K_BOX, BOX_B))

# G04 axis-split identity vs brute force
g04 = True
for (k, x, y, z), c in free.items():
    if k <= K_AXIS:
        if N3_free(k, x, y, z) != c:
            g04 = False
            break
g04 = g04 and all(N3_free(k, x, y, z) == 0
                  for k in range(K_AXIS + 1)
                  for (x, y, z) in [(k + 1, 0, 0), (0, k + 1, 0), (0, 0, k + 2)])
gate("G04", g04, "axis-split identity == brute-force DP, k <= %d, all endpoints"
     " (+ out-of-range zeros)" % K_AXIS)

# G05 torus wrap identity: fold free DP counts mod N and compare with venue transfer counts
tc = torus_counts(nbr, idx[(0, 0, 0)], K_WRAP)
g05 = True
for k in range(K_WRAP + 1):
    folded = [0] * NC
    for (kk, x, y, z), c in free.items():
        if kk == k:
            folded[idx[(x % N, y % N, z % N)]] += c
    if folded != tc[k]:
        g05 = False
        break
gate("G05", g05, "torus wrap identity: venue transfer counts == wrap-folded free DP,"
     " every cell, every k <= %d" % K_WRAP)

# G06 octahedral symmetry of the venue (earns direction-equivalence as venue symmetry)
def relabel(perm_neg):
    m = {}
    for c in cells:
        m[idx[c]] = idx[perm_neg(c)]
    return m

gens = [lambda c: (c[1], c[0], c[2]),
        lambda c: ((-c[0]) % N, c[1], c[2]),
        lambda c: (c[1], c[2], c[0])]
g06 = True
for gnum, g in enumerate(gens):
    m = relabel(g)
    for i in range(NC):
        if sorted(m[j] for j in nbr[i]) != sorted(nbr[m[i]]):
            g06 = False
gate("G06", g06, "octahedral generators (swap, reflect, cycle) preserve the venue"
     " adjacency: all 6 directions are venue-equivalent")

print()
print("=" * 78)
print("SECTION 2 -- WRITER ENERGETICS: earning the per-link amplitudes from the")
print("model's own algebra (configuration-energy changes COMPUTED, in units of dE)")
print("=" * 78)

# configuration instrument (C-71 occupancy count): E_config = dE * (#written grains).
# Delta-counts are computed exact integers on explicit written-set patterns; the
# amplitude table then follows from the model's corrected activation convention:
# flip INTO the metastable well (dN = +1) pays E_b + dE -> amplitude v = u*b;
# flip FROM the metastable well (dN = -1) pays E_b -> amplitude u; a single
# energy-conserving composite (dN = 0) pays its one declared saddle -> amplitude a
# (swept over both readings in E1).
def config_count(written):
    return len(written)

x0 = idx[(0, 0, 0)]
dN_transport, dN_write, dN_erase = [], [], []
# an explicit written trail ending at the tip x_tip (concrete pattern, not symbolic)
trail = [idx[(0, 0, 0)], idx[(1, 0, 0)], idx[(2, 0, 0)], idx[(2, 1, 0)]]
x_tip = trail[-1]
for d in range(6):
    y = nbr[x0][d]
    # transport: written set {x0} -> {y} (erase behind, write ahead, one event)
    dN_transport.append(config_count({y}) - config_count({x0}))
    yt = nbr[x_tip][d]
    if yt not in trail:
        # trail write: S -> S + {y}
        dN_write.append(config_count(set(trail) | {yt}) - config_count(set(trail)))
# erase: the tip un-writes (any trail, tip removed)
dN_erase.append(config_count(set(trail) - {x_tip}) - config_count(set(trail)))
g07 = all(v == 0 for v in dN_transport) and len(dN_transport) == 6
gate("G07", g07, "TRANSPORT step: configuration-count change dN == 0 in EVERY"
     " direction (dE cancels exactly; earns E1's direction-uniform amplitude)",
     "dN per direction = %s" % (dN_transport,))
g08 = (all(v == +1 for v in dN_write) and len(dN_write) == 5
       and all(v == -1 for v in dN_erase))
gate("G08", g08, "TRAIL step: write dN == +1 (every fresh direction), erase dN == -1"
     " (the D-15 nonzero beside G07's zero)",
     "write %s, erase %s" % (dN_write, dN_erase))

# G09 detailed balance from the model's own two-state kernel (per (u,b) sample):
# K2 = [[1-u, u], [v, 1-v]] rows = (from metastable, from stable); stationary pi
# solved exactly; gate pi K2 == pi and pi_meta/pi_stable == v/u == b.
g09 = True
for u in U_SAMPLES:
    for b in B_SAMPLES:
        v = u * b
        pi = (v / (u + v), u / (u + v))
        piK = (pi[0] * (1 - u) + pi[1] * v, pi[0] * u + pi[1] * (1 - v))
        if piK != pi or pi[0] / pi[1] != b or v / u != b:
            g09 = False
gate("G09", g09, "detailed balance COMPUTED from the model's two-state kernel:"
     " stationary pi exact, pi_meta/pi_stable == v/u == b at every (u,b) sample")

print()
print("=" * 78)
print("SECTION 3 -- THE KERNELS ON THE VENUE: measure audit and induced per-link")
print("amplitude, every declared ensemble, every (u,b) sample point")
print("=" * 78)

results = {"E1": [], "E2": [], "E3a": [], "E3b": []}

# ---------------------------------------------------------------- E1 TRANSPORT
g10 = g11 = g12 = g13 = g14 = g15 = g16 = True
e1_mus = set()
for u in U_SAMPLES:
    for b in B_SAMPLES:
        for aname, a in (("a=u", u), ("a=u*b", u * b)):
            W = kernel_pos(nbr, 1 - DEG * a, a)
            rs = set(row_sums(W))
            nonneg = all(vv >= 0 for r in W for vv in r.values())
            if rs != {F(1)} or not nonneg:
                g10 = False
            cs, ms = extract_pos(W, nbr)
            if len(cs) != 1 or len(ms) != 1:
                g11 = False
            c = next(iter(cs))
            m = next(iter(ms))
            e1_mus.add(m)
            if m != MU_C_REF:
                g12 = False
            if not op_identity_pos(W, nbr, c, m):
                g13 = False
            results["E1"].append(dict(u=str(u), b=str(b), variant=aname,
                                      row_sum="1", c=str(c), mu_induced=str(m)))
# measure evolution + divergence witness at one representative sample (exact)
u, b = U_SAMPLES[0], B_SAMPLES[2]
W = kernel_pos(nbr, 1 - DEG * u, u)
vec = [F(0)] * NC
vec[x0] = F(1)
tot_sum = F(0)
for k in range(K_MEAS + 1):
    tot = sum(vec)
    tot_sum += tot
    if tot != 1:
        g14 = False
    if k < K_MEAS:
        vec = apply_forward(W, vec)
if tot_sum != K_MEAS + 1:
    g15 = False
im1 = [1 - sum(W[i].values()) for i in range(NC)]  # (I-W)1 by row sums
if any(vv != 0 for vv in im1):
    g16 = False
gate("G10", g10, "E1 kernel: nonnegative, EVERY row sum exactly 1 -- the transport"
     " writer CONSERVES MEASURE at every (u,b) and both amplitude readings")
gate("G11", g11, "E1 induced amplitude UNIFORM across all cells and all 6 directions"
     " (uniformity EARNED: G06 venue symmetry + G07 computed dE-cancellation)")
gate("G12", g12, "E1 induced per-link amplitude == MU_C_REF = 1/deg_venue for EVERY"
     " (u,b) and BOTH amplitude readings",
     "induced mu set = {%s} == {%s}" % (", ".join(sorted(str(m) for m in e1_mus)), MU_C_REF))
gate("G13", g13, "E1 operator identity entrywise: I - W == c(I - m A) on the venue")
gate("G14", g14, "E1 exact measure evolution from a point source: total == 1 at every"
     " step k <= %d" % K_MEAS)
gate("G15", g15, "E1 propagator partial sums == K+1 (LINEAR divergence: the string sum"
     " sits AT its critical point; massless witness)", "sum = %s" % tot_sum)
gate("G16", g16, "E1 criticality signature: (I - W) 1 == 0 exactly (resolvent singular"
     " AT the writer's own normalization -- the T44-B mu_c signature)")

# ---------------------------------------------------------------- E2 TRAIL WITH RETREAT
g17 = g18 = g19 = g20 = g21 = True
NE = NC * 6
for u in U_SAMPLES:
    for b in B_SAMPLES:
        v = u * b
        W = kernel_edge(nbr, 1 - 5 * v - u, v, u)
        rs = set(row_sums(W))
        if rs != {F(1)}:
            g17 = False
        csum = set(col_sums(W, NE))
        if csum != {F(1)}:
            g18 = False
        cs, mf, mb, tots = extract_edge(W, nbr)
        if len(cs) != 1 or len(mf) != 1 or len(mb) != 1 or tots != {F(1)}:
            g19 = False
        mfv, mbv = next(iter(mf)), next(iter(mb))
        uniform = (mfv == mbv)
        if uniform != (b == 1):
            g20 = False
        if b == 1 and (mfv != MU_C_REF or mbv != MU_C_REF):
            g20 = False
        dr = drift_per_state(W, nbr)
        avg = [F(0)] * 3
        for i in range(NC):
            for d in range(6):
                s = i * 6 + d
                expect = tuple((v - u) * F(DIRS[d][ax]) for ax in range(3))
                if dr[s] != expect:
                    g21 = False
                for ax in range(3):
                    avg[ax] += dr[s][ax]
        if any(vv != 0 for vv in avg):
            g21 = False
        results["E2"].append(dict(u=str(u), b=str(b), m_fresh=str(mfv), m_back=str(mbv),
                                  redistribution_total="1",
                                  uniform=bool(uniform)))
gate("G17", g17, "E2 kernel: EVERY row sum exactly 1 -- the raw two-rate writer with"
     " backtracking kept CONSERVES MEASURE at every dE")
gate("G18", g18, "E2 kernel doubly stochastic (column sums exactly 1): uniform is the"
     " Perron vector both ways; propagator diverges -- no gap opens")
gate("G19", g19, "E2 induced amplitudes: 5 fresh links equal, per-state redistribution"
     " total == 1 EXACTLY at every (u,b) (conserving redistribution)")
gate("G20", g20, "E2 UNIFORMITY iff dE == 0: m_back == m_fresh == 1/deg exactly at"
     " b = 1, and uniformity FAILS at every b != 1 (computed both directions)",
     "at b=1: m = %s; at b=1/2: (m_fresh, m_back) = (%s, %s)"
     % (MU_C_REF, F(1, 2) / (5 * F(1, 2) + 1), F(1) / (5 * F(1, 2) + 1)))
gate("G21", g21, "E2 drift: per-state drift == (v-u)*DIRS[d] exactly (path-relative"
     " persistence bias, nonzero at dE != 0) and stationary-average drift == 0"
     " exactly (NO spatial bias)")

# ---------------------------------------------------------------- E3a TRAIL WITH DECAY (H1 counting)
g22 = g23 = g24 = g25 = g26 = g27 = g28 = g29 = g30 = True
mu_by_b = {}
for u in U_SAMPLES:
    for b in B_SAMPLES:
        v = u * b
        W = kernel_pos(nbr, 1 - DEG * v - u, v)
        rs = set(row_sums(W))
        if rs != {1 - u} or (1 - u) == 1:
            g22 = False
        cs, ms = extract_pos(W, nbr)
        if len(cs) != 1 or len(ms) != 1:
            g23 = False
        c = next(iter(cs))
        m = next(iter(ms))
        if m != b / (DEG * b + 1) or m != v / (DEG * v + u):
            g23 = False
        if not (m < MU_C_REF):
            g24 = False
        if not op_identity_pos(W, nbr, c, m):
            g25 = False
        mu_by_b.setdefault(b, set()).add(m)
        results["E3a"].append(dict(u=str(u), b=str(b), row_sum=str(1 - u),
                                   loss=str(u), mu_induced=str(m),
                                   mass_ratio_muc_over_mu=str(MU_C_REF / m)))
for b, s in mu_by_b.items():
    if len(s) != 1:
        g26 = False   # E_b/f0 must drop out exactly
# measure evolution + convergence witness (representative sample)
u, b = U_SAMPLES[0], B_SAMPLES[3]
v = u * b
W = kernel_pos(nbr, 1 - DEG * v - u, v)
vec = [F(0)] * NC
vec[x0] = F(1)
tot_sum = F(0)
for k in range(K_MEAS + 1):
    tot = sum(vec)
    tot_sum += tot
    if tot != (1 - u) ** k:
        g27 = False
    if k < K_MEAS:
        vec = apply_forward(W, vec)
closed = (1 - (1 - u) ** (K_MEAS + 1)) / u
if tot_sum != closed or not (tot_sum < 1 / u):
    g28 = False
for u2 in U_SAMPLES:
    for b2 in B_SAMPLES:
        m = b2 / (DEG * b2 + 1)
        if MU_C_REF / m != (DEG * b2 + 1) / (DEG * b2):
            g29 = False
mus_sorted = [next(iter(mu_by_b[b2])) for b2 in sorted(B_SAMPLES)]
if not all(mus_sorted[i] < mus_sorted[i + 1] for i in range(len(mus_sorted) - 1)):
    g30 = False
gate("G22", g22, "E3a kernel: EVERY row sum exactly 1-u < 1 -- measure LOST at exactly"
     " the model's own erase probability u (the D-15 nonzero control beside G10/G17)")
gate("G23", g23, "E3a induced per-link amplitude == v/(deg*v+u) == b/(deg*b+1) EXACT:"
     " f0 and E_b cancel; the induced amplitude depends on dE/kT ALONE")
gate("G24", g24, "E3a mu_induced < MU_C_REF STRICTLY at every (u,b): the"
     " non-conserving ensemble lands computably OFF criticality (DONE_WHEN control)",
     "mu(b) = %s" % {str(b2): str(next(iter(mu_by_b[b2]))) for b2 in sorted(B_SAMPLES, reverse=True)})
gate("G25", g25, "E3a operator identity entrywise: I - W == c(I - m A), c = deg*v+u")
gate("G26", g26, "E3a barrier/attempt-clock independence: induced mu identical across"
     " u samples at fixed b (E_b and f0 set the clock, never the criticality)")
gate("G27", g27, "E3a exact measure evolution: total == (1-u)^k at every step k <= %d"
     % K_MEAS)
gate("G28", g28, "E3a propagator partial sums == (1-(1-u)^(K+1))/u exactly, < 1/u:"
     " CONVERGENT (massive) beside E1's linear divergence (massless)",
     "sum = %s < 1/u = %s" % (tot_sum, 1 / u))
gate("G29", g29, "E3a mass ratio mu_c/mu == (deg*b+1)/(deg*b) == 1 + e^{dE/kT}/deg"
     " EXACT as rationals at every sample (the O-58 N2 gap, computed)")
gate("G30", g30, "E3a induced mu strictly increasing in b: the gap closes only as the"
     " erase channel becomes negligible relative to writing")

# ---------------------------------------------------------------- E3b TRAIL WITH DECAY (NB counting)
g31 = g32 = g33 = g34 = True
B = nb_edge_adjacency(nbr)
rsB = set(row_sums(B))
if len(rsB) != 1:
    g31 = False
DEG_NB = next(iter(rsB))
MU_C_NB = F(1, DEG_NB)
for u in U_SAMPLES:
    for b in B_SAMPLES:
        v = u * b
        W = kernel_edge(nbr, 1 - 5 * v - u, v, None)   # retreat channel DROPPED = loss
        rs = set(row_sums(W))
        if rs != {1 - u}:
            g34 = False
        cs, mf, mb, tots = extract_edge(W, nbr)
        if len(cs) != 1 or len(mf) != 1 or mb != set():
            g32 = False
        m = next(iter(mf))
        if m != b / (5 * b + 1) or tots != {5 * b / (5 * b + 1)} or not (next(iter(tots)) < 1):
            g32 = False
        if MU_C_NB / m != 1 + 1 / (5 * b):
            g33 = False
        results["E3b"].append(dict(u=str(u), b=str(b), row_sum=str(1 - u),
                                   mu_induced=str(m),
                                   mass_ratio_mucNB_over_mu=str(MU_C_NB / m)))
gate("G31", g31, "E3b criticality reference EARNED: the venue's directed-edge"
     " (non-backtracking) operator has constant row sums (same Perron instrument)",
     "deg_NB = %s, MU_C_NB = %s" % (DEG_NB, MU_C_NB))
gate("G32", g32, "E3b induced per-link amplitude == b/(5b+1), redistribution total"
     " == 5b/(5b+1) < 1 strictly: off criticality in the NB counting too")
gate("G33", g33, "E3b mass ratio mu_c^NB/mu == 1 + e^{dE/kT}/5 exact: SAME LAW as E3a"
     " with l = 5 available links -- gap = ln(1 + e^{dE/kT}/l), counting-independent")
gate("G34", g34, "E2/E3 closure: E3's row-sum deficit == u == exactly E2's retreat"
     " amplitude -- the measure the extension-only counting loses IS the"
     " backtracking channel H1 keeps")

print()
print("=" * 78)
print("SECTION 4 -- COMPARISON (floats, LABELED, run AFTER all results above)")
print("=" * 78)
from math import log, exp
g35 = True
print("  %-8s %-12s %-14s %-16s %-16s" % ("b", "dE/kT", "mu_E3a", "gap ln(mu_c/mu)", "ln(1+e^x/6)"))
for b in B_SAMPLES:
    x = -log(b)
    m = b / (DEG * b + 1)
    gap_exact_ratio = MU_C_REF / m
    gap = log(gap_exact_ratio)
    law = log(1 + exp(x) / DEG)
    if abs(gap - law) > CMP_TOL:
        g35 = False
    print("  %-8s %-12.6f %-14s %-16.10f %-16.10f" % (b, x, m, gap, law))
gate("G35", g35, "COMPARISON: gap ln(mu_c/mu) agrees with the closed law"
     " ln(1 + e^{dE/kT}/deg) at every sample (float check of the exact rationals)")
m_half = F(1, 2) / (DEG * F(1, 2) + 1)
g36 = (m_half == SEALED_MU_SUB_ROW)
gate("G36", g36, "COMPARISON-TO-SEALED: at dE = kT ln 2 the trail ensemble's induced"
     " amplitude EQUALS the sealed T44-B subcritical sweep row mu = 1/8, whose class"
     " was already computed EXPONENTIAL with OZ-rate agreement",
     "mu_E3a(b=1/2) = %s == sealed row %s" % (m_half, SEALED_MU_SUB_ROW))
print("  NOTE (comparison): at large dE/kT the gap grows as dE/kT - ln(deg): the well")
print("  asymmetry prices the coupling's range linearly, offset by the venue's ln(deg).")
print("  At dE = 0 the minimal trail-ensemble gap is ln(7/6) (H1 counting) or ln(6/5)")
print("  (NB counting): the bare erase channel alone keeps written-trail mediation")
print("  massive.  Persistent-walk remark for E2 at b != 1 (owner: Goldstein 1951,")
print("  Kac 1974; comparison-only): path-persistence classically renormalizes the")
print("  diffusion constant, not the class; the E2 class computation is NOT run here")
print("  and is named open below.")

print()
print("=" * 78)
print("SECTION 5 -- VERDICTS (computed conjunctions of the gates above)")
print("=" * 78)
V1 = g10 and g11 and g12 and g13 and g14 and g15 and g16
V2 = g17 and g18 and g19 and g20 and g21
V3 = g22 and g23 and g24 and g25 and g26 and g27 and g28 and g29 and g30 and g32 and g33 and g34
V4 = g22 and g24 and g32          # DONE_WHEN: non-conserving ensemble computably OFF criticality
V5 = V1 and V2 and V3
print("V1 [%s] TRANSPORT (iv' literal): the energy-conserving writer conserves measure"
      % ("PASS" if V1 else "FAIL"))
print("      IDENTICALLY and induces mu = 1/deg = mu_c at EVERY dE, EVERY barrier,")
print("      both saddle readings: criticality is structural, not tuned.")
print("V2 [%s] TRAIL WITH RETREAT (raw two rates, H1 backtracking kept): measure"
      % ("PASS" if V2 else "FAIL"))
print("      conserved at every dE; the induced redistribution is the UNIFORM 1/deg")
print("      EXACTLY at dE = 0 and non-uniform (path-persistent, zero spatial drift)")
print("      at every dE != 0; no loss channel, propagator diverges.")
print("V3 [%s] TRAIL WITH DECAY (extension-only): measure conservation FAILS by exactly"
      % ("PASS" if V3 else "FAIL"))
print("      the model's erase probability; induced mu = 1/(deg + e^{dE/kT}) lands")
print("      strictly below mu_c at EVERY dE >= 0; mass gap = ln(1 + e^{dE/kT}/l),")
print("      counting-independent in form; f0 and E_b drop out exactly.")
print("V4 [%s] DONE_WHEN CONTROL: an ensemble that does NOT conserve measure lands"
      % ("PASS" if V4 else "FAIL"))
print("      computably OFF criticality -- and it is not an artificial control: it is")
print("      the model's own written-trail mediation with its own decay channel.")
print("V5 [%s] O-58 N2 ON THE WORLD VENUE: in every computed construction, the induced"
      % ("PASS" if V5 else "FAIL"))
print("      amplitude sits AT the measure-conserving normalization 1/deg exactly when")
print("      the writer step loses no measure AND the venue's link symmetry is unbroken")
print("      (both computed); the one measure leak the model owns (erasure) prices the")
print("      gap ln(1 + e^{dE/kT}/l).  Masslessness is measure conservation, not tuning:")
print("      no tuning closes the gap while a loss channel exists, none is needed once")
print("      it does not.")

npass = sum(GATES)
print()
print("GATES: %d PASS, %d FAIL, %d total" % (npass, len(GATES) - npass, len(GATES)))

out = dict(
    lane="LANE_T48_C_WORLD", probe="PROBE 3 -- world venue writer-ensemble computation (O-58 N2)",
    declared=dict(N=N, K_WRAP=K_WRAP, K_AXIS=K_AXIS, K_BOX=K_BOX, BOX_B=BOX_B,
                  K_MEAS=K_MEAS, K_SUM=K_SUM,
                  U_SAMPLES=[str(x) for x in U_SAMPLES],
                  B_SAMPLES=[str(x) for x in B_SAMPLES],
                  MU_C_REF=str(MU_C_REF), MU_C_NB=str(MU_C_NB)),
    gates=dict(total=len(GATES), passed=npass),
    verdicts=dict(V1_transport_critical_at_every_dE=V1,
                  V2_retreat_conserving_uniform_iff_dE0=V2,
                  V3_decay_never_critical_gap_law=V3,
                  V4_done_when_control=V4,
                  V5_N2_world_venue=V5),
    results=results)
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, "t48c_world.RESULT.json"), "w") as f:
    json.dump(out, f, indent=1, sort_keys=True)
print("RESULT written: t48c_world.RESULT.json")
