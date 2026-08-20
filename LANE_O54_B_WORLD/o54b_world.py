"""O-54-B  THE TWO-REGION RELATION -- WORLD TIER.  Interface(A,B) vs EARNED separation.

THE HYPOTHESIS UNDER TEST (stated so it can fail): "the falloff is not Gamma's to give --
pure-Gamma two-region relations are TOPOLOGICAL (separation-independent or contact-only),
and any genuine falloff enters through the MEDIATOR (a coupling term under C-77's openness
clause)."  If the no-mediator sub-experiment shows a separation-dependent interface, the
hypothesis FAILS.

VENUE (generator side, hidden from the analysis where stated): two L x L blocks, A and B,
of census barrier records -- NAND-like floating-gate cells (record identity: flash gate,
E_a ~ 120 kT, owner LANE_GR1_CENSUS numbers.txt; carrier parameters q = 100 e, h = 10 nm,
pitch = 40 nm, eps_r = 3.9, owner LANE_T34_NAND).  Block B is offset by a gap of g pitches;
nearest-cell distance g*pitch, hidden centroid separation (g + L - 1)*pitch.

THE EARNED SEPARATION AXIS (both sub-experiments): T-42-B part (ii)'s instrument,
RECOMPUTED per configuration, never assumed: the analysis receives ONLY the interaction
matrix U, applies the map delta = U^(-1/3) that T-42-B pinned from both sides, runs
classical MDS (Torgerson-Gower, BORROWED), re-earns metric axioms + dimension per gap, and
reads the separation D = |centroid_A - centroid_B| in units of the embedding's own median
nearest-neighbour distance (an EARNED unit; no generator unit enters the axis).  Hidden
positions are consulted ONLY in the final comparison section.

SUB-EXPERIMENT (1)  CONTACT MODEL (no mediator).  Access channels nearest-neighbour only
-- the world-tier interface instrument of C-79 (T42_C: IFACE = number of region-complement
adjacencies, here counted between A and B), on the integer grid, exact.
Expected: contact-only.  D-15: the zero line sits beside (i) adjacent blocks (g=1) with
NON-ZERO interface, and (ii) the same channel-counting instrument fed mediated channels
(sub-experiment 2), which DOES produce a separation-dependent curve -- the control that
would have shown a falloff had one existed.

SUB-EXPERIMENT (2)  MEDIATED MODEL.  The image-monopole coupling
    U(r) = k q^2 (1/r - 1/sqrt(r^2 + 4 h^2))
is added as an EXPLICITLY DECLARED coupling term (openness clause).  Its r^-3 large-r law
is INSERTED PHYSICS (classical image construction, via LANE_T34_NAND) -- nothing about it
is derived here, and nothing about it may be credited to Gamma.  What is measured:
    K(g)      = cross-block coupling matrix U on A x B          (the declared mediator)
    S(g)      = sum of K            (coupling-weighted interface STRENGTH)
    Kmax(g)   = largest entry       (the single strongest channel; the pair law)
    rank_tau  = # singular values > tau * sigma_1, tau in {1e-2, 1e-4, 1e-6}  (DECLARED
                relative thresholds; coupling-weighted interface RANK)
    N_chan    = # entries above the DECLARED absolute threshold U(12*pitch)
Laws are read by finite-difference log-log exponents against the EARNED axis -- no fit
anywhere -- with a two-point out-of-sample prediction check.  REPORTABLE FINDINGS ONLY:
(a) whether the BLOCK-level law inherits the pair exponent or aggregation changes it
(computed, not assumed), and (b) that the falloff enters through the coupling term, not
through Gamma.  CONTROL: shuffled positions (cells reassigned to positions at random,
seeded) must BREAK the law.

GAMMA-SIDE REMARK (cited, not recomputed): the world tier's pure-Gamma metric d_W
(T-42-B part (i), C-78) is a function of the FLIP SET alone -- additivity was EARNED over
all conversion paths -- so it contains no spatial-separation dependence at all; the only
pure-Gamma spatial object available to two regions is the access-channel count tested in
sub-experiment (1).

DISCIPLINE: D-1 no classical gravitational form in any construction step; the comparison
target may be named ONLY in the final comparison section.  D-15 controls beside every zero
and every flat line.  D-22 venue Aut computed.  D-24 audit table at the end.  No literal
verdicts: every PASS/FAIL is gated by a computed boolean.  Exact integers for all contact
counts and channel counts; couplings involve sqrt (irrational), so float64 with stated
tolerances -- the deviation from the exact-arithmetic rule is stated openly, as in T-42-B.

OWNERS: census E_a -- LANE_GR1_CENSUS.  Image-charge electrostatics -- classical (Thomson
image construction) via LANE_T34_NAND.  Classical MDS -- Torgerson-Gower (BORROWED).
delta = U^(-1/3) map + earned-dimension instrument -- LANE_T42_B (EARNED there, reused
with citation).  IFACE-as-adjacency-channels -- LANE_T42_C / C-79.  Topological mutual
information (the corner-tier analogue of a separation-independent two-region term) --
Kitaev-Preskill 2006, Levin-Wen 2006; named for orientation only, not used.
"""
import sys
import numpy as np

def say(*a):
    print(*a); sys.stdout.flush()

def verdict(ok):            # gate: computed boolean -> word. Never called with a literal.
    return "PASS" if bool(ok) else "FAIL"

say("=" * 100)
say("O-54-B   TWO-REGION INTERFACE vs EARNED SEPARATION -- WORLD TIER")
say("=" * 100)

# ------------------------------------------------------------------ carrier constants (T34/GR1)
E_CH  = 1.602176634e-19
EPS0  = 8.8541878128e-12
q, h, pitch, eps_r = 100 * E_CH, 10e-9, 40e-9, 3.9          # LANE_T34_NAND
kq2 = q * q / (4 * np.pi * eps_r * EPS0)

def U_of_r(r_m):
    """Corrected image-monopole law (INSERTED coupling term, declared)."""
    return kq2 * (1.0 / r_m - 1.0 / np.sqrt(r_m * r_m + 4 * h * h))

L = 8                          # block side
NA = NB = L * L                # 64 records per block
GAPS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 28, 32]

def block_positions(g):
    """Hidden integer-grid positions (pitch units). A at x in [0,L), B at x in [L-1+g, ...)."""
    A = [(x, y) for y in range(L) for x in range(L)]
    B = [(x + L - 1 + g, y) for y in range(L) for x in range(L)]
    return A, B

def build_U(pos_units):
    P = np.asarray(pos_units, dtype=float) * pitch
    m = len(P)
    Umat = np.zeros((m, m))
    for i in range(m):
        for j in range(i + 1, m):
            r = float(np.hypot(*(P[i] - P[j])))
            Umat[i, j] = Umat[j, i] = U_of_r(r)
    return Umat

# ------------------------------------------------------------------ earning instrument (T42_B)
def triangle_violations(delta, rtol=1e-9):
    m = len(delta); viol = 0
    mx = np.max(delta)
    for k in range(m):
        excess = delta - (delta[:, [k]] + delta[[k], :])
        viol += int(np.count_nonzero(excess > rtol * mx))
    return viol

def mds(delta):
    m = len(delta)
    J = np.eye(m) - np.ones((m, m)) / m
    B = -0.5 * J @ (delta ** 2) @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]
    return w[order], V[:, order]

def earn(Umat):
    """The analysis side: receives ONLY U. Returns earned quantities."""
    delta = Umat.copy(); np.fill_diagonal(delta, 1.0)
    delta = delta ** (-1.0 / 3.0); np.fill_diagonal(delta, 0.0)
    delta /= np.max(delta)
    tv = triangle_violations(delta)
    ev, V = mds(delta)
    lam1 = ev[0]
    dim = int(np.count_nonzero(ev > 0.01 * lam1))
    negfrac = abs(min(ev.min(), 0.0)) / lam1
    X = V[:, :2] * np.sqrt(np.maximum(ev[:2], 0.0))
    # earned unit: median nearest-neighbour distance in the embedding itself
    m = len(X)
    d2 = np.sum((X[:, None, :] - X[None, :, :]) ** 2, axis=2)
    np.fill_diagonal(d2, np.inf)
    unit = float(np.median(np.sqrt(np.min(d2, axis=1))))
    return dict(tv=tv, ev=ev, dim=dim, negfrac=negfrac, X=X, unit=unit,
                l2_over_l1=ev[1] / ev[0])

def earned_separation(X, unit, idxA, idxB):
    cA = X[idxA].mean(axis=0); cB = X[idxB].mean(axis=0)
    return float(np.linalg.norm(cA - cB)) / unit

def earned_nn_separation(X, unit, idxA, idxB):
    """Earned NEAREST-cell separation: min cross-block distance in the embedding, earned units."""
    d2 = np.sum((X[idxA][:, None, :] - X[idxB][None, :, :]) ** 2, axis=2)
    return float(np.sqrt(np.min(d2))) / unit

# =====================================================================================
# D-22  VENUE GEOMETRY
# =====================================================================================
say("")
say("D-22  VENUE: two L x L blocks on the integer grid, gap g along x.")
gcheck = 8
A0, B0 = block_positions(gcheck)
posU = A0 + B0
Uv = build_U(posU)
span_x = 2 * L - 1 + gcheck   # x runs 0 .. span_x-1 (hidden; generator-side check only)
def apply_perm(tr):
    lookup = {p: i for i, p in enumerate(posU)}
    return [lookup[tr(p)] for p in posU]
trs = [("identity",   lambda p: p),
       ("x-mirror",   lambda p: (span_x - 1 - p[0], p[1])),
       ("y-mirror",   lambda p: (p[0], L - 1 - p[1])),
       ("180-rot",    lambda p: (span_x - 1 - p[0], L - 1 - p[1]))]
n_aut = 0
for nm, tr in trs:
    P = np.array(apply_perm(tr))
    if np.allclose(Uv[np.ix_(P, P)], Uv, rtol=1e-12, atol=0):
        n_aut += 1
aut_ok = (n_aut == 4)
say(f"  permutations preserving U at g={gcheck}: {n_aut}/4 (identity, x-mirror, y-mirror, 180-rot)")
say(f"  => venue Aut contains Z_2 x Z_2; x-mirror SWAPS the blocks: A and B are exchangeable,")
say(f"     so Interface(A,B) is symmetric by venue symmetry .......... {verdict(aut_ok)}")

# =====================================================================================
# STEP 1  EARN THE SEPARATION AXIS (recomputed per gap; positions never given)
# =====================================================================================
say("")
say("=" * 100)
say("STEP 1  THE EARNED SEPARATION AXIS (T42_B instrument, RECOMPUTED per configuration)")
say("-" * 100)
say(f"  per gap: analysis receives only the {2*NA}x{2*NA} matrix U; map delta = U^(-1/3)")
say(f"  (EARNED by T42_B's two-sided p-scan, reused with citation); triangle scan rtol 1e-9;")
say(f"  MDS dim@1%; earned unit = median NN distance in the embedding (float64, sqrt venue).")
say("")
say(f"  TWO earned separations, both read off the embedding: D_cent (between block centroids)")
say(f"  and D_nn (nearest cross-block cells) -- each aggregate below is measured against the")
say(f"  axis that is its own scale (the strongest single channel lives at the nearest cells).")
say("")
say(f"    {'g':>4}{'tri-viol':>9}{'dim@1%':>8}{'|neg|/l1':>10}{'l2/l1':>9}{'D_cent':>9}{'D_nn':>9}   gate")

idxA = list(range(NA)); idxB = list(range(NA, 2 * NA))
D_earn, D_nn, earn_ok_by_g, U_by_g, l2l1 = {}, {}, {}, {}, {}
for g in GAPS:
    A, Bl = block_positions(g)
    Umat = build_U(A + Bl)
    U_by_g[g] = Umat
    e = earn(Umat)
    D_earn[g] = earned_separation(e["X"], e["unit"], idxA, idxB)
    D_nn[g]   = earned_nn_separation(e["X"], e["unit"], idxA, idxB)
    l2l1[g]   = float(e["l2_over_l1"])
    ok = (e["tv"] == 0) and (e["negfrac"] < 0.01) and (e["dim"] == 2)
    earn_ok_by_g[g] = ok
    say(f"    {g:>4}{e['tv']:>9d}{e['dim']:>8d}{e['negfrac']:>10.4f}{e['l2_over_l1']:>9.4f}"
        f"{D_earn[g]:>9.3f}{D_nn[g]:>9.3f}   {verdict(ok)}")

earn_all = all(earn_ok_by_g[g] for g in GAPS)
mono = all(D_earn[a] < D_earn[b] and D_nn[a] < D_nn[b] for a, b in zip(GAPS, GAPS[1:]))
say("")
say(f"  INSTRUMENT SCOPE (computed, stated openly): dim@1% reads 2 only while l2/l1 > 0.01;")
say(f"  at g={GAPS[-1]} the two-block venue reaches l2/l1 = {l2l1[GAPS[-1]]:.4f} -- the sweep ends where")
say(f"  the earned-dimension instrument still certifies the venue, not where the law runs out.")
say("")
say(f"  metric + dim 2 re-earned at every gap ........................ {verdict(earn_all)}")
say(f"  D_earned strictly increasing across the sweep (earned data alone) .... {verdict(mono)}")
say(f"  axis dynamic range: D_earned {D_earn[GAPS[0]]:.2f} -> {D_earn[GAPS[-1]]:.2f}"
    f" (x{D_earn[GAPS[-1]]/D_earn[GAPS[0]]:.1f})")

# D-15 for the earning instrument itself: shuffled-U must fail to embed
rng = np.random.default_rng(7)
Ush = U_by_g[8].copy()
iu = np.triu_indices(2 * NA, 1)
vals = Ush[iu].copy(); rng.shuffle(vals)
Ush = np.zeros_like(Ush); Ush[iu] = vals; Ush += Ush.T
esh = earn(Ush)
instr_ctrl = (esh["negfrac"] > 0.05) or (esh["tv"] > 0)
say(f"  D-15 (instrument): SHUFFLED-ENTRY U at g=8: tri-viol {esh['tv']}, |neg|/l1"
    f" {esh['negfrac']:.3f} -> must fail to embed .... {verdict(instr_ctrl)} (fails, as required)")

# =====================================================================================
# STEP 2  SUB-EXPERIMENT (1): CONTACT MODEL (no mediator)
# =====================================================================================
say("")
say("=" * 100)
say("SUB-EXPERIMENT (1)  CONTACT MODEL -- access channels nearest-neighbour only (no mediator)")
say("-" * 100)
say("  Interface(A,B) = number of A-B nearest-neighbour adjacencies on the integer grid")
say("  (C-79's world-tier IFACE, owner LANE_T42_C, counted between regions; EXACT integers).")
say("")
say(f"    {'g':>4}{'D_earned':>10}{'IFACE_contact':>15}")
iface_contact = {}
for g in GAPS:
    A, Bl = block_positions(g)
    Aset, Bset = set(A), set(Bl)
    cnt = 0
    for (x, y) in A:
        for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            if (x + dx, y + dy) in Bset:
                cnt += 1
    iface_contact[g] = cnt
    say(f"    {g:>4}{D_earn[g]:>10.3f}{cnt:>15d}")

zero_sep   = all(iface_contact[g] == 0 for g in GAPS if g >= 2)
adj_nonzero = iface_contact[1] > 0
axis_varies = (D_earn[GAPS[-1]] / D_earn[1]) > 4.0
say("")
say(f"  RESULT: Interface = {iface_contact[1]} at contact (g=1), exactly 0 for every g >= 2,")
say(f"  while the earned separation varies by x{D_earn[GAPS[-1]]/D_earn[1]:.1f}: the pure-access")
say(f"  two-region relation is CONTACT-ONLY -- a topological (adjacency) quantity, carrying")
say(f"  NO dependence on the earned separation.")
say(f"    zero at every separated gap (exact count) ................. {verdict(zero_sep)}")
say(f"    D-15 control (i): adjacent blocks non-zero ({iface_contact[1]} channels = L) "
    f"....... {verdict(adj_nonzero and iface_contact[1] == L)}")
say(f"    earned axis genuinely varied under the zero ............... {verdict(axis_varies)}")
# contact-radius scan: the zero tracks the DECLARED contact radius, not the separation
say("")
say("  CONTACT-RADIUS SCAN (control: the zero begins exactly where the declared radius ends,")
say("  i.e. the zero is the access model's, not the counting instrument's):")
radius_ok = True
for R_c in (1, 2, 3):
    row = []
    for g in GAPS[:6]:
        A, Bl = block_positions(g)
        PA = np.array(A, float); PB = np.array(Bl, float)
        d2 = np.sum((PA[:, None, :] - PB[None, :, :]) ** 2, axis=2)
        cnt = int(np.count_nonzero(d2 <= R_c * R_c + 1e-12))
        row.append((g, cnt))
        radius_ok = radius_ok and ((cnt > 0) == (g <= R_c))
    say(f"    radius {R_c} pitch: " + "  ".join(f"g={g}:{c}" for g, c in row))
say(f"    non-zero exactly for g <= declared radius, zero beyond .... {verdict(radius_ok)}")
say(f"  (the falloff-capable control for this flat zero line is sub-experiment (2): the SAME")
say(f"   venue and the same earned axis, where channels exist at every separation and the")
say(f"   instrument reads a non-flat law -- see below)")

# =====================================================================================
# STEP 3  SUB-EXPERIMENT (2): MEDIATED MODEL (declared coupling term)
# =====================================================================================
say("")
say("=" * 100)
say("SUB-EXPERIMENT (2)  MEDIATED MODEL -- the image-monopole coupling as a DECLARED term")
say("-" * 100)
say("  INSERTED: U(r) = kq^2(1/r - 1/sqrt(r^2+4h^2)) between every A-cell and B-cell.")
say("  Its large-r r^-3 form is the MEDIATOR'S OWN LAW (classical image construction, owner")
say("  LANE_T34_NAND) -- inserted physics, reportable only as the coupling's contribution.")
say("  Aggregates of the 64x64 cross-block matrix K (float64; sums exact to machine precision,")
say("  SVD tolerance stated by threshold):")
say("")
hdr = (f"    {'g':>4}{'D_earned':>10}{'S=sum K (J)':>14}{'Kmax (J)':>12}"
       f"{'rk@1e-2':>9}{'rk@1e-4':>9}{'rk@1e-6':>9}{'N_chan':>8}")
say(hdr)
TAUS = (1e-2, 1e-4, 1e-6)
theta = U_of_r(12 * pitch)     # DECLARED absolute channel threshold
S, Kmax, ranks, Nchan = {}, {}, {}, {}
for g in GAPS:
    K = U_by_g[g][np.ix_(idxA, idxB)]
    S[g] = float(np.sum(K))
    Kmax[g] = float(np.max(K))
    sv = np.linalg.svd(K, compute_uv=False)
    ranks[g] = tuple(int(np.count_nonzero(sv > t * sv[0])) for t in TAUS)
    Nchan[g] = int(np.count_nonzero(K > theta))
    say(f"    {g:>4}{D_earn[g]:>10.3f}{S[g]:>14.4e}{Kmax[g]:>12.4e}"
        f"{ranks[g][0]:>9d}{ranks[g][1]:>9d}{ranks[g][2]:>9d}{Nchan[g]:>8d}")

say("")
say(f"  N_chan threshold declared as U(12*pitch) = {theta:.3e} J: the count reaches 0 for")
say(f"  g > 12 -- ANY absolute-thresholded channel notion is contact-like with a radius set")
say(f"  by the threshold (an induced observation about thresholding, not a falloff law).")

# ---- finite-difference exponents against the EARNED axes ----------------------------
say("")
say("  FINITE-DIFFERENCE LOG-LOG EXPONENTS (no fit anywhere). Each aggregate against the")
say("  earned axis that is its own scale: S (block aggregate) vs D_cent; Kmax (the single")
say("  strongest channel, which lives at the nearest cells) vs D_nn. The p_pair column is")
say("  the DECLARED law's own analytic slope at the matching scale (D_cent x declared pitch,")
say("  T34 GIVEN-DATA) -- labelling of the inserted physics, not a finding.")
say(f"    {'interval':>14}{'p_eff(S)':>10}{'p_eff(Kmax)':>12}{'p_pair(midpt)':>14}")
def fd_exp(y, axis, g1, g2):
    return -(np.log(y[g2]) - np.log(y[g1])) / (np.log(axis[g2]) - np.log(axis[g1]))
pS, pK, ppair_mid = {}, {}, {}
for g1, g2 in zip(GAPS, GAPS[1:]):
    pS[(g1, g2)] = fd_exp(S, D_earn, g1, g2)
    pK[(g1, g2)] = fd_exp(Kmax, D_nn, g1, g2)
    Dm = 0.5 * (D_earn[g1] + D_earn[g2]) * pitch
    eps = 1e-6
    ppair_mid[(g1, g2)] = -(np.log(U_of_r(Dm * (1 + eps))) - np.log(U_of_r(Dm))) / np.log(1 + eps)
    say(f"    {g1:>6}->{g2:<6}{pS[(g1,g2)]:>10.3f}{pK[(g1,g2)]:>12.3f}{ppair_mid[(g1,g2)]:>14.3f}")

intervals = list(zip(GAPS, GAPS[1:]))
far_pairs = intervals[-3:]
p_far_S = float(np.mean([pS[p] for p in far_pairs]))
p_far_K = float(np.mean([pK[p] for p in far_pairs]))
p_near_S = pS[intervals[0]]
falloff_seen = S[GAPS[-1]] < S[GAPS[0]] / 50.0

# the aggregation correction, MEASURED: residual p_eff(S) - 3 and its own decay exponent
resid = {iv: pS[iv] - 3.0 for iv in intervals}
resid_pos_mono = all(resid[a] > resid[b] > 0 for a, b in zip(intervals, intervals[1:]))
def midD(iv):
    return 0.5 * (D_earn[iv[0]] + D_earn[iv[1]])
q_resid = (np.log(resid[intervals[-4]]) - np.log(resid[intervals[-1]])) \
          / (np.log(midD(intervals[-1])) - np.log(midD(intervals[-4])))
q_ok = 1.5 < q_resid < 2.5
# the strongest single channel: must track the inserted pair law on its own earned axis
pK_vs_pair = max(abs(pK[iv] - ppair_mid[iv]) for iv in far_pairs)
# ppair at the Kmax scale for the far intervals (D_nn axis):
pK_pair_far = []
for g1, g2 in far_pairs:
    Dm = 0.5 * (D_nn[g1] + D_nn[g2]) * pitch
    eps = 1e-6
    pK_pair_far.append(-(np.log(U_of_r(Dm * (1 + eps))) - np.log(U_of_r(Dm))) / np.log(1 + eps))
pK_match = max(abs(pK[iv] - pp) for iv, pp in zip(far_pairs, pK_pair_far)) < 0.02
inherit = (p_far_S - 3.0 < 0.1) and resid_pos_mono and q_ok and pK_match

say("")
say(f"  READINGS (gated):")
say(f"    the mediated interface is NON-FLAT: S drops x{S[GAPS[0]]/S[GAPS[-1]]:.0f} over the sweep")
say(f"    (this is the falloff-capable control beside sub-experiment (1)'s zero) "
    f".... {verdict(falloff_seen)}")
say(f"    FINDING (a) -- WHAT AGGREGATION ACTUALLY GIVES (computed, not assumed):")
say(f"      the strongest single channel Kmax tracks the inserted pair law on its own earned")
say(f"      axis D_nn: far intervals match the declared law's analytic slope to <0.02"
    f" .... {verdict(pK_match)}")
say(f"      the BLOCK aggregate S runs STEEPER than the pair law at finite separation")
say(f"      (contact interval p_eff = {p_near_S:.3f}) and approaches the pair exponent 3 from")
say(f"      above: residual p_eff(S)-3 positive and strictly decreasing at every interval"
    f" .... {verdict(resid_pos_mono)}")
say(f"      the residual's own decay exponent (finite differences on the residuals, far range):")
say(f"      q = {q_resid:.2f} -- a ~1/D^2 aggregation correction (multipole-type, from block")
say(f"      extent), vanishing as separation grows .... {verdict(q_ok)}")
say(f"      far-range block exponent p_eff(S) = {p_far_S:.3f}: the block-level law INHERITS the")
say(f"      inserted pair exponent asymptotically; aggregation changes the exponent only by")
say(f"      the measured, decaying finite-size correction ............ {verdict(inherit)}")
say(f"    coupling-weighted interface RANK falls with separation at every declared threshold")
rank_mono = all(all(ranks[a][k] >= ranks[b][k] for k in range(3)) for a, b in zip(GAPS, GAPS[1:]))
rank_collapse = (ranks[GAPS[-1]][0] * 2 <= ranks[GAPS[0]][0]) and (ranks[GAPS[-1]][0] <= 2)
say(f"    (monotone non-increasing, all thresholds) ................. {verdict(rank_mono)}")
say(f"    near blocks couple through {ranks[GAPS[0]][0]} channels @1e-2; far blocks through"
    f" {ranks[GAPS[-1]][0]}: the")
say(f"    mediated interface collapses toward a single aggregate channel"
    f" .... {verdict(rank_collapse)}")
say(f"    the rank sequence is threshold-dependent (declared taus differ column to column):")
say(f"    reported as the coupling's structure, no universal law claimed")

# ---- out-of-sample checks (no fit: two points determine exponent+amplitude) ---------
say("")
say("  OUT-OF-SAMPLE (two-point, zero fitted parameters):")
g1, g2 = 20, 24
p_oos = fd_exp(S, D_earn, g1, g2)
amp = S[g1] * D_earn[g1] ** p_oos
oos_ok = True
for gt in (28, 32):
    pred = amp * D_earn[gt] ** (-p_oos)
    rel = abs(pred - S[gt]) / S[gt]
    oos_ok = oos_ok and (rel < 0.05)
    say(f"    exponent+amplitude from g=({g1},{g2}) -> predict S at g={gt}:"
        f" {pred:.4e} vs {S[gt]:.4e}  rel err {rel:.4f} .... {verdict(rel < 0.05)}")

# =====================================================================================
# STEP 4  CONTROL: SHUFFLED POSITIONS MUST BREAK THE LAW
# =====================================================================================
say("")
say("=" * 100)
say("STEP 4  D-15 CONTROL -- SHUFFLED POSITIONS (seeded): cells reassigned to positions at")
say("random; the same 128 positions, the same coupling law, block LABELS kept. The earned")
say("separation must collapse and the law must break.")
say("-" * 100)
say(f"    {'g':>4}{'D_earned_sh':>12}{'S_sh (J)':>14}")
rngS = np.random.default_rng(42)
D_sh, S_sh = {}, {}
for g in GAPS:
    perm = rngS.permutation(2 * NA)
    Ush2 = U_by_g[g][np.ix_(perm, perm)]        # cells relabelled = positions shuffled
    e2 = earn(Ush2)
    D_sh[g] = earned_separation(e2["X"], e2["unit"], idxA, idxB)
    K2 = Ush2[np.ix_(idxA, idxB)]
    S_sh[g] = float(np.sum(K2))
    say(f"    {g:>4}{D_sh[g]:>12.3f}{S_sh[g]:>14.4e}")
collapse = max(D_sh.values()) < 0.25 * min(D_earn[g] for g in GAPS if g >= 8)
sh_mono_broken = not all(D_sh[a] < D_sh[b] for a, b in zip(GAPS, GAPS[1:]))
# exponent of S_sh against the TRUE sweep axis D_earn: must not read 3
p_sh = -(np.log(S_sh[32]) - np.log(S_sh[20])) / (np.log(D_earn[32]) - np.log(D_earn[20]))
sh_law_broken = abs(p_sh - 3.0) > 0.5
say("")
say(f"  earned inter-block separation collapses (max {max(D_sh.values()):.2f} vs real"
    f" axis >= {min(D_earn[g] for g in GAPS if g>=8):.2f}) .... {verdict(collapse)}")
say(f"  monotone ordering destroyed ................................. {verdict(sh_mono_broken)}")
say(f"  far-range exponent of S_sh vs sweep axis: {p_sh:.3f} (law's 3 gone)"
    f" .... {verdict(sh_law_broken)}")
shuffle_breaks = collapse and sh_mono_broken and sh_law_broken

# =====================================================================================
# FINAL COMPARISON SECTION (hidden positions and the named target appear ONLY here)
# =====================================================================================
say("")
say("=" * 100)
say("FINAL COMPARISON (hidden generator values and the comparison target, here and only here)")
say("-" * 100)
dev_max = max(abs(D_earn[g] - (g + L - 1)) / (g + L - 1) for g in GAPS)
say(f"  (a) earned vs hidden separation: D_earned agrees with the hidden centroid distance")
say(f"      (g + {L-1} pitches) to max relative deviation {dev_max:.4f} across the sweep"
    f" .... {verdict(dev_max < 0.02)}")
say(f"      (the earning recomputed the axis; nothing was assumed from the generator)")
Sasym_ok = True
say(f"  (b) inheritance against the inserted law's own asymptote S ~ N_A N_B * 2kq^2h^2 / D^3:")
for gt in (24, 32):
    Dm = (gt + L - 1) * pitch
    Sa = NA * NB * 2 * kq2 * h * h / Dm ** 3
    ratio = S[gt] / Sa
    Sasym_ok = Sasym_ok and abs(ratio - 1.0) < 0.1
    say(f"      g={gt}: S/S_asym = {ratio:.4f} .... {verdict(abs(ratio-1)<0.1)}")
say(f"  (c) THE COMPARISON TARGET (named here only): the Newton-form falloff, force ~ r^-2 /")
say(f"      potential ~ r^-1 between separated bodies. What this lane finds is NOT that form:")
say(f"      the mediated block law is the mediator's own r^-3 (image-dipole), INSERTED, and")
say(f"      inherited unchanged by aggregation; the pure-access (Gamma-side) relation is")
say(f"      contact-only, i.e. NO falloff at all. The division of labour reads: Gamma gives")
say(f"      topological/contact interface structure; ANY genuine falloff, including whatever")
say(f"      law the world's mediators actually carry, enters through the coupling term under")
say(f"      the openness clause. Nothing here derives, and nothing here credits to Gamma,")
say(f"      any specific falloff exponent.")

# =====================================================================================
# HYPOTHESIS VERDICT + D-24 AUDIT
# =====================================================================================
say("")
say("=" * 100)
say("HYPOTHESIS VERDICT")
say("-" * 100)
contact_topological = zero_sep and adj_nonzero and axis_varies and radius_ok
mediated_falloff = falloff_seen and inherit and oos_ok
hypothesis_holds = (contact_topological and mediated_falloff and earn_all and mono
                    and instr_ctrl and shuffle_breaks and aut_ok)
say(f"  pure-Gamma (access-only) two-region relation topological/contact-only"
    f" .... {verdict(contact_topological)}")
say(f"  genuine falloff present, and it enters through the DECLARED coupling term"
    f" .... {verdict(mediated_falloff)}")
say(f"  earned axis sound + instrument controls + shuffle control"
    f" .... {verdict(earn_all and mono and instr_ctrl and shuffle_breaks)}")
say(f"  => THE HYPOTHESIS HOLDS on this venue: the falloff is not Gamma's to give;")
say(f"     it is the mediator's, reported as the coupling's contribution"
    f" .... {verdict(hypothesis_holds)}")

say("")
say("=" * 100)
say("D-24 AUDIT TABLE (every concept used)")
say("-" * 100)
audit = [
 ("record identity",      "GIVEN-DATA", "census flash-gate barrier records, owner LANE_GR1_CENSUS"),
 ("carrier parameters",   "GIVEN-DATA", "q=100e, h=10nm, pitch=40nm, eps_r=3.9, owner LANE_T34_NAND"),
 ("coupling U(r)",        "INSERTED",   "image-monopole law, classical image construction via T34; its r^-3 is the mediator's, never Gamma's"),
 ("map delta=U^(-1/3)",   "EARNED (T42_B), cited", "pinned there from both sides; re-validated here per gap (tv=0, dim 2)" if earn_all else "REVALIDATION FAILED"),
 ("classical MDS",        "BORROWED",   "Torgerson-Gower double-centering (the earning instrument)"),
 ("earned separation D",  "EARNED" if earn_all and mono and dev_max < 0.02 else "NOT EARNED",
                          "recomputed per gap from U alone; unit = embedding's own median NN distance"),
 ("interface (contact)",  "EARNED",     "C-79 world-tier IFACE (adjacency channels, owner T42_C) counted between regions; exact integers"),
 ("interface (mediated)", "DECLARED",   "strength/rank/thresholded channels of the cross-block coupling matrix; thresholds declared"),
 ("exponents",            "MEASURED",   "finite-difference log-log vs earned axis; no fit; out-of-sample checked" if oos_ok else "OOS FAILED"),
 ("block-law inheritance","INDUCED" if inherit else "NOT SHOWN",
                          "aggregation over 64x64 grains keeps the far exponent at the inserted pair value (computed)"),
 ("contact-only zero",    "INDUCED" if contact_topological else "NOT SHOWN",
                          "consequence of access structure alone; controls: adjacent nonzero, radius scan, mediated non-flat"),
 ("topological MI (corner analogue)", "NOT USED", "named for orientation; owners Kitaev-Preskill 2006, Levin-Wen 2006"),
 ("Newton-form target",   "COMPARISON ONLY", "named in the final section; no gravitational form in any construction step (D-1)"),
 ("overall energy unit",  "IMPORTED",   "kq2/pitch set scales; every law read on log-log differences, scale-free"),
]
for c, s2, note in audit:
    say(f"  {c:<26}{s2:<26}{note}")

say("")
say(f"EXACT OBJECTS: all contact/adjacency counts, radius-scan counts, N_chan counts (integers on")
say(f"the integer grid). FLOAT OBJECTS: 128x128 U per gap ({len(GAPS)} gaps + shuffles), 64x64 SVDs,")
say(f"MDS eigenproblems; sqrt-irrational venue, tolerances: triangle rtol 1e-9, D4-check rtol 1e-12,")
say(f"thresholds as declared. Deviation from exact arithmetic stated openly (as in T42_B).")
overall = hypothesis_holds and Sasym_ok and rank_mono
say("")
say(f"OVERALL .... {verdict(overall)}")
say("")
say("NEXT STEP (named, per program rule): THREE-region composition -- does the mediated block")
say("law COMPOSE (A-B in the presence of C: superposition vs screening under the image")
say("construction), and does the coupling-weighted interface at contact (g=1) connect to the")
say("C-79 world form IFACE^3 = 216 CAP^2 -- i.e. does the mediated rank at contact reproduce")
say("the access-channel count, tying the coupling's near end to Gamma's topological end?")
