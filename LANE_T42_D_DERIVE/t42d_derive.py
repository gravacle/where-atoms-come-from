#!/usr/bin/env python3
# LANE_T42_D_DERIVE -- C-77 first falsifiable increment, executed honestly.
#
# THE ONE DERIVATION ATTEMPTED (construction uses NO classical gravitational form):
#   From the earned structures only --
#     interface rank IR2 (T42_C, adversarially verified; recomputed here from scratch by F_2 rank),
#     content capacity cap_pauli = 2*AQ - r_in (T42_C, verified),
#     lattice counts AQ/PER/AV/PD (T42_C, earned as counts),
#     polynomial degree by exact finite differences (T42_C world tier, earned instrument),
#     world-tier census access counts (GR1 census model; coupling term under the openness clause)
#   -- derive a scale-level relation between the Gamma boundary-aggregate and the content aggregate.
#
# DERIVED RELATION R1 (stated before checking; checked exactly below):
#   For thick nested regions the interface aggregate IR and the content aggregate C obey
#        deg(IR) = deg(C) - 1
#   with D := deg(C) the EARNED dimension (computed, never imported), equivalently the
#   parameter-free forms
#        corner (exact, all L, all s):  32*C = (IR+10)^2 + 8*(IR+10) - 160
#        world  (exact, all n):         IR^3 = 216 * C^2
#   i.e. IR ~ kappa^(1/D) * C^((D-1)/D) to leading order: the boundary aggregate is one degree
#   below content, and the relation between the two earned aggregates needs no external scale.
#   MECHANISM (theory step, verified): the per-interior-cell density of the Euler count
#   AV + PD - AQ is (1 vertex) + (1 face) - (2 edges) = 0, so IR2's LAW-3 source term receives
#   contributions ONLY from the boundary layer; in the world tier the nested difference
#   n^3 - (n-2)^3 cancels the bulk identically. Bulk cancellation IS the law.
#
# DISCIPLINE: every PASS/FAIL printed from a computed boolean (no literal verdicts);
# exact arithmetic only (F_2 bitmasks, python ints, Fraction); D-15 controls that FAIL;
# D-22 venue lines; classical relations named ONLY in the final comparison section.

import sys
from fractions import Fraction

GATES = []
def gate(name, ok):
    ok = bool(ok)
    GATES.append((name, ok))
    print(("PASS " if ok else "FAIL ") + name)
    return ok

# ---------------- exact F_2 rank (bitmask rows) ----------------
def rank_f2(rows):
    pivots = {}
    r = 0
    for row in rows:
        cur = row
        while cur:
            hb = cur.bit_length() - 1
            if hb in pivots:
                cur ^= pivots[hb]
            else:
                pivots[hb] = cur
                r += 1
                break
    return r

# ---------------- exact finite-difference degree (earned instrument, T42_C) ----------------
def fd_degree(seq):
    d = list(seq); k = 0
    while len(d) >= 2:
        if all(x == d[0] for x in d):
            return k
        d = [d[i+1] - d[i] for i in range(len(d)-1)]
        k += 1
    return None  # unresolved on this window

# ---------------- toric carrier (Kitaev, BORROWED construction; conventions per C-74/T42_C) ----
def toric(L):
    n = 2 * L * L
    def h(i, j): return (i % L) * L + (j % L)            # edge (i,j)-(i,j+1)
    def v(i, j): return L * L + (i % L) * L + (j % L)    # edge (i,j)-(i+1,j)
    stars, plaqs = [], []
    for i in range(L):
        for j in range(L):
            s = (1 << h(i, j)) | (1 << h(i, j - 1)) | (1 << v(i, j)) | (1 << v(i - 1, j))
            p = (1 << h(i, j)) | (1 << h(i + 1, j)) | (1 << v(i, j)) | (1 << v(i, j + 1))
            stars.append(s); plaqs.append(p)
    return n, h, v, stars, plaqs

def region_block(L, h, v, s, oi=0, oj=0):
    # C-74 convention: region = induced edge set of an s x s vertex block (no wrap; s <= L-2 here)
    verts = set((oi + a, oj + b) for a in range(s) for b in range(s))
    mask = 0
    for (i, j) in verts:
        if (i, j + 1) in verts: mask |= (1 << h(i, j))
        if (i + 1, j) in verts: mask |= (1 << v(i, j))
    return verts, mask

def edge_endpoints(L, e):
    if e < L * L:
        i, j = divmod(e, L); return (i, j), (i, (j + 1) % L)
    e -= L * L
    i, j = divmod(e, L); return (i, j), ((i + 1) % L, j)

def popcount(x): return bin(x).count("1")

def region_quantities(L, n, h, v, stars, plaqs, verts, maskR):
    maskAll = (1 << n) - 1
    maskOut = maskAll ^ maskR
    AQ = popcount(maskR)
    PER = 0
    for e in range(n):
        a, b = edge_endpoints(L, e)
        a = (a[0] % L, a[1] % L); b = (b[0] % L, b[1] % L)
        wa = a in vset_mod(verts, L); wb = b in vset_mod(verts, L)
        if wa != wb: PER += 1
    AV = sum(1 for m in stars if m & maskR)
    PD = sum(1 for m in plaqs if m & maskR)
    rM = rank_f2(stars); rP = rank_f2(plaqs)
    r_inX  = rM - rank_f2([m & maskOut for m in stars])
    r_inZ  = rP - rank_f2([m & maskOut for m in plaqs])
    r_outX = rM - rank_f2([m & maskR   for m in stars])
    r_outZ = rP - rank_f2([m & maskR   for m in plaqs])
    r_in  = r_inX + r_inZ
    r_out = r_outX + r_outZ
    IR2 = (rM + rP) - r_in - r_out
    capP = 2 * AQ - r_in
    return dict(AQ=AQ, PER=PER, AV=AV, PD=PD, r_in=r_in, r_out=r_out,
                IR2=IR2, capP=capP, rM=rM, rP=rP)

def vset_mod(verts, L):
    # no caching: keying a cache by id(verts) risks stale hits after garbage
    # collection (error caught in-lane before sealing; see LEDGER.txt E1)
    return set(((i % L), (j % L)) for (i, j) in verts)

print("=" * 78)
print("LANE_T42_D_DERIVE  --  C-77 first increment: ONE derivation, exact checks")
print("=" * 78)

# ---------------- carrier sanity gates ----------------
for L in (6, 12):
    n, h, v, stars, plaqs = toric(L)
    ok_deg = all(sum(1 for m in stars if (m >> e) & 1) == 2 for e in range(n)) and \
             all(sum(1 for m in plaqs if (m >> e) & 1) == 2 for e in range(n))
    acc_s = 0
    for m in stars: acc_s ^= m
    acc_p = 0
    for m in plaqs: acc_p ^= m
    gate("carrier L=%d: each edge in exactly 2 stars and 2 plaqs" % L, ok_deg)
    gate("carrier L=%d: full star and plaq products vanish" % L, acc_s == 0 and acc_p == 0)
    gate("carrier L=%d: rank(stars)=rank(plaqs)=L^2-1" % L,
         rank_f2(stars) == L * L - 1 and rank_f2(plaqs) == L * L - 1)

# ---------------- CORNER TIER: s x s thick regions, L = 6,8,10,12 ----------------
print("\n-- CORNER TIER (pure Gamma: stabilizer boundary structure; exact F_2) --")
corner = {}   # (L,s) -> quantities
Ls = (6, 8, 10, 12)
for L in Ls:
    n, h, v, stars, plaqs = toric(L)
    for s in range(2, L - 1):
        verts, maskR = region_block(L, h, v, s)
        q = region_quantities(L, n, h, v, stars, plaqs, verts, maskR)
        corner[(L, s)] = q

# per-region structural gates (the theory step's exact verification)
ok_euler  = all(q["AV"] + q["PD"] - q["AQ"] == 4 * s - 3 for (L, s), q in corner.items())
ok_per    = all(q["PER"] == 4 * s for (L, s), q in corner.items())
ok_law2   = all(q["IR2"] == 2 * q["PER"] - 10 for (L, s), q in corner.items())
ok_rin    = all(q["r_in"] == (s - 2) ** 2 + (s - 1) ** 2 for (L, s), q in corner.items())
gate("corner: Euler count AV+PD-AQ == 4s-3 on every (L,s)  [bulk density 0: 1+1-2 per cell]", ok_euler)
gate("corner: PER == 4s on every (L,s)", ok_per)
gate("corner: IR2 == 2*PER-10 on every (L,s)  [reproduces T42_C LAW-2, thick scope]", ok_law2)
gate("corner: r_in == (s-2)^2+(s-1)^2 on every (L,s)", ok_rin)

# L-independence of the aggregates (locality of the boundary law)
ok_Lindep = True
for s in range(2, Ls[0] - 1):
    vals = set((corner[(L, s)]["IR2"], corner[(L, s)]["capP"]) for L in Ls if (L, s) in corner)
    if len(vals) != 1: ok_Lindep = False
gate("corner: (IR2, capP) independent of L at fixed s (shared s range)", ok_Lindep)

# translation invariance (in-lane homogeneity evidence; full Aut BORROWED, see D-22 lines)
L = 12
n, h, v, stars, plaqs = toric(L)
verts0, m0 = region_block(L, h, v, 4, 0, 0)
verts1, m1 = region_block(L, h, v, 4, 3, 4)
q0 = region_quantities(L, n, h, v, stars, plaqs, verts0, m0)
q1 = region_quantities(L, n, h, v, stars, plaqs, verts1, m1)
gate("corner D-22 evidence: IR2 and capP invariant under translation (offset (0,0) vs (3,4))",
     q0["IR2"] == q1["IR2"] and q0["capP"] == q1["capP"])
print("  D-22 (cited): venue Aut = 8L^2, C-74 ESTABLISHED (BORROWED, not recomputed here)")

# sequences at L=12, s=2..10
S = list(range(2, 11))
IR_seq  = [corner[(12, s)]["IR2"]  for s in S]
C_seq   = [corner[(12, s)]["capP"] for s in S]
EU_seq  = [corner[(12, s)]["AV"] + corner[(12, s)]["PD"] - corner[(12, s)]["AQ"] for s in S]
print("  s      :", S)
print("  IR2(s) :", IR_seq)
print("  capP(s):", C_seq)
dIR = fd_degree(IR_seq); dC = fd_degree(C_seq); dEU = fd_degree(EU_seq)
print("  computed degrees: deg IR2 = %s, deg capP = %s, deg Euler = %s" % (dIR, dC, dEU))
gate("corner: deg IR2 == 1 (computed by exact finite differences)", dIR == 1)
gate("corner: deg capP == 2 (computed; D_corner := 2 EARNED as content degree)", dC == 2)
gate("corner R1: deg IR2 == deg capP - 1", dIR == dC - 1)

# parameter-free exact relation between the two earned aggregates (no scale parameter s)
ok_pf = all(32 * q["capP"] == (q["IR2"] + 10) ** 2 + 8 * (q["IR2"] + 10) - 160
            for q in corner.values())
gate("corner R1 exact form: 32*capP == (IR2+10)^2 + 8*(IR2+10) - 160 on every (L,s)", ok_pf)

# D-15 control (a): wrong-form candidate IR^D ~ C^D, i.e. IR2^2/capP constant -- must FAIL
ratios = [Fraction(ir * ir, c) for ir, c in zip(IR_seq, C_seq)]
gate("corner D-15 control: IR2^2/capP constant FAILS (witness ratios %s vs %s)"
     % (ratios[0], ratios[-1]), len(set(ratios)) > 1)

# D-15 control (b): scattered regions (same cardinality, no interior) -- law must FAIL
print("\n  scatter control (L=12, deterministic stride spread, |R| = 2s(s-1)):")
IRs_seq, Cs_seq = [], []
for s in S:
    m = 2 * s * (s - 1)
    idxs = sorted(set((k * n) // m for k in range(m)))
    ok_card = (len(idxs) == m)
    maskR = 0
    for e in idxs: maskR |= (1 << e)
    vertsS = set()
    for e in idxs:
        a, b = edge_endpoints(12, e)
        vertsS.add((a[0] % 12, a[1] % 12)); vertsS.add((b[0] % 12, b[1] % 12))
    q = region_quantities(12, n, h, v, stars, plaqs, vertsS, maskR)
    if not ok_card:
        gate("scatter s=%d: cardinality control" % s, False)
    IRs_seq.append(q["IR2"]); Cs_seq.append(q["capP"])
print("  IR2_scat :", IRs_seq)
print("  capP_scat:", Cs_seq)
dIRs = fd_degree(IRs_seq); dCs = fd_degree(Cs_seq)
print("  computed degrees: deg IR2_scat = %s, deg capP_scat = %s" % (dIRs, dCs))
gate("corner D-15 control: boundary-degree law FAILS on scatter (deg IR != deg C - 1)",
     not (dIRs is not None and dCs is not None and dIRs == dCs - 1))

# ---------------- WORLD TIER: census access model (coupling term, openness clause) ----------
print("\n-- WORLD TIER (census barrier-record blocks; access counts; exact integers) --")
def world_counts(nn):
    IFACE = 0; CAP_d1 = 0
    for x in range(nn):
        for y in range(nn):
            for z in range(nn):
                out = 0
                for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
                    X, Y, Z = x+dx, y+dy, z+dz
                    if not (0 <= X < nn and 0 <= Y < nn and 0 <= Z < nn): out += 1
                IFACE += out
                if out: CAP_d1 += 1
    return IFACE, CAP_d1, nn ** 3

Ns = list(range(2, 13))
IF_seq, D1_seq, CT_seq = [], [], []
for nn in Ns:
    IF, D1, CT = world_counts(nn)
    IF_seq.append(IF); D1_seq.append(D1); CT_seq.append(CT)
print("  n        :", Ns)
print("  IFACE(n) :", IF_seq)
print("  CAP_d1(n):", D1_seq)
print("  CAP_tot  :", CT_seq)
gate("world: IFACE == 6n^2 for all n=2..12 (recount matches T42_C established)",
     all(IF_seq[i] == 6 * Ns[i] ** 2 for i in range(len(Ns))))
gate("world: CAP_d1 == n^3-(n-2)^3 for all n (bulk cancellation by nested difference)",
     all(D1_seq[i] == Ns[i] ** 3 - (Ns[i] - 2) ** 3 for i in range(len(Ns))))
dIF = fd_degree(IF_seq); dCT = fd_degree(CT_seq)
print("  computed degrees: deg IFACE = %s, deg CAP_total = %s" % (dIF, dCT))
gate("world: deg IFACE == 2 (computed)", dIF == 2)
gate("world: deg CAP_total == 3 (computed; D_world := 3 EARNED as content degree)", dCT == 3)
gate("world R1: deg IFACE == deg CAP_total - 1", dIF == dCT - 1)
gate("world R1 exact form: IFACE^3 == 216 * CAP_total^2 for every n",
     all(IF_seq[i] ** 3 == 216 * CT_seq[i] ** 2 for i in range(len(Ns))))

# D-15 control (c): wrong-form candidate IFACE^2/CAP_total constant -- must FAIL
wr = [Fraction(IF_seq[i] ** 2, CT_seq[i]) for i in range(len(Ns))]
gate("world D-15 control: IFACE^2/CAP_total constant FAILS (witness %s vs %s)"
     % (wr[0], wr[-1]), len(set(wr)) > 1)

# D-15 control (d): scattered grains (every-other site in (2n)^3) -- law must FAIL
def world_scatter(nn):
    IFACE = 0
    occ = set((2*x, 2*y, 2*z) for x in range(nn) for y in range(nn) for z in range(nn))
    for (x, y, z) in occ:
        for dx, dy, dz in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
            if (x+dx, y+dy, z+dz) not in occ: IFACE += 1
    return IFACE, nn ** 3
Nss = list(range(2, 9))
IFs = []; CTs = []
for nn in Nss:
    a, b = world_scatter(nn); IFs.append(a); CTs.append(b)
dIFs = fd_degree(IFs); dCTs = fd_degree(CTs)
print("  scatter: IFACE_scat:", IFs, " degrees:", dIFs, dCTs)
gate("world D-15 control: boundary-degree law FAILS on scattered grains (deg IF == deg C == 3)",
     not (dIFs is not None and dCTs is not None and dIFs == dCTs - 1))
print("  D-22 (cited): world venue checks (hypercube Aut, grid D4) per T42_B/T42_C (BORROWED)")

# ---------------- summary ----------------
print("\n" + "=" * 78)
npass = sum(1 for _, ok in GATES if ok)
overall = all(ok for _, ok in GATES)
print("GATES: %d/%d passed" % (npass, len(GATES)))
print(("OVERALL PASS" if overall else "OVERALL FAIL") + "  (computed conjunction of %d gates)" % len(GATES))
print("Largest exact objects: toric L=12 (n=288 qubits, F_2 ranks on 144-row bitmask matrices,")
print("  9-point s-sequences); world n=12 (1728 grains); scatter control 512 grains in 16^3.")

print("""
============================ FINAL COMPARISON SECTION ============================
(The only place classical gravitational relations are named. No computation here.)

DERIVED (construction above, no classical form used anywhere):
  R1: For thick nested regions the Gamma boundary aggregate (interface rank) and the
      content aggregate satisfy  deg IR = deg C - 1  with D := deg C the EARNED
      dimension (computed content degree: 2 corner, 3 world); parameter-free exact
      forms  32*C = (IR+10)^2 + 8*(IR+10) - 160  (corner)  and  IR^3 = 216*C^2
      (world), i.e. IR ~ kappa^(1/D) * C^((D-1)/D): sub-extensive boundary scaling of
      the record surface's coupling structure, produced by exact bulk cancellation
      (Euler density 0 per interior cell; nested-difference cancellation), and
      requiring an interior (both scatter controls break the law).

RECOVERY TARGETS AND VERDICTS:
  (a) Area / entropy-boundary law:  S_BH = A / (4 l_P^2)
      OWNERS: Bekenstein (1973), Hawking (1975); holographic bound: 't Hooft (1993),
      Susskind (1995); entanglement-entropy area law: Bombelli-Koul-Lee-Sorkin (1986),
      Srednicki (1993); stabilizer/perimeter instance: Hamma-Ionicioiu-Zanardi (2005),
      Kitaev-Preskill (2006).
      VERDICT: MATCHES IN FORM. The derived aggregate scales as the boundary measure,
      one degree below content -- the form of the area law -- and the matching quantity
      is the interface/entanglement rank, exactly the quantity the entanglement reading
      assigns to horizon entropy. NOT matched: stored record content (volume in both
      tiers; the strong-holographic reading is not delivered, per T42_C -- admissible
      under C-77's different-shape expectation). The corner instance is in substance
      the known stabilizer perimeter law (BORROWED, Hamma-Ionicioiu-Zanardi); OURS is
      the exact constants and scope, the Euler bulk-cancellation route, the two-tier
      degree split with dimension earned as computed degree, and the parameter-free
      IR-C forms. The coefficient is NOT recovered: kappa is venue-dependent (32 vs
      216); no 1/4G analogue exists on this surface.
  (b) Potential / falloff laws (Newton: F ~ 1/r^2, phi ~ 1/r).
      VERDICT: NEITHER-YET from Gamma. T42_B's measured exponent 3 is the world
      model's own image-construction law (measured, not derived from record
      structure); no falloff relation has been derived from Gamma.
  (c) Horizon thermodynamics (Hawking temperature; Bardeen-Carter-Hawking first law).
      VERDICT: NEITHER-YET. No temperature-like or energy-like conjugate pair has
      been earned on the record surface.

C-77 CLAUSE APPLIED: one scale-level relation was derived from record-surface boundary
structure without assuming any classical form in construction, and it MATCHES IN FORM a
known gravitational relation (a). Openness clause invoked: the world-tier instance
runs through the census access/coupling model (a coupling term alongside Gamma), which
the claim explicitly allows. Increment verdict: SATISFIED (form-level, scope as stated).

NEXT STEP (binding, per discipline): derive a falloff relation from Gamma alone --
two-region interface rank as a function of earned separation (T42_A metric / T42_B
positions) -- the candidate for increment (b); and search for an earned conjugate pair
before touching (c).
==================================================================================""")
sys.exit(0 if overall else 1)
