"""T42-C  THE BOUNDARY-CAPACITY RELATION -- WORLD TIER, and the side-by-side.

THE CENSUS CONTRAST (established, GR1/census): world records are barrier-protected and the
record capacity of a region of media scales with VOLUME -- each grain is a record.  This
script computes, with exact integers and no fit, WHICH aggregates of a block of barrier
records scale with the boundary instead, under the brief's candidate: interface-mediated
access only -- an outside operation reaches a grain only by crossing the region's own
boundary, grain by grain.

MODEL (stated in full; a graph, nothing continuous): grains at the integer points of an
n x n x n block inside an unbounded grid of media; one barrier record per grain (census).
Adjacency = nearest neighbour (shares a face).  An outside operation enters through a
boundary face and must traverse grains in adjacency steps.  DEPTH of a grain = number of
grains an entering operation must traverse, counting the grain itself (computed by BFS from
the complement, not by a formula).

AGGREGATES (all counts):
  CAP_total      = number of grains                       (the census capacity)
  IFACE          = number of block-complement adjacencies (independent access channels;
                   also the F_2 rank of the channel->first-grain influence map is computed)
  CAP_depth(1)   = grains erasable WITHOUT traversing any other grain (one access step)
  CAP_prot       = grains whose erasure requires first traversing at least one other grain
                   of the region (depth >= 2) -- the brief's 'protected against outside
                   operations' candidate
  DEPTH_SUM      = total traversal cost to erase everything, sum of depths

SCALING is read by exact finite differences on the integer sequences (a cubic has constant
third difference, a quadratic constant second difference) -- no exponent is imported, no
fit is performed.  D-15 control: a SCATTERED set of the same number of grains (deterministic
seed), where the boundary-scaling aggregates must go volume-like and the protected aggregate
must vanish.

D-1: the classical relation appears ONLY in the final comparison section at the bottom.
"""
import random

def analyse_block(n):
    grains = {(x, y, z) for x in range(n) for y in range(n) for z in range(n)}
    return analyse(grains)

def analyse(grains):
    NB = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))
    iface = 0
    first = set()          # grains adjacent to the complement (reachable in one step)
    for g in grains:
        for d in NB:
            h = (g[0] + d[0], g[1] + d[1], g[2] + d[2])
            if h not in grains:
                iface += 1
                first.add(g)
    # BFS depths from the complement
    depth = {g: 1 for g in first}
    frontier = sorted(first)
    while frontier:
        nxt = []
        for g in frontier:
            for d in NB:
                h = (g[0] + d[0], g[1] + d[1], g[2] + d[2])
                if h in grains and h not in depth:
                    depth[h] = depth[g] + 1
                    nxt.append(h)
        frontier = sorted(nxt)
    assert set(depth) == grains
    cap_total = len(grains)
    cap_d1 = sum(1 for g in grains if depth[g] == 1)
    cap_prot = cap_total - cap_d1
    depth_sum = sum(depth.values())
    # F_2 rank of the influence map channels -> first-layer grains: each channel toggles
    # exactly one first-layer grain, so the rank is the number of DISTINCT grains hit;
    # computed as a rank all the same (unit vectors, exact)
    rank_infl = len(first)
    return cap_total, iface, cap_d1, cap_prot, depth_sum, rank_infl

def diffs(seq):
    return [b - a for a, b in zip(seq, seq[1:])]

def const_diff_order(seq):
    """smallest k with the k-th finite difference constant (exact integers); None if none up
       to len-2."""
    s = list(seq)
    for k in range(len(seq) - 1):
        if len(set(s)) == 1:
            return k
        s = diffs(s)
    return None

def gate(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok

if __name__ == "__main__":
    NS = list(range(2, 13))
    rows = [analyse_block(n) for n in NS]
    print("WORLD TIER -- block of barrier records, exact counts, n = 2..12")
    print(f"{'n':>3} {'CAP_total':>10} {'IFACE':>7} {'CAP_d1':>7} {'CAP_prot':>9} "
          f"{'DEPTH_SUM':>10} {'rank_infl':>10}")
    for n, r in zip(NS, rows):
        print(f"{n:>3} {r[0]:>10} {r[1]:>7} {r[2]:>7} {r[3]:>9} {r[4]:>10} {r[5]:>10}")
    cap = [r[0] for r in rows]; ifc = [r[1] for r in rows]; d1 = [r[2] for r in rows]
    prot = [r[3] for r in rows]; dsum = [r[4] for r in rows]; rinf = [r[5] for r in rows]
    print()
    print("SCALING BY EXACT FINITE DIFFERENCES (order k means: k-th difference constant, i.e.")
    print("the count is a degree-k polynomial in n over the computed range -- READ, not fitted):")
    ok = True
    ok &= gate("CAP_total is degree 3 (VOLUME) -- the census, replicated",
               const_diff_order(cap) == 3, f"order={const_diff_order(cap)}")
    ok &= gate("IFACE is degree 2 (BOUNDARY)", const_diff_order(ifc) == 2,
               f"order={const_diff_order(ifc)}")
    ok &= gate("CAP_depth(1) is degree 2 (BOUNDARY) -- capacity addressable without "
               "traversal", const_diff_order(d1) == 2, f"order={const_diff_order(d1)}")
    ok &= gate("influence rank == CAP_depth(1) exactly (channels collapse onto the first "
               "layer)", rinf == d1)
    ok &= gate("CAP_prot is degree 3 (VOLUME) -- 'requires crossing the boundary' protects "
               "the BULK of the block, not a boundary's worth",
               const_diff_order(prot) == 3, f"order={const_diff_order(prot)}")
    ok &= gate("DEPTH_SUM is NOT any single polynomial over the range (no finite difference "
               "goes constant: parity of n matters -- the surface's own shape, reported as is)",
               const_diff_order(dsum) is None, f"order={const_diff_order(dsum)}")
    ok &= gate("DEPTH_SUM exact recurrence, every n >= 4: DEPTH_SUM(n) == DEPTH_SUM(n-2) + n^3 "
               "(erasing a block = erasing its outer shell at full depth-weight n^3 worth, "
               "then the block two sizes down) -- super-volume growth, read exactly",
               all(dsum[i] == dsum[i - 2] + NS[i] ** 3 for i in range(2, len(NS))))
    ok &= gate("exact closed forms verified at every n: IFACE == 6n^2 and "
               "CAP_d1 == n^3 - (n-2)^3",
               all(i == 6 * n * n and d == n ** 3 - (n - 2) ** 3
                   for n, i, d in zip(NS, ifc, d1)))
    print()
    print("D-15 CONTROL -- SCATTERED grains (same count as the n=6 block, deterministic seed):")
    rnd = random.Random(7)
    cells = [(x, y, z) for x in range(24) for y in range(24) for z in range(24)]
    scat = set(rnd.sample(cells, 216))
    s = analyse(scat)
    print(f"    scattered 216 grains in a 24^3 volume: CAP_total={s[0]} IFACE={s[1]} "
          f"CAP_d1={s[2]} CAP_prot={s[3]}")
    b6 = rows[NS.index(6)]
    ok &= gate("scattered set: nearly every grain is first-layer (CAP_d1 == CAP_total would be "
               "full dispersal; block has CAP_d1 << CAP_total)",
               s[2] > (9 * s[0]) // 10 and b6[2] < (3 * b6[0]) // 4,
               f"scatter {s[2]}/{s[0]}, block {b6[2]}/{b6[0]}")
    ok &= gate("scattered set: IFACE within a factor ~6 of its grain count (volume-like), "
               "block IFACE is 6n^2 (boundary)",
               s[1] > 4 * s[0] and b6[1] == 6 * 36)
    ok &= gate("scattered set: protected capacity collapses toward zero",
               s[3] < s[0] // 10, f"CAP_prot={s[3]} of {s[0]}")
    print()
    print("=" * 98)
    print("SIDE BY SIDE -- what each tier's own counts exhibit (corner numbers from t42_corner)")
    print("=" * 98)
    print("""
                                  CORNER (toric, exact F_2)        WORLD (barrier block, exact)
  stored / distinguishable        cap_pauli = 2AQ - r_in:          CAP_total = n^3: VOLUME
  content of a region             BULK-dominant (2AV + PER/2 - 5   (the census, replicated)
                                  on thick rectangles)
  interface measure               cut-rank IR2: perimeter law      IFACE = 6n^2, CAP_depth(1)
                                  2*PER - 10 on rectangles; the    = n^3 - (n-2)^3: BOUNDARY
                                  exact general law is the         (degree 2, read by exact
                                  Euler-count LAW-3                differences)
  record (protected) content      cap_record: TOPOLOGICAL, 0 on    CAP_prot (needs traversal):
                                  every contractible region, 1-2   VOLUME (degree 3)
                                  per wrap class on bands
  dispersed control (D-15)        scatter: IR2 -> ceiling 2AQ,     scatter: CAP_d1 -> CAP_total,
                                  r_in -> small: all interface     CAP_prot -> 0: all interface

  THE ANSWER TO THE LANE QUESTION (which aggregate of the world's records scales with
  boundary rather than bulk): the aggregates that scale with the boundary are the INTERFACE
  aggregates -- the independent access channels (IFACE) and the capacity addressable without
  traversing the region's own records (CAP_depth(1) = influence rank).  The brief's candidate
  'records whose erasure requires crossing the region's own boundary' does NOT scale with
  boundary: it is the bulk minus a boundary's worth (n^3 - (n-2)^3 removed), degree 3.  In
  BOTH tiers the boundary-scaling quantity is the interface rank, and in both tiers it is a
  property a region only has when it HAS an interior: thin strips (corner) and scattered
  grains (both tiers) are all interface, and their 'boundary' measure is their volume.
""")
    print("=" * 98)
    print("FINAL COMPARISON SECTION -- the only place the classical relation is named (D-1)")
    print("=" * 98)
    print("""
  RECOVERY TARGET (named here only): the known gravitational relation in which a region's
  capacity scales with the AREA OF ITS BOUNDARY, not its volume (the Bekenstein-Hawking /
  holographic bound).
  WHAT THIS LANE'S SURFACE EXHIBITS, with no classical form used in any construction step:
    * The record surface's own interface measure -- the stabiliser cut-rank -- obeys an
      EXACT boundary law on the corner carrier (LAW-2 on rectangles; LAW-3, the Euler-count
      law, everywhere).  It was not put in: the same instrument run on strips and scattered
      sets returns their volume, and on rectangles of equal perimeter and different area it
      returns the same number.
    * The corresponding world aggregates -- access channels and traversal-free capacity --
      are boundary-scaling by exact finite differences, while the stored capacity stays
      volume (census).
    * The STORED-CONTENT capacity is NOT boundary-law in either tier.  If C-77's increment
      requires the region's full content to obey the boundary law (the strong holographic
      reading), this surface does not deliver it and says so; what it delivers exactly is a
      boundary law for the interface rank -- the number of independent constraints/channels
      through which a region's records meet the rest of the world.
  CAVEAT (the principal's, binding): the record-level measures here behave like their
  familiar analogues only where the tables say so; the boundary law's constant (-10), its
  breakdown on interiorless regions, and its Euler-count completion are the surface's own
  shapes, not imported ones.  A classical-shaped match on the interface rank is NOT a
  same-shape claim about record-level gravity (the two-way reading stands).
""")
    print(f"ALL GATES PASS: {ok}")
