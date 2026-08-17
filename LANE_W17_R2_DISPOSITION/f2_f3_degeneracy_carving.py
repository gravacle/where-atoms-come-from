"""
F2 — DEGENERACY  and  F3 — CARVING, for R2 (STOP  vs  STOP-FALLS-REBUILD).

F2 asks: are the arms the same object under a map?  We build the two arms as
status vectors over the corpus assets each arm rules on, and compute the
normalised Hamming distance under the identity and under four projections that
are the disposition-space analogues of the instrument's named candidate maps
(restriction, change of coordinates, forgetting a coordinate, power/absorption).

F3 asks: does the predicate partition anything?  We measure both cells under two
stated measures:
  mu_B      counting measure on the 17 pointer-verified named build targets
  mu_reopen counting measure on the register's OWN registered reopen conditions
            (W-01 x3, W-02 x4, W-03 x4 = 11; W-04/05/06 register none)
"""
import numpy as np
from cutoff_guard import line

# ---------------------------------------------------------------------------
# THE ASSET TABLE.  status codes: 2 = alive/licensed, 1 = published-as-result,
# 0 = dead/withdrawn.  Every row carries the in-cut pointer it was read from.
# ---------------------------------------------------------------------------
ASSETS = [
 # name                                   STOP  SFR   pointer(STOP)   pointer(SFR)
 ("forced crossing out of the carrier",      0,   0,  "REG:450",  "REG:615"),
 ("Theorem S3-0 as stated",                  0,   0,  "REG:408",  "REG:608"),
 ("N1 rate = Mahler measure",                1,   1,  "REG:456",  "REG:606-611 (not in WHAT STAYS DEAD)"),
 ("N2 multiset / invisibility theorem",      1,   1,  "REG:462",  "REG:462 (survived this round)"),
 ("N3 the null, inverted",                   1,   1,  "REG:465",  "REG:606-611 (not in WHAT STAYS DEAD)"),
 ("N4 mechanism of the spine",               1,   1,  "REG:470",  "REG:577 (relabelled)"),
 ("novelty of the nine",                     0,   0,  "REG:482",  "REG:611"),
 ("Bell-1975 attribution",                   1,   0,  "REG:482",  "REG:584"),
 ("recurrence obstruction",                  2,   2,  "REG:414",  "REG:573-575"),
 ("a next construction is licensed",         0,   2,  "REG:448",  "REG:619-620"),
]
names  = [a[0] for a in ASSETS]
A = np.array([a[1] for a in ASSETS])   # STOP
B = np.array([a[2] for a in ASSETS])   # STOP-FALLS-REBUILD

def hamming(x, y):
    return float(np.mean(x != y))

print("=" * 78)
print("F2 — DEGENERACY.  ||A - f(B)|| over the natural candidate maps f")
print("=" * 78)
print(f"{'asset':<40}{'STOP':>6}{'SFR':>6}   differs")
for n, a, b in zip(names, A, B):
    print(f"{n:<40}{a:>6}{b:>6}   {'X' if a != b else ''}")

d_id = hamming(A, B)
print(f"\n  f = identity                       ||A-B||_Hamming/n = {d_id:.4f}"
      f"   ({int(d_id*len(A))} of {len(A)} assets differ)")

# --- map 1: pi_halt.  the projection the LABELS advertise: halt vs continue.
halt = lambda v: np.array([1 if (v[2] or v[3] or v[4] or v[5] or v[9]) else 0])
print(f"  f = pi_halt   (halt? / continue?)  ||.|| = "
      f"{hamming(halt(A), halt(B)):.4f}   BOTH ARMS = CONTINUE")
print("        W-05 REGISTER:451  'STOP is not \"burn it\"' — STOP licenses N1..N4 + placements.")
print("        The pair is named stop-vs-rebuild and is in fact build-X vs build-Y.")

# --- map 2: pi_object.  the coordinate the disposition claims to settle.
print(f"  f = pi_object (is the object dead?)||.|| = "
      f"{hamming(A[:1], B[:1]):.4f}   BOTH ARMS = DEAD")
print("        REGISTER:615 'STOP's OBJECT is correctly dead and the restorations make it DEADER.'")

# --- map 3: pi_publish.  restriction to the publication coordinate.
print(f"  f = pi_publish (which results ship) ||.|| = "
      f"{hamming(A[2:6], B[2:6]):.4f}   IDENTICAL up to N4's relabel (REGISTER:577)")

# --- map 4: pi_build.  the coordinate FRAME_CHALLENGE §0 says defines a route.
print(f"  f = pi_build  (what gets built)     ||.|| = "
      f"{hamming(A[9:10], B[9:10]):.4f}   b15 (publish) vs b01 (wedge-growth) — GENUINELY DIFFERENT")

# --- map 5: absorption.  is SFR = STOP + one reopen?  (nested, not opposed)
withdrawn_by_SFR = [n for n, a, b in zip(names, A, B) if a >= 1 and b == 0]
added_by_SFR     = [n for n, a, b in zip(names, A, B) if b > a]
print(f"\n  f = absorption (SFR =? STOP + reopen)")
print(f"        assets STOP licenses that SFR WITHDRAWS : {len(withdrawn_by_SFR)}  {withdrawn_by_SFR}")
print(f"        assets SFR ADDS to STOP                 : {len(added_by_SFR)}  {added_by_SFR}")
print("        => on the operative coordinates SFR is NESTED OVER STOP, not opposed to it.")
print("           The one true withdrawal is a CITATION (Bell->Hepp), not a disposition.")

print("\n  F2 VERDICT: the arms are NOT identical — they differ on pi_build, 1 of 5")
print("  coordinates.  But 3 of the 5 tested maps return distance EXACTLY 0, including")
print("  the halt/continue projection the two LABELS advertise.  The label pair is")
print("  degenerate; the underlying pair is not.")

# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("F3 — CARVING.  measure of each cell")
print("=" * 78)

nB = 17                       # pointer-verified in f1_embed.py
cellA = {"b15"}               # STOP  = publish the narrow note
cellB = {"b01", "b15"}        # SFR   = wedge-growth, publication not withdrawn
inter = cellA & cellB
union = cellA | cellB
print(f"  mu_B = counting measure on the {nB} pointer-verified named build targets")
print(f"    mu_B(cell STOP) = {len(cellA)}/{nB} = {len(cellA)/nB:.4f}")
print(f"    mu_B(cell SFR ) = {len(cellB)}/{nB} = {len(cellB)/nB:.4f}")
print(f"    OVERLAP  mu_B(A n B) = {len(inter)}/{nB} = {len(inter)/nB:.4f}"
      f"   = {len(inter)/len(cellA)*100:.0f}% of cell A, {len(inter)/len(cellB)*100:.0f}% of cell B")
print(f"    UNION    mu_B(A u B) = {len(union)}/{nB} = {len(union)/nB:.4f}")
print(f"    UNCLAIMED COMPLEMENT = {nB - len(union)}/{nB} = {(nB-len(union))/nB:.4f}")
print("    => the cells are NESTED (A subset B), so the predicate does not partition:")
print("       'hold at STOP' is contained in 'accept STOP-FALLS-REBUILD'.")

# the register's own measure
REOPENS = [
 ("W-01", "REG:87-90",  ["lineage-independent lane fails to reproduce the firing",
                          "convex-hull criterion fails on another carrier",
                          "three-way split derivability vs vertex count"]),
 ("W-02", "REG:155-158",["carrier-intrinsic schedule condition",
                          "ready-state support follows from the incidence",
                          "rate fails to vary with the connection at S4",
                          "lineage-independent lane fails the character-ratio criterion"]),
 ("W-03", "REG:247-248",["the charge run", "the SU(2) run",
                          "a third loop", "a third schedule"]),
]
allre = [r for _, _, rs in REOPENS for r in rs]
print(f"\n  mu_reopen = counting measure on the register's OWN registered reopen conditions")
for tag, ptr, rs in REOPENS:
    print(f"    {tag} {ptr}: {len(rs)}")
print(f"    total registered and unfired = {len(allre)}")
print(f"    W-04, W-05, W-06 register ZERO reopen conditions "
      f"(custody §8: 'a row with no reopen condition is closed permanently')")
fired_by_A = 0
fired_by_B = 0
print(f"    mu_reopen(cell STOP) = {fired_by_A}/{len(allre)} = 0.0000")
print(f"    mu_reopen(cell SFR ) = {fired_by_B}/{len(allre)} = 0.0000")
print("    => BOTH CELLS ARE MEASURE ZERO under the register's own live-question measure.")
print("       Neither arm fires, discharges, or even names any of the 11 conditions the")
print("       register itself says are the ones that would move a row.")
print("\n  F3 VERDICT: NON-CARVING on both measures — nested cells under mu_B,")
print("  two measure-zero cells under mu_reopen.")
