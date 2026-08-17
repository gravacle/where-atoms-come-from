"""
F1 — EMBED.  What parameter space are "STOP" and "STOP-FALLS-REBUILD" points in?

Method.  A disposition in this program is not a scalar.  Every place the corpus
actually USES one, it uses a tuple.  We read the coordinates off the corpus
itself (never invented here), machine-verify each coordinate's VALUES against a
file:line pointer inside the cutoff, and then place the two named arms.

Coordinate O  — OBJECT VERDICT      (is the forced crossing dead?)
Coordinate G  — GROUND              (on what does the verdict rest?)
Coordinate B  — BUILD TARGET        (what gets built or investigated next?)
Coordinate P  — PUBLICATION SET     (which surviving results are published?)
Coordinate A  — SELF-AUDIT STATUS   (has this layer been checked from outside its lineage?)

FRAME_CHALLENGE §0 defines a route decision as one "whose answer determines what
gets built or investigated next" — i.e. the route coordinate is B.  So B is
measured twice: as a single choice (|B| values) and as a subset (2^|B|).
"""
import itertools
from cutoff_guard import CUT, line, TEXT

REPO = "/Users/bgm/MB Work/where-atoms-come-from"

# ---------------------------------------------------------------- B, verified
# (label, pointer, substring that must appear at that pointer)
B_VALUES = [
 ("b00 none / no further build",
  "REGISTER:448", "RECOMMENDATION OF RECORD: STOP"),
 ("b01 wedge-growth sequence V=4k+1 on the dressed algebra",
  "REGISTER:620", "wedge-growth sequence"),
 ("b02 edge-by-edge transport T with T^3 = M_gamma (COR-F)",
  "REGISTER:530", "edge by edge"),
 ("b03 the charge run",
  "REGISTER:247", "the charge run"),
 ("b04 the SU(2) run / gauge group varied independently",
  "REGISTER:247", "the SU(2) run"),
 ("b05 a third loop",
  "REGISTER:247", "a third loop"),
 ("b06 a third schedule",
  "REGISTER:247", "a third schedule"),
 ("b07 audit S1 (the one artifact no adversary read)",
  "REGISTER:292", "S1 IS UNAUDITED"),
 ("b08 the lineage-independent lane",
  "REGISTER:87", "a lineage-independent lane"),
 ("b09 carrier-intrinsic schedule condition",
  "REGISTER:155", "carrier-intrinsic rather than adversarial"),
 ("b10 ready-state support derived from the incidence",
  "REGISTER:156", "shown to follow from the incidence"),
 ("b11 convex-hull criterion on another carrier",
  "REGISTER:88", "fails on a carrier where it should hold"),
 ("b12 three-way split derivability vs vertex count",
  "REGISTER:88", "three-way split's"),
 ("b13 measure the recurrence obstruction directly",
  "REGISTER:574", "The recurrence obstruction is the one thing"),
 ("b14 S5 — look for the free dimensionless coupling (alpha's slot)",
  "FOUNDING_DESIGN_V001.md:108", "S5 — LOOK FOR THE SLOT"),
 ("b15 publish the narrow note with its placements",
  "REGISTER:344", "RECOMMENDATION OF RECORD: **NARROW**"),
 ("b16 does the rate vary with the connection at S4",
  "REGISTER:157", "fails to vary with the connection at S4"),
]

def verify(ptr, needle):
    f, n = ptr.rsplit(":", 1)
    n = int(n)
    if f == "REGISTER":
        return needle in line(n)
    txt = open(f"{REPO}/{f}", encoding="utf-8").read().replace("\r\n", "\n").split("\n")
    return needle in txt[n - 1]

print("F1 — COORDINATE B (BUILD TARGET): NAMED VALUES, EACH POINTER-VERIFIED")
print("-" * 78)
ok = 0
for lab, ptr, needle in B_VALUES:
    good = verify(ptr, needle)
    ok += good
    print(f"  [{'OK ' if good else 'MISS'}] {lab:<58} {ptr}")
print(f"\n  pointer-verified values of B : {ok} / {len(B_VALUES)}")
assert ok == len(B_VALUES), "an enumerated build target failed its pointer"

# ---------------------------------------------------------------- G, verified
G_VALUES = [
 ("g0 the object does not exist (no forced crossing)", "REGISTER:450",
  "does not exist"),
 ("g1 the spine is VACUOUS (IMP-1 disqualifier)",      "REGISTER:408",
  "IS VACUOUS, WHICH IS WORSE THAN FALSE"),
 ("g2 the spine is FALSE with an exhibited counterexample", "REGISTER:571",
  "RESTORED FROM VACUOUS TO FALSE"),
 ("g3 the recurrence obstruction, undented",           "REGISTER:574",
  "The recurrence obstruction is the one thing"),
 ("g4 the disposition rests on an unpointered import", "REGISTER:618",
  "issued on an unpointered import"),
]
print("\nF1 — COORDINATE G (GROUND): NAMED VALUES")
print("-" * 78)
for lab, ptr, needle in G_VALUES:
    print(f"  [{'OK ' if verify(ptr, needle) else 'MISS'}] {lab:<58} {ptr}")

# ---------------------------------------------------------------- P, verified
P_VALUES = [("N1 rate = logarithmic Mahler measure", "REGISTER:456", "N1. THE RATE"),
            ("N2 multiset / invisibility theorem",   "REGISTER:462", "N2. THE MULTISET"),
            ("N3 the null, inverted",                "REGISTER:465", "N3. THE NULL"),
            ("N4 vacuity/fibre-wise-ness of the spine", "REGISTER:470", "N4. THE VACUITY")]
print("\nF1 — COORDINATE P (PUBLICATION SET): NAMED ELEMENTS")
print("-" * 78)
for lab, ptr, needle in P_VALUES:
    print(f"  [{'OK ' if verify(ptr, needle) else 'MISS'}] {lab:<58} {ptr}")

# ---------------------------------------------------------------- O and A
print("\nF1 — COORDINATE O (OBJECT VERDICT): {alive, dead}")
print("      both arms take the SAME value: dead.")
print(f"      W-05 REGISTER:450  ...{line(450).strip()[:66]}")
print(f"      W-06 REGISTER:615  ...{line(615).strip()[:66]}")
print("\nF1 — COORDINATE A (SELF-AUDIT): {checked outside lineage, not checked}")
print(f"      custody §4 grades this corpus adversarially-checked, never independently-corroborated.")
print(f"      W-06 REGISTER:635  ...{line(635).strip()[:66]}")

# ---------------------------------------------------------------- cardinality
nB, nG, nP, nO, nA = len(B_VALUES), len(G_VALUES), 2 ** len(P_VALUES), 2, 2
card_single = nB * nG * nP * nO * nA
card_subset = (2 ** (nB - 1)) * nG * nP * nO * nA   # b00 = the empty subset

print("\n" + "=" * 78)
print("THE SPACE")
print("=" * 78)
print(f"  D = O x G x B x P x A       dimension (independent coordinates) = 5")
print(f"  |O| = {nO}   |G| = {nG}   |B| = {nB} (single choice)   |P| = 2^4 = {nP}   |A| = {nA}")
print(f"  |D| with B a single choice  = {card_single}")
print(f"  |D| with B a subset         = {card_subset}")
print(f"\n  the binary names 2 points.")
print(f"  coverage of D                        = 2/{card_single} = {2/card_single:.3e}")
print(f"  coverage of the ROUTE coordinate B   = 2/{nB} = {2/nB:.4f}")
print(f"  named B-values claimed by NEITHER arm = {nB - 3} of {nB} "
      f"(b00, b01, b15 are the only ones either arm touches)")

print("\nARM COORDINATES")
print("-" * 78)
print("  STOP  (W-05, REGISTER:448-472) =")
print("      O = dead            REGISTER:450")
print("      G = g0 + g1         REGISTER:408, :450")
print("      B = {b15}           REGISTER:456-472  (publish N1..N4 with placements)")
print("      P = {N1,N2,N3,N4}   REGISTER:456-472")
print("      A = not checked     custody §4")
print("  STOP-FALLS-REBUILD (W-06, REGISTER:613-624) =")
print("      O = dead ('DEADER') REGISTER:615")
print("      G = g2 + g4         REGISTER:571, :617-618")
print("      B = {b01} (+ b15 not withdrawn)  REGISTER:619-620")
print("      P = {N1,N2,N3,N4'}  N4 corrected at REGISTER:577, nothing struck")
print("      A = not checked     REGISTER:635 'Discount this layer too.'")
print("\n  The two arms are EQUAL on O, EQUAL on A, EQUAL on P up to one relabel,")
print("  and differ on exactly TWO of five coordinates: G and B.")
