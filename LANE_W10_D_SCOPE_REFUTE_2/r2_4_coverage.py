# W10-D REFUTE-2  LENS 2 = COMPLETENESS.  LEG 4.
#
# A MECHANICAL COVERAGE SCAN, so that "the table omits X" is a checkable statement and not an
# impression.  For each load-bearing claim I read out of REGISTER_V001.md (and, where the
# register points at it, out of the sealed artifact it points at), a DISTINCTIVE TOKEN is
# chosen -- a string that must appear if the claim is being discussed at all -- and searched
# for in TWO places:
#     (1) W10D_SCOPE_TABLE_V001.md          -- the artifact the principal rules on
#     (2) the whole LANE_W10_D_SCOPE/ dir   -- code and captured output included, to be fair
#
# THE TOKEN IS NOT THE CLAIM.  A hit does not prove the claim was scoped; a MISS in both
# columns proves it was not discussed at all, which is what this lens is for.  Every miss
# below was additionally confirmed by reading the table.
#
# NO ARMS AND NO ISOLATION LEDGER: this leg is a byte-level scan, it varies nothing, and it is
# reported as a scan and not as an experiment.

import os
import re

BASE = "/Users/bgm/MB Work/where-atoms-come-from"
TABLE = os.path.join(BASE, "LANE_W10_D_SCOPE", "W10D_SCOPE_TABLE_V001.md")
LANE = os.path.join(BASE, "LANE_W10_D_SCOPE")

table_txt = open(TABLE, encoding="utf-8", errors="replace").read()
lane_txt = ""
for fn in sorted(os.listdir(LANE)):
    p = os.path.join(LANE, fn)
    if os.path.isfile(p):
        lane_txt += open(p, encoding="utf-8", errors="replace").read()

# (label, register/source pointer, tokens -- ANY hit counts)
CLAIMS = [
    ("S3/S2's TWENTY SEALED CORRECTIONS by name",
     "brief's CARRY list; S3 audit :785-800; S2 audit sec8", ["COR-A", "COR-B", "COR-F", "COR-K", "COR-"]),
    ("COR-F: loop transport is NOT diagonal in general",
     "REGISTER W-06 'the thing that decided the spine'", ["edge by edge", "edge-by-edge", "edge transport",
                                                          "edge-transport", "T^3", "T³"]),
    ("the transport CONVENTION as an unledgered stipulation",
     "REGISTER W-06; W-04 'what K1 was' (ii)", ["whole-circuit", "edge tick", "edge-tick", "unledgered stipulation"]),
    ("W-02's trap figure: circuits span 3 DIMENSIONS at N=1 and at N=100",
     "REGISTER W-02; S3:264-271", ["3 dimensions", "three-dimensional", "3-dimensional", "span of {", "circuits grow"]),
    ("the record slot is NON-ABELIAN BY NECESSITY (pure states of an abelian algebra)",
     "REGISTER W-02", ["non-abelian by necessity", "overlap only 0 or 1", "abelian algebra have overlap"]),
    ("P-9 resolved: ||omega_F^N - omega_C^N|| -> 2.000000000000",
     "REGISTER W-02", ["omega_F", "ω_F", "2.000000000000", "orthogonal reduced support"]),
    ("P-9's finite-stage exception: holds exactly on W-01's FIRING LOCUS",
     "S3 sec5.9, cited by REGISTER W-02", ["firing locus", "firing-locus"]),
    ("THE EXCEPTIONAL SET: lambda varies only there; closed form on the primitive locus",
     "S4 sec3.1, carried by REGISTER W-03 and the ERRATUM v W-02", ["exceptional set", "exceptional value",
                                                                    "primitive locus", "one-variable Mahler"]),
    ("the exceptional-value split 527/314/213 (W-03's correction to S4's 638/380/258)",
     "REGISTER W-03 'FURTHER CORRECTIONS SURVIVING ATTACK'", ["527", "314/213", "638"]),
    ("(pi,pi) is a THIRD STRICT SADDLE, not a local minimum; 1-2+2 = +1 vs chi(T^2)=0",
     "REGISTER W-03 'FURTHER CORRECTIONS SURVIVING ATTACK'", ["saddle", "strict saddle", "chi(T^2)"]),
    ("LAWTON's theorem (1983) missing from S4's IMPORT AUDIT",
     "REGISTER W-03; W-05 rediscovery ledger", ["Lawton"]),
    ("the ERRATUM v W-02's substantive content: subtorus average != torus average",
     "REGISTER, ERRATUM AGAINST W-02", ["subtorus", "SUBTORUS", "-0.767014993", "0.767014"]),
    ("W-08's FLOOR half: the floor's contribution is O(1) and K-INDEPENDENT; floor and rate do not race",
     "REGISTER W-08 headline", ["do not race", "independent coordinates", "K-independent", "decay budget"]),
    ("W-08's onset non-uniformity 2,10,47,216,1000,4642 and K_0 ~ t^{-2/3}",
     "REGISTER W-08", ["4642", "K_0", "onset"]),
    ("B5: b1 = 0, gamma_C cannot be designated, THE FORMATION DATUM DOES NOT EXIST",
     "S4 sec4.1, the ten-carrier table lane D quotes", ["B5", "double sphere", "double-covered"]),
    ("B2: the deliberate exception -- both loops bound, NO FLAT HOLONOMY AT ALL",
     "S4 sec4.1, same table", ["B2", "both filled", "both triangles filled"]),
    ("B0a: ring torus, loops DISJOINT -- the three-class arm on the same complex as B0b",
     "S4:575, same table", ["B0a"]),
    ("the two designated loops' LENGTHS, and the circuit clock",
     "REGISTER W-01 'circuit count is carrier-supplied discrete time'; S3 sec3.5",
     ["loop length", "loop lengths", "circuit clock", "edge time", "|gamma_F|"]),
    ("W-03 sec2's foundational sentence 'every vertex phase of s cancels'",
     "REGISTER W-03", ["vertex phase", "phases cancel", "every vertex phase"]),
    ("W-06's Bell R_N transplant, cross-expectation modulus exactly 1 at every N=1..9",
     "REGISTER W-06", ["R_N", "Bell", "4.11e-01", "8.42e-04"]),
    ("W-05's N1 'inherits the entropy theory of algebraic Z^d-actions wholesale'",
     "REGISTER W-05 N1", ["Z^d", "algebraic Z", "entropy theory", "Lind"]),
    ("the CHARGE run: THEOREM S4-1 FAILS at charge != 1; the taxonomy is a charge-1 statement",
     "REGISTER W-03", ["charge", "S4-1"]),
    # --- controls: claims I EXPECT to be present, to show the scan is not rigged ---
    ("[CONTROL] W-01's convex-hull criterion", "REGISTER W-01", ["convex-hull", "convex hull", "hull"]),
    ("[CONTROL] the multiset / invisibility theorem", "REGISTER W-03 / N2", ["multiset"]),
    ("[CONTROL] V = 4k+1 wedge-growth rebuild route", "REGISTER W-06", ["4k+1", "wedge"]),
    ("[CONTROL] the schedule adversary's constants 0.606/0.615/0.588/0.601", "REGISTER W-08",
     ["0.606", "0.615", "0.588"]),
]


# CASE-SENSITIVE for carrier labels, case-insensitive otherwise.
# MY OWN CONFOUND, RECORDED NOT PATCHED: the first version matched case-insensitively
# throughout, and the row "B2: the deliberate exception" therefore printed "yes: B2" against
# the scope table -- a FALSE POSITIVE off the Betti number "b2" (the table writes "b2=1",
# "b2=2" in its topology columns).  `grep -n "B2" W10D_SCOPE_TABLE_V001.md` returns NOTHING.
# A SECOND false positive in the same run: N1's "algebraic Z^d / Lind-Schmidt-Ward" row matched
# on "Lind" inside the table's own phrase "curvature-BLIND".  BOTH false positives ran AGAINST
# my own finding: the first run reported 17 omissions from the table, the corrected run reports
# 19.  Both counts are stated and neither is patched away.
CASE_SENSITIVE = {"B2", "B5", "B0a", "R_N", "K_0", "COR-A", "COR-B", "COR-F", "COR-K", "COR-",
                  "Lawton", "Lind", "T^3"}


def hit(txt, toks):
    for t in toks:
        if t in CASE_SENSITIVE:
            if t in txt:
                return t
        elif t.lower() in txt.lower():
            return t
    return None


print("=" * 112)
print("== 4A  TOKEN-COVERAGE SCAN OF LANE D's SCOPE TABLE AGAINST THE REGISTER'S OWN CLAIMS ==")
print("=" * 112)
print(f"  scope table  : {TABLE}   ({len(table_txt)} bytes)")
print(f"  whole lane   : {LANE}   ({len(lane_txt)} bytes over "
      f"{len([f for f in os.listdir(LANE) if os.path.isfile(os.path.join(LANE,f))])} files)")
print(f"\n  {'claim':66s} {'in TABLE':>12s} {'in LANE dir':>14s}")
miss_table, miss_lane = [], []
for lab, src, toks in CLAIMS:
    ht, hl = hit(table_txt, toks), hit(lane_txt, toks)
    print(f"  {lab[:66]:66s} {(('yes: '+ht) if ht else 'NO'):>12s} "
          f"{(('yes: '+hl) if hl else 'NO'):>14s}")
    if not lab.startswith("[CONTROL]"):
        if not ht:
            miss_table.append(lab)
        if not hl:
            miss_lane.append(lab)

print(f"\n  CONTROLS: every [CONTROL] row must be a hit in both columns, or the scan is broken.")
ctl = [(lab, hit(table_txt, toks), hit(lane_txt, toks)) for lab, _, toks in CLAIMS
       if lab.startswith("[CONTROL]")]
print(f"  controls hitting the table: {sum(1 for _,a,_ in ctl if a)} of {len(ctl)}"
      f"   controls hitting the lane: {sum(1 for _,_,b in ctl if b)} of {len(ctl)}")

print(f"\n  ABSENT FROM THE SCOPE TABLE ({len(miss_table)} of {len(CLAIMS)-len(ctl)} non-control claims):")
for lab in miss_table:
    print(f"    - {lab}")
print(f"\n  ABSENT FROM THE ENTIRE LANE DIRECTORY, CODE AND OUTPUT INCLUDED ({len(miss_lane)}):")
for lab in miss_lane:
    print(f"    - {lab}")

print("\n" + "=" * 112)
print("== 4B  THE SEALED-CORRECTION COUNT, EXACTLY ==")
print("=" * 112)
for nm, path in (("S3 audit COR-A..COR-L", os.path.join(BASE, "S3_THE_CROSSING_AUDIT_V001.md")),
                 ("S2 audit COR-A..COR-H", os.path.join(BASE, "S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md"))):
    t = open(path, encoding="utf-8", errors="replace").read()
    labels = sorted(set(re.findall(r"COR-[A-L]\b", t)))
    print(f"  {nm}: {len(labels)} distinct labels in the sealed artifact -> {labels}")
allc = sorted(set(re.findall(r"COR-[A-L]\b", table_txt)))
print(f"  distinct COR- labels named in LANE D's SCOPE TABLE: {len(allc)} -> {allc}")
allc2 = sorted(set(re.findall(r"COR-[A-L]\b", lane_txt)))
print(f"  distinct COR- labels named ANYWHERE in LANE D's directory: {len(allc2)} -> {allc2}")
reg = open(os.path.join(BASE, "REGISTER_V001.md"), encoding="utf-8", errors="replace").read()
regc = sorted(set(re.findall(r"COR-[A-L]\b", reg)))
print(f"  distinct COR- labels named in REGISTER_V001.md: {len(regc)} -> {regc}")
print("\n  The brief for this round carries four of them forward BY NAME as things not to")
print("  rediscover -- COR-B, COR-E, COR-F, COR-K.  W-06's register row is built on COR-F and")
print("  says in terms that the chain's failure was UNDER-READING it.  Twenty sealed")
print("  corrections; zero named in the scope table.")
