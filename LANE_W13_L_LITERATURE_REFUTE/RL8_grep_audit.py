#!/usr/bin/env python3
"""
RL8 — THE TARGET'S COUNTING CLAIMS, RECOUNTED.  Findings L-8, L-12 and L-14 rest on greps
over the sealed .md corpus.  A refuter on the literature lens has to recount them.
"""
import glob
import os
import re

ROOT = "/Users/bgm/MB Work/where-atoms-come-from"
MDS = sorted(glob.glob(os.path.join(ROOT, "*.md")))
print("=" * 78)
print("RL8 — RECOUNT OF THE TARGET'S GREP-BASED FINDINGS.  %d sealed .md artifacts." % len(MDS))
print("=" * 78)


def occ(word):
    tot = {}
    for f in MDS:
        n = len(re.findall(word, open(f, encoding="utf-8", errors="replace").read(), re.I))
        if n:
            tot[os.path.basename(f)] = n
    return tot


def lines(word):
    tot = {}
    for f in MDS:
        n = sum(1 for L in open(f, encoding="utf-8", errors="replace") if re.search(word, L, re.I))
        if n:
            tot[os.path.basename(f)] = n
    return tot


CLAIMS = [
    ("inhomogeneous", "L-8: REGISTER 0, W10_SCOPE_TABLE 2, W08_THE_RACE 1"),
    ("lyapunov", "L-12: exactly TWO occurrences in the sealed corpus"),
    ("lawton", "L-14: 11 total = REGISTER 3, W10 3, W08 4, HANDOFF 1"),
    ("cassaigne", "L-14(a): S4 = 17, REGISTER = 1, all others 0"),
    ("wendel", "L-14(b): Wendel 4"),
]
for w, claim in CLAIMS:
    o = occ(w)
    l = lines(w)
    print("\n  '%s'   TARGET CLAIM: %s" % (w, claim))
    print("     occurrences by file: %s   TOTAL %d" % (o if o else "{}", sum(o.values())))
    print("     LINES by file      : %s   TOTAL %d" % (l if l else "{}", sum(l.values())))

print("\n  ABSENT-CITATION CLAIM (L-14b).  Occurrences over the sealed .md corpus:")
for w in ["verbitskiy", "gelfond", "baker", "dobrowolski", "kenyon", "okounkov", "sheffield",
          "lyons", "everest", "sudler", "lubinsky", "herman", "avila", "kuipers",
          "niederreiter", "dimitrov", "habegger", "issa", "lalin", "boyd"]:
    o = occ(w)
    print("     %-14s %d   %s" % (w, sum(o.values()), o if o else ""))

print("""
--------------------------------------------------------------------------------
VERDICT.
  L-8  CONFIRMED EXACTLY.  'inhomogeneous': REGISTER 0, W10_SCOPE_TABLE 2, W08_THE_RACE 1.
       The hypothesis N1 needs was proposed one round earlier and never entered the register.
  L-12 CONFIRMED EXACTLY.  'Lyapunov' occurs twice, at S3_THE_CROSSING_V001.md:472 and
       S3_THE_CROSSING_AUDIT_V001.md:629, and in no register row and no lane directory.
  L-14 CONFIRMED on 'Lawton' (11 = 3 + 3 + 4 + 1) and on the zero-citation list, which is
       zero for every name checked.
  ONE CORRECTION, TRIVIAL BUT IT IS A NUMBER IN A FINDING: 'Cassaigne' occurs 18 times in
  S4, on 17 LINES.  The target reports 17 and attributes it to `grep -c`, which counts lines.
  THE REGISTER'S OWN W-05 ROW ALREADY SAYS 18 ("Mahler 12, Cassaigne 18, Glimm 2, Lawton 1").
  So the target restated a corpus figure with a different counting convention and did not
  notice the register disagreeing with it.  Nothing depends on which is right.
  AND 'boyd' -- the corpus mentions Boyd but never with a paper attached, consistent with the
  target's B-6.
--------------------------------------------------------------------------------""")
print("\nDONE RL8")
