#!/usr/bin/env python3
"""LANE W-10 B — LEG F.  CUSTODY: WHAT THE CORPUS ALREADY OWNED ON THIS QUESTION.

Every claim below is a file test, re-runnable.  Nothing here is a computation about
lambda; it is about where the four-class multiset test already lives on disk.
"""
import hashlib
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def sha(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def sealed(lane, fname):
    seals = os.path.join(ROOT, lane, 'SEALS.sha256')
    if not os.path.exists(seals):
        return None
    for line in open(seals):
        h, _, n = line.strip().partition('  ')
        if n == fname:
            return h == sha(os.path.join(ROOT, lane, fname))
    return False


print(__doc__)
print("=" * 96)
print("F.1  A FOUR-CLASS WEIGHT VECTOR HAS ALREADY BEEN RUN THROUGH THE MULTISET TEST,")
print("     IN A SEALED LANE, AND THE REGISTER'S OWN FIGURE IS THAT RUN")
print("=" * 96)
f = os.path.join(ROOT, 'LANE_G_GROUP_REFUTER', 'g0_validate.py')
src = open(f).read()
print("  LANE_G_GROUP_REFUTER/g0_validate.py  sealed and verifying: %s" %
      sealed('LANE_G_GROUP_REFUTER', 'g0_validate.py'))
print("  its OUT.txt sealed and verifying    : %s" % sealed('LANE_G_GROUP_REFUTER', 'g0_validate.OUT.txt'))
m = re.search(r"base = \[.*?\]", src)
c = re.search(r"allcorn = np\.array\(\[.*?\]\)", src, re.S)
print("  the weights it permutes  : %s" % (m.group(0) if m else '??'))
print("  the corners it puts them on: %s" % (' '.join((c.group(0) if c else '??').split())))
print("  -> all four classes occupied.  The line the register quotes at W-03 --")
print("     '24 of 24 permutations invariant, worst spread 2.4e-15' -- is a FOUR-CLASS figure.")
out = open(os.path.join(ROOT, 'LANE_G_GROUP_REFUTER', 'g0_validate.OUT.txt')).read()
for line in out.splitlines():
    if 'multiset' in line or 'pinch' in line:
        print("     OUT: %s" % line.strip())
print("""
  W-09's closing sentence -- 'The corpus has never run a four-class carrier through
  anything' -- is therefore true of the CARRIERS and false of the multiset test, which was
  run on a four-class WEIGHT VECTOR.  Since lambda depends on the carrier only through the
  pushforward, that distinction is immaterial for W-03's theorem and material only for the
  question LEG C asks (what makes the pushforward real), which does need the incidence.""")

print()
print("=" * 96)
print("F.2  AND THE D4-VERSUS-S4 QUESTION WAS ALREADY ASKED, ANSWERED, AND LEFT UNSEALED")
print("=" * 96)
for fn in ('maps_refuter3.py', 'perm_exact.py', 'perm_search.py', 'prec_check.py'):
    path = os.path.join(ROOT, 'LANE_R_MAPS_REFUTER', fn)
    outp = path.replace('.py', '.OUT.txt')
    print("  LANE_R_MAPS_REFUTER/%-16s exists %s   in SEALS.sha256 %-5s   has an OUT.txt %s"
          % (fn, os.path.exists(path), sealed('LANE_R_MAPS_REFUTER', fn), os.path.exists(outp)))
src3 = open(os.path.join(ROOT, 'LANE_R_MAPS_REFUTER', 'maps_refuter3.py')).read()
for key in ("I PREDICTED THIS WOULD FAIL", "3 orbits of 8", "all 24 agree", "K1+", "30 dps"):
    print("  contains %-28s : %s" % (repr(key), key in src3))
print("""
  So the exact question this lane was commissioned on -- 'the Newton polygon gives only
  D4, why do all 24 agree?' -- was posed inside the W-03 audit round, answered 'the
  multiset invariance is a real identity strictly stronger than the polygon symmetry',
  and NEVER REACHED THE REGISTER.  The file carries no output and no seal, so under the
  corpus's own pointer rule its numbers are flagged, not inherited.  This lane's LEG B
  supplies what that file did not: the mechanism (the pointwise identity in both Jensen
  pairings) and the hypothesis (at least three coefficients collinear in C), and therefore
  the boundary at which D4 is all that is left.""")
