#!/usr/bin/env python3
"""LENS 2 (SCOPE), ATTACK 5.  DOES THE LANE'S REPORT SAY WHAT ITS SEALED ARTIFACTS SAY,
AND DID IT FIND EVERYTHING THE CORPUS ALREADY OWNED ON THIS QUESTION?

Two file audits, both re-runnable, no floating-point anywhere.

  5.0  REPRODUCTION.  All five lane scripts re-run and diffed against their sealed
       outputs; SEALS.sha256 verified at bytes.
  5.1  REPORT-VS-ARTIFACT.  Every numeric figure the lane's findings quote, searched for
       in the sealed file the finding cites.
  5.2  CUSTODY THE LANE MISSED.  The lane's own LEG F does custody archaeology and finds
       two items.  There are two more, and one of them is SEALED, which makes it
       inheritable rather than flagged.
"""
import hashlib
import os
import re
import subprocess

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_W10_B_MULTISET"
REPO = "/Users/bgm/MB Work/where-atoms-come-from"


def hdr(t):
    print()
    print('=' * 96)
    print(t)
    print('=' * 96)


print(__doc__)

# =============================================================================== 5.0
hdr("5.0  REPRODUCTION — SEALS AND BYTE-FOR-BYTE RE-RUN")
seals = {}
for line in open(os.path.join(LANE, "SEALS.sha256")):
    h, f = line.split()
    seals[f] = h
okc = 0
for f, h in seals.items():
    d = hashlib.sha256(open(os.path.join(LANE, f), 'rb').read()).hexdigest()
    okc += (d == h)
print("  SEALS.sha256: %d of %d files verify at bytes" % (okc, len(seals)))
for s in ("b1_legA", "b2_legB_complex", "b3_pushforward", "b4_involution_labels", "b5_custody"):
    out = subprocess.run(["python3", s + ".py"], cwd=LANE, capture_output=True, text=True)
    sealed = open(os.path.join(LANE, s + ".OUT.txt")).read()
    print("  %-24s re-run exit %d, output identical to its sealed OUT.txt: %s"
          % (s + ".py", out.returncode, out.stdout == sealed))

# =============================================================================== 5.1
hdr("5.1  REPORT-VS-ARTIFACT — FIGURES QUOTED IN THE FINDINGS, LOOKED FOR ON THE PAGE")
FIG = [
    ("B-07", "b3_pushforward.OUT.txt", "8.0e-16",
     "scalar-branch residual on B0b"),
    ("B-07", "b3_pushforward.OUT.txt", "1.3e-15",
     "scalar-branch residual on B4"),
    ("B-07", "b3_pushforward.OUT.txt", "7.0e-01",
     "edge-branch residual on B0b"),
    ("B-07", "b3_pushforward.OUT.txt", "1.1e+00",
     "edge-branch residual on B4"),
    ("B-06", "b3_pushforward.OUT.txt", "0.69",
     "lower end of the quoted collinearity-defect range"),
    ("B-09", "b4_involution_labels.OUT.txt", "1.88e-13",
     "closed form vs direct |Z_k|^2 -- the finding writes 1.9e-13, a fair rounding"),
    ("B-01", "b1_legA.OUT.txt", "5.55112e-17", "E1 vs E3 on B0b"),
    ("B-03", "b2_legB_complex.OUT.txt", "0.00981", "smallest between-block gap"),
    ("B-11", "b4_involution_labels.OUT.txt", "0.323797", "order-4 spread on GEN"),
    ("B-11", "b4_involution_labels.OUT.txt", "0.003428", "resonant spread on GEN"),
    ("B-13", "b5_custody.OUT.txt", "9.285e-06", "LANE_G's own four-class spread"),
]
for fid, fn, fig, what in FIG:
    txt = open(os.path.join(LANE, fn)).read()
    print("  %-5s %-30s %-12s present in the cited file: %-5s   (%s)"
          % (fid, fn, fig, fig in txt, what))
print()
print("  Rows that fail are of two kinds.  B-09's is a ROUNDING (the page prints 1.88e-13,")
print("  the finding writes 1.9e-13) and is not a defect.  B-06's 0.69 and B-07's four")
print("  figures are not roundings of anything on the page.")
print()
print("  THE FOUR B-07 FIGURES ARE NOT ON THE PAGE THEY CITE.  What the sealed C.3 block")
print("  actually prints:")
for line in open(os.path.join(LANE, "b3_pushforward.OUT.txt")):
    if "residual over n" in line:
        print("     " + line.rstrip())
print("""  The lane's scripts re-run byte-identically (5.0), so this is not irreproducibility:
  it is a REPORT that quotes numbers its own sealed artifact does not contain.  That is
  the defect class the corpus already carries as COR-K against S3 ('published rows not
  reproducible from their parameters') and it is here again, in the report layer.""")
print()
print("  B-01's evidence says the direct ergodic method 'agrees to 1e-6'.  The sealed")
print("  LEG A block prints these E1-E4 deviations:")
for line in open(os.path.join(LANE, "b1_legA.OUT.txt")):
    if "E1-E4" in line:
        print("     " + line.rstrip())
print("  worst is 1.35e-05, an order of magnitude past the quoted 1e-6, and past S4's own")
print("  stated 3.0e-06 simulation error.  The claim is not wrong in substance; the number")
print("  quoted for it is.")

# =============================================================================== 5.2
hdr("5.2  CUSTODY THE LANE'S OWN LEG F MISSED — AND ONE OF THEM IS SEALED")
G2 = os.path.join(REPO, "LANE_G_GROUP_REFUTER", "g2_nonabelian.OUT.txt")
seal_g = os.path.join(REPO, "LANE_G_GROUP_REFUTER", "SEALS.sha256")
sealed_names = set()
if os.path.exists(seal_g):
    for line in open(seal_g):
        parts = line.split()
        if len(parts) == 2:
            sealed_names.add(parts[1])
print("  (i)  THE REGISTRAR'S ATTRIBUTION IS NOT THE REGISTRAR'S.  It is already in the")
print("       corpus, SEALED, from the W-03 audit round:")
print("       LANE_G_GROUP_REFUTER/g2_nonabelian.OUT.txt in that lane's SEALS.sha256: %s"
      % ("g2_nonabelian.OUT.txt" in sealed_names or "g2_nonabelian.py" in sealed_names))
for i, line in enumerate(open(G2)):
    if "NON-NEGATIVE" in line or "multiset theorem and pinch" in line:
        print("       :%d %s" % (i + 1, line.rstrip()))
print("""       So 'the multiset theorem is a theorem about NON-NEGATIVE MASS coefficients'
       is a SEALED corpus statement, not a registrar improvisation.  B-02's headline
       'THE REGISTRAR'S ATTRIBUTION IS WRONG' is correct on the merits and misattributes
       the attribution: it is inherited from a sealed lane, which under the corpus's own
       pointer rule makes it citable and makes the correction a correction OF THE CORPUS,
       not of one registrar sentence.  The lane's LEG F searched LANE_R_MAPS_REFUTER and
       LANE_G_GROUP_REFUTER's g0 and did not read g2.""")
print()
SIB = os.path.join(REPO, "LANE_W10_A_CARRIERS_REFUTE_2", "r2_4_multiset.OUT.txt")
print("  (ii) A SIBLING LANE IN THIS SAME W-10 ROUND REACHED THE SAME CHARACTERISATION")
print("       AND THE SAME COUNTEREXAMPLES, INDEPENDENTLY.  File: %s"
      % os.path.relpath(SIB, REPO))
keep = ("cos D1 = cos D2 = cos D3", "SUFFICIENCY", "NECESSITY", "REFUTE NECESSITY",
        "DEGENERATE: m = log(max", "registrar (real non-negative)")
for i, line in enumerate(open(SIB)):
    if any(k in line for k in keep):
        print("       :%d %s" % (i + 1, line.rstrip()))
print("""       'cos D1 = cos D2 = cos D3' is Lane B's '|phi_a| = |phi_b| = |phi_c|' written
       with cosines -- the same condition, since cos is injective on [0,pi].  Lane B's
       'at least three coefficients collinear' is its geometric form.  TWO INDEPENDENT
       LANES OF THIS ROUND CONVERGED ON ONE NAME, which is corroboration Lane B does not
       claim and could not have known about.
       BUT the sibling also states, and Lane B does not, that the condition is SUFFICIENT
       AND NOT NECESSARY, with its 18 exceptions all in the regime m = log(max|p|).  That
       is exactly the exceptional set ATTACK 3 measures at 91/400 in the generic regime,
       and it is what makes Lane B's 'S4-INVARIANCE <=>' and B-03's 'the invariance group
       is exactly D4' over-strong.""")
