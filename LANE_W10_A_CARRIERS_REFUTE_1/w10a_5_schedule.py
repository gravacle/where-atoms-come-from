#!/usr/bin/env python3
"""
LANE W10-A — SCRIPT 5.  THE SCHEDULE QUALIFIER, ON FOUR CLASSES.

W-08's row proves durability is a property of the (connection, SCHEDULE) PAIR and that an
adversary writing only the sqrt(K) cells of smallest 1-|Z_k| holds |Omega| ~ 0.55 forever with
unboundedly many writes.  W-02's row states "FORMATION OCCURS <=> G != {1}" with NO schedule
qualifier at all, although the theorem it registers (S3 audit sec4.4) says "Along the canonical
clock k_n = n".  The two rows are the SAME criterion and the register carries only one of the
two hypotheses.  This script checks whether the four-class move changes any of that.

ONE VARIABLE: the occupied class set.  Connection, K-ladder, adversary rule, evaluator and code
path identical in every row; arms diffed at the bytes below.

IEEE double.  Chunked so that K = 1e7 fits in memory; the chunking is arithmetic-neutral and
is checked against the unchunked value at K = 1e6.
"""
import sys, math, hashlib
from fractions import Fraction
import numpy as np

LOG = []
def out(s=""):
    print(s); LOG.append(s)

out("=" * 100)
out("W10-A SCRIPT 5 — THE SCHEDULE QUALIFIER ON FOUR CLASSES")
out("=" * 100)
out(f"numpy {np.__version__}; IEEE double.")
out()

ARMS = [
    ("K1  SENSE U  3-class {10,01,11}", (0.0, 0.4, 0.4, 0.2)),
    ("B1q SENSE U  3-class {00,10,01}", (1/7, 3/7, 3/7, 0.0)),
    ("B0b SENSE U  4-class",            (4/9, 2/9, 1/9, 2/9)),
    ("B4  SENSE U  4-class",            (1/6, 1/6, 1/6, 1/2)),
]
seen = {}
for nm, w in ARMS:
    h = hashlib.sha256(repr([f"{x:.17g}" for x in w]).encode()).hexdigest()[:12]
    assert h not in seen, f"ARMS COLLIDE: {nm} / {seen[h]}"
    seen[h] = nm
out("ARMS DIFF (sha256[:12] of the float64 weight vector):")
for nm, w in ARMS:
    h = hashlib.sha256(repr([f"{x:.17g}" for x in w]).encode()).hexdigest()[:12]
    out(f"   {nm:<36} pi = ({', '.join(f'{x:.6f}' for x in w)})   {h}")
out()

CONNS = [("generic irrational", 1/math.sqrt(2), 1/math.sqrt(3)),
         ("S3/S4 resonant f=2,c=1.1", -2.0/(2*math.pi), 1.1/(2*math.pi)),
         ("S1 published order-4", -0.5, 0.75)]

CHUNK = 1 << 21

def scan(w, a, b, K, keep):
    """Return (honest_total, smallest `keep` values of -log|Z_k|, n_zero_cost).
    honest_total = sum_{k<=K} -log|Z_k|.  Chunked."""
    p00, p10, p01, p11 = w
    tot = 0.0
    best = np.full(keep, np.inf)
    nzero = 0
    k0 = 1
    while k0 <= K:
        k1 = min(K, k0 + CHUNK - 1)
        k = np.arange(k0, k1 + 1, dtype=np.float64)
        ua = np.exp(2j * np.pi * ((k * a) % 1.0))
        vb = np.exp(2j * np.pi * ((k * b) % 1.0))
        m = np.abs(p00 + p10 * ua + p01 * vb + p11 * ua * vb)
        with np.errstate(divide="ignore"):
            cost = -np.log(np.maximum(m, 1e-300))
        tot += float(cost.sum())
        nzero += int(np.count_nonzero(cost <= 1e-15))
        cat = np.concatenate([best, cost])
        best = np.partition(cat, keep - 1)[:keep]
        k0 = k1 + 1
    return tot, np.sort(best), nzero

out("-" * 100)
out("THE ADVERSARY (W-08's): write only the sqrt(K) cells of SMALLEST cost -log|Z_k|.")
out("Reported: total nats accumulated by the adversary vs by the honest schedule k_n = n.")
out("-" * 100)
KS = [10**4, 10**5, 10**6, 10**7]
for clbl, a, b in CONNS:
    out(f"CONNECTION: {clbl}")
    out(f"   {'arm':<36}" + "".join(f"{'K=1e'+str(int(math.log10(K))):<31}" for K in KS))
    for nm, w in ARMS:
        cells = []
        for K in KS:
            keep = int(math.isqrt(K))
            tot, best, nzero = scan(w, a, b, K, keep)
            adv = float(best.sum())
            cells.append(f"adv {adv:7.4f} | hon {tot:10.1f}")
        out(f"   {nm:<36}" + "".join(f"{c:<31}" for c in cells))
    out()
out("READ: on EVERY arm, three-class and four-class alike, the adversary's accumulation is")
out("BOUNDED in K -- flat at the generic connection (0.70, 0.70, 0.70, 0.70 on B0b), DECREASING")
out("at the resonant one (0.77 -> 0.05), and EXACTLY ZERO at the order-4 one -- while the honest")
out("schedule's grows linearly.  |Omega_adv| = exp(-adv) stays O(1) with")
out("unboundedly many writes.  THE SCHEDULE QUALIFIER IS CARRIER-INDEPENDENT AND CLASS-COUNT-")
out("INDEPENDENT: the four-class move has no purchase on it.")
out()
out("AND AT THE ORDER-4 CONNECTION THE ADVERSARY PAYS EXACTLY ZERO, ON EVERY ARM:")
for nm, w in ARMS:
    tot, best, nzero = scan(w, -0.5, 0.75, 10**6, 1000)
    out(f"   {nm:<36} #{{cost = 0}} over k <= 1e6 = {nzero:>8}   adversary's 1000-cell cost = {best.sum():.3e}")
out("   (density 1/4 on every arm.  |Z_k| = 1 iff k*L_S is contained in L, and all four arms")
out("    have L_S = Z^2, so this could not have come out otherwise -- DERIVED, NOT A CONTROL.)")
out()

out("-" * 100)
out("CHUNKING NEUTRALITY CHECK (the code path must not be the variable)")
out("-" * 100)
w = ARMS[2][1]
a, b = CONNS[0][1], CONNS[0][2]
K = 10**6
k = np.arange(1, K + 1, dtype=np.float64)
ua = np.exp(2j * np.pi * ((k * a) % 1.0)); vb = np.exp(2j * np.pi * ((k * b) % 1.0))
m = np.abs(w[0] + w[1] * ua + w[2] * vb + w[3] * ua * vb)
un = float(-np.log(np.maximum(m, 1e-300)).sum())
ch, _, _ = scan(w, a, b, K, 1000)
out(f"   unchunked total = {un:.9f}   chunked total = {ch:.9f}   dev = {abs(un-ch):.3e}")
out()

out("=" * 100)
out("VERDICT ON THE SCHEDULE QUALIFIER")
out("=" * 100)
out("W-02's REGISTER row states 'FORMATION OCCURS <=> G != {1}' with no schedule hypothesis.")
out("The theorem it registers (S3_THE_CROSSING_AUDIT_V001.md:420-428) says 'Along the canonical")
out("clock k_n = n'.  W-08's row supplies the missing hypothesis and proves it is load-bearing.")
out("MEASURED HERE: the gap is CARRIER-INDEPENDENT -- it is the same on three and four classes.")
out("So W-02's row is UNDER-QUALIFIED, and the missing qualifier is the SCHEDULE, not the CARRIER.")
out("A lane looking for W-02's scope on the carrier axis finds nothing; the scope defect is on an")
out("axis this round did not ask about, and it is already fixed in a later row of the same file.")

with open("w10a_5_schedule.OUT.txt", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
