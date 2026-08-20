"""V4 -- DOES THE INSTRUMENT EVER CERTIFY A KNOWN-BOUNDED SERIES AS EXTENSIVE?

S6 prints, under its own D-17 sweep, a KNOWN-SATURATING series (chi of ONE record read off f
bath qubits, hard-capped at 1 bit because it is one classical bit) with doubling ratios
1.802 and 1.673 at two of the six settings.  The engine's OWN stated extensivity band is
Q(2N)/Q(N) in [1.6, 2.6].  So on the doubling test alone the known-saturating control lands
INSIDE the extensive band at two settings, while S6's prose says it "does not, at any setting
tested".  This script runs the lane's own classify() on those settings and reports what the
composite EXT? gate actually returns.

It also asks the same question of the exact 1-bit cap: chi(one record : anything) <= 1 bit, so
NO setting of this series can be extensive at any f -- an exact statement the doubling test
cannot see.
"""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from fss_lib import classify
from chi_lib import vN, TIMES
from record_model import Environment

OUT = []
def say(s=""):
    print(s); OUT.append(s)

def frag_chi(f, nq=8, lam=0.8, beta=0.5, times=TIMES):
    """chi(one record : first f bath qubits), record couples to ALL nq sites (S6's construction
       reproduced independently: one +-1 bit driving the whole bath, read on a fragment)."""
    en = tuple((1.0, 1.4, 0.7)[j % 3] for j in range(nq))
    env = Environment(nq=nq, energies=en, beta=beta)
    rB = env.thermal()
    vals = []
    for t in times:
        half = {}
        for s in (+1, -1):
            H = env.HB + lam*s*env.probe
            w, U = np.linalg.eigh(H); ph = np.exp(-1j*w*t)
            Uc = U.conj().T @ rB @ U
            r = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
            half[s] = env._fragment(r, tuple(range(f)))
        vals.append(max(vN(0.5*(half[1]+half[-1])) - 0.5*(vN(half[1])+vN(half[-1])), 0.0))
    v = np.array(vals)
    return v.mean(), v.std(ddof=1)/np.sqrt(len(v))

say("="*112)
say("V4  DOES THE ENGINE'S EXTENSIVITY GATE EVER PASS A SERIES THAT IS EXACTLY CAPPED AT 1 BIT?")
say("="*112)
say()
say("   The series: chi(ONE record : fragment of f bath qubits), f = 1..8.  It is ONE classical")
say("   bit, so chi <= 1 bit at EVERY f, by an exact argument.  It cannot be extensive at any f.")
say()
say("     lam   beta |  chi(1)   chi(4)   chi(8)  | dbl(8/4) | engine category            best dAICc  expo       EXT?")
bad = []
for lam, beta in [(0.4, 0.5), (0.4, 2.0), (0.8, 0.5), (0.8, 2.0), (1.6, 0.5), (1.6, 2.0)]:
    F = np.arange(1, 9, dtype=float)
    m = [frag_chi(int(f), lam=lam, beta=beta) for f in F]
    Q = np.array([x[0] for x in m]); sg = np.array([max(x[1], 1e-6) for x in m])
    r = classify(F, Q, sg, "frag")
    dbl = Q[7]/Q[3]
    say("    %5.2f %5.2f | %7.4f %7.4f %7.4f | %8.3f | %-26s %-5s %5.1f %6.2f %8s"
        % (lam, beta, Q[0], Q[3], Q[7], dbl, r["category"], r["best"], min(r["dAICc"], 999),
           r["expo"], "YES" if r["extensive"] else "no"))
    if 1.6 <= dbl <= 2.6: bad.append((lam, beta, dbl, r["extensive"]))
say()
say("   Settings where the DOUBLING TEST ALONE put this exactly-capped series inside the engine's")
say("   own extensive band [1.6, 2.6]: %d of 6." % len(bad))
for lam, beta, d, e in bad:
    say("     lam=%.2f beta=%.2f  dbl = %.3f  -- composite EXT? gate returned %s"
        % (lam, beta, d, "YES (FALSE POSITIVE)" if e else "no (the other two gate conditions caught it)"))
say()
say("   READ: S6's prose sentence 'the saturating column does not [give 2.000], at any setting")
say("   tested' is not what the numbers say against the engine's own [1.6,2.6] band.  The")
say("   COMPOSITE gate (LIN best by dAICc>=4 AND doubling in band AND |expo-1|<0.25) is what")
say("   rejects these rows.  The doubling ratio on its own is not a sufficient instrument, and")
say("   the lane leans on it verbally in S6.  The EXACT 1-bit cap is what actually settles it.")
say("="*112)
open(LANE+"/VERIFY/v4_instrument.txt", "w").write("\n".join(OUT)+"\n")
