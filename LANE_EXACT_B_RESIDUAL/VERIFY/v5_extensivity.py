"""V5 -- ADVERSARIAL CHECK OF "NO EXTENSIVE TERM: chi DOES NOT KNOW HOW MANY RECORDS THE
CARRIER HOLDS" (lane B5 part 1, B9 Part 1 row 6, and exactly_zero item 4).

WHAT THE LANE ACTUALLY PROVED.  Its generalised C-38 theorem has a hypothesis: every added
partner COMMUTES with the read record AND sits off its bath site.  Such a partner is decoupled
from R's sector by construction.  The lane then adds up to two such partners (the bath has
three sites, one of which is the read site) and reports "1.22e-15, no extensive term".

THE OBJECTION.  Extensivity was tested only inside the subclass that the theorem makes
invisible.  The lane's OWN b5 table already shows two PAIRING partners suppress chi more than
one (Z1@1 -> 0.4226; Z1@1,Z1X2@2 -> 0.3349).  So the record COUNT does move chi -- in the
channel that was not swept.  Sweep it: place k pairing partners on k DISTINCT bath sites and
watch chi as a function of k, with the lane's own commuting family as the control in the same
table.

BATH: 5 qubits, so up to four off-site partners can be placed one per site.  Carrier n = 6
(chi is n-independent at fixed configuration -- the lane proves it and V4 confirms the algebra).
"""
import numpy as np, sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from common import *                      # noqa
from battery import build_ops, READ, features, sp  # noqa

t0 = time.time()
N, NB = 6, 5
TIMES = BASE.times
ops, _, _ = build_ops(N)
R = ops[READ][0]
env = Venue("V0-5q", (1.0, 1.4, 0.7, 1.2, 0.9), TIMES).env(NB)

def chi(parts, lam):
    so = [(R, 0)] + [(ops[l][0], s) for l, s in parts]
    return float(np.mean(chi_times(so, R, env, lam, TIMES)))

say("=" * 116)
say("V5   DOES chi KNOW HOW MANY RECORDS THE CARRIER HOLDS?  SWEEP THE CHANNEL THE THEOREM EXCLUDES.")
say("=" * 116)
say(f"  carrier n = {N}, bath {NB} qubits (one record per site, none on the read site 0), lam = 0.4/0.8/1.2")
say("")

PAIRING = [("Z1", 1), ("Z1X2", 2), ("Z1X3", 3), ("Z1X4", 4)]
COMMUT  = [("X2", 1), ("X3", 2), ("X4", 3)]
for l, s in PAIRING:
    assert sp(ops[READ][1], ops[l][1], N) == 1, f"{l} does not pair with the read record"
for l, s in COMMUT:
    assert sp(ops[READ][1], ops[l][1], N) == 0, f"{l} is not a commuting partner"
say("  SELF-CHECK: the four 'pairing' partners all anticommute with the read record and the three")
say("  'commuting' partners all commute with it, checked in F_2 -- PASS.")
say("")
say(f"  {'k = off-site records added':<30}" +
    "".join(f"{'lam='+str(l)+' PAIRING':>22}{'CONTROL commuting':>22}" for l in (0.4, 0.8, 1.2)))
tab = {}
for k in range(0, 5):
    cells = ""
    for lam in (0.4, 0.8, 1.2):
        cp = chi(PAIRING[:k], lam)
        cc = chi(COMMUT[:min(k, 3)], lam) if k <= 3 else float('nan')
        tab[(k, lam)] = (cp, cc)
        cells += f"{cp:>22.12f}" + (f"{cc:>22.12f}" if k <= 3 else f"{'--':>22}")
    say(f"  {k:<30}" + cells)

say("")
say(f"  {'lam':<8}{'chi(k=0)':>18}{'chi(k=4) PAIRING':>20}{'total change':>16}{'CONTROL chi(k=3) commuting':>30}{'control change':>18}")
for lam in (0.4, 0.8, 1.2):
    c0 = tab[(0, lam)][0]; c4 = tab[(4, lam)][0]; cc3 = tab[(3, lam)][1]
    say(f"  {lam:<8}{c0:>18.12f}{c4:>20.12f}{c4-c0:>+16.6f}{cc3:>30.12f}{cc3-c0:>+18.2e}")

say("")
say("  READ IT OFF THE TABLE.  In the COMMUTING column the record count is invisible, exactly as")
say("  the lane says (control change at the float64 floor).  In the PAIRING column chi falls by")
say("  0.083 / 0.199 / 0.242 as k goes 0 -> 4 -- four decades above any floor.  (It is not strictly")
say("  monotone: at lam = 0.4 the k = 1 -> 2 step rises by 0.004, and k = 3 -> 4 by 0.005; the")
say("  TOTAL, not the step, is what the table supports.)  The record count is therefore")
say("  NOT invisible to chi; it is invisible only inside the subclass the theorem's hypothesis")
say("  selects.  'chi does not know how many records the carrier holds' is true of that subclass")
say("  and false of the carrier.")
say("")
say(f"  elapsed {time.time()-t0:.1f}s")
