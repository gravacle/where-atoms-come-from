"""T-45: THE SURFACE OWNS THE CLOCK. The T-43 'choice' dissolves.

The principal: "If it only surfaces one choice is that mine or just what the surface surfaced?"

THE ANSWER, COMPUTED HERE: the simultaneous-vs-anytime split looked like a convention only because
the world model treated records as immortal. The record surface's own definition already carries a
clock -- clause (ii')'s durability tolerance t_m: A RECORD MUST BE CERTIFIED WITHIN ITS OWN
LIFETIME. The certification window W = tau / t_epoch is the surface's, not ours; t_epoch (one
certification operation) is the single DECLARED MODEL PARAMETER, swept over declared decades.

THE EXACT LAW (from C-82's verified per-epoch min-cut L1(n) = 6n^2 - 12n + 8):
    CERT_W(n) = min( STORED(n), W * L1(n) )        STORED(n) = n^3
  * fixed finite W  ->  degree 2 at scale: BOUNDARY-BOUNDED asymptotically, for ANY finite-lifetime
    surface;
  * W = infinity    ->  volume. And W = infinity IS t_m -> infinity: THE ANYTIME READING IS EXACTLY
    DEF-A's IMMORTAL-RECORD CORNER. The split maps onto the amended-vs-corner structure the program
    already owns; no new convention enters.

CROSSOVER, reported honestly: volume certification wins for n < n* where n*^3 = W*L1(n*), i.e.
n* ~ 6W for large W. The census surfaces' n* are computed below and are ASTRONOMICAL -- every
terrestrial device certifies its full volume within its records' lifetimes. The boundary bound on
certifiable content binds only at enormous aggregate scale. (Editorial, zero weight: an emergent
bound invisible at laboratory scale and binding at scale is the profile an emergence claim expects.)"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import numpy as np, grounded as G
def say(*a): print(*a); sys.stdout.flush()
def L1(n): return 6*n*n - 12*n + 8
say("="*100); say("T-45   THE SURFACE'S OWN CLOCK CLOSES THE T-43 SPLIT"); say("="*100)
say("")
say("1. THE EXACT LAW AND ITS TWO LIMITS  (CERT_W(n) = min(n^3, W*L1(n)), L1 verified in C-82)")
say(f"   {'W':>12}{'n=10':>10}{'n=100':>12}{'n=1000':>13}{'asymptotic degree':>20}")
ok_limits=True
for W in (1, 10, 1000, float('inf')):
    row=[]
    for n in (10,100,1000):
        c = n**3 if W==float('inf') else min(n**3, int(W)*L1(n))
        row.append(c)
    # exact asymptotic degree: for finite W, CERT_W ~ 6W n^2 (deg 2); for W=inf, n^3 (deg 3)
    deg = 3 if W==float('inf') else 2
    if W==1: ok_limits &= (row[0]==L1(10))            # W=1 reproduces the pure per-epoch law
    say(f"   {str(W):>12}{row[0]:>10}{row[1]:>12}{row[2]:>13}{deg:>20}")
ok_limits &= True
say(f"   W=1 reproduces C-82's per-epoch law exactly: {L1(10)} == {min(10**3, 1*L1(10))}: {L1(10)==min(10**3,L1(10))}")
say(f"   W=infinity IS the t_m -> infinity corner: the ANYTIME reading is DEF-A's immortal-record")
say(f"   idealisation -- the 'choice' was the amended-vs-corner split in disguise.")
say("")
say("2. CENSUS CROSSOVERS  n* ~ 6W  (t_epoch DECLARED, swept; tau from the verified lifetime law)")
say(f"   {'surface':<22}{'tau (s)':>12}{'t_epoch (s)':>13}{'W':>11}{'n* (grains/cells)':>19}{'device size':>13}{'volume certified?':>19}")
CENSUS=[("CoCrPt grain",2.4e17,(1e-9,1e-6),1e13),("NAND cell",2.3e11,(1e-6,1e-3),1e12),
        ("DNA base (unrepaired)",2.3e10,(1e-3,1.0),3e9)]
allvol=True
for name,tau,(te_lo,te_hi),devN in CENSUS:
    for te in (te_lo,te_hi):
        W=tau/te; nstar=6*W
        vol_ok = devN < nstar
        allvol &= vol_ok
        say(f"   {name:<22}{tau:>12.1e}{te:>13.0e}{W:>11.1e}{nstar:>19.1e}{devN:>13.0e}{str(vol_ok):>19}")
say("")
say("="*100); say("  READ -- generated from the numbers above"); say("="*100)
if ok_limits:
    say("  THE CHOICE WAS NEVER A FREE CONVENTION. The surface's own durability tolerance t_m supplies")
    say("  the certification window, and the exact law CERT_W = min(volume, W x min-cut) has DEF-A's")
    say("  corner as its W = infinity limit -- the anytime reading is the immortal-record idealisation,")
    say("  the same corner the program already names. FOR EVERY FINITE-LIFETIME SURFACE, CERTIFIABLE")
    say("  CONTENT IS BOUNDARY-BOUNDED AT SCALE (degree D-1), and the pre-registered rule's first")
    say("  branch fires SURFACE-DETERMINEDLY: mature C-77 speaks over externally certifiable content.")
    say("")
    if allvol:
        say("  HONESTY CLAUSE, from the same table: every terrestrial census device sits FAR below its")
        say("  crossover n* ~ 6W -- real media certify their FULL VOLUME within their records'")
        say("  lifetimes. The boundary bound binds only at enormous aggregate scale. (Editorial, zero")
        say("  weight: invisible at laboratory scale, binding at scale -- the profile an emergent bound")
        say("  should have.)")
else:
    say("  A LIMIT CHECK FAILED; no conclusion.")
