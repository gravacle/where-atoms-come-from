# rule_main5.py -- W-19 RULING, part 5.  The floor UNDER THE GAUGE-INVARIANT SYSTEM ALGEBRA.
# Under alg{X_l} the purification jump at |F| = |E| disappears (a classical bit's mutual information
# with anything is bounded by its own entropy), so plateau points = d, not d - 1, and the floor for
# four points drops one step in d: d = 4, hence V >= 5, hence L >= ceil(3*5/2) = 8.
import numpy as np
from rule_verify import *
P("=" * 118)
P("W-19 RULING -- PART 5.  THE FLOOR MOVES AGAIN WHEN THE SYSTEM ALGEBRA IS MADE GAUGE-INVARIANT.")
P("=" * 118)
P("  %-14s %2s %3s %3s %14s   %-34s %-34s" % ("carrier","V","L","d","H_elec(S)","EXT/H  (pts)","CL/H  (pts)"))
rows = [("mg_chain4", mg_chain(4)), ("dbl_chain9", mg_chain(5)), ("tri_chain12", tri_chain12()),
        ("petersen", petersen()), ("heawood", heawood())]
for nm, (V, E) in rows:
    car = Carrier(nm, V, E); L = car.L
    frags, d = rule_A_fragments(V, E, 0)
    vec = car.lift(car.ground(0.50)[1])
    pX = np.abs(hadamard_all(vec, L)) ** 2
    HS = S_ax(vec, L, AX(L, [0]))
    pm = pX.reshape([2] * L).sum(axis=tuple(a for a in range(L) if a != L - 1))
    HX = float(-(pm[pm > 1e-15] * np.log2(pm[pm > 1e-15])).sum())
    ext = [channel_EXT(vec, L, 0, F) / HS for F in frags]
    cl = [channel_CL(pX, L, 0, F) / HX for F in frags]
    pts = lambda v: sum(1 for x in v if abs(x - 1) <= 0.10)
    P("  %-14s %2d %3d %3d %14.9f   %-28s (%d)  %-28s (%d)"
      % (nm, V, L, d, HX, " ".join("%.4f" % x for x in ext), pts(ext),
         " ".join("%.4f" % x for x in cl), pts(cl)))
P("")
P("  READ: EXT gives four points first at L = 12 (simple) / L = 9 (multigraph); alg{X_l} gives four")
P("  points first at L = 8, mg_chain(4), V = 5, d = 4, H_elec(S) = 0.738303126.  petersen moves from")
P("  MARGINAL (3) to EXHIBITED (4) at L = 15 on the same change.  A fourth number for one question.")
open("OUT_rule_main5.txt", "w").write("\n".join(LOG) + "\n")
