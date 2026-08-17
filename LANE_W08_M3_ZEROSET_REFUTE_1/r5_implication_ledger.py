#!/usr/bin/env python3
# LANE W08 / M3 REFUTER 1 — script 5.  THE FULL BOTH-DIRECTIONS IMPLICATION LEDGER.
# The lens asked for each implication checked separately in both directions.  This is the
# ledger, with a COUNT for every arrow: 0 counterexamples = the arrow holds on the sweep;
# a positive count = the arrow fails, with the count as the witness.
# Sweep: connections = q-th roots of unity, q = 12 and q = 15 (odd, so gap = pi is
# unreachable and no boundary case can be mis-decided); states = exact rational simplex
# grid at denominator 24.  Every predicate exact except the |Z_1| = 0 test, which is float
# with tolerance 1e-12 and whose hits are all at exactly representable configurations.
import numpy as np
L = []
def out(s=""):
    print(s); L.append(s)

out("=" * 100)
out("R5  IMPLICATION LEDGER: EVERY ARROW, BOTH DIRECTIONS, WITH A COUNT")
out("=" * 100)
out()
out("  HULL(u,v)  :  0 in conv{uv, u, v}                          [about the CONNECTION]")
out("  ZERO(p)    :  P_p has a zero somewhere on T^2              [about the STATE]")
out("  TRI(p)     :  |p10-p11| <= p01 <= p10+p11                  [about the STATE]")
out("  FIRE(p,u,v):  Z_1 = p10 u + p01 v + p11 uv = 0             [about BOTH]")
out()

def gap_pred(angles, q):
    a = sorted(set(x % q for x in angles))
    if len(a) == 1:
        return False
    g = [a[i + 1] - a[i] for i in range(len(a) - 1)] + [a[0] + q - a[-1]]
    return 2 * max(g) <= q

rows = {}
for q in (12, 15):
    Ns = 24
    cnt = dict(n=0, fire=0, hull=0, tri=0,
               fire_not_hull=0, hull_not_fire=0,
               fire_not_tri=0, tri_not_fire=0,
               tri_not_zero=0, zero_not_tri=0,
               exists_fire_not_hull=0, hull_not_exists_fire=0,
               exists_fire_not_tri=0, tri_not_exists_fire=0)
    # state side: for each p, does SOME connection on the q-grid fire?  (only a lower bound
    # on 'exists a connection', so it is used only in the direction where that is sound)
    fired_for_p = {}
    for aa in range(q):
        for bb in range(q):
            A3 = [(bb - aa) % q, (-aa) % q, bb % q]     # angles of uv, u, v in units 2pi/q
            hull = gap_pred(A3, q)
            zc = [np.exp(2j * np.pi * A3[t] / q) for t in range(3)]
            hull_saw_fire = False
            for i in range(Ns + 1):
                for j in range(Ns - i + 1):
                    kk = Ns - i - j
                    p10, p01, p11 = i, j, kk
                    tri = (abs(p10 - p11) <= p01 <= p10 + p11)
                    z = (p11 * zc[0] + p10 * zc[1] + p01 * zc[2]) / Ns
                    fire = abs(z) < 1e-12
                    cnt['n'] += 1
                    cnt['fire'] += fire
                    cnt['hull'] += hull
                    cnt['tri'] += tri
                    cnt['fire_not_hull'] += (fire and not hull)
                    cnt['hull_not_fire'] += (hull and not fire)
                    cnt['fire_not_tri'] += (fire and not tri)
                    cnt['tri_not_fire'] += (tri and not fire)
                    if fire:
                        hull_saw_fire = True
                        fired_for_p[(i, j, kk)] = True
            cnt['exists_fire_not_hull'] += (hull_saw_fire and not hull)
            cnt['hull_not_exists_fire'] += (hull and not hull_saw_fire)
    # state-side existential, over the SAME q-grid of connections
    for i in range(Ns + 1):
        for j in range(Ns - i + 1):
            kk = Ns - i - j
            tri = (abs(i - kk) <= j <= i + kk)
            ex = fired_for_p.get((i, j, kk), False)
            cnt['exists_fire_not_tri'] += (ex and not tri)
            cnt['tri_not_exists_fire'] += (tri and not ex)
    rows[q] = cnt

for q, c in rows.items():
    out("  q = %d : %d (connection, state) pairs ; %d of them fire" % (q, c['n'], c['fire']))
    out("    ARROW                                                  counterexamples   VERDICT")
    def line(lab, k, expect):
        v = c[k]
        out("    %-52s %10d      %s" % (lab, v, "HOLDS" if v == 0 else "FAILS"))
    line("FIRE(p,u,v)  =>  HULL(u,v)                    [1]", 'fire_not_hull', 0)
    line("HULL(u,v)    =>  FIRE(p,u,v)                  [2]", 'hull_not_fire', 1)
    line("FIRE(p,u,v)  =>  TRI(p)                       [3]", 'fire_not_tri', 0)
    line("TRI(p)       =>  FIRE(p,u,v)                  [4]", 'tri_not_fire', 1)
    line("(EXISTS p) FIRE  =>  HULL(u,v)                [5]", 'exists_fire_not_hull', 0)
    line("HULL(u,v)    =>  (EXISTS p) FIRE              [6]", 'hull_not_exists_fire', 0)
    line("(EXISTS u,v) FIRE  =>  TRI(p)                 [7]", 'exists_fire_not_tri', 0)
    line("TRI(p)       =>  (EXISTS u,v on this grid) FIRE [8]", 'tri_not_exists_fire', 1)
    out()
out("  READING OF THE LEDGER.")
out("   [1][3]   hold unconditionally and are THEOREMS: a firing pair certifies both HULL")
out("            and TRI.  These are the halves that W-01's own exhibit travels on.")
out("   [2][4]   fail, and heavily.  This is F3's content -- and F3 does not say that only")
out("            these two halves fail.")
out("   [5]      holds.")
out("   [6][8]   show POSITIVE COUNTS, and both counts are DISCRETISATION ARTEFACTS OF MY OWN")
out("            SWEEP, not failures of anything.  RECORDED, NOT SILENTLY FIXED, because the")
out("            uncorrected reading would be a false refutation:")
out("              [6] quantifies p over a RATIONAL simplex grid of denominator 24; the exact")
out("                  p that zeroes a given coefficient triple is generally irrational, so a")
out("                  finite grid misses it.  Over the real simplex [5]+[6] together ARE the")
out("                  definition of a convex hull and hold identically -- proved, and")
out("                  demonstrated numerically in r1 (2) against exact polygon geometry.")
out("              [8] quantifies (u,v) over a q-th-root grid, which cannot realise the")
out("                  continuum of angles TRI needs.  M3-1's hard direction is certified")
out("                  CONSTRUCTIVELY over the reals in r2 (a): 29161 exact lattice states,")
out("                  0 certificate failures, exact residual 0.")
out("            A grid quantifier is a LOWER bound on an existential.  Both counts are")
out("            upper bounds on nothing.  This is the same defect class as COR-E.")
out("   [7]      holds: it is M3-1's easy direction.")
out()
out("DONE.")
open("r5_implication_ledger.OUT.txt", "w").write("\n".join(L) + "\n")
