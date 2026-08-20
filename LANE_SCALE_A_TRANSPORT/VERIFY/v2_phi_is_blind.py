"""V2 -- WHAT IS phi_fix ACTUALLY A FUNCTION OF?

Three adversarial questions, each answered by a computation, not an argument.

 Q1  Is phi_fix an ORDINARY quantity?  I recompute it by a completely different route:
     phi = [ sum_E <chi_E,chi_E> ] / [ sum_E chi_E(e)^2 ],  where <chi,chi> = (1/|G|) sum_h
     |chi_E(h)|^2 is Burnside/Plancherel -- the number of orbits of the transport group on
     pairs, i.e. dim End_G(E).  If this matches the lane's isotypic sum on every carrier then
     phi_fix is a normalised SECOND MOMENT OF THE TRANSPORT CHARACTER and nothing else.

 Q2  Can phi_fix respond to the record count AT ALL?  The record count the lane plots is
     m = dim E(-2), and E(-2) carries the TRIVIAL transport rep (chi_{E(-2)}(h) = m for all h).
     So E(-2) contributes m^2 to the numerator and m^2 to the denominator: RATIO EXACTLY 1.
     Printed per eigenspace.  If the ground sector's own contribution is identically 1 on
     every carrier, the statistic is constitutionally blind to the variable being plotted.

 Q3  D-16 IN ITS OWN SHAPE: the lane ranks phi_fix against the RAW record count while the
     carrier dimension varies by a factor of 64 in the same column.  The scale-free variable
     is the record DENSITY m/dim.  Ranked against that instead -- same 35 carriers, same
     numbers -- does the "never the record count" reading survive?
"""
import sys, numpy as np
from collections import defaultdict
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT"
sys.path.insert(0, LANE)
import glib
from carriers import census, isotypic, phi
def say(*a): print(*a); sys.stdout.flush()

def spearman(x, y):
    def rk(v):
        o = np.argsort(v, kind='stable'); r = np.empty(len(v), float); r[o] = np.arange(len(v))
        # average ties
        v = np.asarray(v, float); out = r.copy()
        for u in set(v.tolist()):
            idx = np.where(v == u)[0]
            out[idx] = r[idx].mean()
        return out
    a, b = rk(x), rk(y)
    a = a - a.mean(); b = b - b.mean()
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))

rows = []
say("="*130)
say("V2   WHAT IS phi_fix A FUNCTION OF?   (Q1 method cross-check, Q2 blindness, Q3 D-16 normalisation)")
say("="*130)
say("")
say("  Q1/Q2 TABLE -- phi by two independent routes, and the per-eigenspace contribution.")
say(f"  {'carrier':<13}{'|G|':>4}{'dim':>6}{'RECORDS m':>10}{'|G/Z|':>6}{'phi(isotypic)':>14}"
    f"{'phi(<chi,chi>)':>15}{'diff':>10}   per-E ratio  <chi,chi>/dim^2   [E=-2, -1, 0]")
for G in glib.ladder(64):
    ce = census(G); iso = isotypic(G, ce)
    p_lane = phi(iso, ce['dims'])[2]
    num = 0.0; den = 0.0; per = []
    for v in (-2, -1, 0):
        d = ce['dims'][v]
        if d == 0: per.append(None); continue
        ip = float(np.sum(ce['chis'][v] ** 2) / G.n)     # <chi_E, chi_E>
        num += ip; den += d * d
        per.append(ip / (d * d))
    p_mom = num / den
    ps = ", ".join("--" if x is None else f"{x:.4f}" for x in per)
    say(f"  {G.name:<13}{G.n:>4}{G.n*G.n:>6}{ce['dims'][-2]:>10}{G.n//ce['Z']:>6}{p_lane:>14.10f}"
        f"{p_mom:>15.10f}{p_mom-p_lane:>10.1e}   [{ps}]")
    rows.append(dict(name=G.name, n=G.n, dim=G.n*G.n, m=ce['dims'][-2], GZ=G.n//ce['Z'],
                     ab=G.abelian, phi=p_lane, ground_ratio=per[0]))
say("")
say(f"  MAX |phi(isotypic) - phi(<chi,chi>)| over all 35 carriers: "
    f"{max(abs(r['phi']-r['phi']) for r in rows):.1e}   (both routes printed above, column 'diff')")
gr = sorted({round(r['ground_ratio'], 12) for r in rows})
say(f"  GROUND-SECTOR CONTRIBUTION TO phi, over all 35 carriers: {gr}")
say("     -> the sector whose dimension IS the plotted 'record count' contributes the constant 1.")

say("")
say("="*130)
say("  Q3   D-16: RANK phi_fix AGAINST THE SCALE-FREE RECORD DENSITY m/dim, NOT THE RAW COUNT")
say("="*130)
say(f"  {'carrier':<13}{'abel':>6}{'dim':>6}{'RECORDS m':>10}{'m/dim':>9}{'|G/Z|':>6}{'phi_fix':>9}")
for r in sorted(rows, key=lambda r: -r['m']/r['dim']):
    say(f"  {r['name']:<13}{str(r['ab']):>6}{r['dim']:>6}{r['m']:>10}{r['m']/r['dim']:>9.4f}{r['GZ']:>6}{r['phi']:>9.4f}")
say("")
for lab, sel in (("ALL 35 CARRIERS", rows), ("NON-ABELIAN ONLY (n=20)", [r for r in rows if not r['ab']])):
    ph = [r['phi'] for r in sel]
    say(f"  {lab}")
    say(f"     spearman(phi_fix, RAW record count m)   = {spearman(ph,[r['m'] for r in sel]):+.4f}")
    say(f"     spearman(phi_fix, RECORD DENSITY m/dim) = {spearman(ph,[r['m']/r['dim'] for r in sel]):+.4f}")
    say(f"     spearman(phi_fix, |G/Z|)                = {spearman(ph,[r['GZ'] for r in sel]):+.4f}")
    say(f"     spearman(phi_fix, dim)                  = {spearman(ph,[r['dim'] for r in sel]):+.4f}")
    say(f"     spearman(m/dim, |G/Z|)                  = {spearman([r['m']/r['dim'] for r in sel],[r['GZ'] for r in sel]):+.4f}")
    say("")
say("  PSEUDO-REPLICATION CHECK on the lane's n=20 spearman: how many DISTINCT (phi, m/dim, |G/Z|)")
say("  triples do the 20 non-abelian carriers actually supply?")
na = [r for r in rows if not r['ab']]
tri = sorted({(round(r['phi'],6), round(r['m']/r['dim'],6), r['GZ']) for r in na})
say(f"     20 carriers -> {len(tri)} distinct triples: {tri}")
say(f"     distinct phi values: {sorted({round(r['phi'],4) for r in na})}")
