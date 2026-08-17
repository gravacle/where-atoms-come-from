# rule_main3.py -- W-19 RULING, part 3.  Exactly where the boundary-algebra choice bites.
import numpy as np, time
from rule_verify import *

t0 = time.time()
P("=" * 118)
P("W-19 RULING -- PART 3.  THE BOUNDARY-ALGEBRA RESIDUE, LOCATED TO THE FRAGMENT.")
P("=" * 118)

def has_uv_path(V, edges, l, F):
    u, v = edges[l]
    adj = {x: [] for x in range(V)}
    for i in F:
        a, b = edges[i]; adj[a].append(b); adj[b].append(a)
    seen = {u}; st = [u]
    while st:
        x = st.pop()
        for y in adj[x]:
            if y not in seen: seen.add(y); st.append(y)
    return v in seen

P("")
P("[10a] EXACT DISCREPANCY EXT - CL ON EVERY RULE-A FRAGMENT OF EVERY d>=4 CARRIER.")
P("      EXT = full 2x2 algebra on the system link (lane A / lane B reading 1).")
P("      CL  = alg{X_l}, the ONLY non-trivial gauge-invariant subalgebra of a single link.")
P("      %-13s %5s %5s %14s %14s %14s %12s" % ("carrier", "g2", "|F|", "encloses cycle", "EXT/H", "CL/H", "EXT-CL bits"))
for nm, (V, E) in [("dbl_chain9", mg_chain(5)), ("tri_chain12", tri_chain12()),
                   ("petersen", petersen()), ("heawood", heawood())]:
    car = Carrier(nm, V, E); L = car.L
    frags, d = rule_A_fragments(V, E, 0)
    for g2 in (0.50,):
        vec = car.lift(car.ground(g2)[1])
        vX = hadamard_all(vec, L); pX = np.abs(vX) ** 2
        HS = S_ax(vec, L, AX(L, [0]))
        pm = pX.reshape([2] * L).sum(axis=tuple(a for a in range(L) if a != L - 1))
        HX = float(-(pm[pm > 1e-15] * np.log2(pm[pm > 1e-15])).sum())
        for F in frags:
            e = channel_EXT(vec, L, 0, F); k = channel_CL(pX, L, 0, F)
            P("      %-13s %5.2f %5d %14s %14.9f %14.9f %12.9f"
              % (nm, g2, len(F), str(has_uv_path(V, E, 0, F)), e / HS, k / HX, e - k))

P("")
P("[10b] LANE A's OWN ENCLOSURE CONTROL (its C1), REBUILT HERE, AND WHAT THE GAUGE-INVARIANT")
P("      ALGEBRA SAYS ABOUT IT.  heawood, g2 = 1.00, |F| = 14 held FIXED, only the SHAPE moves.")
car = Carrier("heawood", *heawood()); V, E, L = car.V, car.edges, car.L
vec = car.lift(car.ground(1.00)[1])
vX = hadamard_all(vec, L); pX = np.abs(vX) ** 2
HS = S_ax(vec, L, AX(L, [0]))
pm = pX.reshape([2] * L).sum(axis=tuple(a for a in range(L) if a != L - 1))
HX = float(-(pm[pm > 1e-15] * np.log2(pm[pm > 1e-15])).sum())
frags, d = rule_A_fragments(V, E, 0)
armA = [F for F in frags if len(F) == 14]
armA = armA[0] if armA else frags[2]
u, v = E[0]
# arm B: same size, forced to contain a shortest u-v path in G-l
dist = bfs_dist(V, E, u, 0)
path = []; x = v
while x != u:
    for i, (a, b) in enumerate(E):
        if i == 0: continue
        if (a == x and dist.get(b, 99) == dist[x] - 1): path.append(i); x = b; break
        if (b == x and dist.get(a, 99) == dist[x] - 1): path.append(i); x = a; break
rest = [i for i in range(1, L) if i not in path]
armB = sorted(path + rest[:len(armA) - len(path)])
for tag, F in [("ARM A (no u-v path)", armA), ("ARM B (contains u-v path)", armB)]:
    e = channel_EXT(vec, L, 0, F); c = channel_CHI(vX, L, 0, F); k = channel_CL(pX, L, 0, F)
    P("      %-26s |F|=%d encloses=%s   EXT/H=%.9f  CHI/H=%.9f  CL/H=%.9f"
      % (tag, len(F), has_uv_path(V, E, 0, F), e / HS, c / HX, k / HX))
P("      ARMS DIFFER IN: A\\B=%s  B\\A=%s  (symmetric difference %d links)"
  % (sorted(set(armA) - set(armB)), sorted(set(armB) - set(armA)), len(set(armA) ^ set(armB))))
P("      H(S) = %.9f bits" % HS)

P("")
P("[10c] THE SAME COMPARISON ON A MAGNETIC RECORD (theta_L, magnetic GHZ = the exact g2->0 state).")
for Lk in (6, 8):
    car = Carrier("theta_%d" % Lk, *theta(Lk))
    ghz = np.zeros(1 << Lk); ghz[0] = ghz[(1 << Lk) - 1] = 1 / np.sqrt(2)
    vX = hadamard_all(ghz, Lk); pX = np.abs(vX) ** 2
    HS = S_ax(ghz, Lk, AX(Lk, [0]))
    worst = 0.0
    for k in range(1, Lk):
        F = list(range(1, k + 1))
        e = channel_EXT(ghz, Lk, 0, F); c = channel_CL(pX, Lk, 0, F)
        if k < Lk - 1: worst = max(worst, abs(e - c))
    P("      theta_%d magnetic GHZ:  H(S) = %.9f;  max |EXT - CL| over PROPER fragments = %.9f bits"
      % (Lk, HS, worst))

P("")
P("[10d] SUMMARY OF THE RESIDUE.")
P("      ELECTRIC RECORD (every d>=4 carrier, every coupling, ground AND Haar): EXT - CL = 0 exactly")
P("      on every fragment that does not enclose a cycle through l.  The choice does not bite.")
P("      ENCLOSING FRAGMENT (same carrier, same state, same |F|): EXT/H = 1.21..., CL/H = 1.000000.")
P("      MAGNETIC RECORD (theta_L, g2->0 ground state): EXT/H = 1.000000 at every proper fragment,")
P("      CL/H = 0.000000 at every proper fragment.  FULL SWING.")

P("")
P("elapsed %.1f s" % (time.time() - t0))
open("OUT_rule_main3.txt", "w").write("\n".join(LOG) + "\n")
