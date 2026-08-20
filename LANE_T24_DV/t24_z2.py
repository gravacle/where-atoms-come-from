"""T-24 item 3 -- THE ABELIAN COMPARISON: D(Z_2) ON THE SAME 1x2 TORUS (dim 2^4 = 16).

Same lattice, same machinery, same exact arithmetic as t24_main.py.  The corner-corpus
expectation: the 1x2 torus is a DISTANCE-1 toric code in the vertical dual direction
(the dual vertical cycle crosses a single horizontal edge), so a weight-1 admissible
logical exists and clause (v) must FAIL on the contractible edge.  This is the control
that shows the D(D_4) HOLDS verdict discriminates: same size, same test, opposite answer.
"""
import sys, json, itertools
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV")
from fractions import Fraction
import numpy as np
from t24_lib import *

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T24_DV"
OUT = []
def say(s=""):
    print(s, flush=True); OUT.append(s)

say("=" * 108)
say("T-24 / item 3 -- ABELIAN CONTROL: D(Z_2), 1x2 TORUS (dim 16) -- EXACT")
say("=" * 108)
G = make_Z2()
car = Carrier(G)
n, N = car.n, car.N
MUL = G["MUL"]

# structure
okA = all(np.array_equal(car.permA0(k1)[car.permA0(k2)], car.permA0(int(MUL[k1, k2])))
          for k1 in range(n) for k2 in range(n))
b0, b1 = car.diagB0(), car.diagB1()
okB = all(np.array_equal(b0[car.permA0(k)], b0) and np.array_equal(b1[car.permA1(k)], b1)
          for k in range(n))
say("structure checks (exact): rep %s, [A,B]=0 %s" % (okA, okB))
say("NOTE the Ly=1 abelian degeneracy: both plaquette holonomies reduce to u1 - u0;")
say("the horizontal edges drop out of B entirely, and A_v0(1) = A_v1(1) = X_h0 X_h1 X-star.")

sectors = build_sectors(car)
dims = {}
for s, P in sectors.items():
    tr = sp_trace(P); assert tr == sp_frob_check_projector(P)
    dims[s] = int(tr)
assert sum(dims.values()) == N
eigmap = eigen_projectors(sectors)
eigdims = {k: sum(dims[s] for s in v) for k, v in eigmap.items()}
say("sector dims: %s" % {str(k): v for k, v in sorted(dims.items())})
say("H-eigenspace dims: %s   (note collapsed sectors: A_v0 = A_v1 here)"
    % {("-%d" % k): v for k, v in sorted(eigdims.items())})
eigP = {}
for k, keys in eigmap.items():
    merged = [dict() for _ in range(N)]
    for s in keys:
        P = sectors[s]
        for j in range(N):
            for i, v in P[j].items():
                w = merged[j].get(i, Fraction(0)) + v
                if w: merged[j][i] = w
                else: merged[j].pop(i, None)
    eigP[k] = merged

# record: horizontal Wilson loop Zbar_h = sign(h0)*sign(h1); vertical Zbar_v = sign(u0)
sign = G["irreps"][1][1]
r_H = sign[MUL[car.h0, car.h1]]
r_V0 = sign[car.u0]
Hent = h_entries_exact(car)
say("")
say("RECORDS, clause-verified exactly (D-18):")
recs = {}
for nm, r in [("Zbar_h=sign(h0*h1)", r_H), ("Zbar_v0=sign(u0)", r_V0)]:
    bit = bool(np.all(np.abs(r) == 1))
    durable = not any(i != j and r[i] != r[j] for (i, j), v in Hent.items())
    bal = {}; nt = {}
    for k in sorted(eigP):
        P = eigP[k]
        bal[-k] = sp_trace_diag(P, r)
        acc = Fraction(0)
        for j in range(N):
            for i, v in P[j].items():
                w = P[i].get(j)
                if w: acc += v * w * int(r[i]) * int(r[j])
        nt[-k] = acc
    say("  %-20s (i) %s (ii) %s (iii) nontrivial %s (iv) balanced %s ; Tr(P_E R)=%s"
        % (nm, bit, durable, any(x != 0 for x in nt.values()),
           all(x == 0 for x in bal.values()), {k: str(x) for k, x in bal.items()}))
    recs[nm] = (r, dict(i=bit, ii=durable, iii=any(x != 0 for x in nt.values()),
                        iv=all(x == 0 for x in bal.values())))

say("")
say("C_e (exact integer Gram + Fraction kernel; D-21-clean):")
Ce = {}
for e in ["h0", "h1", "u0", "u1"]:
    basis, grank, _ = edge_admissible_algebra(car, e)
    for M in basis:
        ok, _n = verify_in_commutant(car, e, M)
        assert ok
    Ce[e] = basis
    say("  C_%s: dim %d" % (e, len(basis)))

say("")
say("THE TEST on the contractible edges {h0}, {h1} (same reduction as the D_4 run):")
X_h0_mat = [[Fraction(1) if a == int(MUL[b, 1]) else Fraction(0) for b in range(n)]
            for a in range(n)]
flip_table = {}
for e in ["h0", "h1"]:
    tau = car.edge_perm_right(e, 1)     # X on edge e = translation by the nontrivial element
    tE = {}
    for k in sorted(eigP):
        t = sp_trace_perm(eigP[k], tau)
        assert t.denominator == 1
        tE[-k] = int(t)
    flip_table[e] = tE
    say("  edge %s: t_E = Tr(P_E X_%s) = %s  -> ALL ZERO: %s" % (e, e, tE,
        all(v == 0 for v in tE.values())))
# the flip, exhibited exactly: U = X_h0 is admissible and flips Zbar_h
ok_adm, _ = verify_in_commutant(car, "h0", X_h0_mat)
tau = car.edge_perm_right("h0", 1)
flips = np.array_equal(r_H[tau], -r_H)
say("  U = X_h0 (weight-1, on the CONTRACTIBLE edge h0):")
say("     admissible ([X_h0, H] = 0, exact): %s ;  U^dag Zbar_h U = -Zbar_h (exact): %s"
    % (ok_adm, flips))
assert ok_adm and flips
say("  => CLAUSE (v) FAILS on the abelian 1x2: an admissible weight-1 operation on a single")
say("     contractible region flips the record Zbar_h.  This is the known distance-1 proxy")
say("     artifact (the dual vertical cycle has length Ly = 1), scoped per O-49/D-23.")

say("")
say("positive control bookkeeping (D-15): on this carrier the FLIP ITSELF is the firing")
say("control; the non-admissible analogue X_u0 flips Zbar_v0 while [X_u0,H] != 0:")
lam = car.edge_perm_right("u0", 1)
flips2 = np.array_equal(r_V0[lam], -r_V0)
accH = {}
Ui = np.empty(N, dtype=np.int64); Ui[lam] = np.arange(N)
for (i, j), v in Hent.items():
    k1 = (int(lam[i]), j); accH[k1] = accH.get(k1, Fraction(0)) + v
    k2 = (i, int(Ui[j])); accH[k2] = accH.get(k2, Fraction(0)) - v
cn2 = sum(v * v for v in accH.values())
say("  X_u0 flips Zbar_v0: %s ;  ||[X_u0,H]||_F^2 = %s (nonzero -> not admissible)" % (flips2, cn2))

say("")
say("NON-CONTRACTIBLE single edges on the abelian carrier (contrast with the D_4 nulls):")
say("Zbar_v0 = sign(u0) is balanced (above), so by the flip reduction a record anticommuting")
say("with the ADMISSIBLE diag Z_u0 exists; explicitly, Xbar_v = X_u0 X_u1:")
rho = car.edge_perm_right("u0", 1)[car.edge_perm_right("u1", 1)]
inv_ok = np.array_equal(rho[rho], np.arange(N))
h_inv = True
for (i, j), v in Hent.items():
    if Hent.get((int(rho[i]), int(rho[j]))) != v: h_inv = False
bal_x = {}
for k in sorted(eigP):
    t = sp_trace_perm(eigP[k], rho)
    bal_x[-k] = t
nt_x = Fraction(0)
for k in sorted(eigP):
    P = eigP[k]
    for j in range(N):
        for i, v in P[j].items():
            w = P[int(rho[j])].get(int(rho[i]))
            if w: nt_x += v * w
zdiag = sign[car.u0]
anti = np.array_equal(zdiag[rho], -zdiag)
Zmat = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(-1)]]
z_adm, _ = verify_in_commutant(car, "u0", Zmat)
say("  Xbar_v: involution %s; [R,H]=0 (H conj-invariant under it): %s; Tr(P_E R) = %s;"
    % (inv_ok, h_inv, {k: str(v) for k, v in bal_x.items()}))
say("  Z_u0 admissible (exact): %s;  {Z_u0, Xbar_v} = 0 (exact, rho flips sign(u0)): %s"
    % (z_adm, anti))
say("  => on the abelian 1x2, admissible single-edge flippers exist on the NON-contractible")
say("     edges too (Z_u0 flips Xbar_v), where D(D_4) has none.")

say("")
say("CONTROL VERDICT: same lattice, same test, ABELIAN answer = FAILS (weight-1 admissible")
say("flipper exists, all t_E = 0), NON-ABELIAN answer = HOLDS (t_E != 0 obstruction).")
say("The discriminating difference is exact: for Z_2 the central translation grades every")
say("eigenspace evenly; for D_4 the central r^2 lies in [G,G], conjugation fixes classes,")
say("and the grading is uneven.  Protection on this proxy is a NON-ABELIAN effect.")

with open(LANE + "/t24_z2.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
json.dump(dict(carrier="D(Z_2) 1x2 torus", dim=N,
               eig_dims={str(-k): v for k, v in eigdims.items()},
               Ce_dims={e: len(Ce[e]) for e in Ce},
               tE_X={e: {str(k): v for k, v in flip_table[e].items()} for e in flip_table},
               records={nm: c for nm, (r, c) in recs.items()},
               clause_v="FAILS on contractible edges (X_h0 admissible weight-1 flipper)"),
          open(LANE + "/t24_z2.json", "w"), indent=1)
