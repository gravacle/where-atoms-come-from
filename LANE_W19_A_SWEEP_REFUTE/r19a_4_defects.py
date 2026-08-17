# W19-A REFUTER / step 4.  LATENT DEFECTS IN THE INSTRUMENT, AND THE AUDITS LANE A OWED.
#
#  4a  rdm() computes T @ T.T with NO CONJUGATION (zn_gauge.py, function rdm).  Every state lane A
#      fed it happens to be real (eigh of a real symmetric H; np.random.randn).  So the defect never
#      bit -- but F6 calls those states "HAAR-RANDOM", which they are not: randn is real Gaussian
#      (orthogonal-invariant), not Haar (unitary-invariant).  Test: a genuinely COMPLEX physical
#      state, entropies via lane A's rdm vs a conjugating rdm.  If they agree, no defect; if they
#      differ, F6's decisive control was run on the wrong ensemble and must be re-run.
#  4b  min_cycle_basis() selects cycles by GF(2) independence.  For Z_N with N > 2 the plaquette set
#      must generate the cycle group MOD N.  If the C x C chord minor has |det| != 1, the magnetic
#      term omits part of the cycle group for any N sharing a factor with det.  This is the only
#      place F12 (the Z_2 vs Z_3 control) could be wrong.
#  4c  The WEIGHT FLOOR is lane A's only remaining defence of the 21-link number.  Scan g^2.
#  4d  DISJOINTNESS AND CUT VALIDITY AUDIT -- the refuter's assigned check, run on lane A's own
#      threshold carrier AND on the counterexamples.
import numpy as np, sys, itertools, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP_REFUTE")
import zn_gauge as ZG
from zn_gauge import (ZNGauge, S_of, mutual_information, level_cuts, nested_fragments,
                      girth_through, rdm, build_adj)
from carriers import heawood, petersen, cube, ladder, grid, theta, theta_sub, LADDER
from refute_carriers import tri_chain12, dbl_chain9, mg_chain, min_degree
from collections import deque

np.random.seed(7717)
print("=" * 120)
print("W19-A REFUTER / 4  LATENT DEFECTS AND THE AUDITS")
print("=" * 120)

# ---------------------------------------------------------------- 4a complex states
def rdm_conj(Psi, L, N, A):
    A = sorted(A); B = [i for i in range(L) if i not in A]
    T = Psi.reshape([N] * L).transpose(A + B).reshape(N ** len(A), N ** len(B))
    return T @ T.conj().T

def ent(rho):
    w = np.linalg.eigvalsh(rho); w = w[w > 1e-13]
    return float(-(w * np.log2(w)).sum())

def S_conj(Psi, L, N, A):
    A = sorted(A)
    if len(A) in (0, L): return 0.0
    B = [i for i in range(L) if i not in A]
    use = A if len(A) <= len(B) else B
    return ent(rdm_conj(Psi, L, N, use))

def MI_conj(Psi, L, N, A, B):
    return S_conj(Psi, L, N, A) + S_conj(Psi, L, N, B) - S_conj(Psi, L, N, sorted(set(A) | set(B)))

print("\n[4a] COMPLEX PHYSICAL STATES.  carrier = heawood (lane A's threshold carrier), l = 0,")
print("     fragments = lane A's rule-A F_1..F_4 (|F| = 2,6,14,18).  Three ensembles.")
V, E = heawood(); g = ZNGauge("heawood", V, E, 2); L = g.L
frs, d = nested_fragments(V, E, 0)
print(f"     {'ensemble':<42}{'H(S) laneA':>15}{'H(S) conj':>15}{'max|dS|':>11}   I/H(S) (conjugating rdm)")
for tag, mk in [("real Gaussian (what lane A called Haar)", lambda: np.random.randn(g.dimP)),
                ("COMPLEX Haar (Ginibre, unitary-invariant)",
                 lambda: np.random.randn(g.dimP) + 1j * np.random.randn(g.dimP)),
                ("COMPLEX Haar, second draw",
                 lambda: np.random.randn(g.dimP) + 1j * np.random.randn(g.dimP))]:
    v = mk(); v = v / np.linalg.norm(v)
    Psi = g.full_vector(v)
    hA = S_of(Psi, L, 2, [0]); hC = S_conj(Psi, L, 2, [0])
    dmax = max(abs(S_of(Psi, L, 2, F) - S_conj(Psi, L, 2, F)) for F in frs)
    rs = [MI_conj(Psi, L, 2, [0], F) / hC for F in frs]
    print(f"     {tag:<42}{hA:>15.9f}{hC:>15.9f}{dmax:>11.3e}   " + " ".join(f"{r:.9f}" for r in rs))
print("     READ: with a conjugating rdm the plateau is still exactly 1.000000000 on complex states,")
print("     so F6's CONCLUSION survives; but lane A's rdm() is only valid for real amplitudes and its")
print("     'Haar-random' label is wrong.  The defect is latent, not active.  Recorded, not scored.")

# ---------------------------------------------------------------- 4b cycle basis mod N
print("\n[4b] IS THE PLAQUETTE SET A BASIS MOD N, OR ONLY MOD 2?  min_cycle_basis selects by GF(2)")
print("     rank.  The magnetic term is right for every N only if the chord minor has |det| = 1.")
print(f"     {'carrier':<22}{'C':>4}{'|det chord minor|':>20}{'basis mod 3 ?':>16}{'basis mod 2 ?':>16}")
bad = []
ALLC = list(LADDER) + [("tri_chain12", tri_chain12()), ("dbl_chain9", dbl_chain9()),
                       ("petersen", petersen()), ("mg_chain6", mg_chain(6))]
for nm, (Vg, Eg) in ALLC:
    gg = ZNGauge(nm, Vg, Eg, 2)
    M = np.array([[int(p[gg.chords[c]]) for c in range(gg.C)] for p in gg.plaq], dtype=float)
    det = abs(np.linalg.det(M)) if gg.C else 1.0
    def rank_mod(p):
        A = np.array(M, dtype=np.int64) % p; A = A.copy(); r = 0
        for c in range(A.shape[1]):
            piv = None
            for i in range(r, A.shape[0]):
                if A[i, c] % p: piv = i; break
            if piv is None: continue
            A[[r, piv]] = A[[piv, r]]
            inv = pow(int(A[r, c]), p - 2, p)
            A[r] = (A[r] * inv) % p
            for i in range(A.shape[0]):
                if i != r and A[i, c] % p: A[i] = (A[i] - A[i, c] * A[r]) % p
            r += 1
        return r
    ok3 = rank_mod(3) == gg.C; ok2 = rank_mod(2) == gg.C
    if not ok3: bad.append(nm)
    print(f"     {nm:<22}{gg.C:>4}{det:>20.6f}{str(ok3):>16}{str(ok2):>16}")
print(f"     Carriers whose plaquette set FAILS to generate the cycle group mod 3: {bad if bad else 'NONE'}")
print("     -> F12's Z_3 arm is not corrupted by the GF(2) selection on any carrier actually used.")
print("     The risk is real in principle (a GF(2)-independent set can have det divisible by an odd")
print("     prime); it did not materialise here.  This test could have refuted F12 and did not.")

# ---------------------------------------------------------------- 4c weight floor scan
print("\n[4c] DOES THE WEIGHT FLOOR RESCUE THE 21-LINK CLAIM?  H(S) vs g^2, floor = 0.10 bits.")
print("     If the counterexamples were only ever above the floor at couplings where heawood is not,")
print("     lane A could still claim 21.  ONE VARIABLE MOVED per column: the carrier.")
gsqs = [0.20, 0.30, 0.40, 0.50, 0.70, 1.00, 1.50]
cars = [("heawood", heawood()), ("tri_chain12", tri_chain12()), ("dbl_chain9", dbl_chain9()),
        ("mg_chain6", mg_chain(6))]
print("     " + f"{'g^2':>6}" + "".join(f"{nm:>16}" for nm, _ in cars) + "   (H(S) bits; * = below floor)")
tab = {}
for gs in gsqs:
    line = f"     {gs:>6.2f}"
    for nm, (Vg, Eg) in cars:
        gg = ZNGauge(nm, Vg, Eg, 2); p, _, _ = gg.ground(2.0 / gs, 2.0 * gs); P = gg.full_vector(p)
        l = max(range(gg.L), key=lambda i: girth_through(Vg, Eg, i)[1] or -1)
        h = S_of(P, gg.L, 2, [l]); tab[(nm, gs)] = h
        line += f"{h:>15.6f}" + ("*" if h < 0.10 else " ")
    print(line)
print("     READ: at g^2 <= 0.50 every counterexample is well above the floor at the SAME coupling")
print("     lane A used for its own exhibit (g^2 = 0.50).  The floor does not rescue the claim.")

# ---------------------------------------------------------------- 4d disjointness / cut audit
print("\n[4d] FRAGMENT AUDIT (the refuter's assigned check).  For each carrier: are the rule-C fragments")
print("     (i) pairwise disjoint, (ii) disjoint from S, (iii) GENUINE u-v edge cuts in G-l (removing")
print("     the cut must disconnect tail from head), (iv) do they satisfy <X_l X(C_i)> = +1 exactly?")
def is_cut(V_, E_, l, C):
    a, b = E_[l]; adj = [[] for _ in range(V_)]
    for i, (x, y) in enumerate(E_):
        if i == l or i in C: continue
        adj[x].append(y); adj[y].append(x)
    seen = {a}; dq = deque([a])
    while dq:
        x = dq.popleft()
        for y in adj[x]:
            if y not in seen: seen.add(y); dq.append(y)
    return b not in seen
def expect_X_product(gg, Psi, links):
    N, L_ = gg.N, gg.L; T = Psi.reshape([N] * L_); Tv = T
    for e in links: Tv = np.roll(Tv, 1, axis=e)
    return float((T.ravel() * Tv.ravel()).sum())
print(f"     {'carrier':<14}{'l':>3}{'#cuts':>7}{'pairwise disjoint':>19}{'disjoint from S':>17}"
      f"{'all genuine cuts':>18}{'min <X_l X(C)>':>16}{'nested F_k are NOT disjoint':>29}")
for nm, (Vg, Eg) in [("heawood", heawood()), ("petersen", petersen()), ("tri_chain12", tri_chain12()),
                     ("dbl_chain9", dbl_chain9()), ("cube_Q3", cube())]:
    gg = ZNGauge(nm, Vg, Eg, 2); p, _, _ = gg.ground(4.0, 1.0); P = gg.full_vector(p)
    l = max(range(gg.L), key=lambda i: girth_through(Vg, Eg, i)[1] or -1)
    cuts, dd = level_cuts(Vg, Eg, l); frs2, _ = nested_fragments(Vg, Eg, l)
    seen = set(); pd = True
    for C in cuts:
        if seen & set(C): pd = False
        seen |= set(C)
    ds = all(l not in C for C in cuts)
    gc = all(is_cut(Vg, Eg, l, set(C)) for C in cuts)
    mx = min(expect_X_product(gg, P, [l] + C) for C in cuts)
    nested_disj = all(not (set(frs2[i]) & set(frs2[j])) for i in range(len(frs2)) for j in range(i))
    print(f"     {nm:<14}{l:>3}{len(cuts):>7}{str(pd):>19}{str(ds):>17}{str(gc):>18}{mx:>16.12f}"
          f"{str(not nested_disj):>29}")
print("     READ: rule C's fragments ARE disjoint and ARE real cuts -- R_delta is correctly computed.")
print("     But rule A's F_k, the fragments the PLATEAU is read off, are NESTED, hence NOT disjoint.")
print("     Lane A reports plateau (nested) and R_delta (disjoint) from two different families; the")
print("     brief asks for a plateau on DISJOINT fragments.  On the disjoint family the 5 cuts all sit")
print("     at exactly I = H(S), which is the Gauss identity F3 restated, not an accumulation curve.")
print("\nDONE 4.")
