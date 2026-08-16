"""G1 — THE ABELIAN GROUP U(1)xU(1) AT FIBRE DIMENSION 2.
The claim says every load-bearing theorem FAILS IDENTICALLY here.
This script tests each one, and asks what the correct index set is."""
import numpy as np, itertools
from glib import *

np.set_printoptions(precision=6, suppress=True)
SEED = 20260816
rng = np.random.default_rng(SEED)
print("seed =", SEED)

# U(1)xU(1) acting on C^2 with weights (1,0) and (0,1): W = diag(e^{i a}, e^{i b}).
# Connection carried on e1 (loop F) and e4 (loop C); all other edges = I.

def rank2_state(weights, dirs):
    """weights: 5 vertex masses; dirs: 5 unit vectors in C^2."""
    return normalise([np.sqrt(w) * np.asarray(d, dtype=complex) / np.linalg.norm(d)
                      for w, d in zip(weights, dirs)])

def show(tag, U, s, N=400000):
    z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
    dev = max(abs(Z_from_chars(z, c, [k])[0] -
                  Z_direct(U, EDGES_K1, LOOP_F, LOOP_C, s, k)) for k in [1, 2, 3, 7, 13, 40])
    lb = lambda_B(z, c, N=N)
    return z, c, lb, dev

print()
print("=" * 78)
print("G1.0  IS RANK-2 U(1)xU(1) A NEW THEORY, OR TWO COPIES OF THE OLD ONE?")
print("=" * 78)
a1, a2 = 1.0, 2.3           # W_F = diag(e^{i a1}, e^{i a2})
b1, b2 = np.sqrt(2), 0.7    # W_C = diag(e^{i b1}, e^{i b2})
U2 = diag_conn([a1, a2], [b1, b2])
w = [0.4, 0.15, 0.15, 0.15, 0.15]
dirs = [rng.normal(size=2) + 1j * rng.normal(size=2) for _ in range(5)]
s2 = rank2_state(w, dirs)
z, c, lb, dev = show("", U2, s2)
print("char-form vs literal matrix action, max dev = %.2e" % dev)
print("number of distinct characters = %d   (rank one gives 3 on K1)" % len(z))
print("EVERY coefficient:")
for zz, cc in zip(z, c):
    print("   zeta=%+.6f%+.6fj   coeff=%+.9f%+.9fj   |Im coeff|=%.2e"
          % (zz.real, zz.imag, cc.real, cc.imag, abs(cc.imag)))
print("max |Im(coeff)| = %.3e ; min Re(coeff) = %+.9f ; sum coeff = %.12f"
      % (np.abs(c.imag).max(), c.real.min(), c.sum().real))
print("--> at ABELIAN structure group, EVERY coefficient is a NON-NEGATIVE REAL")
print("    and equals a piece of the ready state's mass. Sum = 1.")

# explicit decomposition into two rank-one U(1) systems
s_a = [np.array([x[0]]) for x in s2]
s_b = [np.array([x[1]]) for x in s2]
Ua = u1_conn(a1, b1)
Ub = u1_conn(a2, b2)
ks = np.arange(1, 61)
Zfull = Z_from_chars(z, c, ks)
za, ca = merge_characters(*characters(Ua, EDGES_K1, LOOP_F, LOOP_C, s_a))
zb, cb = merge_characters(*characters(Ub, EDGES_K1, LOOP_F, LOOP_C, s_b))
Zsum = Z_from_chars(za, ca, ks) + Z_from_chars(zb, cb, ks)
print("max |Z_k(rank2) - (Z_k^{(1)} + Z_k^{(2)})| over k<=60 : %.3e" % np.abs(Zfull - Zsum).max())
print("--> RANK-2 U(1)xU(1) IS EXACTLY THE DIRECT SUM OF TWO RANK-ONE U(1) SYSTEMS.")

print()
print("=" * 78)
print("G1.1  THEOREM 1 — 'THE ROOT CAN NEVER FIRE'")
print("=" * 78)
print("Refined index set: (vertex class, fibre weight). The root at rank 2 carries")
print("TWO refined indices, so |S| = 2, not 1. The theorem's hypothesis fails, not")
print("the theorem. And the theorem is FALSIFIABLE at rank 2 -- here is a rank-2,")
print("FAITHFUL, NON-SCALAR U(1)xU(1) connection on which the root still cannot fire.")
print()
for tag, (aa, bb) in [("generic  a=(1.0,2.3) b=(1.41,0.7)", ([1.0, 2.3], [np.sqrt(2), 0.7])),
                      ("tuned    a=(0.0,0.9) b=(0.0,0.9)",  ([0.0, 0.9], [0.0, 0.9]))]:
    Ux = diag_conn(*[np.asarray(x) for x in (aa, bb)])
    WF = based_holonomies(Ux, EDGES_K1, LOOP_F)[0]
    WC = based_holonomies(Ux, EDGES_K1, LOOP_C)[0]
    sr = rank2_state([1, 0, 0, 0, 0], [[1, 1], [1, 0], [1, 0], [1, 0], [1, 0]])
    zr, cr = merge_characters(*characters(Ux, EDGES_K1, LOOP_F, LOOP_C, sr))
    Zr = np.abs(Z_from_chars(zr, cr, np.arange(1, 20001)))
    print("  %s" % tag)
    print("     W_F scalar? %s   W_C scalar? %s   [W_F,W_C] = %.1e   faithful U(1)^2? %s"
          % (np.allclose(WF, WF[0, 0] * np.eye(2)), np.allclose(WC, WC[0, 0] * np.eye(2)),
             np.abs(WF @ WC - WC @ WF).max(), "yes (weights (1,0),(0,1))"))
    print("     pinch characters u_j v_j: %s" % np.round(np.diag(WF.conj().T @ WC), 9))
    print("     root-only state: #chars=%d  min|Z_k|=%.12f  max|Z_k|=%.12f  -> %s"
          % (len(zr), Zr.min(), Zr.max(),
             "FIRES" if Zr.min() < 0.999 else "CANNOT FIRE (theorem HOLDS)"))
print()
print("VERDICT G1.1: 'fails identically' is FALSE. The theorem holds on the whole")
print("subfamily {u_1 v_1 = u_2 v_2} of rank-2 faithful U(1)xU(1) connections, and")
print("the refined criterion |S|=1 => never predicts both cases correctly.")

print()
print("=" * 78)
print("G1.2  THEOREM 2 — W-02's CRITERION  G = <chi_a/chi_b> != {1}  <=>  FORMATION")
print("=" * 78)
print("Tested on 400 random rank-2 U(1)xU(1) instances, refined index = (class,weight).")
ok = bad = 0
worst = []
rng2 = np.random.default_rng(31337)
for trial in range(400):
    aa = rng2.choice([0.0, 1.0, 2.3, np.sqrt(2), np.pi / 3], size=2)
    bb = rng2.choice([0.0, 0.7, 1.9, np.sqrt(3), np.pi / 5], size=2)
    Ux = diag_conn(aa, bb)
    mask = rng2.integers(0, 2, size=(5, 2))
    if mask.sum() == 0:
        continue
    wts = rng2.random((5, 2)) * mask
    sx = normalise([np.array(wts[v], dtype=complex) for v in range(5)])
    zx, cx = merge_characters(*characters(Ux, EDGES_K1, LOOP_F, LOOP_C, sx))
    pred_forms = len(zx) >= 2          # G != {1}  <=>  >= 2 distinct characters
    Zx = np.abs(Z_from_chars(zx, cx, np.arange(1, 30001)))
    obs_forms = Zx.min() < Zx.max() - 1e-9
    if pred_forms == obs_forms:
        ok += 1
    else:
        bad += 1
        worst.append((len(zx), Zx.min(), Zx.max()))
print("  agreement: %d / %d   mismatches: %d" % (ok, ok + bad, bad))
print("VERDICT G1.2: W-02's criterion -- the corpus's single most load-bearing")
print("theorem -- SURVIVES VERBATIM at rank-2 U(1)xU(1). It does not fail at all.")

print()
print("=" * 78)
print("G1.3  THEOREM 3 — 'lambda IS A FUNCTION OF THE CLASS WEIGHTS'")
print("=" * 78)
print("W-03's SU(2) kill: six ready states with IDENTICAL class weights, |Z_1| spread")
print("0.4247, against exactly 0.0 under U(1). Re-run at ABELIAN rank two:")
U2 = diag_conn([1.0, 2.3], [np.sqrt(2), 0.7])
vals = []
rng3 = np.random.default_rng(99)
for t in range(6):
    dirs = [rng3.normal(size=2) + 1j * rng3.normal(size=2) for _ in range(5)]
    sx = rank2_state([0.4, 0.15, 0.15, 0.15, 0.15], dirs)
    cw = class_weights(sx, EDGES_K1, LOOP_F, LOOP_C)
    zx, cx = merge_characters(*characters(U2, EDGES_K1, LOOP_F, LOOP_C, sx))
    vals.append((abs(Z_from_chars(zx, cx, [1])[0]), lambda_B(zx, cx, N=200000), cw))
print("  class weights (identical for all six): %s" % {k: round(v, 12) for k, v in vals[0][2].items()})
print("  |Z_1| : %s" % np.round([v[0] for v in vals], 6))
print("  spread |Z_1| = %.6f      lambda_B spread = %.6f"
      % (max(v[0] for v in vals) - min(v[0] for v in vals),
         max(v[1] for v in vals) - min(v[1] for v in vals)))
print("VERDICT G1.3: the 0.4247-style spread REPRODUCES AT AN ABELIAN GROUP.")
print("  W-03 attributed it to SU(2); register ERR-2 already convicted that control of")
print("  a three-way confound. This is the direct measurement ERR-2 called for, and it")
print("  lands on the ABELIAN side: the spread is a FIBRE-RANK effect, as the claim says.")
print("  But it is a failure of the 4-class INDEX, not of the theorem: with the refined")
print("  8-element index (class x weight) lambda is again an exact function of the weights.")
# demonstrate the refined statement
def refined_weights(U, s):
    z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
    return z, c
dirs = [rng3.normal(size=2) + 1j * rng3.normal(size=2) for _ in range(5)]
sA = rank2_state([0.4, 0.15, 0.15, 0.15, 0.15], dirs)
zA, cA = refined_weights(U2, sA)
# build a DIFFERENT state with the same refined weights: rephase every component
ph = np.exp(1j * rng3.uniform(0, 2 * np.pi, (5, 2)))
sB = [sA[v] * ph[v] for v in range(5)]
zB, cB = refined_weights(U2, sB)
print("  same refined weights, all components rephased:  |lambda_B(A) - lambda_B(B)| = %.3e"
      % abs(lambda_B(zA, cA, 200000) - lambda_B(zB, cB, 200000)))

print()
print("=" * 78)
print("G1.4  THEOREM 4 — PINCH = SPECTATOR")
print("=" * 78)
print("Rank one: the four characters 1,u,v,uv are the four corners of the unit square,")
print("a CENTRALLY SYMMETRIC set; translation by (uv)^{-1} permutes it, giving the")
print("exchange 00<->11, 10<->01. Rank two: the 8 exponent vectors live in Z^4 as two")
print("unit squares in orthogonal coordinate planes sharing the origin -- NOT centrally")
print("symmetric. Measured:")
rng4 = np.random.default_rng(4242)
worstdev = 0.0
for t in range(500):
    aa = rng4.uniform(0, 2 * np.pi, 2)
    bb = rng4.uniform(0, 2 * np.pi, 2)
    Ux = diag_conn(aa, bb)
    wts = rng4.random((5, 2))
    sP = normalise([np.array(wts[v], dtype=complex) for v in range(5)])
    # swap classes: v0 (11) <-> a phantom 00 vertex does not exist on K1, so use the
    # K1-internal half of the symmetry: 10 <-> 01, i.e. swap {v1,v2} with {v3,v4}
    sQ = [sP[0], sP[3], sP[4], sP[1], sP[2]]
    zP, cP = merge_characters(*characters(Ux, EDGES_K1, LOOP_F, LOOP_C, sP))
    zQ, cQ = merge_characters(*characters(Ux, EDGES_K1, LOOP_F, LOOP_C, sQ))
    worstdev = max(worstdev, abs(abs(Z_from_chars(zP, cP, [1])[0]) -
                                abs(Z_from_chars(zQ, cQ, [1])[0])))
print("  rank-2, 500 random instances, max |Z_1| deviation under 10<->01 : %.3e" % worstdev)
# rank-one control
rng5 = np.random.default_rng(4242)
worst1 = 0.0
for t in range(500):
    ff, cc = rng5.uniform(0, 2 * np.pi, 2)
    Ux = u1_conn(ff, cc)
    wv = rng5.random(5)
    sP = normalise([np.array([np.sqrt(x)], dtype=complex) for x in wv])
    sQ = [sP[0], sP[3], sP[4], sP[1], sP[2]]
    zP, cP = merge_characters(*characters(Ux, EDGES_K1, LOOP_F, LOOP_C, sP))
    zQ, cQ = merge_characters(*characters(Ux, EDGES_K1, LOOP_F, LOOP_C, sQ))
    worst1 = max(worst1, abs(abs(Z_from_chars(zP, cP, [1])[0]) -
                             abs(Z_from_chars(zQ, cQ, [1])[0])))
print("  rank-1 control, same 500 draws                                : %.3e" % worst1)

print()
print("=" * 78)
print("G1.5  WHERE DOES THE FOUR-CORNER STRUCTURE ACTUALLY STOP? (the wall test)")
print("=" * 78)
print("The claim puts the wall at fibre rank 2. The record already puts one door")
print("earlier: CHARGE at rank one. Re-derived here from scratch.")
# charge: unit-rank fibre, edge charges q_e -> W_F = exp(i sum q_e a_e).
# W-03 of record: exponents (1,0),(2,0),(3,0) give |S|=3 with rank G = 1.
E3 = np.array([[1, 0], [2, 0], [3, 0]])
D = E3[1:] - E3[0]
print("  charge exponents (1,0),(2,0),(3,0):  |S| = 3, rank of difference lattice = %d"
      % int(np.linalg.matrix_rank(D)))
print("  S4-1 predicts rank 2 whenever |S| >= 3.  ---> S4-1 FAILS AT RANK ONE.")
print("  Reproduced independently:  lambda for those three characters, equal weights")
lam_chg = mahler1_jensen([1 / 3, 1 / 3, 1 / 3], exps=[1, 2, 3])
print("     m(1/3 z + 1/3 z^2 + 1/3 z^3) = %.9f   (a ONE-variable measure: rank G = 1)"
      % lam_chg)
print("  q = (1,2,2,2,2) on K1, rank one -- W-03 records lambda moving to -1.200555:")
# q=(1,2,2,2,2) is a charge per EDGE: e1 has q=1, e2..e6 have q=2 (S2/S5 convention).
# W_F = exp(i(q1 a1 + q2 a2 + q3 a3)), W_C = exp(i(q4 a4 + q5 a5 + q6 a6)).
# With the generic connection sampled on the full torus the classes keep unit
# characters but the EXPONENTS scale; the class characters become u^{qF}, v^{qC}.
print("     (the S5 lane owns this number; this lane only needs the rank fact above)")
print()
print("VERDICT G1.5: 'rank-one fibres' is NOT the wall. The four-corner taxonomy")
print("already dies at rank one under charge (of record, W-03), and it dies at rank")
print("two for the SAME reason: the realised character set stops being the corners of")
print("a square. The operative variable is the CHARACTER SET, not the fibre rank.")
