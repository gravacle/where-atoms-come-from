# pin_algebra.py -- LANE W20_PRE.  OPERATOR ALGEBRA AND SYMMETRY.  NO STATE IS BUILT.
#
# Three jobs:
#   (A) the gauge-invariant algebra of S, explicitly, with its centre and block structure;
#   (B) the boundary-algebra choice, and which of the three candidate choices COLLAPSE;
#   (C) Aut(tri_chain12), because a symmetry that maps one declared arm onto another makes that
#       arm's null a THEOREM.  This is the could-not-have-failed defect, hunted BEFORE the run.
#
# Pauli strings are handled symbolically as GF(2) symplectic vectors (a|b) meaning X^a Z^b.
# The only matrices built anywhere are 8x8 (three qubits).  No 2^12 vector is ever allocated.
import itertools
import numpy as np

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)
def rule(t=""):
    P("\n" + "=" * 100); P(t); P("=" * 100)

V = 8
E = [(0, 7), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
     (3, 4), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]
L = len(E)
star = [0] * V
for i, (a, b) in enumerate(E): star[a] |= 1 << i; star[b] |= 1 << i
def bits(m, n=L): return [i for i in range(n) if m >> i & 1]
def pop(m): return bin(m).count("1")
def rank2(vs):
    b = []
    for v in vs:
        for x in b: v = min(v, v ^ x)
        if v: b = sorted(b + [v], reverse=True)
    return len(b), b
def dpart(m):
    r = 0
    for i in bits(m):
        a, b = E[i]; r ^= (1 << a) ^ (1 << b)
    return r
CYC = [m for m in range(1 << L) if dpart(m) == 0]
COC = sorted({(lambda am: __import__("functools").reduce(lambda x, v: x ^ star[v], bits(am, V), 0))(am)
              for am in range(1 << V)})

S_LINKS = [1, 2, 3]; S_MASK = sum(1 << i for i in S_LINKS)
SIGMA = [0, 4, 5];   SIG_MASK = sum(1 << i for i in SIGMA)

# ---------------------------------------------------------------- (A) explicit algebra of S
rule("BLOCK A1 -- THE GAUGE-INVARIANT ALGEBRA OF S, BUILT EXPLICITLY ON 3 QUBITS")
I2 = np.eye(2); Xp = np.array([[0., 1.], [1., 0.]]); Zp = np.diag([1., -1.])
def kron(*m):
    out = np.array([[1.]])
    for x in m: out = np.kron(out, x)
    return out
X = [kron(Xp, I2, I2), kron(I2, Xp, I2), kron(I2, I2, Xp)]   # local qubit 0,1,2 = links 1,2,3
W = kron(Zp, Zp, Zp)
P("local qubit 0,1,2  <->  links 1,2,3 = %s" % [E[i] for i in S_LINKS])
P("generators: X_1, X_2, X_3 (electric, each gauge invariant), and W_S = Z_1 Z_2 Z_3 (the Wilson")
P("loop of the triangle; gauge invariant because {1,2,3} is a cycle).")
P("")
P("commutators, printed rather than asserted:")
for i in range(3):
    P("   || X_%d W_S + W_S X_%d ||  = %.9f   (0 means ANTICOMMUTE)"
      % (i + 1, i + 1, np.linalg.norm(X[i] @ W + W @ X[i])))
for i, j in [(0, 1), (0, 2), (1, 2)]:
    P("   || X_%d X_%d - X_%d X_%d ||  = %.9f   (0 means COMMUTE)"
      % (i + 1, j + 1, j + 1, i + 1, np.linalg.norm(X[i] @ X[j] - X[j] @ X[i])))
    A = X[i] @ X[j]
    P("   || (X_%d X_%d) W_S - W_S (X_%d X_%d) || = %.9f  (0 means the PAIR commutes with W_S)"
      % (i + 1, j + 1, i + 1, j + 1, np.linalg.norm(A @ W - W @ A)))
A3 = X[0] @ X[1] @ X[2]
P("   || (X_1X_2X_3) W_S + W_S (X_1X_2X_3) || = %.9f  (0 means the TRIPLE anticommutes)"
  % np.linalg.norm(A3 @ W + W @ A3))
P("")
P(">>> A_S IS NON-ABELIAN.  It carries a Wilson loop.  Both requirements of the brief, verified.")

# span
gens = {"I": np.eye(8)}
basis = []
for a in range(8):
    for b in range(2):
        M = np.eye(8)
        for k in range(3):
            if a >> k & 1: M = M @ X[k]
        if b: M = M @ W
        basis.append(M.reshape(-1))
Bm = np.array(basis)
P("linear span of {X^a W^b : a in F_2^3, b in F_2} has dimension %d (16 expected)"
  % np.linalg.matrix_rank(Bm))
# centre
cen = []
for a in range(8):
    for b in range(2):
        M = np.eye(8)
        for k in range(3):
            if a >> k & 1: M = M @ X[k]
        if b: M = M @ W
        ok = all(np.linalg.norm(M @ G - G @ M) < 1e-12 for G in (X + [W]))
        if ok: cen.append((a, b))
P("centre Z(A_S) = span of the strings %s"
  % ["X^%s W^%d" % (bits(a, 3), b) for a, b in cen])
P("dim Z(A_S) = %d" % len(cen))
P(">>> Z(A_S) = { I, X_1X_2, X_1X_3, X_2X_3 } -- the EVEN-WEIGHT electric strings on S.")
P("    A_S has dim 16 with a 4-dim centre acting on an 8-dim space:  A_S = M_2 + M_2 + M_2 + M_2 .")
P("    Maximum entropy of a state restricted to A_S = log2(4) + log2(2) = 2 + 1 = 3 bits.")

rule("BLOCK A2 -- THE CENTRE OF A_S *IS* THE SURFACE.  EXACT, 2 BITS, DERIVED.")
P("From the Gauss relations (pin_carrier BLOCK 4), ON THE PHYSICAL SECTOR:")
P("    X_1 X_2 = eta_0 * X_0     X_1 X_3 = eta_1 * X_4     X_2 X_3 = eta_2 * X_5")
P("and X_0 X_4 X_5 = eta_0 eta_1 eta_2 .  So:")
P("")
P("    Z(A_S)  =  alg{ X_0, X_4, X_5 } / (X_0 X_4 X_5 = scalar)  =  the SURFACE ALGEBRA.")
P("    dim 4 on both sides.  This is not an analogy; it is an operator identity.")
P("")
P(">>> CONSEQUENCE, AND IT IS THE WHOLE VACUITY MAP:")
P("    THE RECORD ON S IS EXACTLY 3 BITS.  2 OF THEM ARE THE CENTRE, AND THE CENTRE IS THE")
P("    SURFACE, GAUSS-FORCED, FOR EVERY STATE / COUPLING / SECTOR / HISTORY.  Any criterion that")
P("    reads the centre CANNOT FAIL.  The remaining 1 bit lives inside the M_2 blocks.")
P("    Inside a block the surviving generators are X_1 and W_S, and they ANTICOMMUTE:")
P("       block algebra = alg{X_1, W_S} = M_2 , one electric pole and one magnetic pole.")
P("    (X_2 = x12 X_1 and X_3 = x13 X_1 once the central labels x12, x13 are fixed, so the whole")
P("     electric content collapses to X_1 inside a block.)")
P("")
P(">>> THEREFORE, BEFORE ANY STATE EXISTS:  as g^2 -> infinity the state is an X eigenstate,")
P("    <X_1> -> 1, the block is pure, the LIVE bit is EMPTY.  As g^2 -> 0 the state is a W")
P("    eigenstate, <W_S> -> 1, the block is pure, the LIVE bit is EMPTY AGAIN.  Anticommutation")
P("    forbids both being sharp at once, so the live bit is nonzero ONLY IN THE CROSSOVER.")
P("    THE TEST IS LIVE IN A WINDOW AND VACUOUS AT BOTH ENDS, AND THAT IS DERIVED FROM AN")
P("    ANTICOMMUTATOR, NOT MEASURED.")

# ---------------------------------------------------------------- (B) boundary-algebra choices
rule("BLOCK B -- THE BOUNDARY-ALGEBRA CHOICE: THREE CANDIDATES, ONE COLLAPSE")
def gi_strings(region_mask):
    """count gauge-invariant Pauli strings supported in region: X^a (a in region) x Z^b (b in Z_1, b in region)"""
    zs = [z for z in CYC if not (z & ~region_mask)]
    rz, _ = rank2(zs)
    return (1 << pop(region_mask)) * (1 << rz), pop(region_mask), rz
def gauge_group_in(region_mask):
    """cocycles entirely inside the region: these act as SCALARS, so they quotient the algebra"""
    return rank2([c for c in COC if c and not (c & ~region_mask)])[0]
for nm, rm in [("NARROW  S = E[A] = {1,2,3}", S_MASK),
               ("WIDE    S u Sigma = {0,1,2,3,4,5}", S_MASK | SIG_MASK)]:
    tot, nx, nz = gi_strings(rm)
    q = gauge_group_in(rm)
    P("%-34s : %d X-strings x 2^%d Z-strings = %d ; scalars from cocycles inside: 2^%d ; "
      "dim on physical sector = %d" % (nm, 1 << nx, nz, tot, q, tot >> q))
P("")
P("CENTER choice (Casini-Huerta-Rosabal electric centre): A_S^narrow with X_Sigma adjoined as")
P("central elements.  ON THE PHYSICAL SECTOR X_0 = eta_0 X_1X_2 , X_4 = eta_1 X_1X_3 ,")
P("X_5 = eta_2 X_2X_3 -- all three are ALREADY in Z(A_S^narrow).  Adjoining them adds nothing.")
P(">>> CENTER == NARROW ON THE PHYSICAL SECTOR.  The three-way boundary-algebra choice this")
P("    program inherited COLLAPSES TO TWO on this carrier: NARROW(=CENTER) vs WIDE.")
zw = [z for z in CYC if z and not (z & ~(S_MASK | SIG_MASK))]
P("    WIDE's extra content is magnetic: cycles inside S u Sigma are %s"
  % [bits(z) for z in sorted(zw, key=pop)])
P("    i.e. WIDE adds the SECOND triangle W_{3,4,5} = Z_3Z_4Z_5, an operator not in A_S^narrow.")
P("    NARROW is PINNED PRIMARY.  WIDE is ARM A2 and the gap is at least one magnetic bit.")

# ---------------------------------------------------------------- (C) automorphisms
rule("BLOCK C -- Aut(tri_chain12), AND THE COULD-NOT-HAVE-FAILED TRAP INSIDE MY OWN ARMS")
Eset = {frozenset(e) for e in E}
auts = []
for p in itertools.permutations(range(V)):
    if all(frozenset((p[a], p[b])) in Eset for (a, b) in E):
        auts.append(p)
P("|Aut(tri_chain12)| = %d" % len(auts))
P("the automorphisms:")
for p in auts: P("   %s" % str(p))
TRI_V = [frozenset({0, 1, 2}), frozenset({1, 2, 3}), frozenset({4, 5, 6}), frozenset({5, 6, 7})]
P("")
P("orbits of the four triangle regions under Aut:")
orb = {}
for t in TRI_V:
    img = {frozenset(p[v] for v in t) for p in auts}
    orb[t] = img
for t in TRI_V:
    P("   region on vertices %-12s -> orbit %s" % (sorted(t), sorted([sorted(x) for x in orb[t]])))
P("")
P(">>> READ THIS AS A TRAP WARNING, WHICH IS WHY IT IS COMPUTED BEFORE THE RUN:")
P("    ANY two arms related by an element of Aut MUST return identical record content.  A null")
P("    across such a pair is a THEOREM ABOUT THE SYMMETRY, not a measurement, and scoring it as")
P("    'formation does not matter' would be exactly the W-19 defect in a new costume.")
P("")
A_SET = frozenset({0, 1, 2})
P("stabiliser of the declared region A = {0,1,2} (these fix S and Sigma setwise):")
stab = [p for p in auts if {p[v] for v in A_SET} == set(A_SET)]
for p in stab: P("   %s" % str(p))
P("|Stab(A)| = %d" % len(stab))
P("")
P("CHARGE-ROUTE PAIRS, CHECKED AGAINST Stab(A).  A route pair is VACUOUS if some g in Stab(A)")
P("carries one charge set onto the other.")
routes = {"C0": frozenset(), "C1": frozenset({0, 4}), "C2": frozenset({1, 5}),
          "C3": frozenset({4, 5}), "C4": frozenset({0, 1})}
for a, b in [("C1", "C2"), ("C0", "C3"), ("C0", "C4"), ("C3", "C4"), ("C1", "C4")]:
    hit = [p for p in stab if {p[v] for v in routes[a]} == set(routes[b])]
    P("   %s -> %s : %d symmetry maps in Stab(A) -> %s"
      % (a, b, len(hit), "VACUOUS BY SYMMETRY" if hit else "LIVE (no symmetry forces equality)"))
P("")
P("region pairs for the SURFACE-VARIATION replacement test, checked against Aut:")
for a, b in itertools.combinations(range(4), 2):
    ta, tb = TRI_V[a], TRI_V[b]
    hit = [p for p in auts if {p[v] for v in ta} == set(tb)]
    P("   region%s -> region%s : %d automorphisms -> %s"
      % (sorted(ta), sorted(tb), len(hit),
         "VACUOUS BY SYMMETRY" if hit else "LIVE"))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_W20_PRE/OUT_pin_algebra.txt", "w").write("\n".join(LOG) + "\n")
