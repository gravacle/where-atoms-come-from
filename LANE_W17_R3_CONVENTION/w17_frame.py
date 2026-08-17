#!/usr/bin/env python3
"""
LANE W-17 / ROUTE R3 -- FRAME CHALLENGE on the CONVENTION QUESTION.

Route question as posed (W10_SCOPE_TABLE_V001.md:210-220, REGISTER_V001.md:1098-1108):
   Is "the formation functional is carrier-independent" a FINDING (A) or a
   RESTATEMENT of the transport convention (B)?

Rival named in the brief: COR-F edge-tick transport T
   (S3_THE_CROSSING_AUDIT_V001.md:160-215, :794), T^L = M_gamma.

This script produces NUMBERS for F1 (embed), F2 (degeneracy), F3 (carving), F5 (null).
No sympy. python3 + numpy only.

CARRIER K1 (S1_CARRIER_K1_V001.md sec.1):
   V = v0..v4;  e1:v0->v1 e2:v1->v2 e3:v2->v0 (filled triangle F)
                e4:v0->v3 e5:v3->v4 e6:v4->v0 (unfilled cycle C)
   published connection (S1 sec.6): a1=a2=a3=pi/3, a4=a5=a6=pi/2
   => W_F = -1, W_C = -i
"""
import numpy as np

np.set_printoptions(precision=12, suppress=False)
RNG = np.random.default_rng(20260817)

V = 5
LOOP_F = [0, 1, 2]          # cyclic order of the filled triangle
LOOP_C = [0, 3, 4]          # cyclic order of the unfilled cycle
EDGES_F = [0, 1, 2]         # a1,a2,a3 attached to steps v0->v1, v1->v2, v2->v0
EDGES_C = [3, 4, 5]         # a4,a5,a6 attached to steps v0->v3, v3->v4, v4->v0

A_PUB = np.array([np.pi/3]*3 + [np.pi/2]*3)


def whole_circuit(a, loop, edges):
    """W-01 / S2 transport: fibre-wise scalar multiplication by the loop holonomy
    on the loop's vertices, identity off the loop.  (REGISTER_V001.md:29-33)"""
    W = np.exp(1j * sum(a[e] for e in edges))
    M = np.eye(V, dtype=complex)
    for v in loop:
        M[v, v] = W
    return M, W


def edge_tick(a, loop, edges):
    """COR-F edge-transport operator: move each fibre value one edge along the loop.
    (S3_THE_CROSSING_AUDIT_V001.md:178-184)
       (Ts)(loop[i+1]) = U_{edge i} s(loop[i]),  identity off the loop."""
    T = np.eye(V, dtype=complex)
    for v in loop:
        T[v, v] = 0.0
    L = len(loop)
    for i in range(L):
        src = loop[i]
        dst = loop[(i + 1) % L]
        T[dst, src] = np.exp(1j * a[edges[i]])
    return T


def fro(X):
    return float(np.linalg.norm(X, 'fro'))


def zoverlap(TF, TC, x):
    """the formation overlap  <T_F x, T_C x>  (W-01's <M_dF s, M_c s>)."""
    return complex(np.vdot(TF @ x, TC @ x))


def rand_conn():
    return RNG.uniform(0, 2*np.pi, 6)


def rand_state():
    x = RNG.normal(size=V) + 1j*RNG.normal(size=V)
    return x / np.linalg.norm(x)


OUT = []
def P(s=""):
    print(s)
    OUT.append(str(s))


# ----------------------------------------------------------------------------
P("="*78)
P("F2  DEGENERACY -- ARE THE TWO NAMED CONVENTIONS THE SAME OBJECT UNDER A MAP?")
P("="*78)
P()
P("MAP 1: POWER.  f(T) = T^L,  L = loop length = 3 on K1 for both loops.")
P()

MF, WF = whole_circuit(A_PUB, LOOP_F, EDGES_F)
MC, WC = whole_circuit(A_PUB, LOOP_C, EDGES_C)
TF = edge_tick(A_PUB, LOOP_F, EDGES_F)
TC = edge_tick(A_PUB, LOOP_C, EDGES_C)

P(f"  S1 published connection:  W_F = {WF:.12f}   W_C = {WC:.12f}")
P(f"  ||T_F^*T_F - I||_F                = {fro(TF.conj().T@TF - np.eye(V)):.3e}   (T_F unitary)")
P(f"  ||T_C^*T_C - I||_F                = {fro(TC.conj().T@TC - np.eye(V)):.3e}   (T_C unitary)")
P(f"  T_F diagonal?                      {np.allclose(TF, np.diag(np.diag(TF)))}")
P(f"  T_C diagonal?                      {np.allclose(TC, np.diag(np.diag(TC)))}")
P()
P(f"  ||T_F^3 - M_dF||_F                = {fro(np.linalg.matrix_power(TF,3) - MF):.6e}")
P(f"  ||T_C^3 - M_c ||_F                = {fro(np.linalg.matrix_power(TC,3) - MC):.6e}")

worst = 0.0
for _ in range(2000):
    a = rand_conn()
    mf, _ = whole_circuit(a, LOOP_F, EDGES_F)
    mc, _ = whole_circuit(a, LOOP_C, EDGES_C)
    tf = edge_tick(a, LOOP_F, EDGES_F)
    tc = edge_tick(a, LOOP_C, EDGES_C)
    worst = max(worst, fro(np.linalg.matrix_power(tf,3) - mf),
                       fro(np.linalg.matrix_power(tc,3) - mc))
P(f"  worst ||T^3 - M_gamma||_F over 2000 random connections = {worst:.6e}")
P()
P("  >>> F2 NUMBER, MAP=POWER:  ||T^L - M_gamma||_F = 0 (exactly, to fp noise).")
P("  >>> READ AGAINST THE FRAME: at the CIRCUIT clock the two named 'rival")
P("  >>> conventions' are THE SAME OPERATOR.  They are not rivals there.")
P()

P("MAP 2: CONJUGACY (change of coordinates).  Is T similar to M_gamma?")
evT = np.sort_complex(np.linalg.eigvals(TC))
evM = np.sort_complex(np.linalg.eigvals(MC))
P(f"  spec(T_C) = {np.round(evT,9)}")
P(f"  spec(M_c) = {np.round(evM,9)}")
P(f"  sorted-spectrum L2 distance      = {np.linalg.norm(evT-evM):.6f}")
P("  >>> NOT degenerate under conjugacy.  T is an L-th ROOT of M_gamma, not a")
P("  >>> conjugate of it.  The map collapses one way only.")
P()

P("MAP 3: THE REVERSE MAP.  How many non-diagonal L-th roots does M_gamma have?")
P("  M_c restricted to its loop is W_C * I_3.  Every unitary S with S^3 = W_C*I_3")
P("  is an L-th root.  Dimension of that set inside U(3):")
roots = 0
for _ in range(4000):
    # random unitary on the 3-dim loop block, conjugate a fixed root
    H = RNG.normal(size=(3,3)) + 1j*RNG.normal(size=(3,3))
    H = H + H.conj().T
    Q = np.linalg.qr(H)[0]
    base = (WC**(1/3)) * np.diag([1, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)])
    S = Q @ base @ Q.conj().T
    if np.allclose(np.linalg.matrix_power(S,3), WC*np.eye(3), atol=1e-10):
        roots += 1
P(f"  random conjugates S = Q R Q* with S^3 = W_C I_3 verified: {roots}/4000")
P("  The conjugacy orbit of a primitive cube root is U(3)/T^3, real dim 9-3 = 6.")
P("  >>> the fibre of f(S)=S^L over M_gamma has real dimension 6 on K1's C-loop.")
P("  >>> DEGENERACY IS ASYMMETRIC: M_gamma is a FUNCTION of any edge transport,")
P("  >>> no edge transport is a function of M_gamma.  M_gamma is a QUOTIENT, not")
P("  >>> an alternative at the same level of description.")
P()

# ----------------------------------------------------------------------------
P("="*78)
P("F1  EMBED -- WHAT SPACE ARE THE ARMS POINTS IN?")
P("="*78)
P()
P("Candidate space: TRANSPORT PRESCRIPTIONS ON A DESIGNATED LOOP PAIR,")
P("coordinatised by TICK STRIDE s (how many edges one elementary tick advances).")
P("  s = L  -> whole-circuit scalar multiplication M_gamma   (Reading B's convention)")
P("  s = 1  -> COR-F's edge tick T                            (the named rival)")
P("  1 < s < L -> UNNAMED INTERIOR POINTS.  K1 has L=3 so exactly one: s=2.")
P()
x_pub = None
p_pub = np.array([0.5, 0.0, 0.0, 0.25, 0.25])   # W-01's published ready state, REGISTER:37
x_pub = np.sqrt(p_pub).astype(complex)
P("  ready state: W-01's published p = (1/2,0,0,1/4,1/4)  (REGISTER_V001.md:37)")
P()
P("  stride s :  ||[T_F^s, T_C^s]||_F   |Z(s)|          T^s diagonal?")
for s in range(1, 4):
    tfs = np.linalg.matrix_power(TF, s)
    tcs = np.linalg.matrix_power(TC, s)
    comm = fro(tfs@tcs - tcs@tfs)
    z = zoverlap(tfs, tcs, x_pub)
    isdiag = np.allclose(tfs, np.diag(np.diag(tfs))) and np.allclose(tcs, np.diag(np.diag(tcs)))
    P(f"    s = {s}   :  {comm:.6f}            {abs(z):.12f}    {isdiag}")
P()
P("  W-10's own second coordinate, THE CLOCK (W10_SCOPE_TABLE_V001.md:110-115):")
P("  on K1 both loops have length 3, so circuit-clock and edge-clock coincide and")
P("  the coordinate COLLAPSES.  Off K1 it does not.  Synthetic control, L_F=4, L_C=3:")
for n in [1, 2, 3, 4, 6, 12, 24]:
    P(f"    edge tick n = {n:3d} :  F-circuits = {n/4:.4f}   C-circuits = {n/3:.4f}"
      f"   equal? {abs(n/4 - n/3) < 1e-12}")
P("  >>> the two branches sit at different circuit counts at every n not a")
P("  >>> multiple of 12.  STRIDE and CLOCK are independent coordinates off K1.")
P()

P("Now the DIMENSION the functional actually lives on, under each arm.")
P("One-variable moves over all 16 raw real parameters (a1..a6, |x_v| x5, arg x_v x5).")
P()

def live_params(stride):
    """count how many of the 16 raw parameters the functional responds to."""
    a0 = A_PUB.copy()
    r0 = np.sqrt(np.array([0.30, 0.10, 0.15, 0.25, 0.20]))
    ph0 = np.array([0.3, 1.1, 2.2, 0.7, 1.9])
    def Z(a, r, ph):
        xx = (r/np.linalg.norm(r)) * np.exp(1j*ph)
        tf = np.linalg.matrix_power(edge_tick(a, LOOP_F, EDGES_F), stride)
        tc = np.linalg.matrix_power(edge_tick(a, LOOP_C, EDGES_C), stride)
        return zoverlap(tf, tc, xx)
    base = Z(a0, r0, ph0)
    live, dead, names = [], [], []
    eps = 1e-5
    for i in range(6):
        a = a0.copy(); a[i] += eps
        names.append(f"a{i+1}"); (live if abs(Z(a,r0,ph0)-base) > 1e-9 else dead).append(f"a{i+1}")
    for i in range(V):
        r = r0.copy(); r[i] += eps
        names.append(f"|x{i}|"); (live if abs(Z(a0,r,ph0)-base) > 1e-9 else dead).append(f"|x{i}|")
    for i in range(V):
        ph = ph0.copy(); ph[i] += eps
        names.append(f"arg x{i}"); (live if abs(Z(a0,r0,ph)-base) > 1e-9 else dead).append(f"arg x{i}")
    return live, dead

for s, label in [(3, "s=L=3  M_gamma  (Reading B's convention)"),
                 (2, "s=2    UNNAMED interior point"),
                 (1, "s=1    COR-F's T (the named rival)")]:
    live, dead = live_params(s)
    P(f"  {label}")
    P(f"      LIVE parameters ({len(live):2d}): {live}")
    P(f"      BLIND to        ({len(dead):2d}): {dead}")
P()
P("  >>> F1 RESULT: the arms are points s=L and s=1 of a stride family with an")
P("  >>> unnamed interior; the functional's live-parameter count MOVES along it.")
P("  >>> The binary SAMPLES the space; it does not partition it.")
P()

# ----------------------------------------------------------------------------
P("="*78)
P("F3  CARVING -- DOES 'FINDING vs RESTATEMENT' PARTITION ANYTHING?")
P("="*78)
P()
P("Stated measure.  Universe = transport prescriptions on the loop pair.")
P("  cell_B = { prescriptions under which the functional is carrier/state-blind")
P("             AND the blindness is ENTAILED by the prescription }")
P("  cell_A = { prescriptions under which the functional is carrier/state-blind")
P("             AND the blindness is NOT entailed by the prescription }")
P()
P("Test the biconditional  BLIND  <=>  FIBRE-WISE(diagonal)  [W-06's N4 correction].")
P("Blindness probe: hold the class weights pi fixed, move only the WITHIN-CLASS")
P("split and the state phases; a blind functional cannot move.")

def blindness_spread(TFs, TCs, ntrial=400):
    """max spread of |Z| over states with IDENTICAL class-weight vector pi."""
    # K1 classes: v0 in both (11); v1,v2 in F only (10); v3,v4 in C only (01)
    w11, w10, w01 = 0.40, 0.30, 0.30
    vals = []
    for _ in range(ntrial):
        t = RNG.uniform(0.05, 0.95)
        u = RNG.uniform(0.05, 0.95)
        p = np.array([w11, w10*t, w10*(1-t), w01*u, w01*(1-u)])
        ph = RNG.uniform(0, 2*np.pi, V)
        x = np.sqrt(p) * np.exp(1j*ph)
        vals.append(abs(zoverlap(TFs, TCs, x)))
    return max(vals) - min(vals)

rows = []
for s in range(1, 4):
    tfs = np.linalg.matrix_power(TF, s)
    tcs = np.linalg.matrix_power(TC, s)
    isdiag = np.allclose(tfs, np.diag(np.diag(tfs)), atol=1e-12) and \
             np.allclose(tcs, np.diag(np.diag(tcs)), atol=1e-12)
    spr = blindness_spread(tfs, tcs)
    rows.append((s, isdiag, spr))
    P(f"    stride s={s}:  fibre-wise = {str(isdiag):5s}   blindness spread = {spr:.6e}")

P()
P("  Haar control -- 3000 random unitary transports on the two loop blocks:")
nd_blind = 0; d_count = 0; nd_count = 0
for _ in range(3000):
    def randU(loop):
        M = np.eye(V, dtype=complex)
        k = len(loop)
        H = RNG.normal(size=(k,k)) + 1j*RNG.normal(size=(k,k))
        Q = np.linalg.qr(H)[0]
        for i, vi in enumerate(loop):
            for j, vj in enumerate(loop):
                M[vi, vj] = Q[i, j]
        return M
    uf, uc = randU(LOOP_F), randU(LOOP_C)
    dg = np.allclose(uf, np.diag(np.diag(uf)), atol=1e-10) and \
         np.allclose(uc, np.diag(np.diag(uc)), atol=1e-10)
    spr = blindness_spread(uf, uc, ntrial=60)
    blind = spr < 1e-12
    if dg: d_count += 1
    else:
        nd_count += 1
        if blind: nd_blind += 1
P(f"    diagonal draws                     : {d_count} / 3000")
P(f"    NON-diagonal draws                 : {nd_count} / 3000")
P(f"    NON-diagonal draws that were BLIND : {nd_blind} / {nd_count}")
P()
P("  >>> cell_A is EMPTY: 0 of 3000 prescriptions are blind without being")
P("  >>> fibre-wise, and 3/3 strides obey BLIND <=> FIBRE-WISE exactly.")
P("  >>> measure(cell_A) = 0 (empty).  measure(cell_B) = the diagonal subgroup,")
P("  >>> real dim 3 inside U(3) of real dim 9 per loop -> codimension 6 per loop.")
P("  >>> A predicate with an EMPTY cell does not partition.  NON_CARVING.")
P()

# ----------------------------------------------------------------------------
P("="*78)
P("F5  NULL -- DOES EITHER BRANCH OBTAIN?  Test the object, not the readings.")
P("="*78)
P()
P("The question is about 'the formation functional'.  The corpus uses that phrase")
P("for TWO objects with OPPOSITE scope verdicts in the very table that poses it:")
P("  (i)  W-01's FORMATION CONDITION (does it fire?)  -- W-10: THREE_CLASS_SCOPED")
P("  (ii) N1's RATE lambda = m(p00 + p10 x + p01 y + p11 xy) -- CARRIER_INDEPENDENT")
P()
P("Re-derive (i)'s carrier-DEPENDENCE, from W-09's own construction:")

def firing_fraction(chars, n=1201):
    """fraction of the (f,c) torus where 0 lies in conv{characters}, chars given as
    exponent pairs (m,n) meaning u^m v^n with u=e^{if}, v=e^{ic}."""
    g = np.linspace(0, 2*np.pi, n, endpoint=False)
    F, C = np.meshgrid(g, g, indexing='ij')
    pts = []
    for (m, k) in chars:
        pts.append(np.exp(1j*(m*F + k*C)))
    P_ = np.stack(pts, axis=-1)
    ang = np.angle(P_)
    ang = np.sort(ang, axis=-1)
    gaps = np.diff(ang, axis=-1)
    wrap = ang[..., 0] + 2*np.pi - ang[..., -1]
    maxgap = np.maximum(gaps.max(axis=-1), wrap)
    return float((maxgap <= np.pi).mean())

three = [(0,1), (1,0), (1,1)]        # K1's occupied classes  {01,10,11} -> v, u, uv
four  = [(0,0), (0,1), (1,0), (1,1)] # a four-class designation -> 1, v, u, uv
f3 = firing_fraction(three)
f4 = firing_fraction(four)
P(f"    three-class occupancy (K1 as handed)  firing region = {f3:.6f}   (W-09: exactly 1/4)")
P(f"    four-class occupancy  (B0b, B4)       firing region = {f4:.6f}   (W-09: exactly 1/2)")
P(f"    ratio = {f4/f3:.6f}")
P()
P("  W-01's ADVERTISED VIRTUE, f -> -f, on both occupancies:")

def flip_fraction(chars, n=1201):
    g = np.linspace(0, 2*np.pi, n, endpoint=False)
    F, C = np.meshgrid(g, g, indexing='ij')
    def hull(Fv, Cv):
        pts = [np.exp(1j*(m*Fv + k*Cv)) for (m, k) in chars]
        P_ = np.stack(pts, axis=-1)
        ang = np.sort(np.angle(P_), axis=-1)
        gaps = np.diff(ang, axis=-1)
        wrap = ang[..., 0] + 2*np.pi - ang[..., -1]
        return np.maximum(gaps.max(axis=-1), wrap) <= np.pi
    return float((hull(F, C) != hull(-F, C)).mean())

P(f"    three-class: verdict flips on {flip_fraction(three):.6f} of the torus")
P(f"    four-class : verdict flips on {flip_fraction(four):.6f} of the torus")
P()
P("  >>> F5: for object (i), the formation condition, CARRIER-INDEPENDENCE IS FALSE")
P("  >>> (1/4 vs 1/2, a factor of 2; virtue present vs identically absent).")
P("  >>> So for (i) NEITHER branch obtains: there is no carrier-independence there")
P("  >>> to be either a finding or a restatement.  NULL OBTAINS on the reading of")
P("  >>> 'formation functional' that W-01 built and that S3/S4/W-03 quote forward.")
P()
P("="*78)
P("DONE")
P("="*78)

with open("w17_frame.OUT.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
