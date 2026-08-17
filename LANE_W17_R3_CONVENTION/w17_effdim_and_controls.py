#!/usr/bin/env python3
"""
LANE W-17 / ROUTE R3 -- effective dimension under each arm, reproduction of the
pre-cutoff figure the route was posed as needing, and the self-audit of which of
my own tests COULD NOT HAVE FAILED.
"""
import numpy as np
RNG = np.random.default_rng(4242)
V = 5
LOOP_F = [0,1,2]; EDG_F = [0,1,2]
LOOP_C = [0,3,4]; EDG_C = [3,4,5]
A_PUB = np.array([np.pi/3]*3 + [np.pi/2]*3)

def whole_circuit(a, loop, edges):
    W = np.exp(1j*sum(a[e] for e in edges))
    M = np.eye(V, dtype=complex)
    for v in loop: M[v,v] = W
    return M

def edge_tick(a, loop, edges):
    T = np.eye(V, dtype=complex)
    for v in loop: T[v,v] = 0
    L = len(loop)
    for i in range(L):
        T[loop[(i+1)%L], loop[i]] = np.exp(1j*a[edges[i]])
    return T

def fro(X): return float(np.linalg.norm(X,'fro'))

OUT=[]
def P(s=""):
    print(s); OUT.append(str(s))

P("="*78)
P("EFFECTIVE DIMENSION -- WHAT DOES THE FUNCTIONAL ACTUALLY FACTOR THROUGH?")
P("="*78)
P()
MF = whole_circuit(A_PUB, LOOP_F, EDG_F); MC = whole_circuit(A_PUB, LOOP_C, EDG_C)
TF = edge_tick(A_PUB, LOOP_F, EDG_F);     TC = edge_tick(A_PUB, LOOP_C, EDG_C)

def Z(TFo, TCo, x): return complex(np.vdot(TFo@x, TCo@x))

x0 = np.sqrt(np.array([0.30,0.10,0.15,0.25,0.20])) * np.exp(1j*np.array([0.3,1.1,2.2,0.7,1.9]))

P("  compensating move inside one loop's edge phases (a1+e, a2-e): loop holonomy fixed")
for lbl, mk in [("M_gamma (s=L)", lambda a: (whole_circuit(a,LOOP_F,EDG_F), whole_circuit(a,LOOP_C,EDG_C))),
                ("T       (s=1)", lambda a: (edge_tick(a,LOOP_F,EDG_F),  edge_tick(a,LOOP_C,EDG_C)))]:
    worst = 0.0
    for _ in range(500):
        a = A_PUB + RNG.normal(scale=0.4, size=6)
        b = a.copy(); e = RNG.normal(scale=0.3); b[0]+=e; b[1]-=e
        za = Z(*mk(a), x0); zb = Z(*mk(b), x0)
        worst = max(worst, abs(za-zb))
    P(f"    {lbl}: worst |dZ| over 500 draws = {worst:.6e}")
P()
P("  compensating move inside one class's moduli (|x1|^2 + |x2|^2 fixed, both class 10)")
for lbl, (TFo,TCo) in [("M_gamma (s=L)", (MF,MC)), ("T       (s=1)", (TF,TC))]:
    worst = 0.0
    for _ in range(500):
        p = RNG.dirichlet(np.ones(5))
        c10 = p[1]+p[2]; t = RNG.uniform(0.05,0.95)
        q = p.copy(); q[1]=c10*t; q[2]=c10*(1-t)
        ph = RNG.uniform(0,2*np.pi,5)
        za = Z(TFo,TCo, np.sqrt(p)*np.exp(1j*ph)); zb = Z(TFo,TCo, np.sqrt(q)*np.exp(1j*ph))
        worst = max(worst, abs(za-zb))
    P(f"    {lbl}: worst |dZ| over 500 draws = {worst:.6e}")
P()
P("  >>> under M_gamma the functional factors through (class weights, W_F, W_C):")
P("  >>> EFFECTIVE REAL DIMENSION 4 on K1 (2 free class weights + 2 holonomies),")
P("  >>> out of 14 free physical parameters.  Under T nothing is quotiented:")
P("  >>> EFFECTIVE REAL DIMENSION 14.  The binary's two arms differ by 10")
P("  >>> suppressed real dimensions and the question names none of them.")
P()

P("="*78)
P("THE EXPERIMENT W-10's 'NEXT' CLAUSE ASKS FOR -- ALREADY IN W-10's OWN BODY")
P("="*78)
P()
P("  REGISTER_V001.md:1127-1128 / W10:487-489 NEXT: 'run the whole formation")
P("  functional under COR-F's T instead of M_gamma on K1 and on B0b, and see")
P("  whether the incidence becomes visible.  It is the one experiment that would")
P("  tell reading A from reading B.'")
P("  W10_SCOPE_TABLE_V001.md:393 (row 6.7) already reports refuter D-2 having")
P("  run it, on B0b AND K1.  Reproduced here from independent code:")
P()
P(f"    ||[T_F, T_C]||_F on K1              = {fro(TF@TC - TC@TF):.6f}"
  f"    (D-2 printed 2.449490)")
P(f"    ||[M_dF, M_c]||_F on K1             = {fro(MF@MC - MC@MF):.6f}"
  f"    (diagonals always commute)")
P()
p_pub = np.array([0.5,0.0,0.0,0.25,0.25])       # W-01's published ready state
xp = np.sqrt(p_pub).astype(complex)
P("  W-01's HEADLINE, both conventions, on S1's published connection")
P("  (REGISTER_V001.md:35-37: 'IT FIRES ON S1'S OWN PUBLISHED CONNECTION ...")
P("   the overlap is i/2 - i/2 = 0 exactly'):")
P(f"    |<M_dF x, M_c x>|  =  {abs(Z(MF,MC,xp)):.15f}   -> FIRES")
P(f"    |<T_F  x, T_C  x>| =  {abs(Z(TF,TC,xp)):.15f}   -> DOES NOT FIRE")
P()
P("  >>> the corpus's single most-quoted positive result is CONVENTION_SCOPED")
P("  >>> and inverts under the corpus's own sealed alternative.  This is one")
P("  >>> line of code and was available at W-10.")
P()

P("="*78)
P("SELF-AUDIT -- WHICH OF MY TESTS COULD NOT HAVE FAILED?")
P("="*78)
P()
P("  F2 map-1 (||T^L - M_gamma|| = 0): COULD NOT HAVE FAILED AS A DISCOVERY.")
P("    COR-F prints T^3 = diag(W_C,1,1,W_C,W_C) at S3 audit :180-182, which IS")
P("    W-01's M_gamma by definition.  My number is a REPRODUCTION.  What carries")
P("    weight is the DIRECTION it is read in, not the number.  Marked: control.")
P()
P("  F3 biconditional, <= direction (fibre-wise => blind): COULD NOT HAVE FAILED.")
P("    A diagonal unitary acts on |x_v| alone; blindness is immediate.  VOID.")
P("  F3 biconditional, => direction (blind => fibre-wise): COULD have failed and")
P("    did not: 0 of 3000 non-diagonal Haar draws were blind.  CARRIES WEIGHT.")
P()
P("  D-2's 'non-zero exactly when class 11 is occupied' -- the <= half is a")
P("    THEOREM, not a measurement: vertex-disjoint loops have disjoint supports,")
P("    so [T_F,T_C] = 0 identically.  Exhibited:")
LF2=[0,1,2]; LC2=[3,4,5]
V2=6
def et2(loop, edges, a, n):
    T=np.eye(n,dtype=complex)
    for v in loop: T[v,v]=0
    L=len(loop)
    for i in range(L): T[loop[(i+1)%L],loop[i]]=np.exp(1j*a[edges[i]])
    return T
a2 = RNG.uniform(0,2*np.pi,6)
A_=et2(LF2,[0,1,2],a2,V2); B_=et2(LC2,[3,4,5],a2,V2)
P(f"      disjoint loops on 6 vertices: ||[T_F,T_C]||_F = {fro(A_@B_-B_@A_):.3e}")
P("    So B0a's 0.000000 is void as evidence; only the K1/B0b non-zeros carry.")
P()
P("  F1 live-parameter count, F4 lambda spread, F5 firing regions: all could have")
P("    come out otherwise.  CARRY WEIGHT.")
P()
with open("w17_effdim_and_controls.OUT.txt","w") as fh:
    fh.write("\n".join(OUT)+"\n")
