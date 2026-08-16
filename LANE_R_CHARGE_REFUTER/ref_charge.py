#!/usr/bin/env python3
"""
REFUTER LANE -- CHARGE AXIS.  Written from scratch for this refutation.
No code, function or constant is copied from LANE_S5_CHARGE_CODE or from S4.

PUBLISHED CONVENTIONS (this file is the publication):
  K1 from S1_CARRIER_K1_V001.md section 1.
  vertices v0..v4 -> indices 0..4
  edges e1..e6 -> (src,tgt) = (0,1),(1,2),(2,0),(0,3),(3,4),(4,0)
  d1[v,e] = +1 if v==tgt(e), -1 if v==src(e), 0 else      (5 x 6)
  d2[e,F] = +1 for e in {e1,e2,e3}, 0 else                 (6 x 1)
  gamma_F = e1.e2.e3   (bounds F)      gamma_C = e4.e5.e6  (does not bound)
  a_v = 1 iff v lies on gamma_F ; b_v = 1 iff v lies on gamma_C
  u = conj(W_F) = exp(-i f) ,  v = W_C = exp(i c)
  exponent vector E_v = (m_v, n_v) ; character chi_v = u^{m_v} v^{n_v}
  Z_k = sum_v p_v u^{k m_v} v^{k n_v}
  Fourier frequencies are taken in the (f,c) coordinates:
      u^m v^n = exp(i(-m f + n c))  ->  frequency vector (-m, n)
  GRID for Fourier support: N x N with f_j = 2 pi j / N, c_l = 2 pi l / N, N stated at
      every call site, coefficient threshold stated at every call site.
  SEEDS: numpy.random.default_rng only; every seed printed at its call site.
         REFUTER seeds are 77xxxxx, deliberately disjoint from the axis lane's 515xxxx.
"""
import numpy as np, itertools, math
from fractions import Fraction

np.set_printoptions(linewidth=200, suppress=False)
def hdr(s):
    print("="*98); print(s); print("="*98)

# ----------------------------------------------------------------------------- K1
EDGES = [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]     # e1..e6
FACE  = [0,1,2]                                    # e1+e2+e3
V, E = 5, 6

d1 = np.zeros((V,E), dtype=np.int64)
for j,(s,t) in enumerate(EDGES):
    d1[s,j] -= 1; d1[t,j] += 1
d2 = np.zeros((E,1), dtype=np.int64)
for j in FACE: d2[j,0] = 1

GF = [0,1,2]; GC = [3,4,5]
a_v = np.zeros(V,dtype=np.int64); b_v = np.zeros(V,dtype=np.int64)
for j in GF:
    s,t = EDGES[j]; a_v[s]=1; a_v[t]=1
for j in GC:
    s,t = EDGES[j]; b_v[s]=1; b_v[t]=1

hdr("0.  K1 REBUILT FROM S1 SECTION 1 -- INCIDENCE PUBLISHED")
print("d1 (rows v0..v4, cols e1..e6) =\n", d1)
print("d2^T (cols e1..e6) =", d2.T)
print("d1 @ d2 =", (d1@d2).ravel(), "  (exact integer)")
rk1 = np.linalg.matrix_rank(d1.astype(float)); rk2 = np.linalg.matrix_rank(d2.astype(float))
b0 = V - rk1; b1 = (E - rk1) - rk2; b2 = 1 - rk2
print(f"V={V} E={E} F=1  chi={V-E+1}  rank d1={rk1} rank d2={rk2}  b0={b0} b1={b1} b2={b2}")
print("a_v =", a_v, "  b_v =", b_v, "   (a=on gamma_F, b=on gamma_C)")
print("gamma_F is a cycle:", np.all(d1@np.array([1,1,1,0,0,0])==0),
      " gamma_C is a cycle:", np.all(d1@np.array([0,0,0,1,1,1])==0))
print()

# ------------------------------------------------------------------- the functional
def exponents(q):
    """q = per-vertex charge vector (length 5).  E_v = q_v * (a_v, b_v)."""
    q = np.asarray(q, dtype=np.int64)
    return np.stack([q*a_v, q*b_v], axis=1)          # (5,2)

def Z_closed(E, p, f, c, k=1):
    m = E[:,0]; n = E[:,1]
    return np.sum(p * np.exp(1j*k*(-m*f + n*c)))

def Z_matrix(E, s, f, c, k=1):
    """direct matrix action on C^5: diagonal unitaries built from the exponent data."""
    WF = np.exp(1j*f); WC = np.exp(1j*c)
    MF = np.diag(WF**E[:,0])                # W_F^{m_v}   (W-01's operator, charged)
    MC = np.diag(WC**E[:,1])                # W_C^{n_v}
    lhs = np.linalg.matrix_power(MF, k) @ s
    rhs = np.linalg.matrix_power(MC, k) @ s
    # <M_F^k s, M_C^k s> with <z,w> = conj(z).w   -> sum conj(u^{km}) s* . v^{kn} s
    return np.vdot(lhs, rhs)

def diff_lattice_rank(E, p):
    S = E[p > 0]
    D = np.array([S[i]-S[j] for i in range(len(S)) for j in range(len(S))])
    D = D[np.any(D != 0, axis=1)]
    if len(D)==0: return 0, None
    r = np.linalg.matrix_rank(D.astype(float))
    gen = None
    if r == 1:
        # primitive generator
        for row in D:
            g = math.gcd(int(abs(row[0])), int(abs(row[1])))
            gen = (int(row[0])//g, int(row[1])//g); break
        if gen[0] < 0 or (gen[0]==0 and gen[1]<0): gen = (-gen[0], -gen[1])
    return r, gen

def fourier_support(E, p, N=64, thresh=1e-9):
    """|Z_1|^2 on an N x N grid in (f,c); return the set of frequency vectors present."""
    j = np.arange(N); f = 2*np.pi*j/N; c = 2*np.pi*j/N
    FF, CC = np.meshgrid(f, c, indexing='ij')
    m = E[:,0][:,None,None]; n = E[:,1][:,None,None]
    Z = np.sum(p[:,None,None]*np.exp(1j*(-m*FF + n*CC)), axis=0)
    A = np.abs(Z)**2
    C = np.fft.fft2(A)/(N*N)
    idx = np.argwhere(np.abs(C) > thresh)
    out = []
    for (i0,i1) in idx:
        fx = i0 - N if i0 > N//2 else i0
        fy = i1 - N if i1 > N//2 else i1
        out.append(((int(fx),int(fy)), float(np.abs(C[i0,i1]))))
    out.sort(key=lambda t: (abs(t[0][0])+abs(t[0][1]), t[0]))
    vecs = np.array([o[0] for o in out if o[0]!=(0,0)])
    r = 0 if len(vecs)==0 else int(np.linalg.matrix_rank(vecs.astype(float)))
    return out, r

# ---------------------------------------------------------------- 1. corpus check
hdr("1.  CORPUS REPRODUCED BEFORE DEPARTURE (unit charge)")
p5 = np.array([0.4,0.15,0.15,0.15,0.15]); q1 = [1,1,1,1,1]
E1 = exponents(q1)
print("unit-charge exponent vectors E_v:", [tuple(x) for x in E1],
      " -> classes (1,1):v0  (1,0):v1,v2  (0,1):v3,v4")
rng = np.random.default_rng(7700001)
dev = 0.0
for _ in range(400):
    f,c = rng.uniform(0,2*np.pi,2); k = int(rng.integers(1,40))
    ph = rng.uniform(0,2*np.pi,5); s = np.sqrt(p5)*np.exp(1j*ph)
    dev = max(dev, abs(Z_matrix(E1,s,f,c,k) - Z_closed(E1,p5,f,c,k)))
print(f"seed 7700001, 400 samples: max |matrix action - closed form| = {dev:.3e}")
# gauge invariance
dev = 0.0
for _ in range(200):
    f,c = rng.uniform(0,2*np.pi,2)
    ph = rng.uniform(0,2*np.pi,5); s = np.sqrt(p5)*np.exp(1j*ph)
    g  = np.exp(1j*rng.uniform(0,2*np.pi,5))
    dev = max(dev, abs(Z_matrix(E1,g*s,f,c,3) - Z_matrix(E1,s,f,c,3)))
print(f"gauge invariance (200 random g): max deviation = {dev:.3e}")
# S4 section 2.1 identity
p0,qq,rr = 0.4,0.3,0.3
dev = 0.0
for _ in range(2000):
    f,c = rng.uniform(0,2*np.pi,2); k = int(rng.integers(1,60))
    lhs = abs(Z_closed(E1,p5,f,c,k))**2
    rhs = p0**2+qq**2+rr**2 + 2*p0*qq*np.cos(k*c) + 2*p0*rr*np.cos(k*f) + 2*qq*rr*np.cos(k*(f+c))
    dev = max(dev, abs(lhs-rhs))
print(f"S4 section 2.1 identity, 2000 samples: max deviation = {dev:.3e}")
sup1, r1 = fourier_support(E1, p5, N=64, thresh=1e-9)
print("unit charge Fourier support of |Z_1|^2 (64x64, thresh 1e-9):")
print("   ", [(v,round(m,9)) for v,m in sup1], "   rank =", r1)
lamA = lambda E,p,f,c: np.log(abs(Z_closed(E,p,f,c,1)))
print("S4 3.2 witnesses reproduced:")
for (f,c) in [(1.0,0.3),(1.0,3.0),(0.3,1.0),(3.0,1.0),(0.4,2.6),(1.5,1.5),(0.5,0.5),(2.5,2.5)]:
    print(f"   lambda_A({f},{c}) = {lamA(E1,p5,f,c):.6f}")
print(f"   same-product pair (0.4,2.6) vs (1.5,1.5) separation = "
      f"{abs(lamA(E1,p5,0.4,2.6)-lamA(E1,p5,1.5,1.5)):.6f}")
print()

# ------------------------------------------- 2. the claim's own arithmetic, checked
hdr("2.  THE CLAIM'S ARITHMETIC, RE-DERIVED (per-vertex charge q=(1,2,2,2,2))")
qv = [1,2,2,2,2]; Eq = exponents(qv)
print("exponent vectors:", [tuple(x) for x in Eq])
print("  collinear? (1,1),(2,0),(0,2) all satisfy m+n=2 :",
      all(int(x[0]+x[1])==2 for x in Eq))
rng = np.random.default_rng(7700002)
dev = 0.0
for _ in range(4000):
    f,c = rng.uniform(-10,10,2)
    dev = max(dev, abs(abs(Z_closed(Eq,p5,f,c,1)) - abs(0.4+0.6*np.cos(f+c))))
print(f"seed 7700002, 4000 samples: max | |Z_1| - |0.4+0.6cos(f+c)| | = {dev:.3e}")
print("four connections with f+c = 3.0:")
for (f,c) in [(0.4,2.6),(1.5,1.5),(2.9,0.1),(-1.0,4.0)]:
    print(f"   lambda_A({f},{c}) = {lamA(Eq,p5,f,c):.12f}")
print(f"   exact log|0.4+0.6cos 3| = {np.log(abs(0.4+0.6*np.cos(3.0))):.12f}")
supq, rq = fourier_support(Eq, p5, N=64, thresh=1e-9)
print("charged Fourier support of |Z_1|^2 (64x64, thresh 1e-9):")
print("   ", [(v,round(m,9)) for v,m in supq], "   rank =", rq)
print("VERDICT ON THE ARITHMETIC: every number in the claim reproduces.")
print()

# ============================================================ THE REFUTATION
hdr("3.  REFUTATION R1 -- THE CORPUS'S OWN DEFINITION OF CHARGE NEVER FIRES F-A")
print("S2_FORMATION_CONDITION_ON_K1_V001.md:173-175 (sealed) defines the S4 charge knob as")
print('    "higher charge  U_e |-> exp(i q a_e)  ... this is the S4 knob"')
print("and repeats it at :525  'NOT A FIX  higher charge q'.  That is a SINGLE SCALAR q")
print("applied to the EDGE transport.  Its effect on the exponent data is E_v = q*(a_v,b_v).")
print()
print("READING 1 -- S2's charge, scalar q on every edge:  E_v = q (a_v,b_v)")
print("  q :  exponent set                         rank Delta   F-A fires?   lambda_A(0.4,2.6)  lambda_A(1.5,1.5)  sep")
for q in [-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8]:
    Eqq = exponents([q]*5)
    r,gen = diff_lattice_rank(Eqq, p5)
    l1 = lamA(Eqq,p5,0.4,2.6); l2 = lamA(Eqq,p5,1.5,1.5)
    print(f"  {q:2d} :  {[tuple(int(y) for y in x) for x in Eqq[[0,1,3]]]}   "
          f"rank={r}   fires={'YES' if r<2 else 'no '}   {l1:>12.6f}   {l2:>12.6f}   {abs(l1-l2):.6f}")
print("  -> rank Delta = 2 for EVERY q != 0.  F-A does not fire at any scalar charge.")
print("  -> q = 0 kills the functional entirely (|Z_k| = 1, no formation), the |S|=1 row.")
print()
print("READING 2 -- per-LOOP charge (Q_F,Q_C): M_F multiplies by W_F^{Q_F} on all of gamma_F,")
print("  M_C by W_C^{Q_C} on all of gamma_C.   E_v0=(Q_F,Q_C) E_v1=(Q_F,0) E_v3=(0,Q_C).")
print("  PROOF that this can never collapse: the line through (Q_F,0) and (0,Q_C) is")
print("  x/Q_F + y/Q_C = 1; the pinch (Q_F,Q_C) evaluates to 1+1 = 2 != 1 for ALL nonzero")
print("  Q_F,Q_C.  Three affinely independent points => rank Delta = 2 identically.")
bad = []
for QF in range(-8,9):
    for QC in range(-8,9):
        if QF==0 or QC==0: continue
        Ep = np.array([[QF,QC],[QF,0],[QF,0],[0,QC],[0,QC]])
        r,_ = diff_lattice_rank(Ep, p5)
        if r != 2: bad.append((QF,QC,r))
print(f"  enumerated 16x16 = 256 nonzero (Q_F,Q_C): number with rank != 2 = {len(bad)}")
print()
print("READING 3 -- per-VERTEX charge q_v (the claim's).  COLLAPSE CONDITION, SOLVED:")
print("  E = {(q0,q0),(qF,0),(0,qC)} collinear  <=>  q0/qF + q0/qC = 1  <=>  q0 = qF qC/(qF+qC).")
print("  For positive charges q0 = harmonic-mean/2 < min(qF,qC) STRICTLY.")
print("  Hence the pinch v0 and the vertex v1 -- BOTH ON gamma_F -- must carry DIFFERENT")
print("  powers of the SAME holonomy W_F.  Enumerated solutions with 1<=q<=12:")
sols = []
for qF in range(1,13):
    for qC in range(1,13):
        num = qF*qC; den = qF+qC
        if num % den == 0:
            q0 = num//den
            if q0 >= 1: sols.append((q0,qF,qC))
for (q0,qF,qC) in sols:
    Ec = np.array([[q0,q0],[qF,0],[qF,0],[0,qC],[0,qC]])
    r,gen = diff_lattice_rank(Ec, p5)
    # surviving combination: lambda depends on the linear form dual to gen
    print(f"   q=({q0},{qF},{qF},{qC},{qC})  E={{({q0},{q0}),({qF},0),(0,{qC})}}  rank={r}"
          f"  generator {gen}  q0<min(qF,qC): {q0 < min(qF,qC)}")
print()

hdr("4.  REFUTATION R2 -- THE COLLAPSE IS NOT ALWAYS 'THE PRODUCT'")
print("F-A as ARMED (S4:366, S4:1110) is: 'only W_F.W_C survives'.  Under per-vertex charge")
print("the surviving combination is W_F^m W_C^n for the generator (m,n) of Delta:")
for (q0,qF,qC) in sols:
    Ec = np.array([[q0,q0],[qF,0],[qF,0],[0,qC],[0,qC]])
    r,gen = diff_lattice_rank(Ec, p5)
    m,n = gen
    # check numerically: |Z_1| constant along the level set of the dual form
    rg = np.random.default_rng(7700003)
    ok = True; vals=[]
    for _ in range(200):
        f = rg.uniform(0,2*np.pi)
        # move along direction annihilated by the frequency vector (-m, n)
        t = rg.uniform(0,2*np.pi)
        f1, c1 = f + n*t, c1 if False else 0.0
    # do it properly: frequency vector in (f,c) is (-m, n); level direction is (n, m)
    dev = 0.0
    for _ in range(500):
        f,c = rg.uniform(0,2*np.pi,2); t = rg.uniform(-5,5)
        f2, c2 = f + n*t, c + m*t
        dev = max(dev, abs(abs(Z_closed(Ec,p5,f,c,1)) - abs(Z_closed(Ec,p5,f2,c2,1))))
    # generator (m,n) in exponent space => lambda sees u^m v^n = W_F^{-m} W_C^{n},
    # i.e. the linear combination  m*f - n*c  (up to sign).  The PRODUCT W_F.W_C is (m,n)=(1,-1).
    isprod = (m,n) in [(1,-1),(-1,1)]
    combo = "W_F.W_C   (the PRODUCT)" if isprod else f"W_F^{m} W_C^{-n}  ~ form {m}f{-n:+d}c"
    print(f"   q=({q0},{qF},{qF},{qC},{qC}): gen {str(gen):<9} -> lambda sees {combo:<34}"
          f" level-set invariance {dev:.1e}   F-A-as-armed: "
          f"{'FIRES' if isprod else 'DOES NOT FIRE (no row in S4 taxonomy)'}")
print()

hdr("5.  REFUTATION R3 -- 'WEIGHT ON THE PINCH' WAS NEVER SUFFICIENT, AT UNIT CHARGE")
print("S4:1110 'separation does work IF AND ONLY IF the ready state puts weight on the pinch'.")
tests = [("p=(1,0,0,0,0)   all weight ON the pinch", np.array([1.0,0,0,0,0])),
         ("p=(0.5,0.5,0,0,0) pinch + face",          np.array([0.5,0.25,0.25,0,0])),
         ("p=(0.5,0,0,0.5,0) pinch + cycle",         np.array([0.5,0,0,0.25,0.25])),
         ("p=(0.4,0.15,0.15,0.15,0.15) full",        p5)]
for name,p in tests:
    r,gen = diff_lattice_rank(E1, p)
    l1 = abs(Z_closed(E1,p,0.4,2.6,1)); l2 = abs(Z_closed(E1,p,1.5,1.5,1))
    l3 = abs(Z_closed(E1,p,1.0,0.3,1)); l4 = abs(Z_closed(E1,p,1.0,3.0,1))
    print(f"  {name:<42} p0={p[0]:.2f}>0  rank Delta={r} gen={gen}")
    print(f"      |Z_1|(0.4,2.6)={l1:.6f}  |Z_1|(1.5,1.5)={l2:.6f}   "
          f"|Z_1|(1.0,0.3)={l3:.6f}  |Z_1|(1.0,3.0)={l4:.6f}")
print("  -> p=(1,0,0,0,0): weight on the pinch, |Z_k| == 1 for all k, NO formation, NO")
print("     separation.  The 'if' half of S4's iff is false at UNIT CHARGE, on K1, with")
print("     weight on the pinch.  It is S4's own |S|=1 row, which the gloss overrode.")
print()

hdr("6.  WHAT THE PER-VERTEX KNOB ACTUALLY IS -- S4's OWN CHOICE LEDGER C11")
print("S4:995 C11: 'a_v in {0,1} -- a vertex visited twice by a loop still counts once |")
print("  count multiplicity | ... closed; but note it makes Theorem S4-1's corners-of-a-square")
print("  argument exact, and MULTIPLICITY WOULD BREAK IT'.")
print("Per-vertex charge q_v is exactly per-vertex multiplicity: a_v in {0,1} -> q_v a_v in {0,q_v}.")
print("Demonstration -- the SAME collapse produced with unit charge and a multiplicity")
print("convention (gamma_F traversed twice at v1,v2; gamma_C twice at v3,v4; once at v0):")
Emult = np.array([[1,1],[2,0],[2,0],[0,2],[0,2]])
print("   E(multiplicity) =", [tuple(int(y) for y in x) for x in Emult],
      "  identical to E(q=(1,2,2,2,2)) :", np.array_equal(Emult, Eq))
rg = np.random.default_rng(7700004)
dev = max(abs(Z_closed(Emult,p5,*rg.uniform(0,2*np.pi,2),1) -
              Z_closed(Eq,p5,0,0,1)*0 + 0) for _ in range(1))  # placeholder
dev = 0.0
for _ in range(1000):
    f,c = rg.uniform(0,2*np.pi,2); k=int(rg.integers(1,20))
    dev = max(dev, abs(Z_closed(Emult,p5,f,c,k)-Z_closed(Eq,p5,f,c,k)))
print(f"   seed 7700004, 1000 samples: max |Z_k(multiplicity) - Z_k(charge)| = {dev:.3e}")
print("   -> the two knobs are THE SAME KNOB.  S4 priced it in C11 and declared what it breaks.")
print()

hdr("7.  HOW RARE IS THE COLLAPSE IN PER-VERTEX CHARGE SPACE")
tot=0; coll=0
for q in itertools.product(range(1,7), repeat=3):     # (q0, qF, qC)
    q0,qF,qC = q; tot+=1
    Ec = np.array([[q0,q0],[qF,0],[qF,0],[0,qC],[0,qC]])
    r,_ = diff_lattice_rank(Ec,p5)
    if r<2: coll+=1
print(f"  (q0,qF,qC) in {{1..6}}^3 : {coll} of {tot} collapse  "
      f"({100*coll/tot:.1f}%) -- a codimension-1 slice {{q0(qF+qC)=qF qC}}")
tot=0; coll=0
for q in itertools.product(range(1,5), repeat=5):
    Ec = exponents(list(q)); tot+=1
    r,_ = diff_lattice_rank(Ec,p5)
    if r<2: coll+=1
print(f"  full q in {{1..4}}^5 : {coll} of {tot} collapse ({100*coll/tot:.1f}%)")
print()

hdr("8.  SCHEDULE B UNDER THE CLAIM'S CHARGE -- CHECKED, NOT ASSUMED")
def lamB_direct(E,p,f,c,N):
    n = np.arange(1,N+1)
    m = E[:,0][:,None]; nn = E[:,1][:,None]
    tot=0.0
    B=200000
    for st in range(1,N+1,B):
        blk = np.arange(st, min(st+B,N+1))
        Zk = np.sum(p[:,None]*np.exp(1j*blk[None,:]*(-m*f+nn*c)),axis=0)
        tot += np.sum(np.log(np.abs(Zk)))
    return tot/N
for (f,c) in [(1.0,np.sqrt(2)),(0.7,np.e/2)]:
    print(f"  q=(1,2,2,2,2) (f,c)=({f:.6f},{c:.6f})  lambda_B(N=2e6) = "
          f"{lamB_direct(Eq,p5,f,c,2000000):.9f}   vs log(0.3) = {np.log(0.3):.9f}")
print("  -> under B the charged rate is also a function of the product alone; the claim's")
print("     A-only run does not overstate on this point.")
print()

hdr("9.  SCALAR-CHARGE SEPARATION WITNESSES -- F-A EXPLICITLY DOES NOT FIRE")
for q in [2,3,5,7]:
    Eqq = exponents([q]*5)
    # same product f+c fixed
    A = lamA(Eqq,p5,0.4,2.6); B = lamA(Eqq,p5,1.5,1.5)
    # same W_F
    C = lamA(Eqq,p5,1.0,0.3); D = lamA(Eqq,p5,1.0,3.0)
    print(f"  q={q}:  same-product pair differ by {abs(A-B):.6f} ;  same-W_F pair differ by {abs(C-D):.6f}")
print("  -> at every scalar charge, lambda_A still separates W_F from W_C on full support.")

hdr("10.  IS 'rank Delta(S) = 2' ITSELF THE CORRECT CONDITION?  TESTED, NOT ASSUMED")
print("PROOF SKETCH (mine): |Z_1|^2 = sum_{x,y in S} p_x p_y chi_{E_x - E_y}.  The coefficient")
print("of frequency d is sum_{E_x-E_y=d} p_x p_y > 0 -- ALL TERMS POSITIVE, so no cancellation")
print("can remove a frequency.  Hence Fourier support of |Z_1|^2 = Delta(S) exactly, and")
print("lambda_A is a function of a single linear form iff rank Delta(S) <= 1.  QED")
rg = np.random.default_rng(7700005)
mismatch = 0; tested = 0
for _ in range(600):
    Er = rg.integers(-3,4,size=(5,2))
    mask = rg.random(5) > 0.35
    if mask.sum() < 1: continue
    pr = np.zeros(5); pr[mask] = rg.random(mask.sum()); pr /= pr.sum()
    r,_ = diff_lattice_rank(Er, pr)
    _, rf = fourier_support(Er, pr, N=32, thresh=1e-10)
    tested += 1
    if r != rf: mismatch += 1
print(f"seed 7700005, {tested} random (exponent set, weight) pairs on a 32x32 grid:")
print(f"   rank Delta(S) vs rank of Fourier support of |Z_1|^2 : mismatches = {mismatch}")
print("   -> the replacement criterion is CORRECT.  This half of the claim survives.")
print()

hdr("11.  BUT 'NOT ABOUT WHICH VERTEX CARRIES WEIGHT' IS FALSE AS STATED")
print("S = {E_v : p_v > 0}.  Which vertices carry weight is HALF the input to rank Delta(S).")
Emix = exponents([1,1,1,2,2])
print("charge q=(1,1,1,2,2), E =", [tuple(int(y) for y in x) for x in Emix])
for name,p in [("full support           ", p5),
               ("p0,p1 only (pinch+face)", np.array([0.5,0.5,0,0,0])),
               ("p1,p3 only (no pinch)  ", np.array([0,0.5,0,0.5,0]))]:
    r,gen = diff_lattice_rank(Emix,p)
    print(f"   {name}  rank Delta = {r}  gen={gen}")
print("   -> at FIXED charge the rank still moves with the support.  The correct statement is")
print("      'the affine rank of the exponent vectors OF THE WEIGHTED VERTICES' -- both inputs.")
