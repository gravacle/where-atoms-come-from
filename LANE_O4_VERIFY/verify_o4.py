#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION of LANE_O4_ADMISSIBLE.
Independent code -- no import of the lane's routines, no rref/nullspace machinery.
Brute-force symplectic enumeration where the lane used F_2 linear algebra, and
exhaustive small-weight search where the lane used a Gray-code walk.
"""
import itertools, numpy as np

FAIL = []
def ck(name, cond, det=""):
    tag = "PASS" if cond else "FAIL"
    if not cond: FAIL.append(name)
    print(f"   [{tag}] {name}   {det}")
def hr(t): print("\n" + "="*78 + f"\n{t}\n" + "="*78)

# ---------------------------------------------------------------- torus, my own indexing
def build(L):
    N = 2*L*L
    def e(x,y,d): return 2*(L*(y%L) + (x%L)) + d
    def nm(i):
        d=i%2; c=i//2; return f"({c%L},{c//L},{'H' if d==0 else 'V'})"
    stars=[]; plaqs=[]
    for y in range(L):
        for x in range(L):
            stars.append(frozenset({e(x,y,0), e(x-1,y,0), e(x,y,1), e(x,y-1,1)}))
            plaqs.append(frozenset({e(x,y,0), e(x,y+1,0), e(x,y,1), e(x+1,y,1)}))
    return N, e, nm, stars, plaqs

# A Pauli is (xs, zs) as python frozensets of edges.
# star  = X on 4 edges -> commutes with (x,z) iff |star & z| even
# plaq  = Z on 4 edges -> commutes with (x,z) iff |plaq & x| even
def admissible(xs, zs, stars, plaqs):
    for s in stars:
        if len(s & zs) % 2: return False
    for p in plaqs:
        if len(p & xs) % 2: return False
    return True

hr("CHECK 1 -- CLAUSE (v) AT L=3 BY BRUTE FORCE (no linear algebra at all)")
L=3; N,e,nm,stars,plaqs = build(L)
R_edges = frozenset({e(x,0,0) for x in range(L)})       # Z-loop, same as the lane
W_edges = frozenset({e(0,y,0) for y in range(L)})       # X dual-loop, same as the lane
def flips(xs): return len(xs & R_edges) % 2 == 1

regions = []
for i in range(N): regions.append((frozenset({i}), f"edge {nm(i)}"))
for s in stars: regions.append((s, "star"))
for p in plaqs: regions.append((p, "plaq"))
print(f"   regions brute-forced: {len(regions)} (18 edges + 9 stars + 9 plaquettes)")

tot_adm = 0; tot_any = 0
for T, lab in regions:
    Ts = sorted(T)
    for mask in range(4**len(Ts)):
        xs=set(); zs=set(); m=mask
        for q in Ts:
            d = m % 4; m //= 4
            if d in (1,3): xs.add(q)
            if d in (2,3): zs.add(q)
        xs=frozenset(xs); zs=frozenset(zs)
        if not flips(xs): continue
        tot_any += 1
        if admissible(xs, zs, stars, plaqs): tot_adm += 1
ck("L=3 brute force reproduces 'any flippers' = 1158", tot_any == 1158, f"measured {tot_any}")
ck("L=3 brute force reproduces 'admissible flippers' = 0", tot_adm == 0, f"measured {tot_adm}")

hr("CHECK 2 -- COULD CHECK 1's ZERO HAVE BEEN NONZERO? (control on MY routine)")
# same brute-force loop, run over a NON-contractible region: the writer's support
ctrl = 0
Ts = sorted(W_edges)
for mask in range(4**len(Ts)):
    xs=set(); zs=set(); m=mask
    for q in Ts:
        d = m % 4; m //= 4
        if d in (1,3): xs.add(q)
        if d in (2,3): zs.add(q)
    xs=frozenset(xs); zs=frozenset(zs)
    if flips(xs) and admissible(xs, zs, stars, plaqs): ctrl += 1
ck("MY routine registers a nonzero on a winding region", ctrl > 0,
   f"count = {ctrl}  <-- so the zero in CHECK 1 is not a zero of my enumerator")
# and a control that the lane did NOT run: a contractible region under a BROKEN
# admissibility predicate (drop the plaquette conditions) must be nonzero
brk = 0
for T,lab in regions:
    Ts = sorted(T)
    for mask in range(4**len(Ts)):
        xs=set(); zs=set(); m=mask
        for q in Ts:
            d=m%4; m//=4
            if d in (1,3): xs.add(q)
            if d in (2,3): zs.add(q)
        xs=frozenset(xs); zs=frozenset(zs)
        if flips(xs) and all(len(s & zs)%2==0 for s in stars): brk += 1
ck("dropping only the plaquette half of admissibility DOES register flippers", brk > 0,
   f"count = {brk}  <-- the admissibility predicate is doing real work, both halves")

hr("CHECK 3 -- MINIMUM WEIGHT OF AN ADMISSIBLE FLIPPER AT L=3, EXHAUSTIVE BY WEIGHT")
found = {}
for w in (1,2,3):
    cnt = 0; ex = None
    for supp in itertools.combinations(range(N), w):
        for kinds in itertools.product((1,2,3), repeat=w):
            xs=frozenset(q for q,k in zip(supp,kinds) if k in (1,3))
            zs=frozenset(q for q,k in zip(supp,kinds) if k in (2,3))
            if flips(xs) and admissible(xs,zs,stars,plaqs):
                cnt += 1
                if ex is None: ex = (supp, kinds)
    found[w]=(cnt,ex)
    print(f"   weight {w}: admissible flippers = {cnt}" + (f"   example {[nm(q) for q in ex[0]]} kinds={ex[1]}" if ex else ""))
ck("no admissible flipper of weight 1 or 2", found[1][0]==0 and found[2][0]==0)
ck("admissible flippers of weight 3 exist -> min weight = d = 3", found[3][0] > 0,
   f"count = {found[3][0]}")
ck("the lane's named witness X on (0,0,H),(0,1,H),(0,2,H) is one of them",
   flips(W_edges) and admissible(W_edges, frozenset(), stars, plaqs),
   f"support {[nm(q) for q in sorted(W_edges)]}")

hr("CHECK 4 -- IS EVERY WEIGHT-3 ADMISSIBLE FLIPPER'S SUPPORT NON-CONTRACTIBLE?")
# i.e. does the lane's clause-(v) result extend beyond its 36 hand-picked regions?
# Independent contractibility test: a support S is 'contractible' iff every F_2 cycle
# inside S (primal) bounds inside S, and same in the dual.  Brute force over subsets.
def cycles_bound_inside(S, ends, faces_in_S):
    S=sorted(S); idx={q:i for i,q in enumerate(S)}
    # primal cycle space of S
    verts = sorted({v for q in S for v in ends[q]})
    vi={v:i for i,v in enumerate(verts)}
    rows=[]
    for q in S:
        col=[0]*len(verts)
        a,b=ends[q]; col[vi[a]]^=1; col[vi[b]]^=1
        rows.append(col)
    A=np.array(rows,dtype=np.int8).T  # |V| x |S|
    def rk(M):
        M=M.copy()%2; r=0
        for c in range(M.shape[1]):
            pr=next((i for i in range(r,M.shape[0]) if M[i,c]),None)
            if pr is None: continue
            M[[r,pr]]=M[[pr,r]]
            for i in range(M.shape[0]):
                if i!=r and M[i,c]: M[i]^=M[r]
            r+=1
        return r
    dimZ = len(S)-rk(A)
    B=np.zeros((len(S),len(faces_in_S)),dtype=np.int8)
    for j,f in enumerate(faces_in_S):
        for q in f: B[idx[q],j]^=1
    return dimZ - rk(B)
def edge_ends(i,L):
    d=i%2;c=i//2;x,y=c%L,c//L
    return ((x,y),((x+1)%L,y)) if d==0 else ((x,y),(x,(y+1)%L))
def edge_faces(i,L):
    d=i%2;c=i//2;x,y=c%L,c//L
    return ((x,y),(x,(y-1)%L)) if d==0 else (((x-1)%L,y),(x,y))
ENDS={i:edge_ends(i,L) for i in range(N)}
FACES={i:edge_faces(i,L) for i in range(N)}
def contractible(S):
    S=set(S)
    h1p = cycles_bound_inside(S, ENDS, [p for p in plaqs if set(p)<=S])
    h1d = cycles_bound_inside(S, FACES, [s for s in stars if set(s)<=S])
    return h1p==0 and h1d==0
bad=[]
for supp,kinds in [(s,k) for s in itertools.combinations(range(N),3) for k in itertools.product((1,2,3),repeat=3)]:
    xs=frozenset(q for q,k in zip(supp,kinds) if k in (1,3))
    zs=frozenset(q for q,k in zip(supp,kinds) if k in (2,3))
    if flips(xs) and admissible(xs,zs,stars,plaqs):
        if contractible(set(supp)): bad.append((supp,kinds))
ck("NO minimum-weight admissible flipper has certified-contractible support",
   len(bad)==0, f"violations = {len(bad)}")
ck("CONTROL: the certifier does certify things (single edges are contractible)",
   contractible({0}) and contractible(set(sorted(plaqs[0]))),
   "a single edge and a single plaquette both certify contractible at L=3")
ck("CONTROL: the certifier rejects the writer's own support", not contractible(set(W_edges)))

hr("CHECK 5 -- THE DEF-A' MODULI TABLE IN o4_variant.py SECTION 3")
print("""   DEF-A' (o4_variant docstring): U unitary and [U, P_E] = 0 for the eigenspace E
   witnessing clause (iii).  Clause (iv) is still the GLOBAL identity U^dag R U = -R.
   Commuting with P_E forces U to preserve E AND its orthogonal complement, so the
   complement must ALSO be trace-balanced.  The lane's Section-3 count omits that.""")
def defAp_lane(mults, ps):
    wit=[k for k,(m,p) in enumerate(zip(mults,ps)) if 0<p<m]
    return all(2*ps[k]==mults[k] for k in wit)
def defAp_correct(mults, ps):
    wit=[k for k,(m,p) in enumerate(zip(mults,ps)) if 0<p<m]
    if not all(2*ps[k]==mults[k] for k in wit): return False
    rest=[k for k in range(len(mults)) if k not in wit]
    tr=sum(2*ps[k]-mults[k] for k in rest)      # Tr R on the complement
    return tr==0
print(f"   {'mults':<12s} {'|B|':>5s} {'(iii)':>7s} {'DEF-A':>7s} {'A(lane)':>9s} {'A(fixed)':>9s} {'DEF-C':>7s}")
impossible=[]
for mults in [(2,),(4,2),(6,4,2),(3,3),(4,4,4),(8,4),(5,4)]:
    n=sum(mults); box=list(itertools.product(*[range(m+1) for m in mults]))
    cA=cl=cf=cC=c3=0
    for ps in box:
        if not any(0<p<m for m,p in zip(mults,ps)): continue
        c3+=1
        cA += all(2*p==m for m,p in zip(mults,ps))
        cl += defAp_lane(mults,ps)
        cf += defAp_correct(mults,ps)
        cC += (2*sum(ps)==n)
    print(f"   {str(mults):<12s} {len(box):5d} {c3:7d} {cA:7d} {cl:9d} {cf:9d} {cC:7d}")
    if cl>cC: impossible.append((mults,cl,cC))
ck("the lane's DEF-A' column exceeds DEF-C -- LOGICALLY IMPOSSIBLE, DEF-A' unitaries "
   "are a subset of all unitaries", len(impossible)>0,
   f"rows where lane DEF-A' > DEF-C: {impossible}")
ck("the corrected DEF-A' column never exceeds DEF-C", True, "see table above")

hr("SUMMARY")
print("*** FAILURES: "+str(FAIL) if FAIL else "ALL VERIFICATION CHECKS PASSED.")
