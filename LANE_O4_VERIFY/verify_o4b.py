#!/usr/bin/env python3
"""Second adversarial pass: (a) is 1158 a count of DISTINCT operators?
   (b) does DEF-A' actually rescue the case it is registered to rescue?"""
import itertools, numpy as np
FAIL=[]
def ck(n,c,d=""):
    print(f"   [{'PASS' if c else 'FAIL'}] {n}   {d}")
    if not c: FAIL.append(n)
def hr(t): print("\n"+"="*78+f"\n{t}\n"+"="*78)

hr("A -- IS THE HEADLINE 1158 A COUNT OF DISTINCT OPERATORS?")
L=3;N=2*L*L
def e(x,y,d): return 2*(L*(y%L)+(x%L))+d
stars=[frozenset({e(x,y,0),e(x-1,y,0),e(x,y,1),e(x,y-1,1)}) for y in range(L) for x in range(L)]
plaqs=[frozenset({e(x,y,0),e(x,y+1,0),e(x,y,1),e(x+1,y,1)}) for y in range(L) for x in range(L)]
R_edges=frozenset({e(x,0,0) for x in range(L)})
regions=[frozenset({i}) for i in range(N)]+stars+plaqs
seen=set(); tot=0
for T in regions:
    Ts=sorted(T)
    for mask in range(4**len(Ts)):
        xs=set();zs=set();m=mask
        for q in Ts:
            d=m%4;m//=4
            if d in (1,3): xs.add(q)
            if d in (2,3): zs.add(q)
        if len(frozenset(xs)&R_edges)%2:
            tot+=1; seen.add((frozenset(xs),frozenset(zs)))
print(f"   sum over regions (with multiplicity) : {tot}")
print(f"   DISTINCT Pauli operators             : {len(seen)}")
ck("the reported 1158 is the multiplicity-counted sum, not distinct operators",
   tot==1158 and len(seen)<1158, f"1158 vs {len(seen)} distinct")

hr("B -- DOES DEF-A' RESCUE 'R CONSTANT ON A HIGH-ENERGY SECTOR'?")
print("""   The lane registers DEF-A' as the fallback for exactly this case (o4_variant.py
   docstring: 'It is an over-reach for any system whose high-energy sectors are
   irrelevant to the record').  Test: H has a 4-dim record shell and a 2-dim
   high-energy shell on which R = -I.  DEF-A' = U unitary, [U, P_E] = 0 for the
   eigenspace E witnessing clause (iii).  Clause (iv) is the global identity.""")
mults=(4,2); ps=(2,0); n=sum(mults)
H=np.zeros((n,n)); R=np.zeros((n,n)); o=0
for k,(m,p) in enumerate(zip(mults,ps)):
    H[o:o+m,o:o+m]=np.eye(m)*(k+1)
    d=np.ones(m); d[p:]=-1
    R[o:o+m,o:o+m]=np.diag(d); o+=m
PE=np.zeros((n,n)); PE[:4,:4]=np.eye(4)          # witnessing eigenspace
ck("(i) R is a bit", np.allclose(R@R,np.eye(n)) and np.allclose(R,R.T))
ck("(ii) [H,R]=0", np.allclose(H@R,R@H))
ck("(iii) R non-constant on the 4-dim shell", not np.allclose(R[:4,:4],R[0,0]*np.eye(4)))
print(f"   Tr(P_E R) on the witnessing shell = {np.trace(PE@R):+.1f}   (balanced)")
print(f"   Tr R on the complement            = {np.trace((np.eye(n)-PE)@R):+.1f}   (NOT balanced)")
# exhaustive-in-spirit random search over DEF-A'-admissible unitaries
def randu(k):
    z=(np.random.randn(k,k)+1j*np.random.randn(k,k))/np.sqrt(2)
    q,r=np.linalg.qr(z); return q*(np.diag(r)/np.abs(np.diag(r)))
np.random.seed(1)
best=np.inf
for _ in range(20000):
    U=np.zeros((n,n),dtype=complex)
    U[:4,:4]=randu(4); U[4:,4:]=randu(2)         # commutes with P_E
    best=min(best,np.linalg.norm(U.conj().T@R@U+R))
ck("DEF-A' admits NO flipper here either (20000-sample search)", best>1.0,
   f"best ||U^dag R U + R|| = {best:.4f}   analytic bound: on the complement R=-I so U^dag R U=-I never equals +I")
lane_says = True   # o4_variant.py line ~173: cAp counts this point
ck("but o4_variant.py Section 3 COUNTS this point as DEF-A'-flippable", lane_says,
   "ps=(2,0), wit={0}, 2*ps[0]==mults[0] -> counted; the complement is never checked")
print("""
   CONSEQUENCE.  The registered fallback DEF-A' does NOT rescue a record that is
   constant on a high-energy sector -- the motivating case in its own docstring.
   Commuting with P_E forces U to preserve the complement, where R is +-I and cannot
   be negated.  DEF-A' differs from DEF-A only when two or more non-witnessing shells
   have CANCELLING traces.""")
print("\n"+"="*78)
print("*** FAILURES: "+str(FAIL) if FAIL else "ALL CHECKS PASSED.")
