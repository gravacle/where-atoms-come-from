# W-14 leg B — the state side is OPEN. Is it a CHANNEL, or only a KNOB?
# A constraint helps only if it lets the functional SEE something it could not see before.
# The functional sees the state only through pi (W-11's biconditional). So the test is:
# with pi PINNED EXACTLY, can any state freedom reach the functional?
#
# ISOLATION LEDGER. Held fixed: carrier, connection, observable, convention, k-range, code path,
# and pi BY CONSTRUCTION. Moved, one at a time: the within-class distribution; the phases.
# Run on B0b -- a FOUR-CLASS carrier -- so this is not a K1 fact.
import numpy as np
rng=np.random.default_rng(20260821)
def V(i,j): return 3*(j%3)+(i%3)
def H(i,j): return 3*(j%3)+(i%3)
def Wg(i,j): return 9+3*(j%3)+(i%3)
E=[None]*18
for j in range(3):
    for i in range(3):
        E[H(i,j)]=(V(i,j),V(i+1,j)); E[Wg(i,j)]=(V(i,j),V(i,j+1))
gF=[(H(0,0),1),(Wg(1,0),1),(H(0,1),-1),(Wg(0,0),-1)]; gC=[(H(0,0),1),(H(1,0),1),(H(2,0),1)]
def loopv(g):
    s=set()
    for (e,sg) in g: s.update(E[e])
    return s
VF,VC=loopv(gF),loopv(gC); NV=9
CLS=[(int(v in VF),int(v in VC)) for v in range(NV)]
def hol(g,a): return float(sum(s*a[e] for (e,s) in g))
def pi_of(w):
    p={(0,0):0.,(1,0):0.,(0,1):0.,(1,1):0.}
    for v in range(NV): p[CLS[v]]+=w[v]
    return np.array([p[(0,0)],p[(1,0)],p[(0,1)],p[(1,1)]])

a=rng.uniform(0,2*np.pi,18); WF,WC=np.exp(1j*hol(gF,a)),np.exp(1j*hol(gC,a))
def Zk(s,k):
    x=np.array([(WF**k if CLS[v][0] else 1)*s[v] for v in range(NV)])
    y=np.array([(WC**k if CLS[v][1] else 1)*s[v] for v in range(NV)])
    return abs(np.vdot(x,y))

base=np.array([.10,.12,.09,.14,.11,.11,.11,.11,.11]); base/=base.sum()
target=pi_of(base)
print(f"  carrier B0b (FOUR classes). class map {CLS}")
print(f"  pi PINNED at {np.round(target,12)}\n")
states=[]
for _ in range(60):                                   # states with EXACTLY this pi
    w=np.zeros(NV)
    for c in {(0,0),(1,0),(0,1),(1,1)}:
        idx=[v for v in range(NV) if CLS[v]==c]
        tot=target[[(0,0),(1,0),(0,1),(1,1)].index(c)]
        d=rng.dirichlet(np.ones(len(idx)))*tot
        for i,v in enumerate(idx): w[v]=d[i]
    ph=rng.uniform(0,2*np.pi,NV)
    states.append(np.sqrt(w)*np.exp(1j*ph))
assert max(np.abs(pi_of(np.abs(s)**2)-target).max() for s in states) < 1e-12
print("  60 states, identical pi, differing in BOTH within-class weight and phase:")
print(f"  {'k':>3} {'min |Z_k|':>18} {'max |Z_k|':>18} {'spread':>12}")
for k in (1,2,3,7,20):
    v=[Zk(s,k) for s in states]
    print(f"  {k:>3} {min(v):>18.14f} {max(v):>18.14f} {max(v)-min(v):>12.2e}")
print()
print("  ==> EXACTLY BLIND, on a four-class carrier. The state reaches the functional ONLY")
print("      through pi. So a singular constraint on the state is a KNOB -- it selects a pi --")
print("      and NOT a CHANNEL. It cannot make the construction sensitive to anything new.")
print("      The 'wrong side' possibility is therefore CLOSED, for a different reason than the")
print("      connection side: not because the rate cannot move, but because nothing else gets in.")
