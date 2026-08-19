"""T-11 / PF-3B: does clause (v) under DEF-A hold on a SECOND carrier?

PF-3B could not run: its bouquet had X-cosystole 1, so d=1 and the carrier fails clause (v) on
its own merits. T-8 built one that qualifies -- torus 2x2 + a disk capping the wrap {0,2},
[[8,1,2]], NON-MANIFOLD, both distances 2, same dimension.

CONTRACTIBLE REGION on a graph = an edge subset containing NO cycle, i.e. a FOREST. The test is
the one O-4 ran on the incumbent: count operators supported in such a region that FLIP the
record, then count how many of those are ADMISSIBLE ([U,H] = 0). A zero with no non-zero control
beside it is not a measurement."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); Y=1j*(X@Z)
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
           .read().split('say("="*104); say("0.')[0],g)
op=g['op']; L=g['L']; EDGES=g['EDGES']; V=g['nx']*g['ny']
H_A=g['H0']; capZ=op({0:Z,2:Z},L); H_B=H_A-capZ
P4=[I2,X,Y,Z]
def has_cycle(sub):
    par=list(range(V))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for l in sub:
        a,b=EDGES[l]; ra,rb=find(a),find(b)
        if ra==rb: return True
        par[ra]=rb
    return False
def connected(sub):
    if not sub: return False
    vs={x for l in sub for x in EDGES[l]}
    seen={next(iter(vs))}; changed=True
    while changed:
        changed=False
        for l in sub:
            a,b=EDGES[l]
            if (a in seen) ^ (b in seen): seen.add(a); seen.add(b); changed=True
    return seen==vs
# CLAUSE (v) SAYS "A SINGLE CONTRACTIBLE REGION". A scattered forest is contractible but is NOT
# ONE region -- it is several. O-4 measured 0 admissible flippers using geometric disks; using
# arbitrary forests instead gives 10, and the difference is entirely the word SINGLE.
allf=[list(s) for k in range(1,L+1) for s in itertools.combinations(range(L),k) if not has_cycle(s)]
regions=[r for r in allf if connected(r)]
scattered=[r for r in allf if not connected(r)]
say("="*98); say("T-11 / PF-3B   CLAUSE (v) UNDER DEF-A ON A SECOND CARRIER"); say("="*98)
say(f"  forests in total          : {len(allf)}")
say(f"  CONNECTED (a SINGLE region): {len(regions)}, largest {max(len(r) for r in regions)} edges")
say(f"  scattered (several regions): {len(scattered)}  -- excluded by the word SINGLE in clause (v)")
say("")
say(f"  {'carrier':<36}{'record':<8}{'SINGLE: any':>14}{'SINGLE: adm':>14}{'scatt: any':>13}{'scatt: adm':>13}")
for nm,H,R,rn in (("A  toric 2x2  [[8,2,2]]",H_A,g['Zbar'],"Zbar"),
                  ("B  torus+cap  [[8,1,2]] NON-MANIFOLD",H_B,g['Zbar2'],"Zbar2")):
    c1=np.linalg.norm(R@H-H@R)
    anyf=admf=0; sanyf=sadmf=0
    for T in regions:
        for combo in itertools.product(range(4),repeat=len(T)):
            if all(c==0 for c in combo): continue
            A=op({l:P4[c] for l,c in zip(T,combo)},L)
            if np.linalg.norm(A.conj().T@R@A + R) < 1e-8:
                anyf+=1
                if np.linalg.norm(A@H-H@A) < 1e-8: admf+=1
    for T in scattered:
        for combo in itertools.product(range(4),repeat=len(T)):
            if all(c==0 for c in combo): continue
            A=op({l:P4[c] for l,c in zip(T,combo)},L)
            if np.linalg.norm(A.conj().T@R@A + R) < 1e-8:
                sanyf+=1
                if np.linalg.norm(A@H-H@A) < 1e-8: sadmf+=1
    say(f"  {nm:<36}{rn:<8}{anyf:>14}{admf:>14}{sanyf:>13}{sadmf:>13}")
say("")
say("  A zero in the right-hand column is a measurement only because the left-hand column is not")
say("  zero: the same enumeration, with admissibility dropped, finds flippers.")
