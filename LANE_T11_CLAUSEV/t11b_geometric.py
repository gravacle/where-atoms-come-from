"""T-11b: 'CONTRACTIBLE REGION' HAS TWO READINGS AND THEY DISAGREE.

O-4 measured 0 admissible flippers at L=3 using CERTIFIED-CONTRACTIBLE PLAQUETTE BLOCKS -- disks
of the surface. t11_clausev used EDGE SUBGRAPHS with no cycle, and got 8 even after requiring
connectedness. The clause says 'no admissible operation on a single contractible region' and does
not say which. This measures both, on both carriers."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); Y=1j*(X@Z)
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
           .read().split('say("="*104); say("0.')[0],g)
op=g['op']; L=g['L']; PLAQ=g['PLAQ']; EDGES=g['EDGES']
H_A=g['H0']; capZ=op({0:Z,2:Z},L); H_B=H_A-capZ
P4=[I2,X,Y,Z]
# GEOMETRIC reading: a region is a set of plaquettes; its support is the edges they contain.
# A PROPER subset of the plaquettes on a torus is a disk; the full set wraps and is not.
blocks=[]
for k in range(1,len(PLAQ)):
    for s in itertools.combinations(range(len(PLAQ)),k):
        supp=sorted({e for p in s for e in PLAQ[p]})
        blocks.append((s,supp))
say("="*100); say("T-11b   'CONTRACTIBLE REGION' -- TWO READINGS, MEASURED"); say("="*100)
say(f"  GEOMETRIC regions (proper plaquette subsets): {len(blocks)}, support up to {max(len(s) for _,s in blocks)} edges")
say(f"  code distance d = 2 on both carriers")
say("")
say(f"  {'carrier':<36}{'record':<8}{'any-unitary':>14}{'ADMISSIBLE':>13}")
for nm,H,R,rn in (("A  toric 2x2  [[8,2,2]]",H_A,g['Zbar'],"Zbar"),
                  ("B  torus+cap [[8,1,2]] NON-MANIFOLD",H_B,g['Zbar2'],"Zbar2")):
    anyf=admf=0
    for _,supp in blocks:
        for combo in itertools.product(range(4),repeat=len(supp)):
            if all(c==0 for c in combo): continue
            A=op({l:P4[c] for l,c in zip(supp,combo)},L)
            if np.linalg.norm(A.conj().T@R@A + R) < 1e-8:
                anyf+=1
                if np.linalg.norm(A@H-H@A) < 1e-8: admf+=1
    say(f"  {nm:<36}{rn:<8}{anyf:>14}{admf:>13}")
say("")
say("  READ: if the geometric reading also gives a non-zero admissible count, then clause (v) is")
say("  FALSE ON THE INCUMBENT AT THIS SIZE under BOTH readings, and what O-4 measured at L=3 was")
say("  a statement about regions being SMALL RELATIVE TO d, not about contractibility at all.")
