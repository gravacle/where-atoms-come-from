"""H-15b: DO RECORDS SURVIVE ON A CURVED CARRIER?

v1 tried to flip ONE plaquette and got a sector of dimension 0. That is not a null, it is a
constraint: on a torus every link lies in exactly two faces, so prod_p B_p = I and the number of
faces with holonomy -1 must be EVEN. A single flux is forbidden by the carrier itself.

So flip TWO. That is a genuinely curved carrier -- two faces with non-trivial holonomy -- and the
question is whether a record exists there at all, asked of the model rather than assumed."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
           .read().split('say("="*104); say("0.')[0],g)
op=g['op']; L=g['L']; Z=g['Z']; X=g['X']; PLAQ=g['PLAQ']; STAR=g['STAR']
H0=g['H0']; Zbar=g['Zbar']; Zbar2=g['Zbar2']; N=2**L
BP=[op({l:Z for l in p},L) for p in PLAQ]; AV=[op({l:X for l in s},L) for s in STAR]
say("="*100); say("H-15b   DO RECORDS SURVIVE ON A CURVED CARRIER?"); say("="*100)
prod=np.eye(N,dtype=complex)
for B in BP: prod=prod@B
say(f"  carrier constraint: ||prod_p B_p - I|| = {np.linalg.norm(prod-np.eye(N)):.3e}")
say(f"  -> the number of faces with holonomy -1 must be EVEN. A single flux is forbidden.")
say("")
def sector(flip):
    P=np.eye(N,dtype=complex)
    for i,B in enumerate(BP): P=P@((np.eye(N)+(-1 if i in flip else 1)*B)/2)
    for A in AV: P=P@((np.eye(N)+A)/2)
    return P
say(f"  {'sector':<26}{'dim':>5}{'curved faces':>14}{'Zbar non-trivial?':>20}{'Zbar2':>14}")
for flip,lbl in (((),"flat (the code space)"), ((0,1),"two fluxes"), ((0,3),"two fluxes, apart"),
                 ((0,1,2,3),"all four faces")):
    P=sector(set(flip)); d=int(round(np.real(np.trace(P))))
    if d==0: say(f"  {lbl:<26}{d:>5}{len(flip):>14}   empty sector"); continue
    out=[]
    for R in (Zbar,Zbar2):
        M=P@R@P; c=np.trace(M)/d
        out.append("YES" if np.linalg.norm(M-c*P)>1e-9 else "no")
    say(f"  {lbl:<26}{d:>5}{len(flip):>14}{out[0]:>20}{out[1]:>14}")
say("")
say("  AND DOES THE RECORD COUNT CHANGE WITH CURVATURE?")
say(f"  {'sector':<26}{'dim':>5}{'records (2^k)':>16}")
for flip,lbl in (((),"flat"), ((0,1),"two fluxes"), ((0,1,2,3),"all four")):
    P=sector(set(flip)); d=int(round(np.real(np.trace(P))))
    say(f"  {lbl:<26}{d:>5}{d:>16}")
say("")
say("  READ: if the curved sectors carry records of the same count and the same operators, then")
say("  carrier curvature and record content are INDEPENDENT on this carrier -- neither sources nor")
say("  constrains the other -- and the triad's proportionality has coefficient zero here.")
