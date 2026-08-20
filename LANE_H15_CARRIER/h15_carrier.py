"""H-15: CARRIER CURVATURE -- transport along the carrier, the object H-14 fixed.

Our carrier is a lattice gauge theory. The gauge field lives on LINKS, so a link variable IS a
connection variable, and transport around a face IS the plaquette operator B_p. Carrier curvature
is therefore the plaquette holonomy, and it has been defined in this program from the start.

THE CODE SPACE IS THE SECTOR WHERE EVERY PLAQUETTE IS SATISFIED -- zero holonomy around every face.
So the questions are:
  1. is the carrier flat where records live?
  2. do records survive when the carrier is CURVED (a violated plaquette = a flux)?
  3. does record content SOURCE curvature, or require its absence?"""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
           .read().split('say("="*104); say("0.')[0],g)
op=g['op']; L=g['L']; Z=g['Z']; X=g['X']; PLAQ=g['PLAQ']; STAR=g['STAR']
H0=g['H0']; Zbar=g['Zbar']; Zbar2=g['Zbar2']; N=2**L
BP=[op({l:Z for l in p},L) for p in PLAQ]          # transport around each face
AV=[op({l:X for l in s},L) for s in STAR]
say("="*100); say("H-15   CARRIER CURVATURE: TRANSPORT AROUND THE FACES OF THE CARRIER"); say("="*100)
say(f"  carrier: toric 2x2, {L} links, {len(PLAQ)} faces, dim {N}")
say(f"  a link variable IS a connection variable; transport around face p IS B_p")
w,V=np.linalg.eigh(H0); gs=int(np.sum(np.abs(w-w[0])<1e-9)); Pg=V[:,:gs]@V[:,:gs].conj().T
say("")
say("1.  IS THE CARRIER FLAT WHERE RECORDS LIVE?")
say(f"  {'face':<8}{'<B_p> on the code space':>26}{'holonomy':>12}")
for i,B in enumerate(BP):
    val=np.real(np.trace(Pg@B))/gs
    say(f"  p={i:<6}{val:>26.8f}{('FLAT' if abs(val-1)<1e-9 else 'CURVED'):>12}")
say("")
say("2.  DO RECORDS SURVIVE WHERE THE CARRIER IS CURVED?")
say("    A flux sector is reached by flipping one plaquette's holonomy to -1. The sector is an")
say("    eigenspace of H0 with the same dimension, so ask the MODEL whether a record exists there.")
# project onto the sector with B_0 = -1 and all others +1
Pflux=np.eye(N,dtype=complex)
for i,B in enumerate(BP):
    s = -1 if i==0 else 1
    Pflux = Pflux @ ((np.eye(N)+s*B)/2)
for A in AV: Pflux = Pflux @ ((np.eye(N)+A)/2)
dflux=int(round(np.real(np.trace(Pflux))))
say(f"    flux sector dimension: {dflux}   (code space: {gs})")
if dflux>0:
    say(f"    {'observable':<16}{'<O> on the flux sector':>26}{'is it still a record?':>24}")
    for nm,R in (("Zbar",Zbar),("Zbar2",Zbar2)):
        commH=np.linalg.norm(R@H0-H0@R)
        M=Pflux@R@Pflux; c=np.trace(M)/max(dflux,1)
        nontriv=np.linalg.norm(M-c*Pflux)
        say(f"    {nm:<16}{np.real(np.trace(M))/max(dflux,1):>26.8f}"
            f"{('YES, non-trivial' if nontriv>1e-9 else 'trivial here'):>24}")
say("")
say("3.  DOES RECORD CONTENT SOURCE CURVATURE?")
say("    Write a record and see whether any face's holonomy changes.")
say(f"  {'state':<28}" + "".join(f"{'<B_'+str(i)+'>':>10}" for i in range(len(BP))))
rho0=Pg/gs
for nm,R in (("code space, no record set",None),("Zbar = +1",Zbar),("Zbar2 = +1",Zbar2),
             ("both = +1", None)):
    if nm=="both = +1":
        P=Pg@((np.eye(N)+Zbar)/2)@((np.eye(N)+Zbar2)/2)
    elif R is None: P=Pg
    else: P=Pg@((np.eye(N)+R)/2)
    tr=np.real(np.trace(P))
    if tr<1e-9: say(f"  {nm:<28}   empty"); continue
    say(f"  {nm:<28}" + "".join(f"{np.real(np.trace(P@B))/tr:>10.5f}" for B in BP))
say("")
say("  READ: if every holonomy stays +1 whatever record is written, then record content does NOT")
say("  source carrier curvature, and records live strictly in the FLAT sector.")
