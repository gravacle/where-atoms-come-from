"""T-20: RECORD-LEVEL HOLONOMY. The criterion is H-9, registered before this ran.

  H_ABCD = (G_AB G_BD)(G_AC G_CD)^-1      H = I -> flat;  H != I structured -> curvature

WHAT 'TRANSPORT BETWEEN RECORDS' IS HERE. A record is written by its admissible writer (C-13), so
the natural transport from one record's frame to another's is the writer that flips it. Comparing
two chains between the same endpoints is then comparing two ORDERS of writing:

    path 1:  A then B        path 2:  B then A
    H = (U_A U_B)(U_B U_A)^-1

which is I exactly when the two writers commute. So the diagnostic asks: DOES THE ORDER IN WHICH
RECORDS ARE WRITTEN MATTER? That is holonomy in the ordinary sense -- transport around a closed
loop, here the loop A->B->A^-1->B^-1.

A second reading is measured too, because the first could be accused of being trivially group
theory: transport by the FORMATION channel rather than by the writer -- does forming A then B leave
the same state as forming B then A?"""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment, build_writer
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def pl(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
say("="*100); say("T-20   RECORD-LEVEL HOLONOMY:  DOES TRANSPORT AROUND A RECORD LOOP CLOSE?"); say("="*100)

CAR=[]
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
           .read().split('say("="*104); say("0.')[0],g)
CAR.append(("toric 2x2 [[8,2,2]]", g['H0'], [g['Zbar'],g['Zbar2']], [g['Xbar'], None]))
CAR.append(("[[4,2,2]] dim 16", -(pl('XXXX')+pl('ZZZZ')), [pl('ZZII'),pl('ZIZI')], [pl('XXII'),pl('XIXI')]))

say("\n1.  WRITER TRANSPORT -- does the ORDER of writing two records matter?")
say(f"  {'carrier':<24}{'||[U_A,U_B]||':>15}{'||H - I||':>12}{'verdict':>12}")
for nm,H,recs,writers in CAR:
    m=RecordModel(H,[])
    Ws=[]
    for i,R in enumerate(recs):
        U = writers[i] if writers[i] is not None else build_writer(R,m.es)
        Ws.append(U)
    UA,UB=Ws[0],Ws[1]
    if UA is None or UB is None: say(f"  {nm:<24}   writer construction failed"); continue
    comm=float(np.linalg.norm(UA@UB-UB@UA))
    Hol = (UA@UB) @ np.linalg.inv(UB@UA)
    dev = float(np.linalg.norm(Hol-np.eye(H.shape[0])))
    say(f"  {nm:<24}{comm:>15.3e}{dev:>12.3e}{('FLAT' if dev<1e-9 else 'CURVED'):>12}")

say("\n2.  FORMATION TRANSPORT -- does forming A then B leave the same state as B then A?")
say(f"  {'carrier':<24}{'||rho_AB - rho_BA||':>21}{'chi diff':>11}{'verdict':>12}")
env=Environment(nq=3)
for nm,H,recs,_ in CAR:
    m=RecordModel(H,[]); nS=H.shape[0]
    A,B=recs[0],recs[1]
    rAB=m.evolve(A,env,lam=0.8,t=2.0)
    # second leg: evolve the SYSTEM part again under the other coupling
    rS=rAB.reshape(nS,env.dim,nS,env.dim).trace(axis1=1,axis2=3)
    rAB2=m.evolve(B,env,lam=0.8,t=2.0,state0=rS/np.trace(rS))
    rBA=m.evolve(B,env,lam=0.8,t=2.0)
    rS2=rBA.reshape(nS,env.dim,nS,env.dim).trace(axis1=1,axis2=3)
    rBA2=m.evolve(A,env,lam=0.8,t=2.0,state0=rS2/np.trace(rS2))
    d=float(np.linalg.norm(rAB2-rBA2))
    cA=env.holevo(rAB2,A,nS); cB=env.holevo(rBA2,A,nS)
    say(f"  {nm:<24}{d:>21.3e}{abs(cA-cB):>11.3e}{('FLAT' if d<1e-9 else 'CURVED'):>12}")

say("\n3.  WHY, AND WHAT WOULD BE NEEDED FOR CURVATURE")
say("     Our records are commuting logical operators of an ABELIAN stabiliser code, and their")
say("     writers are Pauli operators that either commute or anticommute -- so a loop closes up to")
say("     a SIGN, and a sign is not curvature: it is the same element of an abelian group.")
say("     Curvature needs transport whose composition depends on ORDER in a way no relabelling")
say("     removes -- a NON-ABELIAN connection. That is exactly the class O-8 flagged as untested:")
say("     Levin-Wen and Fibonacci, where the writers are Wilson loops of a non-abelian theory.")
