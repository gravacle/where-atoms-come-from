import sys, numpy as np; sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0],g)
H0=g['H0']; Zbar=g['Zbar']; Zbar2=g['Zbar2']; Xbar=g['Xbar']; Z=g['Z']; op=g['op']; L=g['L']
m=RecordModel(H0,[]); env=Environment(); Pg,k=m.ground_space()
def zop(links): return op({l:Z for l in links},L)
def on_code(A):
    """how A acts on the code space, as a multiple of Zbar or of the identity"""
    M=Pg@A@Pg; Zc=Pg@Zbar@Pg
    a=np.trace(M.conj().T@Zc)/np.trace(Zc.conj().T@Zc)      # component along Zbar
    b=np.trace(M)/k                                          # component along identity
    resid=np.linalg.norm(M - a*Zc - b*Pg)
    return complex(a), complex(b), float(resid)
print("  WHY THE WEIGHT-4 CYCLES FAIL: what does each coupling do ON THE CODE SPACE?\n")
print(f"  {'links':<18}{'||[A,Xbar]||':>14}{'~ Zbar':>10}{'~ I':>8}{'residual':>11}{'chi':>13}")
for links in ([0,2],[4,6],[3,7],[1,5],[3,4,6,7],[2,3,4,5],[0,3,5,6]):
    A=zop(links); a,b,r=on_code(A)
    c=np.linalg.norm(A@Xbar-Xbar@A)
    x=m.formation(Zbar,A,env,lam=0.8,t=4.0)
    print(f"  {str(links):<18}{c:>14.1f}{abs(a):>10.3f}{abs(b):>8.3f}{r:>11.3f}{x:>13.8f}")
print()
print("  Zbar2 vs Xbar commutator :", f"{np.linalg.norm(Zbar2@Xbar-Xbar@Zbar2):.1f}")
print("  Zbar*Zbar2 vs Xbar       :", f"{np.linalg.norm((Zbar@Zbar2)@Xbar-Xbar@(Zbar@Zbar2)):.1f}")
print()
print("  READ: a coupling can anticommute with the writer and still be a DIFFERENT logical")
print("        (Zbar*Zbar2 anticommutes with Xbar but is not Zbar). The bath then learns THAT")
print("        logical, which in a maximally mixed code state says nothing about Zbar.")
