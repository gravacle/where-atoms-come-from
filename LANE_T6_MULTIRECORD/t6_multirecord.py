import sys, numpy as np; sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0],g)
H0=g['H0']; Zbar=g['Zbar']; Zbar2=g['Zbar2']; Z=g['Z']; op=g['op']; L=g['L']; PLAQ=g['PLAQ']
m=RecordModel(H0,[]); env=Environment()
fam=[Zbar,Zbar2]; names=['Zbar','Zbar2']
BP=[op({l:Z for l in p},L) for p in PLAQ]
coup=[Zbar, Zbar2, Zbar@Zbar2, sum(BP)]
cnames=['Zbar','Zbar2','Zbar*Zbar2','sum plaquettes']
print("T-6  MULTI-RECORD DEPENDENCIES\n")
print("1. CHANNEL MAP -- which coupling serves which record")
M=m.channel_map(fam,coup)
print(f"   {'coupling':<18}" + "".join(f"{n:>10}" for n in names))
for i,cn in enumerate(cnames):
    print(f"   {cn:<18}" + "".join(f"{str(bool(M[i,j])):>10}" for j in range(len(fam))))
print()
print("2. CAN ONE BE FORMED WITHOUT DISTURBING THE OTHER? -- measured")
res=m.formation_independence(fam,coup,env,lam=0.8,t=4.0)
print(f"   {'coupling':<18}{'targets':>10}{'chi(Zbar)':>12}{'chi(Zbar2)':>12}{'|d<Zbar>|':>12}{'|d<Zbar2>|':>12}{'indep':>8}")
for r in res:
    print(f"   {cnames[r['coupling']]:<18}{str([names[t] for t in r['targets']]):>10}"
          f"{r['learned'][0]:>12.8f}{r['learned'][1]:>12.8f}{r['moved'][0]:>12.2e}{r['moved'][1]:>12.2e}"
          f"{str(r['independent']):>8}")
