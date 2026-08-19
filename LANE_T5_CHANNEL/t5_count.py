import sys, numpy as np; sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0],g)
H0=g['H0']; Zbar=g['Zbar']; Z=g['Z']; op=g['op']; L=g['L']
m=RecordModel(H0,[]); env=Environment()
Z1=[int(g['toint'](v)) for v in g['nullspace'](g['d1'])]
cyc=sorted({np.bitwise_xor.reduce([Z1[i] for i in range(len(Z1)) if (mm>>i)&1]).item()
            for mm in range(1,1<<len(Z1))} - {0})
def zop(c): return op({l:Z for l in range(L) if (c>>(L-1-l))&1},L)
ok=bad=0; opens=0
for v in cyc:
    A=zop(v); ch=m.channel(Zbar,A); x=m.formation(Zbar,A,env,lam=0.8,t=4.0)
    opens += ch['opens_channel']
    if ch['opens_channel']==(x>1e-10): ok+=1
    else: bad+=1
print(f"  cycles tested            : {len(cyc)}")
print(f"  channel() predicts OPEN  : {opens}")
print(f"  agreement with measured chi : {ok} of {len(cyc)}   mismatches {bad}")
print(f"  -> {'T-5 DONE_WHEN MET' if bad==0 else 'NOT MET'}")
