import sys, time, numpy as np; sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0],g)
Z=g['Z']; X=g['X']; op=g['op']; L=g['L']; Y=1j*(X@Z)
print("  carrier                                    time   min projections   record possible?", flush=True)
for lbl,Ls in (("toric 2x2, no noise", []),
               ("toric 2x2, single-site Z noise", [op({l:Z},L) for l in range(L)]),
               ("toric 2x2, single-site X,Y,Z noise", [op({l:P},L) for l in range(L) for P in (X,Y,Z)])):
    t=time.time(); m=RecordModel(g['H0'],Ls); el=time.time()-t
    k=len(m.projs)
    print(f"  {lbl:<38} {el:>6.1f}s {k:>15}   {'NO -- scalars only' if k<=1 else 'yes'}", flush=True)
