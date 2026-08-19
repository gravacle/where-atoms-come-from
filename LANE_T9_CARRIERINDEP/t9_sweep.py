"""T-9: does every carrier-dependent PROVED row survive a SECOND carrier?

Three carriers, chosen so the differences are structural and not cosmetic (D-12: only a change
to the STABILISER GROUP changes the carrier):
  A  toric 2x2                 [[8,2,2]]  manifold, lattice gauge theory   -- the incumbent
  B  torus + capped wrap       [[8,1,2]]  NON-MANIFOLD, one logical        -- same dim 256
  C  [[4,2,2]]                 [[4,2,2]]  not a lattice at all, dim 16     -- for composition,
                                          which B cannot test with only one record"""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def pl(s,n):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z,'Y':1j*(X@Z)}[c])
    return M
def op_on(d,n):
    M=np.array([[1]],dtype=complex)
    for l in range(n): M=np.kron(M,d.get(l,I2))
    return M

# ---- carrier A: toric 2x2 ----
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
           .read().split('say("="*104); say("0.')[0],g)
A=dict(name="A toric 2x2 [[8,2,2]]", n=8, H=g['H0'], recs=[g['Zbar'],g['Zbar2']],
       stab=[op_on({l:X for l in s},8) for s in g['STAR']]+[op_on({l:Z for l in p},8) for p in g['PLAQ']],
       loc=[op_on({l:P},8) for l in range(8) for P in (X,Z)],
       gi_local=[op_on({l:Z for l in p},8) for p in g['PLAQ']])
# ---- carrier B: torus 2x2 + a disk capping the wrap {0,2}: one extra Z-stabiliser ----
capZ=op_on({0:Z,2:Z},8)
Hb=A['H'] - capZ
# with the wrap capped, Zbar on that wrap is now a STABILISER; the surviving logical is Zbar2
B=dict(name="B torus+cap [[8,1,2]]", n=8, H=Hb, recs=[g['Zbar2']],
       stab=A['stab']+[capZ], loc=A['loc'], gi_local=A['gi_local']+[capZ])
# ---- carrier C: [[4,2,2]] ----
Sc=[pl('XXXX',4),pl('ZZZZ',4)]
C=dict(name="C [[4,2,2]] dim 16", n=4, H=-sum(Sc), recs=[pl('ZZII',4),pl('ZIZI',4)],
       stab=Sc, loc=[op_on({l:P},4) for l in range(4) for P in (X,Z)], gi_local=Sc)

say("="*104); say("T-9   CARRIER INDEPENDENCE -- every carrier-dependent PROVED row on a SECOND carrier"); say("="*104)
P=F=0
def chk(c,lbl,got,want,tol=1e-7):
    global P,F
    ok = (abs(got-want)<tol) if isinstance(want,float) else (got==want)
    P+=ok; F+=(not ok)
    say(f"    [{'PASS' if ok else 'FAIL'}] {lbl:<48} {str(got)[:14]:>15}  want {want}")

for c in (A,B,C):
    n=c['n']; m=RecordModel(c['H'],[]); env=Environment(nq=3)
    Pg,k=m.ground_space()
    say(f"\n  {c['name']}   dim {2**n}   ground-space dim {k}   logicals {len(c['recs'])}")
    # stabiliser groups must differ, or it is the same carrier (D-12)
    R=c['recs'][0]
    chk(c,"ground-space dim = 2^(#logicals)", k, 2**len(c['recs']))
    chk(c,"record commutes with H", float(np.linalg.norm(R@c['H']-c['H']@R)), 0.0)
    # C-17: generic single-site noise admits no record; Z-only does
    mz=RecordModel(c['H'],[op_on({l:Z},n) for l in range(n)])
    mg=RecordModel(c['H'],[op_on({l:P},n) for l in range(n) for P in (X,1j*(X@Z),Z)])
    chk(c,"C-17  Z-only noise: record possible", len(mz.projs)>1, True)
    chk(c,"C-17  generic noise: scalars only", len(mg.projs)==1, True)
    # G-16 + C-18 + F-20 + F-21
    chk(c,"G-16  channel(R, R) opens", m.channel(R,R)['opens_channel'], True)
    gi=sum(c['gi_local'])
    chk(c,"G-16  channel(R, gauge-inv local) shut", m.channel(R,gi)['opens_channel'], False)
    chk(c,"C-18  chi for that coupling", float(m.formation(R,[(x,i) for i,x in enumerate(c['gi_local'])],env,0.8,4.0)), 0.0)
    chk(c,"F-20  chi at t=0 is exactly 0", float(m.formation(R,R,env,0.8,0.0)), 0.0)
    x=m.formation(R,R,env,0.8,4.0); chk(c,"F-20  chi at t=4 is > 0", x>0.1, True)
    chk(c,"F-21  weight-1 coupling gives 0", float(m.formation(R,op_on({0:Z},n),env,0.8,4.0)), 0.0)
    if len(c['recs'])>=2:
        res=m.formation_independence(c['recs'],c['recs'],env,0.8,4.0)
        chk(c,"C-19  records compose independently", all(r['independent'] for r in res) and len(res)==2, True)
    else:
        say(f"    [ -- ] {'C-19  composition':<48} {'n/a':>15}  only one logical")
say("")
say("="*104); say(f"  {P} PASS, {F} FAIL across three carriers")
