import sys, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def pauli(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
H0=-sum(pauli(g) for g in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']); Zbar=pauli('ZZZZZ')
m=RecordModel(H0,[])
print(f"  {'bath':<8}{'whole':>10}{'mean frag':>12}{'spread':>10}{'rel spread':>12}   per-fragment", flush=True)
for q in (3,5,7):
    rng=np.random.default_rng(1); en=tuple(0.6+1.2*rng.random(q))
    env=Environment(nq=q,energies=en,beta=2.0)
    whole,fr=m.redundancy(Zbar,Zbar,env,lam=0.8,t=8.0)
    print(f"  nq={q:<5}{whole:>10.4f}{fr.mean():>12.4f}{fr.max()-fr.min():>10.4f}"
          f"{(fr.max()-fr.min())/max(fr.mean(),1e-12):>12.3f}   [{', '.join(f'{x:.3f}' for x in fr)}]", flush=True)
