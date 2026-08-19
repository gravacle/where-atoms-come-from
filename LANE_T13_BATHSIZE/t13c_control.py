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
print("  CONTROL: is the uneven redundancy STRUCTURAL, or just my random bath energies?\n")
print(f"  {'bath energies':<20}{'nq':>4}{'whole':>9}{'mean frag':>11}{'spread':>9}{'rel spread':>12}")
for lbl, mk in (("RANDOM (as run)", lambda q: tuple(0.6+1.2*np.random.default_rng(1).random(q))),
                ("IDENTICAL", lambda q: tuple([1.0]*q))):
    for q in (3,5,7):
        env=Environment(nq=q, energies=mk(q), beta=2.0)
        whole,fr=m.redundancy(Zbar,Zbar,env,lam=0.8,t=8.0)
        rel=(fr.max()-fr.min())/max(fr.mean(),1e-12)
        print(f"  {lbl:<20}{q:>4}{whole:>9.4f}{fr.mean():>11.4f}{fr.max()-fr.min():>9.4f}{rel:>12.3f}")
print()
print("  If IDENTICAL energies give spread ~0, the unevenness is the energies and redundancy")
print("  DOES even out for equivalent fragments -- which is the case quantum Darwinism assumes.")
print("  If it stays large, the unevenness is structural and different observers genuinely")
print("  get different amounts, which would matter for calling the record OBJECTIVE.")
