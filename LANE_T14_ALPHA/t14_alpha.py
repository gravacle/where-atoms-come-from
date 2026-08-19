"""T-14: ALPHA'S ROLE, RESTATED AFTER PF-6.

Theorem D registered 'a local perturbation splits the degeneracy as eps^d'. PF-6 showed that
holds for a GENERIC single-site perturbation and NOT for a record-commuting one, where the
exponent is n* -- the fewest AVAILABLE local terms whose product acts non-scalar on the code
space. Measured here in-house, on three carriers, with n* computed rather than assumed.

The splitting is w[k-1] - w[0], the spread WITHIN the code space. Using w[k] instead measures
the gap to the first excited state, barely moves with eps, and gave a slope of -0.05 in an
earlier lane."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); Y=1j*(X@Z)
def pl(s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z,'Y':Y}[c])
    return M
def sing(n,letters):
    return [pl(''.join(c if j==q else 'I' for j in range(n))) for q in range(n) for c in letters]
def code_proj(S,n):
    P=np.eye(2**n,dtype=complex)
    for s in S: P=P@((np.eye(2**n)+s)/2)
    return P
def nstar(S,n,ops,R,Pg,k,cap=5):
    """fewest of `ops` whose PRODUCT acts non-scalar on the code space"""
    for w in range(1,cap+1):
        for combo in itertools.combinations(range(len(ops)),w):
            A=np.eye(2**n,dtype=complex)
            for i in combo: A=A@ops[i]
            M=Pg@A@Pg; c=np.trace(M)/k
            if np.linalg.norm(M-c*Pg)>1e-9: return w
    return None
def slope(H0,V,k,epss):
    dE=[]
    for e in epss:
        w=np.linalg.eigvalsh(H0+e*V); dE.append(w[k-1]-w[0])
    dE=np.array(dE)
    if dE.min()<=0: return None,dE
    return float(np.polyfit(np.log(epss),np.log(dE),1)[0]), dE

CARRIERS=[]
g={}; exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
           .read().split('say("="*104); say("0.')[0],g)
CARRIERS.append(("toric 2x2 [[8,2,2]]", 8,
    [g['op']({l:X for l in s},8) for s in g['STAR'][:3]]+[g['op']({l:Z for l in p},8) for p in g['PLAQ'][:3]],
    g['Zbar'], 2))
CARRIERS.append(("[[5,1,3]] non-CSS", 5, [pl(s) for s in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']], pl('ZZZZZ'), 3))
CARRIERS.append(("[[7,1,3]] Steane", 7,
    [pl(s) for s in ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ']], pl('ZZZZZZZ'), 3))

say("="*100); say("T-14   ALPHA'S TWO EXPONENTS, MEASURED"); say("="*100)
epss=np.array([0.02,0.05,0.10,0.20])
say(f"  {'carrier':<22}{'d':>3}{'n* (Z-only)':>13}{'GENERIC slope':>15}{'Z-ONLY slope':>14}{'':>10}")
for nm,n,S,R,d in CARRIERS:
    H0=-sum(S); Pg=code_proj(S,n); k=int(round(np.real(np.trace(Pg))))
    zs=sing(n,'Z'); ns=nstar(S,n,zs,R,Pg,k,cap=min(n,5))
    rng=np.random.default_rng(7)
    Vg=sum(rng.normal()*o for o in sing(n,'XYZ'))
    Vz=sum(zs)
    sg,_=slope(H0,Vg,k,epss); sz,_=slope(H0,Vz,k,epss)
    agree = (sg is not None and abs(sg-d)<0.25)
    say(f"  {nm:<22}{d:>3}{str(ns):>13}{(f'{sg:.4f}' if sg else 'flat'):>15}"
        f"{(f'{sz:.4f}' if sz else 'flat'):>14}{('' if agree else '  <-- generic != d'):>10}")
say("")
say("  READ: the GENERIC exponent is the code distance d. The Z-ONLY exponent is n*, the fewest")
say("  AVAILABLE terms whose product reaches the code space -- which equals d only when the")
say("  available terms happen to include a minimum-weight logical. On the toric they coincide;")
say("  on [[5,1,3]] and Steane they do not, and Theorem D's unqualified form is wrong there.")
