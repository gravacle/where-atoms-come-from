# FW1 was confounded: the gauge arm got a PROPER fragment (3 of 8) and the spin arm got the WHOLE
# complement (8 of 8), so its 2.000000 was the purification identity, not a constraint effect.
# Redone with MATCHED fragment sizes. ONE VARIABLE: local gauge constraint vs global conservation.
# THE QUESTION: does a PROPER, LOCAL fragment determine the system? That is what a SURFACE means.
import numpy as np
I2=np.eye(2); Zp=np.diag([1,-1]).astype(complex)
def op(i,P,n):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
def red(rho,keep,n):
    d=[2]*n; r=rho.reshape(d+d); tr=[i for i in range(n) if i not in keep]
    for k,i in enumerate(sorted(tr,reverse=True)): r=np.trace(r,axis1=i,axis2=i+n-k)
    return r.reshape(2**len(keep),2**len(keep))
def S(r):
    w=np.linalg.eigvalsh(r); return float(-sum(x*np.log2(x) for x in w if x>1e-12))
rng=np.random.default_rng(20260828)
def ratio(P,n,region,frag,trials=3):
    w_,vec=np.linalg.eigh(P); B=vec[:,w_>0.5]; d=B.shape[1]; out=[]
    for _ in range(trials):
        c=rng.normal(size=d)+1j*rng.normal(size=d); psi=B@c; psi/=np.linalg.norm(psi)
        rho=np.outer(psi,psi.conj()); hs=S(red(rho,region,n))
        out.append((hs+S(red(rho,frag,n))-S(red(rho,sorted(region+frag),n)))/hs)
    return d,out

EG=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4)]
nA=len(EG); ZA=[op(i,Zp,nA) for i in range(nA)]
def gauss(v):
    M=np.eye(2**nA,dtype=complex)
    for i,(a,b) in enumerate(EG):
        if a==v or b==v: M=M@ZA[i]
    return M
PA=np.eye(2**nA,dtype=complex)
for v in range(6): PA=PA@((np.eye(2**nA)+gauss(v))/2)

nB=9; Mtot=sum(op(i,Zp,nB) for i in range(nB))
w,v=np.linalg.eigh(Mtot)
PB=np.zeros((2**nB,2**nB),dtype=complex)
for k in range(len(w)):
    if abs(w[k]-1)<1e-9: PB+=np.outer(v[:,k],v[:,k].conj())

print("  MATCHED FRAGMENT SIZES. system = site/link 0 in both arms.")
print(f"  {'fragment':<34}{'|F|':>4}{'GAUGE  I/H(S)':>18}{'GLOBAL-CONSERVED  I/H(S)':>28}")
FRAGS=[("the Gauss cut at vertex 0",[2,7,8]),
       ("three other sites",[1,2,3]),
       ("five other sites",[1,2,3,4,5]),
       ("the whole complement",[1,2,3,4,5,6,7,8])]
for tag,F in FRAGS:
    dA,rA=ratio(PA,nA,[0],F); dB,rB=ratio(PB,nB,[0],F)
    print(f"  {tag:<34}{len(F):>4}{np.mean(rA):>18.6f}{np.mean(rB):>28.6f}")
print()
print("  GAUGE: a PROPER, LOCAL fragment -- the Gauss cut, 3 of 8 -- already gives 1.000000.")
print("         The determining set is a SURFACE.")
print("  GLOBAL: no proper fragment determines the system; only the whole complement does, and then")
print("         only via purification. There is NO surface, because the constraint is not local.")
print()
print("  ==> WHAT GAUGE BUYS IS THAT THE DETERMINING SET IS LOCAL. That is the whole of it, and it")
print("      is exactly what a 'record surface' would need. But PURE gauge freezes the charge")
print("      ([H,Q] = 0 exactly), so gauge alone gives a surface with nothing happening on it.")
