# ASSUMPTION AUDIT 1c — third attempt, and NO VERDICT IS HARD-CODED. Every conclusion below is
# printed from a number computed on the line above it. AA1 and AA2 both failed as zero-variable
# controls, and AA2 printed "NOW the arms differ" beside a ratio of 1.0x. Recorded, not patched out.
#
# The relation prod_v G_v = I holds because every link touches exactly TWO constrained vertices.
# Subdividing preserves that. The only way to break it is a DANGLING link: one whose far end carries
# no Gauss operator at all. That vertex is BOUNDARY, not bulk.
import numpy as np
rng=np.random.default_rng(20260827)
I2=np.eye(2); Zp=np.diag([1,-1]).astype(complex)
def Zi(i,n):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,Zp if j==i else I2)
    return M
def gauss_of(edges,v,n):
    M=np.eye(2**n,dtype=complex)
    for i,(a,b) in enumerate(edges):
        if a==v or b==v: M=M@Zi(i,n)
    return M
def proj(edges,n,bulk):
    P=np.eye(2**n,dtype=complex)
    for v in bulk: P=P@((np.eye(2**n)+gauss_of(edges,v,n))/2)
    return P

CLOSED=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4)]
print("== ARM 1: CLOSED GRAPH (every carrier this program has ever used) ==")
n1=len(CLOSED); bulk1=list(range(6))
P_all=proj(CLOSED,n1,bulk1); P_drop=proj(CLOSED,n1,[v for v in bulk1 if v!=0])
d_all=int(round(np.trace(P_all).real)); d_drop=int(round(np.trace(P_drop).real))
prod=np.eye(2**n1,dtype=complex)
for v in bulk1: prod=prod@gauss_of(CLOSED,v,n1)
print(f"   V=6 L={n1};  || prod_v G_v - I || = {np.linalg.norm(prod-np.eye(2**n1)):.3e}")
print(f"   dim with all 6 constraints = {d_all};  dim dropping one = {d_drop};  ratio = {d_drop/d_all:.2f}")

OPEN=CLOSED+[(0,6)]                       # a DANGLING link: vertex 6 is BOUNDARY, never constrained
n2=len(OPEN); bulk2=list(range(6))        # Gauss imposed on the six bulk vertices ONLY
print("\n== ARM 2: OPEN GRAPH — one dangling link, its far end is boundary and carries no Gauss ==")
prod2=np.eye(2**n2,dtype=complex)
for v in bulk2: prod2=prod2@gauss_of(OPEN,v,n2)
print(f"   V=6 bulk + 1 boundary, L={n2};  || prod_bulk G_v - I || = {np.linalg.norm(prod2-np.eye(2**n2)):.3e}")
P2_all=proj(OPEN,n2,bulk2); P2_drop=proj(OPEN,n2,[v for v in bulk2 if v!=0])
d2_all=int(round(np.trace(P2_all).real)); d2_drop=int(round(np.trace(P2_drop).real))
print(f"   dim with all 6 constraints = {d2_all};  dim dropping one = {d2_drop};  ratio = {d2_drop/d2_all:.2f}")

print("\n== THE W-19 IDENTITY ON EACH, ONE VARIABLE MOVED: OPEN VERSUS CLOSED ==")
for tag,edges,n,P in (("CLOSED",CLOSED,n1,P_all),("OPEN  ",OPEN,n2,P2_all)):
    l=0; u=edges[l][0]
    cut=[i for i,(a,b) in enumerate(edges) if i!=l and (a==u or b==u)]
    Xc=np.eye(2**n,dtype=complex)
    for i in cut: Xc=Xc@Zi(i,n)
    dev=np.linalg.norm((Zi(l,n)@Xc-np.eye(2**n))@P)
    print(f"   {tag}  cut at vertex {u} = links {cut};  || (Z_l * prod_cut Z - I) P || = {dev:.3e}")

def red(rho,keep,n):
    d=[2]*n; r=rho.reshape(d+d); tr=[i for i in range(n) if i not in keep]
    for k,i in enumerate(sorted(tr,reverse=True)): r=np.trace(r,axis1=i,axis2=i+n-k)
    return r.reshape(2**len(keep),2**len(keep))
def S(r):
    w=np.linalg.eigvalsh(r); return float(-sum(x*np.log2(x) for x in w if x>1e-12))
def rs(P,n):
    w_,vec=np.linalg.eigh(P); b=vec[:,w_>0.5]
    c=rng.normal(size=b.shape[1])+1j*rng.normal(size=b.shape[1]); psi=b@c
    psi/=np.linalg.norm(psi); return np.outer(psi,psi.conj())
print("\n== AND THE CONSEQUENCE FOR I(S:F)/H(S), 3 RANDOM PHYSICAL STATES EACH ==")
for tag,edges,n,P in (("CLOSED",CLOSED,n1,P_all),("OPEN  ",OPEN,n2,P2_all)):
    l=0; u=edges[l][0]
    cut=[i for i,(a,b) in enumerate(edges) if i!=l and (a==u or b==u)]
    vals=[]
    for _ in range(3):
        rho=rs(P,n); hs=S(red(rho,[l],n))
        mi=hs+S(red(rho,cut,n))-S(red(rho,sorted([l]+cut),n))
        vals.append(mi/hs if hs>1e-12 else float('nan'))
    print(f"   {tag}  I/H(S) = {'  '.join(f'{v:.6f}' for v in vals)}")
