# IS THE SECTOR-CHANGING STRUCTURE A PROCESS, OR A LABEL?
# W-22 found a three-way split and one operator, X_9, that moves between superselection sectors.
# That is structure. A RECORD PROCESS needs more: the sector must be able to CHANGE, the change must
# be DURABLE, and the sector must be CORRELATED with something -- otherwise it records nothing.
# W-06 killed attempt one on exactly the third clause: "carrier and record are never correlated;
# the record records which of two counterfactual transports was applied -- a CONTROL SETTING."
import numpy as np
I2=np.eye(2); Xp=np.array([[0,1],[1,0]],dtype=complex); Zp=np.diag([1,-1]).astype(complex)
def op(i,P,n):
    M=np.array([[1]],dtype=complex)
    for j in range(n): M=np.kron(M,P if j==i else I2)
    return M
OPEN=[(0,1),(1,2),(2,0),(2,3),(3,4),(4,2),(4,5),(5,0),(0,4),(0,6)]
n=len(OPEN); BULK=list(range(6))
Z=[op(i,Zp,n) for i in range(n)]; X=[op(i,Xp,n) for i in range(n)]
def gauss(v):
    M=np.eye(2**n,dtype=complex)
    for i,(a,b) in enumerate(OPEN):
        if a==v or b==v: M=M@Z[i]
    return M
Q=gauss(6)

# which links lie on a cycle? a dangling link lies on none, so it is in NO plaquette term
import itertools
def on_cycle(i):
    """link i is on a cycle iff removing it leaves its endpoints still connected"""
    a,b=OPEN[i]; adj={}
    for j,(u,v) in enumerate(OPEN):
        if j==i: continue
        adj.setdefault(u,[]).append(v); adj.setdefault(v,[]).append(u)
    seen={a}; stack=[a]
    while stack:
        x=stack.pop()
        for y in adj.get(x,[]):
            if y not in seen: seen.add(y); stack.append(y)
    return b in seen
cyc=[i for i in range(n) if on_cycle(i)]
print(f"  links on a cycle (hence appearing in plaquette terms): {cyc}")
print(f"  links on NO cycle (appearing in NO plaquette term)   : {[i for i in range(n) if i not in cyc]}")

# standard pure-gauge Hamiltonian: magnetic (X around independent cycles) + electric (Z on links)
cycles=[[0,1,2],[3,4,5],[6,7,8]]
def wil(c):
    M=np.eye(2**n,dtype=complex)
    for i in c: M=M@X[i]
    return M
for g2 in (0.3,1.0,3.0):
    H=-(1.0/g2)*sum(wil(c) for c in cycles)-g2*sum(Z)
    print(f"\n  g^2={g2}:  || [H, Q] ||  =  {np.linalg.norm(H@Q-Q@H):.3e}")
print("\n  -> if [H,Q] = 0 the boundary charge is EXACTLY CONSERVED: the sector is fixed at t=0")
print("     and never changes. That is a LABEL, not a process. The structure would be real and inert.")

print("\n== WHAT WOULD UNFREEZE IT ==")
Hm=-(1.0/1.0)*sum(wil(c) for c in cycles)-1.0*sum(Z)-0.7*X[9]   # a term containing the dangling link
print(f"  add a term containing X_9 (the sector-changing operator):")
print(f"    || [H + 0.7*X_9, Q] || = {np.linalg.norm(Hm@Q-Q@Hm):.3e}")
print("  X_9 lies on NO cycle, so NO plaquette term can contain it. In PURE gauge theory there is")
print("  nothing in the Hamiltonian that can move boundary charge. The only way to get such a term")
print("  is MATTER at the boundary vertex -- a field that can create and destroy charge.")
