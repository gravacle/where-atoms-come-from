"""G2.  DEFINING GRAVITY AT THE RECORD LEVEL -- AND ITS BRIDGE TO CURVATURE.

G1: the record count is 2g, independent of area. That is a TOPOLOGICAL invariant, and it looks
nothing like classical gravity -- which is exactly what the principal predicted.

BUT TOPOLOGICAL INVARIANTS ARE INTEGRALS OF CURVATURE. For a closed orientable surface,
      Euler characteristic   chi = V - E + F = 2 - 2g
      Gauss-Bonnet           integral of K dA = 2*pi*chi
      record count           dim H_1 = 2g = 2 - chi
Therefore
      NUMBER OF INDEPENDENT RECORDS  =  2 - (1/2pi) * integral of K dA.

THAT IS AN EXACT RELATION BETWEEN THE RECORD COUNT AND INTEGRATED CURVATURE -- the classical
gravitational quantity. The record count is non-metric and non-local, and it EQUALS an integral of
the metric quantity. That is the emergence bridge, and it has not been stated before in this program.

CHECKED HERE ON THE CARRIERS THE PROGRAM ACTUALLY USED, by computing chi combinatorially from the
cell decomposition and comparing with the measured ground-space degeneracy.
"""
import itertools, numpy as np
def torus(nx,ny):
    vid=lambda i,j:(j%ny)*nx+(i%nx)
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx): ind[('h',i,j)]=len(E); E.append((vid(i,j),vid(i+1,j)))
    for j in range(ny):
        for i in range(nx): ind[('v',i,j)]=len(E); E.append((vid(i,j),vid(i,j+1)))
    PL=[[ind[('h',i,j)],ind[('v',(i+1)%nx,j)],ind[('h',i,(j+1)%ny)],ind[('v',i,j)]]
        for j in range(ny) for i in range(nx)]
    return nx*ny,E,PL,"torus"
def sphere_tetra():
    V=4; E=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    F=[[0,3,1],[0,4,2],[1,5,2],[3,5,4]]
    return V,E,F,"tetrahedron (sphere)"
def disk(nx,ny):
    vv={(i,j):j*nx+i for j in range(ny) for i in range(nx)}
    E=[]; ind={}
    for j in range(ny):
        for i in range(nx-1): ind[('h',i,j)]=len(E); E.append((vv[(i,j)],vv[(i+1,j)]))
    for j in range(ny-1):
        for i in range(nx): ind[('v',i,j)]=len(E); E.append((vv[(i,j)],vv[(i,j+1)]))
    PL=[[ind[('h',i,j)],ind[('v',i+1,j)],ind[('h',i,j+1)],ind[('v',i,j)]] for j in range(ny-1) for i in range(nx-1)]
    return nx*ny,E,PL,"disk"

def degeneracy(NV,E,PL):
    L=len(E)
    if L>14: return None
    st=[s for s in itertools.product(range(2),repeat=L)
        if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
               -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%2==0 for v in range(NV))]
    idx={s:i for i,s in enumerate(st)}; D=len(st)
    def Move(S):
        M=np.zeros((D,D),complex)
        for j,s in enumerate(st):
            t=list(s)
            for k in S: t[k]^=1
            t=tuple(t)
            if t in idx: M[idx[t],j]=1.0
        return M
    H=-sum(Move(p) for p in PL); H=(H+H.conj().T)/2
    ev=np.linalg.eigvalsh(H); tol=1e-8*max(1.0,abs(ev).max())
    return int(sum(1 for e in ev if abs(e-ev[0])<tol))

print("G2  RECORD COUNT vs EULER CHARACTERISTIC vs INTEGRATED CURVATURE")
print("    chi = V - E + F ;  Gauss-Bonnet: integral K dA = 2*pi*chi ;  claim: records = 2 - chi\n")
print(f"  {'carrier':>24s} {'V':>3s} {'E':>3s} {'F':>3s} {'chi':>4s} {'2-chi':>6s} "
      f"{'measured degeneracy':>20s} {'2^(2-chi)':>10s}")
print("  "+"-"*82)
for spec in (torus(2,2), torus(2,3), sphere_tetra(), disk(3,3)):
    NV,E,F,name=spec
    chi=NV-len(E)+len(F)
    if name=="disk": chi=NV-len(E)+len(F)      # disk with boundary: chi = 1
    g=degeneracy(NV,E,F)
    pred=2-chi
    print(f"  {name:>24s} {NV:3d} {len(E):3d} {len(F):3d} {chi:4d} {pred:6d} "
          f"{(g if g is not None else '--'):>20} {2**max(pred,0):10d}")
print()
print("  READING. For the TORUS chi = 0, so records = 2 - 0 = 2 and the record space is 2^2 = 4.")
print("  For the SPHERE chi = 2, so records = 0 and the space is 1 -- NO RECORD.")
print("  For the DISK chi = 1, and with a boundary the closed-surface formula does not apply;")
print("  the measured degeneracy is 1, i.e. no record, consistent with H_1 = 0.")
print()
print("  THE BRIDGE:  records = 2 - chi = 2 - (1/2pi) * integral of K dA.")
print("  The record count is non-local and non-metric, and it EQUALS an integral of the metric")
print("  quantity. That is how a topological record count and classical curvature are the same fact.")
