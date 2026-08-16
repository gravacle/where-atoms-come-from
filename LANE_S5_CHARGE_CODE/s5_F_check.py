import numpy as np
from itertools import product
from s5lib import *
# B0b under regime 3 -- verify by direct simulation
ge=[];ix={}
for i in range(3):
    for j in range(3): ix[('h',i,j)]=len(ge); ge.append((3*i+j,3*((i+1)%3)+j))
for i in range(3):
    for j in range(3): ix[('v',i,j)]=len(ge); ge.append((3*i+j,3*i+(j+1)%3))
gf=[{ix[('h',i,j)]:1, ix[('v',(i+1)%3,j)]:1, ix[('h',i,(j+1)%3)]:-1, ix[('v',i,j)]:-1} for i in range(3) for j in range(3)]
T=CW("T",9,ge,gf); row0={ix[('h',0,0)]:1,ix[('h',1,0)]:1,ix[('h',2,0)]:1}
a,b=T.classes(gf[0],row0); p=np.ones(9)/9
for lab,q in [("R1 unit",np.ones(9,dtype=int)),("R3 alt",np.array([1+(v%2) for v in range(9)]))]:
    E=exponents_from_charge(a,b,q); Es,ps=support_exponents(E,p)
    print(f"B0b {lab}: E={sorted(set(map(tuple,Es)))} rank={len(difference_lattice(Es))} "
          f"quad={mahler_generic(Es,ps,Nx=32768):.9f} direct={lambda_B_direct(Es,ps,1.0,np.sqrt(2),N=2_000_000):.9f}")
print("log(4/9) =", np.log(4/9))
# spread of lambda_B^gen over all per-vertex charges on K1, q in {0,1,2}^5
K1=CW("K1",5,[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)],[{0:1,1:1,2:1}])
a,b=K1.classes({0:1,1:1,2:1},{3:1,4:1,5:1}); pu=np.array([0.4,0.15,0.15,0.15,0.15])
vals=[]; nf=0
for q in product(range(3),repeat=5):
    E=exponents_from_charge(a,b,np.array(q)); Es,ps=support_exponents(E,pu)
    if len(difference_lattice(Es))==0: nf+=1; continue
    vals.append(mahler_generic(Es,ps,Nx=4096))
print(f"K1 per-vertex charge {{0,1,2}}^5: {len(vals)} forming, {nf} never-forming; "
      f"lambda_B^gen min {min(vals):.6f} max {max(vals):.6f} spread {max(vals)-min(vals):.6f}, "
      f"{len(set(round(v,7) for v in vals))} distinct values")
print("unit charge value:", mahler_generic(np.array([[1,1],[1,0],[0,1]]),np.array([0.4,0.3,0.3]),Nx=4096))
# how many carriers move under regime 3
