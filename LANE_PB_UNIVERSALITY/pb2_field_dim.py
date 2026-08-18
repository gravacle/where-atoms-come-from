"""PHASE B / B2-B3.  DOES T3 SURVIVE CHANGING THE FIELD AND THE DIMENSION?

T3 (capacity = m - 1) was proved over GF(2) on 2-complexes. The proof uses only:
  (i) the boundary map is LINEAR over a field;
  (ii) protectability via a cell L is the kernel of ONE linear functional.
Neither mentions the field's characteristic or the dimension. So it should hold over GF(p) for any
prime p, and in any dimension. Test both, since "should" is not "does".

B2  FIELD: Z_N gauge groups for N = 2, 3, 5. The boundary of a region is the net signed winding on
    each link, mod N -- linear over GF(N).
B3  DIMENSION: a 3-complex. Regions are sets of CUBES; their boundaries are FACES; protectability is
    "some face lies on no boundary". Nothing in the proof cared that boundaries were 1-dimensional.

WHERE THE PROOF IS EXPECTED TO STOP, stated in advance: it needs an ABELIAN group, because a
non-abelian holonomy depends on ordering and base point and the boundary map is then not linear.
That limit is assessed, not tested, and is recorded as untested.
"""
import itertools, numpy as np

def rank_gfp(rows,ncols,p):
    M=[r[:] for r in rows]; r=0
    for c in range(ncols):
        piv=None
        for i in range(r,len(M)):
            if M[i][c]%p: piv=i; break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        inv=pow(M[r][c],p-2,p)
        M[r]=[(x*inv)%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%p:
                f=M[i][c]
                M[i]=[(M[i][k]-f*M[r][k])%p for k in range(ncols)]
        r+=1
        if r==len(M): break
    return r

def capacity_gfp(cells,ncell_lower,p,maxcoef=None):
    """cells: list of dicts {lower_cell_index: signed incidence}. capacity over GF(p)."""
    m=len(cells)
    P=list(itertools.product(range(p),repeat=m))[1:]        # nonzero coefficient vectors
    def bd(coef):
        acc={}
        for a,c in zip(coef,cells):
            if a%p==0: continue
            for lk,sg in c.items(): acc[lk]=(acc.get(lk,0)+a*sg)%p
        return set(lk for lk,v in acc.items() if v%p)
    best=0
    for lk in range(ncell_lower):
        free=[coef for coef in P if lk not in bd(coef)]
        if not free: continue
        best=max(best, rank_gfp([list(c) for c in free],m,p))
    return best

print("PHASE B / B2.  T3 OVER GF(p): does the field's characteristic matter?")
print(f"  {'complex':>26s} {'p':>3s} {'m':>3s} {'predicted m-1':>14s} {'measured':>9s}")
print("  "+"-"*62)
# 2x2 patch of plaquettes on a 3x3 vertex grid -> 4 plaquettes, 12 links, oriented incidences
def patch_cells():
    nx=ny=3
    vid={(i,j):j*nx+i for j in range(ny) for i in range(nx)}
    E=[]
    for j in range(ny):
        for i in range(nx-1): E.append(('h',i,j))
    NH=len(E)
    for j in range(ny-1):
        for i in range(nx): E.append(('v',i,j))
    ind={e:k for k,e in enumerate(E)}
    cells=[]
    for j in range(ny-1):
        for i in range(nx-1):
            cells.append({ind[('h',i,j)]:+1, ind[('v',i+1,j)]:+1,
                          ind[('h',i,j+1)]:-1, ind[('v',i,j)]:-1})
    return cells,len(E)
cells,L=patch_cells()
for p in (2,3,5):
    print(f"  {'3x3 patch (4 plaquettes)':>26s} {p:3d} {len(cells):3d} {len(cells)-1:14d} "
          f"{capacity_gfp(cells,L,p):9d}")
tri=[{0:+1,1:+1,2:-1},{0:+1,3:+1,4:-1},{1:+1,3:+1,5:-1}]
for p in (2,3,5):
    print(f"  {'3 triangles, 6 links':>26s} {p:3d} {len(tri):3d} {len(tri)-1:14d} "
          f"{capacity_gfp(tri,6,p):9d}")

print()
print("PHASE B / B3.  T3 IN THREE DIMENSIONS: regions are CUBES, boundaries are FACES.")
def cubes(nx,ny,nz):
    """cubic complex: 3-cells with their 6 bounding faces, signed"""
    def fx(i,j,k): return ('x',i,j,k)
    def fy(i,j,k): return ('y',i,j,k)
    def fz(i,j,k): return ('z',i,j,k)
    F=[]
    for i in range(nx+1):
        for j in range(ny):
            for k in range(nz): F.append(fx(i,j,k))
    for i in range(nx):
        for j in range(ny+1):
            for k in range(nz): F.append(fy(i,j,k))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz+1): F.append(fz(i,j,k))
    ind={f:n for n,f in enumerate(F)}
    cells=[]
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                cells.append({ind[fx(i,j,k)]:-1, ind[fx(i+1,j,k)]:+1,
                              ind[fy(i,j,k)]:-1, ind[fy(i,j+1,k)]:+1,
                              ind[fz(i,j,k)]:-1, ind[fz(i,j,k+1)]:+1})
    return cells,len(F)
print(f"  {'complex':>26s} {'p':>3s} {'m':>3s} {'faces':>6s} {'predicted m-1':>14s} {'measured':>9s}")
print("  "+"-"*72)
for (a,b,c) in ((2,1,1),(2,2,1),(2,2,2)):
    cl,NF=cubes(a,b,c)
    for p in (2,3):
        if p**len(cl) > 20000: continue
        print(f"  {f'{a}x{b}x{c} cubes':>26s} {p:3d} {len(cl):3d} {NF:6d} {len(cl)-1:14d} "
              f"{capacity_gfp(cl,NF,p):9d}")

print()
print("  WHERE THE PROOF STOPS, ASSESSED NOT TESTED:")
print("    T3 needs the boundary map to be LINEAR over a field, which needs an ABELIAN group.")
print("    For a non-abelian group the holonomy of a region depends on ORDERING and BASE POINT,")
print("    so 'the product of the plaquettes' is not a linear function of the region and the")
print("    kernel argument does not apply. WHETHER AN ANALOGUE HOLDS IS UNTESTED AND OPEN.")
