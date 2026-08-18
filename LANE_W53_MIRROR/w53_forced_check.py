"""W-53.  FORCED-OR-NOT CHECK, BEFORE ANY DYNAMICS.

The proposed mirror of W-51: W-51 asked whether geometry can be reconstructed from RECORDS and found
it cannot (the record algebra is identical on block, chain and ring). The mirror asks whether it can
be reconstructed from THE ENVIRONMENT COUPLINGS AND WHICH REGIONS THEY SELECT -- which is where W-51
and W-38 both located it.

W-51's lesson was that a reconstruction can succeed perfectly and mean nothing, because the map was
an exact identity. So the identity check runs FIRST, on paper and then in code, and no dynamics is
built until it passes.

THE SUSPICION, STATED BEFORE TESTING. The selection rule is "a region is protected iff its boundary
avoids the bath". For a single-link bath {L}, the protected set is exactly
        phi({L}) = { S : L not in bd(S) } = ker f_L
and the family { ker f_L : L } determines each functional f_L up to scale, hence every ROW of the
boundary map, hence the incidence structure. If that holds, the mirror test is the SAME identity
wearing different clothes and is void before it starts.

TEST: reconstruct the incidence matrix from the protected-set family ALONE, using no link labels and
no dynamics. If it comes back exactly, the question is circular.
"""
import itertools, numpy as np

def carriers():
    out={}
    nx=ny=3
    vid={(i,j):j*nx+i for j in range(ny) for i in range(nx)}
    E=[]
    for j in range(ny):
        for i in range(nx-1): E.append(('h',i,j))
    for j in range(ny-1):
        for i in range(nx): E.append(('v',i,j))
    ind={e:k for k,e in enumerate(E)}
    PL=[]
    for j in range(ny-1):
        for i in range(nx-1):
            PL.append({ind[('h',i,j)],ind[('v',i+1,j)],ind[('h',i,j+1)],ind[('v',i,j)]})
    out['3x3 block (4 plaq)']=(len(E),PL)
    out['ring of 4 plaquettes']=(8,[{0,1,2,3},{2,3,4,5},{4,5,6,7},{6,7,0,1}])
    out['chain of 4']=(13,[{0,1,2,3},{3,4,5,6},{6,7,8,9},{9,10,11,12}])
    out['bowtie']=(7,[{0,1,2},{2,3,4},{4,5,6},{6,0,1}])
    return out

def bd(S,PL):
    c={}
    for p in S:
        for lk in PL[p]: c[lk]=c.get(lk,0)+1
    return set(lk for lk,v in c.items() if v%2)

def protected_family(L,PL):
    """phi({L}) for every single-link bath: the set of regions whose boundary avoids that link.
       Regions are encoded as bitmasks. NO link labels are exported -- only the SETS."""
    m=len(PL)
    REG=[S for r in range(1,m+1) for S in itertools.combinations(range(m),r)]
    fam=[]
    for lk in range(L):
        fam.append(frozenset(sum(1<<i for i in S) for S in REG if lk not in bd(S,PL)))
    return fam,REG

def reconstruct_incidence(fam,m,REG):
    """From the protected-set family alone, recover for each bath-class which plaquettes it touches.
       f_L(S) = 1 iff L in bd(S). The protected set IS ker f_L, so f_L is recovered by asking which
       single plaquettes are excluded -- no dynamics, no labels, pure set algebra."""
    rows=[]
    for ks in fam:
        row=set()
        for p in range(m):
            if (1<<p) not in ks: row.add(p)      # link L lies on plaquette p
        rows.append(frozenset(row))
    return rows

print("W-53  FORCED-OR-NOT: is the mirror test the same identity as W-51?")
print(f"\n  {'carrier':>22s} {'L':>3s} {'m':>3s} {'incidence recovered exactly?':>30s}")
print("  "+"-"*64)
allok=True
for name,(L,PL) in carriers().items():
    fam,REG=protected_family(L,PL)
    rows=reconstruct_incidence(fam,len(PL),REG)
    truth=[frozenset(p for p in range(len(PL)) if lk in PL[p]) for lk in range(L)]
    ok=(rows==truth)
    allok=allok and ok
    print(f"  {name:>22s} {L:3d} {len(PL):3d} {str(ok):>30s}")
print()
print(f"  incidence recovered exactly on every carrier: {allok}")
print()
print("  DOES IT DISTINGUISH GEOMETRIES? (the property W-51's record algebra LACKED)")
fams={}
for name,(L,PL) in carriers().items():
    fam,_=protected_family(L,PL)
    fams[name]=frozenset(fam)
ks=list(fams)
for i in range(len(ks)):
    for j in range(i+1,len(ks)):
        same = fams[ks[i]]==fams[ks[j]]
        print(f"    {ks[i]:>22s}  vs  {ks[j]:<22s} identical family? {same}")
print()
if allok:
    print("  VERDICT: the protected-set family IS the incidence matrix, row by row, recovered by pure")
    print("  set algebra with no dynamics and no link labels. The mirror test is W-51's identity in")
    print("  different clothing and is VOID BEFORE IT STARTS. No dynamics is built.")
