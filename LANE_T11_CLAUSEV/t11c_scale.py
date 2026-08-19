"""T-11c: is clause (v)'s failure a SMALL-LATTICE ARTEFACT?

At L=2 the carrier has 8 edges and d=2, and a 'region' in the geometric test spans up to 8 of
them -- the whole system. An operation on the whole system flipping the record is not a
counterexample to protection; it is what protection never claimed.

THE TEST: hold the region SMALL relative to the carrier and grow the carrier. If the admissible
count goes to zero as regions become genuinely local, clause (v) is fine and L=2 was too small."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); Y=1j*(X@Z)
P4=[I2,X,Y,Z]
def build(n):
    ind={}; E=[]; k=0
    for j in range(n):
        for i in range(n):
            ind[('h',i,j)]=k; E.append((j*n+i,j*n+(i+1)%n)); k+=1
            ind[('v',i,j)]=k; E.append((j*n+i,((j+1)%n)*n+i)); k+=1
    PLQ=[[ind[('h',i,j)],ind[('v',(i+1)%n,j)],ind[('h',i,(j+1)%n)],ind[('v',i,j)]] for j in range(n) for i in range(n)]
    STAR=[[ind[('h',i,j)],ind[('h',(i-1)%n,j)],ind[('v',i,j)],ind[('v',i,(j-1)%n)]] for j in range(n) for i in range(n)]
    return len(E),ind,E,PLQ,STAR
say("="*96); say("T-11c   IS CLAUSE (v)'s FAILURE A SMALL-LATTICE ARTEFACT?"); say("="*96)
say(f"  {'lattice':<10}{'edges':>7}{'d':>4}{'region':>9}{'region/system':>15}{'any-unitary':>14}{'ADMISSIBLE':>13}")
for n in (2,3):
    L,ind,E,PLQ,STAR = build(n)
    if L>18: continue
    # operators are built only on the region, so the Hilbert space is the region's -- we test
    # the SYMPLECTIC condition instead of building 2^L matrices: an operator flips the Z-type
    # logical iff it anticommutes with it, and is admissible iff it commutes with every stabiliser.
    logical = [ind[('h',i,0)] for i in range(n)]           # a Z-type wrap: the record
    stabs = [('X',s) for s in STAR] + [('Z',p) for p in PLQ]
    def anticomm(kind, supp, other_kind, other_supp):
        if kind==other_kind: return False
        return len(set(supp)&set(other_supp)) % 2 == 1
    for rsize in (1, n):                                   # ONE plaquette, then n of them
        regions=[list(s) for s in itertools.combinations(range(len(PLQ)), rsize)]
        anyf=admf=0
        for rg in regions:
            supp=sorted({e for p in rg for e in PLQ[p]})
            for combo in itertools.product(range(4), repeat=len(supp)):
                if all(c==0 for c in combo): continue
                xs=[l for l,c in zip(supp,combo) if c in (1,2)]
                zs=[l for l,c in zip(supp,combo) if c in (2,3)]
                # flips a Z-type logical iff its X-part overlaps the logical oddly
                if len(set(xs)&set(logical))%2!=1: continue
                anyf+=1
                ok=True
                for kind,ss in stabs:
                    ov = len(set(zs)&set(ss))%2 if kind=='X' else len(set(xs)&set(ss))%2
                    if ov: ok=False; break
                if ok: admf+=1
        frac=len(supp)/L
        say(f"  {'torus '+str(n)+'x'+str(n):<10}{L:>7}{n:>4}{rsize:>9}{frac:>15.2f}{anyf:>14}{admf:>13}")
say("")
say("  READ: if the admissible count falls to 0 once a region is a small FRACTION of the carrier,")
say("  clause (v) is sound and L=2 simply has no room for a region to be local.")
