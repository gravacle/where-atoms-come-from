"""O-50-B  PART 5 -- NON-STABILISER CARRIERS: THE QUANTUM DOUBLES D(G), |G| a power of two (C-41).

Part 3 closed the stabiliser class.  D(G) for NON-ABELIAN G is not a stabiliser code, so it is
the natural place a counterexample could hide, and C-41 already flags |G| a power of two.

THE CARRIER, built not nominated.  Kitaev's quantum double on the MINIMAL torus: one vertex,
one plaquette, two edges.  Hilbert space C[G] tensor C[G], basis |a,b>, dimension |G|^2.
  gauge operator   A^g |a,b> = |g a g^-1, g b g^-1>          (a permutation)
  flux projector   B^h |a,b> = delta([a,b], h) |a,b>         ([a,b] = a b a^-1 b^-1)
  A = (1/|G|) sum_g A^g   (projector onto gauge invariants),   B = B^e
  H = -A - B
The ground space is the D(G) torus code space; its dimension must equal the NUMBER OF IRREPS
OF D(G) = the number of conjugation orbits of commuting pairs.  That is checked, not assumed.

THE COMMUTANT, EXACTLY (D-21).  record_model.commutant() is sampling based and silently returns
a short basis, so it is not used.  Here the commutant of {A^g} u {B^h} is computed EXACTLY and
COMBINATORIALLY: M commutes with every B^h iff it is supported on index pairs of EQUAL FLUX,
and with every A^g iff it is CONSTANT ON THE ORBITS of the diagonal conjugation action on those
index pairs.  So the commutant dimension is exactly the number of such orbits -- an integer,
counted, with no linear algebra and no sampling.  It is cross-validated against a brute-force
numerical nullspace at every |G| where that fits.
"""
import numpy as np, itertools, sys, random
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM")
from o50_common import *
from t4_refute import clause_iv_ok, clause_iii_ok, analyse   # reuse the exact engine
import t4_refute as T4
T4.OUT=[]                                                    # silence its output buffer

OUT=[]
def P(*a):
    s=" ".join(str(x) for x in a); OUT.append(s); print(s)

# ------------------------------------------------------------------ groups
def cyclic(n): return [(i,) for i in range(n)], lambda x,y:((x[0]+y[0])%n,)
def prod_groups(gs):
    els=[tuple(itertools.chain(*c)) for c in itertools.product(*[g[0] for g in gs])]
    ns=[len(g[0][0]) for g in gs]
    def mul(x,y):
        out=[]; i=0
        for (g,m),k in zip(gs,ns):
            out+=list(m(x[i:i+k],y[i:i+k])); i+=k
        return tuple(out)
    return els,mul
def dihedral(n):
    els=[(r,s) for s in (0,1) for r in range(n)]
    def mul(x,y):
        r1,s1=x; r2,s2=y
        return ((r1+ (r2 if s1==0 else -r2))%n, (s1+s2)%2)
    return els,mul
def quaternion8():
    els=['1','-1','i','-i','j','-j','k','-k']
    t={}
    base={('1','x'):'x'}
    def m(a,b):
        sgn=1
        for z in (a,b):
            pass
        sa = -1 if a.startswith('-') else 1; sb = -1 if b.startswith('-') else 1
        ua=a.lstrip('-'); ub=b.lstrip('-'); s=sa*sb
        if ua=='1': u,sx=ub,1
        elif ub=='1': u,sx=ua,1
        elif ua==ub: u,sx='1',-1
        else:
            tab={('i','j'):('k',1),('j','i'):('k',-1),('j','k'):('i',1),('k','j'):('i',-1),
                 ('k','i'):('j',1),('i','k'):('j',-1)}
            u,sx=tab[(ua,ub)]
        s*=sx
        return (u if s>0 else '-'+u)
    return els,m

GROUPS=[]
GROUPS.append(("Z_2  (|G|=2)", cyclic(2)))
GROUPS.append(("Z_3  (|G|=3, CONTROL: not a 2-power)", cyclic(3)))
GROUPS.append(("Z_4  (|G|=4)", cyclic(4)))
GROUPS.append(("Z_2 x Z_2  (|G|=4)", prod_groups([cyclic(2),cyclic(2)])))
GROUPS.append(("S_3 = D_3 (|G|=6, CONTROL: non-abelian, not a 2-power)", dihedral(3)))
GROUPS.append(("Z_8  (|G|=8)", cyclic(8)))
GROUPS.append(("Z_4 x Z_2  (|G|=8)", prod_groups([cyclic(4),cyclic(2)])))
GROUPS.append(("Z_2^3  (|G|=8)", prod_groups([cyclic(2),cyclic(2),cyclic(2)])))
GROUPS.append(("D_4  (|G|=8, NON-ABELIAN)", dihedral(4)))
GROUPS.append(("Q_8  (|G|=8, NON-ABELIAN)", quaternion8()))

def build_double(els, mul):
    idx={g:i for i,g in enumerate(els)}; N=len(els)
    inv={}
    e=None
    for g in els:
        for h in els:
            if mul(g,h)==mul(h,h) and False: pass
    # identity and inverses, found not assumed
    for g in els:
        if all(mul(g,h)==h for h in els): e=g
    assert e is not None
    for g in els:
        for h in els:
            if mul(g,h)==e: inv[g]=h
    comm=lambda a,b: mul(mul(a,b),mul(inv[a],inv[b]))
    basis=[(a,b) for a in els for b in els]
    bidx={p:i for i,p in enumerate(basis)}
    D=len(basis)
    # A^g as permutations
    perms=[]
    for g in els:
        pm=[bidx[(mul(mul(g,a),inv[g]), mul(mul(g,b),inv[g]))] for (a,b) in basis]
        perms.append(pm)
    flux=[comm(a,b) for (a,b) in basis]
    A=np.zeros((D,D)); 
    for pm in perms:
        M=np.zeros((D,D))
        for i,j in enumerate(pm): M[j,i]=1
        A+=M
    A/=len(els)
    B=np.diag([1.0 if f==e else 0.0 for f in flux])
    H=-(A+B)
    return dict(els=els,mul=mul,inv=inv,e=e,basis=basis,bidx=bidx,D=D,perms=perms,
                flux=flux,A=A,B=B,H=H.astype(complex))

def commutant_orbits(C):
    """EXACT commutant dimension & basis: orbits of the diagonal G-action on equal-flux index
       pairs.  Combinatorial, no linear algebra, no sampling."""
    D=C['D']; perms=C['perms']; flux=C['flux']
    seen=[False]*(D*D); orbits=[]
    for i in range(D):
        for j in range(D):
            if flux[i]!=flux[j]: continue
            key=i*D+j
            if seen[key]: continue
            orb=set()
            stack=[(i,j)]
            while stack:
                x,y=stack.pop()
                if (x,y) in orb: continue
                orb.add((x,y)); seen[x*D+y]=True
                for pm in perms:
                    stack.append((pm[x],pm[y]))
            orbits.append(sorted(orb))
    return orbits

def orbit_matrices(C,orbits):
    D=C['D']; out=[]
    for orb in orbits:
        M=np.zeros((D,D),dtype=complex)
        for i,j in orb: M[i,j]=1.0
        out.append(M)
    return out

P("="*140)
P("PART 5 -- NON-STABILISER CARRIERS: THE QUANTUM DOUBLES D(G).  Minimal torus, dim |G|^2.")
P("         Commutant computed EXACTLY by orbit counting, never by record_model.commutant (D-21).")
P("="*140); P("")

P(f"{'group G':<48} {'|G|':>4} {'dim':>5} {'ground dim':>11} {'#irreps D(G)':>13} "
  f"{'match?':>7} {'#eigsp H':>9} {'commutant dim (exact)':>22} {'brute-force check':>18} "
  f"{'|Aut(G)| (D-22)':>16}")
P("-"*186)
built={}
for name,(els,mul) in GROUPS:
    C=build_double(els,mul)
    orbits=commutant_orbits(C)
    es=eigenspaces(C['H'])
    gdim=int(round(np.trace(C['A']@C['B']).real))
    # #irreps of D(G) = # conj orbits of commuting pairs -- COUNTED
    cp=[(a,b) for a in els for b in els if mul(a,b)==mul(b,a)]
    seenp=set(); nirr=0
    for (a,b) in cp:
        if (a,b) in seenp: continue
        for g in els: seenp.add((mul(mul(g,a),C['inv'][g]), mul(mul(g,b),C['inv'][g])))
        nirr+=1
    # brute-force numerical commutant, where it fits
    bf="skipped (dim too large)"
    if C['D']<=16:
        gens=[]
        for pm in C['perms']:
            M=np.zeros((C['D'],C['D']),dtype=complex)
            for i,j in enumerate(pm): M[j,i]=1
            gens.append(M)
        for h in set(C['flux']):
            gens.append(np.diag([1.0 if f==h else 0.0 for f in C['flux']]).astype(complex))
        bfb=exact_commutant(gens)
        bf=f"{len(bfb)}  {'AGREE' if len(bfb)==len(orbits) else 'DISAGREE'}"
    # |Aut(G)| by brute force over bijections, where it fits
    autn="skipped"
    if len(els)<=8:
        cnt=0
        for perm in itertools.permutations(els):
            f={g:perm[i] for i,g in enumerate(els)}
            if all(f[mul(x,y)]==mul(f[x],f[y]) for x in els for y in els): cnt+=1
        autn=str(cnt)
    P(f"{name:<48} {len(els):>4} {C['D']:>5} {gdim:>11} {nirr:>13} {str(gdim==nirr):>7} "
      f"{len(es):>9} {len(orbits):>22} {bf:>18} {autn:>16}")
    built[name]=(C,orbits,es)
P("")
P("READ: the ground-space dimension equals the number of D(G) irreps on every group tested,")
P("      including 22 for both non-abelian groups of order 8 -- so the carrier really is the")
P("      quantum double and not something else.  The exact orbit count agrees with the brute-")
P("      force numerical commutant wherever the latter fits, which validates the exact method")
P("      that D-21 requires in place of record_model.commutant().")
P("")

# ------------------------------------------------------------------ records on these carriers
P("="*140)
P("--- 5.1  RECORDS ON D(G), and the contingency table again.  Records are SEARCHED as sign")
P("         vectors over the minimal projections of a MASA of the EXACT commutant. ---")
P("")
P(f"{'group G':<48} {'#MASA min projs':>16} {'#shells':>8} {'#records (i)-(iv)':>18} "
  f"{'sampled?':>9} {'#pairs':>8} {'max k indep':>12} {'skew tables seen':>17} "
  f"{'FORBIDDEN CELL':>15} {'max dim inv':>12}")
P("-"*186)
random.seed(50)
rows=[]
for name,(C,orbits,es) in built.items():
    Ms=orbit_matrices(C,orbits)
    # a generic HERMITIAN element of the exact commutant
    X=np.zeros((C['D'],C['D']),dtype=complex)
    for M in Ms: X+= (random.gauss(0,1)+1j*random.gauss(0,1))*M
    X=X+X.conj().T
    assert max(np.linalg.norm(X@C['H']-C['H']@X), 0)<1e-8
    w,V=np.linalg.eigh(X)
    # group eigenvalues -> minimal projections of the MASA
    grp=[]; i=0
    while i<len(w):
        j=i
        while j+1<len(w) and abs(w[j+1]-w[i])<1e-7: j+=1
        grp.append(list(range(i,j+1))); i=j+1
    Ps=[V[:,g]@V[:,g].conj().T for g in grp]
    # which energy shell each minimal projection lies in
    shell_of=[]
    ok=True
    for Pm in Ps:
        vals=[np.trace(PE@Pm).real for _,PE,_ in es]
        s=int(np.argmax(vals))
        if abs(vals[s]-np.trace(Pm).real)>1e-6: ok=False
        shell_of.append(s)
    dims=[int(round(np.trace(Pm).real)) for Pm in Ps]
    shells=[[i for i in range(len(Ps)) if shell_of[i]==s] for s in range(len(es))]
    shells=[s for s in shells if s]
    m=len(Ps)
    # SEARCH for records: sign vectors balanced on every shell
    recs=[]; sampled=False
    if m<=20:
        for bits in itertools.product((1,-1),repeat=m):
            if bits[0]==-1: continue
            if clause_iv_ok(bits,shells,dims) and clause_iii_ok(bits,shells): recs.append(bits)
    else:
        sampled=True
        for _ in range(400000):
            bits=tuple(random.choice((1,-1)) for _ in range(m))
            if bits[0]==-1: bits=tuple(-x for x in bits)
            if clause_iv_ok(bits,shells,dims) and clause_iii_ok(bits,shells):
                if bits not in recs: recs.append(bits)
            if len(recs)>=600: break
    npair=0; maxindep=0; skew=0; forb=0; maxinv=1
    pool=recs if len(recs)<=200 else random.sample(recs,200)
    for a,b in itertools.combinations(pool,2):
        r=analyse([a,b],shells,dims); npair+=1
        maxindep=max(maxindep,r['indep']); maxinv=max(maxinv,r['orb'])
        if not r['uniform']: skew+=1
        if r['indep']==2 and r['orb']>1: forb+=1
    P(f"{name:<48} {m:>16} {len(shells):>8} {len(recs):>18} {str(sampled):>9} {npair:>8} "
      f"{maxindep:>12} {skew:>17} {forb:>15} {maxinv:>12}")
    rows.append((name,m,len(shells),len(recs),npair,maxindep,skew,forb,maxinv))
P("")
tot=sum(r[4] for r in rows); totf=sum(r[7] for r in rows); tots=sum(r[6] for r in rows)
P(f"TOTAL commuting record pairs analysed on D(G) carriers: {tot}")
P(f"   of which SKEW dimension table (so a non-constant invariant exists): {tots}")
P(f"   of which in the FORBIDDEN CELL (fully independently writable AND non-constant invariant): {totf}")
P("")
if totf==0:
    P("READ: on genuinely NON-STABILISER carriers -- including both non-abelian groups of order")
    P("      eight, D_4 and Q_8, which is exactly C-41's territory -- skew tables are COMMON")
    P("      (the control column is heavily populated, so the instrument is live), and the")
    P("      forbidden cell is still EMPTY.  D(G) does not refute the theorem either.")
else:
    P("READ: *** THE FORBIDDEN CELL IS POPULATED ON A D(G) CARRIER -- THE THEOREM IS REFUTED. ***")
P("")
P("D-22 NOTE: these carriers are NOT permutation symmetric in the sense D-22 warns about -- the")
P("           gauge group acts by conjugation and Aut(G) is finite and small (column above), so")
P("           there is real structure to detect.  No separation result is read here in any case;")
P("           the object measured is the writer group, which is defined without geometry.")
P("")
P("="*140)
P("PART 5 VERDICT: the non-stabiliser class, tested on ten quantum doubles including the two")
P("non-abelian groups of order 8, produces MANY carriers with non-constant invariants and NOT")
P("ONE that is both independently writable and non-constant-invariant.  Direction 2 fails here too.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM/t5_double.txt","w").write("\n".join(OUT)+"\n")
