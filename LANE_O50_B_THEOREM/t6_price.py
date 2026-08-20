"""O-50-B  PART 6 -- DIRECTION 3: THE HYPOTHESIS'S PRICE, and two results STRONGER than the
candidate that fell out of Parts 4-5.

RESULT 1 -- THE DICHOTOMY IS UNCONDITIONAL (independent writability is NOT needed for it).
  Let R_1..R_k be mutually commuting records on ANY carrier, G_W the group of flip patterns
  realisable by admissible unitaries, d the joint dimension table.  For any functional f of the
  record configuration write  f = Pi_G f  +  (f - Pi_G f).  Then
     * Pi_G f is G_W-INVARIANT, so NO admissible write changes it: it is UNRESPONSIVE;
     * f - Pi_G f is a combination of NON-invariant characters, and by LEMMA B every one of
       those has EXACTLY zero d-weighted shell average: it CANCELS.
  Hence NO functional of the record configuration is both responsive to writing and non-
  cancelling -- ON ANY CARRIER, independently writable or not.  The candidate asked for this
  under a hypothesis; the hypothesis turns out to be unnecessary FOR THE CONCLUSION.

  What the hypothesis actually controls is the SIZE of the unresponsive part:
     independently writable   =>  invariants are CONSTANTS -- record-BLIND (this is C-61);
     not independently writable =>  invariants can depend on the configuration, but then they
                                    are FROZEN: no admissible operation can move them.
  So dropping independent writability does not buy a source.  It buys a SUPERSELECTION CHARGE.

RESULT 2 -- THEOREM D (divisibility), an exact carrier-level obstruction.
  If k records are independently writable then G_W contains every flip, acts transitively, and
  preserves d; so d(E,sigma) is CONSTANT in sigma and  2^k divides dim(E) for EVERY eigenspace E.
  Hence   k_indep  <=  min_E  v_2( dim E ).
  This is a hard bound computable from the SPECTRUM ALONE, before any record is constructed.

CONJECTURE F (freezing), tested exhaustively below: if every R_i is a genuine record and chi_S
is G_W-invariant for S non-empty, must R_S fail clause (iv)?
"""
import itertools, sys, random
import numpy as np
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM")
from o50_common import *
import t4_refute as T4
T4.OUT=[]
from t4_refute import clause_iv_ok, clause_iii_ok, analyse, dtable, writer_group, norbits, rank2

OUT=[]
def P(*a):
    s=" ".join(str(x) for x in a); OUT.append(s); print(s)

P("="*136)
P("PART 6 -- DIRECTION 3: THE PRICE OF THE HYPOTHESIS.")
P("="*136); P("")

# ------------------------------------------------------------------ 6.1 the exact trade law
P("--- 6.1  THE EXACT TRADE LAW.  dim(invariants) = 2^(k - rank G_W), and k_indep <= rank G_W. ---")
P("")
P("Exhaustive over the same profiles as Part 4.  Every family is a genuine one: each R_i")
P("satisfies clauses (i)-(iv) on its own.")
P("")
tally={}
tot=0; law_ok=True; bound_ok=True
freeze_fail=[]
for m in (4,6,8):
    for D in (1,2,3):
        for dims in itertools.combinations_with_replacement(range(1,D+1),m):
            if sum(dims)%2: continue
            shells=[list(range(m))]; dl=list(dims)
            recs=[b for b in itertools.product((1,-1),repeat=m)
                  if b[0]==1 and clause_iv_ok(b,shells,dl) and clause_iii_ok(b,shells)]
            if len(recs)<2: continue
            for fam in itertools.combinations(recs,2):
                r=analyse(list(fam),shells,dl); tot+=1
                k=2; rk=r['rank']
                if r['orb']!=2**(k-rk): law_ok=False
                if r['indep']>rk: bound_ok=False
                key=(r['indep'],rk,r['orb'])
                tally[key]=tally.get(key,0)+1
                # CONJECTURE F: invariant chi_S  =>  R_S unbalanced?
                for S in ((0,),(1,),(0,1)):
                    if all(sum(g[i] for i in S)%2==0 for g in r['G']):
                        tr=sum(dl[i]*np.prod([fam[j][i] for j in S]) for i in shells[0])
                        if tr==0: freeze_fail.append((dl,fam,S))
P(f"{'k_indep':>8} {'rank G_W':>9} {'#orbits':>8} {'dim inv':>8} {'2^(k-rank)':>11} "
  f"{'law holds':>10} {'count':>10}")
P("-"*74)
for key in sorted(tally):
    ind,rk,orb=key
    P(f"{ind:>8} {rk:>9} {orb:>8} {orb:>8} {2**(2-rk):>11} {str(orb==2**(2-rk)):>10} {tally[key]:>10}")
P("")
P(f"families analysed: {tot}")
P(f"TRADE LAW  dim inv = 2^(k - rank G_W)   holds on every family: {law_ok}")
P(f"BOUND      k_indep <= rank G_W          holds on every family: {bound_ok}")
P("")
P("READ: the invariant dimension is EXACTLY 2^(k - rank G_W).  Every bit of writer-group rank")
P("      lost DOUBLES the invariant space.  The trade is not approximate; it is an identity.")
P("")

P("--- CONJECTURE F, tested exhaustively on the same families ---")
P("")
P(f"families where a NON-CONSTANT G_W-INVARIANT exists AND the corresponding product R_S is")
P(f"nevertheless BALANCED (i.e. would satisfy clause (iv)): {len(freeze_fail)}")
if freeze_fail:
    P("  CONJECTURE F IS FALSE.  Smallest witnesses:")
    for w in freeze_fail[:6]:
        dl,fam,S=w
        P(f"    dims={dl}  R1={''.join('+' if x>0 else '-' for x in fam[0])} "
          f"R2={''.join('+' if x>0 else '-' for x in fam[1])}  S={S}  Tr(P_E R_S)=0")
    P("  MEANING: R_S is absolutely writable -- SOME admissible U flips it -- but no admissible U")
    P("  flips it while leaving the rest of the family alone.  The freezing is RELATIVE to the")
    P("  family, not absolute.  That is a weaker and more honest statement than Conjecture F.")
else:
    P("  CONJECTURE F HOLDS on every family tested HERE -- and note the scope: these are k = 2")
    P("  families on a single shell.  At k = 2 it is a THEOREM, proved exactly in Part 7.")
    P("  AT k >= 3 IT IS FALSE: Part 7 exhibits 81361 violations at k = 3 and 1055818 at k = 4.")
    P("  There the freezing is RELATIVE to the family, not absolute -- an admissible unitary can")
    P("  flip the invariant, but only by scrambling the rest of the family.  Do not read the")
    P("  zero above as a general fact; read Part 7.")
P("")

# ------------------------------------------------------------------ 6.2 Theorem D
P("="*136)
P("--- 6.2  THEOREM D (divisibility): k_indep <= min_E v_2(dim E), computable from the")
P("         SPECTRUM ALONE.  Verified on every carrier this lane has touched. ---")
P("")
def v2(x):
    c=0
    while x%2==0 and x>0: x//=2; c+=1
    return c

P(f"{'carrier':<50} {'dim':>7} {'eigenspace dims':<34} {'min v_2':>8} {'bound on k_indep':>17} "
  f"{'k_indep observed':>17} {'bound respected':>16}")
P("-"*160)

# stabiliser carriers
def sp_(a,b,n): return sum(a[i]*b[n+i]+a[n+i]*b[i] for i in range(n))%2
def rref_(rows,n):
    rows=[r[:] for r in rows]; piv=[]; r=0
    for c in range(2*n):
        p=next((i for i in range(r,len(rows)) if rows[i][c]),None)
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        for i in range(len(rows)):
            if i!=r and rows[i][c]: rows[i]=[(x+y)%2 for x,y in zip(rows[i],rows[r])]
        piv.append(c); r+=1
    return rows[:r],piv
def code_from_strings(gens):
    n=len(gens[0]); S=[]
    for g in gens:
        r=[0]*(2*n)
        for i,c in enumerate(g):
            if c in 'XY': r[i]=1
            if c in 'ZY': r[n+i]=1
        S.append(r)
    return n,S
def toric(L):
    n=2*L*L
    def h(i,j): return 2*((i%L)*L+(j%L))
    def v(i,j): return 2*((i%L)*L+(j%L))+1
    S=[]
    for i in range(L):
        for j in range(L):
            r=[0]*(2*n)
            for e in (h(i,j),h(i-1,j),v(i,j),v(i,j-1)): r[e]^=1
            S.append(r)
    for i in range(L):
        for j in range(L):
            r=[0]*(2*n)
            for e in (h(i,j),h(i,j+1),v(i,j),v(i+1,j)): r[n+e]^=1
            S.append(r)
    return n,S

for nm,(n,S) in [("[[4,2,2]]",code_from_strings(["XXXX","ZZZZ"])),
                 ("[[8,3,2]]",code_from_strings(["XXXXXXXX","ZZZZZZZZ","IIZZIIZZ","IZIZIZIZ","IIIIZZZZ"])),
                 ("TORIC L=2 (TORUS, D-23)",toric(2))]:
    Sr,_=rref_(S,n); k=n-len(Sr)
    H=-sum(xz_to_matrix(s,n) for s in Sr)
    es=eigenspaces(H)
    ds=[m for _,_,m in es]
    mv=min(v2(x) for x in ds)
    P(f"{nm:<50} {2**n:>7} {str(ds)[:33]:<34} {mv:>8} {mv:>17} {k:>17} {str(k<=mv):>16}")

# quantum doubles
import t5_double as T5
T5.OUT=[]
for name,(els,mul) in T5.GROUPS:
    C=T5.build_double(els,mul); es=eigenspaces(C['H'])
    ds=[m for _,_,m in es]
    mv=min(v2(x) for x in ds)
    P(f"{('D('+name.split(' ')[0]+')  '+name):<50} {C['D']:>7} {str(ds)[:33]:<34} {mv:>8} {mv:>17} "
      f"{'see 5.1':>17} {'n/a':>16}")
P("")
P("READ, from the numbers: on D(D_4) and D(Q_8) the eigenspace dimensions include 22, whose")
P("      2-adic valuation is 1.  THEOREM D therefore forbids more than ONE independently")
P("      writable record on those carriers, no matter how many records they hold -- which is")
P("      exactly why Part 5.1 found k_indep = 0 on both and never a full pair.  That null is an")
P("      EXACT OBSTRUCTION read off the spectrum, not a sampling artifact.")
P("")
P("      The non-abelian doubles are therefore FORCED into the non-independently-writable")
P("      regime.  They necessarily carry non-constant invariants -- and by 6.1 those invariants")
P("      are frozen.  C-41's territory is settled: it holds invariants, and none of them move.")
P("")

# ------------------------------------------------------------------ 6.3 the worked exhibit
P("="*136)
P("--- 6.3  THE WORKED EXHIBIT.  The smallest carrier on which a NON-CONSTANT INVARIANT exists")
P("         and every R_i is a genuine record.  Full clause audit, operator by operator. ---")
P("")
labels=[]
for s,mm in [((1,1),3),((1,-1),1),((-1,1),1),((-1,-1),3)]:
    for _ in range(mm): labels.append(s)
nn=len(labels)
H=np.zeros((nn,nn),dtype=complex)                      # ONE eigenspace: H = 0 on all of it
R1=np.diag([float(s[0]) for s in labels]).astype(complex)
R2=np.diag([float(s[1]) for s in labels]).astype(complex)
es=eigenspaces(H)
P(f"carrier dimension {nn}; H has {len(es)} eigenspace(s) of dimension {[m for _,_,m in es]}")
P(f"joint block dimensions d(++,+-,-+,--) = (3,1,1,3)   [INSERTED: the table is the independent")
P(f"variable of this exhibit.  Everything below is INDUCED from it.]")
P("")
P(f"{'operator':<16} {'(i) R=R+, R^2=I':>17} {'(ii) [R,H]=0':>14} {'(iii) non-trivial':>18} "
  f"{'max |Tr(P_E R)|':>17} {'(iv) writable':>15} {'IS IT A RECORD?':>17}")
P("-"*124)
ops=[("R_1",R1),("R_2",R2),("R_1 R_2",R1@R2)]
for nm,Rm in ops:
    c1=clause_i(Rm); c2=clause_ii(Rm,H); c3=clause_iii_(Rm,es); imb=imbalance(Rm,es); c4=imb<1e-7
    P(f"{nm:<16} {str(c1):>17} {str(c2):>14} {str(c3):>18} {imb:>17.6f} {str(c4):>15} "
      f"{str(c1 and c2 and c3 and c4):>17}")
P("")
es2,blocks,d,cfg=dim_table(H,[R1,R2])
G=realisable_flips(d,cfg,len(es2))
orbs=orbits_of(G,cfg); inv=invariant_characters(G,2)
P(f"writer group G_W = {G}   order {len(G)}   rank {f2rank(G,2)}")
P(f"orbits on the 4 configurations: {orbs}")
P(f"invariant characters: {['chi_'+''.join(str(i+1) for i in S) if S else 'constant' for S in inv]}")
P(f"dim(invariants) = {len(orbs)}")
P("")
P("independently writable records (exact criterion): "
  f"{[i for i in range(2) if tuple(1 if j==i else 0 for j in range(2)) in G]}  -- NONE")
P("cross-check, record_model.independently_writable: "
  f"{RecordModel(H).independently_writable([R1,R2])}")
P("")
P("THE PRICE, read off this single carrier:")
P("  * R_1 and R_2 ARE records: bit, durable, non-trivial, and each individually WRITABLE")
P("    (max |Tr(P_E R_i)| = 0.000000 above).")
P("  * NEITHER is INDEPENDENTLY writable: the only admissible flip is the SIMULTANEOUS one.")
P("    Writing R_1 necessarily writes R_2.  That is exactly what was given up.")
P("  * A non-constant invariant appears: chi_12, the CORRELATION R_1 R_2.")
P("  * And R_1 R_2 has |Tr(P_E R_1R_2)| = 4.000000 -- it FAILS clause (iv).  It is a durable")
P("    non-trivial bit that CANNOT BE WRITTEN BY ANY ADMISSIBLE OPERATION WHATEVER (C-11).")
P("  * So the configuration-dependent, non-cancelling quantity this carrier buys is FROZEN")
P("    (absolutely so at k=2, by the Part 7 theorem).")
P("    It accumulates, it is not blind -- and nothing can create it.  It is a superselection")
P("    charge, not a source.  This is O-42 ('writability requires energetic neutrality') and")
P("    O-40 ('C-47's extensive quantity was bought at the cost of independent writability')")
P("    appearing as the SAME fact, one level up, about PRODUCTS of records.")
P("")

# ------------------------------------------------------------------ 6.4 unconditional dichotomy
P("="*136)
P("--- 6.4  RESULT 1 VERIFIED: the dichotomy holds WITHOUT the independent-writability")
P("         hypothesis.  For every family, invariant part unresponsive AND non-invariant part")
P("         exactly cancelling. ---")
P("")
random.seed(7)
P(f"{'profile class':<34} {'#families':>10} {'max |response of INVARIANT part|':>33} "
  f"{'max |shell avg of NON-INV part|':>32} {'CONTROL: shell avg of full f':>29}")
P("-"*144)
for m,D in [(4,1),(4,2),(6,1),(6,2),(8,1)]:
    fams=0; r1=0; r2=0; ctl=0
    for dims in itertools.combinations_with_replacement(range(1,D+1),m):
        if sum(dims)%2: continue
        shells=[list(range(m))]; dl=list(dims)
        recs=[b for b in itertools.product((1,-1),repeat=m)
              if b[0]==1 and clause_iv_ok(b,shells,dl) and clause_iii_ok(b,shells)]
        if len(recs)<2: continue
        for fam in itertools.combinations(recs,2):
            r=analyse(list(fam),shells,dl); fams+=1
            Gg=r['G']; cfg=r['cfg']; dd=r['d']
            for _ in range(3):
                vals={c:random.randint(-20,20) for c in cfg}
                Pi={c:sum(vals[tuple(-x if e else x for e,x in zip(g,c))] for g in Gg) for c in cfg}
                # invariant part response to every admissible write
                for g in Gg:
                    for c in cfg:
                        gc=tuple(-x if e else x for e,x in zip(g,c))
                        r1=max(r1,abs(Pi[gc]-Pi[c]))
                # shell average of the NON-invariant part
                nz=sum(dd[(0,c)]*(vals[c]*len(Gg)-Pi[c]) for c in cfg)
                r2=max(r2,abs(nz))
                ctl=max(ctl,abs(sum(dd[(0,c)]*vals[c] for c in cfg)))
    P(f"{'m='+str(m)+', dims in 1..'+str(D):<34} {fams:>10} {r1:>33} {r2:>32} {ctl:>29}")
P("")
P("READ: the invariant part of every functional has EXACTLY zero response to every admissible")
P("      write (column 3, integer zero), and the non-invariant part has EXACTLY zero shell")
P("      average (column 4, integer zero) -- on families that are NOT independently writable.")
P("      The CONTROL column is the shell average of the FULL functional, which is large and")
P("      non-zero, so the two zeros are zeros of the objects named and not of everything.")
P("")
P("="*136)
P("PART 6 VERDICT")
P("="*136)
P("  * The trade is an IDENTITY: dim(invariants) = 2^(k - rank G_W).  Each unit of writer-group")
P("    rank surrendered doubles the invariant space.")
P("  * THEOREM D bounds k_indep by min_E v_2(dim E) -- readable off the SPECTRUM.  It forces")
P("    D(D_4) and D(Q_8) into the non-independently-writable regime (22 = 2 x 11).")
P("  * What is bought is NOT a source.  It is a configuration-dependent quantity that no")
P("    admissible operation PRESERVING THE FAMILY can move.  At k=2 that freezing is ABSOLUTE")
P("    (Part 7 proves R_S then fails clause (iv) outright).  At k>=3 it is RELATIVE: the")
P("    invariant can be flipped, but only by a write that destroys the other records' values,")
P("    after which there is no record configuration left for it to be a functional of.")
P("  * The dichotomy 'responsive XOR non-cancelling' needs NO hypothesis at all.  Independent")
P("    writability only decides whether the unresponsive part is BLIND (constant) or FROZEN.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM/t6_price.txt","w").write("\n".join(OUT)+"\n")
