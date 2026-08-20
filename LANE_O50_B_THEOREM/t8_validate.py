"""O-50-B  PART 8 -- VALIDATION OF THE REDUCTION, and the final statement.

Part 4's exhaustive hunt rests on one reduction: that ANY mutually commuting family of records
on ANY carrier is a family of SIGN VECTORS over the minimal projections of a MASA of the
commutant, so that searching profiles searches carriers.  Part 8 tests that reduction against
DENSE OPERATORS on randomly generated carriers -- random Hamiltonians with prescribed
degeneracies and random commuting Lindblad operators -- by computing the writer group twice,
once combinatorially and once by BUILDING the unitary and verifying it as an operator.
"""
import numpy as np, itertools, sys, random
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM")
from o50_common import *
import t4_refute as T4
T4.OUT=[]
from t4_refute import clause_iv_ok, clause_iii_ok, analyse

OUT=[]
def P(*a):
    s=" ".join(str(x) for x in a); OUT.append(s); print(s)

P("="*140)
P("PART 8 -- VALIDATION OF THE REDUCTION ON RANDOM DENSE CARRIERS (non-stabiliser, non-double)")
P("="*140); P("")

rng=np.random.default_rng(8050)
random.seed(8050)
P(f"{'carrier':<44} {'dim':>5} {'shell dims':<20} {'#MASA projs':>12} {'#records':>9} {'#families':>10} "
  f"{'combinatorial |G_W|':>20} {'operator-verified |G_W|':>24} {'agree':>6} {'forbidden cell':>15}")
P("-"*186)
tot_f=0; tot_forb=0; agree_all=True
for trial in range(14):
    # a random carrier: random unitary basis, prescribed degeneracies, random commuting Ls
    degs=random.choice([[4,4],[6,2],[8],[4,4,4],[2,2,4],[6,6],[8,4],[4,2,2,4]])
    n=sum(degs)
    Q,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    evals=[]; 
    for i,dg in enumerate(degs): evals += [float(i)]*dg
    H=Q@np.diag(evals).astype(complex)@Q.conj().T
    H=(H+H.conj().T)/2
    es=eigenspaces(H)
    # a random HERMITIAN element of the commutant of H, built exactly by shell projection
    Xr=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); Xr=Xr+Xr.conj().T
    Xc=sum(PE@Xr@PE for _,PE,_ in es)
    w,V=np.linalg.eigh(Xc)
    grp=[]; i=0
    while i<len(w):
        j=i
        while j+1<len(w) and abs(w[j+1]-w[i])<1e-8: j+=1
        grp.append(list(range(i,j+1))); i=j+1
    Ps=[V[:,g]@V[:,g].conj().T for g in grp]
    dims=[int(round(np.trace(Pm).real)) for Pm in Ps]
    shell_of=[]
    for Pm in Ps:
        vals=[np.trace(PE@Pm).real for _,PE,_ in es]; shell_of.append(int(np.argmax(vals)))
    shells=[[i for i in range(len(Ps)) if shell_of[i]==s] for s in range(len(es))]
    shells=[s for s in shells if s]
    m=len(Ps)
    if m>18: continue
    recs=[b for b in itertools.product((1,-1),repeat=m)
          if b[0]==1 and clause_iv_ok(b,shells,dims) and clause_iii_ok(b,shells)]
    if len(recs)<2:
        P(f"{'trial '+str(trial)+'  degs='+str(degs):<44} {n:>5} {str([m2 for _,_,m2 in es]):<20} "
          f"{m:>12} {len(recs):>9} {0:>10} {'--':>20} {'--':>24} {'n/a':>6} {'n/a':>15}")
        continue
    pool=recs if len(recs)<=25 else random.sample(recs,25)
    nf=0; nforb=0; agree=True; gc=[]; go=[]
    for fam in itertools.combinations(pool,2):
        r=analyse(list(fam),shells,dims); nf+=1
        Rs=[sum(s*Pm for s,Pm in zip(f,Ps)) for f in fam]
        for R in Rs:
            assert clause_i(R) and clause_ii(R,H) and clause_iv_(R,es)
        es2,blocks,dtab,cfg=dim_table(H,Rs)
        # OPERATOR-VERIFIED writer group: build U for each eps and check it as an operator
        gver=[]
        for mm in range(4):
            eps=tuple((mm>>i)&1 for i in range(2))
            U=build_flip_unitary(blocks,eps,len(es2),cfg)
            if U is None: continue
            ok = (np.linalg.norm(U.conj().T@U-np.eye(n))<1e-8 and np.linalg.norm(U@H-H@U)<1e-7)
            for i,R in enumerate(Rs):
                sgn=-1.0 if eps[i] else 1.0
                ok = ok and np.linalg.norm(U.conj().T@R@U - sgn*R)<1e-7
            if ok: gver.append(eps)
        if sorted(gver)!=sorted(r['G']): agree=False; agree_all=False
        gc.append(len(r['G'])); go.append(len(gver))
        if r['indep']==2 and r['orb']>1: nforb+=1
    tot_f+=nf; tot_forb+=nforb
    P(f"{'trial '+str(trial)+'  degs='+str(degs):<44} {n:>5} {str([m2 for _,_,m2 in es]):<20} "
      f"{m:>12} {len(recs):>9} {nf:>10} {str(sorted(set(gc))):>20} {str(sorted(set(go))):>24} "
      f"{str(agree):>6} {nforb:>15}")
P("")
P(f"TOTAL families operator-verified on random dense carriers: {tot_f}")
P(f"combinatorial writer group == operator-verified writer group on every family: {agree_all}")
P(f"FORBIDDEN CELL population on random dense carriers: {tot_forb}")
P("")
P("READ: on random carriers that are neither stabiliser codes nor quantum doubles, the writer")
P("      group computed from the dimension table alone agrees EXACTLY with the group verified")
P("      by explicitly constructing each unitary and testing three operator norms.  The")
P("      reduction Part 4 rests on is therefore sound on carriers built with no structure at")
P("      all, and the forbidden cell is empty here too.")
P("")
P("D-22: these carriers are generated from a HAAR-random basis, so their automorphism group is")
P("      trivial -- there is no permutation symmetry to hide behind, and equally no geometry to")
P("      detect.  Nothing geometric is read from them; only the writer group, which needs none.")
P("")

# ------------------------------------------------------------------ final statement
P("="*140)
P("THE THEOREM, AS IT SURVIVES O-50-B")
P("="*140)
P("")
P("SETTING.  A carrier (H, {L_k}).  ADMISSIBLE U := unitary with [U,H] = 0 (O-4).  R_1..R_k a")
P("family of mutually commuting RECORDS -- each satisfying clauses (i)-(iv) ON ITS OWN, with NO")
P("even-splitting demand.  V(E,sigma) their joint eigenspaces inside energy shell E, and")
P("d(E,sigma) = dim V(E,sigma).  G_W := the group of flip patterns eps realisable by an")
P("admissible U with U* R_i U = (-1)^{eps_i} R_i for every i.")
P("")
P("LEMMA W  (proved, Part 2).  eps in G_W  <=>  d(E,sigma) = d(E,eps.sigma) for all E, sigma.")
P("          Both directions.  A negative verdict is a NON-EXISTENCE proof over all unitaries.")
P("")
P("LEMMA B  (proved, Part 2).  If chi_S is not G_W-invariant then Tr(P_E R_S) = 0 on every")
P("          shell: every non-invariant character has EXACTLY zero shell average.")
P("")
P("THEOREM 1  (proved, Parts 1-2).  If the records are INDEPENDENTLY WRITABLE then G_W contains")
P("          every coordinate flip, acts SIMPLY TRANSITIVELY on the 2^k configurations, and the")
P("          space of G_W-invariant functionals is EXACTLY ONE-DIMENSIONAL -- the constants.")
P("")
P("THEOREM 2 -- THE CANCELLATION LAW  (proved, Parts 1, 2, 6; NO hypothesis of independent")
P("          writability).  Every functional f of the record configuration splits uniquely as")
P("               f  =  Pi_{G_W} f   +   (f - Pi_{G_W} f)")
P("          where Pi_{G_W} f is UNRESPONSIVE (no admissible family-preserving write changes it)")
P("          and f - Pi_{G_W} f has EXACTLY ZERO shell average (it cancels).  Hence NO functional")
P("          of the record configuration is both RESPONSIVE TO WRITING and NON-CANCELLING.")
P("")
P("THEOREM D  (proved, Part 6).  k independently writable records force 2^k | dim E on every")
P("          eigenspace; so k_indep <= min_E v_2(dim E), a bound read off the SPECTRUM alone.")
P("")
P("TRADE LAW  (exact identity, Part 6).  dim(invariants) = 2^(k - rank G_W).")
P("")
P("WHAT IS FALSE.  The candidate's literal second clause -- 'every non-constant functional is")
P("G_W-odd, hence has mean exactly zero' -- is REFUTED by f(sigma) = #{i : sigma_i = +1}, which")
P("is non-constant, responsive, and has mean k/2.  Non-constant does not imply odd.  The split")
P("in Theorem 2 is the correct repair.")
P("")
P("WHAT IS NOT PROVED, and is the honest open edge.  G_W is built from admissible unitaries that")
P("PRESERVE the family.  Clause (iv) gives each record a writer, but that writer need not")
P("normalise the family: at k >= 3 there are families in which NO admissible unitary flips any")
P("one record while fixing the others (Part 7), so the flip that clause (iv) guarantees must")
P("SCRAMBLE the rest.  After such a write there is no record configuration left for a functional")
P("to be a functional of -- which is why Theorem 2 is stated for family-preserving writes.  A")
P("source built on family-scrambling writes is not excluded by anything proved here.  That is")
P("the one route this lane did not close.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM/t8_validate.txt","w").write("\n".join(OUT)+"\n")
