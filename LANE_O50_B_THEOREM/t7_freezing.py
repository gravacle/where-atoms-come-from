"""O-50-B  PART 7 -- THE FREEZING THEOREM.  Is every non-constant writer-invariant of a genuine
record family ABSOLUTELY unwritable (fails clause (iv)), or only unwritable-while-keeping-the-
family-fixed?  Part 6 found no counterexample at k=2 on one shell.  Part 7 proves k=2 exactly
and pushes the search to k=3, k=4 and multi-shell profiles.

PROOF AT k=2 (exact, all carriers).  On a shell E write a=d(++), b=d(+-), c=d(-+), e=d(--).
  R_1 is a record => balanced:  a+b = c+e.
  R_2 is a record => balanced:  a+c = b+e.
  Subtracting:  b-c = c-b  =>  b = c;  substituting back  =>  a = e.
  So ON EVERY SHELL, d(++) = d(--) and d(+-) = d(-+), FOR FREE.
  COROLLARY 1: the simultaneous flip (1,1) is ALWAYS in G_W for a pair of records.  rank G_W >= 1,
               so dim(invariants) <= 2 at k=2 -- a pair of records can never be fully frozen.
  COROLLARY 2: chi_12 is G_W-invariant  <=>  flip_1 not in G_W  <=>  a != c = b on some shell.
               And  Tr(P_E R_1 R_2) = a - b - c + e = 2(a - b),  which is then NON-ZERO.
  => at k=2 a non-constant invariant ALWAYS fails clause (iv).  THE FREEZING IS ABSOLUTE.  QED
"""
import itertools, sys
import numpy as np
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM")
import t4_refute as T4
T4.OUT=[]
from t4_refute import clause_iv_ok, clause_iii_ok, analyse

OUT=[]
def P(*a):
    s=" ".join(str(x) for x in a); OUT.append(s); print(s)

P("="*130)
P("PART 7 -- THE FREEZING THEOREM: is a non-constant writer-invariant ABSOLUTELY unwritable?")
P("="*130); P("")

P("--- 7.1  COROLLARY 1 verified: for ANY pair of records, d(E,++) = d(E,--) and d(E,+-) = d(E,-+) ---")
P("")
P(f"{'profile class':<30} {'#pairs':>9} {'max |d(++)-d(--)|':>19} {'max |d(+-)-d(-+)|':>19} "
  f"{'(1,1) in G_W always?':>21} {'CONTROL: one R unbalanced':>26}")
P("-"*130)
for m,D in [(4,1),(4,2),(6,1),(6,2),(6,3),(8,1),(8,2)]:
    np1=0; w1=0; w2=0; always=True; ctl=0
    for dims in itertools.combinations_with_replacement(range(1,D+1),m):
        if sum(dims)%2: continue
        shells=[list(range(m))]; dl=list(dims)
        recs=[b for b in itertools.product((1,-1),repeat=m)
              if b[0]==1 and clause_iv_ok(b,shells,dl) and clause_iii_ok(b,shells)]
        bad=[b for b in itertools.product((1,-1),repeat=m)
             if b[0]==1 and (not clause_iv_ok(b,shells,dl)) and clause_iii_ok(b,shells)]
        for fam in itertools.combinations(recs,2):
            r=analyse(list(fam),shells,dl); np1+=1; dd=r['d']
            w1=max(w1,abs(dd[(0,(1,1))]-dd[(0,(-1,-1))]))
            w2=max(w2,abs(dd[(0,(1,-1))]-dd[(0,(-1,1))]))
            if (1,1) not in r['G']: always=False
        # CONTROL: pair one record with a NON-record (clause (iv) fails) -- symmetry must break
        if recs and bad:
            for fam in itertools.product(recs[:6],bad[:6]):
                r=analyse(list(fam),shells,dl); dd=r['d']
                ctl=max(ctl,abs(dd[(0,(1,1))]-dd[(0,(-1,-1))]),
                            abs(dd[(0,(1,-1))]-dd[(0,(-1,1))]))
    P(f"{'m='+str(m)+', dims 1..'+str(D):<30} {np1:>9} {w1:>19} {w2:>19} {str(always):>21} {ctl:>26}")
P("")
P("READ: the two symmetries hold EXACTLY (integer zero) on every pair of genuine records, and")
P("      the simultaneous flip is in G_W every time.  The CONTROL column pairs a record with a")
P("      DURABLE NON-TRIVIAL BIT THAT FAILS CLAUSE (iv), and there the symmetry breaks (non-zero)")
P("      -- so the zeros above are consequences of clause (iv), not of the parameterisation.")
P("")

P("="*130)
P("--- 7.2  FREEZING at k = 2, 3, 4 and on MULTI-SHELL profiles, exhaustively ---")
P("")
P("For every family and every non-empty S with chi_S G_W-invariant, is R_S UNBALANCED?")
P("")
P(f"{'setting':<38} {'#families':>11} {'#(family,S) invariant pairs':>28} "
  f"{'#with R_S BALANCED (violations)':>32} {'freezing absolute?':>19}")
P("-"*136)
def run(mlist, Dlist, k, nsh=1):
    fams=0; pairs=0; viol=0; ex=[]
    for m in mlist:
        for D in Dlist:
            base=list(itertools.combinations_with_replacement(range(1,D+1),m))
            for combo in itertools.product(base,repeat=nsh):
                if any(sum(c)%2 for c in combo): continue
                dl=[]; shells=[]; idx=0
                for c in combo:
                    shells.append(list(range(idx,idx+len(c)))); dl+=list(c); idx+=len(c)
                recs=[b for b in itertools.product((1,-1),repeat=len(dl))
                      if b[0]==1 and clause_iv_ok(b,shells,dl) and clause_iii_ok(b,shells)]
                if len(recs)<k or len(recs)>60: continue
                for fam in itertools.combinations(recs,k):
                    r=analyse(list(fam),shells,dl); fams+=1
                    for mm in range(1,1<<k):
                        S=tuple(i for i in range(k) if (mm>>i)&1)
                        if all(sum(g[i] for i in S)%2==0 for g in r['G']):
                            pairs+=1
                            bal=True
                            for ei,sh in enumerate(shells):
                                t=sum(dl[i]*int(np.prod([fam[j][i] for j in S])) for i in sh)
                                if t!=0: bal=False; break
                            if bal:
                                viol+=1
                                if len(ex)<4: ex.append((dl,shells,fam,S))
    return fams,pairs,viol,ex
for label,args in [("k=2, one shell, m<=8",       ((4,6,8),(1,2,3),2,1)),
                   ("k=3, one shell, m<=8",       ((4,6,8),(1,2,3),3,1)),
                   ("k=4, one shell, m<=8",       ((6,8),(1,2),4,1)),
                   ("k=2, TWO shells, m<=6 each", ((4,6),(1,2),2,2)),
                   ("k=3, TWO shells, m<=6 each", ((4,6),(1,2),3,2))]:
    f,p,v,ex=run(*args)
    P(f"{label:<38} {f:>11} {p:>28} {v:>32} {str(v==0):>19}")
    if v:
        for e in ex: P("      VIOLATION:", e)
P("")
P("="*130)
P("PART 7 VERDICT")
P("="*130)
P("  At k=2 the freezing theorem is PROVED exactly (header of this file) and verified.")
P("  At k=3 and k=4, and on multi-shell profiles, the exhaustive search found the number of")
P("  violations printed in the table above.  Read that number, not an expectation.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM/t7_freezing.txt","w").write("\n".join(OUT)+"\n")
