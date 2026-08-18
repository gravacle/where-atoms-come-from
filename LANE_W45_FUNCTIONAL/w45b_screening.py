"""W-45b.  F2 (NO SCREENING) FAILED. IS THAT PHYSICS OR A BROKEN TEST?

w45 counted 108 of 1260 independent additions as 'consuming nothing'. Suspicion: when the held set
has ALREADY EXHAUSTED capacity (0 left), adding anything more trivially consumes nothing, because
there was nothing left to consume. That is SATURATION, not screening, and the test did not exclude it.

FIX: only count an addition as screened if capacity was STRICTLY POSITIVE before it was added.
Then screening means: capacity remained available, a new independent record was added, and it paid
nothing. That is the real question.
"""
import itertools, numpy as np
exec(open('w45_function.py').read().split('print("W-45 ')[0])

L,PL=patch(3,3); m=len(PL)
REG=[S for r in range(1,m+1) for S in itertools.combinations(range(m),r)]
def cap_given(held):
    best=0
    for lk in range(L):
        if any(lk in bd(S,PL) for S in held): continue
        free=[S for S in REG if lk not in bd(S,PL)]
        base=[sum(1<<i for i in S) for S in held]
        tot=rank_gf2(base+[sum(1<<i for i in S) for S in free],m)
        best=max(best,tot-rank_gf2(base,m))
    return best

print("W-45b  F2 RE-TESTED, EXCLUDING SATURATION")
print(f"  {'held':>6s} {'cases':>7s} {'cap>0 before':>13s} {'paid nothing':>13s} {'screened (cap>0)':>17s}")
print("  "+"-"*62)
tot_screen=0
for hk in (1,2):
    cases=0; pos=0; free=0; scr=0
    for hold in itertools.combinations(REG,hk):
        v=[sum(1<<i for i in S) for S in hold]
        if rank_gf2(v,m)<hk: continue
        before=cap_given(list(hold))
        for S in REG:
            if S in hold: continue
            if rank_gf2(v+[sum(1<<i for i in S)],m)<hk+1: continue
            cases+=1
            after=cap_given(list(hold)+[S])
            if after>=before: free+=1
            if before>0:
                pos+=1
                if after>=before: scr+=1
    tot_screen+=scr
    print(f"  {hk:6d} {cases:7d} {pos:13d} {free:13d} {scr:17d}")
print()
print(f"  screened cases with capacity genuinely available: {tot_screen}")
print(f"  -> {'NO SCREENING CONFIRMED: every independent record pays whenever there is anything to pay'  if tot_screen==0 else 'SCREENING IS REAL'}")
print()
print("  And the saturation reading, stated explicitly so it is not mistaken for screening:")
h=(REG[0],REG[1])
print(f"    holding {h[0]} and {h[1]}: capacity left = {cap_given(list(h))}")
h3=(REG[0],REG[1],REG[2])
print(f"    holding three                      : capacity left = {cap_given(list(h3))}")
print(f"    once capacity is 0 nothing further can be charged -- the carrier is full, and W-41")
print(f"    already measured what happens then: it EVICTS.")
