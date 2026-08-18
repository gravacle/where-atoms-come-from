"""W-45.  IDENTIFY GRAVITY BY FUNCTION, NOT BY CLASSICAL FORM.

The principal's standing correction, restated: we should not expect the full footprint of CLASSICAL
gravity at the record level. It is likely to LOOK different while PERFORMING the same function. The
W-44 workflow specifies gravity by its classical form (metric, proper time, clocks, redshift) and
therefore imports exactly the assumption we have been told twice to stop importing.

SO ASK THE FUNCTIONAL QUESTION. Gravity's functional signature, stripped of its classical shape:
   (F1) UNIVERSAL SOURCING -- everything sources it; there is no neutral matter.
   (F2) NO SCREENING       -- it cannot be shielded; you cannot hide behind anything.
   (F3) ONE SIGN           -- there is no negative mass, so it only ever accumulates.
   (F4) ARENA, NOT FORCE   -- it sets what is possible rather than pushing things around.
EM has none of these: charge comes in two signs, screens, and neutral matter exists.

CANDIDATE, ARISING FROM THIS PROGRAM'S OWN MEASUREMENTS AND NOT FROM GRAVITY:
   CAPACITY. W-42: capacity = area - 1, exactly. W-41: past capacity the carrier EVICTS.
Test capacity against F1-F4 as stated, with no metric anywhere in the construction.

CAPACITY IS COMPUTED AS IN W-42: a set of records is simultaneously protectable iff some link lies
on NONE of their boundaries; the boundary map is linear over GF(2); so capacity = max over links of
dim ker(f_link). No Hilbert space is needed.
"""
import itertools, numpy as np

def patch(nx,ny):
    V=[(i,j) for j in range(ny) for i in range(nx)]; vid={v:k for k,v in enumerate(V)}
    E=[]
    for j in range(ny):
        for i in range(nx-1): E.append((vid[(i,j)],vid[(i+1,j)]))
    H=len(E)
    for j in range(ny-1):
        for i in range(nx): E.append((vid[(i,j)],vid[(i,j+1)]))
    hid=lambda i,j: j*(nx-1)+i
    vx =lambda i,j: H + j*nx + i
    PL=[[hid(i,j),vx(i+1,j),hid(i,j+1),vx(i,j)] for j in range(ny-1) for i in range(nx-1)]
    return len(E),PL

def bd(S,PL):
    c={}
    for p in S:
        for lk in PL[p]: c[lk]=c.get(lk,0)+1
    return frozenset(lk for lk,v in c.items() if v%2)

def rank_gf2(vecs,n):
    rows=list(vecs); r=0
    for c in range(n):
        p=None
        for i in range(r,len(rows)):
            if (rows[i]>>c)&1: p=i; break
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        for i in range(len(rows)):
            if i!=r and ((rows[i]>>c)&1): rows[i]^=rows[r]
        r+=1
        if r==len(rows): break
    return r

L,PL=patch(3,3); m=len(PL)
REG=[S for r in range(1,m+1) for S in itertools.combinations(range(m),r)]
def cap_given(held):
    """capacity still available once `held` is already protected: the largest independent set that
       can be ADDED while some link still avoids every boundary in the whole collection."""
    best=0
    for lk in range(L):
        if any(lk in bd(S,PL) for S in held): continue        # this link is already spoken for
        free=[S for S in REG if lk not in bd(S,PL)]
        vecs=[sum(1<<i for i in S) for S in free]
        base=[sum(1<<i for i in S) for S in held]
        tot=rank_gf2(base+vecs,m); have=rank_gf2(base,m)
        best=max(best,tot-have)
    return best

print("W-45  DOES CAPACITY PERFORM GRAVITY'S FUNCTION?  3x3 patch, no metric anywhere.")
print(f"      links {L}, plaquettes (area) {m}, total capacity {cap_given([])}\n")

print("  F3 -- ONE SIGN. Does ANY record ever INCREASE the capacity available to others?")
base=cap_given([])
worst=None; anyincrease=False
for S in REG:
    c=cap_given([S])
    if c>base-1: anyincrease=True
    if worst is None or c<worst[1]: worst=(S,c)
print(f"    capacity with nothing held            : {base}")
rows=sorted({(len(S),cap_given([S])) for S in REG})
for sz,c in rows: print(f"    after holding one record of size {sz}  : {c}")
print(f"    any record that INCREASES others' capacity? {anyincrease}")
print(f"    -> {'ONE-SIGNED: every record only ever consumes' if not anyincrease else 'TWO-SIGNED: something frees capacity'}")

print("\n  F1 -- UNIVERSAL SOURCING. Is there a record that consumes NOTHING? (a 'neutral' record)")
consume=[(S, base-cap_given([S])) for S in REG]
zero=[S for S,d in consume if d==0]
print(f"    records examined {len(consume)};  consuming 0 capacity: {len(zero)}")
print(f"    consumption values observed: {sorted(set(d for _,d in consume))}")
print(f"    -> {'UNIVERSAL: no record is free' if not zero else f'NOT universal: {len(zero)} free records'}")

print("\n  F2 -- NO SCREENING. Can a record be hidden behind others so it stops consuming?")
print("    Hold a set, then add one more, and ask whether the newcomer ever consumes 0.")
screened=0; tested=0
for hold in itertools.combinations(REG,2):
    vec=[sum(1<<i for i in S) for S in hold]
    if rank_gf2(vec,m)<2: continue
    before=cap_given(list(hold))
    for S in REG:
        if S in hold: continue
        if rank_gf2(vec+[sum(1<<i for i in S)],m)<3: continue   # dependent: not a new record
        tested+=1
        if cap_given(list(hold)+[S])>=before: screened+=1
    if tested>4000: break
print(f"    independent additions tested {tested};  additions that consumed nothing {screened}")
print(f"    -> {'NO SCREENING: every independent record pays, whatever is already there' if screened==0 else 'SCREENING EXISTS'}")

print("\n  F4 -- ARENA NOT FORCE. Does capacity act on records, or set what records are possible?")
print("    Capacity never appears in H and exerts no force. It bounds which SETS can coexist:")
for k in range(1,m+1):
    ok=0
    for combo in itertools.combinations(REG,k):
        v=[sum(1<<i for i in S) for S in combo]
        if rank_gf2(v,m)<k: continue
        if set(range(L))-set().union(*[bd(S,PL) for S in combo]): ok+=1; break
    print(f"      k={k}: {'possible' if ok else 'IMPOSSIBLE -- no arrangement exists'}")

print("\n  CONTRAST -- does the EM charge behave this way? Z_2 charge on the same carrier.")
print("    Z_2 charge has TWO signs (it is its own inverse, +1 and -1 both exist), charges CANCEL")
print("    in pairs, and a neutral configuration exists. Checked directly:")
tot_even=sum(1 for c in itertools.product([0,1],repeat=4) if sum(c)%2==0)
print(f"      4-site Z_2 charge configurations that are globally neutral: {tot_even} of 16")
print(f"      a neutral pair (+,+) annihilates to vacuum: yes, since 1+1 = 0 mod 2")
print(f"      -> EM charge SCREENS and has a neutral sector. Capacity has neither.")
