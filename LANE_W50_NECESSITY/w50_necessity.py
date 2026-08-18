"""W-50.  IS THE GRAVITY-LIKE STRUCTURE NECESSARY, OR MERELY CO-PRESENT?

Proof obligation C. Showing records CAN form this way is not showing this is HOW they form. The
gravity identification (W-45) rests on capacity satisfying four functional markers. So strip the
structure that produced them and see what survives.

WHAT IS STRIPPABLE. capacity = m - 1 (W-42) is a RANK fact and cannot be removed -- any carrier with
m independent records has it. What CAN be removed is the SHARING: whether record boundaries overlap.

  CONNECTED    3x3 patch: 4 plaquettes, boundaries SHARE links (links 2,3,7,10 lie on two)
  DISCONNECTED four separate 2x2 squares: 4 plaquettes, boundaries pairwise DISJOINT

Same number of records, same cycle rank, same physical dimension. Only the sharing differs.
Run W-45's four markers on both, and W-41's eviction test on both.

  If the markers hold in BOTH -> capacity's gravity-like character does not depend on sharing, and
     is forced by counting alone. Then it is necessary but possibly vacuous.
  If they FAIL on the disconnected carrier -> the gravity-like character comes from the sharing,
     which is a structural property a carrier could lack. Then it is NOT necessary, and records
     could form without anything gravity-like.
Either answer is a result. Neither is assumed.
"""
import itertools, numpy as np

def rank_gf2(vecs):
    piv=[]; r=0
    for v in vecs:
        for p in piv: v=min(v,v^p)
        if v: piv.append(v); piv.sort(reverse=True); r+=1
    return r

def connected():
    nx=ny=3
    V=[(i,j) for j in range(ny) for i in range(nx)]; vid={v:k for k,v in enumerate(V)}
    E=[]
    for j in range(ny):
        for i in range(nx-1): E.append((vid[(i,j)],vid[(i+1,j)]))
    NH=len(E)
    for j in range(ny-1):
        for i in range(nx): E.append((vid[(i,j)],vid[(i,j+1)]))
    hid=lambda i,j:j*(nx-1)+i; vx=lambda i,j:NH+j*nx+i
    PL=[[hid(i,j),vx(i+1,j),hid(i,j+1),vx(i,j)] for j in range(ny-1) for i in range(nx-1)]
    return len(E),PL,len(V),1
def disconnected():
    PL=[]; L=0; NV=0
    for q in range(4):
        PL.append([L,L+1,L+2,L+3]); L+=4; NV+=4
    return L,PL,NV,4

def bd(S,PL):
    c={}
    for p in S:
        for lk in PL[p]: c[lk]=c.get(lk,0)+1
    return set(lk for lk,v in c.items() if v%2)

def analyse(name,L,PL,NV,ncomp):
    m=len(PL)
    REG=[S for r in range(1,m+1) for S in itertools.combinations(range(m),r)]
    cyc=L-NV+ncomp
    mult={lk:sum(1 for p in range(m) if lk in PL[p]) for lk in range(L)}
    def cap_given(held):
        best=0
        for lk in range(L):
            if any(lk in bd(S,PL) for S in held): continue
            free=[S for S in REG if lk not in bd(S,PL)]
            base=[sum(1<<i for i in S) for S in held]
            best=max(best, rank_gf2(base+[sum(1<<i for i in S) for S in free])-rank_gf2(base))
        return best
    base=cap_given([])
    print(f"\n  === {name} ===")
    print(f"    links {L}, vertices {NV}, components {ncomp}, plaquettes {m}, cycle rank {cyc}, "
          f"physical dim {2**cyc}")
    print(f"    link multiplicity (how many plaquettes each link lies on): "
          f"{sorted(set(mult.values()))}  -> {'SHARED boundaries' if max(mult.values())>1 else 'DISJOINT boundaries'}")
    print(f"    capacity = {base}   (m-1 = {m-1})")
    # F3 one sign
    inc=any(cap_given([S])>base-1 for S in REG)
    # F1 universal
    cons=sorted(set(base-cap_given([S]) for S in REG))
    # F2 no screening, gated on capacity>0 (the W-45b fix)
    scr=0; tested=0
    for hold in itertools.combinations(REG,2):
        v=[sum(1<<i for i in S) for S in hold]
        if rank_gf2(v)<2: continue
        before=cap_given(list(hold))
        if before<=0: continue
        for S in REG:
            if S in hold: continue
            if rank_gf2(v+[sum(1<<i for i in S)])<3: continue
            tested+=1
            if cap_given(list(hold)+[S])>=before: scr+=1
    # F4 arena: which k are possible
    poss=[]
    for k in range(1,m+1):
        ok=False
        for combo in itertools.combinations(REG,k):
            if rank_gf2([sum(1<<i for i in S) for S in combo])<k: continue
            if set(range(L))-set().union(*[bd(S,PL) for S in combo]): ok=True; break
        poss.append(ok)
    print(f"    F1 universal sourcing : consumption values {cons}   "
          f"{'PASS (no free record)' if 0 not in cons else 'FAIL'}")
    print(f"    F2 no screening       : {scr} screened of {tested} tested   "
          f"{'PASS' if scr==0 else 'FAIL'}")
    print(f"    F3 one sign           : any record frees capacity? {inc}   "
          f"{'PASS' if not inc else 'FAIL'}")
    print(f"    F4 arena              : k possible = {[k+1 for k,o in enumerate(poss) if o]}, "
          f"impossible = {[k+1 for k,o in enumerate(poss) if not o]}   "
          f"{'PASS (a k is impossible)' if not all(poss) else 'FAIL (all k possible)'}")
    # eviction profile (W-41), first-order counting rate = 2*gamma*|bath INTERSECT boundary|
    singles=[(i,) for i in range(m)]
    prof=None
    for lk in range(L):
        v=sorted(len(bd(S,PL)&{lk}) for S in singles)
        if prof is None or v.count(0)>prof[1].count(0): prof=(lk,v)
    print(f"    W-41 eviction (best single-link bath, overlap counts per record): {prof[1]}"
          f"   -> {prof[1].count(0)} protected, {len(prof[1])-prof[1].count(0)} evicted")
    return base,cons,scr,inc,poss

print("W-50  NECESSITY: does the gravity-like character survive removing the sharing?")
a=analyse("CONNECTED  (3x3 patch, boundaries share links)",*connected())
b=analyse("DISCONNECTED (4 separate squares, boundaries disjoint)",*disconnected())
print()
print("  COMPARISON")
print(f"    capacity           connected {a[0]}   disconnected {b[0]}")
print(f"    F1 consumption     connected {a[1]}  disconnected {b[1]}")
print(f"    F2 screened        connected {a[2]}   disconnected {b[2]}")
print(f"    F3 frees capacity  connected {a[3]}  disconnected {b[3]}")
print(f"    F4 k possible      connected {[k+1 for k,o in enumerate(a[4]) if o]}  "
      f"disconnected {[k+1 for k,o in enumerate(b[4]) if o]}")
