"""O-30: DOES TRANSPORT GROW WITH HOW MANY RECORDS ARE ENCLOSED?

H-9's density signature: 'stronger with enclosed record-density'. O-29 measured exactly 50% at every
k on S_3, which is not a coincidence and not a small-group artifact:

  THE PRODUCT OF INDEPENDENT UNIFORM GROUP ELEMENTS IS UNIFORM, IN ANY GROUP.

So enclosing k uniformly-drawn records can NEVER show a density law -- the transport is already fully
mixed after one. A density law can only exist if the enclosed records are of a PARTICULAR KIND. That
is the physically right condition anyway: mass density is a lot of the SAME stuff, not a lot of
uniformly random stuff.

So: enclose k records ALL FROM ONE CONJUGACY CLASS and measure how transport depends on k.
Abelian controls are carried in the same table (D-15)."""
import sys, itertools
def say(*a): print(*a); sys.stdout.flush()
# ---------- groups as permutation groups, generated from generators ----------
def close(gens,d):
    def mul(a,b): return tuple(a[b[i]] for i in range(d))
    e=tuple(range(d)); G={e}; frontier=[e]
    while frontier:
        nf=[]
        for x in frontier:
            for g in gens:
                y=mul(x,g)
                if y not in G: G.add(y); nf.append(y)
        frontier=nf
    return sorted(G)
def ops(d):
    def mul(a,b): return tuple(a[b[i]] for i in range(d))
    def inv(a):
        r=[0]*d
        for i,x in enumerate(a): r[x]=i
        return tuple(r)
    return mul,inv
def cyc(n):  # Z_n
    return close([tuple((i+1)%n for i in range(n))],n), n
GROUPS=[]
for n in (2,4,6):
    G,d=cyc(n); GROUPS.append((f"Z_{n}",G,d,True))
GROUPS.append(("S_3", close([(1,0,2),(1,2,0)],3),3,False))
GROUPS.append(("D_4", close([(1,2,3,0),(1,0,3,2)],4),4,False))
GROUPS.append(("A_4", close([(1,2,0,3),(1,0,3,2)],4),4,False))
GROUPS.append(("S_4", close([(1,0,2,3),(1,2,3,0)],4),4,False))
GROUPS.append(("S_5", close([(1,0,2,3,4),(1,2,3,4,0)],5),5,False))
say("="*104); say("O-30   DOES TRANSPORT GROW WITH ENCLOSED RECORD DENSITY?"); say("="*104)
say("  enclose k records ALL OF ONE KIND (one conjugacy class); T(k) = fraction of records that come")
say("  back CHANGED.  Uniformly-drawn records cannot show a density law: the product of independent")
say("  uniform group elements is uniform in ANY group, which is why O-29 read 50% at every k.")
say("")
KMAX=12
for name,G,d,abelian in GROUPS:
    mul,inv=ops(d)
    # conjugacy classes
    seen=set(); cls=[]
    for g in G:
        if g in seen: continue
        c=sorted({mul(mul(h,g),inv(h)) for h in G}); seen|=set(c); cls.append(c)
    # pick the largest non-central class (the 'kind of record' with the most structure)
    noncentral=[c for c in cls if len(c)>1]
    pick=max(noncentral,key=len) if noncentral else cls[0]
    say("-"*104)
    say(f"  {name:<6} order {len(G):<4} classes {len(cls):<3} {'ABELIAN (control)' if abelian else 'non-abelian'}"
        f"   enclosed record kind: class of size {len(pick)}")
    # exact T(k) by distribution over products of k elements of `pick`
    dist={tuple(range(d)):1.0}
    rows=[]
    for k in range(1,KMAX+1):
        nd={}
        for p,w in dist.items():
            for x in pick: nd[mul(p,x)]=nd.get(mul(p,x),0.0)+w/len(pick)
        dist=nd
        # T(k) = P over (holonomy ~ dist, record ~ uniform G) that h g h^-1 != g
        T=sum(w*sum(1 for g in G if mul(mul(p,g),inv(p))!=g)/len(G) for p,w in dist.items())
        rows.append(T)
    say(f"     {'k enclosed':<12}" + "".join(f"{k:>8}" for k in range(1,KMAX+1)))
    say(f"     {'T(k)':<12}" + "".join(f"{t:>8.4f}" for t in rows))
    mono=all(rows[i+1]>=rows[i]-1e-12 for i in range(len(rows)-1))
    flat=max(rows)-min(rows)<1e-12
    osc=any(rows[i+1]<rows[i]-1e-12 for i in range(len(rows)-1)) and not flat
    say(f"     -> {'IDENTICALLY ZERO -- nothing to rotate' if max(rows)<1e-12 else ('FLAT in k' if flat else ('MONOTONE INCREASING in k' if mono else 'CONVERGES with damped oscillation -- it does NOT grow with k'))}")
    # WHY: the enclosed records of one kind perform a random walk on the subgroup they GENERATE,
    # and that walk MIXES. The limit is therefore an exact, predictable number -- checked, not asserted.
    H=close(pick,d)
    comm=sum(1 for a in H for g in G if mul(mul(a,g),inv(a))==g)
    lim=1.0-comm/(len(H)*len(G))
    say(f"     one kind of record generates a subgroup of order {len(H)} inside order {len(G)}"
        f"  ->  predicted mixing limit {lim:.4f}")
    # The walk is PERIODIC, not ergodic -- S_3 and D_4 have period 2, A_4 has period 3 -- so the
    # limit is the TIME AVERAGE. Average the last 6 values: a whole number of periods for 1,2,3,6.
    tail=rows[-6:]; avg=sum(tail)/len(tail)
    say(f"     time-average of T over k = {KMAX-5}..{KMAX}: {avg:.4f}"
        f"   {'MATCHES the mixing limit' if abs(avg-lim)<5e-3 else 'DOES NOT MATCH'}")
say("-"*104); say("")
say("  READ")
say("  Every abelian control is identically zero at every k: no amount of enclosed record density")
say("  produces any transport at all. The density question does not even arise there.")
say("  On the non-abelian carriers transport is non-zero, but T(k) does NOT grow monotonically with")
say("  the number of enclosed records -- it saturates, or oscillates with the parity of k. So the")
say("  signature H-9 names, 'stronger with enclosed record-density', is NOT PRESENT as a monotone")
say("  law in the counting of same-kind records on any finite group tested.")
say("  A count of records is the wrong density variable. What varies monotonically in a group is not")
say("  how MANY elements are multiplied but WHICH class the product lands in.")
