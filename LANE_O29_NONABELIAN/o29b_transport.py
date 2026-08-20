"""O-29b: TRANSPORT BETWEEN RECORDS -- the object the EMERGENCE question names.

The principal: we are looking for where gravity EMERGES FROM the record surface, not for gravity AT
it. Emergence means the property belongs to the COLLECTION and to no single member. Every prior
measurement in this program was on ONE record's neighbourhood, which is why they could only return
nulls.

So: take one record around another and ask what comes back. In a gauge theory a record transported
around a record h returns CONJUGATED, g -> h g h^-1. That is frame transport between records. It
needs at least TWO records to exist at all, and it is the identity for every abelian group -- which
is why every Z_2 measurement in this program returned zero.

ERRATA, corrected below.
 v1  The READ claimed the change GROWS with the number of enclosed records. The v1 data showed 3
     outcomes at k=1,2,3 -- SATURATED, not growing -- and the claim contradicted the table above it.
 v1  Item 4 formed the unordered product of every group element, which is not well-defined for a
     non-abelian group and was not a measurement.
 v2  The READ called closure 'THE STRONGEST RESULT'. Item 5 returns [1] for Z_2 as well as S_3,
     because 'given all but one the last is determined' is the group inverse z = p^-1 and holds in
     ANY group. It does not discriminate and is not an emergence signature."""
import sys, itertools
def say(*a): print(*a); sys.stdout.flush()
def grp(perms):
    d=len(perms[0])
    def mul(a,b): return tuple(a[b[i]] for i in range(d))
    def inv(a):
        r=[0]*d
        for i,x in enumerate(a): r[x]=i
        return tuple(r)
    return list(perms), mul, inv
S3,mul,inv=grp(list(itertools.permutations(range(3)))); e3=(0,1,2)
Z2,mul2,inv2=grp([(0,1),(1,0)]); e2=(0,1)
def classes(G,m,iv):
    seen=set(); out=[]
    for g in G:
        if g in seen: continue
        c=frozenset(m(m(h,g),iv(h)) for h in G); seen|=c; out.append(sorted(c))
    return out
say("="*100); say("O-29b   TRANSPORT BETWEEN RECORDS ON A NON-ABELIAN SURFACE"); say("="*100)
res={}
for name,G,m,iv,ee in (("Z_2  (every carrier built so far)",Z2,mul2,inv2,e2),
                       ("S_3  (the first non-abelian carrier)",S3,mul,inv,e3)):
    cl=classes(G,m,iv); R={}
    say(""); say("-"*100); say(f"  {name}   order {len(G)},  record classes {len(cl)}"); say("-"*100)
    say("  1. ONE RECORD AROUND ANOTHER   g -> h g h^-1")
    moved=sum(1 for g in G for h in G if m(m(h,g),iv(h))!=g)
    say(f"     configurations where the record COMES BACK CHANGED: {moved} of {len(G)**2}")
    say(f"     -> {'TRANSPORT IS NON-TRIVIAL' if moved else 'IDENTITY -- the record is untouched, always'}")
    R['moved']=moved
    say("  2. DENSITY: does each ADDITIONAL enclosed record contribute?")
    say("     fix a record and a sequence of enclosed records; compare enclosing k against k+1")
    say(f"     {'k -> k+1':>12}{'sequences where the outcome CHANGES':>40}{'of':>8}")
    for k in (1,2,3):
        ch=0; tt=0
        for g in G:
            for seq in itertools.product(G,repeat=k+1):
                p=ee
                for x in seq[:k]: p=m(p,x)
                q=m(p,seq[k]); tt+=1
                if m(m(p,g),iv(p))!=m(m(q,g),iv(q)): ch+=1
        say(f"     {f'{k} -> {k+1}':>12}{ch:>40}{tt:>8}")
        R[f'dens{k}']=(ch,tt)
    say("  3. ORDER: enclose two records, h1 then h2, against h2 then h1")
    dis=sum(1 for g in G for h1 in G for h2 in G
            if m(m(m(h1,h2),g),iv(m(h1,h2)))!=m(m(m(h2,h1),g),iv(m(h2,h1))))
    say(f"     configurations where the ORDER MATTERS: {dis} of {len(G)**3}")
    say(f"     -> {'PATH-DEPENDENT MEMORY' if dis else 'order is irrelevant -- no path memory'}")
    R['order']=dis
    say("  4. CLOSURE: on a closed surface the ORDERED product of all record holonomies must be e.")
    say("     This is a constraint linking the records to EACH OTHER -- H-9's 'curvature cannot")
    say("     appear arbitrarily'. How many placements of k records are admissible?")
    say(f"     {'k records':>12}{'admissible (prod = e)':>24}{'of':>10}{'fraction':>12}")
    for k in (2,3,4):
        ok=0; tt=0
        for seq in itertools.product(G,repeat=k):
            p=ee
            for x in seq: p=m(p,x)
            tt+=1; ok+= (p==ee)
        say(f"     {k:>12}{ok:>24}{tt:>10}{ok/tt:>12.6f}")
        R[f'clos{k}']=(ok,tt)
    say("  5. CONTROL -- given the first k-1 records, how many choices remain for the last?")
    say("     This is the group inverse z = p^-1 and must return 1 in ANY group. It is here as a")
    say("     CONTROL: if a signature returns the same answer for Z_2 it discriminates nothing.")
    det=set()
    for k in (2,3,4):
        for seq in itertools.product(G,repeat=k-1):
            p=ee
            for x in seq: p=m(p,x)
            det.add(len([z for z in G if m(p,z)==ee]))
    say(f"     number of admissible completions, over all prefixes: {sorted(det)}")
    say(f"     -> {'determined -- as it must be in any group; NOT a discriminator' if det=={1} else 'not determined'}")
    R['det']=sorted(det); res[name]=R
say(""); say("="*100); say("  READ"); say("="*100)
z,s3=res["Z_2  (every carrier built so far)"],res["S_3  (the first non-abelian carrier)"]
say("  Transport between records requires TWO records; no single record has it. It is the IDENTITY")
say("  for every abelian group, so on a Z_2 carrier a record taken around another comes back")
say(f"  UNCHANGED ({z['moved']} of 4), the order never matters ({z['order']} of 8), and no additional enclosed record")
say("  changes anything (0 at every k).")
say("")
say("  THAT SINGLE FACT EXPLAINS THE PROGRAM'S ENTIRE NEGATIVE RECORD ON THIS QUESTION. Every")
say("  carrier built here has gauge group Z_2. The nulls were not evidence against emergence; they")
say("  were measured where the effect is identically zero by construction.")
say("")
say(f"  On S_3: the record COMES BACK CHANGED ({s3['moved']} of 36); the ORDER of encirclement MATTERS ({s3['order']} of 216);")
say(f"  and EACH ADDITIONAL enclosed record changes the outcome in half of all sequences")
say(f"  ({s3['dens1'][0]}/{s3['dens1'][1]}, {s3['dens2'][0]}/{s3['dens2'][1]}, {s3['dens3'][0]}/{s3['dens3'][1]}). Path-dependent transport between records is PRESENT.")
say("")
say("  WHICH SIGNATURES ACTUALLY DISCRIMINATE")
say(f"  {'signature':<34}{'Z_2':>12}{'S_3':>14}{'discriminates?':>18}")
for lbl,a,b in (("record comes back changed",f"{z['moved']}/4",f"{s3['moved']}/36"),
                ("order of encirclement matters",f"{z['order']}/8",f"{s3['order']}/216"),
                ("each added record contributes",f"{z['dens1'][0]}/{z['dens1'][1]}",f"{s3['dens1'][0]}/{s3['dens1'][1]}"),
                ("closure fraction",f"{z['clos2'][0]}/{z['clos2'][1]}",f"{s3['clos2'][0]}/{s3['clos2'][1]}"),
                ("last record determined",str(z['det']),str(s3['det']))):
    d="YES" if lbl in ("record comes back changed","order of encirclement matters","each added record contributes") else "no"
    say(f"  {lbl:<34}{a:>12}{b:>14}{d:>18}")
say("")
say("  Closure holds at fraction 1/|G| on BOTH carriers and the last record is determined on BOTH,")
say("  so neither is an emergence signature -- both are properties of being a group at all.")
say("")
say("  WHAT IS NOT HERE: a MONOTONE density law. The reachable outcomes are bounded by the")
say("  conjugacy class, |class| = 3 on S_3, so 'more records -> more transport' cannot be seen on a")
say("  group this small. That is H-4's small-carrier warning, and it is the next thing to measure.")
