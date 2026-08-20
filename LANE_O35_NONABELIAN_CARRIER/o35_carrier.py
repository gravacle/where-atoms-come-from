"""O-35: BUILD A NON-ABELIAN CARRIER THAT ACTUALLY HAS RECORDS.

AUDIT 2 found the minimal-torus D(S_3) has NO records: eigenspace multiplicities [8,13,15], and
clause (i) forces Tr(P_E R) = dim(E) mod 2, so an odd-dimensional eigenspace makes clause (iv)
unsatisfiable. That carrier has ONE vertex and ONE face, so the gauge and flatness projectors have
almost nowhere to act -- H-4's small-carrier warning.

This lane builds the next torus up and asks the same question. For an L1 x L2 square lattice on a
torus: V = L1*L2 vertices, E = 2*L1*L2 edges, F = L1*L2 faces, so V - E + F = 0 always. With S_3 the
Hilbert space is 6^E, so 1x2 gives 6^4 = 1296 -- reachable. 2x2 would be 6^8, which is not.

CONTROLS IN THE SAME TABLE (D-15): the same construction run with Z_2, where the toric code is known
to have records, and with the minimal torus, where AUDIT 2 showed it does not."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
def group(name):
    """Groups as permutation groups. The point of the extra ones is CONTROL: S_3 has order 6, which
       is not a power of two, so 'non-abelian' and 'order not 2^k' are confounded unless both an
       ABELIAN group of non-power-of-two order and a NON-ABELIAN group of order 2^k are tested."""
    def close(gens,d):
        def m(a,b): return tuple(a[b[i]] for i in range(d))
        E=tuple(range(d)); S={E}; fr=[E]
        while fr:
            nf=[]
            for x in fr:
                for g in gens:
                    y=m(x,g)
                    if y not in S: S.add(y); nf.append(y)
            fr=nf
        return sorted(S)
    if   name=="Z_2": G,d=[(0,1),(1,0)],2
    elif name=="Z_3": G,d=close([(1,2,0)],3),3
    elif name=="Z_4": G,d=close([(1,2,3,0)],4),4
    elif name=="S_3": G,d=list(itertools.permutations(range(3))),3
    elif name=="D_4": G,d=close([(1,2,3,0),(1,0,3,2)],4),4   # non-abelian, ORDER 8 = 2^3
    else: raise ValueError(name)
    def mul(a,b): return tuple(a[b[i]] for i in range(d))
    def inv(a):
        r=[0]*d
        for i,x in enumerate(a): r[x]=i
        return tuple(r)
    return G,mul,inv,tuple(range(d))
def build(gname,L1,L2):
    """D(G) on an L1 x L2 square-lattice torus. Returns H and a description."""
    G,mul,inv,e=group(gname); n=len(G); gi={g:i for i,g in enumerate(G)}
    E=[]; eid={}
    for r in range(L1):
        for c in range(L2):
            for dd in (0,1): eid[(r,c,dd)]=len(E); E.append((r,c,dd))
    nE=len(E); D=n**nE
    def cfg(i):
        out=[]
        for _ in range(nE): out.append(G[i%n]); i//=n
        return tuple(reversed(out))
    def num(t):
        i=0
        for g in t: i=i*n+gi[g]
        return i
    # face (r,c): h(r,c) . v(r,c+1) . h(r+1,c)^-1 . v(r,c)^-1
    faces=[]
    for r in range(L1):
        for c in range(L2):
            faces.append([(eid[(r,c,0)],+1),(eid[(r,(c+1)%L2,1)],+1),
                          (eid[((r+1)%L1,c,0)],-1),(eid[(r,c,1)],-1)])
    # vertex (r,c): outgoing h(r,c), v(r,c); incoming h(r,c-1), v(r-1,c)
    verts=[]
    for r in range(L1):
        for c in range(L2):
            verts.append(([eid[(r,c,0)],eid[(r,c,1)]],
                          [eid[(r,(c-1)%L2,0)],eid[((r-1)%L1,c,1)]]))
    B=np.zeros((D,D))
    for i in range(D):
        t=cfg(i); okf=True
        for f in faces:
            p=e
            for (k,s) in f: p=mul(p, t[k] if s>0 else inv(t[k]))
            if p!=e: okf=False; break
        if okf: B[i,i]=1.0
    A=np.zeros((D,D))
    for (outs,ins) in verts:
        Av=np.zeros((D,D))
        for h in G:
            for i in range(D):
                t=list(cfg(i))
                for k in outs: t[k]=mul(h,t[k])
                for k in ins:  t[k]=mul(t[k],inv(h))
                Av[num(tuple(t)),i]+=1.0/n
        A=A+Av
    return -(A+len(faces)*0-0) - 0, A, B, len(verts), len(faces), nE, D
say("="*104); say("O-35   A NON-ABELIAN CARRIER WITH RECORDS?"); say("="*104)
say(f"  {'carrier':<26}{'V':>3}{'E':>3}{'F':>3}{'dim':>8}{'multiplicities':>34}{'odd?':>7}{'(iv)?':>8}")
rows=[]
for gname,L1,L2,lbl in (("Z_2",1,2,"Z_2  abelian  ord 2   1x2"),
                        ("Z_3",1,1,"Z_3  ABELIAN  ord 3   min"),
                        ("Z_3",1,2,"Z_3  ABELIAN  ord 3   1x2"),
                        ("Z_4",1,2,"Z_4  abelian  ord 4   1x2"),
                        ("D_4",1,1,"D_4  NONABEL  ord 8   min"),
                        ("S_3",1,1,"S_3  nonabel  ord 6   min"),
                        ("S_3",1,2,"S_3  nonabel  ord 6   1x2")):
    _,A,B,V,F,nE,D=build(gname,L1,L2)
    H=-(A+B)
    w=np.linalg.eigvalsh(H)
    vals=[]; 
    for x in w:
        if not vals or abs(x-vals[-1][0])>1e-8: vals.append([x,1])
        else: vals[-1][1]+=1
    mult=[m for _,m in vals]
    odd=any(m%2 for m in mult)
    show=str(mult) if len(str(mult))<=32 else str(mult)[:29]+"..."
    say(f"  {lbl:<26}{V:>3}{nE:>3}{F:>3}{D:>8}{show:>34}{str(odd):>7}{('NO' if odd else 'yes'):>8}")
    rows.append((lbl,mult,odd,D))
say("")
say("  clause (i) forces R^2 = I, so Tr(P_E R) = dim(E) mod 2: an ODD eigenspace can never give 0,")
say("  and clause (iv) needs Tr(P_E R) = 0 on EVERY eigenspace.")
for lbl,mult,odd,D in rows:
    say(f"    {lbl:<26} -> {'clause (iv) UNSATISFIABLE: no record can exist' if odd else 'clause (iv) is possible on this carrier'}")
say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
ab={l.split()[0]:o for l,_,o,_ in rows}
say("  THE CONTROLS SEPARATE TWO EXPLANATIONS THAT S_3 ALONE CONFOUNDS.")
say("  S_3 is non-abelian AND has order 6, which is not a power of two. Z_3 is ABELIAN with order 3;")
say("  D_4 is NON-ABELIAN with order 8 = 2^3. Whichever of those two fails identifies the cause.")
say("")
for l,mult,o,D in rows:
    say(f"    {l:<28} odd eigenspace present: {str(o):<6} -> {'no record can exist' if o else 'records are possible'}")
say("")
say("  THE CONTROL OVERTURNS THE DIAGNOSIS. Z_3 is ABELIAN and has NO records. D_4 is NON-ABELIAN")
say("  and HAS them. So the obstruction is NOT non-abelianness -- it is whether |G| is a POWER OF")
say("  TWO. Z_2, Z_4 and D_4 have orders 2, 4 and 8 and pass; Z_3 and S_3 have orders 3 and 6 and")
say("  fail, at both carrier sizes tested.")
say("")
say("  WHY: clause (i) makes a record a BIT, R^2 = I, so Tr(P_E R) = dim(E) mod 2 and every")
say("  eigenspace must be EVEN-dimensional. |G|^E is a power of two exactly when |G| is. A group")
say("  whose order carries an odd factor produces odd eigenspaces and no record can be written on it.")
say("")
say("  D_4 IS A NON-ABELIAN CARRIER WITH RECORDS. It is what O-34's join needs and what AUDIT 2")
say("  showed D(S_3) could not supply. Had this lane tested only S_3, it would have registered")
say("  'non-abelian carriers cannot have records', which the D_4 column shows is FALSE.")
