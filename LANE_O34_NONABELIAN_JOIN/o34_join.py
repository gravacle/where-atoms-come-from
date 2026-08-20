"""O-34: THE NON-ABELIAN JOIN, ON A CARRIER THAT ACTUALLY HAS RECORDS.

O-35 supplied D_4: non-abelian (24 of 64 pairs do not commute, so flux transport is non-trivial) AND
with all eigenspace multiplicities even, so clause (iv) is satisfiable. D(S_3) had the first property
and not the second; every Z_2 carrier had the second and not the first.

model.records() cannot enumerate here -- 64 minimal projections against O-28's cap of 20 -- so a
record is CONSTRUCTED and every clause CHECKED, rather than searched for or assumed.

The natural record is a Wilson loop in a one-dimensional representation: R|g1,g2> = chi(g1)|g1,g2>
for a linear character chi: G -> {+1,-1}. It is a bit by construction. Whether it is DURABLE,
NON-TRIVIAL and WRITABLE is measured.

THE QUESTION THIS LANE EXISTS TO ASK. Clause (ii) requires [H,R] = 0, and H contains the GAUGE
projector. Gauge transport is conjugation. So the clause that makes a record durable may be the same
clause that makes it untransportable -- in which case no record can ever be moved by transport, on
any carrier, and that would explain every null this program has recorded. Measured below, with the
toric code as the control."""
import sys, os, itertools, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel
def say(*a): print(*a); sys.stdout.flush()
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
def carrier(gens,d,label):
    G=close(gens,d); n=len(G); gi={g:i for i,g in enumerate(G)}
    def mul(a,b): return tuple(a[b[i]] for i in range(d))
    def inv(a):
        r=[0]*d
        for i,x in enumerate(a): r[x]=i
        return tuple(r)
    e=tuple(range(d)); D=n*n
    def ket(a,b): return gi[a]*n+gi[b]
    Ah={}
    for h in G:
        Mh=np.zeros((D,D))
        for g1 in G:
            for g2 in G:
                Mh[ket(mul(mul(h,g1),inv(h)),mul(mul(h,g2),inv(h))), ket(g1,g2)]=1.0
        Ah[h]=Mh
    A=sum(Ah.values())/n
    B=np.zeros((D,D))
    for g1 in G:
        for g2 in G:
            if mul(mul(g1,g2),mul(inv(g1),inv(g2)))==e: B[ket(g1,g2),ket(g1,g2)]=1.0
    H=-(A+B)
    # linear characters: all homomorphisms G -> {+1,-1}
    chars=[]
    for bits in itertools.product((1,-1),repeat=n):
        c=dict(zip(G,bits))
        if all(c[mul(a,b)]==c[a]*c[b] for a in G for b in G): chars.append(c)
    return dict(label=label,G=G,n=n,D=D,mul=mul,inv=inv,e=e,ket=ket,Ah=Ah,A=A,B=B,H=H,chars=chars,gi=gi)
def record(C,chi,which):
    D=C['D']; R=np.zeros((D,D))
    for g1 in C['G']:
        for g2 in C['G']:
            R[C['ket'](g1,g2),C['ket'](g1,g2)]=chi[g1 if which==0 else g2]
    return R
say("="*104); say("O-34   THE NON-ABELIAN JOIN, ON A CARRIER WITH RECORDS"); say("="*104)
CS=[carrier([(1,0)],2,"D(Z_2)  abelian  (control)"),
    carrier([(1,2,3,0),(1,0,3,2)],4,"D(D_4)  NON-ABELIAN")]
for C in CS:
    M=RecordModel(C['H'].astype(complex))
    nc=sum(1 for a in C['G'] for b in C['G'] if C['mul'](a,b)!=C['mul'](b,a))
    say(""); say("-"*104)
    say(f"  {C['label']}   |G|={C['n']}  dim {C['D']}  non-commuting pairs {nc}/{C['n']**2}"
        f"  multiplicities {[m for _,_,m in M.es]}")
    say("-"*104)
    nontrivial=[c for c in C['chars'] if any(v==-1 for v in c.values())]
    say(f"  linear characters: {len(C['chars'])} total, {len(nontrivial)} non-trivial")
    if not nontrivial:
        say("  no non-trivial character -- no Wilson-loop record of this form exists here"); continue
    chi=nontrivial[0]
    R=record(C,chi,0)
    say("  CLAUSES, checked on R|g1,g2> = chi(g1)|g1,g2>:")
    say(f"    (i)   R=R^dag, R^2=I         ||R-R^T||={np.linalg.norm(R-R.T):.2e}  ||R^2-I||={np.linalg.norm(R@R-np.eye(C['D'])):.2e}")
    say(f"    (ii)  [H,R]=0                ||[H,R]|| = {np.linalg.norm(C['H']@R-R@C['H']):.2e}")
    iii=any(abs(np.trace(P@R))>1e-9 and abs(abs(np.trace(P@R))-np.trace(P))>1e-9 for _,P,_ in
            [(v,Pj,m) for v,Pj,m in M.es]) or any(
            np.linalg.norm(P@R@P - (np.trace(P@R)/np.trace(P))*P)>1e-9 for _,P,_ in M.es)
    say(f"    (iii) non-constant on an eigenspace: {iii}")
    tr=[float(np.real(np.trace(P@R))) for _,P,_ in M.es]
    say(f"    (iv)  Tr(P_E R) per eigenspace: {[round(t,6) for t in tr]}   all zero: {all(abs(t)<1e-9 for t in tr)}")
    isrec = np.linalg.norm(C['H']@R-R@C['H'])<1e-9 and iii and all(abs(t)<1e-9 for t in tr)
    say(f"    -> {'THIS IS A RECORD by clauses (i)-(iv)' if isrec else 'NOT a record on this carrier'}")
    if not isrec: continue
    say("")
    say("  DOES TRANSPORT MOVE IT?  gauge transport is conjugation by h; the record moves iff [A_h,R] != 0")
    worst=0.0; moved=0
    for h in C['G']:
        c=np.linalg.norm(C['Ah'][h]@R - R@C['Ah'][h])
        worst=max(worst,c); moved += (c>1e-9)
    say(f"    group elements whose transport MOVES the record: {moved} of {C['n']}")
    say(f"    largest ||[A_h, R]|| over the whole group: {worst:.3e}")
    say(f"    -> {'TRANSPORT MOVES THE RECORD' if moved else 'TRANSPORT LEAVES THE RECORD FIXED, for every element of the group'}")
say(""); say("="*104)
say("  THE WILSON-LOOP CHARACTER RECORD WORKS ON Z_2 AND FAILS ON D_4. Even multiplicities make")
say("  clause (iv) POSSIBLE, not automatic. So construct records DIRECTLY on the eigenspaces, where")
say("  (ii), (iii) and (iv) hold by construction, and ask the real question of THEM.")
say("="*104)
rng=np.random.default_rng(7)
for C in CS:
    M=RecordModel(C['H'].astype(complex)); D=C['D']
    say(""); say(f"  {C['label']}   dim {D}   multiplicities {[m for _,_,m in M.es]}")
    if any(m%2 for _,_,m in M.es):
        say("    an eigenspace is odd-dimensional -- no record exists here"); continue

    # Build the eigenspaces from H DIRECTLY rather than relying on the layout of M.es.
    wH,VH=np.linalg.eigh(C['H'])
    blocks=[]
    for i,x in enumerate(wH):
        if blocks and abs(x-blocks[-1][0])<1e-8: blocks[-1][1].append(i)
        else: blocks.append([x,[i]])
    say(f"    eigenspaces from H: {[(round(float(v),4),len(ix)) for v,ix in blocks]}")
    made=0; movers=0; worst=0.0; fixed_all=True
    for trial in range(40):
        R=np.zeros((D,D),dtype=complex)
        for _,ix in blocks:
            U=VH[:,ix]; m=len(ix)
            Q=np.linalg.qr(rng.normal(size=(m,m)))[0]        # a random orthonormal basis OF the eigenspace
            sgn=np.array([1.0]*(m//2)+[-1.0]*(m-m//2))
            W=U@Q
            R += W @ np.diag(sgn) @ W.conj().T
        okI=np.linalg.norm(R@R-np.eye(D))<1e-8 and np.linalg.norm(R-R.conj().T)<1e-8
        okII=np.linalg.norm(C['H']@R-R@C['H'])<1e-8
        okIV=all(abs(float(np.real(np.trace(VH[:,ix].conj().T@R@VH[:,ix]))))<1e-8 for _,ix in blocks)
        okIII=any(np.linalg.norm(VH[:,ix].conj().T@R@VH[:,ix]
                  - (np.trace(VH[:,ix].conj().T@R@VH[:,ix])/len(ix))*np.eye(len(ix)))>1e-8 for _,ix in blocks)
        if not (okI and okII and okIII and okIV): continue
        made+=1
        mv=sum(1 for h in C['G'] if np.linalg.norm(C['Ah'][h]@R-R@C['Ah'][h])>1e-8)
        wr=max(np.linalg.norm(C['Ah'][h]@R-R@C['Ah'][h]) for h in C['G'])
        worst=max(worst,wr); movers += (mv>0); fixed_all &= (mv==0)
    say(f"    records constructed and verified against (i)-(iv): {made} of 40")
    say(f"    of those, records that SOME group element's transport MOVES: {movers}")
    say(f"    largest ||[A_h, R]|| seen: {worst:.3e}")
    say(f"    -> {'TRANSPORT MOVES RECORDS ON THIS CARRIER' if movers else 'TRANSPORT LEAVES EVERY RECORD FIXED'}")
say(""); say("="*104)
say("  IS '40 OF 40 MOVED' STRUCTURE, OR JUST THAT RANDOM OPERATORS COMMUTE WITH NOTHING?")
say("  The decisive question is whether a transport-FIXED record EXISTS AT ALL.")
say("")
say("  ERRATUM: a first attempt decomposed each eigenspace by the eigenvalues of sum_h (A_h + A_h^t)")
say("  and did a subset-sum on those block dimensions. It reported NO transport-fixed record on")
say("  D(Z_2) -- CONTRADICTING the direct measurement four lines above it, which found 0 of 40 moved")
say("  with ||[A_h,R]|| = 0.000e+00 exactly. The bug: for an abelian group A_h IS the identity, so")
say("  the commutant is EVERYTHING, but a symmetrised sum of identities has one eigenvalue and the")
say("  method reported a single coarse block. Replaced with the model's own commutant.")
say("")
say("  A record fixed by transport is exactly a record of the model built with the A_h as the")
say("  theory's operators: clause (ii) then demands [A_h,R] = 0 as well as [H,R] = 0.")
say("="*104)
for C in CS:
    Mall=RecordModel(C['H'].astype(complex))
    Mfix=RecordModel(C['H'].astype(complex), Ls=[C['Ah'][h].astype(complex) for h in C['G']])
    say(""); say(f"  {C['label']}   dim {C['D']}")
    say(f"    minimal projections, records only durable against H      : {len(Mall.projs)}")
    say(f"    minimal projections, ALSO durable against every A_h      : {len(Mfix.projs)}")
    def count(MM,lbl):
        try:
            r=MM.records(); return len(r), f"{len(r)}"
        except RuntimeError as ex:
            return None, f"REFUSED ({ex})"
    na,sa=count(Mall,"all"); nf,sf=count(Mfix,"fixed")
    say(f"    records satisfying (i)-(iv) against H alone              : {sa}")
    say(f"    records ALSO commuting with every A_h -- TRANSPORT-FIXED : {sf}")
    if nf is not None:
        say(f"    -> {'TRANSPORT-FIXED RECORDS EXIST: transport does not have to move a record' if nf>0 else 'NO TRANSPORT-FIXED RECORD EXISTS: every record on this carrier is moved by transport'}")
say(""); say("="*104); say("  READ -- from the numbers above"); say("="*104)