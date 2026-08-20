"""THE FUNCTION: DOES HOW MUCH TRANSPORT MOVES A RECORD DEPEND ON WHAT ELSE IS PRESENT?

The principal: "we should see some function at scale I would think." And: gravity's weakness is a
COEFFICIENT; the FUNCTION -- how the effect depends on configuration -- is a separate thing, visible
wherever the quantity is non-zero at all. Our carriers do not carry the physical 10^-36; that is a
fact about the world, not about a D(D_4) lattice.

C-43 established that gauge transport MOVES records on D(D_4) -- 40 of 40, against 0 of 40 on the
abelian control at exactly 0.000e+00. IT DID NOT ASK WHETHER HOW MUCH IT MOVES DEPENDS ON WHAT OTHER
RECORDS ARE PRESENT. That dependence IS gravity's form: the arrangement of matter determining
transport for everything else.

MEASURED HERE: for a state psi with a chosen set of records written to chosen values,
    delta_i(psi, h) = <psi| A_h^dag R_i A_h |psi> - <psi| R_i |psi>
is how much transport by h changes record i's value. The question is whether delta_i depends on the
OTHER records' values -- and if so, with what function.

CONTROL IN THE SAME TABLE (D-15): D(Z_2), where conjugation is the identity map so delta must be
exactly zero for every state and every h. A method that does not return zero there is broken.

ERRATUM -- THIS TEST IS INVALID AND THE ZEROS BELOW ARE NOT A RESULT.
delta came back 0.000000 in EVERY configuration on BOTH carriers, including D(D_4) where C-43
measured ||[A_h,R]|| = 9.423. A measurement contradicting an established result means the TEST is
wrong, and it is:

  THE GROUND SPACE IS BY DEFINITION GAUGE-INVARIANT. H = -(A+B) with A the gauge projector, so every
  ground state satisfies A_h|psi> = |psi>, and therefore
      <psi| A_h^dag R A_h |psi> = <A_h psi| R |A_h psi> = <psi| R |psi>
  identically, for ANY operator R whatsoever. The operator moves; no gauge-invariant state can
  notice. delta is zero by construction and would have been zero for any R, any carrier, any h.

AND THE DEEPER DEFECT: the minimal torus has ONE VERTEX, so A_h is a GLOBAL GAUGE TRANSFORMATION,
not transport along a path. Transport between places requires more than one place. This is the third
measurement today taken in a venue where the effect cannot appear -- the same class as O-29's abelian
carriers and AUDIT 1's one-qubit bath sites.

WHAT THE VALID VERSION REQUIRES:
  * a carrier with SEVERAL VERTICES, so a Wilson line between two of them is a real path; and
  * an observable that is not gauge-fixed away -- the HOLONOMY around a closed loop, or a state in a
    FLUX sector rather than the flat ground space.
The zeros below are retained as the record of the error, with the self-check that catches it."""
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
    nc=sum(1 for a in G for b in G if mul(a,b)!=mul(b,a))
    return dict(label=label,G=G,n=n,D=D,Ah=Ah,H=-(A+B),nc=nc)
def make_records(C, rng, want):
    """Records CONSTRUCTED on the eigenspaces of H, then every clause CHECKED. Never assumed."""
    D=C['D']; w,V=np.linalg.eigh(C['H'])
    blocks=[]
    for i,x in enumerate(w):
        if blocks and abs(x-blocks[-1][0])<1e-8: blocks[-1][1].append(i)
        else: blocks.append([x,[i]])
    if any(len(ix)%2 for _,ix in blocks): return [],blocks
    out=[]
    while len(out)<want:
        R=np.zeros((D,D),dtype=complex)
        for _,ix in blocks:
            m=len(ix); U=V[:,ix]
            Q=np.linalg.qr(rng.normal(size=(m,m)))[0]
            W=U@Q
            R+= W@np.diag([1.0]*(m//2)+[-1.0]*(m//2))@W.conj().T
        if np.linalg.norm(R@R-np.eye(D))>1e-8: continue
        if np.linalg.norm(R-R.conj().T)>1e-8: continue
        if np.linalg.norm(C['H']@R-R@C['H'])>1e-8: continue
        if any(abs(float(np.real(np.trace(V[:,ix].conj().T@R@V[:,ix]))))>1e-8 for _,ix in blocks): continue
        out.append(R)
    return out,blocks
say("="*104); say("THE FUNCTION: DOES TRANSPORT'S EFFECT ON A RECORD DEPEND ON THE CONFIGURATION?"); say("="*104)
rng=np.random.default_rng(11)
for C in (carrier([(1,0)],2,"D(Z_2)  abelian  (CONTROL)"),
          carrier([(1,2,3,0),(1,0,3,2)],4,"D(D_4)  NON-ABELIAN")):
    recs,blocks=make_records(C,rng,3)
    say(""); say("-"*104)
    say(f"  {C['label']}   dim {C['D']}   non-commuting pairs {C['nc']}/{C['n']**2}"
        f"   eigenspaces {[(round(float(v),3),len(ix)) for v,ix in blocks]}")
    say("-"*104)
    if not recs:
        say("    an eigenspace is odd-dimensional -- no record exists on this carrier (C-41)"); continue
    R0,R1,R2=recs
    say(f"    3 records constructed and verified against clauses (i)-(iv)")
    # commuting? we need a joint eigenbasis to 'write' several records at once
    c01=np.linalg.norm(R0@R1-R1@R0); c02=np.linalg.norm(R0@R2-R2@R0)
    say(f"    ||[R0,R1]|| = {c01:.3e}   ||[R0,R2]|| = {c02:.3e}")
    say("")
    say("  delta_0(psi,h) = <psi|A_h^dag R0 A_h|psi> - <psi|R0|psi>   -- how far transport moves record 0")
    say("  psi is drawn from the GROUND SPACE, then projected into the +-1 eigenspace of R1 (and R2),")
    say("  which is what 'writing' those records to a value means.")
    Pg,kdim=RecordModel(C['H'].astype(complex)).ground_space()
    say(f"    ground space dimension {kdim}")
    def project(P, vecs):
        out=[]
        for v in vecs:
            u=P@v
            nrm=np.linalg.norm(u)
            if nrm>1e-8: out.append(u/nrm)
        return out
    wg,Vg=np.linalg.eigh(Pg); base=[Vg[:,i] for i in range(len(wg)) if wg[i]>0.5]
    rows=[]
    for cfg,desc in (((None,None),"nothing else written"),
                     ((+1,None),"R1 = +1"),
                     ((-1,None),"R1 = -1"),
                     ((+1,+1),"R1 = +1, R2 = +1"),
                     ((+1,-1),"R1 = +1, R2 = -1"),
                     ((-1,-1),"R1 = -1, R2 = -1")):
        vs=base
        if cfg[0] is not None: vs=project((np.eye(C['D'])+cfg[0]*R1)/2, vs)
        if cfg[1] is not None: vs=project((np.eye(C['D'])+cfg[1]*R2)/2, vs)
        if not vs:
            rows.append((desc,0,None,None)); continue
        ds=[]
        for v in vs:
            for h in C['G']:
                a=C['Ah'][h]
                before=float(np.real(v.conj()@(R0@v)))
                after =float(np.real((a@v).conj()@(R0@(a@v))))/max(float(np.real((a@v).conj()@(a@v))),1e-12)
                ds.append(abs(after-before))
        rows.append((desc,len(vs),float(np.mean(ds)),float(np.max(ds))))
    say(f"    {'configuration':<26}{'states':>8}{'mean |delta_0|':>17}{'max |delta_0|':>16}")
    for desc,nv,mn,mx in rows:
        if mn is None: say(f"    {desc:<26}{nv:>8}{'(empty sector)':>17}{'':>16}")
        else: say(f"    {desc:<26}{nv:>8}{mn:>17.6f}{mx:>16.6f}")
    vals=[mn for _,_,mn,_ in rows if mn is not None]
    spread=(max(vals)-min(vals)) if vals else 0.0
    say("")
    say(f"    spread in mean |delta_0| across configurations: {spread:.6f}")
    # THE SELF-CHECK THAT SHOULD HAVE RUN FIRST: is transport acting trivially on these states?
    triv=max(np.linalg.norm(C['Ah'][h]@v - v) for v in base for h in C['G'])
    say(f"    SELF-CHECK  max ||A_h|psi> - |psi>|| over ground states and group elements: {triv:.3e}")
    if triv<1e-8:
        say("    -> TEST INVALID. Transport acts as the IDENTITY on every state used, because the")
        say("       ground space is gauge-invariant by construction. delta is zero for ANY operator")
        say("       whatsoever and NO conclusion may be drawn -- least of all from D(D_4), where")
        say("       C-43 measured ||[A_h,R]|| = 9.423 on the same carrier.")
    else:
        say(f"    -> {'THE EFFECT DEPENDS ON THE CONFIGURATION -- there is a function here' if spread>1e-6 else 'no configuration dependence over the states tested'}")
say(""); say("="*104); say("  READ -- from the numbers above"); say("="*104)
