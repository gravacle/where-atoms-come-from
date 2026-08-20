"""ADVERSARIAL: TRY HARD TO MAKE CROWDING *HELP* A RECORD.

sd2 found the crowding suppression sign-definite -- 27 of 27 pairs at n=4 cost the record chi, none
helped, |sum|/sum|d| = 1.0000 at every coupling. That is the form gravity needs.

BUT THIS PROGRAM MADE EXACTLY THIS MISTAKE TODAY. O-32 registered chi-per-site as 'the first quantity
that knows how much is enclosed'; AUDIT 1 showed it was the Holevo bound, and C-37 was withdrawn. If
the Holevo bound FORCES d > 0 -- capacity can only be consumed, never created -- then sign-definiteness
here is conventional physics and it SATURATES rather than accumulating without bound.

THE HONEST TEST OF A SIGN-DEFINITENESS CLAIM IS TO TRY TO BREAK IT. This lane searches hard for any
configuration with d < 0, sweeping every knob at once: coupling strength over two orders of magnitude,
bath size, bath energies, temperature, time window, and which records are paired. A single negative
falsifies sign-definiteness. Sustained failure to find one, across a wide sweep, is evidence the sign
is FORCED -- which is itself the answer, just not the exciting one.

Also measured: whether the accumulated suppression SATURATES with bath size, which is what a capacity
bound must do and what a gravitational source must NOT do."""
import sys, os, itertools, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel, Environment, symplectic_logicals, xz_to_matrix
def say(*a): print(*a); sys.stdout.flush()
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); I2=np.eye(2)
def paulis(n,word):
    M=np.array([[1]],dtype=complex)
    for c in word: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
say("="*104); say("ADVERSARIAL: CAN CROWDING EVER *HELP* A RECORD?"); say("="*104)
n=4
stab=[[1]*n+[0]*n,[0]*n+[1]*n]
prs=symplectic_logicals(stab,n)
M=RecordModel(-(paulis(n,'X'*n)+paulis(n,'Z'*n))); nS=M.n
flat=[xz_to_matrix(v,n) for pr in prs for v in pr]
say(f"  carrier [[4,2,2]], dim {nS}, {len(flat)} record operators from {len(prs)} conjugate pairs")
say("")
say("1. WIDE SWEEP -- every knob at once. A SINGLE d < 0 falsifies sign-definiteness.")
LAMS=[0.05,0.1,0.2,0.4,0.8,1.2,2.0,3.0,5.0]
BETAS=[0.1,0.5,2.0,8.0]
NBS=[2,3,4]
ENERGY_SETS=[(1.0,1.4,0.7,1.2),(0.1,0.2,0.3,0.4),(3.0,0.1,2.0,0.5),(1.0,1.0,1.0,1.0)]
WINDOWS=[np.linspace(0.2,2.0,9),np.linspace(1.0,13.0,9),np.linspace(10.0,60.0,9)]
tot=0; neg=0; worst=1e9; negex=[]
for lam in LAMS:
    for beta in BETAS:
        for nb in NBS:
            for es in ENERGY_SETS:
                env=Environment(nb, energies=es[:nb], beta=beta)
                for TS in WINDOWS:
                    for i in (0,1):
                        Ri=flat[i]
                        for j in range(len(flat)):
                            if j==i: continue
                            Pj=flat[j]
                            a=float(np.mean([env.holevo(M.evolve([(Ri,0),(Pj,0)],env,lam=lam,t=t),Ri,nS) for t in TS]))
                            b=float(np.mean([env.holevo(M.evolve([(Ri,0),(Pj,1%nb)],env,lam=lam,t=t),Ri,nS) for t in TS]))
                            if max(a,b)<1e-9: continue
                            d=b-a; tot+=1
                            if d<worst: worst=d
                            if d<-1e-9:
                                neg+=1
                                if len(negex)<5: negex.append((lam,beta,nb,es,round(float(TS[0]),2),i,j,d))
say(f"   configurations tested: {tot}")
say(f"   configurations where crowding HELPED (d < 0): {neg}")
say(f"   most negative d found: {worst:.6e}")
if negex:
    say("   examples where it helped:")
    for e in negex: say(f"     lam={e[0]} beta={e[1]} nb={e[2]} E={e[3]} t0={e[4]} rec={e[5]} partner={e[6]} d={e[7]:.6e}")
say(f"   -> {'SIGN-DEFINITENESS FALSIFIED -- crowding can help' if neg else 'NO COUNTEREXAMPLE FOUND across the whole sweep: the sign appears FORCED'}")
say("")
say("2. DOES THE ACCUMULATED SUPPRESSION SATURATE WITH BATH SIZE?")
say("   A capacity bound MUST saturate. A gravitational source must NOT.")
say(f"   {'bath qubits':>13}{'site dim':>10}{'total suppression over pairs':>31}{'per pair':>11}")
prev=None
for nb in (1,2,3,4,5):
    env=Environment(nb, energies=tuple([1.0,1.4,0.7,1.2,0.9][:nb]))
    TS=np.linspace(1.0,13.0,9); tots=[]
    for i in (0,1):
        Ri=flat[i]
        for j in range(len(flat)):
            if j==i: continue
            Pj=flat[j]
            a=float(np.mean([env.holevo(M.evolve([(Ri,0),(Pj,0)],env,lam=0.8,t=t),Ri,nS) for t in TS]))
            b=float(np.mean([env.holevo(M.evolve([(Ri,0),(Pj,min(1,nb-1))],env,lam=0.8,t=t),Ri,nS) for t in TS]))
            if max(a,b)<1e-9: continue
            tots.append(b-a)
    if not tots: say(f"   {nb:>13}{2**nb:>10}{'(no valid pairs)':>31}"); continue
    S=sum(tots)
    say(f"   {nb:>13}{2**nb:>10}{S:>31.6f}{S/len(tots):>11.6f}")
say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
