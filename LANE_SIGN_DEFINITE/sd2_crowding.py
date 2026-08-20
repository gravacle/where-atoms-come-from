"""SIGN-DEFINITENESS OF THE CROWDING SUPPRESSION -- the candidate that could accumulate.

sd1 showed the integer intersection pairing is TWO-SIGNED and CANCELS, with |sum|/sum|.| falling
0.3333 -> 0.2222 -> 0.1698 as records are added. It SCREENS, exactly like electric charge, which is
why it cannot be gravity's source.

The crowding suppression is the other candidate, and its sign WAS genuinely in question: sharing an
environment site could in principle HELP a record as easily as hurt it -- constructive interference
in the joint evolution is not forbidden by anything. If it turns out never to help, it is
sign-definite, it accumulates rather than cancelling, and it is the first quantity in this program
with the FORM gravity requires.

EXCLUDED AS TRIVIAL: norms, squares, and magnitudes-of-differences are non-negative by construction
and prove nothing. The quantity here is a SIGNED DIFFERENCE b - a of two chi values, either of which
could be the larger."""
import sys, os, numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','model'))
from record_model import RecordModel, Environment, symplectic_logicals, xz_to_matrix
def say(*a): print(*a); sys.stdout.flush()
ENERGIES=(1.0,1.4,0.7,1.2,0.9,1.6,1.1,0.8)
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); I2=np.eye(2)
def paulis(n,word):
    M=np.array([[1]],dtype=complex)
    for c in word: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
say("="*104); say("IS THE CROWDING SUPPRESSION SIGN-DEFINITE?"); say("="*104)
say("  signed difference  d = chi(spread) - chi(crowded).  d > 0 means sharing the site COST the")
say("  record; d < 0 would mean sharing HELPED it. Both are a priori possible.")
say("")
TS=np.linspace(1.0,13.0,13); NB=4
say(f"  {'n':>4}{'k':>4}{'lam':>7}{'pairs':>8}{'d > 0 (cost)':>14}{'d < 0 (helped)':>16}{'min d':>11}{'max d':>11}{'|sum|/sum|d|':>14}")
for n in (4,6):
    stab=[[1]*n+[0]*n,[0]*n+[1]*n]
    prs=symplectic_logicals(stab,n); k=len(prs)
    M=RecordModel(-(paulis(n,'X'*n)+paulis(n,'Z'*n))); nS=M.n
    flat=[xz_to_matrix(v,n) for pr in prs for v in pr]
    for lam in (0.4,0.8,1.2):
        env=Environment(NB, energies=ENERGIES[:NB])
        ds=[]
        for i in range(min(len(flat),3)):
            Ri=flat[i]
            for j in range(len(flat)):
                if j==i: continue
                Pj=flat[j]
                a=float(np.mean([env.holevo(M.evolve([(Ri,0),(Pj,0)],env,lam=lam,t=t),Ri,nS) for t in TS]))
                b=float(np.mean([env.holevo(M.evolve([(Ri,0),(Pj,1)],env,lam=lam,t=t),Ri,nS) for t in TS]))
                if max(a,b)<1e-9: continue
                ds.append(b-a)
        if not ds:
            say(f"  {n:>4}{k:>4}{lam:>7.2f}{0:>8}{'(no valid pairs -- no conclusion)':>57}"); continue
        pos=sum(1 for d in ds if d>1e-9); neg=sum(1 for d in ds if d<-1e-9)
        S=sum(ds); SA=sum(abs(d) for d in ds)
        say(f"  {n:>4}{k:>4}{lam:>7.2f}{len(ds):>8}{pos:>14}{neg:>16}{min(ds):>11.6f}{max(ds):>11.6f}"
            f"{(abs(S)/SA if SA else float('nan')):>14.4f}")
say("")
say("  A ratio |sum d| / sum|d| of 1.0000 means EVERY term has the same sign -- it ACCUMULATES.")
say("  A ratio well below 1 means the terms CANCEL, as the integer intersection pairing does")
say("  (0.3333, 0.2222, 0.1698 at k = 2, 4, 6 -- cancelling MORE as records are added, which is")
say("  screening).")
say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
