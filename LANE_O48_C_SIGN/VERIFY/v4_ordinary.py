"""
VERIFY-4.  IS ANY OF THIS MORE THAN THE TEXTBOOK 1D ISING CHAIN / CLASSICAL REPETITION CODE?

(a) The 'searched, never nominated' writer set: solve for it in closed form and compare to the
    exhaustive 4^n search.  Also MEASURE the energy change the lane hard-coded as dE = 0
    (s1_clauses.py line: `dE = 0 if cnt else None`).
(b) I(Z_0..Z_{k-1} ; E) = k-1 : test whether it is just the entropy of the k-1 bond variables
    t_0..t_{k-2} that k consecutive spins determine.
(c) The 'decaying' family that 'keeps coherence 0.60': is that coherence, or one dominant bond?
"""
import itertools, math, random
import numpy as np
from fractions import Fraction

print("="*104)
print("VERIFY-4a   THE ADMISSIBLE WRITER SET IN CLOSED FORM, AND THE ENERGY CHANGE ACTUALLY MEASURED")
print("="*104)
print()
print("  A Pauli (x,z) commutes with Z_i Z_{i+1} iff x_i + x_{i+1} = 0 mod 2.  On a connected chain")
print("  that forces x = 0 (Z-type) or x = 1...1.  Anticommuting with Z_0 needs x_0 = 1, so the")
print("  admissible flippers of Z_0 are EXACTLY {X-or-Y on every site} : 2^n of them, weight n.")
print()
I2=np.eye(2,dtype=complex); X2=np.array([[0,1],[1,0]],dtype=complex)
Y2=np.array([[0,-1j],[1j,0]],dtype=complex); Z2=np.array([[1,0],[0,-1]],dtype=complex)
PA=[I2,X2,Y2,Z2]
def kron(ops):
    o=np.array([[1.0+0j]])
    for a in ops: o=np.kron(o,a)
    return o
print(f"  {'n':>3} {'closed form 2^n':>16} {'exhaustive 4^n search':>22} {'agree?':>7} "
      f"{'min weight':>11} {'MEASURED max |E(psi)-E(W psi)|':>31}")
for n in (3,4,5,6,7):
    J=[i+1 for i in range(n-1)]
    cnt=0; mw=None
    for code in itertools.product(range(4),repeat=n):
        x=[1 if c in (1,2) else 0 for c in code]; z=[1 if c in (2,3) else 0 for c in code]
        ok=all((x[i]+x[i+1])%2==0 for i in range(n-1))
        if not ok: continue
        if x[0]!=1: continue
        cnt+=1; w=sum(1 for c in code if c); mw=w if mw is None else min(mw,w)
    # measure the energy change of the minimum-weight writer on every eigenstate, densely
    dEmax=0.0
    if n<=7:
        Hd=np.zeros((2**n,2**n),dtype=complex)
        for i,j in enumerate(J):
            ops=[I2]*n; ops[i]=Z2; ops[i+1]=Z2; Hd+=j*kron(ops)
        W=kron([X2]*n)
        Ediag=np.real(np.diag(Hd))
        WE=np.real(np.diag(W.conj().T@Hd@W))
        dEmax=float(np.max(np.abs(Ediag-WE)))
    print(f"  {n:>3} {2**n:>16} {cnt:>22} {str(cnt==2**n):>7} {mw:>11} {dEmax:>31.6f}")
print()
print("  READ: the exhaustive 4^n search returns exactly the closed-form set.  The dE column the")
print("  lane printed as 0 was hard-coded (`dE = 0 if cnt else None`), not measured; measured here")
print("  it is 0 -- but it is 0 because [W,H]=0, i.e. by clause (ii), for ANY H and ANY admissible")
print("  W.  'Records flip for free' is the definition of admissible, not a result about records.")
print()

print("="*104)
print("VERIFY-4b   IS I(k records ; E) = k-1 ANYTHING BUT THE k-1 BOND VARIABLES THEY DETERMINE?")
print("="*104)
print()
def spin_table(n):
    xs=np.arange(2**n)
    return np.array([1-2*((xs>>i)&1) for i in range(n)],dtype=np.int8).T
def mi(a,b,N):
    _,ai=np.unique(a,return_inverse=True); _,bi=np.unique(b,return_inverse=True)
    joint=np.zeros((ai.max()+1,bi.max()+1),dtype=np.int64); np.add.at(joint,(ai,bi),1)
    pj=joint/N; pc=pj.sum(1,keepdims=True); pe=pj.sum(0,keepdims=True)
    with np.errstate(divide="ignore",invalid="ignore"): t=pj*np.log2(pj/(pc*pe))
    return float(np.where(joint>0,t,0.0).sum())
print(f"  {'family':<10} {'n':>3} " + " ".join(f"{'k='+str(k):>8}" for k in (1,2,3,4,5))
      + f"   {'= H(t_0..t_{k-2})?':>20}")
for fam,mk in (("superinc",lambda m:[2**i for i in range(m)]),
               ("uniform", lambda m:[1]*m)):
    for n in (11,13):
        J=mk(n-1); S=spin_table(n); N=2**n
        E=np.zeros(N,dtype=np.int64)
        for i in range(n-1): E+=J[i]*S[:,i].astype(np.int64)*S[:,i+1].astype(np.int64)
        cells=[]; match=True
        for k in (1,2,3,4,5):
            code=np.zeros(N,dtype=np.int64)
            for i in range(k): code=code*2+(S[:,i]>0).astype(np.int64)
            v=mi(code,E,N); cells.append(f"{v:>8.5f}")
            # the bond variables the same k spins determine
            bond=np.zeros(N,dtype=np.int64)
            for i in range(k-1): bond=bond*2+((S[:,i]*S[:,i+1])>0).astype(np.int64)
            vb=mi(bond,E,N)
            if abs(v-vb)>1e-12: match=False
        print(f"  {fam:<10} {n:>3} " + " ".join(cells) + f"   {str(match):>20}")
print()
print("  READ: the joint record reading buys EXACTLY what the k-1 bond variables buy, and nothing")
print("  else -- I(Z_0..Z_{k-1};E) = I(t_0..t_{k-2};E) identically.  k spins on an open chain are")
print("  (k-1 bond variables) + (1 global sign); H depends only on bonds; the global sign is the")
print("  Z_2 symmetry.  That is the whole 'I = k-1 EXACT' law.")
print()

print("="*104)
print("VERIFY-4c   THE 'decaying' FAMILY THAT 'KEEPS COHERENCE 0.60'")
print("="*104)
print()
print(f"  {'m':>4} {'J_0':>7} {'sum of the rest':>16} {'J_0 > rest?':>12} {'exact mean|E|':>14} "
      f"{'J_0':>7} {'equal?':>7} {'coh':>9} {'J_0/M':>9}")
def dist(J):
    R=sum(abs(j) for j in J); off=R; cur=[0]*(2*R+1); cur[off]=1; lo=hi=off
    for j in J:
        a=abs(j); nxt=[0]*(2*R+1)
        for i in range(lo,hi+1):
            c=cur[i]
            if c: nxt[i+a]+=c; nxt[i-a]+=c
        lo-=a; hi+=a; cur=nxt
    return {i-off:c for i,c in enumerate(cur) if c}
for m in (8,16,32,64,128):
    J=[max(1,5040//((i+1)**2)) for i in range(m)]
    M=sum(J); rest=M-J[0]
    d=dist(J); tot=sum(d.values())
    ma=Fraction(sum(abs(e)*c for e,c in d.items()),tot)
    print(f"  {m:>4} {J[0]:>7} {rest:>16} {str(J[0]>rest):>12} {float(ma):>14.4f} "
          f"{J[0]:>7} {str(ma==J[0]):>7} {float(ma/M):>9.6f} {J[0]/M:>9.6f}")
print()
print("  READ: for the decaying family J_0 alone exceeds the sum of all other couplings, so the")
print("  sign of E is set by ONE bond and mean|E| equals J_0 EXACTLY.  coh = J_0/M is the fraction")
print("  of the l1 norm sitting in a single term.  That is not accumulation across terms; it is a")
print("  one-term sum with a decaying tail, and it is the reason the same family is not extensive.")
