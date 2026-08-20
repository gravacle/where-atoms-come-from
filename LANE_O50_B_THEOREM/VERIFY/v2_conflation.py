"""ADVERSARIAL VERIFY 2 -- THE AXIS-5 TEST.
Does the lane's "CANCELS" mean (A) "zero when averaged over the whole configuration space /
energy shell", or (B) "a given configuration yields no accumulation"?  The headline
   "no functional of the record configuration can be BOTH responsive AND non-cancelling
    -- which is exactly what a source must be"
only bites if (B).  Test (A) vs (B) with an explicit functional on the CANONICAL carrier.
"""
import numpy as np, itertools
from fractions import Fraction as F

print("="*100)
print("TEST 1.  Is the lane's F5 zero  sum_sigma [f(flip_i sigma) - f(sigma)] = 0  a")
print("         record-theoretic fact, or a bijection identity true for ANY function?")
print("="*100)
rng=np.random.default_rng(7)
worst=0
for k in range(1,8):
    cfg=list(itertools.product((1,-1),repeat=k))
    for trial in range(200):
        f={s:int(rng.integers(-50,50)) for s in cfg}      # ARBITRARY function, NO records anywhere
        for i in range(k):
            tot=sum(f[tuple(-x if j==i else x for j,x in enumerate(s))]-f[s] for s in cfg)
            worst=max(worst,abs(tot))
print("  arbitrary integer functions on {+-1}^k, k=1..7, 200 each, NO carrier, NO clause (i)-(v):")
print("     max |sum_sigma Delta_i f| =",worst)
print("  READ: the F5 zero holds for functions that have nothing to do with records.  It is the")
print("        statement that flip_i is a BIJECTION of a finite set.  It carries no record content.")
print("        The lane's non-zero CONTROL column is produced only by re-weighting with a")
print("        flip-NON-invariant d -- i.e. by breaking the bijection, not by finding a source.")

print()
print("="*100)
print("TEST 2.  A functional of the record configuration that is RESPONSIVE and ACCUMULATES.")
print("         Carrier: the CANONICAL one -- toric code, g tori, k = 2g genuine records.")
print("="*100)
print("  f(sigma) = sum_i sigma_i     (the record 'magnetisation')")
print()
print("   g   k=2g   mean over C   Pi_GW f   value at all-+   |Delta_i f| at every sigma   'cancels' by the lane   accumulates at a CONFIG")
print("  "+"-"*140)
for g in range(1,9):
    k=2*g
    cfg=list(itertools.product((1,-1),repeat=k))
    f=lambda s: sum(s)
    mean=F(sum(f(s) for s in cfg),len(cfg))
    # G_W is the full flip group on the toric code (verified in v1), so Pi_GW f = mean = const
    pi=mean
    allplus=f(tuple([1]*k))
    delt=set(abs(f(tuple(-x if j==0 else x for j,x in enumerate(s)))-f(s)) for s in cfg)
    print("  %2d   %4d   %11s   %7s   %14d   %-27s   %-21s  %s"
          %(g,k,mean,pi,allplus,sorted(delt),"YES (shell avg = 0)","YES: value = %d, grows as 2g"%allplus))
print()
print("  READ: f is a FUNCTIONAL OF THE RECORD CONFIGURATION.  It is RESPONSIVE -- every single")
print("        admissible write changes it by exactly 2.  Its shell average is EXACTLY 0, so by")
print("        THEOREM 2 it is 100% 'cancelling' (Pi_GW f = 0, f - Pi_GW f = f).  And yet at the")
print("        realised configuration it equals 2g and GROWS WITHOUT BOUND with the number of")
print("        records.  RESPONSIVE *and* ACCUMULATING-AT-A-CONFIGURATION, simultaneously.")

print()
print("="*100)
print("TEST 3.  The same on the operator side: 'cancels' == 'traceless per shell' == clause (iv).")
print("="*100)
import sys
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_B_THEOREM/VERIFY")
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def kron(ms):
    M=np.array([[1]],dtype=complex)
    for m in ms: M=np.kron(M,m)
    return M
def pau(n,s):
    d={'X':X,'Z':Z,'I':I2}; return kron([d[c] for c in s])
n=8
S=['XIXIIIII','IZIZIIII']   # the two records found by search in v1 (R_1, R_2)
R1=pau(n,'XIXIIIII'); R2=pau(n,'IZIZIIII')
A=R1+R2                      # the operator whose eigenvalue IS f(sigma)
L=2
def h(i,j): return 2*((i%L)*L+(j%L))
def v(i,j): return 2*((i%L)*L+(j%L))+1
St=[];Pl=[]
for i in range(L):
    for j in range(L):
        s=['I']*n
        for e in (h(i,j),h(i-1,j),v(i,j),v(i,j-1)): s[e]='X'
        St.append(''.join(s))
        p=['I']*n
        for e in (h(i,j),h(i,j+1),v(i,j),v(i+1,j)): p[e]='Z'
        Pl.append(''.join(p))
H=-sum(pau(n,s) for s in St+Pl)
w,V=np.linalg.eigh(H)
shells=[]; i=0
while i<len(w):
    j=i
    while j+1<len(w) and abs(w[j+1]-w[i])<1e-8: j+=1
    shells.append((w[i],V[:,i:j+1])); i=j+1
print("  shell   E        dim     Tr(P_E A)   <sigma|A|sigma> for the joint eigenstates in that shell")
for e,C in shells:
    P_=C@C.conj().T
    tr=np.trace(P_@A).real
    # joint eigenvalues of R1,R2 inside the shell
    vals=set()
    Cc=C
    groups=[Cc]
    for R in (R1,R2):
        ng=[]
        for Cx in groups:
            Rs=Cx.conj().T@R@Cx; ws,Vs=np.linalg.eigh(Rs)
            for s_ in (1,-1):
                idx=[t for t in range(len(ws)) if (ws[t]>0)==(s_>0)]
                if idx: ng.append(Cx@Vs[:,idx])
        groups=ng
    for Cx in groups:
        a=(Cx.conj().T@A@Cx)
        vals.add(round(np.real(np.trace(a)/Cx.shape[1]),6))
    print("  %5.1f  %7.3f  %5d  %11.6f   %s"%(0,e,C.shape[1],tr,sorted(vals)))
print()
print("  READ: Tr(P_E A) = 0 on EVERY shell -- the lane's 'it cancels'.  But the eigenvalues of A")
print("        on the joint record blocks are -2, 0, 0, +2: a definite record configuration gives")
print("        a NON-ZERO, configuration-dependent value.  'Cancels' here means 'traceless on the")
print("        shell', i.e. the expectation in the MAXIMALLY MIXED state.  That is literally")
print("        clause (iv) itself.  THEOREM 2 restates the definition of a record; it does not")
print("        forbid an accumulating configuration-dependent quantity.")
