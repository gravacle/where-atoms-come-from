"""LEG D -- THE DECLARED WEAKNESS: T IS NOT UNIQUE.  SETTLE IT.
   Three parts:
   D1  the exact solution set of  A unitary, A^L = M_gamma  (characterised, then sampled)
   D2  the exact invisibility criterion, proved and verified
   D3  a CANONICITY THEOREM for T under locality  -- which cuts FOR the registrar
"""
import numpy as np, sys
sys.path.insert(0,'.')
from wcore import *
np.set_printoptions(precision=9, linewidth=200)

f, c = 1.0, np.sqrt(2.0)
k = K1(f, c)
pA = np.array([0.40,0.15,0.15,0.15,0.15]); sA = np.sqrt(pA).astype(complex)
pB = np.array([0.40,0.25,0.05,0.02,0.28]); sB = np.sqrt(pB).astype(complex)
sC = np.sqrt(pA)*np.exp(1j*np.array([0.0,0.7,-1.9,2.3,0.4]))
S = [sA,sB,sC]
cls = k.classes(); order=[(0,0),(1,0),(0,1),(1,1)]
w = np.exp(2j*np.pi/3)

def spread(AF, AC, nmax=24):
    """max over n<=nmax of the spread of |<AF^n s, AC^n s>| across the three equal-pi states"""
    a=[s.copy() for s in S]; b=[s.copy() for s in S]; worst=0.0
    for n in range(nmax):
        a=[AF@x for x in a]; b=[AC@x for x in b]
        vals=np.array([abs(np.vdot(a[i],b[i])) for i in range(3)])
        worst=max(worst, np.ptp(vals))
    return worst

def diag_classconst(Q, tol=1e-10):
    off = np.linalg.norm(Q - np.diag(np.diag(Q)))
    if off > tol: return False, off, None
    d = np.diag(Q); bad = 0.0
    for o in order:
        idx = cls[o]
        if len(idx) > 1: bad = max(bad, np.max(np.abs(d[idx]-d[idx[0]])))
    return bad <= tol, off, bad

# ---------------------------------------------------------------- D1
print("D1  THE SOLUTION SET OF  A unitary, A^3 = M_gamma  ON K1")
print("    A commutes with A^3 = M, so A preserves M's eigenspaces.  For gamma_F,")
print("    M_F = diag(W_F,W_F,W_F,1,1): A = A|_{span e0,e1,e2} (+) A|_{span e3,e4},")
print("    with A|_1^3 = W_F I_3 and A|_2^3 = I_2.  Sampled below; 8 free real parameters")
print("    generically (U(3)/U(1)^3 = 6, U(2)/U(1)^2 = 2) plus discrete cube-root labels.\n")

def rand_root(car, which, rng, force_diag=False, force_classconst=False):
    onloop = car.inF if which=='F' else car.inC
    ang = car.f if which=='F' else car.c
    Wr = np.exp(1j*ang/3.0)                    # a principal cube root of W
    idx1 = np.where(onloop)[0]; idx2 = np.where(~onloop)[0]
    A = np.zeros((car.nV,car.nV), complex)
    for idx, base in ((idx1, Wr), (idx2, 1.0+0j)):
        m = len(idx)
        if m == 0: continue
        labs = rng.integers(0,3,m)
        D = base*np.power(w, labs)
        if force_classconst:                    # one label for the whole block
            D = base*np.power(w, rng.integers(0,3))*np.ones(m)
        if force_diag:
            A[np.ix_(idx,idx)] = np.diag(D)
        else:
            X = rng.normal(size=(m,m)) + 1j*rng.normal(size=(m,m))
            V, _ = np.linalg.qr(X)
            A[np.ix_(idx,idx)] = V @ np.diag(D) @ V.conj().T
    return A

rng = np.random.default_rng(11235813)
NS = 400
res = {'generic pair (random cube roots, both loops)':[],
       'both DIAGONAL, class-constant labels':[],
       'both DIAGONAL, per-vertex labels (NOT class-constant)':[]}
for _ in range(NS):
    AF = rand_root(k,'F',rng); AC = rand_root(k,'C',rng)
    assert np.linalg.norm(np.linalg.matrix_power(AF,3)-k.M('F'))<1e-9
    assert np.linalg.norm(np.linalg.matrix_power(AC,3)-k.M('C'))<1e-9
    res['generic pair (random cube roots, both loops)'].append(spread(AF,AC))
for _ in range(NS):
    AF = rand_root(k,'F',rng,force_diag=True,force_classconst=True)
    AC = rand_root(k,'C',rng,force_diag=True,force_classconst=True)
    res['both DIAGONAL, class-constant labels'].append(spread(AF,AC))
for _ in range(NS):
    AF = rand_root(k,'F',rng,force_diag=True); AC = rand_root(k,'C',rng,force_diag=True)
    res['both DIAGONAL, per-vertex labels (NOT class-constant)'].append(spread(AF,AC))
print("    %d draws each, seed 11235813.  'spread' = max over n<=24 of the |Z_n| spread across")
print("    the three ready states that share pi exactly.")
for nm,v in res.items():
    v=np.array(v)
    print("      %-52s  min %.2e  med %.2e  max %.2e   #(<1e-12)=%d/%d" %
          (nm, v.min(), np.median(v), v.max(), int((v<1e-12).sum()), len(v)))

# ---------------------------------------------------------------- named alternatives
print("\n    NAMED members of the solution set, each with A^3 = M exactly:")
def report(nm, AF, AC):
    ok = all(np.linalg.norm(np.linalg.matrix_power(X,3)-M)<1e-9 for X,M in ((AF,k.M('F')),(AC,k.M('C'))))
    Q1 = np.linalg.inv(AF)@AC
    dcc,off,bad = diag_classconst(Q1)
    print("      %-46s A^3=M:%s  Q_1 diag+class-const:%-5s  spread=%.3e" % (nm, ok, dcc, spread(AF,AC)))
TF, TC = k.T('F'), k.T('C')
report("COR-F edge tick T (moves one edge)", TF, TC)
DF = np.diag(np.where(k.inF, np.exp(1j*k.f/3), 1.0).astype(complex))
DC = np.diag(np.where(k.inC, np.exp(1j*k.c/3), 1.0).astype(complex))
report("'smeared holonomy' root diag(W^{1/3}) on loop", DF, DC)
# a diagonal root that is NOT class-constant: put omega on v1 and omega^2 on v2
DF2 = DF@np.diag([1,w,w**2,1,1]); DC2 = DC@np.diag([1,1,1,w,w**2])
report("diagonal root, labels differing within a class", DF2, DC2)
# conjugated tick
V = np.diag([1,1,1,1,1]).astype(complex); V[1,1]=np.exp(1j*0.9); V[2,2]=np.exp(-1j*0.3)
report("V T_F V* / T_C   (V diagonal, commutes with M_F)", V@TF@V.conj().T, TC)
print("      T_C = T_F * D for some diagonal D ?  ", np.linalg.norm(np.linalg.inv(TF)@TC - np.diag(np.diag(np.linalg.inv(TF)@TC))) < 1e-10)

# ---------------------------------------------------------------- D2
print("\nD2  THE EXACT INVISIBILITY CRITERION")
print("    <A_F^n s, A_C^n s> = <s, Q_n s> with Q_n = (A_F^n)^* A_C^n = A_F^{-n} A_C^n.")
print("    <s,Q s> is a function of pi(s) alone for every s  <=>  Q is diagonal with")
print("    class-constant diagonal.  (=>: off-diagonals are read by phase changes that fix")
print("    every |s_v|; within-class variation of the diagonal is read by moving weight")
print("    inside a class.  <=: <s,Q s> = sum_ab q_ab p_ab.)   VERIFIED, 600 random pairs:")
rng2 = np.random.default_rng(271828)
agree = 0; tot = 0; mism = []
for _ in range(600):
    AF = rand_root(k,'F',rng2, force_diag=rng2.random()<0.5, force_classconst=rng2.random()<0.5)
    AC = rand_root(k,'C',rng2, force_diag=rng2.random()<0.5, force_classconst=rng2.random()<0.5)
    allD = True
    QF = np.linalg.inv(AF); Qn = np.eye(5,dtype=complex); AFn=np.eye(5,dtype=complex); ACn=np.eye(5,dtype=complex)
    for n in range(1,13):
        AFn = AFn@AF; ACn = ACn@AC
        Qn = np.linalg.inv(AFn)@ACn
        if not diag_classconst(Qn)[0]: allD = False; break
    sp = spread(AF,AC)
    pred_zero = allD
    obs_zero  = sp < 1e-12
    tot += 1; agree += int(pred_zero == obs_zero)
    if pred_zero != obs_zero: mism.append((allD, sp))
print("      criterion agrees with observation on %d of %d pairs; mismatches: %s" % (agree, tot, mism[:5]))

# ---------------------------------------------------------------- D3
print("\nD3  CANONICITY THEOREM FOR T -- AND IT CUTS FOR THE REGISTRAR, NOT AGAINST IT")
print("    THEOREM (proved here, verified below).  Let gamma be a cycle w_0->...->w_{L-1}->w_0,")
print("    L>=3.  Let A be unitary, identity off gamma, with support inside gamma contained in")
print("    the diagonal together with the loop's own directed edges (each fibre value either")
print("    stays put or moves one edge forward).  Then A is EITHER diagonal OR a pure shift.")
print("    PROOF. On the loop A[j,j]=d_j, A[j+1,j]=t_j, all else 0.  Columns j and j+1 overlap")
print("    only in row j+1, so orthogonality gives conj(t_j) d_{j+1} = 0 for every j.  If some")
print("    t_j != 0 then d_{j+1}=0, so |t_{j+1}|=1 by normalisation, so t_{j+1} != 0; going")
print("    round the cycle kills every d and forces |t_j|=1 for all j -- a pure shift.  Else")
print("    every t_j = 0 and A is diagonal.  QED  -- no intermediate 'partly moving' unitary.")
print("    COROLLARY.  Adding A^L = M_gamma: the shift branch has prod t_j = W_gamma, and if")
print("    each t_j is required to be a function of its OWN edge's transport U_e alone then")
print("    t_j = zeta*U_{e_j} with zeta^L = 1: T IS UNIQUE UP TO AN L-TH ROOT OF UNITY.")
print("    The diagonal branch is the 'smeared holonomy' family diag(W^{1/L}) -- which needs")
print("    the WHOLE circuit to write down, so it is NOT built from edge-local data.")
print("    => among LOCAL transports there are exactly TWO conventions, MOVE and DON'T MOVE,")
print("       and only MOVE is edge-local.  Verification, 4000 random local unitaries:")
rng3 = np.random.default_rng(31415926)
bad = 0; found_mixed = 0
for _ in range(4000):
    d = rng3.normal(size=3)+1j*rng3.normal(size=3); t = rng3.normal(size=3)+1j*rng3.normal(size=3)
    # try to build a unitary with this support pattern by projecting: just test the algebraic condition
    # instead: enumerate the two families and confirm nothing else satisfies unitarity
    A = np.zeros((3,3),complex)
    for j in range(3):
        A[j,j] = d[j]; A[(j+1)%3, j] = t[j]
    # normalise columns then test unitarity
    A = A/np.linalg.norm(A,axis=0)
    if np.linalg.norm(A.conj().T@A-np.eye(3)) < 1e-9:
        nd = np.linalg.norm(np.diag(np.diag(A))); ns = np.linalg.norm(A-np.diag(np.diag(A)))
        if nd > 1e-9 and ns > 1e-9: found_mixed += 1
print("      random support-pattern matrices that are unitary AND mixed (both d and t nonzero): %d" % found_mixed)
print("      (the theorem says 0; a random draw is measure-zero on the unitary locus, so this is")
print("       a consistency check only -- the proof above is the content.)")
# direct constructive check of the two families
ok1 = np.linalg.norm(TC.conj().T@TC-np.eye(5))<1e-12
ok2 = np.linalg.norm(DC.conj().T@DC-np.eye(5))<1e-12
print("      pure-shift family unitary: %s    diagonal family unitary: %s" % (ok1, ok2))
print("      rate check on the diagonal ('smeared') convention -- invisibility holds, but the")
print("      per-tick rate is NOT m(P)/3 either:")
for nm, s in (('sA',sA),('sB',sB),('sC',sC)):
    a=s.copy(); b=s.copy(); acc=0.0; NN=200000
    dF=np.diag(DF); dC=np.diag(DC)
    n=np.arange(1,NN+1)
    # closed form: Z_n = sum_v conj(dF_v)^n dC_v^n |s_v|^2
    q = np.conj(dF)*dC
    Zn = (np.abs(s)**2)[None,:] * np.power(q[None,:], n[:,None])
    Z = Zn.sum(axis=1)
    print("        %s  per-tick rate = %.9f   (m(P) = %.9f,  m(P)/3 = %.9f)" %
          (nm, float(np.mean(np.log(np.abs(Z)))), mahler4(*k.pi(sA)), mahler4(*k.pi(sA))/3))
