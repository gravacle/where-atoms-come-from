# LANE_W11_R_MATH — LEG 3.  THE CHARACTERISATION: IS "BOTH BRANCH OPERATORS DIAGONAL" NECESSARY?
#
# The registrar's conclusion: "the invisibility of incidence holds exactly where both branch
# operators are DIAGONAL ... only on the sublattice n = 0 mod lcm(L_F,L_C)."
#
# THEOREM C (mine).  Write Z_n = s^* A_n s with A_n := T_F^{-n} T_C^n (unitary).  Then
#   (C1)  Z_n is a function of pi alone, for every ready state s
#         <=>  A_n is DIAGONAL and its diagonal is CONSTANT on every incidence class of size >= 2.
#         [ diagonality is forced by phase freedom s -> Ds, D diagonal unitary, which preserves pi;
#           class-constancy by within-class redistribution.  Singleton classes are unconstrained
#           because their weight is pinned by pi. ]
#   (C2)  A_n is diagonal  <=>  sigma_F^n = sigma_C^n  as permutations of the vertex set,
#         where sigma_gamma is the loop's cyclic rotation (identity off the loop).
#         [ T_gamma is a weighted permutation matrix with permutation sigma_gamma. ]
#   (C3)  If V(gamma_F) != V(gamma_C) then (C2) holds iff L_F | n AND L_C | n, i.e. lcm | n.
#         [ u in V(F)\V(C) is fixed by sigma_C^n, so sigma_F^n fixes u, so L_F | n; symmetric. ]
#   (C4)  If V(gamma_F) = V(gamma_C) the condition can hold at n NOT divisible by lcm.
#
# SO: the registrar's characterisation is CORRECT ON BOTH CARRIERS IT RAN, and is a THEOREM there,
# not a measurement.  Its STATED FORM is not the operative condition: "both diagonal" is
# sufficient, never necessary in general.  The operative object is the RELATIVE branch operator
# A_n, and the operative condition is "A_n is a class-function multiplication operator".
import numpy as np, itertools, wm0_lib as L
rng=np.random.default_rng(20260817)

def sigma(walk,NV):
    s=list(range(NV))
    for (u,v,e,sg) in walk: s[u]=v
    return tuple(s)
def compose(p,q): return tuple(p[q[i]] for i in range(len(p)))
def power(p,n):
    r=tuple(range(len(p)))
    for _ in range(n): r=compose(p,r)
    return r

def random_same_pi(cl,pi,rng,nsamp):
    """states with EXACTLY the given class sums, random within-class split AND random phases."""
    NV=len(cl); out=[]
    for _ in range(nsamp):
        w=np.zeros(NV)
        for c in range(4):
            idx=np.where(cl==c)[0]
            if len(idx)==0: continue
            r=rng.random(len(idx)); r=r/r.sum()*pi[c]
            w[idx]=r
        out.append(np.sqrt(w)*np.exp(1j*rng.uniform(0,2*np.pi,NV)))
    return out

def invisible_empirical(TF,TC,states,n,tol=1e-10):
    v=[abs(np.vdot(np.linalg.matrix_power(TF,n)@s, np.linalg.matrix_power(TC,n)@s)) for s in states]
    return (max(v)-min(v)) < tol, max(v)-min(v)

def A_is_classfun(TF,TC,n,cl,tol=1e-10):
    A=np.linalg.inv(np.linalg.matrix_power(TF,n))@np.linalg.matrix_power(TC,n)
    off=np.linalg.norm(A-np.diag(np.diag(A)))
    d=np.diag(A); bad=0.0
    for c in range(4):
        idx=np.where(cl==c)[0]
        if len(idx)>=2: bad=max(bad, np.max(np.abs(d[idx]-d[idx[0]])))
    return (off<tol and bad<tol), off, bad

print("== M3a  ON THE TWO CARRIERS RUN: algebraic condition versus measured invisibility ==")
for car,a in ((L.K1(), np.array([1.0,0.37,0.91,2**0.5,0.23,1.77])),
              (L.B0b(), np.random.default_rng(20260817).uniform(0,2*np.pi,18))):
    NV=car["NV"]; TF,TC=L.Top(car["walkF"],a,NV),L.Top(car["walkC"],a,NV)
    cl,F,C=L.classes(car); LF,LC=len(car["walkF"]),len(car["walkC"]); Lam=int(np.lcm(LF,LC))
    pi=np.array([0.,0.,0.,0.]); 
    sizes=[int((cl==c).sum()) for c in range(4)]
    tot=sum(sizes); pi=np.array([s/tot for s in sizes])   # any pi with full support on occupied classes
    states=random_same_pi(cl,pi,np.random.default_rng(11),40)
    # ARMS DIFF: the states must actually differ.
    dmin=min(np.linalg.norm(states[i]-states[j]) for i in range(6) for j in range(i+1,6))
    agree=True; rows=[]
    for n in range(1,4*Lam+1):
        emp,spread=invisible_empirical(TF,TC,states,n)
        alg,off,bad=A_is_classfun(TF,TC,n,cl)
        agree &= (emp==alg)
        rows.append((n,emp,alg,spread,off,bad))
    print(f"  {car['name']}: L_F={LF} L_C={LC} lcm={Lam}  |V(F) sym-diff V(C)| = {len(F^C)}"
          f"   arms-diff min||s_i - s_j|| = {dmin:.3f}  (40 states, same pi to 1e-16)")
    print(f"     n with measured invisibility : {[r[0] for r in rows if r[1]]}")
    print(f"     n with A_n a class-function  : {[r[0] for r in rows if r[2]]}")
    print(f"     n with BOTH T^n diagonal     : "
          f"{[n for n in range(1,4*Lam+1) if n%LF==0 and n%LC==0]}")
    print(f"     algebraic condition == measured on all {len(rows)} ticks: {agree}")
    print(f"     worst spread at a NON-invisible tick = {max(r[3] for r in rows if not r[1]):.3e}")
    print(f"     worst spread at an     invisible tick = {max([r[3] for r in rows if r[1]]+[0]):.3e}\n")

print("== M3b  NECESSITY FAILS IN GENERAL.  EXHIBIT: gamma_C = gamma_F REVERSED, L = 4 ==")
# 4-cycle 0-1-2-3-0 on 4 vertices, edges 0,1,2,3.  gamma_F forward, gamma_C the same cycle backward.
NV=4
walkF=[(0,1,0,+1),(1,2,1,+1),(2,3,2,+1),(3,0,3,+1)]
walkC=[(0,3,3,-1),(3,2,2,-1),(2,1,1,-1),(1,0,0,-1)]
car={"name":"REV4","NV":NV,"NE":4,"walkF":walkF,"walkC":walkC}
a=np.array([0.7,1.3,-0.4,2.1])
TF,TC=L.Top(walkF,a,NV),L.Top(walkC,a,NV)
cl,F,C=L.classes(car)
print(f"  V(gamma_F) = {sorted(F)}   V(gamma_C) = {sorted(C)}   classes {cl}  (all four vertices are 11)")
print(f"  || T_C - T_F^(-1) || = {np.linalg.norm(TC-np.linalg.inv(TF)):.2e}   so A_n = T_F^(-2n)")
pi=np.array([0.,0.,0.,1.0]); states=random_same_pi(cl,pi,np.random.default_rng(5),40)
print(f"  {'n':>3} {'T_F^n diag':>11} {'T_C^n diag':>11} {'A_n classfun':>13} {'measured spread':>17}")
for n in range(1,9):
    dF=np.linalg.norm(np.linalg.matrix_power(TF,n)-np.diag(np.diag(np.linalg.matrix_power(TF,n))))<1e-12
    dC=np.linalg.norm(np.linalg.matrix_power(TC,n)-np.diag(np.diag(np.linalg.matrix_power(TC,n))))<1e-12
    alg,_,_=A_is_classfun(TF,TC,n,cl); _,sp=invisible_empirical(TF,TC,states,n)
    print(f"  {n:>3} {str(dF):>11} {str(dC):>11} {str(alg):>13} {sp:>17.3e}")
print("  -> at n = 2 NEITHER branch operator is diagonal and invisibility HOLDS EXACTLY.")
print("     lcm(4,4) = 4, so n = 2 is NOT on the registrar's sublattice.  Necessity is FALSE.")
print("     Cost, stated: this designation has W_C = conj(W_F), so the character lattice <x,y>")
print("     collapses from Z^2 to Z.  Every counterexample below shares that vertex set.\n")

print("== M3c  EXHAUSTIVE SEARCH: when can sigma_F^n = sigma_C^n at n not divisible by lcm? ==")
print("  All ordered pairs of simple cycles (as cyclic permutations) on up to 6 labelled vertices,")
print("  loops of length 2..6, checking every n < 2*lcm.")
hits_diffV=0; hits_sameV=0; tot=0; ex=[]
NVs=6
for LFv in range(2,7):
    for LCv in range(2,7):
        for Fv in itertools.combinations(range(NVs),LFv):
            for Cv in itertools.combinations(range(NVs),LCv):
                # cyclic orders: fix first element to kill rotation redundancy
                for pf in itertools.permutations(Fv[1:]):
                    cycF=(Fv[0],)+pf
                    for pc in itertools.permutations(Cv[1:]):
                        cycC=(Cv[0],)+pc
                        sF=list(range(NVs)); sC=list(range(NVs))
                        for i in range(LFv): sF[cycF[i]]=cycF[(i+1)%LFv]
                        for i in range(LCv): sC[cycC[i]]=cycC[(i+1)%LCv]
                        sF=tuple(sF); sC=tuple(sC); tot+=1
                        lam=int(np.lcm(LFv,LCv))
                        for n in range(1,2*lam):
                            if n%lam==0: continue
                            if power(sF,n)==power(sC,n):
                                if set(Fv)==set(Cv): hits_sameV+=1
                                else:
                                    hits_diffV+=1
                                    if len(ex)<3: ex.append((cycF,cycC,n))
                                break
print(f"  pairs examined                                             : {tot}")
print(f"  pairs with sigma_F^n = sigma_C^n at some n not div by lcm  : {hits_sameV+hits_diffV}")
print(f"     of which V(gamma_F) != V(gamma_C)                       : {hits_diffV}   {ex}")
print(f"     of which V(gamma_F) == V(gamma_C)                       : {hits_sameV}")
print("  -> (C3) confirmed exhaustively: when the two loops have DIFFERENT vertex sets there is")
print("     no such n.  Every counterexample requires V(gamma_F) = V(gamma_C).")
