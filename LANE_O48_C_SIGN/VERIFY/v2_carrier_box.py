"""
VERIFY-2.  IS THE 'CARRIER' BOX A MEASUREMENT OF ENERGY CONTENT, OR OF DETERMINISM?

The lane's carrier axis is  all_const & varies : the observable must be a DETERMINISTIC FUNCTION
of the energy.  That is much stronger than "carries energy".  This script re-runs the same sweep
with the honest energetic axis the lane itself adopted for Exhibit 5 -- exact mutual information
I(O;E) > 0 -- and compares the population of the boxes.

It also tests, separately, whether clause (iv) has any power on the INFORMATIONAL axis:
   does any clause-(iv)-balanced observable have I(O;E) > 0 ?
"""
import numpy as np, random

def hadamard(n):
    Hm = np.ones((1,1), dtype=np.int8)
    for _ in range(n): Hm = np.block([[Hm,Hm],[Hm,-Hm]])
    return Hm
def spin_table(n):
    xs = np.arange(2**n)
    return np.array([1-2*((xs>>i)&1) for i in range(n)], dtype=np.int8).T
def couplings(name,m,seed=0):
    rnd=random.Random(seed)
    if name=="uniform":  return [1]*m
    if name=="randpos":  return [rnd.randrange(1,61) for _ in range(m)]
    if name=="superinc": return [2**i for i in range(m)]
    if name=="randsign": return [rnd.choice((1,-1))*rnd.randrange(1,61) for _ in range(m)]

print("="*104)
print("VERIFY-2   THE LANE'S 'ENERGY' BOX COUNTS OBSERVABLES DETERMINED BY E, NOT OBSERVABLES")
print("           THAT CARRY ENERGY.  Exact mutual information, same sweep, same observables.")
print("="*104)
print()
print(f"  {'family':<10} {'n':>3} {'2^n':>6} {'lane: energy-only':>18} {'#obs with I(O;E)>0':>20} "
      f"{'#records w/ I>0':>16} {'max I over records':>19} {'max I over I>0 set':>19}")
for fam in ("superinc","randpos","uniform","randsign"):
    for n in (7,9,11,13):
        J=couplings(fam,n-1); S=spin_table(n)
        E=np.zeros(2**n,dtype=np.int64)
        for i in range(n-1): E+=J[i]*S[:,i].astype(np.int64)*S[:,i+1].astype(np.int64)
        order=np.argsort(E,kind="stable"); Es=E[order]
        starts=np.flatnonzero(np.r_[True,Es[1:]!=Es[:-1]])
        Hm=hadamard(n)[:,order].astype(np.int32)
        bsum=np.add.reduceat(Hm,starts,axis=1)
        bmax=np.maximum.reduceat(Hm,starts,axis=1); bmin=np.minimum.reduceat(Hm,starts,axis=1)
        const=(bmax==bmin); balanced=np.all(bsum==0,axis=1); nonconst=np.any(~const,axis=1)
        is_rec=balanced&nonconst
        carries=np.all(const,axis=1)&(bmax.max(axis=1)!=bmax.min(axis=1))
        sizes=np.add.reduceat(np.ones(2**n,dtype=np.int64),starts); N=2**n
        npl=(sizes[None,:]+bsum)//2; nmi=(sizes[None,:]-bsum)//2
        tp=npl.sum(axis=1).astype(float); tm=nmi.sum(axis=1).astype(float)
        def term(cnt,tot):
            p=cnt/N; pb=sizes[None,:]/N; po=(tot/N)[:,None]
            with np.errstate(divide="ignore",invalid="ignore"):
                t=p*np.log2(p/(pb*po))
            return np.where(cnt>0,t,0.0).sum(axis=1)
        MI=np.maximum(term(npl,tp)+term(nmi,tm),0.0)
        informative=MI>1e-12
        mi_rec = MI[is_rec].max() if is_rec.any() else float('nan')
        print(f"  {fam:<10} {n:>3} {2**n:>6} {int((carries&~is_rec).sum()):>18} "
              f"{int(informative.sum()):>20} {int((is_rec&informative).sum()):>16} "
              f"{mi_rec:>19.6f} {MI[informative].max():>19.6f}")
    print()
print("  READ: for the degenerate families the lane's 'energy-only' box holds 1 observable out of")
print("  8192 while thousands of observables have strictly positive I(O;E).  The box is a")
print("  determinism test, not an energy-content test, and 'the other three boxes are populated'")
print("  is a statement about that determinism test only.")
print()
print("  THE ONE THING THAT IS TRUE, AND IT IS TWO LINES OF ALGEBRA, NOT A SWEEP:")
print("  clause (iv) says Tr(P_E O)=0 on every eigenspace.  For a +-1-valued O that is exactly")
print("  #(O=+1 | E) = #(O=-1 | E) in EVERY block, i.e. P(O=+1|E)=1/2 for every E, i.e. I(O;E)=0.")
print("  The column '#records with I>0' above is 0 for that reason and for no other; it is the")
print("  framework's own stated equivalence for clause (iv) restated, true at every n with no")
print("  computation, and it needs neither 67 million pairs nor an Ising chain.")
