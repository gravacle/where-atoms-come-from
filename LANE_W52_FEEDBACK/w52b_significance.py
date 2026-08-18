"""W-52b.  IS THE SPATIAL EFFECT SIGNIFICANT, AND WHAT TURNS IT ON?

w52 (after fixing the bath and the initialisation) found that loading a record on p0 helps formation
on its NEIGHBOUR p1 more than on the NON-NEIGHBOUR p3 -- but ONLY when g^2 is nonzero:
    g2=0.00   d(neighbour) 0.037   d(far) 0.040   difference -0.003   no spatial structure
    g2=0.05   d(neighbour) 0.089   d(far) 0.049   difference +0.040
Before that is read as anything, it needs error bars: E[<R>^2] is a Monte Carlo average over
trajectories and the differences are small compared with the values.

THIS RUN: independent seeds, standard errors on every quantity, and a g^2 sweep to locate where the
spatial structure switches on. A MECHANISM IS ALSO STATED IN ADVANCE, so it can be checked rather
than invented afterwards: at g^2 = 0 the Hamiltonian is purely magnetic and the plaquette records are
INDEPENDENT degrees of freedom -- nothing can couple them, so no spatial structure is possible. The
electric term acts on LINKS, and neighbouring plaquettes SHARE a link. If that is the mechanism, the
spatial difference must vanish at g^2 = 0 and grow with g^2, and it must track SHARED LINKS rather
than anything else.
"""
import itertools, numpy as np
exec(open('w52_feedback.py').read().split('print()\nprint("  GATE')[0])

def trial(load,g2,gam,seed,NT=500,T=15.0):
    return run(load,g2,gam,T=T,NT=NT,seed=seed)

def stat(load,g2,gam,seeds):
    a=[];b=[]
    for s in seeds:
        o=trial(load,g2,gam,s)
        a.append(o[NEI][0]); b.append(o[FAR][0])
    a=np.array(a); b=np.array(b)
    return a.mean(),a.std(ddof=1)/np.sqrt(len(a)),b.mean(),b.std(ddof=1)/np.sqrt(len(b))

SEEDS=[11,23,37,51,67,83]
print("W-52b  SIGNIFICANCE.  6 independent seeds x 500 trajectories each.")
print(f"  bath = {BATH} (avoids p0, p1 and p3)   NEIGHBOUR=p{NEI} shares 1 link with p0, "
      f"FAR=p{FAR} shares {len(sup[0]&sup[FAR])}")
print()
print(f"  {'g^2':>6s} {'gamma':>6s} {'d(neigh)':>18s} {'d(far)':>18s} {'difference':>18s} {'sigma':>6s}")
print("  "+"-"*80)
for g2 in (0.0,0.02,0.05,0.10,0.20):
    for gam in (0.5,):
        n0,e_n0,f0,e_f0 = stat(False,g2,gam,SEEDS)
        n1,e_n1,f1,e_f1 = stat(True ,g2,gam,SEEDS)
        dn=n1-n0; edn=np.hypot(e_n1,e_n0)
        df=f1-f0; edf=np.hypot(e_f1,e_f0)
        d=dn-df;  ed=np.hypot(edn,edf)
        print(f"  {g2:6.2f} {gam:6.2f} {dn:11.5f} +-{edn:6.5f} {df:11.5f} +-{edf:6.5f} "
              f"{d:11.5f} +-{ed:6.5f} {abs(d)/max(ed,1e-12):6.1f}")
print()
print("  MECHANISM CHECK -- stated before running: the effect requires SHARED LINKS and the electric")
print("  term. At g^2 = 0 the Hamiltonian is purely magnetic and the plaquette records are")
print("  independent, so no spatial structure is possible.")
H0=-MAG
print(f"    ||[W_p0, W_p1]|| = {np.linalg.norm(W[0]@W[1]-W[1]@W[0]):.3e}   (plaquettes commute)")
print(f"    ||[H(g2=0), W_p1]|| = {np.linalg.norm(H0@W[1]-W[1]@H0):.3e}")
for g2 in (0.05,):
    Hg=-MAG-g2*ELEC
    print(f"    ||[H(g2={g2}), W_p1]|| = {np.linalg.norm(Hg@W[1]-W[1]@Hg):.3e}   "
          f"<- the electric term is the only channel coupling plaquettes")
