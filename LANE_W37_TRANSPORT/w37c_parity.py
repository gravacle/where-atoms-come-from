"""W-37c.  WHY DO ODD RINGS CARRY EXACTLY ZERO INFORMATION?

W-37b: even rings reach ~1 bit and obey a clean perimeter law; n=5 and n=7 give EXACTLY 0.0000.
An exact zero is a symmetry, not a small effect. Find it or find the artifact.

TEST 1  Is it the initial condition? Try every start site and both a random and a localised start.
TEST 2  Are the two R-sectors' Hamiltonians actually isospectral / unitarily equivalent on the
        probe? If some unitary V acts only on the gauge factor with V H V^dag = H and V R V^dag = -R,
        the probe can never tell the sectors apart and the zero is exact and structural.
        Detect it directly: compare the probe-position dynamics generated in each sector.
"""
import itertools, numpy as np
exec(open('w37b_scaling.py').read().split('print("W-37b')[0])

def sector_dynamics(n, tau=1.0, tmax=12.0, dt=0.25):
    ST,IDX,D,R,H,_=ring(n,tau)
    Pp=(np.eye(D)+R)/2; Pm=(np.eye(D)-R)/2
    out={}
    for nm,Pr in (("R=+1",Pp),("R=-1",Pm)):
        # probe starts at site 1, projected into this sector
        v=np.zeros(D,complex)
        for i,(p,s) in enumerate(ST):
            if p==1: v[i]=1.0
        v=Pr@v
        if np.linalg.norm(v)<1e-12: out[nm]=None; continue
        v/=np.linalg.norm(v)
        U=expm(-1j*H*dt); traj=[]; psi=v.copy(); t=0.0
        while t<tmax:
            psi=U@psi; t+=dt
            pr=np.zeros(n)
            for i,(p,s) in enumerate(ST): pr[p]+=abs(psi[i])**2
            traj.append(pr)
        out[nm]=np.array(traj)
    return out,ST,D,R,H

print("W-37c  ODD-RING NULL: symmetry or artifact?")
print()
print("  TEST 1 -- does ANY start site give information on an odd ring?")
for n in [5,7]:
    best_over_starts=0.0
    for start in range(n):
        ST,IDX,D,R,H,_=ring(n)
        Pp=(np.eye(D)+R)/2; Pm=(np.eye(D)-R)/2
        g=np.random.default_rng(3)
        mask=np.array([1.0 if p==start else 0.0 for p,_ in ST])
        w=(g.normal(size=D)+1j*g.normal(size=D))*mask
        a=Pp@w; b=Pm@w
        if min(np.linalg.norm(a),np.linalg.norm(b))<1e-12: continue
        a/=np.linalg.norm(a); b/=np.linalg.norm(b); psi=(a+b); psi/=np.linalg.norm(psi)
        U=expm(-1j*H*0.25); t=0.0
        while t<20.0:
            psi=U@psi; t+=0.25
            br=[]
            for Pr in (Pp,Pm):
                v=Pr@psi; pr=float(np.vdot(v,v).real)
                if pr<1e-14: br.append((0.0,None)); continue
                v/=np.sqrt(pr)
                M=np.zeros((n,n),complex)
                for i,(p,s) in enumerate(ST):
                    for jj,(q,t2) in enumerate(ST):
                        if s==t2: M[p,q]+=v[i]*np.conj(v[jj])
                br.append((pr,M))
            avg=sum(p*m for p,m in br if m is not None)
            best_over_starts=max(best_over_starts, vn(avg)-sum(p*vn(m) for p,m in br if m is not None))
    print(f"    n={n}: best I(R:probe) over ALL start sites and times = {best_over_starts:.3e}")

print()
print("  TEST 2 -- are the two sectors' probe dynamics identical? (an exact match = a symmetry)")
print(f"    {'n':>4s} {'max |P_+(t) - P_-(t)| over sites and times':>44s}")
print("    "+"-"*50)
for n in [4,5,6,7,8,9,10]:
    d,ST,D,R,H=sector_dynamics(n)
    if d["R=+1"] is None or d["R=-1"] is None:
        print(f"    {n:4d}   (a sector is empty for a localised start)"); continue
    diff=np.abs(d["R=+1"]-d["R=-1"]).max()
    print(f"    {n:4d} {diff:44.3e}   {'IDENTICAL -> probe cannot distinguish the sectors' if diff<1e-10 else 'different -> probe can read the flux'}")

print()
print("  TEST 3 -- the spectra. For a ring, R=+1 is periodic and R=-1 antiperiodic hopping.")
print("  Odd n makes the two spectra coincide as SETS, so no probe observable can separate them.")
for n in [4,5,6,7,8,9]:
    per =np.sort(np.round(-2*np.cos(2*np.pi*np.arange(n)/n),9))
    anti=np.sort(np.round(-2*np.cos((2*np.arange(n)+1)*np.pi/n),9))
    same=np.allclose(per,anti)
    print(f"    n={n}: periodic {per}  antiperiodic {anti}  -> {'SAME SET' if same else 'different'}")
