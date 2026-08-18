"""W-49c.  DOES THE CARRIER'S SIZE DEPEND ON WHAT IS RECORDED IN IT?

w49b gave the bath a temperature and the structural register thermalised to Gibbs to 1e-9. But that
raised the real tension:
    a THERMAL bath relaxes everything to Gibbs, INCLUDING the records.
    our records survived precisely because the old bath COMMUTED with them.
So "how many records are loaded" is not even a property of a thermal steady state -- it forgets.

RESOLVE IT THE ONLY HONEST WAY: couple the bath to something that COMMUTES with the records. Then
the records are exactly conserved (they persist), while the structural register still thermalises.
Here the bath couples to sigma_x on the STRUCTURE ALONE:  A = I_gauge (x) sigma_x. It cannot touch
the gauge sector at all, so every plaquette flux is conserved.

THEN THE QUESTION IS SHARP AND NOT CIRCULAR:
   within each conserved record sector, what size does the carrier settle to?
   if <Nhat> DIFFERS between record sectors, the content determines the geometry.
   if it is the same in all of them, the carrier does not care what it holds.

NOTHING IN H REFERENCES A RECORD, A BOUNDARY OR A CAPACITY. The only structure-content coupling is
the gating -Nhat (x) (W2 + h.c.), i.e. "the third plaquette exists only if the carrier is large",
which is what it MEANS for the carrier to have that plaquette. The correlation, if any, is a
consequence.

CONTROL: gate strength -> 0 must make every sector identical.
"""
import itertools, numpy as np
exec(open('w49_dynamics.py').read().split('BATH=[0]')[0])
def expm(A):
    nr=np.linalg.norm(A,np.inf)
    k=max(0,int(np.ceil(np.log2(nr)))+1) if nr>0 else 0
    B=A/(2.0**k); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for m in range(1,60):
        T=T@B/m; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(k): X=X@X
    return X
I2=np.eye(2,dtype=complex); sx=np.array([[0,1],[1,0]],complex); Nh=np.array([[0,0],[0,1]],complex)
def K(a,b): return np.kron(a,b)
def Hg(mu,Delta,g2,gate=1.0):
    H =K(-(W[0]+W[0].conj().T),I2)+K(-(W[1]+W[1].conj().T),I2)
    H+=gate*K(-(W[2]+W[2].conj().T),Nh)
    H+=K(-g2*ELEC,I2)+K(IG,mu*Nh)-Delta*K(IG,sx)
    return (H+H.conj().T)/2
def davies(H,A,T,gam=0.2):
    E,U=np.linalg.eigh(H); Ab=U.conj().T@A@U; Ls=[]
    for a in range(len(E)):
        for b in range(len(E)):
            if abs(Ab[a,b])<1e-12: continue
            w=E[a]-E[b]; f=1.0/(1.0+np.exp(w/T))
            if f<1e-14: continue
            M=np.zeros_like(H); M[a,b]=np.sqrt(gam*f)*Ab[a,b]; Ls.append(U@M@U.conj().T)
    return Ls
def liou(mu,Delta,g2,T,gate=1.0):
    H=Hg(mu,Delta,g2,gate); A=K(IG,sx)          # bath touches the STRUCTURE ONLY
    Ls=davies(H,A,T); Id=np.eye(D,dtype=complex)
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for Lk in Ls: M+=np.kron(Lk,Lk.conj())-0.5*(np.kron(Lk.conj().T@Lk,Id)+np.kron(Id,(Lk.conj().T@Lk).T))
    return M,H
NHAT=K(IG,Nh)
R2=K((W[2]+W[2].conj().T)/2.0,I2)                 # the record on the gated plaquette
ev=np.unique(np.round(np.linalg.eigvals(R2).real,6))
print("W-49c  DOES THE CARRIER'S SIZE DEPEND ON WHAT IT HOLDS?")
print(f"  bath couples to structure only; ||[A, plaquette ops]|| = "
      f"{max(np.linalg.norm(K(IG,sx)@K(W[i],I2)-K(W[i],I2)@K(IG,sx)) for i in range(3)):.1e}"
      f"  -> records exactly conserved")
print(f"  record sectors of R2: {ev}")
print()
R2c=lambda g2: np.linalg.norm(R2@Hg(1.0,0.3,g2)-Hg(1.0,0.3,g2)@R2)
print("  GATE: the record must be CONSERVED or the sectors mix and the test is void.")
for g2 in (0.05,0.0):
    print(f"    ||[R2, H]|| at g2={g2:4.2f} : {R2c(g2):.3e}"
          + ("   <- g2*ELEC does not commute with the plaquettes; sectors LEAK" if R2c(g2)>1e-9
             else "   <- record exactly conserved"))
print()
print(f"  {'mu':>5s} {'T':>5s} {'g2':>5s} {'gate':>5s} " + " ".join(f"<N|R2={v:+.0f}>" for v in ev) + "   difference")
print("  "+"-"*66)
for gate in (1.0,0.5,0.0):
    for mu in (1.0,2.0):
        for g2 in (0.05,0.0):
            T=0.5
            M,H=liou(mu,0.3,g2,T,gate)
            outs=[]
            for v in ev:
                P=np.eye(D,dtype=complex)
                for u in ev:
                    if abs(u-v)>1e-9: P=P@((R2-u*np.eye(D))/(v-u))
                r0=P@P.conj().T; r0=r0/np.trace(r0).real
                rr=(expm(M*400.0)@r0.reshape(-1)).reshape(D,D)
                rr=(rr+rr.conj().T)/2; rr=rr/np.trace(rr).real
                outs.append(np.trace(NHAT@rr).real)
            print(f"  {mu:5.2f} {T:5.1f} {g2:5.2f} {gate:5.2f} " + " ".join(f"    {o:9.6f}" for o in outs)
                  + f"   {max(outs)-min(outs):.3e}")
print()
print("  READING: a nonzero difference means the carrier settles to a DIFFERENT SIZE depending on")
print("  what is recorded in it. gate=0 removes the structure-content coupling and must give 0.")
