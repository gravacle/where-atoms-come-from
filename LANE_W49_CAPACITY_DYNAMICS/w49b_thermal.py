"""W-49b.  THE BATH HAS BEEN INFINITE-TEMPERATURE ALL ALONG. GIVE IT A TEMPERATURE.

w49's controls failed three ways and together they diagnose one thing:
  * <Nhat> = 0.50000 at EVERY mu from 0 to 8 -- the energy cost did nothing;
  * the Delta=0 control returned <Nhat> = -0.208, impossible for a projector, because the
    zero-eigenvector extraction picks an arbitrary element of a DEGENERATE steady manifold, which
    need not be a state;
  * 0.5 is exactly the maximally mixed value of the structural register.

CAUSE. Every Lindblad bath in this program -- W-28 onward -- used UNITARY jump operators (Z's,
plaquettes). A unitary jump operator has no preferred energy direction: the channel is unital and
drives everything toward the maximally mixed state. THAT IS AN INFINITE-TEMPERATURE BATH. It is
perfectly adequate for dephasing and pointer-selection questions, which is what the program has
been asking, but AN ENERGY COST CANNOT REGISTER IN IT. Growth-versus-eviction is exactly an
energetic question, so it could not have been asked with the old bath.

FIX. A Davies (thermal) bath built from a coupling operator A, with detailed balance at temperature
T:  L_ab = sqrt(f(E_a - E_b)) |a><b| A_ab , with f(w)/f(-w) = exp(-w/T), taking f(w)=1/(1+exp(w/T)).
Then the steady state is Gibbs-like and mu is meaningful.
STEADY STATE is obtained by EVOLVING a physical initial state to long time, so the result is always
a density matrix -- never by picking a null eigenvector.

CONTROLS THAT MUST NOW FIRE: <Nhat> must FALL as mu rises; must approach the Gibbs value; T -> inf
must return the old behaviour (<Nhat> -> 0.5); and Delta = 0 must freeze the structure.
"""
import itertools, numpy as np
# reuse the carrier, Move/Zop, W, ELEC, K and H_of from w49 (everything before its BATH line)
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

def davies(H,A,T,gam=0.2):
    """Davies jump operators with detailed balance at temperature T."""
    E,U=np.linalg.eigh(H); Ab=U.conj().T@A@U; Ls=[]
    for a in range(len(E)):
        for b in range(len(E)):
            if abs(Ab[a,b])<1e-12: continue
            w=E[a]-E[b]
            f=1.0/(1.0+np.exp(w/T)) if T>0 else (1.0 if w<0 else 0.0)
            if f<1e-14: continue
            Lab=np.zeros_like(H); Lab[a,b]=np.sqrt(gam*f)*Ab[a,b]
            Ls.append(U@Lab@U.conj().T)
    return Ls

def liou_thermal(mu,Delta,g2,T,gam=0.2):
    H=H_of(mu,Delta,g2); A=K(Zop([0])+Zop([1]),I2)+K(IG,sx)   # couples BOTH sectors to the bath
    Ls=davies(H,A,T,gam); Id=np.eye(D,dtype=complex)
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for Lk in Ls:
        M+=np.kron(Lk,Lk.conj())-0.5*(np.kron(Lk.conj().T@Lk,Id)+np.kron(Id,(Lk.conj().T@Lk).T))
    return M,H

def evolve_ss(M,T_end=400.0):
    rho=np.eye(D,dtype=complex)/D
    v=(expm(M*T_end)@rho.reshape(-1)).reshape(D,D)
    v=(v+v.conj().T)/2
    return v/np.trace(v).real

NHAT=K(IG,Nh)
print("W-49b  THERMAL BATH.  <Nhat> = the carrier's size.  Gibbs comparison shown.")
print(f"  {'mu':>6s} {'T':>6s} {'<Nhat>':>10s} {'Gibbs <Nhat>':>13s} {'|diff|':>9s}")
print("  "+"-"*50)
for T in (0.5,2.0,50.0):
    for mu in (0.0,0.5,1.0,2.0,4.0):
        M,H=liou_thermal(mu,0.3,0.05,T)
        r=evolve_ss(M)
        n=np.trace(NHAT@r).real
        E,U=np.linalg.eigh(H); w=np.exp(-(E-E.min())/T); w/=w.sum()
        g=float(np.sum(w*np.real(np.einsum('ai,ij,ja->a',U.conj().T,NHAT,U))))
        print(f"  {mu:6.2f} {T:6.1f} {n:10.5f} {g:13.5f} {abs(n-g):9.2e}")
    print()
print("  CONTROL -- Delta = 0 freezes the structure (no tunnelling): <Nhat> must not move with mu.")
for mu in (0.0,4.0):
    M,H=liou_thermal(mu,0.0,0.05,0.5)
    print(f"    mu={mu:4.1f} Delta=0  <Nhat> = {np.trace(NHAT@evolve_ss(M)).real:.6f}")
print()
print("  CONTROL -- T -> infinity must reproduce the OLD infinite-temperature behaviour, <Nhat>=0.5")
for T in (50.0,500.0,5000.0):
    M,H=liou_thermal(2.0,0.3,0.05,T)
    print(f"    T={T:7.0f}  mu=2.0  <Nhat> = {np.trace(NHAT@evolve_ss(M)).real:.6f}")
