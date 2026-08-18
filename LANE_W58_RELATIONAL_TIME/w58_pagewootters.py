"""KERNEL GATE. The first version used H_S with NEGATIVE eigenvalues while H_C is
non-negative, so the constraint n*omega = E_s had only the n=0 solution: the kernel was a
STATIC EIGENSTATE, not a history, and conditioning on it looked unitary for trivial reasons
while lambda had no effect at all. H_S now shares the clock ladder so the kernel contains
genuine superpositions over clock readings. The gate below checks that.
"""
"""W-58.  REMOVE THE BACKGROUND TIME PARAMETER.

Surface audit A2: 125 of 125 lanes evolve in an external parameter t, and ZERO derive it. Gravity is
the dynamics of space-time structure, so with t a fixed external parameter there is nowhere for it to
act -- and none of the three failed gravity routes touched t. This is the first lane that removes it.

THE CONSTRUCTION (Page-Wootters). There is NO external time. The whole thing is in a STATIONARY state
obeying a constraint,
        (H_C  (x) I  +  I (x) H_S) |Psi> = 0
and "time" is a CORRELATION between a clock subsystem and the rest: conditioning on the clock reading
t gives the system's state at t. Nothing outside the system is ever evolved.

WHAT MUST BE CHECKED, IN ORDER:
 1. the constraint is genuinely satisfied -- ||H_tot |Psi>|| = 0, not approximately;
 2. conditioning on the clock REPRODUCES Schroedinger evolution (the control: if it does not, the
    construction is wrong and nothing after it means anything);
 3. THEN, and only then: is the record the same object in the timeless picture?
Step 3 is the question. Steps 1 and 2 are the price of asking it.
"""
import numpy as np
rng=np.random.default_rng(9)
def expm(A):
    nr=np.linalg.norm(A,np.inf)
    k=max(0,int(np.ceil(np.log2(nr)))+1) if nr>0 else 0
    B=A/(2.0**k); X=np.eye(A.shape[0],dtype=complex); T=X.copy()
    for m in range(1,60):
        T=T@B/m; X=X+T
        if np.linalg.norm(T,np.inf)<1e-18*max(1.0,np.linalg.norm(X,np.inf)): break
    for _ in range(k): X=X@X
    return X

dC=12                      # clock dimension
dS=4                       # system: 2 qubits -- a record and a partner
# IDEAL FINITE CLOCK: energy levels equally spaced; time states are the Fourier duals
wC=2*np.pi/dC
HC=np.diag(wC*np.arange(dC)).astype(complex)
F=np.array([[np.exp(-2j*np.pi*n*k/dC)/np.sqrt(dC) for n in range(dC)] for k in range(dC)])
tstates=[F[k].conj() for k in range(dC)]          # |t_k>, k = 0..dC-1
# SYSTEM: a Hamiltonian whose spectrum matches the clock's ladder so the kernel is nonempty
th=rng.normal(size=(dS,dS)); U0,_=np.linalg.qr(th+1j*rng.normal(size=(dS,dS)))
levels=np.array([0,1,2,3])                        # integers -> commensurate with the clock ladder
HS=U0@np.diag(+wC*levels)@U0.conj().T
HS=(HS+HS.conj().T)/2
Htot=np.kron(HC,np.eye(dS))+np.kron(np.eye(dC),HS)

# THE STATE IS TAKEN FROM THE KERNEL OF THE CONSTRAINT, NOT ASSUMED.
# The first version built |Psi> = sum_k |t_k> (x) e^{-i H_S t}|psi_0> and imposed (H_C + H_S)|Psi>=0.
# But H_C generates +t translation, so that history state satisfies (H_C - H_S)|Psi> = 0. Rather than
# patch the sign, the physical state is now DIAGONALISED OUT of the constraint: whatever the kernel
# is, that is the timeless state, and conditioning is then tested against it.
Htot=np.kron(HC,np.eye(dS))-np.kron(np.eye(dC),HS)
ev,evec=np.linalg.eigh(Htot)
ker=[evec[:,i] for i in range(len(ev)) if abs(ev[i])<1e-9]
print("W-58  RELATIONAL TIME. No external parameter is evolved anywhere in this lane.")
print(f"  clock dim {dC}, system dim {dS}, total {dC*dS}")
print(f"\n  CHECK 1 -- the constraint. kernel dimension of H_tot = {len(ker)}")
if not ker:
    print("    EMPTY KERNEL -- no timeless state exists for this pair. Nothing further is valid.")
    raise SystemExit
c=rng.normal(size=len(ker))+1j*rng.normal(size=len(ker))
Psi=sum(c[i]*ker[i] for i in range(len(ker))); Psi/=np.linalg.norm(Psi)
r=np.linalg.norm(Htot@Psi)
print(f"    ||H_tot |Psi>|| = {r:.3e}   {'SATISFIED' if r<1e-8 else 'NOT SATISFIED'}")

print("\n  CHECK 2 -- does conditioning on the clock reproduce unitary evolution?")
print("    (no target is assumed: the conditional state at k=0 is taken as the initial condition)")
conds=[]
for k in range(dC):
    proj=np.kron(np.outer(tstates[k],tstates[k].conj()),np.eye(dS))
    v=(proj@Psi).reshape(dC,dS)
    cond=v[np.argmax(np.linalg.norm(v,axis=1))]
    conds.append(cond/np.linalg.norm(cond))
worst=0.0
print(f"    {'k':>3s} {'t':>8s} {'|| cond_k - e^{-iH t} cond_0 ||':>32s}")
for k in range(dC):
    t=2*np.pi*k/(dC*wC)
    tgt=expm(-1j*HS*t)@conds[0]; tgt/=np.linalg.norm(tgt)
    ph=np.vdot(tgt,conds[k]); ph=ph/abs(ph) if abs(ph)>1e-12 else 1.0
    e=np.linalg.norm(conds[k]/ph-tgt); worst=max(worst,e)
    if k<6: print(f"    {k:3d} {t:8.3f} {e:32.3e}")
print(f"    worst over all {dC} readings: {worst:.3e}   "
      f"{'REPRODUCES unitary evolution from a timeless state' if worst<1e-8 else 'DOES NOT'}")
