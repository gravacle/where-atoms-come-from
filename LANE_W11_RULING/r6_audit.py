# RULING LEG 6 -- I READ THE FIVE CROSS-REFUTATIONS AS SCEPTICALLY AS THE LANES.
# Only the charges that bear on MY ruling are re-run.  Default: refute.
import numpy as np, rlib
rng = np.random.default_rng(rlib.SEED)
K = rlib.K1()

print("== R6.1  THE STEELMAN-CROSS'S 'EXHAUSTIVE OVER ALL 729 PAIRS' IS NOT EXHAUSTIVE ==")
print("   It claims the gauge-covariant root variety on a 3-cycle is 27 branch choices on T's")
print("   eigenvectors, and enumerates 729 pairs.  Lam.T is covariant, continuous, a root, and")
print("   is a 2-TORUS per loop.  Is it a function of T (i.e. one of the 27)?  [T,U] = 0 iff yes.")
a = rlib.a_generic(K, rng, 1.0, 2**0.5); T = rlib.Tedge(K, K.walkF, a)
ncomm = 0; wcov = 0.0; wroot = 0.0
for _ in range(200):
    ph = rng.uniform(0, 2*np.pi, 3); ph[-1] = -ph[:-1].sum()
    Lam = np.eye(5, dtype=complex)
    for k, v in enumerate(sorted(K.VF)): Lam[v, v] = np.exp(1j*ph[k])
    U = Lam@T
    wroot = max(wroot, np.linalg.norm(np.linalg.matrix_power(U,3)-rlib.Mcirc(K,K.VF,rlib.holon(K.walkF,a))))
    th = rng.uniform(0, 2*np.pi, 5)
    ag = np.array([a[j]+th[t]-th[s] for j,(s,t) in enumerate(K.edges)])
    Ug = Lam@rlib.Tedge(K, K.walkF, ag); g = np.exp(1j*th); s = rng.normal(size=5)+1j*rng.normal(size=5)
    wcov = max(wcov, np.linalg.norm(Ug@(g*s)-g*(U@s)))
    ncomm += (np.linalg.norm(U@T-T@U) < 1e-9)
print(f"   200 draws of Lam.T:  max||U^3 - M|| = {wroot:.2e}   max covariance defect = {wcov:.2e}")
print(f"   commutes with T (=> a function of T, i.e. one of the 27): {ncomm}/200")
print("   -> the covariant root variety is POSITIVE-DIMENSIONAL, not 27 points.  The steelman-")
print("      cross's leg E3 is not exhaustive and its 'EXACTLY ONE passes' is not established.")
print("      Its VERDICT is unharmed -- every member of Lam.T is pi-visible (my R4.3) -- but the")
print("      one leg it offers as EXACT and EXHAUSTIVE is neither.\n")

print("== R6.2  'THE COMPARISON TIME' IS NOT THE OPERATIVE VARIABLE.  Hold n = 1, move the operator ==")
for C in (rlib.K1(), rlib.B0b()):
    a = rlib.a_generic(C, rng, 1.0, 2**0.5)
    WF, WC = rlib.holon(C.walkF,a), rlib.holon(C.walkC,a)
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 24)
    TF,TC = rlib.Tedge(C,C.walkF,a), rlib.Tedge(C,C.walkC,a)
    DF,DC = rlib.Droot(C,C.VF,WF,C.LF), rlib.Droot(C,C.VC,WC,C.LC)
    print(f"  {C.name:4s} TICK HELD AT n = 1.  T: spread {rlib.pi_spread(C,TF,TC,S,[1]):.2e}"
          f"    D: spread {rlib.pi_spread(C,DF,DC,S,[1]):.2e}")
print("  -> one variable moved, the tick fixed, and the answer flips.  The time is ONE coordinate")
print("     of the operative variable, not the whole of it.  (Symmetrically, R1.4 fixed the")
print("     operator at T and moved only the tick, and the answer flips there too.)\n")

print("== R6.3  'FIBRE-WISE AND LOOP-CONSTANT' IS SUFFICIENT, NOT NECESSARY, AND TOO STRONG ==")
print("   (a) CLASS-constant beats LOOP-constant: give gamma_F's two classes DIFFERENT phases.")
for C in (rlib.K1(), rlib.B0b()):
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 24)
    ph = rng.uniform(0, 2*np.pi, 4)
    AF = np.diag([np.exp(1j*ph[C.CLASSES.index(C.cls[v])]) for v in range(C.nv)])
    AC = np.diag([np.exp(1j*rng.uniform(0,2*np.pi)*(C.cls[v][1] or 0.7)) for v in range(C.nv)])
    loopconst = len({round(float(np.angle(AF[v,v])),9) for v in C.VF}) == 1
    print(f"    {C.name:4s} A_F class-constant, LOOP-constant? {loopconst}"
          f"   pi-spread over n<=12 = {rlib.pi_spread(C,AF,AC,S,range(1,13)):.2e}")
print("   (b) NON-fibre-wise pairs that are pi-blind (the CORRELATED locus), U_F = L_F R, U_C = L_C R:")
def haar(n,rng):
    z=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))/np.sqrt(2)
    q,r=np.linalg.qr(z); return q*(np.diag(r)/np.abs(np.diag(r)))
for C in (rlib.K1(), rlib.B0b()):
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 24)
    worst = 0.0; nondiag = 0
    for _ in range(200):
        R = np.zeros((C.nv,C.nv),dtype=complex)
        for c in C.CLASSES:
            vs = C.idx[c]
            if vs: R[np.ix_(vs,vs)] = haar(len(vs),rng)
        dF, dC = rng.uniform(0,2*np.pi,4), rng.uniform(0,2*np.pi,4)
        LF = np.diag([np.exp(1j*dF[C.CLASSES.index(C.cls[v])]) for v in range(C.nv)])
        LC = np.diag([np.exp(1j*dC[C.CLASSES.index(C.cls[v])]) for v in range(C.nv)])
        UF, UC = LF@R, LC@R
        nondiag += (not np.allclose(UF, np.diag(np.diag(UF))))
        worst = max(worst, rlib.pi_spread(C,UF,UC,S,range(1,13)))
    print(f"    {C.name:4s} 200 draws: U_F NON-diagonal in {nondiag}/200;  worst pi-spread n<=12 = {worst:.2e}")
print("   -> the criterion is on the RELATIVE operator Q_n = (U_F^n)^* U_C^n, not on each branch")
print("      operator.  'fibre-wise and loop-constant' names a SUFFICIENT condition as if it were")
print("      the whole content.  Lane M's C1 and lane I's D2 have this right; three lanes do not.\n")

print("== R6.4  THE TRIVIAL-CONNECTION CONTACT POINT.  Measured, and scored NEITHER WAY ==")
print("   FOUNDING_DESIGN :117-118 / S2 :583 pre-register 'no formation at the trivial connection'.")
for C in (rlib.K1(), rlib.B0b()):
    a0 = np.zeros(len(C.edges))
    MF,MC = rlib.Mcirc(C,C.VF,1.0), rlib.Mcirc(C,C.VC,1.0)
    TF,TC = rlib.Tedge(C,C.walkF,a0), rlib.Tedge(C,C.walkC,a0)
    base = rng.dirichlet(np.ones(C.nv)); s = np.sqrt(base)*np.exp(1j*rng.uniform(0,2*np.pi,C.nv))
    zM = [abs(rlib.Z(MF,MC,s,n)) for n in range(1,13)]
    zT = [abs(rlib.Z(TF,TC,s,n)) for n in range(1,13)]
    e = np.zeros(C.nv,dtype=complex); e[sorted(C.VF)[0]] = 1.0
    print(f"  {C.name:4s} a = 0:  min_n |Z^M_n| = {min(zM):.12f}   min_n |Z^T_n| = {min(zT):.6f}")
    print(f"       and at a = 0, T is a PERMUTATION: |<e0, T_F e0>| = {abs(np.vdot(e,TF@e)):.3f}")
print("  -> the contact point EXCLUDES every tick that moves amplitude and ADMITS exactly the")
print("     fibre-wise class.  It therefore does discriminate -- but its admitted set IS the")
print("     stipulated category, and it also excludes the corpus's OWN ledgered alternative")
print("     (S2 audit A2's 'a real parameter t with a Hamiltonian').  READS TWO WAYS.  SEE RULING.")
