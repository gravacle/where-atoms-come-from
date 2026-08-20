"""ADVERSARIAL VERIFY of C-69/C-70 and model/grounded.py (2026-08-20 review).

Checks, in order:
  A. Independent tau for the T-29 grain by DIRECT TIME EVOLUTION of the population
     difference under expm(L t) -- no eigendecomposition, no Rayleigh quotient.
  B. The Law-2 label: is tau = exp(E_b/kT)/f0 what the lane actually computed?
     Compute both for all six surfaces; report the ratio 2*cosh(dE/2kT).
  C. Algebraic-identity check: are Law 1 and Law 2 identities of the construction
     (steady state (gu-gd)/(gu+gd) = tanh(dE/2kT) and record rate = gu+gd, exactly,
     BY CONSTRUCTION of detailed-balance rates)?
  D. clause_ii on a ROTATING observable (R = sx, H = -(dE/2)sz, no dissipation):
     does the shipped instrument certify as durable an observable that C-75 says
     is NOT durable?
  E. The zircon 'decline': does the model decline non-thermal physics, or only
     float64 overflow? Feed it an absurd thermal barrier just under the overflow.
  F. The CMB 'decline': what does the Law-1 steady-state extraction actually
     return if called on the bath-free surface (the lane never calls it)?
  G. Azobenzene reality check: model tau_cis vs the measured cis lifetime scale,
     under both barrier conventions.
"""
import sys, os, numpy as np

def expm(A, terms=60):
    """Scaling-and-squaring Taylor expm; adequate here (||L t|| is O(1))."""
    A = np.asarray(A, dtype=complex)
    s = max(0, int(np.ceil(np.log2(max(1.0, np.linalg.norm(A, 1))))) + 4)
    B = A / (2**s)
    E = np.eye(A.shape[0], dtype=complex); term = np.eye(A.shape[0], dtype=complex)
    for k in range(1, terms):
        term = term @ B / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'model'))
import grounded as G

eV = 1.602176634e-19
sz = np.array([[1,0],[0,-1]], dtype=complex)
sp = np.array([[0,1],[0,0]], dtype=complex)
sm = sp.conj().T

def build(dE, E_b, T, f0):
    H  = -(dE/2)*sz
    gd = f0*np.exp(-(E_b+dE/2)/(G.KB*T))
    gu = f0*np.exp(-(E_b-dE/2)/(G.KB*T))
    return H, [np.sqrt(gd)*sm, np.sqrt(gu)*sp], gd, gu

print("="*100)
print("A. INDEPENDENT tau BY DIRECT TIME EVOLUTION (T-29 grain)")
T29 = dict(T=300.0, K_u=2.0e5, V=1.26e-24, f0=1.0e9)
E_b = T29['K_u']*T29['V']; dE = 3.0*G.KB*T29['T']
H, Ls, gd, gu = build(dE, E_b, T29['T'], T29['f0'])
L = G.liouvillian(H, Ls)
rho0 = np.array([[0,0],[0,1]], dtype=complex)        # prepared 'down' (the written state)
z_ss = np.tanh(dE/(2*G.KB*T29['T']))
taus = []
tau_claim = G.lifetime(H, Ls, sz)
for t in [0.2*tau_claim, 0.5*tau_claim, 1.0*tau_claim, 2.0*tau_claim]:
    rho_t = expm(L*t) @ rho0.reshape(-1,1,order='F')
    z_t = float(np.real(np.trace(rho_t.reshape(2,2,order='F') @ sz)))
    # <sz>(t) - <sz>_ss = (z0 - z_ss) exp(-t/tau)  =>  tau = -t / ln(...)
    tau_fit = -t/np.log((z_t - z_ss)/(-1.0 - z_ss))
    taus.append(tau_fit)
    print(f"   t = {t:.3e} s   <sz>(t) = {z_t:+.6f}   tau_fit = {tau_fit:.6e} s")
tau_neel = 1.0/(gu+gd)
print(f"   lifetime() (Rayleigh quotient)  = {tau_claim:.6e} s")
print(f"   1/(gu+gd)                       = {tau_neel:.6e} s")
worst = max(abs(x-tau_neel)/tau_neel for x in taus)
print(f"   worst rel deviation of time-evolution fits from 1/(gu+gd): {worst:.2e}")
print(f"   -> {'CONFIRMED: single-exponential at rate gu+gd' if worst < 1e-8 else 'DISAGREES'}")

print()
print("="*100)
print("B. IS THE LAW-2 LABEL  tau = exp(E_b/kT)/f0  WHAT WAS COMPUTED?")
SURFACES = [
 ("CoCrPt HDD grain",  3.0*G.KB*300,  60.8*G.KB*300, 300.0, 1e9),
 ("NAND floating gate",0.30*eV,       1.60*eV,       300.0, 1e13),
 ("DNA base tautomer", 0.35*eV,       1.30*eV,       310.0, 1e13),
 ("Fe(II) crossover",  2.10*G.KB*300, 0.85*eV,       300.0, 1e12),
 ("Azobenzene",        0.60*eV,       1.05*eV,       300.0, 1e13),
 ("Alanine",           1.0e-13*eV,    1.40*eV,       300.0, 1e13),
]
print(f"   {'surface':<20}{'label formula (s)':>18}{'computed 1/(gu+gd) (s)':>24}{'ratio = 2cosh(dE/2kT)':>24}")
for name, dE_, Eb_, T_, f0_ in SURFACES:
    lab = np.exp(Eb_/(G.KB*T_))/f0_
    _,_,gd_,gu_ = build(dE_, Eb_, T_, f0_)
    comp = 1.0/(gu_+gd_)
    print(f"   {name:<20}{lab:>18.4e}{comp:>24.4e}{lab/comp:>24.4e}")

print()
print("="*100)
print("C. ARE THE TWO LAWS IDENTITIES OF THE CONSTRUCTION?")
print("   steady state of the 2-level GKSL:  p_up/p_dn = gu/gd  (rate balance), and the rates are")
print("   BUILT as gu/gd = exp(dE/kT), so <sz>_ss = (gu-gd)/(gu+gd) = tanh(dE/2kT) EXACTLY.")
print("   record mode: L_adj(sz) = (gu-gd)I - (gu+gd)sz, and <sz, I> = 0, so the Rayleigh")
print("   quotient is gu+gd EXACTLY. Symbolic check on one surface:")
H2, Ls2, gd2, gu2 = build(0.35*eV, 1.30*eV, 310.0, 1e13)
Lad = G.liouvillian(H2, Ls2).conj().T
v = sz.reshape(-1,1,order='F')/np.sqrt(2)
q = complex((v.conj().T @ Lad @ v)[0,0])
print(f"   Rayleigh quotient = {q:.6e};  -(gu+gd) = {-(gu2+gd2):.6e};  match: {abs(q+gu2+gd2)/(gu2+gd2):.1e}")
print("   -> Both 'agreements' are algebraic identities; the six mechanisms enter ONLY as the")
print("      scalars (dE, E_b, T, f0) fed to the IDENTICAL 2x2 model. No measured value is compared.")

print()
print("="*100)
print("D. clause_ii ON A ROTATING OBSERVABLE (the C-75 correction, applied to the shipped instrument)")
Hrot = -(0.5*eV)*sz; Lrot = []          # closed system, R = sx rotates at dE/hbar
c = G.clause_ii(Hrot, Lrot, np.array([[0,1],[1,0]],dtype=complex), t_m=3.156e8)
print(f"   R = sigma_x, H = -(dE/2) sigma_z, no dissipation: rotation frequency {0.5*eV/G.HBAR:.2e} rad/s")
print(f"   clause_ii returns rate = {c['rate']:.3e} /s, durable = {c['durable']}")
print(f"   -> {'HOLE: the instrument certifies a rotating (non-durable per C-75) observable as durable' if c['durable'] else 'ok'}")

print()
print("="*100)
print("E. DOES THE MODEL DECLINE NON-THERMAL PHYSICS, OR ONLY FLOAT OVERFLOW?")
for x in [500.0, 700.0, 708.0, 710.0]:
    Eb_ = x*G.KB*300
    _,_,gd_,gu_ = build(0.0, Eb_, 300.0, 1e13)
    tau = (1.0/(gu_+gd_)) if (gu_+gd_) > 0 else float('inf')
    print(f"   E_b/kT = {x:6.1f}:  tau = {tau:.3e} s   ({'returned a number' if np.isfinite(tau) else 'declined (underflow)'})")
print("   -> the 'decline' boundary is the float64 exponent range (~709), not any thermal-vs-")
print("      tunneling physics; a barrier of 708 kT (tau >> age of universe x 10^280) still returns")
print("      a number. Nothing in the model asks whether the system is thermally activated.")

print()
print("="*100)
print("F. THE CMB CONTROL: WHAT DOES THE LAW-1 EXTRACTION RETURN IF ACTUALLY CALLED?")
Lc = G.liouvillian(0.0*sz, [])
wl, Vl = np.linalg.eig(Lc)
j = int(np.argmin(np.abs(wl))); rho = Vl[:, j].reshape(2,2,order='F')
tr = np.trace(rho)
print(f"   eigenvalues all zero: {np.allclose(wl,0)};  chosen 'steady state' trace = {tr:.3f}")
if abs(tr) > 1e-12:
    rho = rho/tr
    print(f"   normalised, <sz> = {float(np.real(np.trace(rho@sz))):+.3f}  <- an ARBITRARY number")
print(f"   G.lifetime(H=0, Ls=[], R=sz) = {G.lifetime(0.0*sz, [], sz)}")
print("   -> the T-33 script never calls the extraction here; the 'model declines' framing is the")
print("      SCRIPT declining. lifetime() returns inf (defensible); the steady-state code, if run,")
print("      returns an arbitrary eigenvector, not a refusal.")

print()
print("="*100)
print("G. AZOBENZENE REALITY CHECK")
_,_,gdA,guA = build(0.60*eV, 1.05*eV, 300.0, 1e13)
tau_model = 1.0/(guA+gdA)
tau_from_cis_barrier = np.exp(1.05*eV/(G.KB*300.0))/1e13
print(f"   model tau (barrier measured from the state MIDPOINT): {tau_model:.3e} s")
print(f"   tau if E_b is the cis-side activation energy (as measured): {tau_from_cis_barrier:.3e} s")
print("   measured thermal cis->trans lifetime in the dark: ~1e4-1e6 s (hours to days)")
print(f"   -> the model's 0.4 s is off from the real record's lifetime by ~5-6 ORDERS; the midpoint")
print("      convention silently halves the escape barrier by dE/2 = 0.30 eV.")
