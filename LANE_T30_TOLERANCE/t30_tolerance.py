"""T-30: DOES THE GRAVITY-STRENGTH EXCLUSION SURVIVE AT FINITE TOLERANCE?

C-52 -- the Z[i] quantisation theorem -- is the load-bearing exclusion behind two days of negative
results: a product of record operators on an abelian stabiliser carrier is a single SIGNED PAULI, so
every Gaussian-integer polynomial in record operators has trace ratio in Z[i], which has no element
of modulus strictly between 0 and 1. Hence "a 1e-36 residual has no slot".

BUT IT WAS PROVED AT TOLERANCE ZERO, over EXACT signed Paulis in an EXACT code space -- DEF-A's
corner. A real record is a SLOW MODE with a finite lifetime, not an exact Pauli. The censuses found
ZERO real records in that corner. Exclusions proved at tolerance zero do not automatically survive at
finite tolerance, and the census named this the computation that decides whether the exclusion was a
theorem about the WORLD or about the corner.

THE MEASUREMENT. Take the record of a metastable carrier -- the CoCrPt grain of T-29, and a
tunable family interpolating from exact to metastable -- and compute the same trace ratios. Ask:
are they QUANTISED in Z[i], or CONTINUOUS? And if continuous, how small can they get?

CONTROL IN THE SAME TABLE (D-15): the exact stabiliser carrier, where C-52 says the ratios must land
exactly on Z[i] and the minimum non-zero modulus must be exactly 1. If the instrument does not
reproduce that, it cannot be trusted on the metastable side."""
import sys, os, itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import grounded as G
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def word(n,s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
say("="*100); say("T-30   IS THE Z[i] QUANTISATION AN ARTIFACT OF EXACTNESS?"); say("="*100)

# ---------------------------------------------------------------- CONTROL: the exact carrier
say("")
say("1. CONTROL — THE EXACT STABILISER CARRIER. C-52 says every ratio lands on Z[i] and the")
say("   minimum non-zero modulus is exactly 1. If this fails, the instrument is broken.")
n=4
Hs=-(word(n,'X'*n)+word(n,'Z'*n))
w,V=np.linalg.eigh(Hs); ix=[i for i in range(len(w)) if abs(w[i]-w.min())<1e-9]
Pg=V[:,ix]@V[:,ix].conj().T; k=len(ix)
recs=[word(n,'ZZII'), word(n,'IZZI'), word(n,'XXII')]
# ERRATUM: a first version summed only NON-IDENTITY logicals. Every such term has vanishing
# normalised trace on the code space, so all 4000 ratios were exactly 0 and the instrument measured
# nothing. C-52's non-zero values come from terms landing ON the stabiliser group -- squares and the
# identity -- so the polynomial family must contain them.
rng=np.random.default_rng(3); mods=[]
for _ in range(4000):
    M=int(rng.integers(-2,3))*np.eye(2**n,dtype=complex)          # the identity term C-52 needs
    for R in recs:
        M = M + int(rng.integers(-2,3))*R + int(rng.integers(-1,2))*(R@R)
    for a,b in itertools.combinations(recs,2): M = M + int(rng.integers(-2,3))*(a@b)
    r=complex(np.trace(Pg@M)/k)
    if abs(r)>1e-12: mods.append(abs(r))
say(f"   4000 integer-coefficient polynomials: {len(mods)} non-zero")
say(f"   minimum non-zero modulus: {min(mods):.10f}    (C-52 predicts exactly 1)")
say(f"   any value strictly inside (0,1)? {any(1e-12<m<1-1e-9 for m in mods)}")
say(f"   -> {'CONTROL PASSES: the quantisation is reproduced' if abs(min(mods)-round(min(mods)))<1e-9 else 'CONTROL FAILS -- instrument broken, no conclusion'}")

# ---------------------------------------------------------------- the metastable record
say("")
say("2. THE METASTABLE RECORD — the T-29 grain, and a family tuned from exact to metastable.")
say("   The record is the SLOW MODE of the Liouvillian, not an exact Pauli. Its projector onto the")
say("   long-lived manifold is the object the trace ratio must be taken against.")
T=300.0; K_u=2.0e5; V0=1.26e-24; f0=1e9
say(f"   {'E_b/kT':>9}{'dE/kT':>9}{'tau (s)':>13}{'record rate (1/s)':>20}{'<R>_ss':>16}{'in Z[i]?':>10}")
sz=np.array([[1,0],[0,-1]],dtype=complex); sp=np.array([[0,1],[0,0]],dtype=complex); sm=sp.conj().T
vals=[]
# ERRATUM: a first version varied only the BARRIER E_b and reported the interval as "populated
# continuously by the barrier height" while its own table printed the SAME value six times.
# <R>_ss = tanh(dE/2kT) depends on the SPLITTING, not the barrier. Both are varied here, and the
# splitting is the one that moves the value.
for scale, dEk in ((1.0,3.0),(1.0,1.0),(1.0,0.3),(1.0,0.1),(1.0,0.01),(1.0,1e-4),(1.0,1e-8),
                   (0.5,3.0),(0.25,3.0)):
    E_b=K_u*V0*scale; dE=dEk*G.KB*T
    H=-(dE/2)*sz
    gd=f0*np.exp(-(E_b+dE/2)/(G.KB*T)); gu=f0*np.exp(-(E_b-dE/2)/(G.KB*T))
    Ls=[np.sqrt(gd)*sm, np.sqrt(gu)*sp]
    c=G.clause_ii(H,Ls,sz,np.inf)
    Lv=G.liouvillian(H,Ls); wl,Vl=np.linalg.eig(Lv)
    j=int(np.argmin(np.abs(wl))); rho=Vl[:,j].reshape(2,2,order='F'); rho=rho/np.trace(rho)
    ratio=float(np.real(np.trace(rho@sz)))
    vals.append(ratio)
    say(f"   {E_b/(G.KB*T):>9.1f}{dEk:>9.0e}{c['tau']:>13.3e}{c['rate']:>20.4e}{ratio:>16.9f}"
        f"{str(abs(ratio-round(ratio))<1e-9):>10}")
say("")
say(f"   closed form tanh(dE/2kT) at dE = 1e-8 kT: {np.tanh(1e-8/2):.3e}  — the value goes to ZERO")
say(f"   CONTINUOUSLY as the splitting shrinks, taking EVERY value in (0,1) on the way.")
say("")
inside=[v for v in vals if 1e-12<abs(v)<1-1e-9]
say(f"   values strictly inside (0,1) in modulus: {len(inside)} of {len(vals)}")
if inside: say(f"   smallest such: {min(abs(v) for v in inside):.12f}")
say("")
say("="*100); say("  READ — generated from the numbers above"); say("="*100)
if abs(min(mods)-round(min(mods)))>1e-9:
    say("  THE CONTROL FAILED. No conclusion may be drawn.")
elif inside:
    say("  THE QUANTISATION IS AN ARTIFACT OF EXACTNESS. On the exact carrier every integer-coefficient")
    say(f"  polynomial lands on Z[i] with minimum non-zero modulus {min(mods):.6f}; on the metastable")
    say(f"  record the same quantity takes {len(inside)} values strictly inside (0,1), the smallest")
    say(f"  {min(abs(v) for v in inside):.3e}. THE INTERVAL C-52 CALLS EMPTY IS NOT EMPTY FOR A REAL")
    say("  RECORD. It is <R>_ss = tanh(dE/2kT), a CONTINUOUS function of the SPLITTING -- not of the")
    say("  barrier, which moves the lifetime over 25 orders of magnitude and leaves the value fixed.")
    say("  Every value in (0,1) is reachable, and arbitrarily small ones by arbitrarily small dE.")
    say("")
    say("  CONSEQUENCE. C-52's physical conclusion -- that a gravity-strength residual has NO SLOT --")
    say("  was a theorem about DEF-A's corner, not about the world. The exclusion does not inherit.")
else:
    say("  The quantisation SURVIVES on the metastable record: no value fell strictly inside (0,1).")
    say("  C-52's exclusion inherits, and the gravity-strength residual has no slot in the world either.")
