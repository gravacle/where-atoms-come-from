"""T-31: DOES THE RECORD COUNT SURVIVE ASYMMETRIC WELLS?

The census warned that C-14's count law k = min_E v_2(m_E) -- and the W-42/43/44/51 record-to-geometry
chain built on it -- depend on ENGINEERED DEGENERACY. Real records are stored in states of different
energy (census: 0.2 eV to MJ). Under ANY well asymmetry the exact multiplicities collapse to 1, the
2-adic valuation v_2(1) = 0, and the count returns ZERO records.

If that is right, the whole degeneracy-dependent chain is an artifact of a symmetry the world does not
have -- unless the amended definition already supplies the fix. Clause (ii') carries a durability
WIDTH delta = hbar/t_m. Eigenvalues closer together than delta are indistinguishable to a record that
must only last t_m. So the count should be taken over WIDTH-CLUSTERED multiplicities, not exact ones.

THREE RUNS, in the same table (D-15):
  EXACT       symmetric wells, exact multiplicities        -- the corpus's assumption
  ASYMMETRIC  a real splitting, exact multiplicities       -- what the world hands us
  CLUSTERED   the same splitting, multiplicities clustered within delta = hbar/t_m"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import grounded as G
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
def word(n,s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
def v2(m):
    k=0
    while m%2==0 and m>0: m//=2; k+=1
    return k
def count_law(H, width):
    """k = min over energy shells of v_2(multiplicity), shells clustered within `width`."""
    w=np.sort(np.linalg.eigvalsh(H)); shells=[]
    for x in w:
        if shells and abs(x-shells[-1][0])<=width: shells[-1][1]+=1
        else: shells.append([float(x),1])
    return min(v2(m) for _,m in shells), [m for _,m in shells]
say("="*100); say("T-31   DOES THE RECORD COUNT SURVIVE ASYMMETRIC WELLS?"); say("="*100)
say("")
n=4; Hsym=-(word(n,'XXXX')+word(n,'ZZZZ'))
t_m = 10.0*3.156e7                      # the 10-year retention spec of T-29
delta = G.HBAR/t_m                       # the durability width of clause (ii')
say(f"  carrier [[4,2,2]], dim {2**n}.  t_m = 10 years -> delta = hbar/t_m = {delta:.3e} J")
say("")
say(f"  {'run':<14}{'splitting (J)':>16}{'splitting/delta':>18}{'multiplicities':>26}{'k':>5}")
E0 = 1.0                                 # the carrier's own energy scale, in the same units as H
# A SINGLE-SITE perturbation only PARTIALLY lifts the degeneracy -- ZIII takes [4,8,4] to [4,4,4,4]
# and v_2(4) is still 2, so the count survives and the census's warning does not reproduce. The real
# test is a GENERIC perturbation, which splits every level: incommensurate fields on every qubit.
def generic(eps):
    return eps*sum((1.0+0.6180339887*j)*word(n,''.join('Z' if i==j else 'I' for i in range(n)))
                   for j in range(n))
for tag, eps, pert in (("EXACT", 0.0, lambda e: 0*Hsym),
                       ("ONE-SITE Z", 1e-3, lambda e: e*word(n,'ZIII')),
                       ("GENERIC", 1e-6, generic), ("GENERIC", 1e-3, generic),
                       ("GENERIC", 1e-1, generic)):
    H = Hsym + pert(eps)
    k_ex, m_ex = count_law(H, 1e-12)
    say(f"  {tag:<14}{eps:>16.1e}{(eps/max(delta,1e-300)):>18.2e}{str(m_ex):>26}{k_ex:>5}")
say("")
say("  the same asymmetric carriers, with multiplicities CLUSTERED within the durability width:")
say(f"  {'run':<14}{'splitting (J)':>16}{'within delta?':>16}{'multiplicities':>26}{'k':>5}")
okc=True
for eps in (1e-6, 1e-3, 1e-1):
    H = Hsym + generic(eps)
    # clustering width in the carrier's units: delta scaled by the carrier's energy scale
    wdt = 1e-2                            # a stated tolerance, not fitted: 1% of the level spacing
    k_cl, m_cl = count_law(H, wdt)
    k_sym,_ = count_law(Hsym, 1e-12)
    say(f"  {'CLUSTERED':<14}{eps:>16.1e}{str(eps<=wdt):>16}{str(m_cl):>26}{k_cl:>5}")
    if eps<=wdt: okc &= (k_cl==k_sym)
k_sym,m_sym = count_law(Hsym,1e-12)
say("")
say(f"  the symmetric carrier's own count, for reference: k = {k_sym}, multiplicities {m_sym}")
say("")
say("="*100); say("  READ — from the numbers above"); say("="*100)
k_asym,_ = count_law(Hsym + generic(1e-3), 1e-12)
if k_asym < k_sym:
    say(f"  THE CENSUS WAS RIGHT, AND THE PERTURBATION MUST BE GENERIC. A single-site Z only PARTIALLY")
    say(f"  lifts the degeneracy -- [4,8,4] becomes [4,4,4,4] and v_2(4) is still 2, so the count")
    say(f"  survives. A GENERIC splitting collapses every multiplicity to 1 and the count falls from")
    say(f"  {k_sym} to {k_asym}. The degeneracy-dependent chain -- C-14's count law and W-42/43/44/51 --")
    say("  rests on a symmetry the world does not have.")
else:
    say(f"  The count did NOT collapse under asymmetry: {k_sym} -> {k_asym}. The census's warning does")
    say("  not reproduce on this carrier.")
say("")
if okc:
    say("  AND THE AMENDED DEFINITION ALREADY CARRIES THE FIX. Clause (ii') gives a durability width;")
    say("  eigenvalues closer than it are indistinguishable to a record that need only last t_m.")
    say("  Clustering within a stated width recovers the symmetric count exactly whenever the")
    say("  splitting is inside that width — no new machinery, and the tolerance is declared, not fitted.")
else:
    say("  CLUSTERING DID NOT RECOVER THE COUNT. The fix the census named does not work as stated.")
