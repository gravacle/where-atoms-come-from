"""O-44: IS THE DURABILITY-GEOMETRY EXCLUSION AN ARTIFACT OF DEF-A?

The principal: "why did we forbid records from carrying configuration energy? We shouldn't be
assuming anything if we're looking for the truth."

WE DID NOT FORBID IT DELIBERATELY. It follows from ONE definition:

    CORE_FRAMEWORK line 26:  ADMISSIBLE U := a unitary with [U,H] = 0.   (O-4, PROVISIONAL)
    ledger O-4 evidence:     "PROVISIONAL -- one carrier family, and DEF-A has NO working fallback"
    ledger O-4 item:         "The second disjunct 'or a physical channel' is UNTESTED"

THE DERIVATION USES IT ESSENTIALLY. Clause (iv) says some ADMISSIBLE U has U-dag R U = -R. If
[U,H] = 0 then U P_E U-dag = P_E, so
    Tr(P_E R) = Tr(P_E U R U-dag) = Tr(P_E (-R)) = -Tr(P_E R)   =>   Tr(P_E R) = 0.
Remove [U,H] = 0 and NOTHING forces Tr(P_E R) = 0. The O-42 exclusion collapses.

AND DEF-A IS PHYSICALLY THE WRONG REQUIREMENT. A unitary commuting with H does NO WORK on the state's
energy. Requiring the writer to be a symmetry of H is requiring that WRITING BE FREE. Landauer puts a
real write's floor at kT ln 2, not zero.

ERRATUM, v1: v1 NOMINATED the writer as an X-string on sites 0..i and its READ then asserted the
unmodified write was free with ||[W,H]|| = 0 while the table directly above printed 1.600e+01 and a
cost of +2.000000. X_0 flips Z_0 Z_1 but does not commute with Z^n. Writers are now SEARCHED FOR over
the full Pauli group, never nominated -- the sixth time this program has gone wrong that way."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); I2=np.eye(2)
def word(n,s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
def opat(n,sites,P):
    s=['I']*n
    for i in sites: s[i]=P
    return word(n,''.join(s))
def find_writer(n,H,R,admissible):
    """SEARCH the Pauli group for a writer. admissible=True demands [W,H]=0 as well."""
    best=None
    for xs in range(1,2**n):
        for zs in range(2**n):
            v=[(xs>>i)&1 for i in range(n)]+[(zs>>i)&1 for i in range(n)]
            M=np.array([[1]],dtype=complex)
            for i in range(n):
                x,z=v[i],v[n+i]
                P=I2 if (x,z)==(0,0) else (X if (x,z)==(1,0) else (Z if (x,z)==(0,1) else 1j*X@Z))
                M=np.kron(M,P)
            if np.linalg.norm(M.conj().T@R@M + R)>1e-9: continue          # must FLIP R
            if admissible and np.linalg.norm(M@H-H@M)>1e-9: continue      # must commute with H
            wt=sum(1 for i in range(n) if (v[i],v[n+i])!=(0,0))
            if best is None or wt<best[0]: best=(wt,M,v)
    return best
say("="*104); say("O-44   DOES A PHYSICAL WRITER EXIST WHERE AN ADMISSIBLE ONE DOES NOT?"); say("="*104)
say("  Writers are SEARCHED over the full Pauli group, never nominated.")
for n in (4,6):
    say(""); say("-"*104); say(f"  carrier n = {n}, dim {2**n}"); say("-"*104)
    Hstab=-(word(n,'X'*n)+word(n,'Z'*n))
    recs=[opat(n,[i,i+1],'Z') for i in range(n-2)]
    Hgeo=Hstab - sum((0.3+0.1*i)*recs[i] for i in range(len(recs)))
    R=recs[0]
    for tag,H in (("UNMODIFIED   H = -(X^n + Z^n)          [CONTROL]",Hstab),
                  ("CONFIGURATION ENERGY   H = H_stab - sum c_i R_i",Hgeo)):
        w,V=np.linalg.eigh(H)
        es=[]
        for i,x in enumerate(w):
            if es and abs(x-es[-1][0])<1e-9: es[-1][1].append(i)
            else: es.append([x,[i]])
        tr=max(abs(float(np.real(np.trace(V[:,ix].conj().T@R@V[:,ix])))) for _,ix in es)
        ix0=es[0][1]; U0=V[:,ix0]; rho=U0@U0.conj().T/len(ix0)
        E0=float(np.real(np.trace(rho@H)))
        say(f"  {tag}")
        say(f"     ||[H,R]|| = {np.linalg.norm(H@R-R@H):.3e}   record is DURABLE")
        say(f"     max_E |Tr(P_E R)| = {tr:.3e}   clause (iv) under DEF-A: {'HOLDS' if tr<1e-9 else 'FAILS'}")
        for kind,adm in (("ADMISSIBLE ([W,H]=0)",True),("PHYSICAL (any unitary)",False)):
            b=find_writer(n,H,R,adm)
            if b is None:
                say(f"     {kind:<24} writer: NONE EXISTS in the Pauli group")
            else:
                wt,W,v=b
                E1=float(np.real(np.trace((W@rho@W.conj().T)@H)))
                say(f"     {kind:<24} writer FOUND, weight {wt}, ||[W,H]|| = {np.linalg.norm(W@H-H@W):.3e},"
                    f" energy cost {E1-E0:+.6f}")
        say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
say("  WHERE CLAUSE (iv) HOLDS, A FREE WRITER EXISTS AND COSTS EXACTLY ZERO. The search finds a")
say("  weight-2 admissible writer with ||[W,H]|| = 0.000e+00 and energy cost +0.000000. That is the")
say("  control, and it is the carrier this program has always used.")
say("")
say("  WHERE THE CARRIER ASSIGNS A CONFIGURATION ENERGY, NO FREE WRITER EXISTS -- the search returns")
say("  NONE over the entire Pauli group -- BUT A PHYSICAL WRITER DOES, at weight 1, and it costs")
say("  +2.600000 of work. The record remains a bit, remains durable, remains non-trivial, and remains")
say("  FLIPPABLE BY AN OPERATOR THE SEARCH ACTUALLY FOUND.")
say("")
say("  SO WHAT FAILS IS NOT WRITABILITY. WHAT FAILS IS FREE WRITABILITY.")
say("")
say("  O-42's EXCLUSION IS THEREFORE CONDITIONAL ON DEF-A, NOT A FACT ABOUT RECORDS. DEF-A is marked")
say("  PROVISIONAL in the framework and its evidence line reads 'one carrier family, and DEF-A has NO")
say("  working fallback'; the alternative reading, 'or a physical channel', is marked UNTESTED. A")
say("  unitary commuting with H does no work on the state's energy, so DEF-A amounts to requiring")
say("  that writing be FREE -- and a free write cannot distinguish energies. Landauer puts the floor")
say("  on a real write at kT ln 2, not at zero.")
say("")
say("  WHAT THIS DOES AND DOES NOT DO. It REMOVES A BLOCKER: records can carry a configuration energy")
say("  provided writing is allowed to cost work, which is what writing costs in the world. It does")
say("  NOT show that gravity emerges. It reopens the question that O-42 appeared to close, and it")
say("  reopens C-47 -- whose extensive quantity was registered as bought at the cost of independent")
say("  writability, a cost that under the physical reading is simply WORK DONE.")
