"""O-46: IS THERE A GAP BETWEEN 'NOT A FUNCTION OF H' AND 'UNCORRELATED WITH H'?

The principal: "Are we sure that we aren't importing an assumption and restriction here that is not
applicable for the record level?"

BOTH READINGS OF 'ADMISSIBLE' IMPORT. DEF-A ([U,H]=0) imports that operations must conserve what H
generates. O-44's physical-channel fix imports thermodynamics -- energy, a reservoir, a temperature,
Landauer's kT ln 2, and an external agent doing work. Neither is given at the record level, and there
is no 'outside' there to do the work. Worse: 'admissible' decides what counts as WRITING, which is a
fact about the PROCESS, and finding the process is the program's stated goal. Defining it in advance
assumes the answer.

BUT ONE THING IS CHECKABLE WITHOUT SETTLING THAT, AND IT IS INTERNAL TO THE FRAMEWORK'S OWN REASONING.

  CORE_FRAMEWORK on clause (iii): "If R = f(H) its value is fixed by the energy: knowing the energy
  tells you the record, so it isn't independent information."   -> R is NOT A FUNCTION OF H.

  What clause (iv) enforces once derived: Tr(P_E R) = 0 on EVERY eigenspace of H.
                                          -> R is COMPLETELY UNCORRELATED WITH H.

THOSE ARE DIFFERENT CONDITIONS. A record can be PARTIALLY correlated with H without being a function
of it. This lane measures whether the gap between them is non-empty, and whether operators in the gap
carry a configuration energy.

ERRATUM, v1: v1 built H from all-distinct single-qubit fields, giving a NON-DEGENERATE spectrum
(multiplicities all 1). The framework's own P-1 says clause (iii) implies H is degenerate, so that
carrier can hold NO RECORD AT ALL and the test could only return zero -- which it did, 0 of 256 and
0 of 200000. A venue where the effect cannot appear. H is now built DEGENERATE, and the degeneracy is
printed before any count is read.

NOTHING HERE IMPORTS. No energy is called work, no reservoir, no temperature, no agent. The only
objects are H, the operators, and counting. Tr(P_E R) is read as 'the record's imbalance within a
spectral block of H', not as an energy."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); I2=np.eye(2)
def word(n,s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
def pauli(n,v):
    M=np.array([[1]],dtype=complex)
    for i in range(n):
        x,z=v[i],v[n+i]
        P=I2 if (x,z)==(0,0) else (X if (x,z)==(1,0) else (Z if (x,z)==(0,1) else 1j*X@Z))
        M=np.kron(M,P)
    return M
say("="*104); say("O-46   THE GAP BETWEEN 'NOT A FUNCTION OF H' AND 'UNCORRELATED WITH H'"); say("="*104)
say("  Both are stated on the SAME objects: H and R. No energy, no reservoir, no agent, no cost.")
say("  Tr(P_E R) is read as the record's IMBALANCE within a spectral block of H.")
say("")
for n in (3,4):
    say("-"*104); say(f"  carrier: {n} qubits, dim {2**n}"); say("-"*104)
    # A carrier with genuine spectral structure -- NOT permutation-symmetric (D-22)
    # DEGENERATE by construction (P-1: clause (iii) requires it), and NOT permutation-symmetric (D-22):
    # a chain of ZZ couplings with distinct strengths.
    H=sum((1.0+0.5*j)*word(n,''.join('Z' if i in (j,j+1) else 'I' for i in range(n))) for j in range(n-1))
    w,V=np.linalg.eigh(H)
    es=[]
    for i,x in enumerate(w):
        if es and abs(x-es[-1][0])<1e-8: es[-1][1].append(i)
        else: es.append([x,[i]])
    say(f"    eigenvalues: {len(es)} distinct, multiplicities {[len(ix) for _,ix in es]}")
    if all(len(ix)==1 for _,ix in es):
        say("    H IS NON-DEGENERATE -- by P-1 no record can exist here at all. NO TEST IS POSSIBLE."); say(""); continue
    # enumerate every Hermitian +-1 operator that is a signed sum of eigenprojectors of H
    # (these are exactly the R with [H,R] = 0 that are diagonal in H's eigenbasis)
    tot=0; i_ii=0; iii_ok=0; iv_ok=0; gap=[]
    # to allow R NOT constant on an eigenspace, split each degenerate eigenspace too
    dims=[len(ix) for _,ix in es]
    if any(d>1 for d in dims):
        say("    (degenerate blocks present -- signs may vary inside a block)")
    cells=[]
    for k,(val,ix) in enumerate(es):
        for j in ix: cells.append((k,j))
    if len(cells)>14:
        say(f"    {len(cells)} cells -- enumerating a random 200000 sample of sign patterns")
        rng=np.random.default_rng(5); pats=[rng.integers(0,2,len(cells))*2-1 for _ in range(200000)]
    else:
        pats=[np.array(p) for p in itertools.product((1,-1),repeat=len(cells))]
    for sgn in pats:
        R=np.zeros((2**n,2**n),dtype=complex)
        for c,(k,j) in enumerate(cells):
            v=V[:,j:j+1]; R+= sgn[c]*(v@v.conj().T)
        tot+=1
        # (i) bit  (ii) [H,R]=0 -- both automatic by construction; verify once
        if tot==1:
            assert np.linalg.norm(R@R-np.eye(2**n))<1e-8 and np.linalg.norm(H@R-R@H)<1e-8
        i_ii+=1
        # (iii) NOT a function of H: R is not constant on SOME eigenspace
        notfn=any(len(set(int(sgn[c]) for c,(k,j) in enumerate(cells) if k==kk))>1 for kk in range(len(es)))
        if not notfn: continue
        iii_ok+=1
        trs=[abs(sum(int(sgn[c]) for c,(k,j) in enumerate(cells) if k==kk)) for kk in range(len(es))]
        if max(trs)<1e-9: iv_ok+=1
        else: gap.append(max(trs))
    say(f"    operators tested (bit + commutes with H)          : {i_ii}")
    say(f"    ALSO not a function of H   -- clause (iii)        : {iii_ok}")
    say(f"    ALSO uncorrelated with H   -- clause (iv) derived : {iv_ok}")
    say(f"    IN THE GAP (iii yes, iv no)                       : {len(gap)}")
    if gap:
        say(f"    their imbalance max_E |Tr(P_E R)| ranges over     : {min(gap)} .. {max(gap)}")
    say(f"    -> {'THE GAP IS NON-EMPTY' if gap else 'the gap is EMPTY -- the two conditions coincide here'}")
    say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)
say("  THE GAP IS LARGE. Of the operators that are a bit, commute with H, and are NOT a function of H")
say("  -- clauses (i), (ii) and (iii) -- clause (iv) as derived keeps only 16 of 240 at n=3 and 772 of")
say("  199204 at n=4. It DISCARDS 93% and 99.6% of them.")
say("")
say("  WHAT IT DISCARDS ARE OPERATORS CORRELATED WITH H WITHOUT BEING A FUNCTION OF H. Their block")
say("  imbalance |Tr(P_E R)| is 2 throughout: not a function of H, and not uncorrelated with it either.")
say("  The framework's stated reason for clause (iii) -- 'if R = f(H) its value is fixed by the energy'")
say("  -- rules out only the FUNCTION case. The derived form of clause (iv) rules out ALL correlation.")
say("")
say("  BUT THE GAP TURNS OUT NOT TO BE THE PLACE TO LOOK. O-47 shows configuration energy and free")
say("  energy-conserving writers coexist WITHOUT entering this gap at all: on H = Z1 Z2 both records")
say("  satisfy clause (iv) exactly, with Tr(P_E R) = [0,0], and the energy sits in their CORRELATION,")
say("  which clause (iv) never constrains. No clause has to be weakened and nothing imported.")
say("")
say("  So this lane's measurement stands as a fact about the two conditions -- they are genuinely")
say("  different and the difference is large -- and it is NOT the route to a configuration energy.")
