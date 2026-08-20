"""O-47: CONFIGURATION ENERGY *AND* ENERGY-CONSERVING WRITERS. BOTH, AT ONCE.

The principal: "Why can't we have configuration energy AND defining a writer as an energy-conserving
symmetry?"

WE CAN. The apparent conflict came from applying clause (iv) to the wrong object.

Clause (iv) constrains EACH RECORD: Tr(P_E R_i) = 0 on every eigenspace -- each record's value must be
balanced at every energy. It says NOTHING about the JOINT configuration. Two records can each be
individually balanced while their CORRELATION is definite on each eigenspace, and the correlation is
then exactly what carries the energy.

THE CONSTRUCTION.  H = c * Z_1 Z_2.
  eigenspace E_+ = span{|++>, |-->}   Tr(P_+ Z_1) = 1 + (-1) = 0,  Tr(P_+ Z_2) = 0
  eigenspace E_- = span{|+->, |-+>}   Tr(P_- Z_1) = 1 + (-1) = 0,  Tr(P_- Z_2) = 0
  so CLAUSE (iv) HOLDS FOR EACH RECORD -- each is balanced at every energy.
  but Z_1 Z_2 = +1 on E_+ and -1 on E_-, so THE PAIR CORRELATION IS DEFINITE AT EACH ENERGY.

So each record can be flipped by an energy-conserving symmetry, and the PAIR still costs energy to
rearrange. That is an INTERACTION ENERGY BETWEEN RECORDS coexisting with free single-record writes.

Everything is verified, nothing nominated: writers are SEARCHED over the full Pauli group.
CONTROLS IN THE SAME TABLE (D-15): a carrier with no configuration energy, where the same quantities
must come out zero."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex); I2=np.eye(2)
def word(n,s):
    M=np.array([[1]],dtype=complex)
    for c in s: M=np.kron(M,{'I':I2,'X':X,'Z':Z}[c])
    return M
def paulis(n):
    for xs in range(2**n):
        for zs in range(2**n):
            v=[(xs>>i)&1 for i in range(n)]+[(zs>>i)&1 for i in range(n)]
            M=np.array([[1]],dtype=complex)
            for i in range(n):
                x,z=v[i],v[n+i]
                P=I2 if (x,z)==(0,0) else (X if (x,z)==(1,0) else (Z if (x,z)==(0,1) else 1j*X@Z))
                M=np.kron(M,P)
            yield v,M
def eigblocks(H,d):
    w,V=np.linalg.eigh(H); es=[]
    for i,x in enumerate(w):
        if es and abs(x-es[-1][0])<1e-9: es[-1][1].append(i)
        else: es.append([x,[i]])
    return w,V,es
say("="*104); say("O-47   CONFIGURATION ENERGY AND ENERGY-CONSERVING WRITERS, TOGETHER"); say("="*104)
n=2
for tag,H,expect in (("H = 1.0 * Z1 Z2      [interaction energy]", 1.0*word(2,'ZZ'), True),
                     ("H = 0                [CONTROL: no energy]", 0.0*word(2,'ZZ'), False)):
    say(""); say("-"*104); say(f"  {tag}"); say("-"*104)
    w,V,es=eigblocks(H,2**n)
    say(f"    spectrum: {[ (round(float(v),4), len(ix)) for v,ix in es ]}")
    R1=word(2,'ZI'); R2=word(2,'IZ'); P=word(2,'ZZ')
    for name,R in (("R1 = Z1",R1),("R2 = Z2",R2)):
        trs=[float(np.real(np.trace(V[:,ix].conj().T@R@V[:,ix]))) for _,ix in es]
        nonconst=any(np.linalg.norm(V[:,ix].conj().T@R@V[:,ix]
                     - (np.trace(V[:,ix].conj().T@R@V[:,ix])/len(ix))*np.eye(len(ix)))>1e-9 for _,ix in es)
        say(f"    {name}:  ||[H,R]|| = {np.linalg.norm(H@R-R@H):.1e}   Tr(P_E R) per block = {[round(t,6) for t in trs]}")
        say(f"             clause (iii) non-constant on some block: {nonconst}"
            f"   clause (iv) {'HOLDS' if max(abs(t) for t in trs)<1e-9 else 'FAILS'}")
    trsP=[float(np.real(np.trace(V[:,ix].conj().T@P@V[:,ix])))/len(ix) for _,ix in es]
    say(f"    PAIR CORRELATION R1*R2 = Z1Z2:  value per block = {[round(t,6) for t in trsP]}")
    say(f"             {'DEFINITE and DIFFERENT across blocks -- IT CARRIES THE ENERGY' if len(set(round(t,6) for t in trsP))>1 else 'same on every block -- carries no energy'}")
    say("")
    say("    ADMISSIBLE WRITERS, searched over the full Pauli group ([W,H]=0 AND W flips the record):")
    for name,R in (("R1",R1),("R2",R2)):
        found=None
        for v,W in paulis(n):
            if np.linalg.norm(W.conj().T@R@W + R)>1e-9: continue
            if np.linalg.norm(W@H-H@W)>1e-9: continue
            wt=sum(1 for i in range(n) if (v[i],v[n+i])!=(0,0))
            if found is None or wt<found[0]: found=(wt,W,v)
        if found is None:
            say(f"       {name}: NO admissible writer exists")
        else:
            wt,W,v=found
            ix0=es[0][1]; U0=V[:,ix0]; rho=U0@U0.conj().T/len(ix0)
            dE=float(np.real(np.trace((W@rho@W.conj().T)@H)))-float(np.real(np.trace(rho@H)))
            flipsP=np.linalg.norm(W.conj().T@P@W - P)>1e-9
            say(f"       {name}: FOUND weight {wt}, ||[W,H]|| = {np.linalg.norm(W@H-H@W):.1e},"
                f" energy change {dE:+.6f}, changes the PAIR? {flipsP}")
    say("")
    say("    COST OF CHANGING THE PAIR CORRELATION (search for ANY Pauli that flips Z1Z2):")
    best=None
    for v,W in paulis(n):
        if np.linalg.norm(W.conj().T@P@W + P)>1e-9: continue
        ix0=es[0][1]; U0=V[:,ix0]; rho=U0@U0.conj().T/len(ix0)
        dE=float(np.real(np.trace((W@rho@W.conj().T)@H)))-float(np.real(np.trace(rho@H)))
        adm=np.linalg.norm(W@H-H@W)
        wt=sum(1 for i in range(n) if (v[i],v[n+i])!=(0,0))
        if best is None or (adm<1e-9, -wt)>(best[0]<1e-9,-best[1]): best=(adm,wt,dE)
    if best is None: say("       no Pauli flips the pair correlation")
    else:
        adm,wt,dE=best
        say(f"       best flipper: weight {wt}, ||[W,H]|| = {adm:.1e}"
            f" ({'ADMISSIBLE' if adm<1e-9 else 'NOT admissible -- costs work'}), energy change {dE:+.6f}")
say(""); say("="*104); say("  READ -- from the numbers above"); say("="*104)
say("  BOTH, AT ONCE, AND WITH NOTHING IMPORTED.")
say("")
say("  On H = Z1 Z2 each record is BALANCED AT EVERY ENERGY -- Tr(P_E R) = [0, 0] for both -- so")
say("  CLAUSE (iv) HOLDS FOR EACH, and the search finds a weight-2 ADMISSIBLE writer for each with")
say("  ||[W,H]|| = 0 and energy change exactly +0.000000. Single records flip for FREE.")
say("")
say("  AND THE PAIR CORRELATION Z1 Z2 IS DEFINITE AND DIFFERENT ACROSS BLOCKS, -1 and +1. It is what")
say("  carries the energy. The free single-record writers do NOT change it, and the best operator")
say("  that does change it is NOT admissible and costs +2.000000.")
say("")
say("  THE CONTROL SHOWS THE INSTRUMENT WORKS. On H = 0 the pair correlation is the same on every")
say("  block, it carries no energy, and the admissible writers -- now weight 1 -- DO change it, for")
say("  free, because there is no energy to protect.")
say("")
say("  SO THE EXCLUSION IN O-42 CAME FROM APPLYING CLAUSE (iv) TO THE WRONG OBJECT. Clause (iv)")
say("  constrains EACH RECORD: Tr(P_E R_i) = 0. It says NOTHING about the JOINT CONFIGURATION. The")
say("  configuration energy lives in the CORRELATIONS BETWEEN records, which clause (iv) never")
say("  touches -- and an interaction energy between records is exactly what a source needs.")
say("")
say("  AND NOTHING WAS IMPORTED TO GET HERE. No reservoir, no temperature, no Landauer bound, no")
say("  external agent, no relaxation of any clause. Only H, the records, and the five clauses as")
say("  written. O-44's physical-channel rescue is not needed: DEF-A can stand exactly as it is.")
