"""PF-3 part B, rewritten. v1 used SINGLE LINKS as regions; the toric test that established
the result used regions of up to 24 edges, and a single link cannot discriminate. v1 also
used the theta graph, which has d = 1 -- a carrier where clause (v) genuinely fails, so it
cannot test a definition either.

Carrier here: BOUQUET OF TWO TRIANGLES -- a pinch-point NON-MANIFOLD (G-10's class,
structurally different from the torus), 6 links, dim 64, cycles of length 3.
CONTRACTIBLE REGION on a graph = an edge subset containing NO cycle, i.e. a FOREST."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=complex); Z=np.array([[1,0],[0,-1]],dtype=complex)
E=[(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]; V=5; nL=len(E)
def op(d):
    M=np.array([[1]],dtype=complex)
    for l in range(nL): M=np.kron(M,d.get(l,I2))
    return M
STAR=[[l for l,(a,b) in enumerate(E) if a==v or b==v] for v in range(V)]
H=-sum(op({l:X for l in s}) for s in STAR)
w,Vec=np.linalg.eigh(H); gs=int(np.sum(np.abs(w-w[0])<1e-9))
say("="*98); say("PF-3B  CLAUSE (v) ON A NON-MANIFOLD CARRIER (bouquet of two triangles)"); say("="*98)
say(f"  V={V}  links={nL}  dim={2**nL}  ground degeneracy={gs}  (dim H_1 = 2 => expect 4)")
# the record: Z along one triangle (a cycle, hence a logical); its cycle length is 3
CYC=[0,1,2]
R=op({l:Z for l in CYC})
say(f"  record R = Z on the cycle {CYC}   ||[R,H]|| = {np.linalg.norm(R@H-H@R):.3e}   "
    f"{'PASS' if np.linalg.norm(R@H-H@R)<1e-9 else 'FAIL'}")
# contractible regions = FORESTS (edge subsets containing no cycle)
def has_cycle(sub):
    par=list(range(V))
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for l in sub:
        a,b=E[l]; ra,rb=find(a),find(b)
        if ra==rb: return True
        par[ra]=rb
    return False
regions=[list(s) for k in range(1,nL+1) for s in itertools.combinations(range(nL),k) if not has_cycle(s)]
say(f"  contractible regions (forests) of every size: {len(regions)}   "
    f"largest {max(len(r) for r in regions)} edges")
P1=[I2,X,Z,1j*(X@Z)]
anyf=admf=0; ex=None
for T in regions:
    for combo in itertools.product(range(4),repeat=len(T)):
        if all(c==0 for c in combo): continue
        A=op({l:P1[c] for l,c in zip(T,combo)})
        if np.linalg.norm(A.conj().T@R@A + R) < 1e-8:
            anyf+=1
            if np.linalg.norm(A@H-H@A) < 1e-8: admf+=1
            elif ex is None: ex=(T,combo)
say("")
say(f"  {'ANY-unitary flippers inside a contractible region':<56}{anyf:>10}   <- POSITIVE CONTROL")
say(f"  {'ADMISSIBLE flippers ([U,H]=0) inside a contractible region':<56}{admf:>10}")
say("")
if anyf==0:
    say("  VOID: the control did not fire -- no operator of any kind flips R in these regions,")
    say("        so this carrier cannot discriminate the definition. Not a result.")
elif admf==0:
    say("  CLAUSE (v) HOLDS UNDER DEF-A AND FAILS UNDER 'ANY UNITARY' -- the same behaviour as the")
    say("  torus, now on a NON-MANIFOLD carrier. O-4's definition is not an artefact of one family.")
else:
    say(f"  CLAUSE (v) FAILS UNDER DEF-A TOO ({admf} admissible flippers) -- O-4's definition does NOT")
    say("  transfer to this carrier. That is a finding against the anchor.")
sys.exit(0)
