# W-25 — THE CLOSED-CURVE CARRIER, built in the PHYSICAL BASIS directly.
# Gauss operators are diagonal in the electric basis, so the physical sector is a combinatorial
# condition on link values, not a matrix product. Enumerate, filter, then build operators there.
#
# CARRIER: wheel W_4. Hub 0 = interior. Rim 1,2,3,4 = A CLOSED CURVE = the boundary.
#   spokes 0:(0,1) 1:(0,2) 2:(0,3) 3:(0,4)   rim 4:(1,2) 5:(2,3) 6:(3,4) 7:(4,1)
import numpy as np, itertools
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(2,3),(3,4),(4,1)]
SPOKE=[0,1,2,3]; RIM=[4,5,6,7]; L=len(E); V=5

def on_cycle(i):
    a,b=E[i]; adj={}
    for j,(u,v) in enumerate(E):
        if j==i: continue
        adj.setdefault(u,[]).append(v); adj.setdefault(v,[]).append(u)
    seen={a}; st=[a]
    while st:
        x=st.pop()
        for y in adj.get(x,[]):
            if y not in seen: seen.add(y); st.append(y)
    return b in seen
print("  THE TEST EVERY PREVIOUS CARRIER FAILED:")
for i in range(L):
    print(f"    link {i} = {E[i]}  on a cycle: {str(on_cycle(i)):>5}   "
          f"{'RIM (boundary)' if i in RIM else 'spoke (interior)'}")
print(f"  -> boundary links on a cycle: {sum(on_cycle(i) for i in RIM)} of {len(RIM)}. "
      f"A flux CAN thread this boundary.\n")

def physical(N):
    """states are tuples of link values in Z_N; Gauss at v: sum out - sum in = 0 mod N"""
    keep=[]
    for s in itertools.product(range(N),repeat=L):
        ok=True
        for v in range(V):
            t=0
            for i,(a,b) in enumerate(E):
                if a==v: t+=s[i]
                if b==v: t-=s[i]
            if t%N: ok=False; break
        if ok: keep.append(s)
    return keep

def loop_op(states,idx,moves,N):
    """DEFECT RECORDED: the first version applied a plaquette as two sequential shifts, which
    passes through a NON-PHYSICAL intermediate and raised KeyError. A loop must be applied as ONE
    simultaneous move. moves = [(link, +1 or -1)] traversing the loop with orientation."""
    D=len(states); M=np.zeros((D,D),dtype=complex)
    for j,s in enumerate(states):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        k=idx.get(tuple(t))
        if k is None: return None                 # not gauge invariant: leaves the physical sector
        M[k,j]=1.0
    return M
def elec_op(states,links,N):
    """sum over named links of the clock operator: diagonal, entries omega^k"""
    w=np.exp(2j*np.pi/N)
    return np.diag([sum(w**s[l] for l in links) for s in states])

for N in (2,3,4):
    st=physical(N); idx={s:j for j,s in enumerate(st)}; D=len(st)
    # rim traversed 1->2->3->4->1: links 4,5,6 forward, link 7 is (4,1) so forward too
    Wrim=loop_op(st,idx,[(4,+1),(5,+1),(6,+1),(7,+1)],N)
    # triangle 0->1->2->0: link 0 (0,1) forward, link 4 (1,2) forward, link 1 (0,2) BACKWARD
    Wtri=loop_op(st,idx,[(0,+1),(4,+1),(1,-1)],N)
    ok=lambda M: "leaves physical sector" if M is None else \
        ("NONTRIVIAL" if np.linalg.norm(M-np.eye(D))>1e-9 else "= identity")
    print(f"  N={N}:  full {N**L:>6}   physical {D:>5}   rim loop: {ok(Wrim):<22} triangle: {ok(Wtri)}")
    if Wrim is not None and Wtri is not None:
        # is the rim loop the product of the four interior triangles? (Stokes on the disc)
        T=[loop_op(st,idx,[(k,+1),(4+k,+1),((k+1)%4,-1)],N) for k in range(4)]
        if all(x is not None for x in T):
            prod=T[0]
            for x in T[1:]: prod=prod@x
            print(f"         rim loop == product of the 4 interior triangles?  "
                  f"{np.linalg.norm(prod-Wrim) < 1e-9}   (discrete Stokes)")
print()
print("  The rim loop is a nontrivial gauge-invariant operator on the physical sector.")
print("  THAT is what a boundary is: a closed curve carrying a flux. No previous carrier had one.")
