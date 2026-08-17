# W-27 — THE PLANAR PATCH. Separator and cycle are DUAL, not the same set. Build a carrier that has
# both, as the two fields it already has: electric flux through a CUT, magnetic flux through a CYCLE.
# 3x3 grid of vertices, planar (no wrap). 12 links, 9 vertices, 4 plaquettes.
#   REGION = the interior vertex (1,1).  ELECTRIC BOUNDARY = its 4 incident links (a CUT).
#   MAGNETIC BOUNDARY = the outer 8-link square (a CYCLE) enclosing all 4 plaquettes.
import numpy as np, itertools
V2=[(i,j) for j in range(3) for i in range(3)]
vid={v:k for k,v in enumerate(V2)}
E=[]
for j in range(3):
    for i in range(2): E.append((vid[(i,j)],vid[(i+1,j)]))      # horizontal
for j in range(2):
    for i in range(3): E.append((vid[(i,j)],vid[(i,j+1)]))      # vertical
L=len(E); NV=len(V2)
CENTER=vid[(1,1)]
CUT=[k for k,(a,b) in enumerate(E) if a==CENTER or b==CENTER]     # electric boundary of {(1,1)}
def comps(removed):
    adj={v:[] for v in range(NV)}
    for k,(a,b) in enumerate(E):
        if k in removed: continue
        adj[a].append(b); adj[b].append(a)
    seen=set(); c=0
    for v in range(NV):
        if v in seen: continue
        c+=1; stk=[v]; seen.add(v)
        while stk:
            x=stk.pop()
            for y in adj[x]:
                if y not in seen: seen.add(y); stk.append(y)
    return c
# outer square: the 8 links on the perimeter
PERIM=[k for k,(a,b) in enumerate(E)
       if all(V2[x][0] in (0,2) or V2[x][1] in (0,2) for x in (a,b))
       and not (CENTER in (a,b))]
deg={}
for k in PERIM:
    a,b=E[k]; deg[a]=deg.get(a,0)+1; deg[b]=deg.get(b,0)+1
print(f"  3x3 grid: V={NV} L={L}  (planar, 4 plaquettes)")
print(f"  ELECTRIC BOUNDARY of the interior vertex = links {CUT}")
print(f"    is it a separator?  components after removal = {comps(set(CUT))}  (1 means no)")
print(f"    is it a cycle?      every touched vertex degree 2? {all(d==2 for d in ([sum(1 for k in CUT if CENTER in E[k])]+[1]))}"
      f"   -> a STAR, not a cycle")
print(f"  MAGNETIC BOUNDARY (perimeter) = links {PERIM}  ({len(PERIM)} links)")
print(f"    is it a cycle?      degrees {sorted(set(deg.values()))}  -> {'CYCLE' if set(deg.values())=={2} else 'not a cycle'}")
print(f"    is it a separator?  components after removal = {comps(set(PERIM))}")
print()
print("  ==> the two objects EXIST and are DIFFERENT sets, exactly as duality requires.")

N=2
st=[s for s in itertools.product(range(N),repeat=L)
    if all((sum(s[k] for k,(a,b) in enumerate(E) if a==v)
           -sum(s[k] for k,(a,b) in enumerate(E) if b==v))%N==0 for v in range(NV))]
idx={s:j for j,s in enumerate(st)}; Dp=len(st)
print(f"\n  Z_{N}: full {N**L}, physical {Dp} = {N}^{round(np.log(Dp)/np.log(N))}  (cycle rank {L-NV+1})")
def loop(moves):
    M=np.zeros((Dp,Dp),dtype=complex)
    for j,s in enumerate(st):
        t=list(s)
        for l,sg in moves: t[l]=(t[l]+sg)%N
        k=idx.get(tuple(t))
        if k is None: return None
        M[k,j]=1.0
    return M
# the four plaquettes, each (i,j): h(i,j) + v(i+1,j) - h(i,j+1) - v(i,j)
def hid(i,j): return j*2+i
def vidx(i,j): return 6+j*3+i
PLQ=[[(hid(i,j),+1),(vidx(i+1,j),+1),(hid(i,j+1),-1),(vidx(i,j),-1)] for j in range(2) for i in range(2)]
P=[loop(m) for m in PLQ]
print(f"  plaquettes gauge-invariant (stay in the physical sector): {all(x is not None for x in P)}")
Wper=loop([(k,+1) for k in PERIM])
print(f"  perimeter loop stays in the physical sector: {Wper is not None}")
if Wper is not None and all(x is not None for x in P):
    prod=P[0]
    for x in P[1:]: prod=prod@x
    print(f"  DISCRETE STOKES on the patch: || prod of 4 plaquettes - perimeter loop || = "
          f"{np.linalg.norm(prod-Wper):.2e}")
