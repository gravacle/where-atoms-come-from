# WHAT DISTINGUISHES "INTERIOR" FROM "SURFACE" ON A GRAPH?
# In the continuum a surface is codimension-1: it SEPARATES. There is no geometry here, so the only
# intrinsic meaning available is: A SURFACE IS A SEPARATOR -- a set of links whose removal
# disconnects. Test the wheel's rim and spokes against that, and against W-24's requirement that a
# boundary be a CYCLE (so a flux can thread it).
import itertools
def wheel(n): return [(0,k+1) for k in range(n)]+[(k+1,(k+1)%n+1) for k in range(n)], n+1
def components(E,V,removed):
    adj={v:[] for v in range(V)}
    for i,(a,b) in enumerate(E):
        if i in removed: continue
        adj[a].append(b); adj[b].append(a)
    seen=set(); c=0
    for v in range(V):
        if v in seen: continue
        c+=1; st=[v]; seen.add(v)
        while st:
            x=st.pop()
            for y in adj[x]:
                if y not in seen: seen.add(y); st.append(y)
    return c
def is_cycle(E,links):
    """do these links form a single closed cycle? every touched vertex has degree exactly 2"""
    deg={}
    for i in links:
        a,b=E[i]; deg[a]=deg.get(a,0)+1; deg[b]=deg.get(b,0)+1
    return all(d==2 for d in deg.values())

for n in (4,5,6):
    E,V=wheel(n); SPK=set(range(n)); RIM=set(range(n,2*n))
    print(f"  wheel n={n}:  V={V} L={len(E)}   components with nothing removed: {components(E,V,set())}")
    print(f"    remove the RIM    -> components {components(E,V,RIM)}   is a cycle: {is_cycle(E,RIM)}")
    print(f"    remove the SPOKES -> components {components(E,V,SPK)}   is a cycle: {is_cycle(E,SPK)}")
print()
print("  A SURFACE MUST BE BOTH: a SEPARATOR (removing it disconnects) and a CYCLE (a flux can")
print("  thread it). The rim is a cycle and separates NOTHING. The spokes separate and are NOT a")
print("  cycle. THE WHEEL HAS NEITHER OBJECT. W-24 required only the cycle half and I built to that.")
print()
print("  And note what this says about the labels: the edge boundary of the region {hub} is the")
print("  SPOKE SET, not the rim. I have been calling the rim the boundary because it is drawn as a")
print("  circle. That is a picture, not a property of the graph.")
