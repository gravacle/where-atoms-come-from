# W-14 leg C — the ONE state rule that is not a stipulation, and what it buys.
# Every ready state in this corpus was CHOSEN. The exception is SENSE U, the uniform state, where
# pi is determined by the carrier's own class SIZES and nothing is picked. Under it lambda becomes
# a function of the carrier alone -- so the carrier does enter the rate. The question is HOW MUCH
# of the carrier enters, and S4's own published table already answers it.
import numpy as np
def m_jensen(p,n=1<<22):
    a,b,c,d=p; t=2*np.pi*np.arange(n)/n; ct=np.cos(t)
    A=np.sqrt(np.maximum(a*a+b*b+2*a*b*ct,0)); B=np.sqrt(np.maximum(c*c+d*d+2*c*d*ct,0))
    return float(np.log(np.maximum(A,B)+1e-300).mean())

# S4:511-515 topology and S4:574-583 class multisets + published SENSE U lambda, quoted verbatim
ROWS=[
 ("B3  horn torus (octahedron, poles ident.)", dict(V=5,chi=1,b1=1,b2=1), {"01":2,"10":2,"11":1}, -0.756573586),
 ("B1  K1 (the pinch, as handed)",             dict(V=5,chi=0,b1=1,b2=0), {"01":2,"10":2,"11":1}, -0.756573586),
 ("B2  K1, both triangles filled",             dict(V=5,chi=1,b1=0,b2=0), {"01":2,"10":2,"11":1}, -0.756573586),
 ("B1s K1, every edge subdivided",             dict(V=11,chi=0,b1=1,b2=0),{"01":5,"10":5,"11":1}, -0.724759919),
 ("B1q K1-bridged + spectator vertex",         dict(V=7,chi=0,b1=1,b2=0), {"00":1,"01":3,"10":3}, -0.741029583),
 ("B0b ring torus 3x3 grid, loops meet",       dict(V=9,chi=0,b1=2,b2=1), {"00":4,"01":1,"10":2,"11":2}, -0.810930216),
 ("B4  spindle (two spheres, two glue pts)",   dict(V=6,chi=2,b1=1,b2=2), {"00":1,"01":1,"10":1,"11":3}, -0.693147181),
]
print("== C1  UNDER THE UNIFORM STATE, lambda IS A FUNCTION OF THE CARRIER. RECOMPUTED. ==")
print(f"  {'carrier':<42} {'V':>2} {'chi':>4} {'b1':>3} {'b2':>3}  {'pi (uniform)':<28} {'lambda':>14} {'S4 published':>14}")
out=[]
for name,topo,mult,pub in ROWS:
    V=topo["V"]; pi=np.array([mult.get(k,0)/V for k in ("00","10","01","11")])
    lam=m_jensen(pi)
    out.append((name,topo,tuple(sorted(mult.values())),lam))
    print(f"  {name:<42} {V:>2} {topo['chi']:>4} {topo['b1']:>3} {topo['b2']:>3}  "
          f"{str(np.round(pi,4)):<28} {lam:>14.9f} {pub:>14.9f}")
print()
print("== C2  AND IT SEES THE CARRIER ONLY THROUGH THE MULTISET OF CLASS SIZES ==")
print("  S4's own table contains the exhibit and never remarks on it. B3, B1 and B2 share the")
print("  class multiset {01:2, 10:2, 11:1} and differ in EVERY topological invariant on the page:")
for name,topo,ms,lam in out[:3]:
    print(f"    {name:<42} chi={topo['chi']} b1={topo['b1']} b2={topo['b2']}   lambda = {lam:.12f}")
lams=[l for _,_,_,l in out[:3]]
print(f"  spread across the three = {max(lams)-min(lams):.3e}")
print()
print("  a horn torus, K1, and K1 with both triangles filled -- chi = 1, 0, 1; b1 = 1, 1, 0;")
print("  b2 = 1, 0, 0 -- give the IDENTICAL rate, because their class-size multisets coincide.")
print()
print("  ==> The uniform state is the only non-stipulated state rule, and under it the carrier")
print("      enters the rate ONLY through the MULTISET OF CLASS SIZES. That is the whole of the")
print("      channel. Nothing about chi, b1, b2, the pinch, the faces or the 2-cells survives it.")
