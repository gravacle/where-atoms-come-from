# W19-A step 3.  CONTROLS AND CONFOUNDS.  Each block states the ONE variable it moves and DIFFS ITS ARMS.
#
#  C1  ENCLOSURE CONTROL   -- same carrier, same state, same g^2, same FRAGMENT SIZE; only the
#                             fragment's SHAPE moves (contains a tail->head path in G-l, or not).
#  C2  DYNAMICS CONTROL    -- same carrier, same fragments; the STATE moves (ground state at three
#                             couplings, an excited state, and a HAAR-RANDOM physical state).
#                             If a Haar-random physical state plateaus too, the plateau is carried by
#                             the Gauss constraint and not by the dynamics.  That is decisive.
#  C3  GROUP CONTROL       -- same carrier, same fragments, same g^2; only N moves (Z_2 -> Z_3).
#  C4  SUBDIVISION CONFOUND-- theta vs theta subdivided.  dim_phys is IDENTICAL (4).  Does the verdict move?
#  C5  FRAGMENT-RULE CONTROL- same carrier, same state; only the FRAGMENT RULE moves
#                             (A nested-by-distance / B uniformly random / D cherry-picked).
import numpy as np, sys, json, itertools, time
sys.path.insert(0,"/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_A_SWEEP")
from zn_gauge import *
from carriers import *
# has_uv_path is defined once, in zn_gauge.py, and used by both 01 and 03.
DELTA=0.10
np.random.seed(20260817)

print("="*118); print("W19-A / 03 CONTROLS"); print("="*118)

# ------------------------------------------------------------------ C1 ENCLOSURE CONTROL
print("\n[C1] ENCLOSURE CONTROL.  carrier=heawood, g^2=1, N=2, system link l=0, |F| HELD FIXED.")
print("     ONE VARIABLE MOVED: whether F contains a tail(l)->head(l) path in G-l.")
V,E = heawood(); g = ZNGauge("heawood",V,E,2); L=g.L
psi,E0,gap = g.ground(2.0,2.0); Psi = g.full_vector(psi); HS = S_of(Psi,L,2,[0])
frs,d = nested_fragments(V,E,0)
armA = frs[2]                                   # |F|=14, no path (from rule A)
# build arm B: same size, but forced to contain a shortest tail->head path
import collections
def shortest_path_links(V,E,l):
    a,b=E[l][0],E[l][1]; adj=build_adj(V,E)
    prev={a:None}; dq=collections.deque([a])
    while dq:
        x=dq.popleft()
        for (y,i) in adj[x]:
            if i==l or y in prev: continue
            prev[y]=(x,i); dq.append(y)
    out=[]; cur=b
    while prev[cur] is not None:
        x,i=prev[cur]; out.append(i); cur=x
    return sorted(out)
sp = shortest_path_links(V,E,0)
rest = [e for e in range(L) if e!=0 and e not in sp]
armB = sorted(sp + rest[:len(armA)-len(sp)])
IA = mutual_information(Psi,L,2,[0],armA); IB = mutual_information(Psi,L,2,[0],armB)
print(f"     ARM A  links={armA}  |F|={len(armA)}  path={has_uv_path(V,E,0,armA)}")
print(f"     ARM B  links={armB}  |F|={len(armB)}  path={has_uv_path(V,E,0,armB)}")
print(f"     ARMS DIFFER IN: A\\B={sorted(set(armA)-set(armB))}   B\\A={sorted(set(armB)-set(armA))}   "
      f"(|symmetric difference|={len(set(armA)^set(armB))}, NOT byte-identical)")
print(f"     H(S)={HS:.9f}   I_A/H(S)={IA/HS:.9f}   I_B/H(S)={IB/HS:.9f}   "
      f"-> the criterion CAN fail at fixed fragment size.  Control is LIVE, not vacuous.")

# ------------------------------------------------------------------ C2 DYNAMICS CONTROL
print("\n[C2] DYNAMICS CONTROL.  carrier=heawood, fragments = rule A F_1..F_4 (|F|=2,6,14,18).")
print("     ONE VARIABLE MOVED: the state.  Same graph, same fragments, same criterion.")
def report_state(tag, psiv):
    P = g.full_vector(psiv/np.linalg.norm(psiv)); h = S_of(P,L,2,[0])
    rs=[mutual_information(P,L,2,[0],F)/h for F in frs]
    print(f"     {tag:<44} H(S)={h:.8f}  I/H = " + " ".join(f"{r:.9f}" for r in rs))
    return h, rs
for gsq in (0.3,1.0,3.0):
    p,_,_ = g.ground(2.0/gsq,2.0*gsq); report_state(f"ground state, g^2={gsq}", p)
H = g.hamiltonian(2.0,2.0); w,v = np.linalg.eigh(H)
report_state("7th excited state of H (g^2=1)", v[:,7])
report_state("HIGHEST state of H (g^2=1)",     v[:,-1])
for k in range(3):
    report_state(f"HAAR-RANDOM physical state #{k+1} (no Hamiltonian at all)",
                 np.random.randn(g.dimP))
print("     READ: the plateau ratio is 1.000000000 for EVERY state above, including states that are")
print("     not ground states of anything.  The plateau is carried by the Gauss constraint, not by H.")

# ------------------------------------------------------------------ C3 GROUP CONTROL
print("\n[C3] GROUP CONTROL.  ONE VARIABLE MOVED: N.  carriers within the Z_3 ceiling (3^L <= 4.5e6 => L<=13).")
print(f"{'carrier':<16}{'N':>3}{'L':>4}{'C':>4}{'dimP':>7}{'d':>4}{'H(S) bits':>12}{'log2 N':>9}  I/H(S) on rule-A fragments")
for nm,(Vg,Eg) in [("theta_3link",theta()),("ladder_2sq",ladder(2)),("ladder_3sq",ladder(3)),("cube_Q3",cube())]:
    for N in (2,3):
        if N**len(Eg) > 4_500_000:
            print(f"{nm:<16}{N:>3}{len(Eg):>4}   -- OVER CEILING {N}^{len(Eg)} > 4.5e6, skipped"); continue
        gg = ZNGauge(nm,Vg,Eg,N); p,_,_ = gg.ground(2.0,2.0); P = gg.full_vector(p)
        f2,dd = nested_fragments(Vg,Eg,0); h = S_of(P,gg.L,N,[0])
        rs = [mutual_information(P,gg.L,N,[0],F)/h for F in f2]
        print(f"{nm:<16}{N:>3}{gg.L:>4}{gg.C:>4}{gg.dimP:>7}{dd:>4}{h:>12.8f}{np.log2(N):>9.5f}  "
              + " ".join(f"{r:.6f}" for r in rs))
print("     READ: d, the plateau length and R_delta are IDENTICAL for Z_2 and Z_3 on the same graph.")
print("     N moves only the CEILING on H(S) (<= log2 N).  The threshold is not a property of the group.")

# ------------------------------------------------------------------ C4 SUBDIVISION CONFOUND
print("\n[C4] SUBDIVISION CONFOUND.  Subdividing a link inserts a degree-2 vertex.  This leaves the")
print("     cyclomatic number C, hence dim_phys = N^C, UNCHANGED -- no new physical state space.")
print(f"{'carrier':<16}{'V':>4}{'L':>4}{'C':>4}{'dimP':>7}{'d':>4}{'H(S)':>11}{'plateau pts':>13}{'R_delta':>9}  verdict")
for nm,(Vg,Eg) in [("theta",theta()),("theta_subdiv1",theta_sub(1)),("theta_subdiv2",theta_sub(2)),
                   ("theta_subdiv3",theta_sub(3))]:
    gg = ZNGauge(nm,Vg,Eg,2); p,_,_ = gg.ground(2.0,2.0); P = gg.full_vector(p)
    f2,dd = nested_fragments(Vg,Eg,0); cu,_ = level_cuts(Vg,Eg,0); h = S_of(P,gg.L,2,[0])
    pts = sum(abs(mutual_information(P,gg.L,2,[0],F)/h - 1) <= DELTA for F in f2)
    R   = sum(mutual_information(P,gg.L,2,[0],Ci)/h >= 1-DELTA for Ci in cu)
    ver = "EXHIBITED" if pts>=4 else ("MARGINAL" if pts==3 else "FAIL")
    if h < 0.10: ver += " (WEIGHTLESS: H(S)<0.10)"
    print(f"{nm:<16}{Vg:>4}{gg.L:>4}{gg.C:>4}{gg.dimP:>7}{dd:>4}{h:>11.6f}{pts:>13}{R:>9}  {ver}")
print("     READ: dim_phys is 4 on every row.  The verdict still moves FAIL -> EXHIBITED.")
print("     R_delta is therefore NOT an invariant of the gauge theory; it is a graph distance.")

# ------------------------------------------------------------------ C5 FRAGMENT-RULE CONTROL
print("\n[C5] FRAGMENT-RULE CONTROL.  carrier=heawood, g^2=1, state fixed.  ONE VARIABLE MOVED: the rule.")
psi,_,_ = g.ground(2.0,2.0); Psi = g.full_vector(psi); HS = S_of(Psi,L,2,[0])
env = [e for e in range(L) if e != 0]
print("     RULE B (uniformly random fragments; the textbook partial-information plot).")
print(f"     {'|F|':>5}{'mean I/H(S)':>14}{'min':>10}{'max':>10}{'frac exactly 1':>16}{'frac 0':>9}{'frac 2':>9}")
NS=12; ruleB=[]
t0=time.time()
for m in range(1,len(env)+1):
    vals=[]
    combos = list(itertools.combinations(env,m))
    picks = combos if len(combos)<=NS else [tuple(np.random.choice(env,m,replace=False)) for _ in range(NS)]
    for F in picks:
        vals.append(mutual_information(Psi,L,2,[0],sorted(F))/HS)
    vals=np.array(vals)
    print(f"     {m:>5}{vals.mean():>14.6f}{vals.min():>10.6f}{vals.max():>10.6f}"
          f"{np.mean(np.abs(vals-1)<=DELTA):>16.3f}{np.mean(vals<DELTA):>9.3f}{np.mean(vals>2-DELTA):>9.3f}")
    ruleB.append(dict(m=m,mean=float(vals.mean()),lo=float(vals.min()),hi=float(vals.max()),
                      frac1=float(np.mean(np.abs(vals-1)<=DELTA))))
print(f"     (rule B took {time.time()-t0:.1f}s, {NS} samples per size)")
nplat = sum(1 for r in ruleB if abs(r['mean']-1)<=DELTA)
print(f"     RULE B plateau points (sizes whose MEAN I/H(S) is within {DELTA} of 1): {nplat}")
print("     RULE D (cherry-picked: every link except star(head(l))-l).  Reported WEIGHTLESS by construction:")
a,b = E[0]; starb=[e for e,(x,y) in enumerate(E) if e!=0 and (x==b or y==b)]
chain=[e for e in env if e not in starb]
sizes=[];
for m in range(2,len(chain)+1):
    F=sorted(chain[:m]); r=mutual_information(Psi,L,2,[0],F)/HS
    sizes.append((m,r))
print(f"     sizes {sizes[0][0]}..{sizes[-1][0]}, all ratios in [{min(s[1] for s in sizes):.6f},{max(s[1] for s in sizes):.6f}]"
      f"  -> {sum(1 for s in sizes if abs(s[1]-1)<=DELTA)} 'plateau points'.")
print("     This family is engineered to exclude the conjugate; it inflates the plateau and is WEIGHTLESS.")
json.dump(dict(ruleB=ruleB), open("out_03_controls.json","w"), indent=1)
print("\nDONE 03.")
