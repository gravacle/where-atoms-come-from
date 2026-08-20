"""O-29: A CARRIER WHERE CURVATURE CAN ACTUALLY APPEAR.

Every carrier so far has gauge group Z_2, whose holonomies are +-1 -- a SIGN, which cannot rotate a
frame. H-9 asks for a frame mismatch, so those nulls were measured where the effect has no room to
exist.

AND Z_N WOULD NOT HELP. On the minimal torus cell structure (one vertex, two edges, one face) the
plaquette holonomy is the GROUP COMMUTATOR g1 g2 g1^-1 g2^-1, which is trivial for ANY abelian group
-- U(1) and Z_N included. Only a NON-ABELIAN group gives non-trivial carrier curvature on a torus.

Smallest non-abelian group: S_3, order 6. Hilbert space |G|^edges = 36 on the minimal torus."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
# ---- S_3 as permutations of (0,1,2) ----
G=list(itertools.permutations(range(3)))
def mul(a,b): return tuple(a[b[i]] for i in range(3))
def inv(a):
    r=[0]*3
    for i,x in enumerate(a): r[x]=i
    return tuple(r)
e=(0,1,2)
idx={g:i for i,g in enumerate(G)}; n=len(G)
say("="*100); say("O-29   NON-ABELIAN CARRIER: S_3 ON THE MINIMAL TORUS"); say("="*100)
say(f"  group S_3, order {n}; minimal torus cell structure: 1 vertex, 2 edges, 1 face")
say(f"  Hilbert space |G|^edges = {n*n}")
nonab=[(a,b) for a in G for b in G if mul(a,b)!=mul(b,a)]
say(f"  non-commuting pairs: {len(nonab)} of {n*n}  -> the group commutator is NOT always trivial")
say("")
# ---- the carrier: basis |g1,g2>; plaquette holonomy = commutator ----
def comm(a,b): return mul(mul(a,b),mul(inv(a),inv(b)))
say("1.  CARRIER CURVATURE IS A GROUP ELEMENT, NOT A SIGN")
classes={}
for a,b in itertools.product(G,G): classes.setdefault(comm(a,b),[]).append((a,b))
say(f"  {'holonomy g1g2g1^-1g2^-1':<28}{'# configurations':>18}{'flat?':>8}")
for h,cfg in sorted(classes.items(), key=lambda kv:-len(kv[1])):
    say(f"  {str(h):<28}{len(cfg):>18}{('FLAT' if h==e else 'CURVED'):>8}")
say("")
say("2.  THE FLAT SECTOR IS THE RECORD SPACE")
flat=[(a,b) for a,b in itertools.product(G,G) if comm(a,b)==e]
say(f"  flat configurations (commuting pairs): {len(flat)} of {n*n}")
# gauge transformations act by simultaneous conjugation
def conj(h,p): return (mul(mul(h,p[0]),inv(h)), mul(mul(h,p[1]),inv(h)))
orb=[]; seen=set()
for p in flat:
    if p in seen: continue
    o={conj(h,p) for h in G}; seen|=o; orb.append(sorted(o))
say(f"  modulo gauge (simultaneous conjugation): {len(orb)} classes")
say(f"  -> the record space has {len(orb)} elements = Hom(pi_1(T^2), S_3)/conj")
say("")
say("3.  DOES CURVATURE CHANGE WHAT THE CARRIER CAN HOLD?")
cur=[(a,b) for a,b in itertools.product(G,G) if comm(a,b)!=e]
corb=[]; seen=set()
for p in cur:
    if p in seen: continue
    o={conj(h,p) for h in G}; seen|=o; corb.append(sorted(o))
say(f"  {'sector':<22}{'configurations':>16}{'gauge classes':>15}")
say(f"  {'FLAT (curvature = e)':<22}{len(flat):>16}{len(orb):>15}")
say(f"  {'CURVED':<22}{len(cur):>16}{len(corb):>15}")
say("")
say("4.  IS THE HOLONOMY PATH-DEPENDENT?  (H-9's actual question)")
say("    transport around the face one way vs the other: g1g2g1^-1g2^-1 against g2g1g2^-1g1^-1")
diff=0; tot=0
for a,b in itertools.product(G,G):
    h1=comm(a,b); h2=comm(b,a); tot+=1
    if h1!=h2: diff+=1
say(f"    configurations where the two orders DISAGREE: {diff} of {tot}")
say(f"    -> {'PATH-DEPENDENT: the two chains give different frame elements' if diff else 'path-independent'}")
say("")
say("  READ: on a Z_2 carrier every holonomy is +-1 and both orders always agree, so H-9's mismatch")
say("  cannot occur. Here it can. This is the first carrier in the program where a non-zero answer")
say("  to the criterion has room to exist.")
