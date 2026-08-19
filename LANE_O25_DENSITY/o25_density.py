"""O-25: does record DENSITY force curvature?

Before sweeping density, a structural check. In a stabiliser code the logical operators obey a
SYMPLECTIC algebra: conjugate pairs (X_i, Z_i) anticommute, and every other pair COMMUTES. In
particular X_i and X_j always commute, whatever their supports and however densely packed.

If that holds, then overlap is irrelevant among Pauli writers and density cannot force anything --
the curvature O-24 found came from X_0 * S_1, where S_1 is a PHASE gate and not a Pauli.

Tested on every k>=2 stabiliser code available, with all logicals COMPUTED by symplectic
Gram-Schmidt (never nominated), across a range of overlap."""
import sys, itertools, numpy as np
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import symplectic_logicals, xz_to_matrix
def say(*a): print(*a); sys.stdout.flush()
def s2xz(s,n): return [1 if c in 'XY' else 0 for c in s]+[1 if c in 'ZY' else 0 for c in s]
def supp(v,n): return {i for i in range(n) if v[i] or v[n+i]}
say("="*100); say("O-25   CAN DENSITY FORCE CURVATURE AMONG PAULI WRITERS?"); say("="*100)
CODES=[("[[4,2,2]]", 4, ['XXXX','ZZZZ']),
       ("[[8,3,2]]", 8, ['XXXXXXXX','ZZZZZZZZ','ZZZZIIII','ZZIIZZII','ZIZIZIZI']),
       ("[[6,4,2]]", 6, ['XXXXXX','ZZZZZZ']),
       ("toric 2x2 [[8,2,2]]", 8, ['XXIXIXII','IXXIXIXI','ZZZZIIII','IIZZZZII'])]
say(f"  {'code':<22}{'n':>3}{'records':>9}{'writer overlap':>16}{'max ||[X_i,X_j]||':>20}{'verdict':>10}")
allflat=True
for nm,n,gens in CODES:
    try:
        pairs=symplectic_logicals([s2xz(g,n) for g in gens],n)
    except Exception as e:
        say(f"  {nm:<22}{n:>3}   construction failed: {e}"); continue
    if len(pairs)<2:
        say(f"  {nm:<22}{n:>3}{len(pairs):>9}   fewer than two records -- skipped"); continue
    Xs=[xz_to_matrix(b,n) for a,b in pairs]
    Sm=[xz_to_matrix(s2xz(g,n),n) for g in gens]; H=-sum(Sm)
    if not all(np.linalg.norm(A@H-H@A)<1e-9 for A in Xs):
        say(f"  {nm:<22}{n:>3}   VERIFICATION FAILED -- skipped"); continue
    ov=max(len(supp(pairs[i][1],n)&supp(pairs[j][1],n))
           for i,j in itertools.combinations(range(len(pairs)),2))
    mc=max(float(np.linalg.norm(Xs[i]@Xs[j]-Xs[j]@Xs[i]))
           for i,j in itertools.combinations(range(len(pairs)),2))
    allflat &= (mc<1e-9)
    say(f"  {nm:<22}{n:>3}{len(pairs):>9}{ov:>16}{mc:>20.3e}{('FLAT' if mc<1e-9 else 'CURVED'):>10}")
say("")
say(f"  every code flat: {allflat}")
say("")
say("  READ")
say("    Writers overlap on up to several sites and STILL commute exactly. In a stabiliser code the")
say("    logical algebra is symplectic: conjugate pairs anticommute, every other pair commutes, and")
say("    support has nothing to do with it. So DENSITY CANNOT FORCE CURVATURE among Pauli writers,")
say("    and O-25's premise is wrong.")
say("")
say("    O-24's curvature came from X_0 * S_1 with S_1 a PHASE gate -- not a Pauli, and not a")
say("    minimal writer. Curvature requires leaving the stabiliser formalism, which is exactly what")
say("    minimality forbids.")
