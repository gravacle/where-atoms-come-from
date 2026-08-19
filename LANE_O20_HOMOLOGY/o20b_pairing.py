"""O-20b: chi is NOT a function of the homology class alone. Non-boundaries SPLIT --
[4,6] and [0,2] give 0.908, [3,7] and [1,5] give 0.000, all four non-boundaries.

dim H_1 = 2, so there are THREE non-trivial classes. The discriminator must be which class,
relative to the record's. The candidate: the INTERSECTION PAIRING of the coupling's cycle with
the WRITER's support -- |c ∩ supp(Xbar)| mod 2 -- which is exactly the second half of the pair
(H_1, <,>) the program registered as Gamma. H_1 ALONE does not decide it; H_1 WITH ITS FORM does.

This needs no time evolution: the pairing is GF(2) arithmetic, checked against F-15's criterion."""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
supp = lambda v: [l for l in range(L) if (v>>(L-1-l))&1]
sX, sZ = supp(Xc1), supp(Zc1)
say("="*98); say("O-20b  WHAT DECIDES chi: THE CLASS, OR THE PAIRING WITH THE RECORD?"); say("="*98)
say(f"  record  Zbar supported on links {sZ}")
say(f"  writer  Xbar supported on links {sX}")
say("")
say(f"  {'coupling':<14}{'|c ∩ supp Xbar|':>18}{'pairing mod 2':>15}{'||[Z_c, Xbar]||':>18}{'chi measured':>15}")
MEASURED = {(4,6):0.90811968, (3,7):0.00000000, (1,5):0.00000000, (0,2):0.90811968}
ok=True
for links,chi in MEASURED.items():
    A = op({l:Z for l in links}, L)
    inter = len(set(links) & set(sX)); par = inter % 2
    comm = np.linalg.norm(A@Xbar - Xbar@A)
    pred = (par == 1)
    got  = (chi > 1e-10)
    ok &= (pred == got)
    say(f"  {str(list(links)):<14}{inter:>18}{par:>15}{comm:>18.4f}{chi:>15.8f}"
        f"   {'ok' if pred==got else 'MISMATCH'}")
say("")
say(f"  chi > 0  IFF  the coupling's cycle pairs ODDLY with the writer : {'CONFIRMED' if ok else 'FALSIFIED'}")
say("")
say("  READ: H_1 alone does not decide it -- all four couplings above are non-boundaries.")
say("        What decides it is the INTERSECTION PAIRING with the record's conjugate.")
say("        That is the pair (H_1, <,>), not H_1. The form is doing the work.")
