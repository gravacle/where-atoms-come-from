"""VERIFY check 1 (independent): toric d_W by DIRECT low-weight enumeration, L=2,3.
No affine solver, no coset span: enumerate ALL x-masks of weight <= 2L+1, keep those that
commute with every plaquette (cycles), classify by crossing parity with the two cut systems,
record min weight per class.  Also: NO admissible writer of nonzero class below the claimed
minimum exists (exhaustive over the weight shell), and the claimed minima are achieved."""
import itertools, sys

def run(L):
    n = 2*L*L
    def h(i,j): return (i%L)*L + (j%L)
    def v(i,j): return L*L + (i%L)*L + (j%L)
    plaqs = []
    for i in range(L):
        for j in range(L):
            m = 0
            for e in (h(i,j), h(i+1,j), v(i,j), v(i,j+1)):
                m |= 1<<e
            plaqs.append(m)
    # cut systems: row loops (h-edges of row i), column loops (v-edges of col j)
    rowloop0 = 0
    for j in range(L): rowloop0 |= 1<<h(0,j)
    colloop0 = 0
    for i in range(L): colloop0 |= 1<<v(i,0)
    best = {}
    maxw = 2*L + 1
    for w in range(0, maxw+1):
        for comb in itertools.combinations(range(n), w):
            m = 0
            for e in comb: m |= 1<<e
            if any(bin(m & p).count("1") % 2 for p in plaqs):
                continue  # not a cycle -> not admissible (anticommutes with a plaquette)
            a = bin(m & colloop0).count("1") % 2   # winding in horizontal direction
            b = bin(m & rowloop0).count("1") % 2
            key = (a,b)
            if key not in best:
                best[key] = w
    return best

ok = True
for L in (2,3):
    best = run(L)
    want = {(0,0):0, (1,0):L, (0,1):L, (1,1):2*L}
    match = (best == want)
    ok &= match
    print(("PASS" if match else "FAIL"), "independent enumeration L=%d: min weight per class %s (claimed %s)" % (L, best, want))
print("V1 OVERALL:", "PASS" if ok else "FAIL")
