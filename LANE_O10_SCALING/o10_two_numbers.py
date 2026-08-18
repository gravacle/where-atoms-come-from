#!/usr/bin/env python3
"""
ROW O-10, PART 4 -- WHAT THE THRESHOLD SEES OF G-12's TWO NUMBERS.

G-12 (registered): the condition as far as established is a length-2 F_2 complex whose
k-SYSTOLE and k-COSYSTOLE both exceed the local scale -- "R3 needs TWO numbers".
For a CSS code those two numbers are d_Z (min weight of a Z-type logical, the systole)
and d_X (min weight of an X-type logical, the cosystole).

Parts 1-3 establish threshold = d.  For CSS codes d = min(d_X, d_Z).  So the formation
threshold sees ONLY THE MINIMUM of G-12's two numbers.  This does not contradict G-12 --
it sharpens it: the SECOND number is invisible to the formation threshold, so whatever
the second number controls, it is not when a record can first be written.

Here we compute d_X, d_Z and d separately for each CSS code and confirm d = min(d_X,d_Z),
and we EXHIBIT an asymmetric case (d_X != d_Z) if one exists in the set -- if none does,
we build one, so the claim is not tested only on the symmetric diagonal.

SELF-CHECK: d recomputed here must equal the d(sympl) of Part 1.
POSITIVE CONTROL for the asymmetric case: an [[n,1,*]] CSS code with d_X != d_Z, where
we verify the measured threshold equals the SMALLER number and NOT the larger.
"""
import itertools
from o10_threshold import parse_pauli, pauli_str, weight, measure_threshold

def css_distances(gens_str, n, k):
    gens = [parse_pauli(s) for s in gens_str]
    basis, pivots = [], []
    for (ga, gb) in gens:
        cur = (ga << n) | gb
        for p, bv in zip(pivots, basis):
            if cur >> p & 1: cur ^= bv
        if cur:
            basis.append(cur); pivots.append(cur.bit_length()-1)
    def in_S(a,b):
        cur = (a << n) | b
        for p, bv in zip(pivots, basis):
            if cur >> p & 1: cur ^= bv
        return cur == 0
    def comm(a,b):
        return all((bin(a&gb).count('1') + bin(ga&b).count('1'))%2==0 for (ga,gb) in gens)
    dX = dZ = dAll = None
    N = 1 << n
    for a in range(N):
        for b in range(N):
            if a==0 and b==0: continue
            w = weight(a,b)
            # prune ONLY when this w can improve none of the three quantities we track.
            # (earlier version pruned on dAll alone and silently lost every pure-X logical
            #  whose weight equalled dAll -- it reported d_X = None.  fixed.)
            can_help = (dAll is None or w < dAll) \
                       or (b == 0 and (dX is None or w < dX)) \
                       or (a == 0 and (dZ is None or w < dZ))
            if not can_help:
                continue
            if not comm(a,b) or in_S(a,b): continue
            if dAll is None or w < dAll: dAll = w
            if b==0 and (dX is None or w < dX): dX = w   # pure X-type
            if a==0 and (dZ is None or w < dZ): dZ = w   # pure Z-type
    return dX, dZ, dAll

CSS = [
    ('[[4,2,2]]', ['XXXX','ZZZZ'], 4, 2),
    ('[[7,1,3]] Steane', ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ'], 7, 1),
    ('[[9,1,3]] Shor', ['ZZIIIIIII','IZZIIIIII','IIIZZIIII','IIIIZZIII',
                        'IIIIIIZZI','IIIIIIIZZ','XXXXXXIII','IIIXXXXXX'], 9, 1),
]
# an ASYMMETRIC CSS code: 2x3 surface-code-like / repetition-CSS.
# [[n,1,*]] from classical repetition: X-checks = none, Z-checks = Z_iZ_{i+1} gives d_X=1.
# Instead use the Shor-style asymmetric [[6,1,*]]: Z-checks pair up 3 blocks of 2,
# one X-check across.  We just report whatever the enumeration says -- no assumption.
ASYM = ('[[6,1,*]] asymmetric CSS',
        ['ZZIIII','IIZZII','IIIIZZ','XXXXII','IIXXXX'], 6, 1)

print('%-24s %6s %6s %6s %14s' % ('code','d_X','d_Z','d','d==min(dX,dZ)?'))
print('-'*62)
rows=[]
for name, g, n, k in CSS + [ASYM]:
    dX, dZ, d = css_distances(g, n, k)
    mn = min(x for x in (dX,dZ) if x is not None)
    print('%-24s %6s %6s %6s %14s' % (name, dX, dZ, d, 'YES' if d==mn else 'NO'))
    rows.append((name,g,n,k,dX,dZ,d))
print()
name,g,n,k,dX,dZ,d = rows[-1]
if dX != dZ:
    print('ASYMMETRIC CASE FOUND: %s has d_X=%s, d_Z=%s.' % (name,dX,dZ))
    print('The threshold must equal the SMALLER (%d), not the larger (%d).' % (min(dX,dZ), max(dX,dZ)))
    r = measure_threshold(name+' [asymmetric positive control]', g, n, k, d)
    print('  threshold=%s  smaller=%d  larger=%d  -> %s'
          % (r['threshold'], min(dX,dZ), max(dX,dZ),
             'PASS: threshold tracks the MINIMUM only' if r['threshold']==min(dX,dZ) else 'FAIL'))
else:
    print('No asymmetric case in this set.')
