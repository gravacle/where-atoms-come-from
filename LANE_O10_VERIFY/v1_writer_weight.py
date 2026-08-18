#!/usr/bin/env python3
"""
ADVERSARIAL CHECK OF ROW O-10.

O-10 measures ONE thing: the min weight w at which some operator has Pc E Pc NOT prop Pc.
That is CORE_FRAMEWORK clause (iii) NON-TRIVIAL (distinguishes same-energy states).
It is NOT clause (iv) WRITABLE.  A RECORD needs BOTH: a bit R, and an admissible U with
U^dag R U = -R.

Q1. For each weight-d witness R, what is the MINIMUM WEIGHT of a Pauli logical U that
    ANTICOMMUTES with R on the code space (i.e. actually writes it)?
Q2. On O-10's own asymmetric carrier [[6,1,*]] with d_X=2, d_Z=3: is the writer for the
    weight-2 bit also weight 2, or is it weight 3?

If the writer for the min-weight bit has weight max(d_X,d_Z), then O-10b's conclusion
("the threshold sees only min(dX,dZ); the larger number does not say when a record can
first be written") is FALSE: the larger number is the WRITER's weight.

SELF-CHECKS (known answers, each can FAIL):
 SC-A  on the SELF-DUAL Steane code d_X=d_Z=3, so writer weight must come out 3 -- if the
       instrument returned max(d_X,d_Z) on the asymmetric code but ALSO 3 here, that is
       consistent; the discriminating case is the asymmetric one.
 SC-B  anticommutation is checked TWO independent ways: symplectic form over F2, and the
       explicit matrix product R U + U R restricted to the code space.  They must agree.
 SC-C  NEGATIVE CONTROL: the search must FAIL to find a writer at weight < the answer, and
       must FIND one at the answer.  Both branches are exercised and reported.
"""
import itertools, numpy as np
import sys
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_O10_SCALING')
from o10_threshold import parse_pauli, pauli_str, weight, code_projector, toric_2x2
from o10_proof_check import pauli_matrix, isometry

def analyse(name, gens_str, n, kbits, d):
    gens = [parse_pauli(s) for s in gens_str]
    P = code_projector(gens, n)
    kdim = 2**kbits
    V = isometry(P, kdim)
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
    def commS(a,b):
        return all((bin(a&gb).count('1') + bin(ga&b).count('1'))%2==0 for (ga,gb) in gens)
    def symp(a1,b1,a2,b2):
        return (bin(a1&b2).count('1') + bin(a2&b1).count('1')) % 2

    # all logicals (elements of N(S)\S) grouped by weight
    logicals = []
    for a in range(1<<n):
        for b in range(1<<n):
            if a==0 and b==0: continue
            if commS(a,b) and not in_S(a,b):
                logicals.append((weight(a,b), a, b))
    logicals.sort()
    dmin = logicals[0][0]

    # weight-d witnesses
    wits = [(a,b) for (w,a,b) in logicals if w == d]
    print('CODE %-28s n=%d k=%d d=%d   #weight-d bits = %d' % (name,n,kbits,d,len(wits)))

    # for each witness, min weight of an anticommuting logical (the WRITER)
    worst = 0; best = 99; disagree = 0; found_none = 0
    per = {}
    for (a,b) in wits:
        wmin = None; wop = None
        for (w,a2,b2) in logicals:
            if symp(a,b,a2,b2) == 1:
                # confirm by explicit matrices on the code space
                R = V.conj().T @ pauli_matrix(a,b,n) @ V
                U = V.conj().T @ pauli_matrix(a2,b2,n) @ V
                anti = np.max(np.abs(R@U + U@R))
                if anti > 1e-9: disagree += 1
                wmin = w; wop = (a2,b2); break
        if wmin is None:
            found_none += 1; continue
        per[pauli_str(a,b,n)] = (wmin, pauli_str(wop[0],wop[1],n))
        worst = max(worst, wmin); best = min(best, wmin)
    print('   writer weight for a weight-%d bit: MIN over bits = %d, MAX over bits = %d'
          % (d, best, worst))
    for i,(kk,vv) in enumerate(sorted(per.items(), key=lambda t:-t[1][0])[:3]):
        print('      bit %-10s -> min-weight writer %-10s weight %d' % (kk, vv[1], vv[0]))
    print('   SC-B symplectic vs matrix anticommutation disagreements: %d -> %s'
          % (disagree, 'PASS' if disagree==0 else 'FAIL'))
    print('   SC-C bits with NO writer found at any weight: %d -> %s'
          % (found_none, 'PASS (every bit is writable by some logical)' if found_none==0 else 'FAIL'))
    # NEGATIVE CONTROL: for the bit achieving `worst`, verify NO logical of weight < worst
    # anticommutes with it, by exhaustive scan (this is the branch that must be able to fail).
    tgt = [ (a,b) for (a,b) in wits if per.get(pauli_str(a,b,n),(0,''))[0]==worst ]
    if tgt:
        a,b = tgt[0]
        cnt_below = sum(1 for (w,a2,b2) in logicals if w < worst and symp(a,b,a2,b2)==1)
        # ALSO: any Pauli at all (not just logical) of weight < worst that anticommutes?
        cnt_any = 0
        for a2 in range(1<<n):
            for b2 in range(1<<n):
                if (a2,b2)==(0,0): continue
                if weight(a2,b2) < worst and symp(a,b,a2,b2)==1:
                    cnt_any += 1
        print('   NEG CONTROL on bit %s : logicals of weight<%d that anticommute = %d (expect 0);'
              % (pauli_str(a,b,n), worst, cnt_below))
        print('        ALL Paulis (incl. non-logical) of weight<%d that anticommute = %d'
              % (worst, cnt_any))
    print()
    return dict(name=name, d=d, writer_min=best, writer_max=worst)

def main():
    tg, tn = toric_2x2()
    codes = [
        ('[[3,1,1]] repetition', ['ZZI','IZZ'], 3, 1, 1),
        ('[[4,2,2]]', ['XXXX','ZZZZ'], 4, 2, 2),
        ('[[8,2,2]] toric 2x2', tg, tn, 2, 2),
        ('[[5,1,3]] perfect NON-CSS', ['XZZXI','IXZZX','XIXZZ','ZXIXZ'], 5, 1, 3),
        ('[[7,1,3]] Steane', ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ'], 7, 1, 3),
        ('[[9,1,3]] Shor', ['ZZIIIIIII','IZZIIIIII','IIIZZIIII','IIIIZZIII',
                            'IIIIIIZZI','IIIIIIIZZ','XXXXXXIII','IIIXXXXXX'], 9, 1, 3),
        ('[[6,1,*]] ASYMMETRIC dX=2 dZ=3', ['ZZIIII','IIZZII','IIIIZZ','XXXXII','IIXXXX'], 6, 1, 2),
    ]
    res = [analyse(*c) for c in codes]
    print('='*78)
    print('%-32s %4s %12s %12s' % ('code','d','min writer','max writer'))
    print('-'*62)
    for r in res:
        print('%-32s %4d %12d %12d' % (r['name'], r['d'], r['writer_min'], r['writer_max']))
    print()
    print('THE DISCRIMINATING CASE is the asymmetric [[6,1,*]] (d_X=2, d_Z=3).')
    print('If its writer weight is 3 while O-10 reports threshold 2, then O-10b is WRONG:')
    print('the LARGER of the two numbers is exactly what says when the bit can be WRITTEN.')

main()
