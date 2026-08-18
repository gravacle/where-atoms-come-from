#!/usr/bin/env python3
"""
Is 'weight >= d' SUFFICIENT?  Part 3's printed VERDICT says F-13's "coupling of
weight >= d" is EXACT.  Exact = necessary AND sufficient.  Test sufficiency.
NAMED COUNTEREXAMPLES wanted: operators of weight >= d with residual EXACTLY 0.
POSITIVE CONTROL: the same routine on a support that DOES carry a logical must fire.
"""
import itertools, numpy as np, sys
sys.path.insert(0,'/Users/bgm/MB Work/where-atoms-come-from/LANE_O10_SCALING')
from o10_threshold import parse_pauli, pauli_str, weight, code_projector, toric_2x2
from o10_proof_check import pauli_matrix
from o10_arbitrary_coupling import embed, resid
RNG = np.random.default_rng(7)

def probe(name, gens_str, n, kb, d):
    gens=[parse_pauli(s) for s in gens_str]; P=code_projector(gens,n); kd=2**kb
    dead=[]; live=[]
    for sup in itertools.combinations(range(n), d):
        M = RNG.normal(size=(1<<d,1<<d)) + 1j*RNG.normal(size=(1<<d,1<<d))
        M = (M+M.conj().T)/2
        r = resid(P, embed(M, sup, n), kd)
        (live if r>1e-9 else dead).append((sup, r))
    print('%-28s d=%d : supports of size d = %d | FIRE %d | SILENT %d'
          % (name, d, len(dead)+len(live), len(live), len(dead)))
    if dead:
        s,r = dead[0]
        print('    COUNTEREXAMPLE: a random dense HERMITIAN coupling of weight %d on qubits %s'
              % (d, list(s)))
        print('    gives residual %.3e -- weight >= d but NO record formed.' % r)
    if live:
        s,r = live[0]
        print('    POSITIVE CONTROL: same routine on qubits %s gives residual %.6f' % (list(s), r))
    # also: an arbitrarily HIGH weight operator that is silent -- a stabiliser element
    return len(dead), len(live)

tg,tn = toric_2x2()
probe('[[7,1,3]] Steane', ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ'],7,1,3)
probe('[[9,1,3]] Shor', ['ZZIIIIIII','IZZIIIIII','IIIZZIIII','IIIIZZIII',
                         'IIIIIIZZI','IIIIIIIZZ','XXXXXXIII','IIIXXXXXX'],9,1,3)
probe('[[8,2,2]] toric', tg, tn, 2, 2)
probe('[[5,1,3]] perfect', ['XZZXI','IXZZX','XIXZZ','ZXIXZ'],5,1,3)

print()
print('MAXIMAL-WEIGHT SILENT COUPLINGS (weight n >> d):')
for nm,g,n,kb in [('Steane',['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ'],7,1),
                  ('Shor',['ZZIIIIIII','IZZIIIIII','IIIZZIIII','IIIIZZIII',
                           'IIIIIIZZI','IIIIIIIZZ','XXXXXXIII','IIIXXXXXX'],9,1)]:
    gens=[parse_pauli(s) for s in g]; P=code_projector(gens,n)
    # product of all stabiliser generators = a high-weight element of S: acts as lambda*Pc
    a=b=0
    for (ga,gb) in gens: a^=ga; b^=gb
    E=pauli_matrix(a,b,n); A=P@E@P; lam=np.trace(A)/(2**kb)
    print('  %-8s stabiliser product %s  weight %d  residual %.3e -> SILENT'
          % (nm, pauli_str(a,b,n), weight(a,b), np.max(np.abs(A-lam*P))))
print()
print('CONCLUSION: weight >= d is NECESSARY (proved) but NOT SUFFICIENT.')
print('Part 3\'s printed verdict "F-13\'s coupling of weight >= d is EXACT" is an OVERSTATEMENT.')
