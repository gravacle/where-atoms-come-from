#!/usr/bin/env python3
"""
MUTATION TEST: can O-10's instrument produce a FAIL / a 'match: NO'?
A control that cannot fail is not a control.
Four mutations, each with a KNOWN correct answer:
  M1 wrong d_known supplied           -> 'match: NO' expected
  M2 wrong k_expect supplied          -> SC2 and SC5 FAIL expected
  M3 a NON-COMMUTING 'stabiliser' set -> SC1/SC2/SC3 must break (Pc is not a projector
     onto a joint +1 eigenspace / is zero)
  M4 phase_vector deliberately broken -> SC1/SC3 FAIL expected
Also: M5 asks whether Part 2's 'balanced (tr R = 0)' check CAN fail for any stabiliser
code -- proved below that it cannot, so it is not a control.
"""
import sys, io, contextlib, numpy as np
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_O10_SCALING')
import o10_threshold as T

def cap(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            r = fn(*a, **k)
        except Exception as e:
            return buf.getvalue() + '\nEXCEPTION: %r' % (e,), None
    return buf.getvalue(), r

print('--- M1: tell it Steane has d=2 (it is 3). Expect match: NO ---')
out, r = cap(T.measure_threshold, '[[7,1,3]] Steane MUTATED d_known=2',
             ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ'], 7, 1, 2)
print([l for l in out.splitlines() if 'MEASURED THRESHOLD' in l][0])
print('  -> instrument CAN report a mismatch:', 'match: NO' in out)

print()
print('--- M2: tell it Steane has k=2 (it is 1). Expect SC2 and SC5 FAIL ---')
out, r = cap(T.measure_threshold, '[[7,1,3]] Steane MUTATED k=2',
             ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ'], 7, 2, 3)
for l in out.splitlines():
    if 'SC2' in l or 'SC5' in l or 'MEASURED THRESHOLD' in l: print(l)
print('  -> self-checks CAN fail:', 'FAIL' in out)

print()
print('--- M3: non-commuting generator set (XIIIIII, ZIIIIII). Expect breakage ---')
out, r = cap(T.measure_threshold, 'NONCOMMUTING (mutation)', ['XIIIIII','ZIIIIII'], 7, 5, 1)
for l in out.splitlines():
    if 'SC1' in l or 'SC2' in l or 'SC3' in l: print(l)
print('  -> detected:', 'FAIL' in out)

print()
print('--- M4: break phase_vector (drop the Z sign). Expect SC1/SC3 FAIL ---')
orig = T.phase_vector
T.phase_vector = lambda b, dim: np.ones(dim)
out, r = cap(T.measure_threshold, 'PHASE BUG (mutation)',
             ['IIIXXXX','IXXIIXX','XIXIXIX','IIIZZZZ','IZZIIZZ','ZIZIZIZ'], 7, 1, 3)
T.phase_vector = orig
for l in out.splitlines():
    if l.strip().startswith('SC1') or l.strip().startswith('SC2') or l.strip().startswith('SC3'):
        print(l)
print('  -> detected:', 'FAIL' in out or 'EXCEPTION' in out)

print()
print('--- M5: CAN Part 2\'s "balanced (tr R = 0)" check ever FAIL? ---')
print('  R = V^dag E V with E a HERMITIAN Pauli commuting with S, E not in <+-S>.')
print('  tr R = tr(E Pc) = (1/|S|) sum_{s in S} tr(E s).  tr of a Pauli string is 0 unless')
print('  the string is +-I; E s = +-I iff E = +-s iff E in <+-1, S>, which is EXCLUDED.')
print('  Therefore tr R = 0 identically, for EVERY weight-d witness of EVERY stabiliser code.')
print('  Likewise R^dag = V^dag E^dag V = R (E Hermitian by construction in pauli_matrix),')
print('  and R^2 = V^dag E Pc E V = V^dag E^2 Pc V = I  (E commutes with Pc, E^2 = I).')
print('  => "BIT clause (i) PASS on all 119 witnesses" and "balanced PASS" are THEOREMS with')
print('     NO FAILURE MODE.  They verify the numerics, they are NOT evidence.')
print()
print('  numeric confirmation that these are forced: scan ALL logicals of the [[5,1,3]] code')
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_O10_SCALING')
from o10_proof_check import pauli_matrix, isometry
from o10_threshold import parse_pauli, weight, code_projector
gens = [parse_pauli(s) for s in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']]
P = code_projector(gens,5); V = isometry(P,2)
basis,pivots=[],[]
for (ga,gb) in gens:
    cur=(ga<<5)|gb
    for p,bv in zip(pivots,basis):
        if cur>>p&1: cur^=bv
    if cur: basis.append(cur); pivots.append(cur.bit_length()-1)
def in_S(a,b):
    cur=(a<<5)|b
    for p,bv in zip(pivots,basis):
        if cur>>p&1: cur^=bv
    return cur==0
def commS(a,b): return all((bin(a&gb).count('1')+bin(ga&b).count('1'))%2==0 for (ga,gb) in gens)
mx=0.0; cnt=0
for a in range(32):
    for b in range(32):
        if (a,b)==(0,0): continue
        if commS(a,b) and not in_S(a,b):
            R=V.conj().T@pauli_matrix(a,b,5)@V
            mx=max(mx, abs(np.trace(R)), np.max(np.abs(R-R.conj().T)), np.max(np.abs(R@R-np.eye(2))))
            cnt+=1
print('  ALL %d logicals of [[5,1,3]] (every weight, not just weight-d): max deviation from'%cnt)
print('  (tr R = 0, R = R^dag, R^2 = I) = %.3e  -- forced at EVERY weight, so the check does'%mx)
print('  not single out weight d at all.')
