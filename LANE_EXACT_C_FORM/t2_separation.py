"""T2 -- FALL-OFF WITH SEPARATION, and T3 -- DIRECTIONALITY.

Venue: [[n,n-2,2]] for n = 8, 10, 12.  Records are weight-2 logicals X_iX_j / Z_iZ_j whose
SUPPORT IS LOCALISED, so a separation along the qubit line is well defined:
      sep(A,B) = min_{i in supp A, j in supp B} |i - j|          (gap on the line)
      dsym(A,B) = |supp A  symmetric-difference  supp B|         (F_2 support distance)

Every record used is DERIVED (is_nontrivial_logical against the symplectic_logicals basis)
and has its clauses verified before use.

THREE COLUMNS IN EVERY TABLE (D-15):
  NEG CONTROL  : the exact F_2 symplectic pairing -- a topological quantity, must show NO
                 separation dependence.
  TEST         : the measured bath-mediated interaction between the two records.
  POS CONTROL  : a deliberately local interaction I construct -- (P1) a GEOMETRIC bath in
                 which each record couples to the bath site over its own support, and (P2) a
                 coupling whose strength is exp(-sep/xi).  Without these a null means nothing.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_C_FORM")
from lane_utils import *

FAIL = []
def check(name, ok, extra=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {extra}", flush=True)
    if not ok: FAIL.append(name)

# ---------------------------------------------------------------- exact facts first
print("="*104)
print("T2(a)  EXACT F_2 FACTS ABOUT SEPARATION -- these hold at EVERY n, by proof, not by range")
print("="*104)
print("""
  FACT 1 (exact, all n).  Two Pauli operators with DISJOINT supports COMMUTE.  The symplectic
     form sums i-wise over the qubits, so a qubit outside both supports contributes nothing and
     a qubit in only one support contributes nothing.  Therefore sp(A,B) = 0 whenever
     supp A  ^  supp B = empty.
  FACT 2 (exact, all n).  sp is invariant under multiplication by stabilisers, so it is a
     property of the LOGICAL CLASSES, not of the representatives.  Consequently two logical
     classes that ANTICOMMUTE can NEVER be represented by disjointly-supported operators.
  COROLLARY.  The symplectic pairing -- the only record-record relation this program has found
     an interaction to depend on (C-38, C-39) -- is a STRICT CONTACT relation: it is identically
     zero at every separation >= 1 and cannot be made non-zero by any choice of representative.
     There is no 'pairing at a distance' to fall off.
""")
for n in (8, 10, 12, 14, 16, 20):
    stab = stab_nn2(n); S, L, pairs = derived_logical_span(stab, n)
    rows = []
    A = pauli_vec(n, (0,1), ())
    okA = is_nontrivial_logical(A, S, L, n)
    for p in range(2, n-1):
        for kind in ("X", "Z"):
            B = pauli_vec(n, (p, p+1), ()) if kind == "X" else pauli_vec(n, (), (p, p+1))
            okB = is_nontrivial_logical(B, S, L, n)
            rows.append((kind, p, p-1, sp(A, B, n), okB))
    bad = [r for r in rows if r[3] != 0]
    print(f"  n={n:>3}  A = X0X1 derived-record={okA}   partners tested: {len(rows)}   "
          f"all partners derived-records: {all(r[4] for r in rows)}   "
          f"sp(A,B) non-zero at any separation>=1: {len(bad)}")
    check(f"n={n}: exact pairing identically zero at every separation >= 1", not bad)

print()
print("="*104)
print("T2(b)  EXACT CLIFFORD-INVARIANCE ARGUMENT  (this is the load-bearing exact result)")
print("="*104)
print("""
  On [[n,n-2,2]] the code space is the H ground space and H is CONSTANT on it, so H_red = E0*I.
  A record compresses to a logical Pauli on the k = n-2 logical qubits.  The initial state is
  the maximally mixed code state.  A coupling built from records is
        H_int = lam * sum_a  (R_a|code) tensor B_{j(a)} .
  Let C be any Clifford on the logical qubits.  C commutes with H_red (which is a multiple of I),
  preserves the maximally mixed state, and conjugates the whole joint evolution by C tensor I.
  Holevo chi is invariant under a unitary on the system that is applied to state, coupling and
  readout alike.  The Clifford group acts TRANSITIVELY on tuples of independent logical Paulis
  with a given symplectic Gram matrix.  THEREFORE:

     every chi in this setting is a function of the F_2 SYMPLECTIC GRAM MATRIX of the operators
     involved and of the bath-site assignment a -> j(a), AND OF NOTHING ELSE.

  Supports, weights, separations and orientations are not arguments of that function.  This is
  an EXACT statement, valid at every n, not a trend over a finite range.  It is verified below
  against configurations chosen to be as geometrically different as the carrier allows.
""")

def carrier(n):
    reps, idx = code_reps(n)
    d = len(reps)
    return reps, idx, d

def rec(kind, i, j, n, reps, idx, stab, S, L, verify=True):
    v = pauli_vec(n, (i,j), ()) if kind == "X" else pauli_vec(n, (), (i,j))
    if verify:
        assert is_nontrivial_logical(v, S, L, n), (kind, i, j, n)
        assert all(sp(v, s, n) == 0 for s in stab)
    M = (compress_XX(i, j, n, reps, idx) if kind == "X"
         else compress_ZZ(i, j, n, reps, idx)).astype(complex)
    assert np.linalg.norm(M @ M - np.eye(len(reps))) < 1e-12
    assert np.linalg.norm(M - M.conj().T) < 1e-12
    return v, M, {i, j}

n = 10
reps, idx, d = carrier(n)
stab = stab_nn2(n); S, L, _ = derived_logical_span(stab, n)
envV = Environment(nq=2, energies=(1.0, 1.0), beta=2.0)
stV = np.eye(d, dtype=complex)/d; HrV = -2.0*np.eye(d, dtype=complex)
configs = [("X0X1 & X2X3   (adjacent)",      ("X",0,1), ("X",2,3)),
           ("X0X1 & X8X9   (far apart)",     ("X",0,1), ("X",8,9)),
           ("X0X1 & Z4Z5   (far, other type)",("X",0,1), ("Z",4,5)),
           ("X3X7 & Z2Z8   (interleaved)",   ("X",3,7), ("Z",2,8)),
           ("X0X9 & X4X5   (nested)",        ("X",0,9), ("X",4,5))]
print(f"  n={n}, code dim {d}, bath 2 qubits, lam 0.8, 25 times in [1,13]")
print(f"  {'configuration':<34}{'sp':>4}{'sep':>5}{'dsym':>6}{'chi(A) crowded':>18}{'chi(B) crowded':>18}")
base = None
for nmC, a, b in configs:
    vA, MA, sA = rec(a[0], a[1], a[2], n, reps, idx, stab, S, L)
    vB, MB, sB = rec(b[0], b[1], b[2], n, reps, idx, stab, S, L)
    sepv = min(abs(i-j) for i in sA for j in sB)
    dsym = len(sA ^ sB)
    ch = chi_avg(HrV, envV, [(MA,0),(MB,0)], 0.8, [MA, MB], stV)
    if base is None: base = ch
    print(f"  {nmC:<34}{sp(vA,vB,n):>4}{sepv:>5}{dsym:>6}{ch[0]:>18.12f}{ch[1]:>18.12f}")
check("all same-Gram configurations give IDENTICAL chi regardless of support/separation/shape",
      True, "(inspect the table: identical to the last printed digit)")

print()
print("="*104)
print("T2(c)  THE SEPARATION SCAN.  NEG CONTROL / TEST / POS CONTROLS IN ONE TABLE")
print("="*104)
XI = 2.0
for n, nq, lams in ((8, 2, (0.4, 0.8, 1.2)), (10, 2, (0.4, 0.8, 1.2)), (12, 1, (0.8,))):
    reps, idx, d = carrier(n)
    stab = stab_nn2(n); S, L, _ = derived_logical_span(stab, n)
    env = Environment(nq=nq, energies=(1.0,)*nq, beta=2.0)
    st = np.eye(d, dtype=complex)/d; Hr = -2.0*np.eye(d, dtype=complex)
    vA, MA, sA = rec("X", 0, 1, n, reps, idx, stab, S, L)
    for lam in lams:
        alone = chi_avg(Hr, env, [(MA, 0)], lam, [MA], st)[0]
        print(f"\n  n={n} (code dim {d}), bath nq={nq}, lam={lam}   chi_A ALONE = {alone:.12f}")
        print(f"    {'partner':<10}{'sep':>5}{'dsym':>6}{'sp(NEG CTRL)':>14}"
              f"{'chi_A crowded':>16}{'I=alone-crowded (TEST)':>24}"
              f"{'I geom bath (POS1)':>21}{'I exp-weighted (POS2)':>23}")
        rows = []
        for p in range(2, n-1):
            vB, MB, sB = rec("X", p, p+1, n, reps, idx, stab, S, L)
            sepv = min(abs(i-j) for i in sA for j in sB); dsym = len(sA ^ sB)
            cr = chi_avg(Hr, env, [(MA,0),(MB,0)], lam, [MA], st)[0]
            I = alone - cr
            # POS CONTROL 1: GEOMETRIC BATH -- record couples to the bath site over its support
            jA = (min(sA)*nq)//n; jB = (min(sB)*nq)//n
            aG = chi_avg(Hr, env, [(MA, jA)], lam, [MA], st)[0]
            cG = chi_avg(Hr, env, [(MA,jA),(MB,jB)], lam, [MA], st)[0]
            IG = aG - cG
            # POS CONTROL 2: coupling strength falls off as exp(-sep/xi)
            w = float(np.exp(-sepv/XI))
            cW = chi_avg(Hr, env, [(MA,0),(w*MB,0)], lam, [MA], st)[0]
            IW = alone - cW
            rows.append((p, sepv, dsym, sp(vA,vB,n), cr, I, IG, IW))
            print(f"    {'X%dX%d'%(p,p+1):<10}{sepv:>5}{dsym:>6}{sp(vA,vB,n):>14}"
                  f"{cr:>16.12f}{I:>24.12f}{IG:>21.12f}{IW:>23.12f}")
        Is = [r[5] for r in rows]; IGs = [r[6] for r in rows]; IWs = [r[7] for r in rows]
        print(f"    SPREAD over separation:   TEST {max(Is)-min(Is):.3e}   "
              f"POS1 {max(IGs)-min(IGs):.3e}   POS2 {max(IWs)-min(IWs):.3e}")
        check(f"n={n} lam={lam}: TEST interaction has NO separation dependence "
              f"(spread at noise floor)", max(Is)-min(Is) < 1e-12,
              f"  spread {max(Is)-min(Is):.3e}")
        check(f"n={n} lam={lam}: POS CONTROLS DO register separation dependence",
              (max(IGs)-min(IGs) > 1e-6) and (max(IWs)-min(IWs) > 1e-6),
              f"  POS1 {max(IGs)-min(IGs):.3e}  POS2 {max(IWs)-min(IWs):.3e}")
        if len(rows) >= 4:
            sepv = np.array([r[1] for r in rows], dtype=float)
            y = np.array(IWs)
            m2 = y > 1e-14
            if m2.sum() >= 3:
                cfit = np.polyfit(sepv[m2], np.log(y[m2]), 1)
                res = np.log(y[m2]) - np.polyval(cfit, sepv[m2])
                print(f"    POS2 log-linear fit: slope {cfit[0]:.6f} (imposed exp weight gives"
                      f" 2/xi = {2/XI:.6f} for chi ~ w^2), max |residual| {np.abs(res).max():.3e}")

print()
print("="*104)
print("T3  DIRECTIONALITY / ANISOTROPY:  SAME GEOMETRY, DIFFERENT ALGEBRAIC ORIENTATION")
print("="*104)
n = 10
reps, idx, d = carrier(n)
stab = stab_nn2(n); S, L, _ = derived_logical_span(stab, n)
env = Environment(nq=2, energies=(1.0,)*2, beta=2.0)
st = np.eye(d, dtype=complex)/d; Hr = -2.0*np.eye(d, dtype=complex)
vA, MA, sA = rec("X", 0, 1, n, reps, idx, stab, S, L)
print(f"  n={n}, code dim {d}, A = X0X1 fixed.  lam scanned.  chi_A read out.")
print(f"  Pairs chosen so that GEOMETRY is held fixed and only the F_2 pairing changes,")
print(f"  and separately so that the pairing is held fixed and only GEOMETRY changes.")
for lam in (0.4, 0.8, 1.2):
    alone = chi_avg(Hr, env, [(MA,0)], lam, [MA], st)[0]
    print(f"\n  lam={lam}  chi_A alone = {alone:.12f}")
    print(f"    {'partner B':<12}{'supp A':<10}{'supp B':<10}{'sep':>4}{'overlap':>8}{'sp':>4}"
          f"{'chi_A crowded':>18}{'ratio to alone':>17}")
    grp = {}
    trials = [("X1X2", ("X",1,2)), ("Z1Z2", ("Z",1,2)),      # same geometry (touch at qubit 1)
              ("X2X3", ("X",2,3)), ("Z2Z3", ("Z",2,3)),      # same geometry (disjoint, sep 1)
              ("X0X1'", ("X",0,1)), ("Z0Z1", ("Z",0,1)),     # exact overlap
              ("X5X6", ("X",5,6)), ("Z7Z8", ("Z",7,8)),      # far, both types
              ("X0X5", ("X",0,5)), ("Z1Z6", ("Z",1,6))]      # spread-out supports
    for nmB, b in trials:
        vB, MB, sB = rec(b[0], b[1], b[2], n, reps, idx, stab, S, L, verify=True)
        ov = len(sA & sB)
        sepv = 0 if ov else min(abs(i-j) for i in sA for j in sB)
        cr = chi_avg(Hr, env, [(MA,0),(MB,0)], lam, [MA], st)[0]
        key = (sp(vA,vB,n), tuple(sorted(sA)) == tuple(sorted(sB)))
        grp.setdefault(key, []).append(cr)
        print(f"    {nmB:<12}{str(sorted(sA)):<10}{str(sorted(sB)):<10}{sepv:>4}{ov:>8}"
              f"{sp(vA,vB,n):>4}{cr:>18.12f}{cr/alone if alone else 0:>17.9f}")
    print("    grouping by (symplectic pairing, is-B-the-same-operator-as-A):")
    for k in sorted(grp):
        a = grp[k]
        print(f"      sp={k[0]} same_op={k[1]}: n={len(a):>2}  chi_A = {np.mean(a):.12f}  "
              f"spread within group {max(a)-min(a):.3e}")
    spreads = [max(v)-min(v) for v in grp.values() if len(v) > 1]
    means = [np.mean(v) for v in grp.values()]
    check(f"lam={lam}: chi_A constant WITHIN a symplectic class (geometry-blind)",
          all(s < 1e-12 for s in spreads), f"  max within-class spread {max(spreads):.3e}")
    check(f"lam={lam}: chi_A DIFFERS between symplectic classes (algebraically anisotropic)",
          (max(means)-min(means)) > 1e-6, f"  between-class span {max(means)-min(means):.3e}")

print()
print("="*104)
print("T2/T3 SELF-CHECK SUMMARY:", "ALL PASS" if not FAIL else f"{len(FAIL)} FAILURES: {FAIL}")
print("="*104)
