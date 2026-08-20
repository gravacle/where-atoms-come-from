"""
VERIFY-1.  IS EXHIBIT 4 (the 2^n diagonal-algebra sweep, s5_grading.py) MEASURING ANYTHING?

The lane's classifier, copied verbatim from s5_grading.py:
      nonconst  = any(~const_per_block)          <- clause (iii)
      balanced  = all(Tr(P_E Z_S) == 0)          <- clause (iv)
      is_rec    = balanced & nonconst
      all_const = all(const_per_block)           == NOT nonconst
      carries   = all_const & varies
so   is_rec & carries = (balanced & nonconst) & ((NOT nonconst) & varies) = EMPTY, identically,
for ANY Hamiltonian, any couplings, any n, with or without clause (iv), with or without the
flip symmetry.  This script tests that reading by REMOVING clause (iv) from the record axis and
by destroying every symmetry of H, and asking whether the both-box ever fills.
"""
import itertools, random
import numpy as np

def hadamard(n):
    Hm = np.ones((1,1), dtype=np.int8)
    for _ in range(n): Hm = np.block([[Hm,Hm],[Hm,-Hm]])
    return Hm

def spin_table(n):
    xs = np.arange(2**n)
    return np.array([1-2*((xs>>i)&1) for i in range(n)], dtype=np.int8).T

def boxes(E, n, use_clause_iv=True):
    """Exactly the lane's classifier.  use_clause_iv=False drops `balanced` from the record axis."""
    order = np.argsort(E, kind="stable"); Es = E[order]
    starts = np.flatnonzero(np.r_[True, Es[1:] != Es[:-1]])
    Hm = hadamard(n)[:, order].astype(np.int32)
    bsum = np.add.reduceat(Hm, starts, axis=1)
    bmax = np.maximum.reduceat(Hm, starts, axis=1)
    bmin = np.minimum.reduceat(Hm, starts, axis=1)
    const_per_block = (bmax == bmin)
    balanced = np.all(bsum == 0, axis=1)
    nonconst = np.any(~const_per_block, axis=1)
    is_rec = (balanced & nonconst) if use_clause_iv else nonconst
    all_const = np.all(const_per_block, axis=1)
    varies = (bmax.max(axis=1) != bmax.min(axis=1))
    carries = all_const & varies
    return (int(np.sum(is_rec & carries)), int(np.sum(is_rec & ~carries)),
            int(np.sum(carries & ~is_rec)), int(np.sum(~is_rec & ~carries)))

def E_ising(J, n, h=0.0, seed=None):
    S = spin_table(n); E = np.zeros(2**n, dtype=np.int64)
    for i in range(len(J)): E += J[i]*S[:,i].astype(np.int64)*S[:,i+1].astype(np.int64)
    if h: E += int(h)*S[:,0].astype(np.int64)
    return E

print("="*100)
print("VERIFY-1   THE BOTH-BOX IS EMPTY BY DEFINITION, NOT BY MEASUREMENT")
print("="*100)
print()
print("  The two axes are complementary predicates on the SAME quantity:")
print("     RECORD  requires  clause (iii): NON-constant on SOME eigenspace")
print("     CARRIER requires               CONSTANT on EVERY eigenspace")
print("  so no observable of any kind can be in both, whatever H is.")
print()
print("  TEST 1.  Re-run the lane's own classifier with clause (iv) DELETED from the record axis.")
print("           If clause (iv) were 'the separator', deleting it should let the both-box fill.")
print()
print(f"  {'H':<34} {'n':>3} {'both(lane)':>11} {'rec':>6} {'en':>5} {'none':>6} "
      f"{'both WITHOUT clause(iv)':>24} {'rec':>6}")
rows=[]
rnd = random.Random(11)
for n in (5,7,9,11,13):
    for name, J, h in (("ising uniform J=1",[1]*(n-1),0),
                       ("ising superinc J=2^i",[2**i for i in range(n-1)],0),
                       ("ising randpos + FIELD h=3",[rnd.randrange(1,61) for _ in range(n-1)],3)):
        E = E_ising(J,n,h)
        a = boxes(E,n,True); b = boxes(E,n,False)
        rows.append((a[0],b[0]))
        print(f"  {name:<34} {n:>3} {a[0]:>11} {a[1]:>6} {a[2]:>5} {a[3]:>6} {b[0]:>24} {b[1]:>6}")
print()
print("  TEST 2.  A DIAGONAL H WITH NO SYMMETRY AT ALL -- energies drawn i.i.d. at random, so the")
print("           spectrum is non-degenerate, no global flip exists, and clause (iv) FAILS for")
print("           every observable.  #records must be 0.  Does the both-box fill?")
print()
print(f"  {'H':<34} {'n':>3} {'#balanced (clause iv)':>22} {'both(lane)':>11} "
      f"{'both WITHOUT clause(iv)':>24} {'rec w/o (iv)':>13}")
for n in (5,7,9,11):
    E = np.array(random.Random(n).sample(range(10**7), 2**n), dtype=np.int64)
    order = np.argsort(E); Es=E[order]
    starts = np.flatnonzero(np.r_[True, Es[1:]!=Es[:-1]])
    Hm = hadamard(n)[:,order].astype(np.int32)
    nbal = int(np.sum(np.all(np.add.reduceat(Hm,starts,axis=1)==0,axis=1)))
    a = boxes(E,n,True); b = boxes(E,n,False)
    print(f"  {'random non-degenerate diagonal H':<34} {n:>3} {nbal:>22} {a[0]:>11} "
          f"{b[0]:>24} {b[1]:>13}")
print()
mx_lane = max(r[0] for r in rows); mx_noiv = max(r[1] for r in rows)
print(f"  READ: both-box with clause (iv) = {mx_lane} (max over rows).")
print(f"        both-box with clause (iv) DELETED = {mx_noiv} (max over rows).")
print("        Identical. Clause (iv) does NO work in emptying that box; clause (iii) alone")
print("        empties it, because 'carrier' was defined as the negation of clause (iii).")
print("        The 'BREAKING CONTROL' (field h) cannot break a box that is empty by definition,")
print("        and the h=3 rows in the lane's own output show it stays 0 exactly as it must.")
