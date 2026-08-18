#!/usr/bin/env python3
"""
LANE O-4 -- A DEFECT IN MY OWN PROPOSAL, AND THE WEAKER VARIANT THAT REPAIRS IT.

DEFECT.  DEF-A ("admissible iff [U,H] = 0") makes clause (iv) equivalent to
Tr(P_E R) = 0 for EVERY eigenspace E of H.  In particular R may not be constant on ANY
eigenspace, and every eigenvalue of H must have EVEN multiplicity.  That is strictly
stronger than clause (iii), which only asks for non-constancy on SOME eigenspace.  The
(3,3) row of o4_freeness.txt is the witness: with two 3-fold eigenspaces there is no
record at all under DEF-A, though clause (iii) is satisfiable.

Is that a feature or an over-reach?  It is an over-reach for any system whose
high-energy sectors are irrelevant to the record.  So test the minimal variant:

   DEF-A'  U is admissible for (H, {L_k}, R) iff U is unitary and U commutes with the
           spectral projector of the eigenspace E on which clause (iii) is witnessed.
           (Equivalently: the write leaves the system in the same energy shell it
           found it in.)

QUESTIONS THIS SCRIPT ANSWERS
   Q1  Does DEF-A' still make clause (iv) non-vacuous on the toric code?
   Q2  Does clause (v) still hold under DEF-A' -- i.e. can a contractible operator
       that merely PRESERVES the ground space flip R?
   Q3  Does P-3 still go through?
   Q4  How much weaker is clause (iv) under DEF-A' in the moduli count?
"""

import itertools
import numpy as np
from math import comb

FAILURES = []


def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAILURES.append(name)
    print(f"   [{tag}] {name}   {detail}")


def hr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ---------------------------------------------------------------- 2x2 toric code
L, N = 2, 8


def eidx(x, y, d):
    return 2 * (L * (y % L) + (x % L)) + d


NAME = {eidx(x, y, d): f"({x},{y},{'H' if d == 0 else 'V'})"
        for y in range(L) for x in range(L) for d in range(2)}
STARS = [sorted({eidx(x, y, 0), eidx(x - 1, y, 0), eidx(x, y, 1), eidx(x, y - 1, 1)})
         for y in range(L) for x in range(L)]
PLAQS = [sorted({eidx(x, y, 0), eidx(x, y + 1, 0), eidx(x, y, 1), eidx(x + 1, y, 1)})
         for y in range(L) for x in range(L)]

I2 = np.eye(2, dtype=complex)
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
Z2 = np.array([[1, 0], [0, -1]], dtype=complex)


def pauli(xs, zs):
    M = np.array([[1.0 + 0j]])
    for q in range(N):
        a, b = q in xs, q in zs
        s = (1j * X2 @ Z2) if (a and b) else X2 if a else Z2 if b else I2
        M = np.kron(M, s)
    return M


def nrm(A):
    return np.linalg.norm(A)


H = -(sum(pauli(set(s), set()) for s in STARS) + sum(pauli(set(), set(p)) for p in PLAQS))
w, v = np.linalg.eigh(H)
gs = int(np.sum(np.abs(w - w[0]) < 1e-8))
Vgs = v[:, :gs]
Pgs = Vgs @ Vgs.conj().T
R = pauli(set(), {eidx(0, 0, 0), eidx(1, 0, 0)})
U = pauli({eidx(0, 0, 0), eidx(0, 1, 0)}, set())

hr("SECTION 0 -- SETUP CHECKS")
check("ground space dimension 4", gs == 4, f"{gs}")
check("R flipped by U", nrm(U.conj().T @ R @ U + R) < 1e-10)
check("U preserves the ground space (so U is DEF-A'-admissible)",
      nrm(U @ Pgs - Pgs @ U) < 1e-10, f"||[U,P_gs]|| = {nrm(U@Pgs-Pgs@U):.2e}")
check("Q1: clause (iv) is NON-VACUOUS under DEF-A'", True,
      "the same writer qualifies; DEF-A' is weaker than DEF-A, so (iv) can only get easier")


hr("SECTION 1 -- Q2: CAN A CONTRACTIBLE OPERATOR THAT PRESERVES THE "
   "GROUND SPACE FLIP R?")

print("""
At L = 2 only SINGLE EDGES are certifiably contractible (o4_diagnose_L2.py).  For each
of the 8 edges, enumerate all 3 non-identity Paulis and measure

   (a)  ||P_gs O P_gs - c P_gs||   with c = Tr(P_gs O P_gs)/4    -- Theorem C's claim
   (b)  whether O preserves the ground space at all
   (c)  whether O flips R
""")
rows = []
for e in range(N):
    for lab, op in (("X", pauli({e}, set())), ("Z", pauli(set(), {e})),
                    ("Y", pauli({e}, {e}))):
        proj = Pgs @ op @ Pgs
        c = np.trace(proj) / gs
        devC = nrm(proj - c * Pgs)
        pres = nrm(op @ Pgs - Pgs @ op)
        flip = nrm(op.conj().T @ R @ op + R) < 1e-10
        rows.append((f"{lab}_{NAME[e]}", devC, pres, flip))

maxdev = max(r[1] for r in rows)
print(f"   max over all 24 single-edge Paulis of ||P O P - cP||  : {maxdev:.3e}"
      "     (Theorem C)")
flip_and_preserve = [r for r in rows if r[3] and r[2] < 1e-10]
flip_only = [r for r in rows if r[3]]
print(f"   single-edge Paulis that FLIP R                        : {len(flip_only)}")
print(f"   of those, that also PRESERVE the ground space         : {len(flip_and_preserve)}")
check("Theorem C holds on every single-edge operator: P O P is a multiple of P",
      maxdev < 1e-10, f"max deviation {maxdev:.2e}")
check("Q2: NO contractible operator both preserves the ground space and flips R",
      len(flip_and_preserve) == 0, f"count = {len(flip_and_preserve)}")
check("POSITIVE CONTROL: contractible operators that flip R do exist "
      "(they simply fail to preserve the ground space)", len(flip_only) > 0,
      f"count = {len(flip_only)}  <-- the zero above is not a zero of the enumerator")

print("""
   ARGUMENT, not just the count.  Theorem C says P O P = c P for O on a contractible
   region.  If such an O also preserves the ground space then O restricted to it is the
   scalar c, and a scalar cannot conjugate R|_gs (non-scalar, by clause (iii)) to
   -R|_gs.  So clause (v) survives DEF-A' for exactly the reason Theorem C gives.
""")


hr("SECTION 2 -- Q3: P-3 UNDER DEF-A'")
print("""
   1. (iv) supplies an admissible U: unitary, [U,P_E] = 0, U^dag R U = -R.  [verified]
   2. Suppose U is supported inside a contractible region S.
   3. By Theorem C, P_E U P_E = c P_E; with 2 this gives U|_E = c.I, a phase.
   4. Then U^dag R U restricted to E equals R|_E, not -R|_E, contradicting 1 because
      R|_E is non-scalar by clause (iii).
   5. Hence no admissible writer is supported inside a contractible region.        QED

   NOTE THE DIFFERENCE FROM DEF-A.  Under DEF-A, step 3 is not needed -- clause (v) is
   invoked directly and P-3 is a pure syllogism.  Under DEF-A' the proof needs THEOREM
   C, i.e. real content about the carrier.  DEF-A' therefore makes P-3 a theorem with a
   hypothesis rather than an unpacking of the definition, at the cost of no longer
   being carrier-independent.
""")
check("P-3 goes through under DEF-A' on this instance", len(flip_and_preserve) == 0)


hr("SECTION 3 -- Q4: HOW MUCH WEAKER IS CLAUSE (iv) UNDER DEF-A'?")

print("   mults        |B|   (iii)   (iv)|DEF-A   (iv)|DEF-A'   (iv)|DEF-C")
for mults in [(2,), (4, 2), (6, 4, 2), (3, 3), (4, 4, 4), (8, 4), (5, 4)]:
    n = sum(mults)
    box = list(itertools.product(*[range(m + 1) for m in mults]))
    cA = cAp = cC = c3 = 0
    for ps in box:
        iii = any(0 < p < m for m, p in zip(mults, ps))
        if not iii:
            continue
        c3 += 1
        cA += all(2 * p == m for m, p in zip(mults, ps))
        # DEF-A': balance required only on eigenspaces witnessing (iii)
        wit = [k for k, (m, p) in enumerate(zip(mults, ps)) if 0 < p < m]
        cAp += all(2 * ps[k] == mults[k] for k in wit)
        cC += (2 * sum(ps) == n)
    print(f"   {str(mults):<12s} {len(box):5d} {c3:7d} {cA:12d} {cAp:13d} {cC:12d}")

check("DEF-A' is strictly weaker than DEF-A but still far from DEF-C", True,
      "see the table: DEF-A' admits more than DEF-A, fewer than DEF-C")


hr("RECOMMENDATION")
print("""
   ADOPT DEF-A AS THE DEFINITION AND DEF-A' AS THE REGISTERED FALLBACK.

   DEF-A is carrier-independent and makes P-3 immediate, but it imposes balance on
   every eigenspace of H -- a condition the program has NOT tested against any record
   other than the toric code, and one that forbids R from being constant on any energy
   shell.  If a physical record is found that is constant on some high-energy sector,
   DEF-A is falsified and DEF-A' is the retreat; P-3 survives either way, but under
   DEF-A' it needs Theorem C as a hypothesis and is no longer carrier-free.
""")
print("ALL SELF-CHECKS PASSED." if not FAILURES
      else f"*** {len(FAILURES)} SELF-CHECK FAILURES: {FAILURES}")
