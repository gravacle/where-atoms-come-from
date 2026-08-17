#!/usr/bin/env python3
"""
LANE W-17 / ROUTE R3 -- F4 PRESUPPOSITION.

The sentence the question assumes, in one sentence:

   "There is one object called 'the formation functional', and 'carrier-independent'
    names one property of it."

Tested here as a claim in its own right, independently of the A/B choice.

Two failures are measured:
  (1) 'the formation functional' names TWO objects with OPPOSITE scope verdicts
      in W-10's own table -- the FORMATION CONDITION (THREE_CLASS_SCOPED) and
      the RATE lambda (CARRIER_INDEPENDENT and CONVENTION_SCOPED).
  (2) 'carrier-independent' names TWO inequivalent predicates:
        FORM-independence  : the formula carries no carrier symbol.
        VALUE-independence : the value does not change with the carrier.
      Reading A is a VALUE claim ("the physics does not depend on the complex").
      Reading B is a FORM claim ("only pi CAN enter").
      They are not two readings of one predicate.

Class multisets quoted verbatim:
  S4_THE_MEASUREMENT_V001.md:574-575  (B0a, B0b)
  W10_SCOPE_TABLE_V001.md:130-133     (B4-SQUARE, B4-TRIPENT)
"""
import numpy as np

def mahler(p00, p10, p01, p11, n=4096):
    """m(p00 + p10 x + p01 y + p11 xy) by 2-torus quadrature."""
    t = 2*np.pi*(np.arange(n) + 0.5)/n
    X, Y = np.meshgrid(np.exp(1j*t), np.exp(1j*t), indexing='ij')
    P = p00 + p10*X + p01*Y + p11*X*Y
    A = np.abs(P)
    A = np.where(A < 1e-300, 1e-300, A)
    return float(np.log(A).mean())

OUT = []
def P(s=""):
    print(s); OUT.append(str(s))

P("="*78)
P("F4  PRESUPPOSITION -- VALUE-INDEPENDENCE OF THE RATE, ACROSS CARRIERS")
P("="*78)
P()
P("  carrier / designation           class multiset            lambda        published")
rows = [
    ("K1 as handed  (3-class)",  (0.0, 2/5, 2/5, 1/5), "{10:2,01:2,11:1}/5", None),
    ("B0a ring torus, disjoint",  (2/9, 4/9, 3/9, 0.0), "{00:2,10:4,01:3}/9", -0.747659833),
    ("B0b ring torus, meeting",   (4/9, 2/9, 1/9, 2/9), "{00:4,10:2,01:1,11:2}/9", -0.810930216),
    ("B4-SQUARE  (S4's B4)",      (1/6, 1/6, 1/6, 3/6), "{00:1,10:1,01:1,11:3}/6", np.log(1/2)),
    ("B4-TRIPENT (W-10's B4)",    (1/6, 1/6, 2/6, 2/6), "{00:1,10:1,01:2,11:2}/6", np.log(1/3)),
]
vals = []
for name, p, ms, pub in rows:
    lam = mahler(*p)
    vals.append(lam)
    pubs = f"{pub:.9f}" if pub is not None else "     --      "
    P(f"  {name:30s} {ms:26s} {lam:12.9f}  {pubs}")
P()
P(f"  SPREAD of lambda across the corpus's own carriers = {max(vals)-min(vals):.9f}")
P(f"  ratio max/min |lambda|                            = {max(np.abs(vals))/min(np.abs(vals)):.6f}")
P()
P("  >>> The RATE is FORM-independent (no carrier symbol in the formula) and")
P("  >>> massively VALUE-dependent (spread 0.40 nats across five carriers the")
P("  >>> corpus itself publishes; a factor 1.59 on |lambda|).")
P("  >>> Reading A asserts VALUE-independence ('the physics does not depend on")
P("  >>> the complex') and is FALSE on the corpus's own table.")
P("  >>> Reading B asserts FORM-independence and is a one-line THEOREM, not a")
P("  >>> reading: fibre-wise multiplication acts on |x_v|^2 alone.")
P()
P("  Zero-variable check on the two B4 rows (W-10 N-1): SAME published")
P("  parameters V E F chi b0 b1 b2 gauge inv curv flat, DIFFERENT lambda:")
P(f"    B4-SQUARE  lambda = {vals[3]:.9f}     B4-TRIPENT lambda = {vals[4]:.9f}")
P(f"    |difference| = {abs(vals[3]-vals[4]):.9f}  -- the carrier is not even")
P("    determined by its own published invariants, so 'carrier-independent'")
P("    has no well-defined argument to be independent OF.")
P()
P("="*78)
P("F4  PRESUPPOSITION -- IS THE STIPULATION THE ONE THE QUESTION NAMES?")
P("="*78)
P()
P("  The question names TRANSPORT as the stipulated thing.  F2 measured")
P("  ||T^L - M_gamma||_F = 0: at the circuit clock the two named conventions are")
P("  the same operator.  The stipulation that does the work is therefore not")
P("  'which transport' but 'what counts as one tick' -- the CLOCK.")
P()
P("  W-10 already had the correct name, at W10_SCOPE_TABLE_V001.md:110-115:")
P("    \"the clock. W-01's 'circuit count is carrier-supplied discrete time' is")
P("     LOOP-LENGTH-SCOPED ... on B0b (lengths 4 and 3) T_F^12 = M_dF^3 while")
P("     T_C^12 = M_c^4 -- the two branches are at different circuit counts")
P("     forever, at every positive edge count.\"")
P()
P("  >>> PRESUPPOSITION_FALSE on two counts, both checkable pre-cutoff.")
P()

with open("w17_presupposition.OUT.txt","w") as fh:
    fh.write("\n".join(OUT)+"\n")
