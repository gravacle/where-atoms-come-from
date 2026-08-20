"""E-11  THE LANE'S FINAL CLASSIFICATION LEDGER.

Every scalar this lane could build, with its EXACT verdict, the control that would have caught a
mistake, and -- the column that decides whether a result survives the weakness objection -- whether
the verdict rests on an EXACT ARGUMENT valid at all N or only on exact values at the N tested.

Three status columns, and they are not the same thing:
  CLASS         (Z) exactly zero / (NZ) exactly non-zero
  GAUGE         does the quantity survive replacing each record representative a by a xor s?
                A quantity that does not is a property of the representative, not of the record.
  CARRIER-FREE  is the quantity computable from (H,{L_k}) alone, or does it need the locality
                structure (the tensor factorisation into sites)?  The model's own docstring records
                that clauses (i)-(iv) are carrier-free and clause (v) is not; the same split runs
                straight through this ledger and is its single most consequential feature.
"""
import sys

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()


ROWS = [
    # quantity, class, gauge-inv, carrier-free, basis, largest N, script
    ("sp_F2(R_a,R_a)                    [CTRL-Z]", "Z", "yes", "yes",
     "PROOF: sp(a,a) = 2(x.z) = 0 mod 2, all n", "n=64", "e1"),
    ("sp_F2 of a conjugate pair        [CTRL-NZ]", "NZ", "yes", "yes",
     "EXACT = 1 at every n tested (C-34)", "n=64", "e1"),
    ("sp_F2 pairing matrix, rank", "NZ", "yes", "yes",
     "EXACT full rank 2k at every n tested (C-34)", "n=64", "e1"),
    ("[R_a,R_b] exact matrix", "Z/NZ", "yes", "yes",
     "PROOF: 0 iff sp=0; else ||.||_F^2 = 4*2^n exactly", "n=16 matrix / all n symbolic", "e1"),
    ("I_unsigned / I_signed / I_overlap", "NZ", "NO", "no",
     "EXACT non-zero, but the value MOVES under a -> a xor s", "n=64", "e1,e5"),
    ("[[R_a,R_b],R_c] associator", "Z/NZ", "yes", "yes",
     "PROOF: !=0 iff sp(a,b)=1 AND sp(a^b,c)=1", "n=64", "e2"),
    ("Jacobi cyclic sum                 [CTRL-Z]", "Z", "yes", "yes",
     "PROOF: identity in any associative algebra", "n=64", "e2"),
    ("tau = Tr(Pi Ra Rb Rc)/Tr(Pi)", "Z/NZ", "NO (sign)", "yes",
     "EXACT: 0 unless closed; else a 4th root of unity", "n=10", "e2,e3"),
    ("Im tau (Bargmann invariant)", "Z/NZ", "NO (sign)", "yes",
     "PROOF: !=0 iff sp_ab+sp_bc+sp_ac is ODD", "n=10", "e2"),
    ("tau^2, K = (Ra Rb Rc)^2", "NZ", "yes", "yes",
     "EXACT = (-1)^(sp_ab+sp_bc+sp_ac): pairwise-determined", "n=12", "e2,e3"),
    ("EVERY gauge-inv word, length <= 8", "Z/NZ", "yes", "yes",
     "EXACT: 0 of 1848 words separate same-signature triples", "n=8", "e3"),
    ("any Z[i]-coefficient record polynomial", "quantised", "-", "yes",
     "PROOF: the Pi-trace ratio lies in Z[i]; no 0<|v|<1", "n=20", "e6"),
    ("[A_h,R] on D(G), G abelian        [CTRL-Z]", "Z", "yes", "yes",
     "PROOF: A_h = I exactly for every abelian G, all |G|", "|G|=8", "e4"),
    ("[A_h,R] on D(D_4) and D(Q_8)     [CTRL-NZ]", "NZ", "yes", "yes",
     "EXACT: ||[A_h,R]||_F^2 = 352/9 for the witness record", "|G|=8, dim 64", "e4"),
    ("commutant gap on D(G), non-abelian", "NZ", "yes", "yes",
     "EXACT: 1384 - 736 = 648 on both order-8 non-abelian groups", "|G|=8", "e4"),
    ("cross-region sp / comm / assoc / tau", "Z", "yes", "yes",
     "PROOF: disjoint support; 0 at EVERY separation d >= 1", "n=40, m=6 blocks", "e5"),
    ("J = minimal-crossing number of a pair", "NZ", "yes", "NO -- needs sites",
     "EXACT: sp=0 yet J=2 on 360 of 32385 pairs at n=6 (exhaustive)", "n=16", "e6,e7"),
    ("T = minimal three-way meeting number", "Z*", "yes", "NO -- needs sites",
     "UNINFORMATIVE: its own positive control also returned 0 (D-15)", "n=12", "e8"),
    ("GAP = Jsum - (J_ab+J_bc+J_ac)", "Z at n=4, NZ from n=6", "yes", "NO -- needs sites",
     "EXACT: 0 on all 3375 triples at n=4 (exhaustive); GAP=2 witness at n=6", "n=16", "e9,e10"),
    ("GAP, determinacy by pairwise data", "NZ", "yes", "NO -- needs sites",
     "EXACT: group (sp,J)=(0,0,0,0,0,0) holds GAP in {0,2,4} -- NOT pairwise-determined", "n=16", "e10"),
    ("GAP across disjoint regions       [CTRL-Z]", "Z", "yes", "NO -- needs sites",
     "EXACT: 0 across regions while the on-block control gives {0,2} in the same table", "n=12, m=2", "e10"),
    ("GAP instrument control (diagonal)  [CTRL-NZ]", "NZ", "yes", "NO -- needs sites",
     "EXACT: registers {0,2} already at n=4 where free GAP is 0 -- instrument is not blind", "n=14", "e9"),
]

say("=" * 132)
say("E-11  FINAL CLASSIFICATION LEDGER -- EVERY SCALAR THIS LANE BUILT, DECIDED EXACTLY")
say("=" * 132)
say("  %-44s %-22s %-12s %-18s %-14s"
    % ("quantity", "CLASS", "GAUGE-INV", "CARRIER-FREE", "largest N"))
say("  " + "-" * 128)
for q, c, g, cf, b, N, s in ROWS:
    say("  %-44s %-22s %-12s %-18s %-14s" % (q, c, g, cf, N))
say("")
say("  BASIS OF EACH VERDICT  ('PROOF' = an exact argument valid at ALL N and so immune to the")
say("  weakness objection; 'EXACT' = exact values at the N tested, and no further)")
say("  " + "-" * 128)
for q, c, g, cf, b, N, s in ROWS:
    say("  %-44s %-8s %s" % (q, "[" + s + "]", b))

say("")
say("=" * 132)
say("  THE SPLIT THIS LEDGER MAKES")
say("=" * 132)
say("  CARRIER-FREE COLUMN = yes.  Every scalar computable from (H,{L_k}) alone came out either")
say("    exactly zero, or representative-dependent (so not a record property at all), or an exact")
say("    function of the two-body F_2 pairing.  1848 gauge-invariant words in three records, to")
say("    length 8, separate NOTHING that the pairing does not already separate.  And every such")
say("    scalar is QUANTISED into Z[i]: it is exactly 0 or it is O(1).  There is no room in this")
say("    algebra for a small-but-non-zero residual, at any n.")
say("")
say("  CARRIER-FREE COLUMN = no.  The moment the site decomposition is admitted -- the same")
say("    locality data clause (v) needs and the four other clauses do not -- exactly-non-zero")
say("    quantities appear that the pairing cannot see: J at the two-body level, and GAP at the")
say("    three-body level.  GAP is EXACTLY ZERO on all 3375 triples at n=4 and EXACTLY NON-ZERO")
say("    from n=6, is NOT a function of the six pairwise invariants, and is EXACTLY ZERO across")
say("    disjoint regions at every separation.")
say("=" * 132)

with open(LANE + "/e11_ledger.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")
