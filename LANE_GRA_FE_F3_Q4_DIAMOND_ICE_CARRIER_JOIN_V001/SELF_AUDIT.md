# Q4DICJ builder self-audit

## Result

The new theorem is a graph/support join.  It proves that the translation-
invariant completion of the incidence graph between two consecutive q4 count
fronts is exactly the standard diamond net: an `A3`/FCC Bravais lattice, two
cosets, degree four, tetrahedral bonds, and girth six.  It also proves that the
finite nonnegative q4 slabs locally exhaust that net.

Under the explicitly supplied `Q4-CARRIER/EDGE-LIFT`, the existing `d_*=2`
F3 diamond-ice theorem can therefore use q4 append incidence as its support
shape.  The existing ice bijection, sixth-order coefficient `63/8`, scalar
diagonal (`V_6=0`), and pure-kinetic `mu=0` convention match are inherited;
they are not recalculated by adding a new physical interaction.

## Exact graph checks

1. The two affine `Z^4` hyperplanes differ by one in total count and are
   connected only by the four append vectors.
2. Projection onto the regular tetrahedral frame maps a fixed-total
   hyperplane to the root lattice `A3`, explicitly the FCC parity lattice.
3. The next hyperplane maps to the standard `(1/4,1/4,1/4)` diamond basis
   coset in conventional cubic coordinates.
4. Distinct same-part vertices share at most one opposite-part neighbor, so
   four-cycles are impossible.
5. A six-cycle is exhibited.  Closure and nonbacktracking prove that every
   simple six-cycle uses three distinct append labels and is a standard
   diamond hexagon.
6. A radius-`r` ball based at a count vector with every coordinate at least
   `r` never touches the nonnegative boundary, proving exact local
   exhaustion.

The verifier independently checks all of these statements, finds the expected
twelve undirected diamond hexagons through one vertex, constructs a safe
`L=4` periodic quotient, and re-sums all 720 alternating-hexagon orders to
obtain `63/8`.

## Principal overclaim guards

1. **Negative integer counts are not physical.**  They occur only in the
   mathematical translation completion.  Every finite-radius neighborhood is
   recovered inside a sufficiently deep nonnegative slab.
2. **Front alternatives are not sites.**  The q4 theorem supplies Hilbert
   basis labels and histories.  `Q4-CARRIER` is the unproved physical step
   which realizes those labels as coexisting carrier modes.
3. **Append edges are not automatically matter links.**  `EDGE-LIFT` is the
   unproved step which binds each incidence key to one binary F3 link while
   retaining all lineage ports.
4. **Carrier flips do not reverse record formation.**  `X_e` changes the
   lifted link occupation and acts as identity on the authenticated history
   factor.
5. **Finite slabs have a boundary.**  Their child degree is the number of
   positive coordinates, not always four.  The global regular-support F3
   theorem therefore needs an infinite controlled definition or a separately
   supplied periodic completion.
6. **A mathematical quotient is not a physical derivation.**  Safe periodic
   diamond quotients exist, but BQ4RSW does not identify distinct count/front
   lineages periodically.
7. **The phase result is leading and conditional.**  The Hamiltonian match is
   exact through sixth order on a fixed plaquette-complete support.  The U(1)
   assignment is inherited public numerical evidence for that effective
   model, not a volume-uniform all-orders theorem for the F3 parent.
8. **The degree sector matters.**  The same support at `d_*=1` has the
   negative ordered comparator.  Diamond incidence alone does not force the
   U(1) phase.
9. **No support force is proved.**  Physical binding and a restoring support
   basin remain open under the existing support-force theorem.
10. **No visible or tensor promotion occurs.**  The phase is spin one.  No
    visible-current map, alpha value, rank-two constraint, helicity-two pole,
    universal stress coupling, gravity, or `G` follows.

## Reproduction

Run:

```text
python3 LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/verify_q4_diamond_ice_carrier_join.py
```

Expected result: `Q4DICJ verification: PASS (21/21)`.

## Builder disposition

`READY_FOR_INDEPENDENT_HOSTILE_AUDIT`
