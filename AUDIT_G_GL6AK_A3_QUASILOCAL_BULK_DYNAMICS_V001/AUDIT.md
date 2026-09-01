# Independent hostile pre-freeze audit — GL6AK A3 quasi-local bulk dynamics V001

**Target:** `LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001/`  
**Reviewed mutable theorem SHA-256:** `38cb58ef9fc52e1252e0b0d3415c54488c0471c0ddf25a35b7adf5aba41bccc9`  
**Disposition:** `CLEAN_PREFREEZE__AUTHOR_BYTES_NOT_FROZEN_OR_EDITED`

## Independent reconstruction

The author replay returns `6304/6304` and the mutable structure check returns
`26/26`.  The separate hostile replay imports neither the author program nor
its ledger.  It reconstructs the infinite edge relation, finite translated
collars, shell ownership, finite influence matrices, distance filtration,
factorial tails, Følner averages, and the six-pair irreducible projectors.

The oriented shared-edge convention is complete.  When `a>b`, the edge from
`(x,a)` to `(x+e_a-e_b,b)` is the same unordered edge as the canonical
`b<a` representative anchored at the other parent.  Every active-link site
therefore has exactly three same-parent and three same-child partners.  A
cell touches six internal and twelve shared-child terms.  Translation and
`S4` preserve the underlying two edge relations; an `S4` permutation may
merely re-anchor the canonical shared-edge representative.

For any finite collared patch, the displayed translation sends every parent
to one common `S_N` strict interior and preserves literal child equality.
This establishes finite local ancestry only.  The theorem expressly rejects
one infinite terminal record or simultaneous global authentication.

## Boundary comparison and bulk dynamics

Assigning every pair term to a minimum-radius endpoint cell gives at most
eighteen terms per owner and hence the conservative shell bound

\[
 N_r\le18|B_r\setminus B_{r-1}|\le18(2r+1)^3.
\]

Rerunning the positive Duhamel--Jacobi recursion from a finite support
indicator gives the stated double support union bound without decomposing an
arbitrary operator.  For a two-site pair term, the endpoint sum supplies
`4J ||A|| |X|`.  Shell counting and Duhamel then give
`72J ||A|| |X|/hbar`.  Since the hostile-audited dressing-complete constant is

\[
 \lambda_{\rm F3}=4J\Delta_L/\hbar=24J/\hbar,
\]

termwise integration gives exactly the factor `3` and tail index `+1` in

\[
 3\|A\||X|\sum_{r\ge R}(2r+1)^3
 T_{r-r_X+1}(\lambda_{\rm F3}|t|).
\]

The factorial tail beats the cubic shell uniformly on compact time
intervals.  Comparison with an included complete ball proves the claimed
independence for locally complete exhaustions formed only by restricting the
inherited interaction.  New or arbitrarily strong boundary laws are correctly
excluded.  Treating the onsite product dynamics exactly in the interaction
picture preserves both support and pair-term norm.

## Invariance and stationary response

The `A3` translations and `S4` relabelings preserve both the graph and every
coefficient.  Time averaging requires weak-* cluster nets, not sequential
compactness; the theorem says only “cluster point” and is sound.  The explicit
`Z3` cube is Følner.  Translation averaging preserves stationarity because
the actions commute, and finite `S4` averaging preserves both stationarity
and full translation invariance.

In the invariant state's GNS representation, the Liouvillian spectral theorem
gives a finite positive matrix-valued measure.  With the theorem's convention
`U(t)=exp(itL)`, its coordinate is angular frequency and

\[
 F_{AB}(t)=\int e^{it\nu}\,\mu_{AB}(d\nu).
\]

Stationarity gives
`F_AB(-t)-F_BA(t)=omega([tau_t(M_A),M_B])`, so the retarded sign and index
ordering are correct.  The six-edge permutation representation of `S4` is
multiplicity-free `A1+E+T2`; Schur's lemma therefore yields three positive
scalar measures and no off-block response.

## Strict ceiling

The theorem derives a mathematical homogeneous thermodynamic completion and
stationary response object from finitely authenticated local ancestors.  It
does not turn the completion into a globally authenticated record, select a
state, establish finite-volume state convergence, call a character physical
momentum, or infer a pole, gapless mode, hydrodynamics, physical cone, metric,
Ricci response, gravity, or `G`.

**Hostile pre-freeze verdict: CLEAN.**  A distinct post-freeze custody audit
must pin and replay the final author bytes before canonical promotion.
