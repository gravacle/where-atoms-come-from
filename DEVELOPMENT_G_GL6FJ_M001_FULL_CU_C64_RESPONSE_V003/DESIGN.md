# Design and exact scope

For 24 real P/C quadrature channels, a real symmetric response has

\[
\dim \operatorname{Sym}(24)=24\cdot25/2=300.
\]

The catalog measures the 24 unit directions and the 276 normalized pair
sums.  Their quadratic projections reconstruct every response entry with no
unmeasured row and no redundant production probe.

At `L=4`, `m001=(0,0,1)` is distinct from its signed lift `-m001=(0,0,-1)`.
The engine verifies signed conjugate transport on every node/pair, retains
the parent/child centering phase, and verifies a 24-real unaliased chart.
This is precisely why `m001` is useful: together with the already measured
`m123` response it is the least expensive exact orbit pair known to retain
the GL6CR rank-eight algebraic screen.  That algebraic fact does not itself
assert that the measured physical response satisfies the Ward contraction.

The source-zero schedule is fixed across all 300 rays.  Each raw checkpoint
contains `stay_contact`, `writer_contact`, `stay_spectral`,
`writer_spectral`, `stay_writer_interference`, and `disconnected` exactly
once.  The postprocessor reconstructs every owner before any aggregation and
checks their matrix recombination.
