# Independent post-freeze hostile audit — GL6AI F3 relational-influence envelope V001

**Target:** `LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001/`  
**Frozen theorem SHA-256:** `a51e802f6ba148e5f9848e95f41a80073795b24b7eaf87e36c0766b0856aa494`  
**Frozen MANIFEST-file SHA-256:** `fc50cad54dca00aab1c30d7c12ef07147df1242f94483f63955185695073f706`  
**Frozen SEAL-file SHA-256:** `12daa03d45cd653db24622ae8b3d8166015291534b3295e3c426eb37180fc918`  
**Disposition:** `PASS__EXACT_DEGREE_SQUARE_SPLIT__LINK_DEGREE_Q_PLUS2_LE6__DRESSING_COMPLETE_LAMBDA48_EXACT__RETAINED_SOURCE_K_DUHAMEL_TYPED__LINK_TO_CELL_DISTANCE_DESCENT__UNIFORM_ANALYTIC_EXPONENTIAL_TAIL__NO_EXACT_SPEED_LORENTZ_RICCI_GRAVITY_OR_G_PROMOTION`

## Frozen custody and replay

The requested theorem, manifest-file, and seal-file pins match exactly.  This
audit pins every frozen author file and copies the frozen dependency ledger
byte for byte.  In particular, that ledger pins the frozen `GL6AH` theorem and
seal and its independent audit and audit seal at the requested values.  No
author byte was edited.

The frozen packet replay returns `47/47`, and the underlying author exact
census returns `123777/123777`.  The separate audit reconstruction imports no
author program or ledger.  It rebuilds finite FPSS slabs, exhausts every link
occupation on `N=0,1`, samples larger slabs with exact rational parameters,
exhausts link/cell distances through `N=4`, and checks influence powers and
tail coefficients directly.

## Exact split and topology

For one endpoint `v`, idempotence gives

\[
 (d_v-d_\star)^2=d_\star^2+(1-2d_\star)
 \sum_{e\ni v}n_e+2\sum_{e<f\ni v}n_en_f.
\]

Every physical link has two endpoints, yielding the onsite coefficient
`Delta+2 Ud(1-2 d_star)`.  Two distinct links in the simple incidence graph
share at most one endpoint, so each line-graph pair receives exactly
`2 Ud n_e n_f` once.  The onsite evolution is a product of one-link unitaries;
conjugation therefore preserves each pair's two-link support and norm
`J=2|Ud|`.

For `e=(m,c)`, three other links share `m`; exactly `q(c)-1` other links share
`c`.  These sets are disjoint, so `deg(e)=q(c)+2<=6`, with equality on an
all-positive child.  The bound is exact and independent of slab size.

## Dressing-complete influence constant

For every pair incident to an already reached link, the Duhamel-Jacobi step
has one advancing and one non-advancing channel.  Encoding the former by an
off-diagonal `J` and the latter by diagonal mass `J deg(e)` gives row mass

\[
 J\deg(e)+\sum_{f\sim e}J=2J\deg(e)\le2J\Delta_L.
\]

The commutator inequality supplies the outer factor `2/hbar`.  Therefore

\[
 {2\over\hbar}(2J\Delta_L)
 ={4J\Delta_L\over\hbar}
 ={48|U_d|\over\hbar}.
\]

Diagonal matrix steps do not move an endpoint; any entry connecting links at
distance `d` still needs at least `d` off-diagonal steps.  The audit's exact
sparse powers verify this distance filtration and the row-power ceiling.  No
support-cardinality factor is hidden in this rate: a finite output `Y` enters
only through the initial indicator vector, hence the displayed sum over
`f in Y`, or the later explicit factor `|Y|` after using the minimum distance.

## Retained source comparator and distance descent

The frozen comparator is normalized `beta_s=1` minus `beta_s=0` on the same
active-link Hilbert factor, with every `beta_-s` fixed.  Its exact Hamiltonian
difference is `V_s=-hX_s`, of norm `h`; the Duhamel identity contributes
`1/hbar`, and the commutator bound contributes `2`, giving the certified
`2h/hbar` source prefactor.  This is neither success filtering nor a deletion
switch, and it does not identify sham/KEEP with formed/BREAK full-instrument
ancestry.  Observables that directly read the orthogonal route record are
outside this fixed-sector active-output statement.

Display (AI16) has a missing backslash before `qquad`.  This is a non-material
frozen typesetting typo: the same display still states `V_s=-hX_s` and
`||V_s||=h`, the surrounding prose states both facts, and all downstream
formulas and exact checks use them correctly.  Author bytes remain untouched.

A same-parent line-graph step leaves the parent-cell label fixed.  A
same-child step joins two authenticated adjacent parents.  Consequently any
link path, after repeated labels are deleted, projects to a genuine cell walk,
so `d_link>=d_cell=||m-n||_1/2`.  Taking the minimum over a finite multi-cell
output gives exactly the distance used in the frozen source bound.  Since the
degree and matrix norm ceilings do not depend on `N`, the result is uniform for
each fixed finite source/output block.

## Tail and strict ceiling

Termwise integration proves
`integral_0^t T_d(lambda u)du=lambda^-1 T_(d+1)(lambda t)` for
`lambda>0`.  For `z=e^mu>1`, the coefficientwise inequality
`1<=z^(r-d)` for every `r>=d` proves
`T_d(x)<=exp(zx-mu d)`.  At `Ud=0` the pair generator vanishes and distinct
links factorize, as stated.

The analytic tail is generally nonzero for every positive time.  The theorem
therefore establishes only an exponentially quasi-local relational-influence
upper envelope.  No Lorentz or common physical cone, exact finite-speed
support, stationary bulk mode, continuum or infrared propagator, Ward/Bianchi
closure, Ricci/Einstein response, gravity, or Newton's `G` is inferred.

**Audit verdict: PASS.**
