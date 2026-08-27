# Same-parent multiscale causal-measure theorem

**Claim class:** exact arbitrary-scale reduction-error theorem, exact finite
support/order-coherence theorem, exact compatible-measure projective-limit
theorem, and sharply conditional Lorentzian reconstruction corollary

**Not claimed:** that a finite graph sequence produces a manifold, that the
limiting cone is Lorentzian without a realizability test, that a calibrated
spatial volume is spacetime four-volume, that records cause the limiting
geometry, that gravitationally dressed carriers are compact, or that Einstein
dynamics, universal coupling, \(G\), or \(\Lambda\) have been derived

## 1. Same-parent multiresolution object

Let \(n=0,1,2,\ldots\) index increasing physical resolution.  Thus \(n>m\)
means that level \(n\) is finer than level \(m\).  Let \(E_n\) be a finite set
of event cells and let

\[
 q_{mn}:E_n\twoheadrightarrow E_m,
 \qquad q_{\ell n}=q_{\ell m}\circ q_{mn}
 \quad(n>m>\ell)                                      \tag{1}
\]

be frozen surjective cell maps.  They are physical resolution maps, not maps
chosen after comparing response outcomes.

There is one complete parent process \({\cal K}_\star\) and one frozen
resolution rule \(\Pi_{p,n}\) for each admitted probe family \(p\),

\[
 {\cal K}_{p,n}=\Pi_{p,n}{\cal K}_\star .             \tag{2}
\]

Equation (2) is the **same-parent** premise.  Separately fitted processes that
happen to have similar graphs do not satisfy it.  Every \({\cal K}_{p,n}\)
contains the complete registered intervention law, noise, memory, higher
response or a certified tail, declared external ports, and source/boundary
ledger required by its mission.  A first derivative alone is not a complete
process unless a physical completeness theorem says that it is.

Let \({\cal B}_{m,m+1}\) be the adjacent physical reduction from level
\(m+1\) to level \(m\), and write

\[
 {\cal B}_{m\leftarrow n}
 ={\cal B}_{m,m+1}\circ\cdots\circ{\cal B}_{n-1,n}.
                                                               \tag{3}
\]

A separately implemented direct reduction is denoted
\({\cal B}^{\rm dir}_{mn}\).  The process metrics \(d_n\) are frozen,
dimensionless or physically calibrated, and operationally dominate the complete
registered laws.  Memory poles, state dimension, tail freedom, and scorer
complexity may not be enlarged after a failed scale.

## 2. Theorem M1 -- arbitrary-scale predictive-error composition

For one fixed probe family or the joint complete process, define adjacent
population residuals

\[
 \delta_k=d_k\!\left({\cal K}_k,
       {\cal B}_{k,k+1}{\cal K}_{k+1}\right).          \tag{4}
\]

Assume the adjacent reductions are Lipschitz in the same metrics,

\[
 d_k({\cal B}_{k,k+1}A,{\cal B}_{k,k+1}B)
 \le L_k d_{k+1}(A,B),\qquad L_k<\infty .              \tag{5}
\]

Then for every \(n>m\),

\[
\boxed{
 d_m\!\left({\cal K}_m,
        {\cal B}_{m\leftarrow n}{\cal K}_n\right)
 \le
 \sum_{r=m}^{n-1}
 \left(\prod_{s=m}^{r-1}L_s\right)\delta_r .}
                                                               \tag{6}
\]

The empty product in (6) is one.  If

\[
 \zeta_{mn}=d_m\!\left(
 {\cal B}^{\rm dir}_{mn}{\cal K}_n,
 {\cal B}_{m\leftarrow n}{\cal K}_n\right),          \tag{7}
\]

then

\[
 d_m\!\left({\cal K}_m,
        {\cal B}^{\rm dir}_{mn}{\cal K}_n\right)
 \le
 \sum_{r=m}^{n-1}
 \left(\prod_{s=m}^{r-1}L_s\right)\delta_r+\zeta_{mn}.
                                                               \tag{8}
\]

### Proof

For one adjacent step, (6) is (4).  Insert
\({\cal B}_{m,m+1}{\cal K}_{m+1}\), apply the triangle inequality,
then (5), and repeat at the next scale.  This gives the weighted telescoping
sum (6).  Inserting \({\cal B}_{m\leftarrow n}{\cal K}_n\) and applying the
triangle inequality gives (8).  QED.

M1 generalizes the sealed two-block G5 bound.  It does not infer a continuum
from a small endpoint residual: every adjacent residual and the direct-versus-
sequential residual remain visible.

## 3. Common probe fronts

For each complete probe process, let \(S_{p,n}\subseteq E_n\times E_n\) be its
exact nonzero interventional influence relation, including the reflexive pairs.
At every scale assume:

1. **closed probe domain:** the declared family contains every sector to which
   the words *maximal* and *universal front* refer;
2. **complete ports:** the intervention/read family separates every admitted
   process difference and includes every external bypass;
3. **common maximal envelope:** there is a relation \(C_n\) such that
   \(S_{p,n}\subseteq C_n\) for every admitted \(p\), while at least two
   physically independent, separately actuated and separately read families
   satisfy \(S_{p_1,n}=S_{p_2,n}=C_n\); and
4. **independence through scale:** common support is not produced by one shared
   coarse actuator, preprocessing map, detector, or selection rule.

Two saturating probes establish replication only inside the closed declared
domain.  They do not prove that an omitted faster sector does not exist.

## 4. Theorem M2 -- support-faithful direct/sequential order coherence

For a relation \(R\subseteq E_n\times E_n\), define its existential
pushforward

\[
 (q_{mn})_*R
 =\{(a,b):\exists x,y, q_{mn}(x)=a, q_{mn}(y)=b, (x,y)\in R\}.
                                                               \tag{9}
\]

For each probe family, let \({\mathfrak D}_{p,n}\) be an admitted domain of
complete level-\(n\) process objects.  Require it to contain the actual process
and to be **path closed** under every registered adjacent and direct block:

\[
 {\cal K}_{p,n}\in{\mathfrak D}_{p,n},\qquad
 {\cal B}:{\mathfrak D}_{p,n}\longrightarrow
              {\mathfrak D}_{p,m}.                     \tag{10a}
\]

Thus \({\mathfrak D}_{p,m}\) contains not only the actual
\({\cal K}_{p,m}\), but every predicted intermediate such as
\({\cal B}_{m\leftarrow n}{\cal K}_{p,n}\) that can be fed to a later
registered block.  Assume every such block is a map-level **support
homomorphism** on this whole domain:

\[
 \operatorname{supp}({\cal B}X)
       =(q_{mn})_*\operatorname{supp}(X)
 \quad\text{for every }X\in{\mathfrak D}_{p,n}.         \tag{10b}
\]

Also require the actual processes to obey
\(S_{p,m}=(q_{mn})_*S_{p,n}\).  Equations (10a)--(10b), rather than a
pointwise check only on each actual \({\cal K}_{p,n}\), are the
**support-faithfulness** premise.  Equivalently, one may verify every registered
composite sequential route itself as support-faithful.  The gate must be tested
for the complete interventional family.  Positive injections and reads on
entrywise nonnegative complete responses are one sufficient realization; signed
averaging is not sufficient in general.

Assume the fine envelope is a partial order (or that a declared strongly
connected causal block quotient has first been taken).  Also assume **causal
lumpability**:

\[
 C_m=(q_{mn})_*C_n                                      \tag{11}
\]

is reflexive, transitive, and antisymmetric at every registered scale.  In
elementary form, transitivity requires that whenever some fine representative
of block \(A\) influences a representative of \(B\), and some possibly
different representative of \(B\) influences a representative of \(C\), there
is a fine influence from a representative of \(A\) to one of \(C\).  The
analogous two-block-cycle condition supplies antisymmetry.

Then:

1. every scale has the same declared maximal front for the two independent
   saturating probes, and all other admitted probes remain inside it;
2. every \(C_n\) is a partial order;
3. support from a direct fine-to-coarse block equals support from every
   sequential block path; and
4. for \(n>m>\ell\),
   \[
    (q_{\ell n})_*C_n
       =(q_{\ell m})_*\bigl((q_{mn})_*C_n\bigr)=C_\ell .
                                                               \tag{12}
   \]

### Proof

Existential relation pushforward is functorial:
\((r\circ q)_*R=r_*(q_*R)\).  Path closure in (10a) permits (10b) to be
applied inductively to every predicted intermediate, giving

\[
 \operatorname{supp}({\cal B}_{m\leftarrow n}X)
 =(q_{mn})_*\operatorname{supp}(X).
\]

A registered direct block has the same support by (10b), and actual-process
support compatibility gives the common fronts at every scale.  Applying
functoriality to each saturating support proves (12).  Equation (11) supplies
the partial-order properties that existential pushforward does not supply by
itself.  QED.

M2 proves support and order coherence, not equality of the complete direct and
sequential kernels.  Complete-law agreement is controlled separately by M1 and
the frozen operational metric.

## 5. Theorem M3 -- calibrated volume compatibility and projective limit

For the joint causal-measure conclusion, let the event cells \(E_n\) also be the
atoms of a finite event-region algebra.  The cell map induces the nonnegative
aggregation matrix \(A_{mn}\), with exactly one unit entry in each fine-cell
column.  Hence

\[
 A_{\ell n}=A_{\ell m}A_{mn}.                          \tag{13}
\]

At every scale, an independently calibrated extensive observable and a
separately established constitutive identification give

\[
 \mathsf{M}_n v_n=b_n,\qquad \operatorname{rank}(\mathsf{M}_n)=|E_n|,
 \qquad v_n>0.                                         \tag{14}
\]

Thus \(v_n\) is the unique positive finite-volume assignment at that scale.
A packet that calibrates only spatial atoms may use a separate spatial region
algebra, but then M3 supplies a separate spatial measure rather than a measure on
the causal event-cell limit.  Joining the two requires an explicit compatible
event-cell/product map and, for M4, the independent four-volume identity.
Add the cross-scale physical compatibility test

\[
 \boxed{v_m=A_{mn}v_n.}                                \tag{15}
\]

Then all direct and sequential volume assignments agree exactly.  If adjacent
compatibility has signed \(\ell^1\) residual

\[
 \eta_k=\|v_k-A_{k,k+1}v_{k+1}\|_1,                  \tag{16}
\]

then

\[
 \|v_m-A_{mn}v_n\|_1\le\sum_{k=m}^{n-1}\eta_k.       \tag{17}
\]

### Proof

Exact direct/sequential agreement follows from (13) and repeated use of (15).
For (17), insert \(A_{m,m+1}v_{m+1}\), apply the triangle inequality, and
iterate.  Pushforward by a partition aggregation is \(\ell^1\)-contractive on
signed measures, so no Lipschitz factor greater than one appears.  QED.

Suppose now that the scale family is countable, every cell map is surjective,
and the finite measures \(\mu_n\) determined by \(v_n\) obey

\[
 (q_{mn})_*\mu_n=\mu_m,\qquad
 0<\mu_n(E_n)=M_{\rm tot}<\infty .                    \tag{17a}
\]

Thus (15) holds exactly with one finite common normalization.  The inverse limit

\[
 X=\varprojlim(E_n,q_{mn})                              \tag{18}
\]

is nonempty and compact in the cylinder topology.  The compatible finite
measures define a unique regular Borel measure \(\mu_\infty\) on \(X\).  If
every \(C_n\) satisfies M2, then

\[
 x\preceq_\infty y
 \quad\Longleftrightarrow\quad
 x_n\,C_n\,y_n\ \text{for every }n                    \tag{19}
\]

is a closed partial order.  Reflexivity and transitivity hold coordinatewise;
antisymmetry at every finite scale makes equal all coordinates and hence the
inverse-limit points.  Each \(E_n\) is finite discrete, so each \(C_n\) is
closed; (19) is an intersection of inverse images of these closed relations and
is therefore closed.  The normalized measures \(\mu_n/M_{\rm tot}\) form a
consistent finite-cylinder family.  The standard extension theorem gives a
unique probability measure on the cylinder sigma-algebra, hence a unique regular
Borel measure after restoring \(M_{\rm tot}\).

Therefore M1--M3 can earn an exact **same-parent projective causal-measure
space**.  They do not imply that \(X\) is connected, locally Euclidean,
four-dimensional, differentiable, or Lorentzian.  A binary refinement system,
for example, can have perfect order/measure compatibility and a Cantor-like or
branching limit.

## 6. Corollary M4 -- sharply conditional Lorentzian causal metric-measure limit

The projective result reaches the existing G3 bridge only if all of the following
are independently established for one fixed source/background sector:

1. **manifold realization:** the same-parent cells have a nested realization
   generating the Borel sets of one connected \(C^2\) four-manifold \(M\), with
   mesh tending to zero in an independently constructed operational atlas;
2. **complete order convergence:** the lifted common-front relations converge,
   with inner and outer control on compact subsets away from the null boundary,
   to the complete chronological relation \(I^+\), not sparse adjacency or a
   selection-rule support; the limit is future- and past-distinguishing or
   strongly causal;
3. **Lorentzian tangent-cone convergence:** local blow-ups of the outer fronts
   converge to one smooth, time-oriented, nondegenerate quadratic cone of
   signature \((-+++)\), common to a physically complete domain of clocks,
   matter, light, and the candidate tensor response;
4. **process autonomy:** the M1 adjacent and direct/sequential complete-law
   errors, including noise, memory, higher kernels, and certified tails, meet a
   frozen convergence bound with operational domination;
5. **measure convergence:** the compatible cell measures converge on continuum
   Borel/Alexandrov regions to a smooth strictly positive locally finite density
   \(\mu\); and
6. **four-volume identity:** the extensive calibration is independently shown
   to be the absolutely normalized spacetime metric four-volume.  A spatial
   layer volume alone is insufficient unless a separately calibrated proper-time
   or lapse law proves the conversion to four-volume.

Under these premises, the common chronological relation fixes one conformal
Lorentzian class \([g]\).  For any representative \(g_0\), write

\[
 d\mu=f\,dV_{g_0},\qquad f>0.
\]

In four dimensions,

\[
 \boxed{g=f^{1/2}g_0,\qquad dV_g=f\,dV_{g_0}=d\mu.}    \tag{20}
\]

Thus \((I^+,\mu)\) reconstructs one Lorentzian metric-measure structure up to
diffeomorphism.  If the measure normalization is known only up to one constant,
one global metric scale remains.

M4 is a conditional reconstruction theorem.  Its manifold, complete chronology,
quadratic-cone, and metric-volume premises are not consequences of finite front
and volume compatibility.  If the reconstructed metric is \(C^2\), its
Levi-Civita connection and Riemann tensor exist mathematically.  A physical
connection still requires probe-transport agreement and bounds on torsion,
nonmetricity, dispersion, and anisotropy.  No curvature-error estimate follows
without at least \(C^2\) metric convergence and a stability theorem.

## 7. Dynamical-gravity port boundary

Gravitational dressing forbids treating every physical carrier as a compact
local subsystem.  At every scale, charge/flux, edge/dressing, homogeneous or
radiative modes, initial-cap data, and lateral worldtube ports must be retained
when they carry event contrast.  Discarding them invalidates same-parent and
support-completeness premises.

The sealed gravitational-separator lane supplies only a restricted approximate
alternative: in the D=4 charge-matched localized-matter/common-dressing source
class, declared bounded linear one-field reads have exterior event contrast
\(O(\kappa^2)\), with a finite-read tail bound only after its field, detector,
and remainder constants are supplied.  A fixed residual
\(C_j\kappa^2\) does not vanish merely by taking more resolution levels.
Therefore it cannot establish the exact projective or continuum theorem unless
all distributed ports are retained exactly or an independent higher-order
argument proves a scale-vanishing gravitational error.

## 8. Exact dynamics and ancestry ceiling

M1--M4 establish, at most, autonomous Lorentzian metric **kinematics** in a
fixed source/background sector.  They do not imply:

- a positive protected massless spin-2 pole or two radiative helicities;
- universal coupling to complete stress-energy;
- diffeomorphism Ward identities or nonlinear self-coupling;
- the Einstein--Hilbert action, horizon equilibrium, \(G\), or \(\Lambda\); or
- that record formation caused the causal-measure limit.

Common propagation cones are not universal coupling strengths.  Einstein
dynamics still requires all premises of the existing G4A spin-2/EFT route or
the alternative G4B horizon-equilibrium route.  Same-parent provenance is also
weaker than record ancestry: a common cause can generate both the record and the
scale family.  Lifecycle-qualified KEEP/BREAK/reprepare mediation and a complete
source/work/heat/boundary ledger remain required for the record-origin claim.
