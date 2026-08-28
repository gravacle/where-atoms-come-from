# K5 processed-clock common-node-potential cycle theorem

**Theorem ID:** GRA-CLOCK-K5-CPC-V001  
**Date:** 2026-08-27  
**Claim class:** exact graph-algebra theorem plus covariance-limited descriptive
compatibility diagnostic on official processed pairwise clock outputs

**Status:**
PAIRWISE_ESTIMATED_NOT_COMMON_FIT_DERIVED__STATISTICAL_INDEPENDENCE_NOT_ESTABLISHED__K5_CUT4_CYCLE6_EXACT__MARGINAL_INTERVAL_COMPATIBLE__NO_GRAVITY_OR_LINEAGE_PROMOTION

## 1. Frozen source and exact semantics

The source is the CC BY 4.0 Zenodo process-data deposit for Zheng et al.'s
five-ensemble strontium clock experiment, DOI 10.5281/zenodo.8184043. All four
official CSVs were reacquired from the Zenodo API and reproduce the MD5,
SHA-256, byte-size, and shape ledger frozen in the SPAG second-pass lane.
Fig4.csv is the sole scored table: it contains ten pairwise inferred relative
clock-height differences and ten reported marginal $1\sigma$ uncertainties in
centimetres.

The article gives the critical provenance. The five clocks produced ten
simultaneous pairwise comparisons. In the analysis underlying these outputs,
each pair was re-analyzed individually over the same 14-run raw data, received
pair-specific systematic corrections, and was then weighted-averaged. The
reported inferred height is

\[
 h_{ij}=\frac{\delta f_{ij}}{f}\frac{c^2}{g},            \tag{K5-01}
\]

where the authors explicitly assume general relativity and use one independently
measured local $g$. Thus the ten values are separately estimated pairwise
outputs, not ten algebraically generated differences from one fitted five-node
vector. They are **not statistically independent**: pairs share clocks, raw
runs, and systematic channels, and the deposited files contain no joint
covariance matrix.

The non-derivation preflight is also visible in the deposited decimals. For the
triangle $2\to3\to4\to2$,

\[
 h_{23}+h_{34}-h_{24}=0.43+0.19-0.53=0.09\ \mathrm{cm}. \tag{K5-02}
\]

Under conventional nearest rounding, three values rounded independently to
$0.01\ \mathrm{cm}$ can accumulate at most $0.015\ \mathrm{cm}$ of
half-last-digit closure error. Therefore these displayed values cannot be an
exactly closed common-node vector subjected only to independent nearest
rounding at the displayed precision. This is a statement about that explicit
rounding model, not about unpublished digits or measurement uncertainty. The
article-level pairwise-analysis provenance independently authorizes the
non-tautological algebraic closure diagnostic; neither fact establishes
statistical independence.

## 2. Exact K5 cut/cycle theorem

Orient every edge $i\to j$ for $i<j$, and let
$B\in\mathbb Z^{10\times5}$ be the complete-graph incidence matrix with a
+1 at node $i$ and -1 at node $j$. A pair vector
$y\in\mathbb R^{10}$ comes from one five-node scalar $z$ exactly when

\[
 y=Bz.                                                    \tag{K5-03}
\]

Exact rational elimination gives

\[
 \operatorname{rank}B=4,\qquad
 \dim\operatorname{im}B=4,\qquad
 \dim\ker B^{\mathsf T}=10-4=6.                         \tag{K5-04}
\]

The cut projector and cycle projector are

\[
 P_{\rm cut}=B(B^{\mathsf T}B)^+B^{\mathsf T}
              ={1\over5}BB^{\mathsf T},\qquad
 P_{\rm cyc}=I-P_{\rm cut}.                             \tag{K5-05}
\]

They have exact ranks 4 and 6. Six star-tree triangles
$(1,i,j,1)$, $2\le i<j\le5$, form a full-rank cycle basis $C$, and

\[
 CB=0,\qquad y\in\operatorname{im}B\iff Cy=0
             \iff P_{\rm cyc}y=0.                       \tag{K5-06}
\]

This is a common-**node-scalar** theorem. A gravitational scalar potential must
pass it downstream, but so must every ordinary conservative node-difference
field. Passing it is necessary, never sufficient, for a common gravitational
potential or metric.

## 3. Complete descriptive cycle ledger

The verifier enumerates every simple K5 cycle uniquely up to rotation and
reversal: 10 triangles, 15 quadrilaterals, and 12 pentagons, for 37 total.
For a directed cycle $c$, it computes

\[
 r_c=c^{\mathsf T}y,
 \qquad s_c=\sum_e |c_e|\sigma_e.                        \tag{K5-07}
\]

The maximum absolute observed cycle residual is $0.29\ \mathrm{cm}$. The
largest ratio $|r_c|/s_c$ is $0.32926829268292684$. Every one of the 37
cycle residuals can therefore be canceled inside its componentwise reported
marginal $1\sigma$ box. The full values, orientations, exact integer weights,
exact rational residuals, and marginal envelopes are frozen in RESULT.json.

For orientation only, the unweighted Euclidean projection gives the
sum-zero-gauge node vector

\[
 \hat z=(0.548,\ 0.306,\ -0.066,\ -0.244,\ -0.544)\ \mathrm{cm}. \tag{K5-08}
\]

The projected cycle component has Euclidean norm
$0.1557562197795003\ \mathrm{cm}$ and maximum edge component
$0.092\ \mathrm{cm}$. This is not covariance-weighted and is not used as a
chi-square statistic.

## 4. Covariance-honest global compatibility bound

Because no ten-edge covariance was deposited, neither independent Gaussian
errors nor a conventional chi-square is admissible. The exact replacement is
the deterministic marginal-box problem

\[
 \rho_*=\min_z\max_e\frac{|y_e-(Bz)_e|}{\sigma_e}.        \tag{K5-09}
\]

Exact rational vertex enumeration, with node 5 fixed to zero as a gauge, gives

\[
\rho_*={27\over82}=0.32926829268292684<1.              \tag{K5-10}
\]

The exact lower certificate is the cycle $1\to3\to5\to4\to1$:

\[
 {|0.66+0.57-0.23-0.73|\over0.20+0.20+0.19+0.23}
 ={0.27\over0.82}={27\over82}.                         \tag{K5-10a}
\]

Every candidate node vector has zero sum around this cycle, so the triangle
inequality forces $\rho\ge27/82$. Conversely, the exact node-5-zero vector

\[
 z=\left({4503\over4100},{3587\over4100},{2067\over4100},
          {2399\over8200},0\right)\ \mathrm{cm}         \tag{K5-10b}
\]

has all ten standardized absolute residuals at most $27/82$. Thus the stated
optimum follows from matching rational lower and upper certificates and does
not depend solely on trusting the vertex-enumeration implementation.

Therefore the ten reported marginal intervals have a **simultaneous geometric
intersection** with the K5 cut space: one five-node scalar can place every edge
at an absolute residual no larger than $0.330$ times that edge's reported
marginal standard-uncertainty width. This is a covariance-limited deterministic
compatibility statement. It is not a joint 68% coverage statement, confidence
region, $p$-value, significance test, or probability assigned to the
common-node-scalar hypothesis.

## 5. Exact scientific result and ceiling

**Exact result (K5-11).** The ten processed pairwise redshift-derived heights
were not algebraically forced to close, and they are compatible with one
five-node scalar height vector inside every reported marginal $1\sigma$
interval. Under the source's already-assumed GR mapping and common $g$, that
height vector is proportional to a relative gravitational-potential vector;
the present calculation does not establish that physical interpretation
independently.

Because Eq. (K5-01) is one common scalar multiplication, zero-cycle structure
is the same for the inferred heights and their underlying processed pairwise
redshifts. But the mapping already assumes general relativity and a common local
$g$. The result is therefore a processed-output consistency check, **not** an
independent confirmation of GR, gravity, a common spacetime metric, or local
$g$. The unowned shared covariance also prevents a probabilistic compatibility
claim.

This lane does not test record formation, record lineage, beta_TM, RGRL,
Gravity Formation Theory, or gravity emergence. It does not own run-wise
ten-pair vectors, raw ellipse fits, covariance, a matched intervention, or a
new gravity measurement. No canonical MODEL, URM, or experiment register is
modified.
