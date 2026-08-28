# HUST ToS equal-configuration round-trip history theorem

**Theorem ID:** `GRA-HUST18-TOS-RTHR-V001`  
**Date:** 2026-08-27  
**Claim class:** exact extraction from official figure-level period summaries;
exact linear weight/rank/endpoint-overlap theorem; descriptive apparatus-history
diagnostic with no covariance or lineage promotion

**Status:**
`EQUAL_CONFIGURATION_ENDPOINT_RETURNS_EXTRACTED__PRESENT_MINUS_BACKGROUND_ORIENTATIONS_SEPARATED__EXACT_REUSED_ENDPOINT_WEIGHT_AND_OVERLAP_MATRICES__DESCRIPTIVE_ONLY_NO_LINEAGE_OR_COVERAGE_CLAIM`

## 1. Frozen source and question

The sole empirical source is the official HUST Figure-2 workbook already pinned
in the audited dual-method forward lane. Sheets `a` and `b`, cells `B3:D22`,
contain ten near and ten far three-day period summaries with the source masses
present and absent, respectively. They are raw-like figure intermediates, not
the original angle stream.

The question is deliberately narrower than the published gravity contrast. If
the schedule goes near--far--near (`N-F-N`) or far--near--far (`F-N-F`), how
much does the transformed response change between the two equal-configuration
endpoints, and how does that return change differ between the source-present
and source-absent panels?

This is an apparatus/source-history diagnostic. The schedule does not contain
authenticated record formation, sham formation, KEEP/BREAK, or independently
randomized lineage. No lineage label is assigned retrospectively.

## 2. Exact endpoint-return observable

For each printed period define

\[
 y=(2\pi/T)^2\quad[\mathrm{s}^{-2}].                    \tag{HR01}
\]

For panel \(p\in\{P,B\}\), where `P` is source-present and `B` is the
source-absent background, let \(y^p_{N,i}\) and \(y^p_{F,i}\) be the ten
chronologically ordered equal-configuration endpoint series. Exact workbook
times verify that \(F_i\) lies between \(N_i,N_{i+1}\), and \(N_{i+1}\) lies
between \(F_i,F_{i+1}\). Define

\[
 r^p_{N,i}=y^p_{N,i+1}-y^p_{N,i},\qquad
 r^p_{F,i}=y^p_{F,i+1}-y^p_{F,i},\qquad i=0,\ldots,8. \tag{HR02}
\]

The opposite-configuration middle point certifies the round trip but has zero
weight in this equal-endpoint observable. The panel differentials are

\[
 d_{N,i}=r^P_{N,i}-r^B_{N,i},\qquad
 d_{F,i}=r^P_{F,i}-r^B_{F,i}.                         \tag{HR03}
\]

Here the shared index \(i\) denotes the same sequence ordinal in two separate
panels. It is not evidence that the two loops were simultaneous, randomized, or
an authenticated matched experimental pair. Reordering the loopwise pairing
would leave the two differential means unchanged, but it would change the
individual \(d_i\), \(c_i\), and \(h_i\) values and the RMS reported below.

For descriptive orientation separation only, define

\[
 c_i={d_{N,i}+d_{F,i}\over2},\qquad
 h_i={d_{N,i}-d_{F,i}\over2}.                         \tag{HR04}
\]

The names `c` and `h` do not assign mechanisms. In particular, \(h_i\) is not a
record-history charge and \(c_i\) is not proven to be drift.

Every endpoint return follows an intervening opposite-configuration excursion,
but the release contains no matched equal-duration **no-excursion** trajectory.
Consequently, the return cannot identify hysteresis or memory against ordinary
time drift. It is a history-confound diagnostic, not a causal history effect.

## 3. Exact linear and endpoint-overlap theorem

Order the transformed endpoint vector as

\[
 x=(P_N^{10},P_F^{10},B_N^{10},B_F^{10})\in\mathbb R^{40}. \tag{HR05}
\]

Let \(D\in\mathbb Z^{9\times10}\) be the first-difference matrix

\[
 D_{ij}=-\delta_{ij}+\delta_{i+1,j}.                   \tag{HR06}
\]

Then all 36 panel returns and all 18 panel differentials are

\[
 R x,\quad R=\operatorname{diag}(D,D,D,D),             \tag{HR07}
\]

\[
 W x,\quad
 W=\begin{pmatrix}D&0&-D&0\\0&D&0&-D\end{pmatrix}.   \tag{HR08}
\]

Exact rational elimination gives

\[
 \operatorname{rank}R=36,\qquad \operatorname{rank}W=18. \tag{HR09}
\]

With \(L=DD^{\mathsf T}\), \(L\) has diagonal `2`, first off-diagonal `-1`,
and all other entries zero. Therefore

\[
 RR^{\mathsf T}=\operatorname{diag}(L,L,L,L),\qquad
 WW^{\mathsf T}=\operatorname{diag}(2L,2L).            \tag{HR10}
\]

These exact integer Gram matrices expose endpoint reuse: adjacent returns share
one endpoint with the opposite sign. They are **not** empirical covariance
matrices. The source releases no row covariance from which such a covariance
could be calculated.

The arithmetic mean of each nine-return orientation telescopes:

\[
 {1\over9}\sum_{i=0}^{8}r^p_{c,i}
 ={y^p_{c,9}-y^p_{c,0}\over9}.                         \tag{HR11}
\]

Thus nine overlapping returns do not constitute nine independent replications
of their mean.

## 4. Exact descriptive extraction

In (mathrm{s}^{-2}), the four return means are

| Panel/orientation | Mean return |
|---|---:|
| source-present `N-F-N` | \(3.647457758732247\times10^{-10}\) |
| source-present `F-N-F` | \(3.496735157677661\times10^{-10}\) |
| source-absent background `N-F-N` | \(5.498090526271760\times10^{-10}\) |
| source-absent background `F-N-F` | \(5.324953327068938\times10^{-10}\) |

Consequently,

\[
 \bar d_N=-1.8506327675395133\times10^{-10}\ \mathrm{s}^{-2},
 \qquad
 \bar d_F=-1.8282181693912770\times10^{-10}\ \mathrm{s}^{-2}, \tag{HR12}
\]

and

\[
 \bar c=-1.8394254684653953\times10^{-10}\ \mathrm{s}^{-2},
 \qquad
 \bar h=-1.1207299074118054\times10^{-12}\ \mathrm{s}^{-2}. \tag{HR13}
\]

The full orientation-mean difference is
\(-2.241459814823611\times10^{-12}\,\mathrm{s}^{-2}\), while the descriptive
RMS of the nine \(h_i\) values is
\(6.524819132094064\times10^{-11}\,\mathrm{s}^{-2}\).

These digits are retained so that a byte-for-byte numerical replay can be
checked. The workbook period cells lie on a nominal \(10^{-5}\,\mathrm{s}\)
grid; digits beyond that source resolution are computational digits, not measurement
precision or an uncertainty statement. In addition, the loopwise \(h_i\) RMS
depends on the explicitly declared ordinal cross-panel pairing; it is not the
RMS of matched experimental trials.

Both orientation means therefore share a much larger common panel difference
than their mean orientation difference at this released resolution. That is a
numerical description, not evidence that the orientation-odd term is zero or
that the common term has any particular cause.

## 5. Scientific ceiling

No standard error, confidence interval, \(p\)-value, coverage statement, or
physical null is assigned. The hard reasons are:

1. adjacent returns reuse endpoints exactly as in (HR10);
2. the workbook supplies no row covariance or event-level angle stream;
3. the source-present and background panels are separate acquisitions, and no
   matched no-excursion control separates excursion history from ordinary time
   drift;
4. printed period precision and unowned fibre, mechanical, thermal, controller,
   and source-motion histories remain in the observable; and
5. same-ordinal loop durations are not identical—the largest printed `N-F-N`
   present/background duration difference is `0.4445027 day`—and the primary
   residual is intentionally not time-normalized; the ordinal correspondence
   is not an authenticated experimental match; and
6. the source-present/source-absent meaning is inherited from the pinned
   article-level custody; the worksheet headers themselves encode time, near,
   and far, but do not independently encode panel matching.

Therefore the strongest claim is

\[
 \boxed{\text{official HUST summaries support a reproducible equal-configuration
 return diagnostic and exact endpoint-overlap ledger}.}       \tag{HR14}
\]

This packet does **not** estimate `beta_TM`, instantiate a record-lineage
intervention, identify causal memory, demonstrate a lineage gravitational
charge, test Gravity Formation Theory, establish gravity emergence, execute
full `GC16`, or measure a new \(G\).
