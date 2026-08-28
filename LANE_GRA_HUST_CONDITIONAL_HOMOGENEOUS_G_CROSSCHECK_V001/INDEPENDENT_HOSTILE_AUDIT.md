# Independent hostile audit: HUST conditional homogeneous-kernel G cross-check

**Audit ID:** GRA-HUST-CHGC-V001-IHA

**Date:** 2026-08-27

**Disposition:**
`ACCEPT_CONDITIONAL_QUOTIENT_AND_PUBLIC_PACKET_NONIDENTIFICATION__NOT_NEW_G__NOT_GFT_EVIDENCE`

## Independence and method

This audit was performed by an agent that did not construct CHGC. I treated
`THEOREM.md`, `RESULT.json`, both lane executables, and all documentary claims
as untrusted. I independently reparsed the two byte-pinned parent results,
re-extracted the official Supplement with Poppler, visually inspected the
load-bearing PDF equations and tables, reconstructed the correction factors
from Supplementary Table 1, and recalculated all three AAF quotients, all seven
ToS anchors and slopes, every published local uncertainty component, both
mixed-normalization diagnostics, all post-comparator identities, and both ToS
figure-level anchors. I also inspected the calculator AST rather than trusting
its prose about comparator quarantine.

The executable replay is
`hostile_audit_conditional_homogeneous_g.py`; its sealed transcript records
`255/255` hostile checks passed after final resealing.

## Material defect found and repaired

The original lane used mechanical correction factors inherited from the
forward parent with digits reconstructed from component fields. The official
source directly owns only the displayed Supplementary-Table-1 corrections:
AAF-I/II \(455.40(1.95)\) ppm, AAF-III \(25.74(0.08)\) ppm, and ToS fibres 1--4
\(0.47(0.08)\), \(7.13(1.19)\), \(0.32(0.05)\), and \(0.27(0.08)\) ppm. For
example, the inherited AAF-I/II factor was \(1.000455396097372\), while the
directly displayed factor is \(1.00045540\). The extra digits were not
source-owned precision.

I repaired the primary calculation to use the displayed central corrections
and their displayed one-standard-deviation uncertainties. This changes the
last few digits of the conditional outputs but no scientific conclusion. I
also separated the local processed-kernel forward evaluated with the displayed
factor from the upstream component-reconstructed forward, so the exact
kernel-ratio identity no longer mixes incompatible correction precisions.

## Formula and sign findings

Visual inspection of Supplement pages 3--5 confirms the equations used by the
lane:

\[
G_{\mathrm{AAF}}=\frac{\alpha_t}{K_{\mathrm{full}}}
\left(1+\frac{K}{K_m}\frac{I_m}{I}\right),
\]

and

\[
G_{\mathrm{ToS}}=\frac{\Delta\omega^2}{K_{\mathrm{full}}}
\left(1-\frac{\Delta K}{I\Delta\omega^2}
+\frac{I_mK^2}{IK_m^2}\right).
\]

Thus `c_f=-Delta_K/(I*Delta_omega2)` and the lane's bracket
\(1+\delta_m+c_f\) has the correct signs. Reconstructing the positive
magnetic factors from Table-1 component central values lands within half of the
last displayed \(0.01\)-ppm unit for every campaign/fibre. Units also close:
response divided by a \(\mathrm{kg\,m^{-3}}\) normalized kernel has SI units
\(\mathrm{m^3\,kg^{-1}\,s^{-2}}\).

## Numerator and processing semantics

The released-file response replays are computationally non-circular with
respect to an accepted \(G\) and the processed source coefficient. They are not
statistically or model independent numerators. The official Supplement states
that source-mass gravitational nonlinearity is corrected synchronously in the
ToS \(\Delta\omega^2\), that AAF campaign \(\alpha_t\) is corrected for air
density, and that the released one-second AAF stream follows 20-kHz acquisition
and averaging. The repaired theorem and machine result now say
"non-circularity," not "independence."

The AAF figure harmonic remains unbound to a campaign-average kernel and is
correctly withheld from a fourth quotient. The ToS figure diagnostic inherits
the parent lane's repeated-fibre-1 binding; CHGC does not claim to establish a
new campaign binding from the figure alone.

## Comparator quarantine

The complete pinned parent JSON objects must be parsed before primary
calculation, so byte-level language such as "opened only after" was too strong.
The audited property is key-level nonuse: AST traversal finds none of
`processed_coefficient_kg_m-3`,
`published_G_summary_SI_comparison_only`, or `recomputed_G_SI` in primary
extraction or computation. Comparator keys are selected only in
`attach_post_comparators()`, which is called after the primary object is
complete. No accepted or CODATA \(G\) field or payload is present.

## Uncertainty and joint-coverage ceiling

For every row I independently verified

\[
u_{G,N}=|G|u_N/|N|,\qquad
u_{G,K}=|G|u_K/|K|,\qquad
u_{G,f}=|G|u_f/|f|,
\]

and the labelled zero-covariance partial RSS. The magnetic uncertainty is an
absolute uncertainty on the dimensionless correction, so the factor component
\(|G|u_f/|f|\) is correct. Repeated fibres, shared clocks, shared campaigns,
and shared processing can induce covariance that the public packet does not
supply. These per-row marginal diagnostics are therefore neither joint
statistical coverage nor physical-\(G\) intervals; the lane forms no combined
cross-row \(G\).

## Normalization and interval claims

The quantity called `normalization_collision` in the JSON is the deliberately
incomplete arithmetic mix of a homogeneous-core numerator with full reported
inertia. It demonstrates the numerical consequence of unmatched normalization;
it does not construct a second physically realized HUST apparatus map. The
load-bearing nonidentification is structural: scalar inertia and the relevant
\(m=2\) gravitational numerator are different functionals, and the public
packet supplies neither the missing density/CMM map nor a law relating them.

The defensible conclusion is correspondingly narrow. The public packet owns no
independently identified deterministic compact interval for the nominal-kernel
conditional family because it gives no bounded admissible domain or covariance
law for \(r_{\mathrm{norm}}\) and gives no signed row-level ToS \(c_f\). This is not
a theorem of mathematical unboundedness under every external apparatus domain,
and it does not deny the authors' published processed-model \(G\) summaries and
standard uncertainties. Those remain post-calculation comparators.

## Final claim ceiling

After repair, CHGC establishes deterministic conditional homogeneous-kernel
quotients, exact algebraic ToS affine families, correctly scoped local
uncertainty propagation, and public-packet nonidentification of the matched
full-kernel family. It does not establish a new or independent \(G\), a compact
physical coverage interval, GC16, lineage, beta_TM, RGRL, GR, a common metric,
or Gravity Formation Theory/emergence evidence.
