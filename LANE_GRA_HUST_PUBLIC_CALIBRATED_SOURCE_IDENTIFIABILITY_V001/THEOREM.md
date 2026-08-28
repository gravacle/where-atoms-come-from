# HUST public calibrated source-identifiability theorem V001

Date: 2026-08-27  
Status: **PASS — strongest public calibrated partial source model; point identification remains one normalized-kernel scalar per released row**

## Question

How much of the HUST-2018 physical source and response model can be reconstructed from public primary materials, without importing an accepted or CODATA value of \(G\), and what is the exact smallest missing object if a fully independent point value still cannot be formed?

## Inputs and prohibition

This lane hash-pins the sealed homogeneous-kernel and conditional-forward lanes, the official Nature supplement and error-budget table, and a public mirror of the publisher-produced primary article PDF. The primary calculation uses:

1. independently integrated homogeneous Newtonian kernels \(K_i^{\rm hom}\);
2. published response summaries \(y_i\) (AAF \(\alpha_t\), ToS \(\Delta\omega^2\));
3. the central correction operators and one-standard-deviation uncertainties published in Extended Data Tables 3 and 5 and Supplementary Tables 1–4.

No accepted, recommended, historical or CODATA value of \(G\) enters a kernel, correction, fit or quotient. Authors' processed source coefficients and campaign \(G\) values are attached only after the primary calculations as quarantined comparators.

## The calibrated physical source theorem

Let \(c_{ij}^{\rm mass}\) denote the published correction to \(G\), in dimensionless units, for the coating, clamp, ferrule and other explicitly listed pendulum bodies/maps omitted from the ideal homogeneous block. Under the paper's displayed correction convention, define

\[
K_i^{\rm partial}
=
\frac{K_i^{\rm hom}}
{1+\sum_j c_{ij}^{\rm mass}}.
\]

This mapping preserves numerator/inertia ownership: it converts the ideal-body quotient into the publicly calibrated physical-body quotient rather than inserting missing inertia in the denominator alone.

The public corrections sum to 22.29–22.89 ppm for AAF and 88.70–92.05 ppm for ToS. The resulting independently normalized partial kernels are:

- AAF: 6926.5053, 6926.5415 and 6926.5456 kg m\(^{-3}\);
- ToS: 24912.9336, 24912.1787, 24912.1359, 24911.9706, 24912.0217, 25003.1249 and 25003.0413 kg m\(^{-3}\).

The ToS Extended Data Table 3 also publishes the previously missing **signed** anelastic corrections: \(-6.01\), \(-8.38\), \(-5.68\) and \(-6.92\) ppm by fibre. Therefore the ToS central forward is no longer a free \(c_f\) line. The public partial forward families are

\[
G_i^{\rm AAF}(r_i)
=
\frac{\alpha_i f_{m,i}}
{K_i^{\rm partial}+r_i},
\qquad
G_i^{\rm ToS}(r_i)
=
\frac{\Delta\omega_i^2\,[1+c_{i}^{\rm anel}+c_i^{\rm mag}]}
{K_i^{\rm partial}+r_i}.
\]

Here \(r_i\) is the normalized harmonic remainder of the physical apparatus.
The authors' processed coefficient makes one comparator value numerically
inferable, but the public packet does not independently reconstruct that value
from released physical maps or supply an independently owned deterministic
admissible set.

At \(r_i=0\), the AAF partial forwards are \(6.6743851\), \(6.6741755\) and \(6.6744092\) \(\times10^{-11}\) m\(^3\) kg\(^{-1}\) s\(^{-2}\). The seven ToS partial forwards range from \(6.6739974\) to \(6.6742407\) \(\times10^{-11}\) m\(^3\) kg\(^{-1}\) s\(^{-2}\). These are conditional public-source calculations, not adopted values of \(G\).

## Transfer-function closure

For AAF, the public air-density, half-second averaging and ten-second numerical-derivative corrections can be inverted. The directly recomputed acquisition corrections are 2.570214 and 2058.706957 ppm for AAF-I and 1.142316 and 914.353638 ppm for AAF-II/III, reproducing the displayed 2.57/2058.71 and 1.14/914.35 ppm values. Writing

\[
T_i=(1+c_i^{\rm air})(1+c_i^{\rm avg})(1+c_i^{\rm diff}),
\]

the deprocessed convention \(\alpha_i/T_i\) over \(K_i^{\rm partial}/T_i\) gives exactly the same forward quotient. This is an algebraic custody check, not a recreation of raw encoder data.

For ToS, the published thermoelastic and gravitational-nonlinearity corrections can likewise be inverted and reapplied exactly under the declared multiplicative composition convention. This is an algebraic identity, not proof that the convention uniquely reconstructs the historical raw-processing operator. The latter correction remains authors-model-mediated because it was calculated from the source gravitational potential; the public packet does not contain campaign raw time series for a source-model-free refit.

## Minimum independently owned remainder proof

At the released harmonic, every undisclosed spatial mass field enters the forward law only through the scalar linear functional \(K_i[\rho]\). After all disclosed components have been evaluated, their difference is one scalar \(r_i\). Thus:

- one independently owned \(r_i\) per released row is **sufficient** for row-wise point evaluation;
- at least one distinguishing scalar is **necessary** whenever two admissible undisclosed maps produce different harmonic kernels but identical independently released fields.

Supplementary Tables 2 and 3 already publish the authors' processed kernels,
so an authors-model comparator \(r_i\) is numerically inferable. It is not an
independent source reconstruction because it is the model output being checked.
The ten row values are sufficient coordinates; shared maps may correlate or
constrain them, so the theorem does not assert ten independent physical degrees
of freedom.

The full measured density maps, run-specific sphere orientations, individual 3D CMM coordinates (including the two AAF shear coordinates), attachment/coating maps, and AAF shelf/deformation/compensation maps are sufficient physical data to calculate these remainders, but they are not public. For a genuinely independent raw-response reanalysis, campaign-bound raw samples, the correction ledger and raw design/covariance matrix are additionally required.

Quoted standard uncertainties and isolated upper bounds do not supply a joint deterministic domain for \(r_i\). Therefore the public packet owns no independently reconstructed deterministic or coverage-certified compact numerical \(G\) interval for these families. An authors-model conventional display band can be formed under additional assumptions; this is not an assertion that the physical remainder is mathematically unbounded.

## Comparator diagnostic — not an input

Only after the partial kernels are fixed, comparison with the authors' processed coefficients shows residual normalized-kernel differences of:

- AAF: \(-18.86\) to \(-29.95\) ppm, versus \(-41.15\) to \(-52.84\) ppm before calibration;
- ToS: \(+0.57\) to \(-12.11\) ppm, versus \(-91.48\) to \(-103.35\) ppm before calibration.

This demonstrates that the public central physical corrections account for most of the ToS homogeneous-to-processed kernel gap and roughly half of the AAF gap. It is a post-calculation localization diagnostic, not independent validation of \(K_i^{\rm partial}\).

## Covariance result

The paper's category-level correlation rules reconstruct the AAF standard-uncertainty covariance and reproduce a combined 11.616 ppm from the rounded table entries (reported 11.61 ppm). The analogous hierarchical ToS reconstruction gives 11.637 ppm (reported 11.64 ppm). The same-fibre shared-background components can be inferred only from rounded totals and combined uncertainties. These matrices encode the authors' uncertainty-combination assumptions; they are not empirical raw-data covariance or coverage theorems.

## Strict ceiling

This theorem supplies a calibrated public partial Newtonian source model and the exact public-data identifiability class. It does **not** supply a full finite-element apparatus reconstruction, a source-model-free raw numerator reanalysis, an independently reconstructed deterministic or coverage-certified compact interval, or independent evidence for GR, RGRL, record lineage, \(\beta_{TM}\), gravity emergence, non-Newtonian gravity or a common metric.
