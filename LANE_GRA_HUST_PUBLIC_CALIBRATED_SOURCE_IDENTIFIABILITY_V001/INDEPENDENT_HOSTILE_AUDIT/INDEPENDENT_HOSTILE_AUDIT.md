# Independent hostile audit of the HUST calibrated-source identifiability lane

**Audited lane:** `LANE_GRA_HUST_PUBLIC_CALIBRATED_SOURCE_IDENTIFIABILITY_V001`  
**Audit date:** 2026-08-27  
**Disposition:** `PASS__M1_N1_CORE_REPAIRS_INCORPORATED__PUBLICATION_SAFE_WITH_INTRINSIC_CEILINGS`  
**Builder/core files:** treated as frozen and not edited

## 1. Bottom line

The central correction signs, units, forward equations, calibrated partial
kernels, comparator gaps, and category-covariance arithmetic independently
recompute. No accepted, historical, recommended, or CODATA value of $G$ is
needed to reproduce any primary partial quotient. In particular, the primary
article prints the ToS fibre-anelastic corrections with explicit negative
signs, while the magnetic-damper corrections are positive; the frozen lane
maps them correctly.

The repaired core now incorporates the two qualifications raised by the first
hostile pass:

1. The minimal missing statistic is **one independently owned physical-harmonic
   remainder per released row**, not merely “one unpublished number.” The
   authors' processed coefficients already make a numerical remainder publicly
   inferable. They cannot be used as an independent source reconstruction
   because they are the authors' model output being cross-checked.
2. The ToS “deprocessing” is an exact algebraic inverse only under the lane's
   declared composition convention. The source owns the correction values and
   says they were applied synchronously, but it does not publish the raw
   correction operator or time series. This identity is not recovery of a raw
   numerator.

Neither qualification changes the reported conditional partial forwards. Both
are now explicit in `THEOREM.md`, `RESULT.json`, and the executable analyzer,
so the builder core is publication-safe unchanged within its stated ceilings.
This audit remains independent evidence for that disposition; it is no longer
needed as a semantic repair to the core.

## 2. Independent source custody

All three primary objects were freshly reacquired over HTTPS into a temporary
directory. Each was byte-for-byte identical to the frozen lane object:

| source | bytes | fresh and frozen SHA-256 |
|---|---:|---|
| official Springer Nature supplement | 2,711,453 | `5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb` |
| official Nature Table 1 HTML | 193,307 | `23436d4be7600a7a9dffa02cc4167a20b6eea032a181e77899bb57bb90aa02e9` |
| public university mirror of primary article | 7,999,835 | `40756ec0fb8f00c1fde31020b294521a3b220a196bef884a2ea5f3534d77dfaa` |

The primary article remains a nonpublisher mirror, as the builder disclosed.
The exact stable hash and publisher-produced DOI/title/tables support identity;
the host itself is not upgraded to publisher custody.

The seven sealed upstream theorem/result/executable dependencies also match
their pinned hashes.

## 3. Official correction signs, units, and mapping

An independent visual transcription of primary Extended Data Tables 3 and 5
and Supplementary Tables 2–4 was compared field-by-field with
`CALIBRATION_FIELDS.json`: **80/80 correction fields agree**. The official
Nature Table 1 HTML was parsed directly rather than copied from the builder:
**42/42 error-budget vectors agree**.

Extended Data Table 3 labels every entry as a correction to the value of $G$
in ppm. For ToS it gives:

| fibre | anelastic correction (ppm) | magnetic correction (ppm) | dynamic sum (ppm) |
|---|---:|---:|---:|
| F1 | -6.01 | +0.47 | -5.54 |
| F2 | -8.38 | +7.13 | -1.25 |
| F3 | -5.68 | +0.32 | -5.36 |
| F4 | -6.92 | +0.27 | -6.65 |

The sign is therefore source-owned, not inferred by fitting an authors' $G$
value. Supplementary Section 1 places the anelastic and magnetic terms in the
same correction bracket, so the frozen forward

\[
G_i=\frac{\Delta\omega_i^2
\left[1+(c_{i,\mathrm{anel}}+c_{i,\mathrm{mag}})10^{-6}\right]}
{K_i^{\mathrm{partial}}+r_i}
\]

has the correct sign and units.

For a correction $c_i^{\rm mass}$ tabulated as a quotient-level correction to
$G$,

\[
G_{\rm corrected}=G_{\rm hom}(1+c_i^{\rm mass})
=\frac{Y_i}{K_i^{\rm hom}/(1+c_i^{\rm mass})}.
\]

Thus $K_i^{\rm partial}=K_i^{\rm hom}/(1+c_i^{\rm mass})$ is the correct
quotient-level mapping. It avoids the forbidden operation of adding mass only
to an inertia denominator without its numerator. The limitation is physical
ownership: these component corrections are still authors-model-mediated
calibrations transplanted onto an independently integrated homogeneous
functional. The output is correctly called a *calibrated partial model*, not a
raw-map reconstruction.

## 4. Independent numerical forwards and corrected gaps

The four mass corrections sum independently to:

- AAF: 22.40, 22.89, and 22.29 ppm;
- ToS: 90.33, 90.33, 92.05, 91.24, 91.24, 88.70, and 88.70 ppm.

Using only the pinned homogeneous kernels, primary-table response summaries,
and independently transcribed correction entries reproduces every reported
partial kernel and primary quotient:

| row | independently recomputed $K^{\rm partial}$ (kg m$^{-3}$) | independently recomputed conditional $G$ (SI) | processed-minus-partial gap (ppm) |
|---|---:|---:|---:|
| AAF-I | 6926.505285 | 6.674385095e-11 | -22.1302 |
| AAF-II | 6926.541459 | 6.674175468e-11 | -29.9513 |
| AAF-III | 6926.545615 | 6.674409224e-11 | -18.8572 |
| ToS F1 first | 24912.933578 | 6.674134876e-11 | -2.9534 |
| ToS F1 repeat | 24912.178731 | 6.674204640e-11 | -2.3575 |
| ToS F2 | 24912.135886 | 6.674240737e-11 | +0.5666 |
| ToS F3 first | 24911.970555 | 6.674201402e-11 | -10.8604 |
| ToS F3 repeat | 24912.021729 | 6.674183678e-11 | -12.1118 |
| ToS F4 first | 25003.124942 | 6.673997378e-11 | -2.9973 |
| ToS F4 repeat | 25003.041260 | 6.674079707e-11 | -3.6499 |

The many displayed digits are computational replay digits. They are not
physical significant figures independently owned beyond the resolution and
uncertainty of the public tables.

## 5. Processed response versus raw response

### AAF

Supplementary Table 3 explicitly says the published $\alpha_t$ has already
been corrected for air density. Supplementary Section 4 gives the averaging and
twice-differentiated sinc attenuation formulae. Independent evaluation gives:

- AAF-I: 2.570214 ppm averaging and 2058.706957 ppm differentiation;
- AAF-II/III: 1.142316 ppm averaging and 914.353638 ppm differentiation.

These reproduce the displayed 2.57/2058.71 and 1.14/914.35 ppm entries. Dividing
both $\alpha_t$ and its convention-matched kernel by the same transfer factor
necessarily preserves the quotient. That is a valid algebraic custody check,
not an independent fit to encoder samples.

### ToS

The article states that thermoelastic, fibre-nonlinearity, and source-gravity
nonlinearity were corrected synchronously in determining
$\Delta\omega^2$. The frozen lane chooses the product
$(1+c_{\rm thermal})(1+c_{\rm grav})$, divides by it, and multiplies it back.
That identity is exact by construction. The source does not publish enough raw
processing detail to prove that this product is the unique historical operator.

For scale, using an additive rather than multiplicative ppm composition changes
the chosen factor by at most 0.010084 ppm across the seven rows. It does not
enter the primary forward, which starts from the already corrected published
response. The issue is epistemic wording, not a numerical defect in the
reported partial $G$ values.

The source-gravity nonlinearity correction is itself calculated from the
authors' apparatus potential. Consequently even algebraically “deprocessed”
ToS data are not a source-model-free numerator.

## 6. Accepted-$G$ and normalization attack

An entirely separate reconstruction produced the partial forwards without any
accepted, CODATA, historical, or authors-derived $G$ number. The only inputs
were response summaries, homogeneous kernels, and correction factors. Authors'
campaign $G$ values and processed kernels were accessed only afterward to
form comparator residuals.

The primary source PDFs naturally mention CODATA and authors' final values, and
the loaded conditional JSON contains comparator fields. Mere presence is not a
numerical dependency. The independent reconstruction demonstrates that the
same primary rows are obtained with no such value.

No normalization collision was found in the frozen mapping. However, applying
the correction scalars does not turn them into independent density/CAD data;
the calibrated partial source remains hybrid and conditional.

## 7. Category covariance arithmetic

Starting directly from official Table 1:

- AAF statistical terms were placed on the diagonal and each non-statistical
  category was treated as fully correlated across campaigns, exactly as the
  Supplement specifies. Published-total inverse-variance weights give
  **11.616239 ppm**, reproducing the reported 11.61 ppm after table rounding.
- For ToS, cross-fibre same-category systematics were treated as fully
  correlated. The covariance of each repeated same-fibre pair was inferred
  from its two rounded row totals and rounded combined-fibre uncertainty. The
  resulting hierarchy gives **11.637465 ppm**, reproducing 11.64 ppm.

All independently reconstructed AAF, ToS run-level, and ToS fibre-level
matrices match the frozen matrices entry-by-entry and are positive definite.

This does not make the covariance empirical or unique. The same-fibre shared
background components are inferred from rounded aggregate outputs, and the
100% correlations are authors' combination assumptions. These matrices support
an arithmetic reproduction of the paper's standard uncertainty, not a raw-data
covariance, confidence interval, or coverage theorem.

## 8. Minimal remainder theorem: what is proved

At one released harmonic, write

\[
Y_i=G\left(K_i^{\rm partial}+r_i\right).
\]

Once the released processed central response $Y_i$ is treated as fixed, one
scalar $r_i$ is sufficient to evaluate that row. It is necessary in the
following precise sense: if two admissible undisclosed maps agree on every
independently released field but have different harmonic kernels, a rule using
only those released fields cannot point-identify both. Some statistic that
distinguishes their scalar harmonic values is required.

Three boundaries matter:

1. **Ownership.** Supplementary Tables 2 and 3 already publish the authors'
   processed kernels, so the numerical comparator
   $r_i^{\rm author}=K_i^{\rm processed}-K_i^{\rm partial}$ is public. What is
   missing is an independently owned and reproducible physical-harmonic value.
   Calling the scalar simply “unpublished” would be false.
2. **Joint dimension.** Ten row values are sufficient coordinates. The theorem
   does not prove ten independent physical degrees of freedom; shared maps may
   correlate or constrain the ten values.
3. **Raw response.** One $r_i$ is not sufficient for a source-model-free
   reanalysis. Raw samples, event/correction custody, the design matrix, and
   joint covariance are additional required objects.

With those qualifications, the fixed-harmonic minimal-statistic theorem holds.

## 9. Can the published standard uncertainties make a compact domain?

There are two different questions:

- If one accepts the authors' processed coefficient and a chosen probabilistic
  interpretation, its central value plus/minus a standard uncertainty is a
  finite conventional display band. Mapping such a band through the monotone
  quotient gives a compact conventional $G$ band.
- That band is **not** a deterministic admissible set and is not an independent
  cross-check. Standard uncertainties are not hard support bounds; the public
  release supplies neither a joint distribution for the missing maps nor an
  independently reconstructed central $r_i$.

Therefore the frozen theorem's “no compact domain” statement is correct only
when read as **no independently owned deterministic identified domain or
coverage-certified interval**. It must not be paraphrased as “no finite
authors-model uncertainty display can be written.”

## 10. Repair closure and publication disposition

### M1 — closed in the repaired core

The scalar remainder is numerically public through the authors' processed
coefficient. The missing object is the same scalar with independent physical
custody. The repaired theorem now says this directly, renames its proof
“Minimum independently owned remainder proof,” and distinguishes sufficient row
coordinates from independent physical degrees of freedom. The executable schema
is now `minimal_independent_remainder_theorem`; the obsolete
`minimal_unreported_parameter_theorem` key is absent. M1 is therefore
incorporated in core language and executable output.

### N1 — closed in the repaired core

ToS invert/reapply is exact under the chosen product convention, not an
independent reversal of the historical raw processing operator. The repaired
theorem and executable strict ceiling now say both parts explicitly. This does
not affect the primary corrected-response quotient, and N1 is incorporated in
core language and executable output.

### N2 — computational digits

Extra digits are useful for replay but are not measured significant figures.

### Final disposition

No material arithmetic, sign, unit, custody, accepted-$G$ leakage,
normalization, or covariance-closure defect was found. The lane passes as a
**calibrated public partial Newtonian source family and fixed-harmonic
identifiability theorem**. M1 and N1 are closed in the repaired frozen core.
`THEOREM.md`, `RESULT.json`, and the analyzer are publication-safe unchanged;
their intrinsic ceilings still apply. The result is not a new measurement of
$G$, a raw source reconstruction, a deterministic interval, a coverage
theorem, or evidence for GR, RGRL, GFT, record lineage, gravity emergence, or a
common metric.
