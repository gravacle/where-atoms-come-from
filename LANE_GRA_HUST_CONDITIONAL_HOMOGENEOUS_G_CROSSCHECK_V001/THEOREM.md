# HUST conditional homogeneous-kernel G cross-check and identifiability theorem

**Theorem ID:** GRA-HUST-CHGC-V001

**Date:** 2026-08-27

**Claim class:** deterministic high-precision conditional response/kernel
quotients at serialized public precision; exact algebraic affine ToS correction
family; first-order propagation of only released response, nominal-kernel, and
published mechanical-correction uncertainties; mixed-normalization diagnostic;
source-model identifiability result; strict key-level post-calculation
comparator quarantine.

**Status:**
AAF_RNORM0_CONDITIONAL_QUOTIENTS_COMPUTED__TOS_RNORM0_CF0_ANCHORS_AND_AFFINE_CF_FAMILIES_COMPUTED__FIGURE_LEVEL_RESPONSE_NUMERATOR_CUSTODY_SEPARATED_FROM_CORRECTED_CAMPAIGN_SUMMARIES__NO_PUBLIC_INDEPENDENT_DETERMINISTIC_G_INTERVAL__MATCHED_FULL_NUMERATOR_AND_TOS_DELTAK_REMAIN_MISSING__NO_NEW_G

**Not claimed:** a new or independent measurement of \(G\); an independently
identified physical point or deterministic compact interval for \(G\) from the
nominal-kernel family (the authors' processed-model summaries are acknowledged);
that the real HUST pendulum equals the homogeneous
cuboid; that the AAF placement is uniquely public; that the Table-2 ToS
response is source-model free; that the Table-3 AAF response is raw; that the
processed source coefficients are independent inputs; full GC16; RGRL or
Gravity Formation Theory confirmation; or a record-lineage gravitational
charge.

## 1. Exact question and pinned parents

The audited HUST forward lane proved that the official release contains real
response information in two ownership channels: ToS stiffness response and
AAF forcing response. The separately audited nominal-kernel lane then
reconstructed conditional homogeneous source coefficients without a measured
\(G\), response value, or processed source coefficient. This lane asks the
narrow remaining no-lab question:

> What numerical \(G\)-quotients follow if the independently calculated
> homogeneous normalized kernel is conditionally substituted for the unknown
> full-apparatus kernel, and exactly what remains unidentified?

The load-bearing dependencies are byte-pinned in DEPENDENCIES.sha256. No
accepted or CODATA \(G\) value is imported. The primary calculation selects only
nominal kernels and public-input sensitivities from NSKR, and response
summaries plus displayed mechanical corrections from DMGF and the official
Supplement. The pinned parent JSON objects are parsed before calculation, but
the authors' processed-source and derived-\(G\) keys are not selected until
after every primary quotient, affine family, uncertainty component, and
mixed-normalization diagnostic is computed. This is key-level computational
nonuse, not statistical independence of the response and source models.

## 2. Full-apparatus remainder parameterization

Let the full normalized source coefficient be

\[
 K_{\mathrm{full}}=K_{\mathrm{hom}}+r_{\mathrm{norm}}.       \tag{C01}
\]

Here \(K_{\mathrm{hom}}\) is the pairwise-centred homogeneous functional from
NSKR. The scalar \(r_{\mathrm{norm}}\) includes the normalized effect of
clamps, coatings, ferrules, density multipoles, unreported CMM or shear
coordinates, shape, and their matched inertia. The unprocessed public mass
fields do not identify or bound its central value.

For AAF the exact displayed response equation becomes

\[
 G(r_{\mathrm{norm}})
 =\frac{\alpha_t f_m}{K_{\mathrm{hom}}+r_{\mathrm{norm}}},
 \qquad
 f_m=1+\frac{K}{K_m}\frac{I_m}{I}.                           \tag{C02}
\]

For ToS define

\[
 c_f=-\frac{\Delta K}{I\Delta\omega^2},
 \qquad
 f_T=1+\delta_m+c_f,
 \qquad
 \delta_m=\frac{I_mK^2}{I K_m^2}.                            \tag{C03}
\]

Then

\[
 G(r_{\mathrm{norm}},c_f)
 =\frac{\Delta\omega^2}{K_{\mathrm{hom}}+r_{\mathrm{norm}}}
   (1+\delta_m+c_f).                                         \tag{C04}
\]

Equations (C02)--(C04), rather than any one numerical substitution, are the
identified conditional families. Setting \(r_{\mathrm{norm}}=0\) is an
explicit homogeneous normalized-kernel hypothesis, not a conclusion from
public data. For ToS, setting \(c_f=0\) is a second explicit anchor convention
because the pinned public fields do not supply the signed row-level
\(\Delta K\).

## 3. AAF conditional quotients

Using the three published campaign-average angular accelerations, which are
response quantities distinct from the processed source coefficients, and the
displayed magnetic-damper factors gives:

| campaign | \(\alpha_t\) (nrad s\(^{-2}\)) | \(K_{\mathrm{hom}}\) (kg m\(^{-3}\)) | \(G_{\mathrm{hom}}^{(0)}=G(r_{\mathrm{norm}}=0)\) (SI) | partial zero-covariance RSS |
|---|---:|---:|---:|---:|
| AAF-I | 462.0912 | 6926.660438859 | \(6.674235591786\times10^{-11}\) | \(7.6127367\times10^{-16}\) |
| AAF-II | 462.0791 | 6926.700007763 | \(6.674022699179\times10^{-11}\) | \(7.4547848\times10^{-16}\) |
| AAF-III | 462.2941 | 6926.700007763 | \(6.674260454935\times10^{-11}\) | \(7.1857101\times10^{-16}\) |

These are deterministic 50-digit decimal evaluations of (C02), rounded to the
serialized output precision, under \(r_{\mathrm{norm}}=0\). The equations are
exact; the displayed decimal outputs are not being relabelled exact real
numbers. They are not estimates of the full-apparatus \(G\), because the
hypothesis fixes precisely the unmeasured central remainder.

The AAF campaign \(\alpha_t\) values are not raw encoder outcomes:
Supplementary Table 3 says they have already been corrected for the
air-density effect. The official workbook supplies a more direct two-hour,
one-second response stream from which the forward lane extracts
\(461.993464795\ {\mathrm{nrad\,s^{-2}}}\) without using \(G\) or a source
coefficient. That segment is representative, not a campaign-average
\(\alpha_t\), and is not lawfully bound to one of the three campaign kernels.
No fourth \(G\)-quotient is formed from it.

### Theorem CHGC-1 -- AAF conditional computability

At displayed public precision, the AAF response summaries and independently
reconstructed homogeneous kernels determine the three table entries
conditional on \(r_{\mathrm{norm}}=0\). They do not identify \(G\) when
\(r_{\mathrm{norm}}\) is free. The minimal missing source field is the matched
full-pendulum \(m=2\) numerator divided by the same full inertia, equivalently
an independent value of \(r_{\mathrm{norm}}\) with covariance.

## 4. ToS anchors and the open affine correction

At \(r_{\mathrm{norm}}=c_f=0\), the seven magnetic-only anchors are:

| run | \(G_{\mathrm{hom}}^{(0,0)}\) (SI) | partial zero-covariance RSS |
|---|---:|---:|
| TOS-I-F1 first | \(6.673572163641\times10^{-11}\) | \(9.0668309\times10^{-16}\) |
| TOS-I-F1 repeat | \(6.673641921634\times10^{-11}\) | \(9.3727511\times10^{-16}\) |
| TOS-I-F2 | \(6.673682354302\times10^{-11}\) | \(2.1337984\times10^{-15}\) |
| TOS-I-F3 first | \(6.673630409563\times10^{-11}\) | \(9.9057307\times10^{-16}\) |
| TOS-I-F3 repeat | \(6.673612686798\times10^{-11}\) | \(8.9577182\times10^{-16}\) |
| TOS-II-F4 first | \(6.673451627184\times10^{-11}\) | \(1.1333983\times10^{-15}\) |
| TOS-II-F4 repeat | \(6.673533949745\times10^{-11}\) | \(1.1238585\times10^{-15}\) |

For each row the executable also reports the exact algebraic affine form at
\(r_{\mathrm{norm}}=0\),

\[
 G(c_f)=G_{\mathrm{hom}}^{(0,0)}
       +\frac{\Delta\omega^2}{K_{\mathrm{hom}}}c_f.           \tag{C05}
\]

There is no independently identified deterministic compact interval for (C05)
from the pinned public packet: a quoted standard uncertainty is not a hard
bound, and the signed row-level \(c_f\) is absent from the pinned forward
fields. This does not deny that the paper publishes authors' processed-model
\(G\) summaries and standard uncertainties; those are quarantined comparators,
not bounds on this independently reconstructed parameter family.

The workbook period summaries provide a useful file-level numerator check for
the repeated fibre-1 row. A-B-A extraction and the independent
common-quadratic diagnostic give magnetic-only anchors
\(6.673623904485\times10^{-11}\) and
\(6.673641568497\times10^{-11}\), versus
\(6.673641921634\times10^{-11}\) from the printed Table-2 response. The first
uses overlapping three-day summaries and has no independent coverage claim;
the second is a drift diagnostic.

None is a source-model-free raw numerator. The Supplement states that the
gravitational nonlinearity from the source masses was corrected synchronously
when determining \(\Delta\omega^2\). Event-level periods and the signed
nonlinearity or correction ledger needed to undo that processing are not
public.

### Theorem CHGC-2 -- ToS non-identification

Public data determine the seven \(r_{\mathrm{norm}}=c_f=0\) anchors and the
slopes in (C05), but not a ToS conditional \(G\) point unless both
\(r_{\mathrm{norm}}\) and \(c_f\) are chosen. The minimal missing fields are
an independently reconstructed matched full stiffness numerator and the
signed row-level \(\Delta K/(I\Delta\omega^2)\); a fully independent numerator
further requires uncorrected response or event data and its
source-nonlinearity ledger.

## 5. Available uncertainty propagation only

For a conditional quotient \(G=Nf/K\), this lane propagates

\[
 u_{G,N}=|G|\frac{u_N}{|N|},\qquad
 u_{G,K}=|G|\frac{u_K}{|K|},\qquad
 u_{G,f}=|G|\frac{u_f}{|f|},\qquad
 u_{G,\mathrm{RSS}}
 =\sqrt{u_{G,N}^2+u_{G,K}^2+u_{G,f}^2}.                     \tag{C06}
\]

The last expression is evaluated with three reported components: response,
nominal-kernel public-input sensitivity, and the Supplementary-Table-1
magnetic-correction standard uncertainty. It is explicitly conditional on
zero cross covariance. RESULT.json separately reports the linear sum and the
effect of the nominal lane's local axis-box diagnostic. These are first-order
local diagnostics, not coverage intervals. No uncertainty is invented for the
central \(r_{\mathrm{norm}}\), ToS \(c_f\), missing campaign transfer, AAF
shear coordinates, or unknown covariances. Their absence is why the partial
RSS values cannot be promoted to uncertainties on a physical \(G\)
measurement.

## 6. Mixed-normalization diagnostic and the identifiability obstruction

The public full inertia does not repair the missing numerator. NSKR already
calculates the deliberately incomplete coefficient

\[
 K_{\mathrm{mix}}
 =\frac{N_{\mathrm{homogeneous\ core}}/G}{I_{\mathrm{full}}}, \tag{C07}
\]

which deliberately combines the homogeneous-core numerator with the reported
full inertia while omitting the unreported mass from the relevant gravitational
multipole. Substituting \(K_{\mathrm{mix}}\) changes the AAF quotients by
\(+1631.537953\) ppm and the ToS anchors by \(+152.258414\) ppm for ToS-I
or \(+153.469498\) ppm for ToS-II.

Both numbers use public arithmetic; neither is the full source coefficient.
The mixed value is not asserted to be a physically admissible second HUST
apparatus map. It is a numerical stress diagnostic showing what follows from
an unmatched numerator/denominator normalization. The load-bearing
identifiability result is instead structural: the public packet supplies a
scalar inertia functional but neither the missing density/CMM map nor a law
that determines the source-coupled \(m=2\) numerator from that scalar. Those
functionals are not interchangeable; supplying total inertia alone cannot
complete the matched gravitational numerator.

### Theorem CHGC-3 -- no independently identified public physical set

Because the unprocessed public fields provide neither the matched full
numerator nor a relation fixing it from the disclosed scalar inertia, they do
not identify a unique full-apparatus \(G\) in this reconstruction. The public
packet also specifies no bounded admissible domain or covariance law for the
missing maps and corrections, so it owns no deterministic compact interval for
this conditional family to propagate. This is a public-packet numerator
obstruction, not a claim of mathematical unboundedness under every conceivable
externally imposed apparatus domain, not a denial of the authors' published
processed-model summaries, not a lack of algebra, and not a request for an
accepted value of \(G\).

## 7. Post-calculation comparators and central remainder

After the primary packet is frozen, the executable opens the authors'
processed coefficients and derived rows.

- For AAF, processed-minus-homogeneous source gaps are
  \(-44.5292\), \(-52.8401\), and \(-41.1463\) ppm. The primary conditional
  quotients differ from the authors' rows by \(-44.7085\), \(-52.7841\), and
  \(-41.1332\) ppm. Using the displayed Table-1 correction on both sides, the
  executable verifies the algebraic identity
  \(G_{\mathrm{hom}}/G_{\mathrm{processed}}
  =K_{\mathrm{processed}}/K_{\mathrm{hom}}\) to below
  \(5\times10^{-16}\) dimensionlessly. The upstream forward retains a separate
  component-reconstructed mechanical factor only as a comparator; its extra
  digits are not source-owned by Table 1. Thus inserting the processed
  coefficient merely evaluates the processed-model forward; it does not
  independently reconstruct the missing mass map.
- For ToS, the processed-minus-homogeneous source gaps are \(-91.48\) to
  \(-103.34\) ppm. With the homogeneous kernel, matching the already-derived
  rows would require total bracket corrections of \(+83.11\) to
  \(+97.89\) ppm. With the processed source coefficient, the remaining
  inferred dynamic bracket is \(-5.46\) to \(-8.37\) ppm, reproducing the
  forward lane's unowned correction bracket.

These comparisons locate scalar central remainders; they do not identify the
underlying clamp, coating, CMM, or density map or independently measure
\(\Delta K\). They are never fed back into the primary conditional quotient.

## 8. Scientific disposition

This closes the no-lab arithmetic question sharply:

\[
\boxed{
\begin{gathered}
\text{independent nominal kernel + released response}
\longrightarrow \text{explicit conditional quotient or family},\\
\text{missing matched numerator and ToS correction}
\longrightarrow \text{no independently identified public physical }G
\text{ point or deterministic compact interval}.
\end{gathered}}                                               \tag{C08}
\]

The result materially advances track 3: the program now owns numerical
conditional homogeneous cross-checks, exact missing-parameter locations, and
a normalization stress diagnostic. It is not a new or independent measurement of \(G\).
It does not weaken or confirm
Gravity Formation Theory. The sharp no-lab ceiling is the absence of an
independent full mass-multipole or CMM model matched to full inertia, plus the
signed ToS correction and raw correction custody. Further fitting of the
processed \(G\) rows cannot supply those fields.
