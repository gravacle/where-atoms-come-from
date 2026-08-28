# SPAG Lane A: strongest no-laboratory public-data substitute

**Artifact ID:** `GRA-SPAG-PUBLIC-SUBSTITUTE-V001`  
**Date:** 2026-08-27  
**Status:** `SEALED_PUBLIC_COMPONENT_AND_PLANNING_RESULT__NO_PUBLIC_LANE_A_ESTIMAND`  
**Normative theory:** adopted complete-source-matched Lane A  
**Shared MODEL/URM:** not edited

## 1. Result

The strongest honest no-laboratory result is **not** a retrospective lineage
measurement. It is a joined public-observable inventory, an exact
identifiability obstruction, a raw-native detector/noise reconstruction, and an
optimistic apparatus-planning envelope.

No public packet inspected here contains independently randomized target and
dummy lineage redistributions (L_T,L_D), all eight same-parent
(M\times L_T\times L_D) cells, and a gravity response. Therefore none can
estimate

\[
 \beta_{TM}={1\over8}\sum_{t,d,m\in\{-1,+1\}}tm\,q_{t,d,m}.
 \tag{PD01}
\]

The public data establish ordinary gravity-response and component feasibility;
they do not establish, bound, or refute a Lane-A lineage response. Pooling
different experiments does not manufacture the missing randomized cells because
the packets have different physical parents, apparatuses, response units, and
transfer functions.

## 2. Usable public observables

| Packet | Public observables actually usable here | Exact use | Why it cannot score Lane A |
|---|---|---|---|
| Page--Geilker | 10 decay-derived decisions, deterministically assigned mass configurations, 10 torsion equilibrium shifts | branch-following endpoint; exact association and design-rank audit | decision and mass satisfy (X=M); no independent (L_T,L_D) |
| Fuchs/Zenodo 10300430 r4 | 8 geometries, 24 complex lock-in traces, native magnitude/phase/frequency/time streams | byte-custodied detector response, cadence, half-trace stability, and native residual spectral diagnostics | no lineage factors; demodulator identities, complete SI transfer, null generator, and covariance are absent |
| NIST/BIPM 2026 | paper-level source geometry, (31.1979\,\mathrm{nN\,m}) copper torque, two readout modes, Table-15 Type-A uncertainties | apparatus scale, configuration/systematics warning, and best-case eight-cell planning algebra | no event-level ancestry factors or eight-cell acquisition |
| Panda/Zenodo 10995225 | public four-file inventory and an atom-interferometric gravity endpoint | optional later compositionally different probe and preserved response holdout | no (L_T,L_D), no same physical parent with a torsion packet, no lawful software join |

The official open-access NIST/BIPM PDF is stored byte-for-byte at
`SOURCE/nist_bipm_2026.pdf`, SHA-256
`c79552d62f4d4f4e85cfbbb00f135c1d985b596d9cdcde9bee57cfe4618f33dc`.
The analysis verifies its 31 pages and extracts Table 15 directly from those
pinned bytes. Primary source: [Schlamminger et al., *Metrologia* 63, 025012
(2026)](https://doi.org/10.1088/1681-7575/ae570f); [official NIST
PDF](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=961075).

## 3. Exact public-data identifiability result

### 3.1 Page--Geilker

Even under the deliberately generous recoding (T=X=M) and fixed (D=+1),
the historical packet occupies only two of the eight SPAG cells. Its saturated
eight-column design has exact rank (2/8). In particular (TM=M^2=1), so the
desired interaction aliases an intercept-like term; other lineage and mass
columns also collapse. The actual experiment did not implement authenticated
SPAG lineage redistribution, making the physical ceiling stricter than this
generous algebraic proxy.

### 3.2 Fuchs, NIST/BIPM, and Panda

Each contains an ordinary source/geometry variation but no observed randomized
(L_T,L_D). Treating a missing factor as a constant makes its interaction
column duplicate an ordinary source or intercept column; inventing labels after
exposure has no causal meaning. Thus (\beta_{TM}) is not a function of any one
packet's observed likelihood.

### 3.3 Cross-root pooling

Let (a) index apparatus/root. A lawful combined model must admit root-specific
response maps and nuisance terms (h_a). Because no root crosses the lineage
factors, a transformation of an unobserved lineage coefficient can always be
absorbed into (h_a) without changing any observed prediction. Pooling values
in centimetres, volts, (mathrm{nN,m}), and (mathrm{nm,s^{-2}}) cannot break
that invariance without a separately validated common transfer and a common
intervention. No such object exists here.

Therefore:

\[
 \boxed{\text{no packet admitted here, or its cross-root pool, identifies }
 \beta_{TM}.}
 \tag{PD02}
\]

This is a proof of a data-design ceiling, not evidence that the physical effect
is zero.

## 4. Defensible quantitative work

### 4.1 Fuchs native detector/noise envelope

The pinned reconstruction contains 303,397,806 source bytes, eight geometries,
and 24 traces spanning 5.77--11.97 h. Its physical demodulator map remains
unknown. The residual control-band Welch-ASD ranges are therefore reported only
in native volts per square-root hertz:

| Neutral demodulator | minimum | median | maximum |
|---|---:|---:|---:|
| 1 | (4.024\times10^{-4}) | (3.175\times10^{-3}) | (1.320\times10^{-2}) |
| 2 | (2.027\times10^{-4}) | (5.286\times10^{-4}) | (3.090\times10^{-3}) |
| 3 | (6.893\times10^{-3}) | (9.377\times10^{-3}) | (1.795\times10^{-2}) |

The second-half/first-half harmonic ratios span 0.322--2.745 across channels.
These numbers can commission a raw-stream adapter. They are not a force-noise
model or achieved power because the deposit lacks the channel map, transfer
calibration, correlated null/sham process, independent-unit definition, and
full covariance.

### 4.2 NIST/BIPM best-case planning calculation

The official Table 15 gives Type-A, $k=1$, uncertainties from
$u_\Delta=0.0002$ to $0.0006\,\mathrm{nN\,m}$ for the full mass-position
torque difference $\Delta N$. If one independent difference were measured in
each of the four $(T,D)$ strata, then

\[
 \hat\beta_{TM}={1\over8}\sum_{t,d}t\,\widehat{\Delta N}_{td},
 \qquad
 u(\hat\beta_{TM})={u_\Delta\over4}.
 \tag{PD03}
\]

Using an explicitly conventional two-sided $\alpha=0.01$ and target power
$0.90$, a one-primary-contrast Gaussian planning boundary is

\[
 \delta_{90}
 =\bigl(z_{0.995}+z_{0.90}\bigr){u_\Delta\over4},
 \tag{PD04}
\]

which gives

\[
 \delta_{90}=0.000193\text{--}0.000579\ \mathrm{nN\,m}
 = (6.18\text{--}18.5)\times10^{-6}
 \times 31.1979\ \mathrm{nN\,m}.
 \tag{PD05}
\]

This is an **optimistic lower-bound planning envelope**, not a prospective
SPAG limit. It assumes the published Type-A performance transfers unchanged to
every new route stratum, independent equal-variance mass differences, no added
route or lineage systematics, and no full-family multiplicity penalty. The paper's three
copper orientations already span $0.0137\,\mathrm{nN\,m}$ (free) and
$0.0134\,\mathrm{nN\,m}$ (servo), while paired free-minus-servo offsets are
$0.0014$--$0.0020\,\mathrm{nN\,m}$. Those are not SPAG effects; they show why
geometry, configuration, covariance, and systematic ownership must precede a
claimed detection limit.

## 5. Frozen Lane-A decision rule

### 5.1 Retrospective public packets

The rule is categorical:

1. Require an authenticated, randomized, same-parent
   (M\times L_T\times L_D) eight-cell packet.
2. If either lineage factor, same-parent custody, or any required support cell
   is absent, return `PUBLIC_DATA_NO_LANE_A_SCORE`.
3. Do not repair the packet with post-hoc lineage labels, cross-root
   pseudo-cells, missing-as-zero coding, or ordinary source-motion labels.

Every public packet in this artifact returns `PUBLIC_DATA_NO_LANE_A_SCORE`.

### 5.2 Prospective acquisition

The physical prediction under the adopted complete-source-matched Lane-A
premise is (\beta_{TM}^{\rm phys}=0). Before the first scored response, freeze
the complete apparatus/collateral ledger, independent units and covariance,
all verdict-bearing coefficients, a familywise error rate no larger than 0.01,
power of at least 0.90 at the independently declared sensitivity boundary
(eta_q), and all equivalence bands.

With simultaneous confidence set (C_{TM}), form

\[
 {cal I}_{TM}
 =C_{TM}\oplus[-(\epsilon_{\rm coll}+\epsilon_{\rm geom}),
                 +(\epsilon_{\rm coll}+\epsilon_{\rm geom})].
 \tag{PD06}
\]

- Any failed lineage, route, source, calibration, covariance, sham, collateral,
  or custody gate gives `NO_ANCESTRY_RESULT`.
- Run A with (0\notin{cal I}_{TM}) and every control passed gives only
  `RUN_A_ANCESTRY_CORRELATED_RESIDUAL_CANDIDATE__NO_CONFIRMATION`.
- Run A with ({\cal I}_{TM}\subset[-\eta_q,+\eta_q]) gives only
  `RUN_A_EXPLORATORY_DECLARED_APPARATUS_BOUND`.
- Overlap gives `INCONCLUSIVE`; optional continuation is forbidden.
- A separately hashed held-out Run B that reproduces a nonzero result without
  an independently owned source bucket gives at most
  `REPRODUCED_ANCESTRY_CORRELATED_GRAVITY_RESIDUAL__SOURCE_BUCKET_UNRESOLVED`.
- A complete held-out null gives
  `BOUNDED_NULL__DECLARED_APPARATUS_COLUMN_ONLY`.
- An independently calibrated ordinary owner gives
  `COLLATERAL_OWNED__NO_LINEAGE_GRAVITY_RESULT`.

No off-shell RGRL-C rank claim is inferred from any on-shell residual.

## 6. Exact scientific ceiling

The public material can test or establish:

- branch-following and ordinary gravity-response endpoints;
- raw detector response, cadence, stability, and native noise diagnostics;
- paper-level force/torque scale and best-case factorial planning algebra; and
- separate torsion and atom-probe component feasibility.

It cannot test or establish:

- a causal record-lineage redistribution effect;
- a nonzero or bounded Lane-A (\beta_{TM});
- empirical RGRL-C or Gravity Formation Theory confirmation;
- a same-parent common-freefall ancestry response; or
- a lineage source functional or microscopic derivation of (G).

The no-laboratory lane is therefore complete at its honest ceiling. The next
empirical increment is not another retrospective relabelling exercise; it is a
new same-parent eight-cell acquisition with the frozen Lane-A controls.

## 7. Reproduction

Run:

```text
python3 LANE_GRA_SPAG_PUBLIC_DATA_SUBSTITUTE_V001/analyze_public_substitute.py
```

The zero-argument program verifies every pinned local dependency, extracts the
NIST/BIPM table from the exact PDF bytes, recomputes the Page--Geilker ranks,
summarizes the Fuchs native diagnostics, preserves the Panda response holdout,
and rewrites `RESULT.json` deterministically. It returns
`SPAG_PUBLIC_DATA_SUBSTITUTE: PASS` or refuses.
