# SPAG public-data second pass: clock and femtonewton-force components

**Artifact ID:** `GRA-SPAG-PUBLIC-DATA-SECOND-PASS-V001`  
**Search freeze:** `2026-08-28T00:15:24Z`  
**Status:** `SEALED_BOUNDED_SECOND_PASS__NO_NEW_LINEAGE_ESTIMAND__TWO_NEW_COMPONENT_DATASETS`  
**Parent:** `GRA-SPAG-PUBLIC-SUBSTITUTE-V001`  
**Shared MODEL/register files:** not edited

## 1. Result

This nonduplicative second pass found **two useful public component datasets**:

1. a blinded five-ensemble strontium-clock gravitational-redshift experiment
   with processed numerical data; and
2. a magnetically levitated sensor's femtonewton-scale Newtonian gravity-drive
   calibration with publisher source data.

Neither dataset contains an authenticated record-lineage intervention. Neither
contains randomized same-parent `M x L_T x L_D` support, and neither can estimate

\[
 \beta_{TM}={1\over 8}\sum_{t,d,m\in\{-1,+1\}}tm\,q_{t,d,m}.
 \tag{SP2-01}
\]

The search therefore produces a bounded negative result for the missing
lineage estimand, not a physical null. It also advances the no-laboratory lane
by locating two experimentally real component precedents for clock and force
readout. They may inform apparatus selection and data-release requirements;
they may not be pooled with other roots to create nonexistent factorial cells.

The preserved Panda response holdout was not opened, downloaded, or scored.

## 2. Frozen search contract

The exact repository universe, query strings, inclusion rules, source URLs,
checksums, and field profiles are frozen in `SEARCH_CUSTODY.json`.

The original search runner did **not** freeze the ordered returned-hit lists or
search-engine ranking snapshots for all 28 queries.  The query count and strings
are therefore reproducible, and the two retained roots are independently
re-checkable, but another runner cannot reproduce the exact screened result set
byte-for-byte.  Accordingly, "bounded negative result" here means only that no
qualifying root was retained by this documented pass.  It does not mean that the
listed repositories were exhaustively enumerated or that an unretained public
root cannot exist.

### 2.1 Primary inclusion gate

A dataset could directly score the frozen SPAG contrast only if it supplied:

- native randomized `M`, `L_T`, and `L_D` assignments;
- all eight support cells within the same physical parent;
- authenticated KEEP/BREAK or equivalent lineage custody;
- a gravity response joined to those assignments;
- source/probe geometry and SI calibration; and
- covariance, systematic, event, and analysis custody sufficient to reproduce
  the scored likelihood.

### 2.2 Secondary component gate

A nonduplicative dataset could be retained as a component candidate if it
supplied a real laboratory gravity-sensitive clock, force, torsion, or atom
response; public numerical data; a primary publisher or repository root; and
enough disclosed acquisition detail to materially inform a frozen protocol.

Theory-only records, simulations without an experimental response, paper-only
claims without public numerical data, uncontrolled geophysical maps, duplicate
roots, and records requiring invented lineage labels were excluded.

## 3. Candidate C1: miniature strontium clock network

### 3.1 Source and exact public content

Zheng et al. performed a blinded laboratory gravitational-redshift experiment
with five equally spaced spatial ensembles of strontium-87 spanning one
centimetre. The paper reports 14 blinded runs over three weeks, each lasting
one to four hours, and ten simultaneous pairwise clock comparisons. It reports

\[
 {\Delta f/f\over \Delta h}
 =[-12.4\pm0.7_{\rm stat}\pm2.5_{\rm sys}]\times10^{-19}/{\rm cm},
 \tag{SP2-02}
\]

consistent with the expected `-10.9e-19/cm`. The authors state that covariance
between pairwise comparisons sharing a clock was included in the error
estimation, and that the local gravitational acceleration was independently
measured as `-9.803 m/s^2`.

Primary sources:

- [article, DOI 10.1038/s41467-023-40629-8](https://www.nature.com/articles/s41467-023-40629-8)
- [Zenodo process data, DOI 10.5281/zenodo.8184043](https://zenodo.org/records/8184043)

The four deposited CSV files were inspected as spreadsheets. The deposit owns:

- `Fig1.csv`: five ensemble excitation-probability columns with
  `714, 719, 717, 716, 716` nonblank numerical values;
- `Fig2.csv`: 17 columns of processed systematic diagnostics, including atom
  number, detuning, Zeeman, black-body-radiation, trap-depth, height, and
  uncertainty fields;
- `Fig3.csv`: four height/redshift aggregates and 14 run-level gradient and
  uncertainty pairs; and
- `Fig4.csv`: all ten inferred pairwise height differences and uncertainties.

The Zenodo MD5 values and independently observed SHA-256 values are frozen in
`SEARCH_CUSTODY.json`. The files are not recopied into this lane; the article is
CC BY 4.0, while the Zenodo record's displayed rights field did not itself
state a data license at the freeze.

At independent hostile re-audit, all four files were downloaded again to a
temporary directory.  Their byte sizes, published MD5 values, observed SHA-256
values, table shapes, and numerical-column counts reproduced exactly.  No
downloaded source bytes were retained in this repository.

### 3.2 What it materially supplies

This is a public processed-data precedent for a compositionally distinct clock
probe, blinded acquisition, common-oscillator rejection, millimetre-scale
gravitational-potential resolution, and explicit shared-clock covariance in the
published analysis. It is useful for clock-endpoint sensitivity and release
planning.

### 3.3 Exact missing fields

The deposit does **not** supply:

- `M`, `L_T`, or `L_D`, a KEEP/BREAK operation, or any ancestry assignment;
- the eight same-parent support cells;
- per-run ten-pair frequency vectors or the per-run covariance matrices;
- a join from the probability rows to run, block, time, spin, blind-offset,
  phase, or detuning identifiers sufficient to regenerate the redshift fit;
- the experimental-control, analysis, or simulation code (available only on
  request); or
- a manipulable source-mass geometry for a matched lineage intervention.

The published statement that covariance was used is not a deposited covariance
object. Consequently this root cannot be transplanted into the frozen SPAG
likelihood or scored as a common-freefall lineage column.

## 4. Candidate C2: levitated femtonewton gravity-drive calibration

### 4.1 Source and exact public content

Yin et al. used a magnetically levitated force sensor to search for a symmetron
fifth force. Extended Data Figure 3 is a separate ordinary-gravity calibration
with a rotating tungsten disk. The publisher reports ten runs of 1000 seconds
with the disk at rest and rotating, integration over `15.35--15.36 Hz`, and
thermal/no-drive and gravity-drive baselines. The publisher source workbook
contains:

- 91 frequency/PSD rows from `15.310` to `15.400 Hz` for thermal noise;
- 91 frequency/PSD rows over the same grid for the gravity drive; and
- theoretical Newtonian interval `2.19--2.57 fN`, measured average `2.33 fN`,
  and measured standard deviation `0.33 fN`.

Primary sources:

- [article, DOI 10.1038/s41550-024-02465-8](https://www.nature.com/articles/s41550-024-02465-8)
- [publisher source workbook for Extended Data Figure 3](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41550-024-02465-8/MediaObjects/41550_2024_2465_MOESM6_ESM.xlsx)
- [supplementary information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41550-024-02465-8/MediaObjects/41550_2024_2465_MOESM1_ESM.pdf)

The publisher workbook was inspected as an XLSX and observed at 15,485 bytes,
SHA-256 `853a2f209d77c2b124d8319f16064944d33d6ba4a73cef7acb8228f9c847fb7b`.
It is not redistributed here because the article page states exclusive
publisher rights.

The independent hostile re-audit downloaded the workbook again to temporary
storage and reproduced both the hash and the three-sheet profile.  The 91-row
frequency grids and the four force-summary values were read directly from the
workbook, not inferred from the plotted figure.

### 4.2 What it materially supplies

This is an experimentally real femtonewton-scale ordinary-gravity drive and
noise precedent. It informs force-sensor frequency placement, integration-band
design, run duration, repetition count, and the scale at which a small-source
calibration has been demonstrated.

### 4.3 Exact missing fields

The public workbook does **not** supply:

- `M`, `L_T`, `L_D`, lineage custody, or eight-cell same-parent support;
- the ten individual force estimates underlying the reported mean and standard
  deviation;
- complete tungsten-disk geometry, density uncertainties, or geometric
  covariance in the workbook;
- the complete voltage-to-force SI transfer chain and its covariance;
- run-wise systematics, timestamps, or analysis code; or
- an independently generated Newtonian prediction that does not already use
  the standard gravitational law and its value of `G`.

The workbook therefore cannot supply an independent `G` cross-check. Its
`2.19--2.57 fN` band is a theory-produced comparison target, not a new
measurement of `G`. The paper's main source-film Newtonian amplitude
(`0.04 fN`, reported in the supplement) is distinct from the tungsten-disk
calibration and must not be conflated with it.

## 5. Bounded near-miss disposition

The second pass also checked official or primary routes for torsion balances,
atom interferometers, clocks, and mission data. The following did not create a
new admissible root:

- the AION prototype, Jaffe atom-source experiment, Fuchs, Panda, and the 2025
  semiclassical-gravity torsion experiment were already represented in the
  repository's prior work;
- official Eot-Wash pages surfaced publications but no new public event-level
  response packet with the required geometry/calibration/covariance custody;
- MICROSCOPE searches surfaced final calibrated publications but no newly
  identified official public raw mission packet satisfying this contract; and
- a public Kiskoros geophysical gravity/torsion survey was excluded because it
  is an uncontrolled field map, not a same-parent laboratory intervention.

No exact-lineage query returned a dataset with native randomized KEEP/BREAK or
`M x L_T x L_D` assignments. This statement is restricted to the frozen search
contract and date, with the returned-hit-log limitation stated in Section 2.
It is not a claim that no qualifying dataset existed, exists, or could later be
released.

## 6. Scientific disposition

The exact result is

\[
 \boxed{\text{no newly identified public root estimates }\beta_{TM};
 \text{ two new roots materially inform clock/force components.}}
 \tag{SP2-03}
\]

No lineage labels were manufactured. No incompatible physical roots were
pooled. No missing covariance was replaced by independence. No unavailable
data were called zero. No absence of a qualifying packet was called a null
result. The preserved Panda holdout remained closed: its four response-bearing
filenames remain absent from the holdout directory, while only the already
frozen metadata and paper/protocol custody objects are present.

The best next no-laboratory action is documentary rather than mechanical:
retain these two roots in the component-source inventory and use their missing
field lists when drafting a prospective release schema. A scored gravity
formation experiment still requires a new or newly released same-parent
lineage acquisition.
