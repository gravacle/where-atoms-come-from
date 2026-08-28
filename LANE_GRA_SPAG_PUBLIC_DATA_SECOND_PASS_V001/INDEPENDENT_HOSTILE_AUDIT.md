# Independent hostile audit

**Artifact:** `GRA-SPAG-PUBLIC-DATA-SECOND-PASS-V001`  
**Audit freeze:** `2026-08-28T00:23:04Z`  
**Verdict:** `ACCEPT_WITH_EXPLICIT_SEARCH_REPRODUCIBILITY_CEILING`

## 1. What was independently checked

The audit re-downloaded all five public files to temporary storage from the
frozen Zenodo and Springer Nature URLs.  Nothing was copied into the repository.
Every claimed byte size and source digest reproduced:

| Object | Bytes | Reproduced digest |
|---|---:|---|
| `Fig1.csv` | 25,550 | SHA-256 `5ff66b2229f9c3b57d0fd2fa27e38aa773e9bfb597441c22393f9a65a36fed61`; MD5 `a829e2044f0ef2dd450435d7b790e8c7` |
| `Fig2.csv` | 2,876 | SHA-256 `2df463bc8840f1b366d6acd9e361bf58322b844b72f6309641fc66e60735bc30`; MD5 `8edb5d957e81a07b07e9d350afe43e3a` |
| `Fig3.csv` | 958 | SHA-256 `a4384040d5678893b4df15ea0be2980661023305f5009de876179fab1dea632f`; MD5 `257eab6f29cf3480c536cac73ed3998a` |
| `Fig4.csv` | 369 | SHA-256 `3c14df355b5cf1e6dcf138cf3b3de750f59ad270091e7874c70ad05204fa988d`; MD5 `53987b91b94d75d844f865ac0a778e75` |
| `Source Data Extended Data Fig. 3` XLSX | 15,485 | SHA-256 `853a2f209d77c2b124d8319f16064944d33d6ba4a73cef7acb8228f9c847fb7b` |

The CSV dimensions and numerical counts also reproduced.  In particular,
`Fig1.csv` has 721 rows, five columns, and numerical counts
`714, 719, 717, 716, 716`; the other CSV shapes are `12 x 17`, `16 x 6`, and
`3 x 20`.  The workbook contains exactly the three profiled sheets.  Its first
two sheets each contain 91 PSD rows from `15.310` through `15.400 Hz`, and its
summary row is exactly `2.57, 2.19, 2.33, 0.33 fN` in the publisher's stated
upper-theory, lower-theory, measured-mean, and measured-standard-deviation
order.

The primary article confirms the clock design, 14 blinded runs, five ensembles,
ten simultaneous comparisons, the reported gradient, the use of shared-clock
covariance, and the fact that only process data—not code—were deposited.  The
force article confirms ten 1000-second runs, the two PSD baselines, the
`15.35--15.36 Hz` integration band, and the distinction between a measured
average and a geometry-dependent Newtonian theory band.  Its supplement
separately confirms that the main thin-film source has a `0.04 fN` Newtonian
amplitude, so the packet correctly refuses to conflate that value with the
tungsten-disk calibration.

## 2. Estimand and physical-null challenge

Neither retained root contains native `M`, `L_T`, or `L_D`, a KEEP/BREAK
lineage intervention, or all eight factorial cells in a single parent.  The
clock root lacks deposited run-wise ten-pair vectors and covariance matrices;
the force root lacks individual force estimates and the complete SI transfer,
geometry, and covariance chain.  Thus neither root identifies

\[
  \beta_{TM}=\frac18\sum_{t,d,m}tm\,q_{t,d,m}.
\]

Combining the roots would manufacture cells across different apparatuses,
parents, units, and transfer functions.  The lane correctly forbids that move.
It also consistently describes the failure as a **data-design ceiling**, not a
zero lineage effect, a physical null, or evidence against gravity formation.

## 3. Query-accounting defect and repair

The 28 frozen query strings are present, nonblank, and unique.  The original
packet did not, however, preserve an ordered hit list or ranking snapshot for
each query.  Consequently an auditor can reproduce the queries and inspect the
two retained roots, but cannot reproduce the full screened result set
byte-for-byte.  This was a real custody defect in the phrase "bounded negative
result" if read as an exhaustive enumeration.

The repair makes that ceiling explicit in the report, result, custody object,
and verifier: the admissible claim is **no qualifying root was retained by this
documented pass**, not that every item in the named repositories was enumerated
or that no qualifying public dataset exists.  This limitation does not weaken
the component-dataset characterization or the conclusion that the two retained
roots cannot score `beta_TM`.

## 4. Panda holdout

The independent audit compared the local Zenodo metadata inventory with its
four published MD5 values and directly checked the holdout directory for the
four response-bearing filenames.  None is present.  The metadata file, protocol
freeze, and sign amendment are now separately dependency-pinned, and the
verifier performs the filename and inventory checks itself.  This establishes
filesystem custody at audit time; as always, it cannot prove a person's private
viewing history outside the repository.  No Panda response value was used or
reported by this audit.

## 5. Final disposition

After the two repairs, no material numerical, source-custody, estimand, pooling,
or physical-null defect remains.  The lane is accepted only at its stated
component-precedent and documented-search ceiling.  It supplies no SPAG score,
no empirical lineage evidence, no Gravity Formation Theory confirmation, and no
independent measurement or derivation of `G`.
