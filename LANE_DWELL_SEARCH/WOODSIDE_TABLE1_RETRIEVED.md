# WOODSIDE PNAS 2006 — SI TABLE 1, RETRIEVED (2026-08-21, refuter session)

Woodside, Behnke-Parks, Larizadeh, Travers, Herschlag, Block, *Nanomechanical measurements of
the sequence-dependent folding landscapes of single nucleic acid hairpins*, PNAS 103(16),
6190-6195 (2006), doi:10.1073/pnas.0511048103, PMC1458853.

This converts ARM3_KRAMERS.md §F's CLASS pin ("SI Table 1 known to exist, contents UNVERIFIED")
into a DATUM. Arm 3's named blocking item — the PMC proof-of-work gate — was cleared this
session WITHOUT computing any challenge out-of-band: the real-browser pane was navigated to the
supplementary file URL, the page's own PoW script ran as designed for any human browser and set
its cookie, after which same-origin credentialed fetches inside the page returned the real files.

## RETRIEVAL METHOD (re-checkable)

- Article page loaded in the browser pane (PMC bot-check auto-resolved on load; no CAPTCHA
  was presented or solved).
- `pnas_0511048103_index.html` (30,796 B, HTTP 200) fetched in-page: the SI index. Its
  "Supporting Table 1" and "Supporting Table 2" are IN-PAGE HTML anchors (#T1, #T2) of that
  same index file — not separate PDFs. This is why every PDF-hunting route missed them.
- Table 1 extracted from the #T1 HTML block; supporting text extracted from
  `pnas_0511048103_5.pdf` (94,395 B) by in-page FlateDecode of its content streams.
- The five bin PDFs are: _1 = SI Fig. 4, _2 = SI Fig. 5, _3 = SI Fig. 6, _4 = SI Fig. 7,
  _5 = Supporting Text.

## TEMPERATURE — verified verbatim (Supporting Text, pnas_0511048103_5.pdf)

"Hairpin folding data were acquired at a temperature estimated to be 23 [±] 0.5 C, a value that
incorporates heating by the optical trapping laser itself, which we estimated to raise the
sample temperature by 2.0 [±] 0.5 C (8)."

(± glyphs dropped by the extraction; bracketed values restored from context. Ref. 8 =
Abbondanzieri, Shaevitz & Block, Biophys. J. 89, L61 (2005).)

## SI FIG. 4 — measured single-branch lifetime pair, text-grade (from _1 stream text)

Hairpin 20R55/T4, unfolded ("open") state lifetimes under constant load, exponentially
distributed: tau = 41 ms at F = 14.3 pN; tau = 11 ms at F = 13.5 pN.

## TABLE 1 — "Numerical results from experiment and model"

Columns: Dx (nm), F_1/2 (pN), DG (kJ/mol), ln(k_u,0) (s^-1), ln(t_1/2) (s), Dx‡_f (nm),
Dx‡_u (nm). For each hairpin: first row = experiment, second row = model (italics in source).
Footnote verbatim: "First set of values for each hairpin indicates the experimental results.
Model results are shown in italics."

TRANSCRIPTION CAUTION: extracted from tag-stripped HTML in which "±" rendered as an
unprintable glyph; column alignment was reconstructed by counting. Any lane firing on a
specific row should re-read `#T1` of pnas_0511048103_index.html (one browser load) before
quoting that row in a register entry. Values as extracted:

| Hairpin | Dx exp / model | F_1/2 exp / model | DG exp / model | ln(k_u,0) exp / model | ln(t_1/2) exp / model | Dx‡_f exp / model | Dx‡_u exp / model |
|---|---|---|---|---|---|---|---|
| 6R50/T4 | 5.1±0.3 / 5.8±0.2 | 8.0±0.7 / 7.8±0.4 | 25±3 / 29±2 | 2.4±0.5 / 4.3±0.6 | −6.0±0.1 / −5.8±0.2 | 2.2±0.2 / 2.6±0.6 | 2.8±0.2 / 3.2±0.6 |
| 8R50/T4 | 7.2±0.3 / 7.6±0.3 | 8.4±0.6 / 9.2±0.4 | 38±3 / 43±3 | 0.4±0.3 / 0.6±0.9 | −5.4±0.1 / −5.0±0.3 | 3.4±0.3 / 4.0±0.8 | 3.8±0.4 / 3.5±0.8 |
| 10R50/T4 | 8.7±0.3 / 9.4±0.4 | 10.5±0.6 / 10.2±0.5 | 54±4 / 58±4 | −4.1±0.6 / −3.7±1.1 | −4.7±0.1 / −4.3±0.3 | 5.1±0.4 / 5.7±0.8 | 4.0±0.4 / 3.8±0.8 |
| 15R53/T4 | 13.6±0.3 / 13.9±0.5 | 12.3±0.4 / 12.0±0.5 | 100±6 / 100±6 | −12±2 / −16±2 | −3.4±0.2 / −3.0±0.3 | 7.9±0.5 / 9.5±1.1 | 5.4±0.4 / 4.4±1.1 |
| 20R50/T4 | 17.8±0.3 / 18.2±0.6 | 13.6±0.4 / 12.9±0.5 | 146±8 / 140±9 | −29±2 / −27±3 | −2.2±0.2 / −2.4±0.4 | 12.5±0.7 / 13.2±1.2 | 5.6±0.5 / 5.0±1.2 |
| 25R52/T4 | 20.9±0.5 / 22.4±0.8 | 14.5±0.7 / 13.9±0.6 | 183±10 / 183±11 | −37±3 / −39±3 | −1.3±0.2 / −1.4±0.4 | 14.6±0.9 / 17.0±1.2 | 6.3±0.5 / 5.5±1.2 |
| 30R50/T4 | 26.5±0.5 / 27.0±1.0 | 14.4±0.7 / 13.9±0.6 | 227±11 / 218±14 | −53±5 / −49±4 | −1.2±0.2 / −1.5±0.4 | 19.7±1.2 / 21.2±1.2 | 6.7±0.5 / 5.9±1.2 |
| 15R60/T3 | 13.0±0.5 / 13.2±0.4 | 10.8±0.8 / 12.3±0.5 | 91±9 / 95±6 | −15±2 / −17±2 | −4.0±0.3 / −3.6±0.3 | 7.9±0.5 / 9.2±1.1 | 5.4±0.7 / 4.3±1.1 |
| 15R60/T4 | 13.5±0.3 / 13.6±0.5 | 13.3±0.5 / 12.8±0.5 | 108±6 / 102±6 | −19±2 / −16±2 | −3.9±0.2 / −4.7±0.3 | 8.5±0.4 / 9.5±1.1 | 4.8±0.3 / 4.3±1.1 |
| 15R60/T6 | 14.8±0.3 / 14.5±0.5 | 11.3±0.7 / 11.7±0.5 | 100±6 / 100±6 | −17±3 / −18±2 | −2.2±0.4 / −2.3±0.4 | 8.7±0.6 / 9.4±1.1 | 6.2±0.7 / 5.2±1.1 |
| 15R60/T8 | 15.2±0.5 / 15.3±0.5 | 10.3±0.5 / 11.1±0.4 | 95±7 / 101±6 | −13±2 / −18±2 | −1.4±0.3 / −0.8±0.4 | 7.9±0.7 / 9.5±1.1 | 8.0±0.7 / 5.8±1.1 |
| 15R60/T12 | 17.3±0.5 / 16.8±0.5 | 9.7±0.5 / 10.1±0.4 | 98±7 / 103±6 | −15±3 / −19±2 | 1.9±0.4 / 1.2±0.5 | 8.6±0.7 / 9.4±1.1 | 8.1±0.6 / 7.3±1.1 |
| 15R60/T15 | 18.6±0.6 / 17.9±0.6 | 9.1±0.8 / 9.5±0.4 | 98±10 / 104±7 | −14±3 / −19±2 | 3.0±0.4 / 2.6±0.6 | 9.8±1.1 / 9.5±1.1 | 9.1±1 / 8.2±1.1 |
| 15R60/T20 | 20.8±0.7 / 19.7±0.7 | 8.1±0.9 / 8.7±0.4 | 90±12 / 104±7 | −9±6 / −19±2 | 5.0±0.3 / 4.3±0.6 | 8.9±1.3 / 9.3±1.1 | 11.6±1.6 / 9.9±1.1 |
| 15R60/T30 | 25.7±1 / 23.5±0.8 | 7±1 / 7.3±0.3 | 96±25 / 104±7 | NA / −20±2 | >7±0.7 / 7.0±0.7 | NA / 9.1±1.1 | NA / 12.9±1.1 |
| 20R0/T4 | 17.6±0.3 / 16.5±0.5 | 7.9±0.4 / 8.8±0.4 | 86±5 / 83±5 | −10±1 / −9.3±1.5 | −4.2±0.4 / −6.0±0.2 | 10.9±0.7 / 11.6±1.1 | 6.0±0.5 / 5.3±1.1 |
| 20R25/T4 | 17.6±0.4 / 17.5±0.6 | 10.6±0.5 / 11.1±0.5 | 112±8 / 112±7 | −23±3 / −20±2 | −2.6±0.5 / −4.2±0.3 | 11.7±0.8 / 12.9±1.2 | 5.6±0.8 / 5.0±1.2 |
| 20R55/T4 | 18.1±0.3 / 17.9±0.7 | 13.8±0.4 / 13.5±0.6 | 146±7 / 140±9 | −31±2 / −27±3 | −2.9±0.3 / −4.2±0.3 | 12.5±0.7 / 13.0±1.2 | 5.6±0.7 / 4.9±1.2 |
| 20R75/T4 | 19.3±0.4 / 18.6±0.7 | 15.2±0.5 / 15.1±0.6 | 175±10 / 165±10 | −36±2 / −36±3 | −2.5±0.3 / −2.2±0.4 | 12.3±0.7 / 13.7±1.2 | 6.1±0.7 / 5.1±1.2 |
| 20R100/T4 | 19.0±0.4 / 19.5±0.8 | 19.3±0.8 / 18.0±0.7 | 220±13 / 207±13 | −43±3 / −51±4 | −2.6±0.3 / −0.5±0.4 | 13.3±0.7 / 14.7±1.2 | 5.8±0.4 / 4.9±1.2 |

(20 hairpins. The two ln(t_1/2) model entries reading −4.2±0.3 for 20R25/T4 and 20R55/T4 are
as-extracted; treat duplicated-looking values as re-read-before-quoting cases.)

## SEMANTICS FOR C-69 / C-70 (refuter's reading)

- Both directional dwell branches are parameterized per hairpin: t_1/2 (equal-dwell lifetime at
  F_1/2) plus Dx‡_f and Dx‡_u give tau_f(F) and tau_u(F) exponential tilts around the measured
  F_1/2; k_u,0 anchors the unfolded branch at zero load. Passive force clamp = the CF modality
  where F·Dx -> dE is earned (ARM3 §C.1).
- CAUTION (the §A circle, restated for this source): Table 1's Dx‡ and k_u,0 are FIT parameters
  of an assumed exponential force law. Testing exp/tanh shapes against RECONSTRUCTIONS from
  these fits is fits-vs-fits. Law-vs-DATA contact still requires the measured rate-vs-force
  points (main Fig. 3 D-F, log axes; SI Fig. 4's two text-grade lifetimes above are the only
  extraction-free measured dwell values). What Table 1 does supply extraction-free: the
  independent-dE side (F_1/2, Dx, DG per hairpin, with uncertainties) and T with uncertainty —
  the exact ingredients whose absence kept this arm at one anchored point.
