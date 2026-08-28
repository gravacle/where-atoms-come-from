# Official response-field semantics used by CHGC

The pinned official Supplementary Information is
41586_2018_431_MOESM1_ESM.pdf, SHA-256
5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb.
It was already visually audited by both parent lanes. This record isolates the
load-bearing response meanings for the cross-check.

1. Supplement pages 3--5 give the ToS and AAF equations reproduced as (C02)
   and (C04) in the theorem.
2. Supplementary Table 1 (printed page 18) reports the magnetic-damper
   correction \(\Delta G/G\) and its one-standard-deviation uncertainty:
   AAF-I/II \(455.40(1.95)\) ppm, AAF-III \(25.74(0.08)\) ppm; ToS fibres
   1--4 respectively \(0.47(0.08)\), \(7.13(1.19)\), \(0.32(0.05)\), and
   \(0.27(0.08)\) ppm. CHGC uses these displayed central corrections and
   propagates their displayed uncertainty fields. Extra digits reconstructed
   from component fields upstream are not promoted to source-owned precision.
3. Supplement page 8 states that the source-mass gravitational nonlinearity is
   corrected synchronously when determining \(\Delta\omega^2\). The Table-2
   response is therefore a corrected response summary, not a source-model-free
   raw numerator.
4. Supplement page 9 describes the AAF half-second averaging and numerical
   differentiation corrections. The released one-second stream is a raw-like
   processed intermediate rather than the original 20-kHz encoder record.
5. Supplementary Table 2 (printed page 19) labels \(\Delta C_g/I\), mean
   \(\Delta\omega^2\), and derived \(G\), with one-standard-deviation
   uncertainties.
6. Supplementary Table 3 (printed page 20) labels the processed multipole
   coefficients, campaign-average \(\alpha_t\), and derived \(G\). Its caption
   says \(\alpha_t\) has been corrected for the air-density effect.
7. The forward lane non-circularly extracts a ToS contrast from released
   three-day period summaries and an AAF source harmonic from a released
   one-second figure stream without using an accepted \(G\) or processed
   source coefficient. Those are released-file-level numerator checks, not
   proof of statistical independence or that every upstream correction is
   source-model independent.

Consequently CHGC distinguishes three levels: released-file-level response
extraction, corrected campaign response summaries, and already-derived
response/source quotients. Only the first two enter conditional calculations;
processed source-coefficient and derived-\(G\) keys are selected for comparison
only after the primary packet is complete. The complete pinned parent JSON is
necessarily parsed first; quarantine is key-level nonuse, not byte-level
non-reading.
