# PDF source-reading record

The pinned official Supplementary Information PDF was read with both text and
visual workflows, as required for formula and table custody.

- File: `SOURCE/41586_2018_431_MOESM1_ESM.pdf`
- SHA-256: `5b61d5c831be98c46e47fcc32f1ade0a680b4af6354d2bc34859d94b22279ffb`
- Text extraction: Poppler `pdftotext -layout`
- Visual rendering: Poppler `pdftoppm`, 170 dpi for PDF pages 4--6 and 150 dpi
  for PDF pages 19--21
- Visual inspection: completed at original rendered resolution

The visual reading fixed the following source meanings:

1. PDF page 4, printed Supplement page 3: ToS torque expansion,
   \(C_g\) definition, near/far frequency equation, and magnetic-damper
   correction.
2. PDF pages 5--6, printed Supplement pages 4--5: AAF multipole definition,
   \(P_{g,l,m}=8\pi mQ_{lm}q_{lm}/[I(2l+1)]\), and AAF response equation.
3. PDF page 19, printed Supplement page 18: full-pendulum inertias in
   Supplementary Table 1.
4. PDF page 20, printed Supplement page 19: processed ToS
   \(\Delta C_g/I\) values in Supplementary Table 2. These are quarantined as
   post-calculation comparators.
5. PDF page 21, printed Supplement page 20: processed AAF
   \(\sum_{l=2}^{10}P_{g,l,2}\) values in Supplementary Table 3. These are
   quarantined as post-calculation comparators.

The Extended Data geometry itself was visually read from the pinned official
300-dpi JPEGs for Tables 1, 2 and 4 and Figure 3. The main error-budget fields
were read from the pinned official Nature HTML Table 1. No OCR-derived value is
silently substituted for a visually verified field.

That visual check also establishes an important negative fact: Extended Data
Figure 3 and Table 4 publish labelled pair separations and the nominal layout,
not individual three-dimensional CMM coordinates for every source centre
relative to the rotation axis. Accordingly, the pairwise-centred AAF placement
used by this lane is recorded as an explicit conditional premise; it is not
represented as a uniquely published coordinate reconstruction.

The official Extended Data Table 4 caption was separately checked live at
`https://www.nature.com/articles/s41586-018-0431-5/tables/5` for the published
temperature statements. That page is source-located in `SOURCE_CUSTODY.json`
but is not separately hash-pinned in this lane.
