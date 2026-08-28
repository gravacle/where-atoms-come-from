# Self-audit

**Lane:** `LANE_GRA_HUST_NOMINAL_SOURCE_KERNEL_RECONSTRUCTION_V001`

**Disposition:** `PASS_AFTER_INDEPENDENT_HOSTILE_REPAIR`

1. **Primary sources only?** Yes. The HUST/Nature supplement, Extended Data
   images, and main Nature uncertainty table are official publisher objects.
2. **PDF visually inspected?** Yes. Formula pages 4--6 and table pages 19--21
   were rendered with Poppler and inspected at original rendered resolution.
3. **Published processed coefficients used as inputs?** No. They live in a
   separate JSON object which the executable opens only after both kernel
   families and their sensitivities are complete.
4. **Any measured or accepted value of \(G\) used?** No. No value of \(G\)
   enters any source calculation, scan, fit, normalization, or comparison.
5. **AAF constructed before ToS?** Yes. The result and theorem preserve the
   requested order.
6. **Uniform spheres called exact without a domain?** No. Exactness is
   conditional on homogeneity, sphericity, and disjoint support. Positive
   surface clearances are calculated.
7. **Measured nonsphericity ignored silently?** No. It is transcribed and
   explicitly assigned to the missing full-apparatus remainder.
8. **The AAF four-distance mismatch projected away?** No. All four reported
   horizontal/vertical coordinate separations are retained in the nominal
   pairwise-centred realization.
9. **Is that realization uniquely fixed by the public fields?** No. An
   independent collision preserves all four pair separations and the overall
   centroid while changing the coefficient. The two remaining shear degrees
   of freedom, unreported off-axis coordinates, and residual CMM coordinates
   are explicit conditional/identifiability ceilings.
10. **Full pendulum inertia mixed with a core-only numerator?** No. The code
    reports that mix only as a labelled forbidden diagnostic and proves why it
    is non-identifying.
11. **Numerical convergence checked?** Yes. Independent cubature/azimuth
    resolutions agree at 10^-10 kg m^-3 or better. This is not called a
    certified quadrature error bound: exactness belongs to the conditional
    integral, analytic derivative, and shell-theorem reduction, while the
    displayed coefficients are converged numerical evaluations.
12. **Analytic ToS derivative used?** Yes. Equation (H04) avoids a finite-angle
    differentiation artifact.
13. **Campaign temperature transport complete?** No. AAF-I transports only
    the upper horizontal separation using the public coefficient; its other
    fields remain at the 23.7 C reference. TOS-I combines dimensions at the
    20.2 C reference with run separations tabulated at 20.1/20.3 C. These are
    explicitly partial public transports.
14. **Public sensitivities reconstructed?** Yes. AAF dimensions, masses,
    horizontal and vertical distance classes reproduce the official values;
    ToS dimensions and distance classes also reproduce them.
15. **Standard uncertainties called exact bounds?** No. RSS and L1 summaries
    have explicit coverage ceilings.
16. **Central remainder set to zero or physically localized?** No. The
    processed-minus-nominal scalar is quantified only after the comparator is
    opened; it is not promoted to a unique spatial mass-map remainder.
17. **Published processed coefficient promoted to an independent result?** No.
    It is consistently labelled a post-calculation comparator.
18. **Comparator quarantine proven historically?** No. Executable ordering,
    literal scans, and mutation tests establish code-level input quarantine;
    no claim about unknowable cognitive history is made.
19. **New \(G\) claimed?** No.
20. **RGRL/GFT confirmation claimed?** No.
21. **Complete conserved stress claimed?** No. Drive and support stress
    ownership remains open.
22. **GC16 called complete?** No. The missing transfer, event, covariance,
    stress, and remainder fields are itemized.
23. **Canonical/shared files edited?** No. Every authored or pinned object is
    confined to this isolated lane.
24. **Next step a rescue Hamiltonian or model expansion?** No. It is a bounded
    request for the missing apparatus mass maps and calibration/event packet
    required to execute the already frozen GC16 architecture.

The independent hostile audit reconstructed both kernel families without
importing the production module and passed 94/94 checks after the placement,
numerical-exactness, temperature-transport, remainder, and custody language was
repaired. The lane remains unpromoted to any canonical register.
