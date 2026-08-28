# Prospective protocol for a calibrated finite-apparatus (G) cross-check

**Protocol ID:** `GRA-GC-CFAGC-PROT-V001`

**Date:** 2026-08-27

**Governing theorem:** `GRA-GC-CFAGC-V001`

## 1. Purpose and claim limit

This protocol turns independently calibrated source mass, finite geometry,
trajectory, mechanical response, and readout into the forward map (GC16) and
the identified set (GC26).  It can cross-check the realized Einstein-endpoint
coefficient (G).  It does not infer a lineage charge, convert scalar gamma into
stress-energy, or test RGRL unless a separately derived and independently
calibrated lineage-dependent complete stress column is later supplied.

## 2. Freeze before examining the gravity channel

The following objects must be versioned and frozen prospectively:

1. A positive absolute scan domain ([p_{\min},p_{\max}]) for (p=Gs), chosen
   from apparatus stability and sensitivity rather than an accepted or hidden
   value of (G).
2. Source and detector Lagrangian mass measures, including a global source
   scale set (s\in[s_-,s_+]), obtained by nongravitational mass metrology.
3. Geometry and source-trajectory maps with their metrology covariance.
4. The torsion and auxiliary-mode calibration
   ((I,\kappa,b,I_x,\kappa_x,b_x,\lambda)), measured with an independent
   nongravitational actuator and source modulation disabled.
5. Readout gain, delay, filters, sample timing, and their calibration
   covariance.
6. Scored frequencies/configurations, null configurations, held-out rows,
   nuisance templates, remainder bounds, covariance estimator, and acceptance
   thresholds.
7. Whether calibration uncertainty is propagated at a fixed predeclared
   reference or retained as an explicit nuisance in the implicit forward map.

No scan bound, covariance reference, nuisance basis, remainder class, or row
selection may be retuned after looking at the inferred gravity amplitude.

## 3. Physical source completeness

Inventory source bodies, detector bodies, motors, supports, suspension,
controller reactions, electromagnetic fields, thermal loads, residual gas,
and any other stress needed for a conserved complete apparatus source.  A
prescribed accelerated mass trajectory alone is not a conserved
(T^{\mu\nu}).  At Newtonian order, use the measured mass-density kernel in
(GC05)--(GC09).  Bound velocity, stress-active-mass, retardation, support,
and higher weak-field corrections as physical remainders before the response
inverse.  Explicitly verify the full finite-apparatus zero-mode/range condition
GI25; the balanced density quadrupole by itself does not prove it.

## 4. Acquisition

Acquire phase-referenced complex response over multiple source orientations,
trajectory amplitudes, distances, and frequencies.  Include:

- balanced source modulation with fixed source mass and zero center-of-mass
  trajectory;
- source-modulation-off and symmetry-null rows;
- independent environmental witness channels;
- predeclared held-out orientations/frequencies; and
- calibration data sufficient to resolve the auxiliary mechanical mode.

Raw observations and calibration records must be retained.  A null residual is
not post-solution readout noise merely because it correlates with the source;
classify it by where it physically enters the equations.

## 5. Forward construction with single ownership

For every scored row:

1. Evaluate the exact discrete or extended-source integrals for (a_n) and
   (k_{g,n}); do not fit either geometry coefficient to the gravity data.
2. Construct (d_{\theta,n}) and (d_{x,n}) from the independent mechanical
   calibration.
3. Place (-pk_{g,n}) and (-\lambda^2/d_{x,n}) only in the dressed operator.
4. Place the metric-mediated source torque (pa_n) only in the inhomogeneous
   source column.
5. Place physical force/torque remainders ((r_{\theta,n},r_{x,n})) before the
   retarded inverse, homogeneous data (d_{h,n}) in the solution data column,
   and readout bias/noise ((B\eta,\epsilon)) after the physical solution.
6. Apply the measured transfer (C_n) once.

The resulting prediction must be algebraically reducible both from the full
coupled system (GC12) and from its Schur complement (GC13)--(GC16).

## 6. Domain and null checks

Check (GC17) on the full admitted (p) and geometry domain.  For the two-mode
static determinant, checking both product endpoints suffices because it is
affine in (p), but the first principal minor and every additional retained mode
must also be nonsingular.  Confirm in code that:

- the full coupled solve equals the Schur-complement solve;
- matched zero source, remainder, and homogeneous data give zero physical
  response for every admitted (p);
- source-scale rescaling leaves (F(G,s)=F(Gq,s/q)); and
- the scored signal difference is not contained in the predeclared nuisance
  column space.

## 7. Covariance and inference

Stack real and imaginary lock-in components.  Combine observation and
independent calibration covariance exactly once as in (GC23).  If the
covariance varies materially with (p), evaluate its quadratic form and log
determinant at each (p); otherwise document the prospective fixed
linearization.  Profile only predeclared nuisance columns.

Report, in this order:

1. the identified or nominal profile set for (p=Gs);
2. the independent source-scale set ([s_-,s_+]);
3. the quotient set
   (G\in[p_-/s_+,p_+/s_-]);
4. changes under every declared physical-remainder bound;
5. null-row and held-out-row diagnostics without refitting (p); and
6. the shift caused by replacing the full dressed operator with each
   prospectively declared incomplete operator.

A nominal (\Delta\chi^2) interval is not an exact coverage theorem.  If a
bounded remainder or calibration set is authoritative, report (GC26) rather
than collapsing it to a point estimate.

## 8. Success, failure, and future lineage comparison

The ordinary-source cross-check succeeds only if one common (p) predicts all
scored rows, survives held-out configurations and declared remainder sets, and
maps through independent mass calibration to a nonempty (G) set.  It fails if
the result depends on post hoc row selection, an uncalibrated source scale, a
bare rather than dressed operator, duplicated systematics, or omitted support
physics beyond the frozen remainder bound.

A later RGRL test must add a separately derived, conserved, independently
calibrated lineage stress contrast and perform the matched KEEP/BREAK
intervention.  This protocol supplies the ordinary gravitational transfer and
(G) calibration side only; it cannot manufacture the missing lineage-to-
stress type join.
