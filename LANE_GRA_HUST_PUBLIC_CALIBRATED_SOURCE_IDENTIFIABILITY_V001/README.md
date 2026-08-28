# HUST public calibrated source-identifiability lane

This isolated no-lab lane advances the physical \(G\) cross-check using public HUST primary materials. It adds the published coating, clamp, ferrule, attachment, signed ToS anelastic, air and acquisition-transfer corrections to the sealed independently integrated homogeneous kernels.

Run:

```bash
python3 analyze_hust_public_calibrated_source.py --write
python3 verify_hust_public_calibrated_source.py
```

The central result is not a new adopted value of \(G\). It is a much narrower, physically calibrated family \(G_i(r_i)\), together with a proof that one independently owned normalized harmonic remainder per released row is sufficient for row-wise point evaluation. The authors' processed coefficients make comparator remainders numerically inferable; they do not provide an independent physical reconstruction, and the ten row coordinates are not claimed to be ten independent physical degrees of freedom.

See `THEOREM.md`, `RESULT.json`, `CALIBRATION_FIELDS.json` and `SOURCE_CUSTODY.json` for equations, computed values, provenance and ceilings.
