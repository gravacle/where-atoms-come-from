# FY native-support m=1 complete-H6 response lane

This no-lab successor to FX retains the native location of every physical
source insertion before the Feshbach reduction, proves exact recovery of the
frozen homogeneous source, and then evaluates one nonzero FO quotient
momentum (`m=1`, with `m=29` as its conjugate).

Run the cheap ownership and H2 gates:

```bash
python3 -B LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001/derive_native_support_m1_response.py
python3 -B LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001/verify_exact_m1_h2_lift.py
```

Run the exact complete-H6 enumeration and sampled finite response:

```bash
python3 -B LANE_GRA_FY_F3_Q4_NATIVE_SUPPORT_M1_COMPLETE_H6_RESPONSE_V001/derive_native_support_m1_response.py --full
```

The full run is deliberately expensive because it recomputes every diagonal
closed word through H6 with its tagged source location.  `m=1` is an exact
cyclic-translation label on the finite FO graph.  It is not, by itself, a
proof of continuum locality, a Ward identity, a massless tensor pole, RGRL-B,
gravity, or Newton's constant.
