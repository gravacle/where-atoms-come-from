# Repair2 adversarial matrix

| Case | Mutation | Required result |
|---|---|---|
| exact scalar list | observations = `[0]` | reproduce Repair1 escape; Repair2 refuses non-object |
| mixed list | three typed members plus scalar | refuse before inspecting any member fields |
| empty object | one lifecycle member = `{}` | refuse exact-key failure |
| unknown key | add one unregistered member key | refuse exact-key failure |
| missing identity | empty event identity | refuse identity failure |
| nonfinite value | encode `NaN` in one value | strict JSON/payload refusal |
| cross-event lifecycle | closure carries another event ID | refuse lifecycle identity join |

The typed three-stage lifecycle is also evaluated as an uncounted positive control so
the repair cannot pass merely by refusing every actual platform. All prior 47 cases
run under Repair2 unchanged, giving 54 counted cases total.
