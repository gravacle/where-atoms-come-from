import sys, math; sys.path.insert(0, __file__.rsplit('/',1)[0])
from rx import *
K = K1(); cls = K.classes(); lam0 = lambdaB_exact(class_weights(K,{v:0.2 for v in K.V}))
print("ALL SIX ELEMENTARY EDGE COLLAPSES OF K1, source state SENSE U (p_v = 1/5)")
print("lambda(source) = %.12f\n" % lam0)
names = ['e1(v0~v1)','e2(v1~v2)','e3(v2~v0)','e4(v0~v3)','e5(v3~v4)','e6(v4~v0)']
for j,(s,t) in enumerate(K.E):
    mc = (cls[s][0]|cls[t][0], cls[s][1]|cls[t][1]); compat = cls[s]==cls[t]
    w = {k:0.0 for k in [(0,0),(1,0),(0,1),(1,1)]}
    for v in K.V: w[mc if v in (s,t) else cls[v]] += 0.2
    w2 = {k:0.0 for k in [(0,0),(1,0),(0,1),(1,1)]}; seen=False
    for v in K.V:
        if v in (s,t):
            if seen: continue
            seen=True; w2[mc]+=0.25
        else: w2[cls[v]] += 0.25
    print("%-11s %s+%s->%s  class-compatible=%-5s  PUSHFORWARD gap %.3e   SENSE-U-ON-TARGET gap %.3e"
          % (names[j], cls[s], cls[t], mc, compat, abs(lambdaB_exact(w)-lam0), abs(lambdaB_exact(w2)-lam0)))
print("\nspanning-tree collapse (all 5 vertices -> one, class (1,1)): |Z_k| = 1, lambda = 0 exactly,")
print("S = {(1,1)}, rank G = 0, formation False.")
