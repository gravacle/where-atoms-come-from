"""ADVERSARIAL CHECK: does C-69's stated closed-form law reproduce C-69's own predicted number?
Grounding says: tau = exp(E_b/kT)/f_0, E_b = 60.8 kT, f0 = 1e9 Hz  ->  claims tau = 5.40e16 s."""
import numpy as np
KB=1.380649e-23; T=300.0; f0=1e9
K_u=2.0e5; V=1.26e-24
E_b=K_u*V                     # J
x=E_b/(KB*T)
dE=3.0*KB*T
tau_stated_formula = np.exp(x)/f0                     # the formula IN THE GROUNDING
gd=f0*np.exp(-(E_b+dE/2)/(KB*T)); gu=f0*np.exp(-(E_b-dE/2)/(KB*T))
tau_actually_computed = 1.0/(gu+gd)                   # what the lane computes and reports
print(f"E_b/kT                      = {x:.4f}")
print(f"tau from the STATED formula exp(E_b/kT)/f0 = {tau_stated_formula:.4e} s")
print(f"tau the lane ACTUALLY computes 1/(gu+gd)   = {tau_actually_computed:.4e} s")
print(f"ratio = {tau_stated_formula/tau_actually_computed:.4f}   (= 2*cosh(dE/2kT) = {2*np.cosh(dE/(2*KB*T)):.4f})")
print(f"grounding's claimed number: 5.40e16 s -> matches the {'ACTUAL' if abs(tau_actually_computed-5.40e16)/5.40e16<0.01 else '???'} computation,")
print(f"NOT the stated formula, which gives {tau_stated_formula:.2e} s — off by {tau_stated_formula/5.40e16:.1f}x.")
# and on the t33 surfaces with large dE the gap between stated law and computed check grows:
eV=1.602176634e-19
for name,dE_,E_b_,T_,f0_ in [("NAND",0.30*eV,1.60*eV,300.,1e13),("DNA",0.35*eV,1.30*eV,310.,1e13),("azobenzene",0.60*eV,1.05*eV,300.,1e13)]:
    g1=f0_*np.exp(-(E_b_+dE_/2)/(KB*T_)); g2=f0_*np.exp(-(E_b_-dE_/2)/(KB*T_))
    print(f"{name:12s} stated-law tau = {np.exp(E_b_/(KB*T_))/f0_:.3e} s   computed tau = {1/(g1+g2):.3e} s   ratio {np.exp(E_b_/(KB*T_))/f0_*(g1+g2):.1f}x")
