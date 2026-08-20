#!/usr/bin/env python3
# LANE_T39_A_ENUM probe 1: ALPHA-THROUGH-E_B amplification numbers.
# tau = f0^-1 exp(E_b/kT)  =>  d ln tau / d ln alpha = (E_b/kT) * s,  s = d ln E_b / d ln alpha
# Mechanism scalings (m_e c^2 held fixed; per-atom energies; owners in the lane note):
#   s = 2 : chemical bond / diffusion activation / band offset / electrostatic (Hartree = alpha^2 m_e c^2)
#   s = 4 : shape (magnetostatic dipolar) anisotropy (mu_B^2/a0^3 = alpha^4 m_e c^2 / 4, exact)
#   s = 6 : uniaxial magnetocrystalline anisotropy, 2nd order in spin-orbit (xi ~ Z^2 alpha^4 m_e c^2,
#           K ~ xi^2/W with W ~ Hartree ~ alpha^2  =>  alpha^6)  [van Vleck 1937; Bruno 1989]
#   s = 10: cubic magnetocrystalline anisotropy, 4th order in SOC: xi^4/W^3 => alpha^10 (listed, unused:
#           census magnetite barrier is SHAPE anisotropy)
import math

census = [
    # (name, E_b/kT at operating T, s, mechanism note)
    ("HDD CoCrPt grain",      61,  6, "uniaxial MAE via SOC^2 (K_uV=2.5e-19 J, census #1)"),
    ("SD magnetite TRM",     780,  4, "shape/magnetostatic anisotropy (census #10)"),
    ("Zircon Pb closure",    220,  2, "diffusion activation 5.7 eV (Cherniak-Watson; census #3)"),
    ("Flash gate (Arrhenius detrap ~1.1 eV path)", 43, 2, "thermally activated detrapping"),
    ("Flash gate (full 3.1 eV barrier)", 120, 2, "SiO2 band offset (census #11)"),
    ("DNA depurination",      50,  2, "glycosidic 1.3 eV at 310 K (Lindahl; census #7)"),
    ("AgBr latent image",     45,  2, "Ag_n cluster chemistry, ln(tau*f0)=ln(3e7*1e12)"),
    ("CMOS latch (powered)",  50,  2, "electrostatic barrier 40-60 kT (census #5)"),
    ("Generic 60 kT chemical",60,  2, "the brief's example"),
    ("Generic 60 kT uniaxial-SOC", 60, 6, "same barrier, magnetic mechanism"),
]

print(f"{'mechanism':<46}{'E_b/kT':>7}{'s':>4}{'A=dlnTau/dlnA':>15}{'x for 1% dA/A':>16}{'x for 0.1%':>12}")
for name, b, s, note in census:
    A = b * s
    print(f"{name:<46}{b:>7}{s:>4}{A:>15}{math.exp(0.01*A):>16.3g}{math.exp(0.001*A):>12.4g}")

print()
# ratio signature for coexisting geological records (same Delta alpha/alpha epoch):
A_mag = 780*4; A_zr = 220*2
print(f"coexisting magnetite-vs-zircon log-lifetime-shift ratio: {A_mag}/{A_zr} = {A_mag/A_zr:.2f} (parameter: barriers only)")

# flash tunneling channel alpha-flatness: exponent kappa*d = sqrt(2 m Phi)/hbar * d
# Phi ~ alpha^2 m_e c^2 => kappa ~ alpha * (m_e c/hbar); d = N * a0, a0 ~ hbar/(alpha m_e c)
# => kappa*d ~ N * alpha^0  : alpha-independent at fixed oxide atom count N.
print("tunneling-limited retention exponent kappa*d ~ N (atom layers): d ln(kappa d)/d ln alpha = 0")

# composed with C-83: W = tau/t_epoch; if t_epoch is set by an attempt/phonon frequency (weak alpha dep.),
# d ln W/d ln alpha ~ (E_b/kT)*s, and the certifiability crossover n* ~ 6W shifts by the same factor.
print("C-83 composition: d ln n*/d ln alpha = d ln W/d ln alpha ~ (E_b/kT)*s (same amplification)")
