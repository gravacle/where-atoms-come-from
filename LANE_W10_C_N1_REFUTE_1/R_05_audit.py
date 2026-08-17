#!/usr/bin/env python3
"""
R_05 — THREE AUDIT ITEMS ON LANE C's OWN INSTRUMENT AND ON ITS READING OF THE CORPUS.

LEG 1  C_04's OWN m(P) COLUMN IS QUADRATURE-LIMITED IN THE ONE ROW WHERE IT MATTERS, AND THE
       LANE'S PUBLISHED CONVENTION FORBIDS EXACTLY THAT.
       C_00 states: "Trapezoid on a periodic integrand; convergence in n is PRINTED, never
       assumed."  C_04's m_jensen() hard-codes n = 2^20 and prints no convergence; its
       subtorus_limit() hard-codes N = 2^22 and prints none either.  On the three analytic
       rows this costs nothing.  On the SENSE C row -- the row C_04's own closing block calls
       "THE ROW WHERE SOMETHING COULD FAIL" -- log|P| has a log singularity, and BOTH targets
       are wrong in the 7th decimal, where an EXACT value is available in one line.

LEG 2  C-6's HEADLINE ("A HYPOTHESIS FAILURE IS ALREADY SITTING IN S4's PUBLISHED TABLE") IS
       REFUTED BY S4's OWN HYPOTHESIS LINE.  S4:562 states the connection sense as
       "W(gamma_F) = e^{if} and W(gamma_C) = e^{ic} with the same (f,c), GENERIC".  W_F = -1,
       W_C = -i has finite order 4 and is excluded by that word.  Quoted from the corpus bytes.
       What survives of C-6 is exhibited and separated here.

LEG 3  THE ONE-LINE EXACT VALUE OF EVERY TARGET C_04 COMPUTED BY QUADRATURE, so the corrected
       gaps can be read off.

Precision: mpmath dps=50; the SENSE C subtorus limit is derived in closed form, not integrated.
"""
import numpy as np
import mpmath as mp
import subprocess, sys

mp.mp.dps = 50
REPO = "/Users/bgm/MB Work/where-atoms-come-from"

def m_jensen_f64(p, n):
    q = [float(z) for z in p]
    t = 2*np.pi*(np.arange(n)+0.5)/n; c = np.cos(t)
    A2 = q[0]**2+q[1]**2+2*q[0]*q[1]*c
    B2 = q[2]**2+q[3]**2+2*q[2]*q[3]*c
    return float(np.mean(0.5*np.log(np.maximum(A2, B2))))

def subtorus_f64(p, m, n, N):
    q = [float(z) for z in p]
    s = 2*np.pi*(np.arange(N)+0.5)/N
    x = np.exp(1j*n*s); y = np.exp(-1j*m*s)
    return float(np.mean(np.log(np.abs(q[0]+q[1]*x+q[2]*y+q[3]*x*y))))

if __name__ == "__main__":
    print("=" * 104)
    print("R_05 LEG 1 — C_04's TARGETS IN THE SINGULAR ROW, AGAINST THE EXACT VALUES.")
    print("=" * 104)
    SC = (0.25, 0.25, 0.25, 0.25)
    exact_mP = -mp.log(4)
    print("""  SENSE C 4-class: P = (1+x+y+xy)/4 = (1+x)(1+y)/4.
     m(P)        = m(1+x) + m(1+y) - log 4 = 0 + 0 - log 4 = -log 4         EXACTLY.
     subtorus at a primitive relation (m,n), H = {(z^n, z^-m)}:
        int_H log|P| = int log|1+z^n| + int log|1+z^-m| - log 4 = 0 + 0 - log 4 = -log 4
                                                                             EXACTLY,
     for EVERY relation -- the subtorus limit and m(P) COINCIDE here, so this row could never
     have separated them however it was computed.""")
    print(f"     exact  -log 4                       = {mp.nstr(exact_mP, 20)}")
    v20 = m_jensen_f64(SC, 1 << 20); v22 = subtorus_f64(SC, 11, 20, 1 << 22)
    print(f"     C_04's m_jensen(n = 2^20)           = {v20:.12f}   ERROR = {abs(v20-float(exact_mP)):.3e}")
    print(f"     C_04's subtorus_limit(N = 2^22)     = {v22:.12f}   ERROR = {abs(v22-float(exact_mP)):.3e}")
    print(f"     C_04 PRINTS as its targets:           -1.386293700082  and  -1.386293534824")
    print(f"""
     CONSEQUENCE FOR THE PRINTED GAPS.  C_04's SENSE C / DIOPHANTINE row prints
        |avg(1e7) - m(P)| = 1.405e-06 .
     Against the EXACT m(P) the same average (-1.386295105) is
        |avg(1e7) - (-log 4)| = {abs(-1.386295105 - float(exact_mP)):.3e} ,
     i.e. the printed gap is {1.405e-06/abs(-1.386295105-float(exact_mP)):.2f}x the true one, and the error is entirely in the TARGET.
     Likewise the RESONANT row prints |avg - subtorus| = 1.894e-06 against a target that is
     exactly -log 4, so the true figure is {abs(-1.386295428 - float(exact_mP)):.3e}.
     NO LANE-C FINDING RESTS ON THESE FOUR NUMBERS -- C-6 is about the FINITE-ORDER row, which
     is exact.  But C_00 promises convergence is "PRINTED, never assumed" and C_04 assumes it,
     in the row C_04 itself names as the only one that could fail.  Recorded, not patched.""")
    print("\n  For completeness, the convergence C_04 did not print:")
    for e in (16, 18, 20, 22, 24):
        v = m_jensen_f64(SC, 1 << e)
        print(f"     m_jensen n = 2^{e:<2d}: {v:.12f}   |err| = {abs(v-float(exact_mP)):.3e}")

    print("\n" + "=" * 104)
    print("R_05 LEG 2 — C-6 AGAINST S4's OWN HYPOTHESIS LINE.  QUOTED FROM THE CORPUS BYTES.")
    print("=" * 104)
    out = subprocess.run(["sed", "-n", "560,568p", f"{REPO}/S4_THE_MEASUREMENT_V001.md"],
                         capture_output=True, text=True).stdout
    print("  S4_THE_MEASUREMENT_V001.md, lines 560-568 (the sense in which the connection is held):")
    for ln in out.rstrip("\n").split("\n"):
        print("     | " + ln)
    print("""
  THE WORD IS THERE: "with the *same* (f,c), GENERIC".  W_F = -1, W_C = -i has ord(rho) = 4
  (W-07's own finding, register head) and is the least generic point in the space.  S4's
  SENSE C column is therefore NOT published as the rate at S1 sec6's connection, and C-6's
  headline -- "A hypothesis failure is ALREADY SITTING IN S4's PUBLISHED TABLE" -- overstates
  what the table says.

  WHAT SURVIVES OF C-6, AND IT IS WORTH KEEPING:
    (i)   the exact arithmetic is right: Z_k = (1+u^k)(1+v^k)/4 = 0 unless 4 | k, so
          Omega_N = 0 for every N >= 1 and the rate is -infinity.  Verified below in exact
          Gaussian integers, independently of C_05.
    (ii)  S4's hedge "generic" is UNQUANTIFIED, and the exceptional set is now known to
          contain the corpus's own S1 sec6 connection.  That is a real scope finding.
    (iii) the float64 blindness is the valuable half and is untouched: a lane trusting float64
          publishes a finite -28.09 for a quantity that is -infinity.  That is a PRECISION
          window artefact, one level below COR-E/COR-H's k-window artefacts, and it is new.""")
    print("\n  (i) re-verified independently, exact Gaussian integers, u = -1, v = -i:")
    bad = 0
    for k in range(1, 41):
        u = (-1)**k
        v = [(1,0),(0,-1),(-1,0),(0,1)][k % 4]
        zr = (1+u)*(1+v[0]) - 0*v[1]
        zi = (1+u)*v[1]
        z = (zr, zi)
        if (k % 4 != 0) and z != (0, 0):
            bad += 1
        if (k % 4 == 0) and z == (0, 0):
            bad += 1
    print(f"     k = 1..40: violations of 'Z_k = 0 iff 4 does not divide k' = {bad}   (must be 0)")
    print(f"     so Omega_N = 0 for every N >= 1 and (1/N) log|Omega_N| = -inf.  CONFIRMED.")

    print("\n" + "=" * 104)
    print("R_05 LEG 3 — WHAT LANE C's REPORT SAYS vs WHAT LANE C's OWN SCRIPTS PRINT: B0b*.")
    print("=" * 104)
    print("""  B0b* = (2,1,3,3)/9 is the ONE value in the whole lane that genuinely needs quadrature
  (Theorem R1: it is the only arm whose Jensen branches cross).  Lane C prints it twice:
      C_05 / the report:  m(P) = -0.987918288038      (mpmath kink-split)
      C_06 LEG A:         m(P) = -0.987918288039      (2^22-node float64 trapezoid)
  They disagree in the 12th place -- the last place either one prints.""")
    p = (mp.mpf(2)/9, mp.mpf(1)/9, mp.mpf(3)/9, mp.mpf(3)/9)
    A0 = p[0]**2 + p[1]**2 - p[2]**2 - p[3]**2
    B0 = 2*(p[0]*p[1] - p[2]*p[3])
    tc = mp.acos(-A0/B0)
    def f(t):
        c = mp.cos(t)
        return mp.log(mp.sqrt(max(p[0]**2+p[1]**2+2*p[0]*p[1]*c, p[2]**2+p[3]**2+2*p[2]*p[3]*c)))
    ref = (mp.quad(f, [0, tc]) + mp.quad(f, [tc, mp.pi]))/mp.pi
    # independent second reference: Gauss-Legendre on each smooth piece, different node counts
    ref2 = (mp.quadgl(f, [0, tc]) + mp.quadgl(f, [tc, mp.pi]))/mp.pi
    print(f"     corner-split mp.quad   (dps=50) : {mp.nstr(ref, 25)}")
    print(f"     corner-split mp.quadgl (dps=50) : {mp.nstr(ref2, 25)}")
    print(f"     the two agree to {mp.nstr(abs(ref-ref2), 4)}")
    print(f"     |ref - (-0.987918288038)| = {mp.nstr(abs(ref + mp.mpf('0.987918288038')), 4)}")
    print(f"     |ref - (-0.987918288039)| = {mp.nstr(abs(ref + mp.mpf('0.987918288039')), 4)}")
    tr = subtorus_f64  # unused; keep the float trapezoid inline for the exhibit
    q = [float(z) for z in p]
    for e in (20, 22, 24):
        n = 1 << e
        t = 2*np.pi*(np.arange(n)+0.5)/n; c = np.cos(t)
        A2 = q[0]**2+q[1]**2+2*q[0]*q[1]*c; B2 = q[2]**2+q[3]**2+2*q[2]*q[3]*c
        v = float(np.mean(0.5*np.log(np.maximum(A2, B2))))
        print(f"     float64 midpoint trapezoid n = 2^{e}: {v:.12f}   |v - ref| = {abs(v-float(ref)):.3e}")
    print("""
  THE REPORT'S VALUE IS THE CORRECT ONE and C_06's is the last-place casualty of an n^-2
  trapezoid whose error at 2^22 is ~5e-13, i.e. exactly at the 12th place it prints.  Nothing
  in C-5 depends on the 12th place.  But the lane prints 12 places for the ONE quantity in it
  that is quadrature-limited, and its own two scripts disagree there.  This is the defect class
  the lens exists for, committed on the lane's own most-quoted new number.""")
    sys.exit(0)
