#!/usr/bin/env python3
"""
R_ladder — THE MONOTONE LADDER.  The one construction this lane adds, in one place.

GOAL.  Given an irrational target gamma in (0,1) and a rung sequence k_1 | k_2 | k_3 | ...,
produce a single real alpha with

        frac(k_j alpha) = gamma + eps_j,     0 < eps_j < 2 k_j / k_{j+1}     for every j.

STRICT POSITIVITY OF eps_j IS THE WHOLE POINT.  It is what makes the two-dimensional pair
(alpha, beta) provably NON-RESONANT in R_05 (see the proof there), and it is obtained for free
by building alpha as an INCREASING limit.

CONSTRUCTION.  alpha_j := (n_j + gamma)/k_j.  Given n_j and k_{j+1} = M_j k_j (M_j >= 2 an
integer), put
        W := M_j n_j + (M_j - 1) gamma,      n_{j+1} := floor(W) + 1.
Then n_{j+1} - W = 1 - frac(W) lies STRICTLY in (0,1), because gamma is irrational and M_j >= 2
makes (M_j - 1) gamma irrational, so W is irrational and frac(W) != 0.  Hence

        0 < alpha_{j+1} - alpha_j = (n_{j+1} - W)/k_{j+1} < 1/k_{j+1},

so (alpha_j) is STRICTLY INCREASING, alpha := lim alpha_j exists, and
        eps_j := k_j(alpha - alpha_j) = k_j sum_{i>=j} (alpha_{i+1}-alpha_i)
               in ( 0 , k_j sum_{i>=j} 1/k_{i+1} ) subset ( 0 , 2 k_j/k_{j+1} )
whenever M_i >= 2 for all i.  And frac(k_j alpha) = frac(n_j + gamma + eps_j) = gamma + eps_j
as long as gamma + eps_j < 1.

WHY THE DIPS ARE NECESSARILY SPARSE, STATED BEFORE THE NUMBERS.  A dip of depth D in the
running average at N = k_j needs log(1/eps_j) ~ D k_j, i.e. k_{j+1} ~ k_j e^{D k_j}.  So two
dips of depth ~1 cannot both sit inside seven decades unless the first is at very small k.
THIS IS A PROPERTY OF THE PHENOMENON, NOT OF THE COMPUTATION, and it is the reason four
decades of numerics can never by themselves settle a convergence question of this shape.

Arithmetic: gamma is supplied as an exact Fraction (a truncation of theta* to PREC digits,
from R_lib's pure-integer routines).  Every n_j and eps_j below is exact in that arithmetic;
the only error is the truncation of theta*, which is reported.
"""
import math
from fractions import Fraction


def log_int(x):
    """log of a Python int of arbitrary size."""
    if x <= 0:
        raise ValueError("log_int: non-positive")
    b = x.bit_length()
    if b <= 900:
        return math.log(x)
    sh = b - 900
    return math.log(x >> sh) + sh * math.log(2.0)


def log_fraction(fr):
    return log_int(fr.numerator) - log_int(fr.denominator)


def build_ladder(gamma, k1, n1, Ms):
    """Return (ks, ns, alphas, alpha_final, eps).
       Ms is the list of multipliers M_1, ..., M_{J-1}; k_{j+1} = M_j k_j.
       alpha_final := alpha_J (the last rung).  eps[j] = k_j(alpha_J - alpha_j) for j < J,
       and eps[J-1] = 0 by construction -- rung J is the ladder's truncation and is stated
       as such wherever it is used."""
    ks = [k1]
    ns = [n1]
    for M in Ms:
        W = M * ns[-1] + (M - 1) * gamma
        n_next = W.numerator // W.denominator + 1        # floor(W) + 1
        ns.append(n_next)
        ks.append(ks[-1] * M)
    alphas = [Fraction(n, 1) / k + gamma / k for n, k in zip(ns, ks)]
    alpha = alphas[-1]
    eps = [k * (alpha - a) for k, a in zip(ks, alphas)]
    return ks, ns, alphas, alpha, eps
