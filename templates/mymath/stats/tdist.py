#!/usr/bin/python3

import math
from mymath import integrate
import mymath




# Above around 45 points, the root-finding algorithm starts to return
# complex values, though all roots should be real. Therefore, above this
# value, switch to Simpson's rule instead of Gauss-Legendre quadrature.
__gl_quad_limit = 40
__exp_gamma_limit = 37
def cdf(x, nu, points = 30) :
    """
The cumulative distribution function of Student's t distribution. Approaches
0 as x approaches -infinity, and 1 as x approaches infinity.
x: The point to evaluate.
nu: The number of degrees of freedom.
points: The resolution of the quadrature.
"""
    func = lambda x_ : (1 + x_ ** 2 / nu) ** ((-nu - 1) / 2)
    pf = 0
    if nu > __exp_gamma_limit :
        pf = math.exp(math.lgamma((nu + 1) / 2) - math.lgamma(nu / 2))
    else :
        pf = math.gamma((nu + 1) / 2) / math.gamma(nu / 2)
    if points < __gl_quad_limit :
        return pf / math.sqrt(math.pi * nu) * \
               integrate.gauss_legendre_integral(func, 0, x, points)    

    else :
        return integrate.simpsons_integral(func, 0, x, points) *\
               pf / math.sqrt(math.pi * nu)
# Initial guess comes from the fact that the CDF for 1 dof is less than any
# other CDF. The CDF for 1 dof is simply 1/2 + atan(x) / pi. The prefactor of
# 1.1 is there to make sure the upper bound will always be higher than the
# lower ensuring the root is always bounded.
def statistic(degs, acc = 0.95, cycles = 100, conv = 1e-5, points = 30) :
    """
Calculates the two-tailed t statistic for a given distribution at a given
confidence.
degs: The degrees of freedom.
acc: The amount of confidence.
cycles: The maximum number of cycles for the root finding algorithm.
conv: The absolute convergence.
points: The resolution of the quadrature.
"""
    if not math.isfinite(degs) :    #Limit case. Approaches normality.
        return zstat(acc, cycles, conv)
    if degs >= 1 :
        return mymath.solve(lambda x : cdf(x, degs, points) -
                 cdf(-x, degs, points) - acc, 0,
                 1.1 * math.tan(math.pi * acc - math.pi / 2), cycles, conv)
    elif degs > 0 : # Can't use the above method due to assumptions. Use Newton.
        otacc = 0.5 + acc / 2
        # 0 = mx + b = pdf(0) x - otacc
        # x = otacc / pdf(0)
        x0 = otacc * math.sqrt(math.pi * degs) * math.gamma(degs / 2) / \
            math.gamma((degs + 1) / 2)
        x1 = -1
        count = 0
        while abs(x1 - x0) > conv and count < cycles :
            x1 = x0
            x0 = x0 - t_cdf(x0, degs, points) / (math.gamma((degs + 1) / 2) /
                                                (math.gamma(degs / 2) *
                                                 math.sqrt(degs * math.pi) *
                                                 (1 + x0 ** 2 / degs) **
                                                 ((degs + 1) / 2)))
            count += 1
        if abs(x1 - x0) > 4 * conv :
            raise ArithmeticError(f"Did not converge! x0 = {x0} x1 = {x1} degs = {degs} acc = {acc} cycles = {cycles} conv = {conv} points = {points}")
        else :
            return x1


def confidence(data, acc = 0.95) :
    m = mean(data)
    d = math.sqrt(svariance(data))
    degs = len(data) - 1
    t = tstat(degs, acc)
    return m, t * d / math.sqrt(len(data))
