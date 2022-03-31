#!/usr/bin/python3

import numpy as np
import math
import mymath.integrate as integrate
import mymath

class Distribution :
    
    def __init__(self) :
        self.mean = None
        self.median = None
        self.mode = None
        self.variance = None
        self.entropy = None
        self.skew = None
        self.exkurtosis = None
        pass
    def pdf(self, x) :
        pass
    def cdf(self, x) :
        pass
    def moment(self, order) :
        pass
    def cmoment(self, order) :
        pass
    def test(self, x) :
        pass
    def stat(self, conf) :
        pass
    def quantile(self, x) :
        pass

def cdf(pfd, points = 30) :
    return lambda x: integrate.gauss_laguerre(lambda t: pdf(t) * math.exp(-t),
                                              points, reverse = True) + \
                                integrate.gauss_legendre_integral(pdf, 0, x,
                                                                  points)

def mean_cont(pdf, points = 30) :
    return integrate.gauss_laguerre(lambda t: pdf(t) * math.exp(-t), points, 1,
                                    True) + integrate.gauss_legendre(lambda t:
                                            pdf(t) * math.exp(t), points, 1)
def moment_cont(pdf, order, points = 30) :
    return integrate.gauss_laguerre(lambda t: pdf(t) * math.exp(-t), 
                                    points, order, reverse = True) + \
                integrate.gauss_laguerre(lambda t: pdf(t) * math.exp(t),
                                         points, order)
def cmoment_cont(pdf, order, points = 30) :
    m = mean(pdf, points)
    return integrate.gauss_laguerre(lambda t: pdf(t + m) * math.exp(-t),
                                    points, order, True) + \
                        integrate.gauss_laguerre(lambda t: pdf(t + m) *
                                                 math.exp(t), points, order)
def quantile_cont(cdf_, a, b, cycles = 100, conv = 1e-5) :
    return lambda x: mymath.solve(lambda t: cdf_(t) - x, a, b, cycles, conv)

def median_cont(cdf_, a, b, cycles = 100, conv = 1e-5) :
    return mymath.solve(lambda x: cdf_(x) - 0.5, a, b, cycles, conv)

def mode_cont(der_pdf, a, b, cycles = 100, conv = 1e-5) :
    return mymath.solve(der_pdf, a, b, cycles, conv)

def entropy_cont(pdf, points = 30) :
    return -gauss_hermite_integral(lambda x: pdf(x) * math.log(pdf(x)),
                                  points)

def _gthresh(n, a = 0.95, points = 30) :
    t = tstat(a / (2 * n), n - 2, points)
    return (n - 1) * t / math.sqrt(n) / math.sqrt(n - 2 + t ** 2)

def grubbs_test(data, a = 0.95, points = 30) :
    g = _gthresh(len(data), a, points)
    m = mean(data)
    s = math.sqrt(svariance(data))
    inds = []
    highest = max(map(lambda x : abs(x - m), data)) / s
    for i in range(len(data)) :
        test = abs(data[i] - m) / s
        if test > g :
            inds += [(i, data[i], test)]
    return (highest > g), highest, abs(highest - m) / s, g, inds









