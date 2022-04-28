#!/usr/bin/python3

import math
from mymath import integrate
import mymath
from . import stats


class NormalDistribution(stats.Distribution) :
    def __init__(self, m, sig) :
        self.mean = m
        self.median = m
        self.mode = m
        self.variance = sig ** 2
        self.__sig = sig
        self.skew = 0
        self.exkurtosis = 0
        self.entropy = 1 / 2 * math.log(2 * math.pi * math.exp(1) * sig ** 2)
    def pdf(self, x) :
        return math.exp(-(x - self.mean) ** 2 / (2 * self.__sig ** 2)) / \
               (math.sqrt(2 * math.pi) * self.__sig)
    def cdf(self, x) :
        return 1 / 2 * (1 + math.erf((x - self.mean) /
                                     (self.__sig * math.sqrt(2))))
    def moment(self, order, points = 30) :
        return math.sqrt(2) * self.__sig * integrate.gauss_hermite_integral(
            lambda x : (x * math.sqrt(2) * self.__sig + self.mean) ** order,
            points)
    def cmoment(self, order, points = 30) :
        return math.sqrt(2) * self.__sig * integrate.gauss_hermite_integral(
            lambda x: (x * math.sqrt(2) * self.__sig) ** order,
            points)
    def test(self, x) :
        return (x - self.mean) / self.__sig
    def stat(self, conf, cycles = 100, conv = 1e-5) :
        guess = math.tan(math.pi * acc - math.pi / 2)
        return mymath.solve(lambda x: normal_cdf(0, 1, x, -x) - acc, 0,
                        guess * 1.1, cycles, conv)
