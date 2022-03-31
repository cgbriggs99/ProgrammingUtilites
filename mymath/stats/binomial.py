#!/usr/bin/python3

from mymath import stats
import math
from mymath import integrate
import numpy as np
import mymath

__ibf_quad_limit = 45
# Incomplete beta function.
def ibf(x, a, b, points = 30) :
    if a == 1 :
        return 1 - (1 - x) ** b
    if b == 1 :
        return x ** a
    if x == 0 :
        return 0
    if x == 1 :
        return 1
    if points < __ibf_quad_limit :
        return integrate.gauss_legendre_integral(lambda t: t ** (a - 1) *
                                                 (1 - t) ** (b - 1),
                                               0, x,
                                               points) * math.exp(
                                                   math.lgamma(a + b) -
                                                   math.lgamma(a) -
                                                   math.lgamma(b))
    else :
        return simpsons_integral(lambda t: t ** (a - 1) (1 - t) ** (b - 1),
                                 0, x, points) *\
                                 math.exp(math.lgamma(a + b) - math.lgamma(a) -
                                          math.lgamma(b))

class BinomialDistribution(stats.Distribution) :
    __max_size = 2 ** 16
    __binomials = {}

    def __clear_cache() :
        __binomials.clear()

    def __set_size(size) :
        __max_size = size

    def __is_full() :
        return len(__binomials) >= __max_size

    def choose(n, k) :
        if n == 0 :
            return 1
        elif n <= k or k == 0 :
            return 1
        elif n > k or k < 0 or n < 0 :
            return 0
        elif (n, k) in __binomials :
            return __binomials[(n, k)]
        elif (n, n - k) in __binomials :
            return __binomials[(n, n - k)]
        else :
            ret = choose(n - 1, k) + choose(n - 1, k - 1)
            if not __is_full() :
                __binomials[(n, k)] = ret
            return ret
    @staticmethod
    def logchoose(n, k) :
        return sum(math.log(i) for i in range(1, n + 1)) - \
               sum(math.log(i) for i in range(1, k + 1)) - \
               sum(math.log(i) for i in range(1, n - k + 1))
    
    def __init__(self, p, n) :
        self.__p = p
        self.__n = n
        self.mean = n * p
        self.median = round(n * p)
        self.mode = round((n + 1) * p)
        self.variance = n * p * (1 - p)
        self.skew = (1 - 2 * p) / math.sqrt(self.variance)
        self.exkurtosis = (1 - 6 * p * (1 - p)) / self.variance
        self.entropy = sum(choose(n, k) * p ** k * (1 - p) ** (n - k) *
                           (logchoose(n, k) + k * math.log(p) +
                            (n - k) * math.log(1 - p)) for k in range(n + 1))
    def pdf(self, k) :
        return choose(self.__n, k) * self.__p ** k * \
               (1 - self.__p) ** (self.__n - k)
    def cdf(self, k) :
        if self.__n == k :
            return 1
        return ibf(1 - self.__p, self.__n - k, 1 + k)
    def moment(self, order) :
        return sum(choose(self.__n, k) * self.__p ** k *
            (1 - self.__p) ** (self.__n - k) * k ** order
            for k in range(self.__n + 1))
    def cmoment(self, order) :
        return sum(choose(self.__n, k) * self.__p ** k *
                   (1 - self.__p) ** (self.__n - k) * (k - self.mean) ** order
                   for k in range(self.__n + 1))
    def quantile(self, k, cycles = 100, conv = 1e-5) :
        try :
            return mymath.solve(lambda x: self.cdf(x) - k, 0, self.mode,
                                cycles, conv)
        except Exception :
            return mymath.solve(lambda x: self.cdf(x) - k, n, self.mode,
                                cycles, conv)
    
