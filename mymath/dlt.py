#!/usr/bin/python3

import math
import numpy as np

def dlt(data, tstep) :
    return lambda s : sum(data[i] * math.exp(s * i * tstep) for i in range(len(data)))


__coefs = {}
def laguerre_coefs(n) :
    if n == 0 :
        return [1]
    elif n == 1 :
        return [-1, 1]
    elif n in __coefs :
        return __coefs[n]
    else :
        l1 = laguerre_coefs(n - 1)
        l2 = laguerre_coefs(n - 2)
        ret = [((-l1[i] if i < len(l1) else 0) +
               ((2 * n - 1) * l1[i - 1] if i > 0 and i < len(l1) + 1 else 0) -
               ((n - 1) * l2[i - 2] if i > 1 and i < len(l2) + 2 else 0)) / n
               for i in range(len(l1) + 1)]
        __coefs[n] = ret
        return ret

def laguerre_poly(n, x) :
    coefs = laguerre_coefs(n)
    return sum(x ** (len(coefs) - i - 1) * coefs[i] for i in range(len(coefs)))

# If more points than this, use Simpson's Rule.
__quad_thresh = 32
def laplace(fun, points = 30) :
    if points <= __quad_thresh :
        nodes = np.roots(laguerre_coefs(points))
        weights = [x / ((points + 1) ** 2 * laguerre_poly(points + 1, x) ** 2)
                        for x in nodes]
        return lambda s : (sum(weights[i] * fun(nodes[i] / s)
                               for i in range(len(weights))) / s if s != 0 else
                           sum(weights[i] * fun(nodes[i]) * math.exp(nodes[i])
                               for i in range(len(weights))))
    else :
        raise Exception

#def ilaplace(fun, points = 30) :
#    if points <= __quad_thresh :
#        nodes = np.roots(
