#!/usr/bin/python3

import numpy as np
import math

# Dictionaries
__lcoefs_dict = {}
__jcoefs_dict = {}
__hcoefs_dict = {}
__lacoefs_dict = {}
__glacoefs_dict = {}
__tcoefs_dict = {}
__ucoefs_dict = {}

#Legendre polynomial coefficients
def _lcoefs(n, redo = False) :
    if n == 0 :
        return [1]
    elif n == 1 :
        return [1, 0]
    elif n in __lcoefs_dict and not redo :
        return __lcoefs_dict[n]
    else :
        l1 = _lcoefs(n - 1)
        l2 = _lcoefs(n - 2)
        ret = [(((2 * n - 1) * l1[i] if i < len(l1) else 0) - ((n - 1) *
                l2[i - 2] if i < len(l2) + 2 and i > 1 else 0)) / n for i in
               range(len(l1) + 1)]
        if not _is_cache_full() :
            __lcoefs_dict[n] = ret
        return ret

def _clear_cached_data() :
    """
Clear caches of polynomial coefficients.
"""
    __lcoefs_dict.clear()
    __jcoefs_dict.clear()
    __hcoefs_dict.clear()
    __lacoefs_dict.clear()
    __glacoefs_dict.clear()
    __tcoefs_dict.clear()
    __ucoefs_dict.clear()

__total_cache_entries = 2 ** 16
def _set_total_cache_entries(n = 2 ** 16) :
    """
Sets the total entries in the dicts. Any above this value will be dropped.
"""
    __total_cache_entries = n

def _is_cache_full() :
    return (len(__lcoefs_dict) + len(__jcoefs_dict) +
            len(__hcoefs_dict) + len(__lacoefs_dict) + len(__glacoefs_dict) +
            len(__tcoefs_dict) + len(__ucoefs_dict))\
            >= __total_cache_entries

def _lderiv(x, n) :
    coefs = _lcoefs(n)
    s = 0
    for i in range(len(coefs) - 1) :
        if i == len(coefs) - 2 :
            s += coefs[i]
        else :
            s += x ** (len(coefs) - i - 2) * (len(coefs) - i - 1) * coefs[i]
    return s

def gauss_legendre_integral(f, low = -1, high = 1, points = 30) :
    nodes = np.roots(_lcoefs(points))
    weights = map(lambda x: 2 / ((1 - x ** 2) * _lderiv(x, points) ** 2), nodes)
    return (high - low) / 2 * sum(w * f(x * (high - low) / 2 + (high + low) / 2)
                                  for x, w in zip(nodes, weights))

def _jcoefs(a, b, n) :
    if n == 0 :
        return [1]
    if n == 1 :
        return [(a + b + 2) / 2, a + 1 - (a + b + 2) / 2]
    if (a, b, n) in __jcoefs_dict :
        return __jcoefs_dict[(a, b, n)]
    else :
        j1 = _jcoefs(a, b, n - 1)
        j2 = _jcoefs(a, b, n - 2)
        ret = [((2 * n + a + b - 1) * (2 * n + a + b) *
               j1[i] / (2 * n * (n + a + b)) if i < len(j1) else 0) +
               ((2 * n + a + b - 1) * (a ** 2 - b ** 2) * j1[i - 1] / (2 * n *
                    (n + a + b) * (2 * n + a + b - 2))
                if i > 0 and i <= len(j1) else 0) -
               (2 * (n + a - 1) * (n + b - 1) * (2 * n + a + b) *
                j2[i - 2] / (2 * n * (n + a + b) * (2 * n + a + b - 2)) if
                i > 1 and i - 2 < len(j2) else 0) for i in range(len(j1) + 1)]
        if not _is_cache_full() :
            __jcoefs_dict[(a, b, n)] = ret
        return ret
def _jpoly(x, a, b, n) :
    coefs = _jcoefs(a, b, n)
    return sum(x ** (len(coefs) - i - 1) * coefs[i] for i in range(len(coefs)))

def _jderiv(x, a, b, n) :
    coefs = _jcoefs(a, b, n)
    s = 0
    for i in range(len(coefs) - 1) :
        if i == len(coefs) - 2 :
            s += coefs[i]
        else :
            s += x ** (len(coefs) - i - 2) * (len(coefs) - i - 1) * coefs[i]
    return s

def gauss_jacobi_integral(f, a, b, low = -1, high = 1, points = 10) :
    nodes = np.roots(_jcoefs(a, b, points))
    n = points
    weights = map(lambda x : -(2 * n + a + b + 2) / (n + a + b + 1) *
                  math.exp(math.lgamma(n + a + 1) + math.lgamma(n + b + 1) -
                           math.lgamma(n + a + b + 1) - math.lgamma(n + 2)) *
                  2 ** (a + b) / (_jderiv(x, a, b, n) * _jpoly(x, a, b, n + 1)),
                  nodes)
    return (high - low) / 2 * sum(w * f((high - low) / 2 * x + (high + low) / 2)
                                  for x, w in zip(nodes, weights))

def simpsons_integral(f, low, high, points) :
    return (high - low) / (6 * points) * sum(f(x) +
                                             4 * f(x +
                                                      0.5 * (high - low)
                                                      / points) +
                                             f(x + (high - low) / points) for
                                             x in np.linspace(low, high, points,
                                                              endpoint = False))


def _hcoefs(n) :
    if n == 0 :
        return [1]
    elif n == 1 :
        return [2, 0]
    elif n in __hcoefs_dict :
        return __hcoefs_dict[n]
    else :
        h1 = _hcoefs(n - 1)
        h2 = _hcoefs(n - 2)
        ret = [(2 * h1[i] if i < len(h1) else 0) - (2 * (n - 1) * h2[i - 1]
                                                    if (i - 1) >= 0 and
                                                    (i - 1) < len(h2) else 0)
               for i in range(max(len(h1), len(h2)) + 1)]
        if not _is_cache_full() :
            __hcoefs_dict[n] = ret
        return ret
def _hpoly(x, n) :
    coefs = _hcoefs(n)
    return sum(x ** (len(coefs) - i - 1) * coefs[i] for i in range(len(coefs)))

def gauss_hermite_integral(f, points) :
    nodes = np.roots(_hcoefs(points))
    n = points
    weights = map(lambda x: (2 ** (n - 1) * math.gamma(n - 1) *
                             math.sqrt(math.pi)) / (n ** 2 *
                                                    _hpoly(x, n - 1) ** 2),
                  nodes)
    return sum(w * f(x) for x, w in zip(nodes, weights))

def _lacoefs(n) :
    if n == 0 :
        return [1]
    elif n == 1 :
        return [-1, 1]
    elif n in __lacoefs_dict :
        return __lacoefs_dict[n]
    else :
        l1 = _lacoefs(n - 1)
        l2 = _lacoefs(n - 2)
        ret = [((2 * (n - 1) + 1) / n * l1[i - 1] if i - 1 >= 0 and
                i - 1 < len(l1) else 0) - (l1[i] / n if i < len(l1) else 0) -
               ((n - 1) * l2[i - 1] if i - 1 >= 0 and i - 1 < len(l2) else 0)
               for i in range(max(len(l1), len(l2)) + 1)]
        if not _is_cache_full() :
            __lacoefs_dict[n] = ret
        return ret

def _glacoefs(n, a) :
    if n == 0 :
        return [1]
    elif n == 1 :
        return [-1, 1 + a]
    elif n in __glacoefs_dict :
        return __glacoefs_dict[n]
    else :
        l1 = _glacoefs(n - 1)
        l2 = _lacoefs(n - 2)
        ret = [((2 * (n - 1) + 1 + a) / n * l1[i - 1] if i - 1 >= 0 and
                i - 1 < len(l1) else 0) - (l1[i] / n if i < len(l1) else 0) -
               ((n - 1 + a) * l2[i - 1] if i - 1 >= 0 and
                i - 1 < len(l2) else 0)
               for i in range(max(len(l1), len(l2)) + 1)]
        if not _is_cache_full() :
            __glacoefs_dict[n] = ret
        return ret

def _lapoly(x, n) :
    coefs = _lacoefs(n)
    return sum(x ** (len(coefs) - i - 1) * coefs[i] for i in range(len(coefs)))

def _glapoly(x, n, a) :
    coefs = _glacoefs(n, a)
    return sum(x ** (len(coefs) - i - 1) * coefs[i] for i in range(len(coefs)))

def gauss_laguerre_integral(f, points, alpha = 0, reverse = False) :
    """
Gauss-Laguerre quadrature with option for generalization.
If reverse is True, integrate from -infinity to low. If false, integrate from
low to infinity.
"""
    if alpha == 0 :
        nodes = np.roots(_lacoefs(points))
        n = points
        weights = map(lambda x: x / ((n + 1) ** 2 * _lapoly(x, n + 1) ** 2),
                      nodes)
    else :
        nodes = np.roots(_glapoly(points, alpha))
        n = points
        weights = map(lambda x: x * math.exp(math.lgamma(n + alpha + 1) -
                                             math.lgamma(n + 1)) /
                      ((n + 1) ** 2 * _glapoly(x, n, alpha) ** 2), nodes)
    if reverse :
        return (-1) ** (alpha + 1) * sum(w * f(-x)
                                         for x, w in zip(nodes, weights))
    else :
        return sum(w * f(x) for x, w in zip(nodes, weights))
    
def _tcoefs(n) :
    if n == 0 :
        return [1]
    elif n == 1 :
        return [1, 0]
    elif n in __tcoefs_dict :
        return __tcoefs_dict[1]
    else :
        t1 = _tcoefs(n - 1)
        t2 = _tcoefs(n - 2)
        ret = [(2 * t1[i] if i < len(t1) else 0) - (t2[i - 1] if i - 1 >= 0 and
                                                    i - 1 < len(t2) else 0)
               for i in range(max(len(t1), len(t2)) + 1)]
        if _is_cache_full() :
            __tcoefs_dict[n] = ret
        return ret

def _ucoefs(n) :
    if n == 0 :
        return [1]
    elif n == 1 :
        return [2, 0]
    elif n in __ucoefs_dict :
        return __ucoefs_dict[1]
    else :
        u1 = _ucoefs(n - 1)
        u2 = _ucoefs(n - 2)
        ret = [(2 * u1[i] if i < len(u1) else 0) - (u2[i - 1] if i - 1 >= 0 and
                                                    i - 1 < len(u2) else 0)
               for i in range(max(len(u1), len(u2)) + 1)]
        if _is_cache_full() :
            __ucoefs_dict[n] = ret
        return ret

def _tpoly(x, n) :
    coefs = _tcoefs(n)
    return sum(x ** (len(coefs) - i - 1) * coefs[i] for i in range(len(coefs)))

def _upoly(x, n) :
    coefs = _ucoefs(n)
    return sum(x ** (len(coefs) - i - 1) * coefs[i] for i in range(len(coefs)))

def gauss_chebyshev1_integral(f, low, high, points) :
    return (high - low) / 2 * sum(math.pi / points *
                                  f(x * (high - low) / 2 + (high + low) / 2)
                                  for x in map(lambda i: math.cos((2 * i - 1) /
                                                                  (2 * points) *
                                                                  math.pi),
                                               range(points)))
def gauss_chebyshev2_integral(f, low, high, points) :
    return (high - low) / 2 * sum(w * f(x * (high - low) / 2 + (high + low) / 2)
                                  for x, w in map(lambda i:
                                                  (math.cos(i / (points + 1) *
                                                            math.pi),
                                                   math.pi / (points + 1) *
                                                   math.sin(i / (n + 1) *
                                                            math.pi) ** 2),
                                                  range(points)))
    
