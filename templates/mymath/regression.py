#!/usr/bin/python3

import numpy as np

def poly(coefs, x) :
    s = 0
    p = 1
    for c in coefs :
        s += p * c
        p *= x * x
    return s

def polygrad(coefs, xs, ys) :
    return np.array([-sum(xs[j] ** (2 * i) * (poly(coefs, xs[j]) - ys[j]) for j in range(len(xs))) for i in range(len(coefs))])

def least_squares(xs_, ys_, order, step, conv, gconv) :
    xs = np.array(xs_)
    ys = np.array(ys_)
    w0 = np.zeros(order)
    w1 = np.array([1 for i in range(order)])
    g0 = polygrad(w0, xs, ys)
    g1 = polygrad(w1, xs, ys)

    while np.linalg.norm(w1 - w0) > conv and np.linalg.norm(g1 - g0) > gconv :
        w1 = w0
        w0 = w0 + step * g0
        g1 = g0
        g0 = polygrad(w0, xs, ys)
    return w0

def rsquared(data, fit) :
    m = mean(data)
    stot = sum((d - m) ** 2 for d in data)
    sres = sum((d - f) ** 2 for d, f in zip(data, fit))
    return 1 - sres / stot

def covariance(x, y, weights = None) :
    mx = mean(x, weights)
    my = mean(y, weights)
    if weights == None :
        return sum((x_ - mx) * (y_ - my) for x_, y_ in zip(x, y)) / len(x)
    else :
        return sum((x_ - mx) * (y_ - my) * w for x_, y_, w in
                   zip(x, y, weights)) / sum(weights)
    
def linfit(xs, ys) :
    m = covariance(xs, ys) / pvariance(xs)
    b = mean(ys) - m * mean(xs)
    return m, b
