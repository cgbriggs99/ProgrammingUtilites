#!/usr/bin/python3

import numpy as np
import math

def ellipse_dist(c, r1_, r2_ = [0, 0], pts = 50) :
    r1 = np.array(r1_)
    r2 = np.array(r2_)
    while r1[1] >= math.pi * 2 :
        r1[1] -= 2 * math.pi
    while r2[1] >= math.pi * 2 :
        r2[1] -= 2 * math.pi
    while r1[1] < 0 :
        r1[1] += 2 * math.pi
    while r2[1] < 0 :
        r2[1] += 2 * math.pi
    if r2[1] - r1[1] > math.pi :
        r1, r2 = r2, r1
    f = lambda t : c * math.hypot(math.sinh(r1[0] + (r2[0] - r1[0]) * t), math.sin(r1[1] + (r2[1] - r1[1]) * t)) * np.linalg.norm(r2 - r1)
    s = 0
    for t in np.linspace(0, 1, num = pts, endpoint = False) :
        s += (f(t) + 4 * f(t + 0.5 / pts) + f(t + 1 / pts)) / (6 * pts)
    return s

def cart_to_ellipse(c, rx, y = None, conv = 0.00000000001) :
    if y == None :
        r = np.array(rx)
    else :
        r = np.array([rx, y])
    r1 = np.array([math.acosh(np.linalg.norm(r)), math.atan2(r[1], r[0])])
    r0 = r1 + 1
    grad = np.array([1, 1])
    # Use gradient descent on (x - x0)^2 + (y - y0)^2 to find mu and nu.
    while np.linalg.norm(r1 - r0) > conv and np.linalg.norm(grad) > conv :
        r0 = r1
        gradf = np.array([c * math.sinh(r0[0]) * math.cos(r0[1]),
                          -c * math.cosh(r0[0]) * math.sin(r0[1])])
        gradg = np.array([c * math.cosh(r0[0]) * math.sin(r0[1]),
                          c * math.sinh(r0[0]) * math.cos(r0[1])])
        grad = (c * math.cosh(r0[0]) * math.cos(r0[1]) - r[0]) * gradf + \
               (c * math.sinh(r0[0]) * math.sin(r0[1]) - r[1]) * gradg
        r1 = r0 - 0.1 * grad
    return r1

def ellipse_to_cart(c, rmu, nu = None) :
    if nu == None :
        return (c * math.cosh(rmu[0]) * math.cos(rmu[1]),
                c * math.sinh(rmu[0]) * math.sin(rmu[1]))
    else :
        return (c * math.cosh(rmu) * math.cos(nu), c * math.sinh(rmu) \
                * math.sin(nu))

def cart_to_pol(rx, y = None) :
    if y == None :
        return (np.linalg.norm(rx), math.atan2(rx[1], rx[0]))
    else :
        return (math.hypot(rx, y), math.atan2(y, rx))

def pol_to_cart(rr, theta = None) :
    if theta == None :
        return (rr[0] * math.cos(rr[1]), rr[0] * math.sin(rr[1]))
    else :
        return (rr * math.cos(theta), rr * math.sin(theta))

def pol_to_ellipse(c, rr, theta = None) :
    x = pol_to_cart(rr, theta)
    return cart_to_ellipse(c, np.array(x))

def ellipse_to_pol(c, rm, n = None) :
    x = ellipse_to_cart(c, rm, n)
    return cart_to_ellipse(np.array(x))

def polar_dist(r1, r2) :
    return math.sqrt(r1[0] ** 2 + r2[0] ** 2 -
                     2 * r1[0] * r2[0] * math.cos(r1[1] - r2[1]))

def cart_dist(r1, r2) :
    return np.linalg.norm(np.array(r1) - np.array(r2))

ellip_area = lambda c, m, n : c ** 2 * math.hypot(math.sinh(m), math.sin(n))
ellip_line = lambda c, m, n : c ** 2 * math.hypot(math.sinh(m), math.sin(n)) * \
             np.array([1, 1])

polar_area = lambda r : r
polar_line = lambda r : np.array([1, r ** 2])
