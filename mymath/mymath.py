#!/usr/bin/python3

def solve(f, a, b, cycles = 100, conv = 1e-5) :
    """
Finds a zero for a function. The zero must be bounded by the a and b arguments.
f: A callable object which takes a single value and outputs a single value.
a, b: The bounds for the root. The root must be between these two values, and
    f(a) * f(b) > 0, which tries to isolate the single root.
cycles: The maximum number of cycles to run, in order to prevent hangups due to
    bad convergence.
conv: The convergence criterion. Checked with the absolute difference in the
    current bounds to determine wether the algorithm has converged.
    """
    fa = f(a)
    fb = f(b)
    if fa * fb > 0 :
        raise Exception(f"Root not bounded! a = {a}, b = {b}, fa = {fa}, fb = {fb}")
    if abs(fa) < abs(fb) :
        a, b, fa, fb = b, a, fb, fa
    c = a
    fc = fa
    flag = True
    count = 0
    while fb != 0 and abs(b - a) > conv and count < cycles :
        if fa != fb and fb != fc and fa != fc :
            s = (a * fb * fc) / ((fa - fb) * (fa - fc)) + (fa * b * fc) / ((fb - fa) * (fb - fc)) + (fa * fb * c) / ((fc - fa) * (fc - fb))
        else :
            s = b + fb * (a - b) / (fa - fb)
        if ((3 * a + b) / 4 > b and (s < b or s > (3 * a + b) / 4)) or \
            ((3 * a + b) / 4 < b and (s > b or s < (3 * a + b) / 4)) or \
            (flag and abs(s - b) >= abs(b - c) / 2) or \
            (not flag and abs(s - b) >= abs(c - d) / 2) or \
            (flag and abs(b - c) < conv) or \
            (not flag and abs(c - d) < conv) :
            s = (a + b) / 2
            flag = True
        else :
            flag = False
        fs = f(s)
        d = c
        c = b
        fc = fb
        if fa * fs < 0 :
            b = s
            fb = fs
        else :
            a = s
            fa = fs
        if abs(fa) < abs(fb) :
            a, b, fa, fb = b, a, fb, fa
        count += 1
    return b
