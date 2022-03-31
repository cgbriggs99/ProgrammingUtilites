#!/usr/bin/python3

import numpy as np

def rk4(f, xstart, ystart, xend, points = 100) :
    step = (xend - xstart) / points
    y = ystart

    for i in range(points) :
        q1 = f(xstart + i * step, y)
        q2 = f(xstart + (i + 0.5) * step, y + step * q1 / 2)
        q3 = f(xstart + (i + 0.5) * step, y + step * q2 / 2)
        q4 = f(xstart + (i + 1) * step, y + step * q3)
        y += step / 6 * (q1 + 2 * q2 + 2 * q3 + q4)
    return y

def rk_1d(f, xstart, ystart, xend, coefs = [1 / 6, 1 / 3, 1 / 3, 1 / 6],
          offs = [0, 0.5, 0.5, 1], points = 100) :
    assert(len(coefs) == len(offs))
    step = (xend - xstart) / points
    y = ystart
    norm = sum(coefs)
    qs = [0 for i in range(len(offs))]
    for i in range(points) :
        qs[0] = f(xstart + (i + offs[0]) * step, y)
        for i in range(1, len(qs)) :
            qs[i] = f(xstart + (i + offs[i]) * step,
                      y + offs[i] * step * qs[i - 1])
        y += step / norm * sum(c * q for c, q in zip(coefs, qs))
    return y

class DiffFunc :
    def __init__(self, func, vnames) :
        self.__func = func
        self.__vnames = vnames
    def getvnames(self) :
        return self.__vnames
    def __call__(self, *args) :
        return self.__func(*args)

class DiffOp :
    def __init__(self, fname, vnames = []) :
        self.__func = fname
        self.__vnames = vnames
    def getfname(self) :
        return self.__fname
    def getvnames(self) :
        return self.__vnames

# Assumed to equal zero.
class DiffEq :
    def __init__(self, vnames) :
        self.__value = 0
        self.__vnames = vnames
    def getvnames(self) :
        return self.__vnames
    def __add__(self, other) :
        if type(other) is DiffFunc :
            out = DiffEq(self.getvnames())
            out.__value = ('+', self.__value, other.__value)
        elif type(other) is DiffFunc or type(other) is DiffOp:
            out = DiffEq(self.getvnames())
            out.__value = ('+', self.__value, other)
        else :
            if type(self.__value) is not tuple :
                out = DiffEq(self.getvnames())
                out.__value = self.__value + other
            else :
                out = DiffEq(self.getvnames())
                out.__value = ('+', self.__value, other)
        return out
    def __sub__(self, other) :
        if type(other) is DiffFunc :
            out = DiffEq(self.getvnames())
            out.__value = ('-', self.__value, other.__value)
        elif type(other) is DiffFunc or type(other) is DiffOp:
            out = DiffEq(self.getvnames())
            out.__value = ('-', self.__value, other)
        else :
            if type(self.__value) is not tuple :
                out = DiffEq(self.getvnames())
                out.__value = self.__value - other
            else :
                out = DiffEq(self.getvnames())
                out.__value = ('-', self.__value, other)
        return out

    def __mul__(self, other) :
        if type(other) is DiffFunc :
            out = DiffEq(self.getvnames())
            out.__value = ('*', self.__value, other.__value)
        elif type(other) is DiffFunc or type(other) is DiffOp:
            out = DiffEq(self.getvnames())
            out.__value = ('*', self.__value, other)
        else :
            if type(self.__value) is not tuple :
                out = DiffEq(self.getvnames())
                out.__value = self.__value * other
            else :
                out = DiffEq(self.getvnames())
                out.__value = ('*', self.__value, other)
        return out
def fem(

