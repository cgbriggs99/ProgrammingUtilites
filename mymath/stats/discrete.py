#!/usr/bin/python3

from . import stats
import math

def mean(data, weights = None) :
    """
Calculate the mean of the data set. Can have weights passed to calculate a
weighted average.

data: The data to use to find the average.
weights: An iterable of weights, or None to use an unweighted average. The
weights do not need to add to 1. The sum will be normalized by the total of
the weights.
"""
    if weights == None :
        return sum(data) / len(data)
    else :
        return sum(weights[i] * data[i] for i in range(len(data))) / \
               sum(weights);
    
def median(data, weights = None, do_qs = False) :
    """
Find the median of a set of data. Can optionally compute the quartiles.

data: The data to use to find the median.
do_qs: True calculates the quartiles, returning in the order
Q2 (median), Q1, Q3. False (default) just calculates the median.
"""
    if weights == None :
        s = sorted(data)
        if len(s) % 2 == 1 :
            q2 = s[len(s) // 2]
            if do_qs :
                q1 = median(s[:len(s) // 2], None, False)
                q3 = median(s[len(s) // 2 + 1:], None, False)
                return q2, q1, q3
            else :
                return q2
        else :
            q2 = (s[len(s) // 2 - 1] + s[len(s) // 2]) / 2
            if do_qs :
                q1 = median(s[:len(s) // 2], None, False)
                q3 = median(s[len(s) // 2:], None, False)
                return q2, q1, q3
            else :
                return q2
    else :
        norm = sum(weights)
        pairs = list(zip(data, weights))
        pairs = sorted(pairs, key = lambda x: x[0])
        s0 = 0
        s1 = 0
        q1 = 0
        q2 = 0
        q3 = 0
        for i in range(len(pairs)) :
            s0 = s1
            s1 += pairs[i][1] / norm
            if s1 == 0.25 :
                q1 = s1
            elif s1 > 0.25 and s0 < 0.25 :
                q1 = (0.25 - s0) * (pairs[i][0] - pairs[i - 1][0]) /\
                     (s1 - s0) + pairs[i - 1][0]
            elif s1 == 0.5 :
                q2 = s1
            elif s1 > 0.5 and s0 < 0.5 :
                q2 = (0.5 - s0) * (pairs[i][0] - pairs[i - 1][0]) /\
                     (s1 - s0) + pairs[i - 1][0]
            elif s1 == 0.75 :
                q3 = s1
            elif s1 > 0.75 and s0 < 0.75 :
                q3 = (0.75 - s0) * (pairs[i][0] - pairs[i - 1][0]) /\
                     (s1 - s0) + pairs[i - 1][0]
        if do_qs :
            return q2, q1, q3
        else :
            return q2

def mode(data) :
    counts = {}
    for d in data :
        if d in counts :
            counts[d] += 1
        else :
            counts[d] = 1
    return max(counts)

def svariance(data, weights = None) :
    """
Find the sample variance of a set of data. Arguments are the same as in mean.
See stats.mean for more info. stats.variance is an alias for this function.
"""
    m = mean(data, weights)
    if weights == None :
        return sum((d - m) ** 2 for d in data) / (len(data) - 1)
    else :
        return sum((d - m) ** 2 * w for d, w in zip(data, weights)) / sum(weights)

def pvariance(data, weights = None) :
    """
Find the population variance of a set of data. See stats.mean for arguments.
"""
    m = mean(data, weights)
    if weights == None :
        return sum((d - m) ** 2 for d in data) / len(data)
    else :
        return sum((d - m) ** 2 * w for d, w in
                   zip(data, weights)) / sum(weights)

variance = svariance

def sstddev(data, weights = None) :
    """
Find the sample standard deviation of a set of data. See stats.mean for
arguments. stats.stddev is an alias for this function.
"""
    return math.sqrt(svariance(data, weights))

def pstddev(data, weights = None) :
    """
Find the population standard deviation of a set of data. See stats.mean for
arguments.
"""
    return math.sqrt(pvariance(data, weights))

stddev = sstddev

def moment(data, order, weights = None) :
    if weights == None :
        return sum(d ** order for d in data) / len(data)
    else :
        return sum(d ** order * w for d, w in zip(data, weights)) / \
               sum(weights)

def cmoment(data, order, weights = None) :
    m = mean(data, weights)
    if weights == None :
        return sum((d - m) ** order for d in data) / (len(data) - 1)
    else :
        return sum((d - m) ** order * w for d, w in zip(data, weights)) / \
               sum(weights)

def entropy(data, weights = None) :
    if weights is None :
        return math.log(len(data))
    else :
        norm = sum(weights)
        return -sum(w * math.log(w / norm) / norm for w in weights)

class Discrete(stats.Distribution) :
    def __init__(self, data, weights = None) :
        self.__data = data
        self.__weights = weights
        self.mean = mean(data, weights)
        self.median = median(data, weights)
        self.mode = mode(data)
        self.variance = variance(data, weights)
        self.skew = cmoment(data, 3, weights) / self.variance ** (3 / 2)
        self.exkurtosis = cmoment(data, 4, weights) / self.variance ** 2 - 3
        self.entropy = entropy(data, weights)
        self.range = max(data) - min(data)
    def pdf(self, x) :
        if self.__weights is None :
            return 1 / len(self.__data)
        if x in self.__data :
            return self.__weights[self.__data.index(x)]
        else :
            highi = 0
            lowi = 0
            low = math.inf
            high = -math.inf
            for i in range(len(self.__data)) :
                if self.__data[i] < x and self.__data[i] > low :
                    low = self.__data[i]
                    lowi = i
                elif self.__data[i] > x and self.__data[i] < high :
                    high = self.__data[i]
                    highi = i
            return (self.__weights[highi] - self.__weights[lowi]) /\
                   (self.__data[highi] - self.__data[lowi]) *\
                   (x - self.__data[lowi]) + self.__weights[lowi]
    def cdf(self, x) :
        if self.__weights is None :
            return sum(1 for d in self.__data if d <= x) / len(self.__data)
        else :
            norm = sum(self.__weights)
            return sum(w / norm for w, d in
                       zip(self.__weights, self.__data) if d <= x)
    def moment(self, order) :
        return moment(self.__data, order, self.__weights)
    def cmoment(self, order) :
        return cmoment(self.__data, order, self.__weights)
    def quantile(self, x) :
        if self.__weights is None :
            weights = [1 / (len(self.__data) + 1) for i in range(len(self.__data))]
            norm = 1
        else :
            weights = self.__weights
            norm = sum(weights)
        pairs = list(zip(self.__data, weights))
        pairs = sorted(pairs, key = lambda x: x[0])
        s0 = 0
        s1 = 0
        for i in range(len(pairs)) :
            s0 = s1
            s1 += pairs[i][1] / norm
            if s1 == x :
                return pairs[i][0]
            elif s1 > x and s0 < x and i != 0:
                return ((x - s0) * pairs[i][0] + (s1 - x) * pairs[i - 1][0]) /\
                    (s1 - s0)
            elif s1 > x and i == 0 :
                return pairs[0][0]
        return pairs[-1][0]
    def test(self, x) :
        return 2 * (x - self.median) / (self.quantile(0.75) -
                                        self.quantile(0.25))
    def stat(self, conf) :
        return 1.25
