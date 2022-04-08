
#ifndef __CM_CONNORMATH_H__
#define __CM_CONNORMATH_H__

#include "roots.h"

#define EQUAL(a, b, diff) ((a - b > 0)? (b - a < diff): (a - b < diff))

// Fails if equal to either
extern int betweenxx(double bound1, double bound2, double value);

// Fails if equal to first, passes if equal to second.
extern int betweenxi(double bound1, double bound2, double value);

// Fails if equal to second, passes if equal to first.
extern int betweenix(double bound1, double bound2, double value);

// Passes if equal to either.
extern int betweenii(double bound1, double bound2, double value);

// Uses the absolute difference over the average of the absolute values.
// Good for when the difference is outside the range of significance.
extern int relnear(double a, double b, double maxdiff);

// Uses the raw difference.
extern int absnear(double a, double b, double maxdiff);

#endif
