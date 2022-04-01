
#ifndef __CM_CONNORMATH_H__
#define __CM_CONNORMATH_H__

#define EQUAL(a, b, diff) ((a - b > 0)? (b - a < diff): (a - b < diff))

extern int between(double bound1, double bound2, double value);

#endif
