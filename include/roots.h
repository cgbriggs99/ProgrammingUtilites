#ifndef __CM_ROOTS_H__
#define __CM_ROOTS_H__

typedef double (*roots_func_t)(double x, void *args);

extern double brent_dekker(double a, double b, roots_func_t func, void *args);


#endif
