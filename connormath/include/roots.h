#ifndef __CM_ROOTS_H__
#define __CM_ROOTS_H__

typedef double (*roots_func_t)(double x, void *args);

// Brent-Dekker method for finding roots. a and b need to bracket the root.
extern double brent_dekker(double x0, double x1, double convergence, double tolerance, roots_func_t func, void *args);

extern double newton(double x0, double convergence, roots_func_t func, void *fargs, roots_func_t deriv, void *derargs);

extern double secant(double x0, double x1, double convergence, roots_func_t func, void *fargs);

extern double regula_falsi(double x0, double x1, double convergence, roots_func_t func, void *args);

extern int poly_roots(double *coefs, int len, double *out, double conv);

#endif
