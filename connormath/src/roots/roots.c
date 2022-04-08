
#include "connormath.h"
#include <math.h>
#include <errno.h>

double brent_dekker(double x0, double x1, double convergence, double tolerance, roots_func_t func, void *args) {
  double a = x0, b = x1, c, d, s, fa, fb, fc, fs;
  int err;
  char flag;
  d = 0;

  fa = func(a, args);
  if(isnan(fa)) {
    return (NAN);
  }
  fb = func(b, args);
  if(isnan(fb)) {
    return (NAN);
  }

  if(fa == 0) {
    return (a);
  }
  if(fb == 0) {
    return (b);
  }
  
  if(fa * fb >= 0) {
    // Domain error.
    errno = EDOM;
    return (NAN);
  }

  if(fabs(fa) < fabs(fb)) {
    double temp = a;
    a = b;
    b = temp;
    temp = fa;
    fa = fb;
    fb = temp;
  }

  c = a;
  fc = fa;
  flag = 1;

  while(fb != 0 || !relnear(a, b, convergence)) {
    if(fa != fc && fb != fc) {
      s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
	(fa * b * fc) / ((fb - fa) * (fb - fc)) +
	(fa * fb * c) / ((fc - fa) * (fc - fb));
    } else {
      s = b - fb * (a - b) / (fa - fb);
    }
    
    if(!between((3 * a + b) / 4, b, s) ||
       (flag && (fabs(s - b) / 2 >= fabs(b - c) / 2) ||
	(fabs(b - c) < fabs(tolerance))) ||
       (!flag && (fabs(s - b) / 2 >= fabs(c - d) / 2) ||
	(fabs(c - d) < fabs(tolerance)))) {
      s = (a + b) / 2;
      flag = 1;
    } else {
      flag = 0;
    }
    fs = func(s, args);
    if(isnan(fc)) {
      return (NAN);
    }
    d = c;
    c = b;
    fc = fb;

    if(fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if(fabs(fa) < fabs(fb)) {
      double temp = a;
      a = b;
      b = temp;
      temp = fa;
      fa = fb;
      fb = temp;
    }
  }
  return (b);
}
    
double newton(double x0, double convergence, roots_func_t func, void *fargs, roots_func_t deriv, void *derargs) {
  double a = x0, b;

  b = x0 - func(x0, fargs) / deriv(x0, derargs);

  while(relnear(a, b, convergence)) {
    a = b;
    b = a - func(a, fargs) / deriv(a, derargs);
  }

  return (b);
}

double secant(double x0, double x1, double convergence, roots_func_t func, void *fargs) {
  double a = x0, b = x1, fa, fb;
  fa = func(a, fargs);
  fb = func(b, fargs);

  while(relnear(a, b, convergence)) {
    double temp = b;
    b = b - fb * (b - a) / (fb - fa);
    a = temp;
    fa = fb;
    fb = func(b, fargs);
  }
  return (b);
}

double regula_falsi(double x0, double x1, double convergence, roots_func_t func, void *args) {
  double a = x0, b = x1, c, fa, fb, fc;
  fa = func(a, args);
  fb = func(b, args);

  if(fa * fb >= 0) {
    return (NAN);
  }

  if(fabs(fa) < fabs(fb)) {
    double temp = a;
    a = b;
    b = temp;
    temp = fa;
    fa = fb;
    fb = temp;
  }

  while(relnear(a, b, convergence)) {
    c = b - fb * (b - a) / (fb - fa);
    fc = func(c, args);

    if(fa * fc >= 0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }

    if(fabs(fa) < fabs(fb)) {
      double temp = a;
      a = b;
      b = temp;
      temp = fa;
      fa = fb;
      fb = temp;
    }
  }
  return (b);
}
  
