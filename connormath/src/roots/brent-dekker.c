
#include "connormath.h"
#include <math.h>

double brent_dekker(double x0, double x1, double convergence, double tolerance, roots_func_t func, void *args) {
  double a = x0, b = x1, c, d, s, fa, fb, fc, fs;
  char flag;
  d = 0;

  fa = func(a, args);
  fb = func(b, args);

  if(fa == 0) {
    return (a);
  }
  if(fb == 0) {
    return (b);
  }
  
  if(fa * fb >= 0) {
    return (NAN);
  }

  if(abs(fa) < abs(fb)) {
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

  while(fb != 0 || abs(a - b) > abs(convergence)) {
    if(fa != fc && fb != fc) {
      s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
	(fa * b * fc) / ((fb - fa) * (fb - fc)) +
	(fa * fb * c) / ((fc - fa) * (fc - fb));
    } else {
      s = b - fb * (a - b) / (fa - fb);
    }
    
    if(!between((3 * a + b) / 4, b, s) ||
       (flag && (abs(s - b) / 2 >= abs(b - c) / 2) ||
	(abs(b - c) < abs(tolerance))) ||
       (!flag && (abs(s - b) / 2 >= abs(c - d) / 2) ||
	(abs(c - d) < abs(tolerance)))) {
      s = (a + b) / 2;
      flag = 1;
    } else {
      flag = 0;
    }
    fs = func(s, args);
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
    if(abs(fa) < abs(fb)) {
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

  while(abs(a - b) > convergence) {
    a = b;
    b = a - func(a, fargs) / deriv(a, derargs);
  }

  return (b);
}

double secant(double x0, double x1, double convergence, roots_func_t func, void *fargs) {
  double a = x0, b = x1, fa, fb;
  fa = func(a, fargs);
  fb = func(b, fargs);

  while(abs(x0 - x1) > convergence) {
    double temp = b;
    b = b - fb * (b - a) / (fb - fa);
    a = temp;
    fa = fb;
    fb = func(b, fargs);
  }
  return (b);
}

double regula_falsi(double x0, double x1, double convergence, roots_func_t func, void *args) {
  double a = x0, b = x1, c, fa, fb;
  fa = func(a, args);
  fb = func(b, args);

  if(fa * fb >= 0) {
    return (NAN);
  }

  if(abs(fa) < abs(fb)) {
    double temp = a;
    a = b;
    b = temp;
    temp = fa;
    fa = fb;
    fb = temp;
  }

  while(abs(a - b) > convergence) {
    c = b - fb * (b - a) / (fb - fa);
    fc = func(c, args);

    if(fa * fc >= 0) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }

    if(abs(fa) < abs(fb)) {
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
  
