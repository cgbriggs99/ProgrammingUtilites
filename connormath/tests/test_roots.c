#include "connormath.h"
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>

static double pcoefs[] = { -9694845, 0, 4508102925, 0,
			  -347123925225, 0, 10529425731825, 0,
			  -166966608033225, 0, 1591748329916745, 0,
			  -9888133564634325, 0, 42051732851796525, 0,
			  -126155198555389575, 0, 271274904083157975, 0,
			  -419762220002360235, 0, 463373879223384675, 0,
			  -355924863751295475, 0, 180700315442965395, 0,
			  -54496920530418135, 0, 7391536347803839};
static int len = 31;

double func1(double x, void *coefs) {
  double sum = 0, xp = 1;
  for(int i = 0; i < len; i++) {
    sum += ((double *) coefs)[i] * xp;
    xp *= x;
  }
  return (sum);
}

double func2(double x, void *args) {
  return (sin(x));
}

double func2_deriv(double x, void *args) {
  return (cos(x));
}

double func2_dderiv(double x, void *args) {
  return (-sin(x));
}

double func3(double x, void *args) {
  return (1 / x);
}

double func3_deriv(double x, void *args) {
  return (-1 / x / x);
}

double func4(double x, void *args) {
  return (sin(x) + 2);
}

double func4_deriv(double x, void *args) {
  return (cos(x));
}

#define CONV 1e-6
#define TOL 1e-4
int test_brent_dekker(void) {
  double ret0, ret1, ret2, ret3, ret4;
  int err0, err1, err2, err3, err4, retval = 0;
  // Test first two if statements
  errno = 0;
  ret0 = brent_dekker(0, 1, CONV, TOL, func2, NULL);
  err0 = errno;
  errno = 0;
  ret1 = brent_dekker(1, 0, CONV, TOL, func2, NULL);
  err1 = errno;

  // Test bracketing.
  errno = 0;
  ret2 = brent_dekker(1, 2, CONV, TOL, func2, NULL);
  err2 = errno;

  // Test convergence.
  errno = 0;
  ret3 = brent_dekker(-1, 0.5, CONV, TOL, func2, NULL);
  err3 = errno;

  // Test failing.
  errno = 0;
  ret4 = brent_dekker(-10, 20, CONV, TOL, func3, NULL);
  err4 = errno;

  errno = err0;
  assert(ret0 == 0);
  errno = err1;
  assert(ret1 == 0);
  errno = err2;
  assert(isnan(ret2) && err2 == EDOM);
  errno = err3;
  assert(absnear(ret3, 0, CONV));
  errno = err4;
  assert(isnan(ret4) && err4 == ERANGE);
  errno = 0;
  
  return (retval);
}

int test_newton(void) {
  double ret0, ret1, ret2;
  int err0, err1, err2;

  // Reset the error.
  errno = 0;
  // Test a failing derivative.
  ret0 = newton(0, CONV, func2_deriv, NULL, func2_dderiv, NULL);
  err0 = errno;

  // Test convergence
  errno = 0;
  ret1 = newton(1, CONV, func2, NULL, func2_deriv, NULL);
  err1 = errno;

  // Test failure to converge.
  errno = 0;
  ret2 = newton(1, CONV, func4, NULL, func4_deriv, NULL);
  err2 = errno;

  errno = err0;
  assert(isnan(ret0) && err0 == EDOM);
  errno = err1;
  assert(relnear(ret1, 0, CONV) && err1 == 0);
  errno = err2;
  assert(isnan(ret2) && err2 == ERANGE);
  errno = 0;
  return (0);
}

int test_secant(void) {
  double ret0, ret1;
  double err0, err1;

  //Test convergence.
  errno = 0;
  ret0 = secant(-1, 0.5, CONV, func2, NULL);
  err0 = errno;

  //Test error.
  errno = 0;
  ret1 = secant(1, 2, CONV, func3, NULL);
  err1 = errno;

  errno= err0;
  assert(relnear(ret0, 0, CONV) && err0 == 0);
  errno = err1;
  assert(isnan(ret1) && err1 == ERANGE);
  errno = 0;

  return (0);
}

int test_regula_falsi(void) {
  double ret0, ret1, ret2, ret3;
  int err0, err1, err2, err3;

  // Test bracketing.
  errno = 0;
  ret0 = regula_falsi(1, 2, CONV, func2, NULL);
  err0 = errno;

  // Test swapping and convergence.
  errno = 0;
  ret1 = regula_falsi(-0.5, 1, CONV, func2, NULL);
  err1 = errno;

  // Test failure.
  errno = 0;
  ret2 = regula_falsi(-10, 20, CONV, func3, NULL);
  err2 = errno;

  errno = err0;
  assert(isnan(ret0) && err0 == EDOM);
  errno = err1;
  assert(relnear(ret1, 0, CONV) && err1 == 0);
  errno = err2;
  assert(isnan(ret2) && err2 == ERANGE);
  errno = 0;
  
  return (0);
}

int contains_all(double *coll1, int len1, double *coll2, int len2, double conv) {
  if(len1 != len2) {
    return (0);
  }
  for(int i = 0; i < len1; i++) {
    int found = 0;
    for(int j = 0; j < len2; j++) {
      if(relnear(coll1[i], coll2[j], conv)) {
	found = 1;
	break;
      }
    }
    if(!found) {
      return (0);
    }
  }
  return (1);
}

int test_poly_roots(void) {
  double coefs0[] = {1}, coefs1[] = {1, 1},
    coefs2[] = {0, 1, 1}, coefs3[] = {-1, 0, 1},
    coefs4[] = {1, 2, 1}, coefs5[] = {2, 3, 1},
    coefs6[] = {1, 1, 0}, coefs7[] = {-2, -1, 2, 1},
    coefs8[] = {2, 5, 4, 1}, coefs9[] = {1, 3, 3, 1},
    coefs10[] = {1, 0, 0, 1}, coefs11[] = {24, 50, 35, 10, 1},
    coefs12[] = {18, 39, 29, 9, 1}, coefs13[] = {8, 20, 18, 7, 1},
    coefs14[] = {1, 4, 6, 4, 1}, coefs15[] = {2, 3, 3, 3, 1},
    coefs16[] = {1, 2, 2, 2, 1}, coefs17[] = {6, 6, -5, -5, 1, 1},
    coefs18[] = {6, 11, 12, 12, 6, 1}, coefs19[] = {3, 7, 8, 8, 5, 1},
    coefs20[] = {6, 2, 9, 3, 3, 1}, coefs21[] = {1, 0, -2, 0, 1},
    coefs22[] = {1, 0, 0, 0, 0, 0, 1}, coefs23[] = {1, 5, 10, 10, 5, 1};
  double res0[] = {}, res1[] = {-1}, res2[] = {0, -1}, res3[] = {1, -1},
    res4[] = {-1, -1}, res5[] = {-1, -2}, res6[] = {-1}, res7[] = {-1, 1, -2},
    res8[] = {-1, -1, -2}, res9[] = {-1, -1, -1}, res10[] = {-1},
    res11[] = {-1, -2, -3, -4}, res12[] = {-1, -2, -3, -3},
    res13[] = {-2, -2, -2, -1}, res14[] = {-1, -1, -1, -1},
    res15[] = {-1, -2}, res16[] = {-1, -1},
    res17[] = {-1, -sqrt(2), sqrt(2), -sqrt(3), sqrt(3)},
    res18[] = {-1, -2, -3}, res19[] = {-3, -1, -1}, res20[] = {-3},
    res21[] = {-1, -1, 1, 1}, res22[] = {}, res23[] = {-1, -1, -1, -1, -1};
  double *ret0, *ret1, *ret2, *ret3, *ret4, *ret5, *ret6, *ret7, *ret8, *ret9,
    *ret10, *ret11, *ret12, *ret13, *ret14, *ret15, *ret16, *ret17, *ret18,
    *ret19, *ret20, *ret21, *ret22, *ret23;
  int err0, err1, err2, err3, err4, err5, err6, err7, err8, err9,
    err10, err11, err12, err13, err14, err15, err16, err17, err18, err19,
    err20, err21, err22, err23;
  int herr0, herr1, herr2, herr3, herr4, herr5, herr6, herr7, herr8, herr9,
    herr10, herr11, herr12, herr13, herr14, herr15, herr16, herr17, herr18,
    herr19, herr20, herr21, herr22, herr23;

  // Test for constant.
  errno = 0;
  ret0 = calloc(0, sizeof(double));
  herr0 = poly_roots(coefs0, 1, ret0, CONV);
  err0 = errno;

  // Test for 
  errno = 0;
  ret1 = calloc(1, sizeof(double));
  herr1 = poly_roots(coefs1, 2, ret1, CONV);
  err1 = errno;

  // Test for 
  errno = 0;
  ret2 = calloc(2, sizeof(double));
  herr2 = poly_roots(coefs2, 3, ret2, CONV);
  err2 = errno;

  // Test for 
  errno = 0;
  ret3 = calloc(2, sizeof(double));
  herr3 = poly_roots(coefs3, 3, ret3, CONV);
  err3 = errno;

  // Test for 
  errno = 0;
  ret4 = calloc(2, sizeof(double));
  herr4 = poly_roots(coefs4, 3, ret4, CONV);
  err4 = errno;

  // Test for 
  errno = 0;
  ret5 = calloc(2, sizeof(double));
  herr5 = poly_roots(coefs5, 3, ret5, CONV);
  err5 = errno;

  // Test for 
  errno = 0;
  ret6 = calloc(2, sizeof(double));
  herr6 = poly_roots(coefs6, 3, ret6, CONV);
  err6 = errno;

  // Test for 
  errno = 0;
  ret7 = calloc(3, sizeof(double));
  herr7 = poly_roots(coefs7, 4, ret7, CONV);
  err7 = errno;

  // Test for 
  errno = 0;
  ret8 = calloc(3, sizeof(double));
  herr8 = poly_roots(coefs8, 4, ret8, CONV);
  err8 = errno;

  // Test for 
  errno = 0;
  ret9 = calloc(3, sizeof(double));
  herr9 = poly_roots(coefs9, 4, ret9, CONV);
  err9 = errno;

  // Test for 
  errno = 0;
  ret10 = calloc(3, sizeof(double));
  herr10 = poly_roots(coefs10, 4, ret10, CONV);
  err10 = errno;

  // Test for 
  errno = 0;
  ret11 = calloc(4, sizeof(double));
  herr11 = poly_roots(coefs11, 5, ret11, CONV);
  err11 = errno;

  // Test for 
  errno = 0;
  ret12 = calloc(4, sizeof(double));
  herr12 = poly_roots(coefs12, 5, ret12, CONV);
  err12 = errno;

  // Test for 
  errno = 0;
  ret13 = calloc(4, sizeof(double));
  herr13 = poly_roots(coefs13, 5, ret13, CONV);
  err13 = errno;

  // Test for 
  errno = 0;
  ret14 = calloc(4, sizeof(double));
  herr14 = poly_roots(coefs14, 5, ret14, CONV);
  err14 = errno;

  // Test for 
  errno = 0;
  ret15 = calloc(4, sizeof(double));
  herr15 = poly_roots(coefs15, 5, ret15, CONV);
  err15 = errno;

  // Test for 
  errno = 0;
  ret16 = calloc(4, sizeof(double));
  herr16 = poly_roots(coefs16, 5, ret16, CONV);
  err16 = errno;

  // Test for 
  errno = 0;
  ret17 = calloc(5, sizeof(double));
  herr17 = poly_roots(coefs17, 6, ret17, CONV);
  err17 = errno;

  // Test for 
  errno = 0;
  ret18 = calloc(5, sizeof(double));
  herr18 = poly_roots(coefs18, 6, ret18, CONV);
  err18 = errno;

  // Test for 
  errno = 0;
  ret19 = calloc(5, sizeof(double));
  herr19 = poly_roots(coefs19, 6, ret19, CONV);
  err19 = errno;

  // Test for 
  errno = 0;
  ret20 = calloc(5, sizeof(double));
  herr20 = poly_roots(coefs20, 6, ret20, CONV);
  err20 = errno;

  // Test for 
  errno = 0;
  ret21 = calloc(4, sizeof(double));
  herr21 = poly_roots(coefs21, 5, ret21, CONV);
  err21 = errno;

  // Test for 
  errno = 0;
  ret22 = calloc(6, sizeof(double));
  herr22 = poly_roots(coefs22, 7, ret22, CONV);
  err22 = errno;

  // Test for 
  errno = 0;
  ret23 = calloc(5, sizeof(double));
  herr23 = poly_roots(coefs23, 6, ret23, CONV);
  err23 = errno;

  errno = err0;
  assert(herr0 == sizeof(res0) / sizeof(double) && err0 == EDOM);
  
#define CHECK(n) errno = err ## n; assert(herr ## n == sizeof(res ## n) / sizeof(double) && err ## n == 0 && contains_all(res ## n, sizeof(res ## n) / sizeof(double), ret ## n, herr ## n, CONV)); free(ret ## n)

  CHECK(1);
  CHECK(2);
  CHECK(3);
  CHECK(4);
  CHECK(5);
  CHECK(6);
  CHECK(7);
  CHECK(8);
  CHECK(9);
  CHECK(10);
  CHECK(11);
  CHECK(12);
  CHECK(13);
  CHECK(14);
  CHECK(15);
  CHECK(16);
  CHECK(17);
  CHECK(18);
  CHECK(19);
  CHECK(20);
  CHECK(21);
  CHECK(22);
  errno = err23;
  for(int i = 0; i < herr23; i++) {
    assert(fabs(1 + ret23[i]) < 0.001);
  }
  free(ret23);

  free(ret0);
  return (0);
}
  
  
int main(void) {
  test_brent_dekker();
  printf("Passed Brent-Dekker\n");
  test_newton();
  printf("Passed Newton\n");
  test_secant();
  printf("Passed Secant\n");
  test_regula_falsi();
  printf("Passed Regula Falsi\n");
  test_poly_roots();
  printf("Passed Polynomial\n");
  printf("Passed all tests.\n");

  return (0);
}
