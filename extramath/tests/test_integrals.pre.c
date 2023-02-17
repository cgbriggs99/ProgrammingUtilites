#include "${CMAKE_SOURCE_DIR}/extramath/include/extramath.h"
#include <${LAPACKE_HEADER}>
#include <math.h>
#include <stdio.h>

#define GOOD_ENOUGH 0.01
double finite_integrand(double x, const void *extra) {
  return x;
}

#define ALPHA 1
#define BETA 1
double jacobi_integrand(double x, const void *extra) {
  return x / (x - 1) / (x + 1);
}

double chebychev_integrand(double x, const void *extra) {
  return x * sqrt(1 - x * x);
}

double laguerre_integrand(double x, const void *extra) {
  return 1;
}

double gen_laguerre_integrand(double x, const void *extra) {
  return 1;
}

double hermite_integrand(double x, const void *extra) {
  return 1;
}

#define POINTS 20
int test_integrals(void) {
  int fails = 0;
  double ret;
  double expected1 = 0, expected2 = 1, expected3 = 0,
    expected4 = 1 / M_2_SQRTPI;

  ret = lriemannint(finite_integrand, -1, 1, POINTS, __INTEGRAL_BOTH__, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed lriemann, computed %lf.\n", ret);
    fails++;
  }

  ret = rriemannint(finite_integrand, -1, 1, POINTS, __INTEGRAL_BOTH__, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed rriemann, computed %lf.\n", ret);
    fails++;
  }

  ret = mriemannint(finite_integrand, -1, 1, POINTS, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed mriemann, computed %lf.\n", ret);
    fails++;
  }

  ret = trapezoidint(finite_integrand, -1, 1, POINTS, __INTEGRAL_BOTH__, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed trapezoid, computed %lf.\n", ret);
    fails++;
  }

  ret = simpsonint(finite_integrand, -1, 1, POINTS, __INTEGRAL_BOTH__, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed simpson, computed %lf.\n", ret);
    fails++;
  }

  ret = gausslegendreint(finite_integrand, -1, 1, POINTS, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed legendre, computed %lf.\n", ret);
    fails++;
  }
  
  ret = gaussjacobiint(jacobi_integrand, -1, 1, ALPHA, BETA, POINTS, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed jacobi, computed %lf.\n", ret);
    fails++;
  }

  ret = gausschebychevint1(chebychev_integrand, -1, 1, POINTS, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed chebychev1, computed %lf.\n", ret);
    fails++;
  }

  ret = gausschebychevint2(chebychev_integrand, -1, 1, POINTS, NULL);
  if(fabs(ret - expected3) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed chebychev2, computed %lf.\n", ret);
    fails++;
  }

  ret = gausslaguerreint(laguerre_integrand, 0, POINTS, NULL);
  if(fabs(ret - expected2) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed laguerre, computed %lf.\n", ret);
    fails++;
  }

  ret = gaussgenlaguerreint(laguerre_integrand, 0, ALPHA, POINTS, NULL);
  if(fabs(ret - expected2) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed general laguerre, computed %lf.\n", ret);
    fails++;
  }

  ret = gausshermiteint(hermite_integrand, POINTS, NULL);
  if(fabs(ret - expected4) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed hermite, computed %lf.\n", ret);
    fails++;
  }

  ret = tanhsinhint(finite_integrand, -1, 1, POINTS, NULL);
  if(fabs(ret - expected1) > GOOD_ENOUGH) {
    fprintf(stderr, "Failed tanhsinh, computed %lf.\n", ret);
    fails++;
  }
  
  return fails;
}



int main(void) {
  int fails = 0;
  fails += test_integrals();

  ${MKL_CLEAR};

  return fails;
}
