 
#include "../include/extramath.h"
#include "test.h"
#include <stdio.h>
#include <stdlib.h>

int test_poly(void) {
  int warns = 0;

  double coefs[] = {-1, -1, 1}, roots[2];
  double root1 = (1 + sqrt(5)) / 2, root2 = (1 - sqrt(5)) / 2;

  fprintf(stderr, "Testing general polynomials.\n");
  
  ASSERT_WARN(NEAR(polynomial(0.5, 3, coefs), -1.25), warns);
  polyroots(coefs, 2, roots);
  ASSERT_WARN_MSG((NEAR(roots[0], root1) && NEAR(roots[1], root2)) ||
	      (NEAR(roots[0], root2) && NEAR(roots[1], root1)),
		  warns, "Got roots %lf and %lf, expected roots %lf and %lf.\n",
		  roots[0], roots[1], root1, root2);
  ASSERT_WARN(NEAR(fbinom(10, 5), 252), warns);

  return warns;
}

int test_hermite(void) {
  int n = 10;
  double x = 1;
  double coefs[] = {-30240, 0, 302400, 0, -403200, 0,
    161280, 0, -23040, 0, 1024};
  double coef_test[11];
  int warns = 0;

  fprintf(stderr, "Testing Hermite polynomials.\n");

  ASSERT_WARN(NEAR(hermite(n, x), 8224), warns);
  ASSERT_WARN(NEAR(hermitederiv(n, x, 1), -214400), warns);
  ASSERT_WARN(hermitecof(n, coef_test) == 11, warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs[i], coef_test[i]);
  }
  return warns;
}

int test_laguerre(void) {
  int n = 10;
  double x = 1;
  double coefs[] = {1.0, -10.0, 22.5, -20.0, 8.75,
    -2.100027557319224, 0.2916666666666667,
    -0.023809523809523808, 0.0011160714285714285,
    -2.755731922398589e-05, 2.755731922398589e-07};
  double coef_test[11];
  int warns = 0;

  fprintf(stderr, "Testing Laguerre polynomials.\n");

  ASSERT_WARN(NEAR(laguerre(n, x), 168919.0 / 403200.0), warns);
  ASSERT_WARN(NEAR(laguerrederiv(n, x, 1), 396271.0 / 362880.0), warns);
  ASSERT_WARN(laguerrecof(n, coef_test) == 11, warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs[i], coef_test[i]);
  }
  return warns;
}

int test_assoclaguerre(void) {
  int n = 10;
  double x = 1, alpha = 0.5;
  double coefs[] = {3.70014, -24.6676, 44.4017, -33.8298, 13.156,
    -2.87041, 0.368001, -0.0280382, 0.00123698, -0.0000289352, 2.75573e-7};
  double coef_test[11];
  int warns = 0;

  fprintf(stderr, "Testing Associated Laguerre polynomials.\n");

  ASSERT_WARN(NEAR(assoclaguerre(n, alpha, x), 0.231181752307495), warns);
  ASSERT_WARN(NEAR(assoclaguerrederiv(n, alpha, x, 1), 2.93974), warns);
  ASSERT_WARN(assoclaguerrecof(n, alpha, coef_test) == 11, warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs[i], coef_test[i]);
  }
  return warns;
}

int test_legendre(void) {
  int n = 10;
  double x = 1;
  double coefs[] = {-63.0/256, 0, 3465.0/256, 0, -30030.0/256,
    0, 90090.0/256, 0, -109395.0/256, 0, 46189.0/256};
  double coef_test[11];
  int warns = 0;

  fprintf(stderr, "Testing Legendre polynomials.\n");

  ASSERT_WARN(NEAR(legendre(n, x), 1), warns);
  ASSERT_WARN(NEAR(legendrederiv(n, x, 1), 55), warns);
  ASSERT_WARN(NEAR(assoclegendre(n, 5, 0.5), 30086.169706116175), warns);
  ASSERT_WARN(legendrecof(n, coef_test) == 11, warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs[i], coef_test[i]);
  }
  return warns;
}

int test_chebychevt(void) {
  int n = 10;
  double x = 1;
  double coefs[] = {-1, 0, 50, 0, -400, 0, 1120, 0, -1280, 0, 512};
  double coef_test[11];
  int warns = 0;

  fprintf(stderr, "Testing Chebychev T polynomials.\n");

  ASSERT_WARN(NEAR(chebychevt(n, x), 1), warns);
  ASSERT_WARN_MSG(NEAR(chebychevtderiv(n, x, 1), 100), warns,
		  "Got %lf, expected %lf.\n",
		  chebychevtderiv(n, x, 1), 100);
  ASSERT_WARN(chebychevtcof(n, coef_test) == 11, warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs[i], coef_test[i]);
  }
  return warns;
}

int test_chebychevu(void) {
  int n = 10;
  double x = 1;
  double coefs[] = {-1, 0, 60, 0, -560, 0, 1792, 0, -2304, 0, 1024};
  double coef_test[11];
  int warns = 0;

  fprintf(stderr, "Testing Chebychev U polynomials.\n");

  ASSERT_WARN(NEAR(chebychevu(n, x), 11), warns);
  ASSERT_WARN(NEAR(chebychevuderiv(n, x, 1), 440), warns);
  ASSERT_WARN(chebychevucof(n, coef_test) == 11, warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs[i], coef_test[i]);
  }
  return warns;
}

int test_jacobi(void) {
  int n = 10;
  double x = 0, a = 1, b = 2;
  double coefs[] = {11, -770, 17325, -184800, 1099560, -3958416,
    8953560, -12790800, 11191950, -5471620, 1144066},
    coefs2[] = {66, -3080, 51975, -443520, 2199120, -6785856,
    13430340, -17054400, 13430340, -5969040, 1144066};
  double coef_test[11];
  int warns = 0;

  fprintf(stderr, "Testing Jacobi polynomials.\n");

  ASSERT_WARN_MSG(NEAR(jacobi(n, a, b, x), -231.0/512), warns,
		  "Got %lf, expected %lf.\n",
		  jacobi(n, a, b, x), -231.0/512);
  ASSERT_WARN_MSG(NEAR(jacobideriv(n, a, b, x, 1), -1155.0/256), warns,
		  "Got %lf, expected %lf.\n",
		  jacobideriv(n, a, b, x, 1), -1155.0/256);
  ASSERT_WARN(altjacobicof(n, a, b, coef_test) == 11, warns);
  ASSERT_WARN(NEAR(jacobi(n, a, b, x), jacobi(n, b, a, -x)), warns);
  ASSERT_WARN(NEAR(jacobideriv(n, a, b, x, 1), -jacobideriv(n, b, a, -x, 1)), warns);
  ASSERT_WARN(NEAR(jacobi(n + 1, a, b, x), -jacobi(n + 1, b, a, -x)), warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs[i], coef_test[i]);
  }
  ASSERT_WARN(altjacobicof(n, b, a, coef_test) == 11, warns);

  for(int i = 0; i < 11; i++) {
    ASSERT_WARN_MSG(NEAR(coefs2[i], coef_test[i]), warns,
		    "Coefficients not equal, expected %lf, got %lf.\n",
		    coefs2[i], coef_test[i]);
  }
  return warns;
}

int main(void) {
  int warns = 0, errs = 0, ret;

  ret = test_poly();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }
  ret = test_hermite();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }
  ret = test_laguerre();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }
  ret = test_assoclaguerre();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }
  ret = test_legendre();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }
  ret = test_chebychevt();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }
  ret = test_chebychevu();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }
  ret = test_jacobi();
  if(ret != -1) {
    warns += ret;
  } else {
    errs++;
  }

  return warns + errs;
}
