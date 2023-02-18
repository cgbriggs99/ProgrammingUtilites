
#include "extramath.h"
#include <math.h>
#include "test.h"

int test_js() {

  float x1f = 0, x2f = 1;
  float n1f = 2, n2f = -1, n3f = 1.5;
  double x1 = 0, x2 = 1;
  double n1 = 2, n2 = -1, n3 = 1.5;
  long double x1l = 0, x2l = 1;
  long double n1l = 2, n2l = -1, n3l = 1.5;
  float _Complex cx1f = 0, cx2f = 1, cx3f = I;
  float _Complex cn1f = 2, cn2f = -1, cn3f = 1.5;
  double _Complex cx1 = 0, cx2 = 1, cx3 = I;
  double _Complex cn1 = 2, cn2 = -1, cn3 = -1.5;
  long double _Complex cx1l = 0, cx2l = 1, cx3l = I;
  long double _Complex cn1l = 2, cn2l = -1, cn3l = 1.5;
  int warns = 0;
  int in1 = 0, in2 = 1;

  perror("Starting tests\n");
  ASSERT_WARN_MSG(CNEAR(cj0(cx1), 1), warns, "cj0(%lf + i%lf) = %lf + i%lf\n",
		  creal(cx1), cimag(cx1), creal(cj0(cx1)), cimag(cj0(cx1)))
  else {
    perror("Passed test 1.\n");
  }
  ASSERT_WARN_MSG(CNEAR(cj0(cx2), 0.765198), warns,
		  "cj0(%lf + i%lf) = %lf + i%lf\n",
		  creal(cx2), cimag(cx2), creal(cj0(cx2)), cimag(cj0(cx2)))
  else {
    perror("Passed test 2.\n");
  }
  ASSERT_WARN_MSG(CNEAR(cj0(cx3), 1.26607), warns,
		  "cj0(%lf + i%lf) = %lf + i%lf",
		  creal(cx3), cimag(cx3), creal(cj0(cx3)), cimag(cj0(cx3)))
  else {
    perror("Passed test 3.\n");
  }
  ASSERT_WARN_MSG(CNEAR(cj1(cx1), 0), warns, "cj1(%lf + i%lf) = %lf + i%lf",
		  creal(cx1), cimag(cx1), creal(cj1(cx1)), cimag(cj1(cx1)))
  else {
    perror("Passed test 4.\n");
  }
  ASSERT_WARN_MSG(CNEAR(cj1(cx2), 0.440051), warns,
		  "cj1(%lf + i%lf) = %lf + i%lf",
		  creal(cx2), cimag(cx2), creal(cj1(cx2)), cimag(cj1(cx2)))
  else {
    perror("Passed test 5.\n");
  }
  ASSERT_WARN_MSG(CNEAR(cj1(cx3), 0.565159 * I), warns,
		  "cj1(%lf + i%lf) = %lf + i%lf",
		  creal(cx3), cimag(cx3), creal(cj1(cx3)), cimag(cj1(cx3)))
  else {
    perror("Passed test 6.\n");
  }
  
  ASSERT_WARN_MSG(CNEAR(cjnu(cn1, cx1), 0), warns,
		  "cjnu(%lf + i%lf, %lf + i%lf) = %lf + i%lf",
		  creal(cn1), cimag(cn1), creal(cx1),
		  cimag(cx1), creal(cjnu(cn1, cx1)), cimag(cjnu(cn1, cx1)))
  else {
    perror("Passed test 7.\n");
  }
  ASSERT_WARN_MSG(CNEAR(cjnu(cn2, cx1), 0), warns,
		  "cjnu(%lf + i%lf, %lf + i%lf) = %lf + i%lf",
		  creal(cn2), cimag(cn2), creal(cx1),
		  cimag(cx1), creal(cjnu(cn2, cx1)), cimag(cjnu(cn2, cx1)))
  else {
    perror("Passed test 8.\n");
  }
  ASSERT_WARN_MSG(isinf(creal(cjnu(cn3, cx1))), warns,
		  "cjnu(%lf + i%lf, %lf + i%lf) = %lf + i%lf",
		  creal(cn3), cimag(cn3), creal(cx1),
		  cimag(cx1), creal(cjnu(cn3, cx1)), cimag(cjnu(cn3, cx1)))
  else {
    perror("Passed test 9.\n");
  }
  ASSERT_WARN_MSG(CNEAR(cjnu(cn1, cx2), 0.114903), warns,
		  "cjnu(%lf + i%lf, %lf + i%lf) = %lf + i%lf",
		  creal(cn1), cimag(cn1), creal(cx2),
		  cimag(cx2), creal(cjnu(cn1, cx2)), cimag(cjnu(cn1, cx2)))
  else {
    perror("Passed test 10.\n");
  }

  return warns;
}
  
int main(void) {
  int warns = 0, errs = 0;
  int res;

  res = test_js();
  if(res != -1) {
    warns += res;
  } else {
    errs++;
  }

  return warns + errs;
}
