#include "../include/extramath.h"
#include <math.h>
#include <stdio.h>
#include "test.h"


int test_gamma(void) {
  int warns = 0;

  double _Complex cx1 = 10.4; // Some random value.
  double x1 = 10.4;

  ASSERT_WARN(NEAR(ctgamma(cx1), tgamma(x1)), warns);
  ASSERT_WARN(NEAR(clgamma(cx1), lgamma(x1)), warns);
  ASSERT_WARN_MSG(NEAR(inGamma(x1, 1), 900608.86352224914), warns,
		  "Got %lf, expected %lf.\n",
		  inGamma(x1, 1), 900608.86352224914);

  return warns;
}

int main(void) {
  int warns = 0, errs = 0, res;

  res = test_gamma();
  if(res == -1) {
    errs++;
  } else {
    warns += res;
  }

  return warns + errs;
}
