#include "connormath.h"
#include "test.h"
#include <math.h>
#include <stdio.h>

static double coefs[] = { -9694845, 0, 4508102925, 0,
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
    sum += coefs[i] * xp;
    xp *= pos;
  }
  return (sum);
}

double func2(double x, void *coefs) {
  return (sin(x));
}

#define CONV 1e-6
#define TOL 1e-4
int test_brent_dekker(void) {
  double ret0, ret1, ret2, ret3;
  int err0, err1, err2, err3, retval = 0;
  // Test first two if statements
  ret0 = brent_dekker(0, 1, CONV, TOL, func2, NULL);
  err0 = errno;
  ret1 = brent_dekker(1, 0, CONV, TOL, func2, NULL);
  err1 = errno;

  // Test bracketing.
  ret2 = brent_dekker(1, 2, CONV, TOL, func2, NULL);
  err2 = errno;

  // Test convergence.
  ret3 = brent_dekker(-1, 0.5, CONV, TOL, func2, NULL);
  err3 = errno;

  if(ret0 != 0) {
    retval = -1;
    errno = err0;
    perror("First test failed.");
  }
  if(ret1 != 0) {
    retval = -1;
    errno = err1;
    perror("Second test failed.");
  }
  if(!isnan(ret2) || err2 != EDOM) {
    retval = -1;
    errno = err2;
    perror("Third test failed.");
  }
  if(!absnear(ret3, 0, CONV)) {
    retval = -1;
    errno = err3;
    perror("Fourth test failed.");
  }
  return (retval);
}

test_func_t funcs[] = {
		       test_brent_dekker,
		       NULL
};

test_struct_t init_tests() {
  return (setup_tests(funcs));
}
