#include "connormath.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

int test_utils(void) {

  assert(!betweenxx(0, 1, 1));
  assert(!betweenxx(0, 1, 0));
  assert(betweenxx(0, 1, nexttoward(1, 0)));
  assert(betweenxx(0, 1, nexttoward(0, 1)));
  assert(betweenxx(0, 1, 0.5));
  assert(!betweenxx(0, 1, 2));
  assert(!betweenxx(0, 1, -1));

  assert(betweenxi(0, 1, 1));
  assert(!betweenxi(0, 1, 0));
  assert(betweenxi(0, 1, nexttoward(1, 0)));
  assert(betweenxi(0, 1, nexttoward(0, 1)));
  assert(betweenxi(0, 1, 0.5));
  assert(!betweenxi(0, 1, 2));
  assert(!betweenxi(0, 1, -1));

  assert(!betweenix(0, 1, 1));
  assert(betweenix(0, 1, 0));
  assert(betweenix(0, 1, nexttoward(1, 0)));
  assert(betweenix(0, 1, nexttoward(0, 1)));
  assert(betweenix(0, 1, 0.5));
  assert(!betweenix(0, 1, 2));
  assert(!betweenix(0, 1, -1));

  assert(betweenii(0, 1, 1));
  assert(betweenii(0, 1, 0));
  assert(betweenii(0, 1, nexttoward(1, 0)));
  assert(betweenii(0, 1, nexttoward(0, 1)));
  assert(betweenii(0, 1, 0.5));
  assert(!betweenii(0, 1, 2));
  assert(!betweenii(0, 1, -1));

  assert(relnear(1e-12, -1e-13, 1e-7));
  assert(!relnear(1, 0, 1e-7));
  assert(relnear(1, 1, 1e-7));

  assert(absnear(1e-12, -1e-13, 1e-7));
  assert(!absnear(1, 0, 1e-7));
  assert(absnear(1, 1, 1e-7));

  return (0);
}

int main(void) {
  test_utils();
  printf("All tests passed!\n");

  return (0);
}
