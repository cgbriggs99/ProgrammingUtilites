/*
 * beta.c
 *
 *  Created on: Jan 25, 2023
 *      Author: connor
 */

#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>
#include <float.h>
#include <ctype.h>


EXTRAMATH_FUNDEF(beta, (__TYPENAME__ __a, __TYPENAME__ __b)) {
	return __FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__a) +
			__FNAMESRC__(lgamma)(__b) - __FNAMESRC__(lgamma)(__a + __b));
}

EXTRAMATH_FUNDEF(betainc, (__TYPENAME__ __z, __TYPENAME__ __a, __TYPENAME__ __b)) {
  if(__a == 0 || ((int) __a < 0 && __FNAMESRC__(isnatural)(-__a) &&
		  !(__FNAMESRC__(isnatural)(__b) && (int) __b < (int) -__a))) {
    return INFINITY;
  } else if((int) __a < 0 && __FNAMESRC__(isnatural)(-__a) && __FNAMESRC__(isnatural)(__b) && (int) __b <= (int) -__a) {
    if(__z == 0) {
      return 0;  // The bounds of integration are the same, so zero not infinity.
    }
    __TYPENAME__ sum = 0, coef = 1, pz = 1;
    for(int i = 0; i <= (int) __b - 1; i++) {
      sum += coef * pz / (i - __a);
      coef *= (1 - __b + i);
      coef /= i + 1;
      pz *= __z;
    }
    return sum;
  } else {
    return __FNAMESRC__(pow)(__z, __a) * __FNAMESRC__(hyper2f1)(__a, 1 - __b, __a + 1, __z) / __a;
  }
}


// This one has some definitions, but they are not easy.
EXTRAMATH_FUNDEF(betareg, (__TYPENAME__ __z, __TYPENAME__ __a, __TYPENAME__ __b)) {
	if(__z == 0) {
		return 0;
	}
        return __FNAMESRC__(betainc)(__z, __a, __b) / __FNAMESRC__(beta)(__a, __b);
}

