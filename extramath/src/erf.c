/*
 * erf.c
 *
 *  Created on: Jan 25, 2023
 *      Author: connor
 */


#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>
#include <float.h>

#ifdef __IS_COMPLEX__

EXTRAMATH_FUNDEF(erf, (__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = __z;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz / (2 * i + 1)); i++) {
		sum += coef * pz / (2 * i + 1);
		coef /= -(i + 1);
		pz *= __z * __z;
	}
	return M_2_SQRTPI * sum;
}

#endif

EXTRAMATH_FUNDEF(ierf, (__TYPENAME__ __z)) {
	if(__z == 0) {
		return 0;
#ifdef __IS_COMPLEX__
	} else if(__z == 1 || __z == -1) {
		return NAN;
#else
	} else if(__z >= 1 || __z <= -1) {
		return NAN;
#endif
	} else {
		__TYPENAME__ x0, x1;
#ifdef	__IS_COMPLEX__
		x0 = I;
#else
		x0 = 0;
#endif
		x1 = 1;
		while(!__FNAMESRC__(absconv)(x0, x1 - x0)) {
			x1 = x0;
			x0 = x1 - (__FNAMESRC__(erf)(x1) - __z) / (M_2_SQRTPI * __FNAMESRC__(exp)(-x1 * x1));
		}
		return x0;
	}
}


EXTRAMATH_FUNDEF(ierfc, (__TYPENAME__ __z)) {
	return __FNAMESRC__(ierf)(__z - 1);
}


