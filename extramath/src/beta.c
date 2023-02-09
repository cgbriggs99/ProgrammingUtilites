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



#ifndef INTPOINTS
#define INTPOINTS 100
#endif

static __TYPENAME__ betakern(__SCALARTYPE__ __t, const void *data) {
	const __TYPENAME__ *params = (const __TYPENAME__ *) data;

	return (__FNAMESRC__(pow)(1 - __t * params[2], params[0]) * __FNAMESRC__(pow)(1 + __t * params[2], params[1]));
}

// This one has some definitions, but they are not easy.
EXTRAMATH_FUNDEF(betainc, (__TYPENAME__ __z, __TYPENAME__ __a, __TYPENAME__ __b)) {
	if(__z == 0) {
		return 0;
	}
	__TYPENAME__  params[3] = {__a, __b, __z};
	return __FNAMESRC__(trapezoidint)(betakern, 0, 1, INTPOINTS, __INTEGRAL_BOTH__, (void *) params);
}


// This one has some definitions, but they are not easy.
EXTRAMATH_FUNDEF(betareg, (__TYPENAME__ __z, __TYPENAME__ __a, __TYPENAME__ __b)) {
	if(__z == 0) {
		return 0;
	}
	__TYPENAME__ params[3] = {__a, __b, __z};
	return __FNAMESRC__(trapezoidint)(betakern, 0, 1, INTPOINTS, __INTEGRAL_BOTH__, (void *) params) /
			__FNAMESRC__(beta)(__a, __b);
}

