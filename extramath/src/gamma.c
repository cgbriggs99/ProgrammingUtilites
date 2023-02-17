/*
 * integrated.c
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>
#include <float.h>

#ifdef __IS_COMPLEX__

EXTRAMATH_FUNDEF(lgamma, (__TYPENAME__ __z)) {

}

// Use Lanczos approximation.
EXTRAMATH_FUNDEF(tgamma, (__TYPENAME__ __z)) {

  
}

#endif

// Should work for complex once tgamma is defined.
EXTRAMATH_FUNDEF(inGamma, (__TYPENAME__ __a, __TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1 / __a, pz = __FNAMESRC__(pow)(__z,__a);

	for(int i = 0; !__FNAMESRC__(absconv)(sum, pz * coef); i++) {
		sum += coef * pz;
		coef *= (__a + i) * (-__z) / ((__a + i + 1) * (i + 1));
	}
	return (__FNAMESRC__(tgamma)(__a) - sum);
}

// Should work for complex.
EXTRAMATH_FUNDEF(ingamma, (__TYPENAME__ __a, __TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1 / __a, pz = __FNAMESRC__(pow)(__z,__a);

	for(int i = 0; !__FNAMESRC__(absconv)(pz * coef, sum); i++) {
		sum += coef * pz;
		coef *= (__a + i) * (-__z) / ((__a + i + 1) * (i + 1));
	}
	return (sum);
}

EXTRAMATH_FUNDEF(ingamma2, (__TYPENAME__ __a, __TYPENAME__ __z1, __TYPENAME__ __z2)) {
	return (__FNAMESRC__(ingamma)(__a, __z1) - __FNAMESRC__(ingamma)(__a, __z2));
}

EXTRAMATH_FUNDEF(inGammareg, (__TYPENAME__ __a, __TYPENAME__ __z)) {
	return (__FNAMESRC__(inGamma)(__a, __z) / __FNAMESRC__(tgamma)(__a));
}

EXTRAMATH_FUNDEF(ingammareg, (__TYPENAME__ __a, __TYPENAME__ __z)) {
	return (__FNAMESRC__(ingamma)(__a, __z) / __FNAMESRC__(tgamma)(__a));
}

EXTRAMATH_FUNDEF(ingamma2reg, (__TYPENAME__ __a, __TYPENAME__ __z1, __TYPENAME__ __z2)) {
	return (__FNAMESRC__(ingamma2)(__a, __z1, __z2) / __FNAMESRC__(tgamma)(__a));
}

EXTRAMATH_FUNDEF(subfact, (__TYPENAME__  __z)) {
	return (__FNAMESRC__(inGamma)(__z + 1, -1) / M_E);
}

// Should work for complex, does not rely on gamma functions.
EXTRAMATH_FUNDEF(polygamma, (__TYPENAME__ __z)) {
#	ifndef __IS_COMPLEX__
	// Can't take value when Gamma(z) < 0.
	if(__z < 0 && __FNAMESRC__(fmod)(-__z, 2) < 1) {
		return NAN;
	} else if(__z <= 0 && __FNAMESRC__(fmod)(-__z, 1) == 0) {
		return INFINITY;
	}
#	else
	if(__FNAMESRC__(real)(__z) < 0 && __FNAMESRC__(imag)(__z) == 0 &&
			__FNAMESRC_SCAL__(fmod)(-__FNAMESRC__(real)(__z), 2) < 1) {
		return NAN;
	} else if (__FNAMESRC__(real)(__z) <= 0 && __FNAMESRC__(imag)(__z) == 0 &&
				__FNAMESRC_SCAL__(fmod)(-__FNAMESRC__(real)(__z), 2) == 0) {
		return INFINITY;
	}
#	endif
	__TYPENAME__ sum = 0, coef = 1 / __z;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef); i++) {
		coef = 1.0 / (i + 1) / (i + __z + 1);
		sum += coef;
	}

	return (-1.0 / __z + __z * sum - MASCHERONI);
}

EXTRAMATH_FUNDEF(polygamman, (int __n, __TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

#	ifndef __IS_COMPLEX__
	// Can't take value when Gamma(z) < 0.
	if(__z < 0 && __FNAMESRC__(fmod)(-__z, 2) < 1) {
		return NAN;
	} else if(__z <= 0 && __FNAMESRC__(fmod)(-__z, 1) == 0) {
		return INFINITY;
	} else if(__n == -1) {	// This is the zeroth derivative of ln(Gamma).
		return __FNAMESRC__(lgamma)(__z);
	} else if(__n == 0) {
		return __FNAMESRC__(polygamma)(__z);
	} else if(__n > 0) {
		coef = 1 / __FNAMESRC__(pow)(__z, __n + 1);
		for(int i = 0; __FNAMESRC__(fabs)(coef) > __FNAMESRC__(epsilon)(sum); i++) {
			sum += coef;
			coef = 1 / __FNAMESRC__(pow)(i + 1 + __z, __n + 1);
		}
		return (((__n % 2)? 1: -1) * __FNAMESRC__(tgamma)(__n + 1) * sum);
	} else {
		coef = 1 / __FNAMESRC__(tgamma)(1 - __n);
		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
			sum += coef * pz * __FNAMESRC__(riemannzeta)(i + 1);
			coef *= (i + 1) / (i - __n + 1);
			pz *= -__z;
		}
		return ((__FNAMESRC__(polygamma)(-__n) + MASCHERONI - __FNAMESRC__(log)(__z)) *
				__FNAMESRC__(pow)(__z, -__n - 1) / __FNAMESRC_SCAL__(tgamma)(-__n) +
				__FNAMESRC__(pow)(__z, -__n) * sum);
	}
#	else
	if(__FNAMESRC__(real)(__z) < 0 && __FNAMESRC__(imag)(__z) == 0 &&
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(-__z), 2) < 1) {
		return NAN;
	} else if(__FNAMESRC__(real)(__z) <= 0 && __FNAMESRC__(imag)(__z) == 0 &&
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(-__z), 1) == 0) {
		return INFINITY;
	} else if(__n == -1) {	// This is the zeroth derivative of ln(Gamma).
		return __FNAMESRC__(lgamma)(__z);
	} else if(__n == 0) {
		return __FNAMESRC__(polygamma)(__z);
	} else if(__n > 0) {
		coef = 1 / __FNAMESRC__(pow)(__z, __n + 1);
		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef); i++) {
			sum += coef;
			coef = 1 / __FNAMESRC__(pow)(i + 1 + __z, __n + 1);
		}
		return (((__n % 2)? 1: -1) * __FNAMESRC__(tgamma)(__n + 1) * sum);
	} else {
		coef = 1 / __FNAMESRC__(tgamma)(1 - __n);
		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
			sum += coef * pz * __FNAMESRC__(riemannzeta)(i + 1);
			coef *= (i + 1) / (i - __n + 1);
			pz *= -__z;
		}
		return ((__FNAMESRC__(polygamma)(-__n) + MASCHERONI - __FNAMESRC__(log)(__z)) *
				__FNAMESRC__(pow)(__z, -__n - 1) / __FNAMESRC_SCAL__(tgamma)(-__n) +
				__FNAMESRC__(pow)(__z, -__n) * sum);
	}
#	endif
}

EXTRAMATH_FUNDEF(polygammanu, (__TYPENAME__  __nu, __TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

#	ifndef __IS_COMPLEX__
	// Can't take value when Gamma(z) < 0.
	if(__z < 0 && __FNAMESRC__(fmod)(-__z, 2) < 1) {
		return NAN;
	} else if(__z <= 0 && __FNAMESRC__(fmod)(-__z, 1) == 0) {
		return INFINITY;
	} else if(__nu == -1) {	// This is the zeroth derivative of ln(Gamma).
		return __FNAMESRC__(lgamma)(__z);
	} else if(__nu == 0) {
		return __FNAMESRC__(polygamma)(__z);
	} else if(__nu > 0 && __FNAMESRC__(fmod)(__nu, 1) == 0) {
		coef = 1 / __FNAMESRC__(pow)(__z, __nu + 1);
		for(int i = 0; __FNAMESRC__(fabs)(coef) > __FNAMESRC__(epsilon)(sum); i++) {
			sum += coef;
			coef = 1 / __FNAMESRC__(pow)(i + 1 + __z, __nu + 1);
		}
		return (((__FNAMESRC__(fmod)(__nu, 2) == 1)? 1: -1) * __FNAMESRC__(tgamma)(__nu + 1) * sum);
	} else {
		coef = 1 / __FNAMESRC__(tgamma)(1 - __nu);
		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
			sum += coef * pz * __FNAMESRC__(riemannzeta)(i + 1);
			coef *= (i + 1) / (i - __nu + 1);
			pz *= -__z;
		}
		return ((__FNAMESRC__(polygamma)(-__nu) + MASCHERONI - __FNAMESRC__(log)(__z)) *
				__FNAMESRC__(pow)(__z, -__nu - 1) / __FNAMESRC_SCAL__(tgamma)(-__nu) +
				__FNAMESRC__(pow)(__z, -__nu) * sum);
	}
#	else
	if(__FNAMESRC__(real)(__z) < 0 && __FNAMESRC__(imag)(__z) == 0 &&
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(-__z), 2) < 1) {
		return NAN;
	} else if(__FNAMESRC__(real)(__z) <= 0 && __FNAMESRC__(imag)(__z) == 0 &&
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(-__z), 1) == 0) {
		return INFINITY;
	} else if(__nu == -1) {	// This is the zeroth derivative of ln(Gamma).
		return __FNAMESRC__(lgamma)(__z);
	} else if(__nu == 0) {
		return __FNAMESRC__(polygamma)(__z);
	} else if(__FNAMESRC__(real)(__nu) > 0 && __FNAMESRC__(imag)(__nu) == 0 &&
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) == 0) {
		coef = 1 / __FNAMESRC__(pow)(__z, __nu + 1);
		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef); i++) {
			sum += coef;
			coef = 1 / __FNAMESRC__(pow)(i + 1 + __z, __nu + 1);
		}
		return (((__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 2) == 1)? 1: -1) * __FNAMESRC__(tgamma)(__nu + 1) * sum);
	} else {
		coef = 1 / __FNAMESRC__(tgamma)(1 - __nu);
		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
			sum += coef * pz * __FNAMESRC__(riemannzeta)(i + 1);
			coef *= (i + 1) / (i - __nu + 1);
			pz *= -__z;
		}
		return ((__FNAMESRC__(polygamma)(-__nu) + MASCHERONI - __FNAMESRC__(log)(__z)) *
				__FNAMESRC__(pow)(__z, -__nu - 1) / __FNAMESRC_SCAL__(tgamma)(-__nu) +
				__FNAMESRC__(pow)(__z, -__nu) * sum);
	}
#	endif
}
