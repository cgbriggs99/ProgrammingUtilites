/*
 * integrated.c
 *
 *  Created on: Feb 3, 2023
 *      Author: connor
 */

#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>

EXTRAMATH_FUNDEF(fresnels,(__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = M_PI / 6, pz = 1;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
		sum += coef * pz;
		pz *= -__z * __z * __z * __z;
		coef *= (4 * i + 3) * M_PI * M_PI;
		coef /= (4 * i + 7) * 4 * (2 * i + 2) * (2 * i + 3);
	}

	return __z * __z * __z * sum;
}

EXTRAMATH_FUNDEF(fresnelc,(__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
		sum += coef * pz;
		pz *= -__z * __z * __z * __z;
		coef *= (4 * i + 1) * M_PI * M_PI;
		coef /= (4 * i + 5) * 4 * (2 * i + 1) * (2 * i + 2);
	}

	return __z * sum;
}

EXTRAMATH_FUNDEF(integrale,(__TYPENAME__ __nu, __TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz / (1 + i - __nu)); i++) {
		if(i != __nu) {
			sum += coef * pz / (1 + i - __nu);
		}
		coef /= -(i + 1);
		pz *= __z;
	}
#	ifdef __IS_COMPLEX__
	if(__FNAMESRC__(imag)(__nu) == 0 && __FNAMESRC__(real)(__nu) >= 1 &&
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) == 0) {
		return __FNAMESRC__(pow)(-__z, __nu - 1) / __FNAMESRC__(tgamma)(__nu) *
				(__FNAMESRC__(dilog)(__nu) - __FNAMESRC__(log)(__z)) - sum;
	}
#	else
	if(__nu == 1) {
		if(__z >= 0) {
			return __FNAMESRC__(pow)(-__z, __nu - 1) / __FNAMESRC__(tgamma)(__nu) *
					(__FNAMESRC__(dilog)(__nu) - __FNAMESRC__(log)(__z)) - sum;
		} else {
			return NAN;
		}
	}
#	endif
	return __FNAMESRC__(tgamma)(1 - __nu) * __FNAMESRC__(pow)(__z, __nu - 1) - sum;
}

EXTRAMATH_FUNDEF(integralei,(__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;
	for(int i = 1; !__FNAMESRC__(absconv)(sum, coef * pz / i); i++) {
		sum += coef * pz / i;
		coef /= i + 1;
		pz *= __z;
	}
#	ifndef __IS_COMPLEX__
	return __FNAMESRC__(log)(__FNAMESRC_PREF__(abs)(__z)) + sum + MASCHERONI;
#	else
	if(__FNAMESRC__(imag)(__z) == 0 && __FNAMESRC__(real)(__z) < 0) {
		return __FNAMESRC_SCAL__(log)(__FNAMESRC_PREF__(abs)(__z)) + sum + MASCHERONI;
	} else {
		return __FNAMESRC__(log)(__z) + sum + MASCHERONI;
	}
#	endif
}

EXTRAMATH_FUNDEF(logint,(__TYPENAME__ __z)) {
	return __FNAMESRC__(integralei)(__FNAMESRC__(log)(__z));
}

EXTRAMATH_FUNDEF(sinint, (__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = __z;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz / (2 * i + 1)); i++) {
		sum += coef * pz / (2 * i + 1);
		coef /= -(2 * i + 2) * (2 * i + 3);
		pz *= __z * __z;
	}
	return sum;
}

EXTRAMATH_FUNDEF(cosint, (__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 0.5, pz = 1;
	for(int i = 1; !__FNAMESRC__(absconv)(sum, coef * pz / i); i++) {
		sum += coef * pz / i;
		coef /= (2 * i + 1) * (2 * i + 2);
		pz *= -__z * __z;
	}
#	ifndef __IS_COMPLEX__
	return __FNAMESRC__(log)(__FNAMESRC_PREF__(abs)(__z)) + sum + MASCHERONI;
#	else
	if(__FNAMESRC__(imag)(__z) == 0 && __FNAMESRC__(real)(__z) < 0) {
		return __FNAMESRC_SCAL__(log)(__FNAMESRC_PREF__(abs)(__z)) + sum + MASCHERONI;
	} else {
		return __FNAMESRC__(log)(__z) + sum + MASCHERONI;
	}
#	endif
}

EXTRAMATH_FUNDEF(sinhint, (__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = __z;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz / (2 * i + 1)); i++) {
		sum += coef * pz / (2 * i + 1);
		coef /= (2 * i + 2) * (2 * i + 3);
		pz *= __z * __z;
	}
	return sum;
}

EXTRAMATH_FUNDEF(coshint, (__TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 0.5, pz = 1;
	for(int i = 1; !__FNAMESRC__(absconv)(sum, coef * pz / i); i++) {
		sum += coef * pz / i;
		coef /= (2 * i + 1) * (2 * i + 2);
		pz *= __z * __z;
	}
#	ifndef __IS_COMPLEX__
	return __FNAMESRC__(log)(__FNAMESRC_PREF__(abs)(__z)) + sum + MASCHERONI;
#	else
	if(__FNAMESRC__(imag)(__z) == 0 && __FNAMESRC__(real)(__z) < 0) {
		return __FNAMESRC_SCAL__(log)(__FNAMESRC_PREF__(abs)(__z)) + sum + MASCHERONI;
	} else {
		return __FNAMESRC__(log)(__z) + sum + MASCHERONI;
	}
#	endif
}
