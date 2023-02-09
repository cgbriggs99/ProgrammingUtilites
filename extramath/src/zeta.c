/*
 * zeta.c
 *
 *  Created on: Feb 4, 2023
 *      Author: connor
 */

#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>

EXTRAMATH_FUNDEF(riemannzeta, (__TYPENAME__ __s)) {
#ifdef __IS_COMPLEX__
	if(__FNAMESRC__(real)(__s) > 1) {
#else
	if(__s > 1) {
#endif
		__TYPENAME__ sum1 = 0, sum2 = 1;

		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			sum2 = sum1;
			sum1 += 1 / __FNAMESRC__(pow)(2 * i + 1, __s);
		}
		return sum1 / (1 - __FNAMESRC__(pow)(2, 1 - __s));
	}
#ifdef __IS_COMPLEX__
	if(__FNAMESRC__(real)(__s) > 0 && __s != 1) {
#else
	if(__s > 0 && __s != 1) {
#endif
		__TYPENAME__ sum1 = 0, sum2 = 1;
		int alt = 1;

		for(int i = 1; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			sum2 = sum1;
			sum1 += 1 / __FNAMESRC__(pow)(i, __s) * alt;
			alt *= -1;
		}
		return sum1 / (1 - __FNAMESRC__(pow)(2, 1 - __s));
	} else if(__s == 1) {
		return INFINITY;
	} else if(__s == 0) {
		return -0.5;
	} else {
		return 2 * __FNAMESRC__(pow)(2 * M_PI, __s - 1) * __FNAMESRC__(tgamma)(1 - __s) *
				__FNAMESRC__(sin)(M_PI_2 * __s) * __FNAMESRC__(riemannzeta)(1 - __s);
	}
}

EXTRAMATH_FUNDEF(riemannzeta2, (__TYPENAME__ __s, __TYPENAME__ __a)) {
	if(__a == 1) {
		return __FNAMESRC__(riemannzeta)(__s);
	} else if(__s == 1) {
		return INFINITY;
	}
#ifdef __IS_COMPLEX__
	if(__FNAMESRC__(real)(__s) > 1 && !__FNAMESRC__(isnatural)(-__a) && __a != 0) {
		__TYPENAME__ sum1 = 0, sum2 = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			sum2 = sum1;
			sum1 += __FNAMESRC__(pow)((i + __a) * (i + __a), -__s / 2);
		}
		return sum1;
	} else if(__FNAMESRC__(real)(__s) > 1 && (__FNAMESRC__(isnatural)(-__a) || __a == 0)) {
		__TYPENAME__ sum1 = 0, sum2 = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			if(i == -__a) {
				continue;
			}
			sum2 = sum1;
			sum1 += __FNAMESRC__(pow)((i + __a) * (i + __a), -__s / 2);
		}
		return sum1;
	} else if(__FNAMESRC__(imag)(__a) == 0 && __FNAMESRC__(real)(__a) > -1) {
		__TYPENAME__ sum1 = 0, sum2 = 1, innersum = 0, coef = 1;
		for(int n = 0; !__FNAMESRC__(absconv)(sum1, sum2 - sum1); n++) {
			sum2 = sum1;
			innersum = 0;
			coef = 1;
			for(int i = 0; i <= n; i++) {
				innersum += coef * __FNAMESRC__(pow)(__a + i, 1 - __s);
				coef *= -(n - i);
				coef /= (i + 1);
			}
			sum1 += innersum / (n + 1);
		}
		return sum1 / (__s - 1);
	} else if(__FNAMESRC__(imag)(__a) == 0) {
		int shift = 1 - (int) __FNAMESRC__(real)(__a);
		__TYPENAME__ sum = 0;
		for(int i = 0; i < shift - 1; i++) {
			sum += 1 / __FNAMESRC__(pow)(__FNAMESRC_PREF__(abs)(__a + i), __s);
		}
		return __FNAMESRC__(riemannzeta2)(__s, __a + shift) + sum;
	} else {
		// TODO: Try to find a definition for __s < 1 and __a not real.
		return NAN;
	}
#else
	if(__s > 1 && !__FNAMESRC__(isnatural)(-__a) && __a != 0) {
		__TYPENAME__ sum1 = 0, sum2 = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			sum2 = sum1;
			sum1 += __FNAMESRC__(pow)(__FNAMESRC_PREF__(abs)(i + __a), -__s / 2);
		}
		return sum1;
	} else if(__s > 1 && (!__FNAMESRC__(isnatural)(-__a) || __a == 0)) {
		__TYPENAME__ sum1 = 0, sum2 = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			sum2 = sum1;
			sum1 += __FNAMESRC__(pow)(__FNAMESRC_PREF__(abs)(i + __a), -__s / 2);
		}
		return sum1;
	} else if(__a > -1) {
		__TYPENAME__ sum1 = 0, sum2 = 1, innersum = 0, coef = 1;
		for(int n = 0; !__FNAMESRC__(absconv)(sum1, sum2 - sum1); n++) {
			sum2 = sum1;
			innersum = 0;
			coef = 1;
			for(int i = 0; i <= n; i++) {
				innersum += coef * __FNAMESRC__(pow)(__a + i, 1 - __s);
				coef *= -(n - i);
				coef /= (i + 1);
			}
			sum1 += innersum / (n + 1);
		}
		return sum1 / (__s - 1);
	} else {
		int shift = 1 - (int) __a;
		__TYPENAME__ sum = 0;
		for(int i = 0; i < shift - 1; i++) {
			sum += 1 / __FNAMESRC__(pow)(__FNAMESRC_PREF__(abs)(__a + i), __s);
		}
		return __FNAMESRC__(riemannzeta2)(__s, __a + shift) + sum;
	}
#endif
}

EXTRAMATH_FUNDEF(lerchtrans, (__TYPENAME__ __z, __TYPENAME__ __s, __TYPENAME__ __a)) {
#ifdef __IS_COMPLEX__
	if(__FNAMESRC_PREF__(abs)(__z) < 1 ||
			(__FNAMESRC_PREF__(abs)(__z) == 1 && __FNAMESRC__(real)(__s) > 1)) {
		__TYPENAME__ sum1 = 0, sum2 = 1, pz = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			if(i != __a) {
				sum2 = sum1;
				sum1 += pz / __FNAMESRC__(pow)((__a + i) * (__a + i), __s / 2);
			}
			pz *= __z;
		}
		return sum1;
	} else if(__z == 1) {
		return __FNAMESRC__(riemannzeta2)(__s, __a);
	} else if(__FNAMESRC_PREF__(abs)(__z) > 1) {
		return INFINITY;
	} else {
		return NAN;
	}
#else
	if(__FNAMESRC_PREF__(abs)(__z) < 1 ||
			(__FNAMESRC_PREF__(abs)(__z) == 1 && __s > 1)) {
		__TYPENAME__ sum1 = 0, sum2 = 1, pz = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum1 - sum2); i++) {
			if(i != __a) {
				sum2 = sum1;
				sum1 += pz / __FNAMESRC__(pow)((__a + i) * (__a + i), __s / 2);
			}
			pz *= __z;
		}
		return sum1;
	} else if(__z == 1) {
		return __FNAMESRC__(riemannzeta2)(__s, __a);
	} else if(__FNAMESRC_PREF__(abs)(__z) > 1) {
		return INFINITY;
	} else {
		return NAN;
	}
#endif
}

EXTRAMATH_FUNDEF(dilog, (__TYPENAME__ __z)) {
	if(__FNAMESRC_PREF__(abs)(__z) < 1) {
		__TYPENAME__ sum = 0, pz = __z;
		for(int i = 1; !__FNAMESRC__(absconv)(sum, pz / (i * i)); i++) {
			sum += pz / (i * i);
			pz *= __z;
		}
		return sum;
#	ifndef __IS_COMPLEX__
	} else if(__z > 1) {	// This gives a complex number.
		return NAN;
#	endif
	} else {
		__TYPENAME__ sum = 0, pz = __z;
		for(int i = 1; !__FNAMESRC__(absconv)(sum, 1 / (pz * i * i)); i++) {
			sum += 1 / (pz * i * i);
			pz *= __z;
		}
		return -0.5 * __FNAMESRC__(log)(-__z) - M_PI * M_PI / 6 - sum;
	}
}

EXTRAMATH_FUNDEF(polylog2, (__TYPENAME__ __z, __TYPENAME__ __nu)) {
	if(__FNAMESRC_PREF__(abs)(__z) == 1) {
		return __FNAMESRC__(tgamma)(1 - __nu) * __FNAMESRC__(pow)(-__FNAMESRC__(log)(__z), __nu - 1) +
				__FNAMESRC__(riemannzeta)(__z);
	} else if(__FNAMESRC_PREF__(abs)(__z) < 1) {
		__TYPENAME__ sum = 0, pz = __z;
		for(int i = 1; !__FNAMESRC__(absconv)(sum, pz / __FNAMESRC__(pow)(i, __nu)); i++) {
			sum += pz / __FNAMESRC__(pow)(i, __nu);
			pz *= __z;
		}
		return sum;
	} else if(__nu == 2) {
		return __FNAMESRC__(dilog)(__z);
#	ifndef __IS_COMPLEX__
	} else if(__nu == 0) {
		return 1 / (1 - __z);
	} else {	// This gives a complex number.
		return NAN;
	}
#	else
	} else if(__FNAMESRC__(real)(__nu) < 0){
		__TYPENAME__ sum1 = 0, sum2 = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum1, sum2 - sum1); i++) {
			sum2 = sum1;
			sum1 += __FNAMESRC__(pow)(i + __FNAMESRC__(log)(-__z) / (2 * M_PI * I) + 0.5, __nu - 1);
		}
		return __FNAMESRC__(pow)(2 * M_PI, __nu) / __FNAMESRC__(tgamma)(__nu) *
				__FNAMESRC__(exp)(M_PI * I * __nu / 2) * sum1 - __FNAMESRC__(exp)(M_PI * I * __nu) *
				__FNAMESRC__(polylog2)(1 / __z, __nu);
	} else if(__nu == 0) {
		return 1 / (1 - __z);
	} else {
		// TODO: formula for when Re(__nu) > 0 and |__z| > 1.
		return NAN;
	}
#	endif
}
