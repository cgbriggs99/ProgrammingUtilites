/*
 * hypergeom.c
 *
 *  Created on: Feb 3, 2023
 *      Author: connor
 */

#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>


EXTRAMATH_FUNDEF(hyper1f1, (__TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ __z)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
		sum += coef * pz;
		coef *= __a + i;
		coef /= (__b + i) * (i + 1);
		pz *= __z;
	}

	return sum;
}


#define __LIMIT_POINTS__ 8
#define __LIMIT_DIST__ 128
#define __PHASE__ 0.1
EXTRAMATH_FUNDEF(hyperu, (__TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ __z)) {
#	ifdef __IS_COMPLEX__
	if(__FNAMESRC__(imag)(__b) == 0 && __FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__b), 1) == 0) {
		__TYPENAME__ sum = 0, root = __FNAMESRC_SCAL__(cos)(2 * M_PI / __LIMIT_POINTS__ + __PHASE__) +
				I * __FNAMESRC_SCAL__(sin)(2 * M_PI / __LIMIT_POINTS__ + __PHASE__),
				direction = __FNAMESRC__(epsilon)(__b);
		for(int i = 0; i < __LIMIT_POINTS__; i++) {
			sum += __FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__b - 1 + direction * __LIMIT_DIST__) -
					__FNAMESRC__(lgamma)(__a)) *
					__FNAMESRC__(pow)(__z, 1 - __b - direction * __LIMIT_DIST__) *
					__FNAMESRC__(hyper1f1)(__a - __b + 1 - direction * __LIMIT_DIST__,
							2 - __b - direction * __LIMIT_DIST__, __z) +
					__FNAMESRC__(exp)(__FNAMESRC__(lgamma)(1 - __b - direction * __LIMIT_DIST__) -
							__FNAMESRC__(lgamma)(__a - __b - direction * __LIMIT_DIST__ + 1)) *
					__FNAMESRC__(hyper1f1)(__a, __b + direction * __LIMIT_DIST__, __z);
			direction *= root;
		}
		return sum / __LIMIT_POINTS__;
	} else {
		return __FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__b - 1) -
				__FNAMESRC__(lgamma)(__a)) *
				__FNAMESRC__(pow)(__z, 1 - __b) *
				__FNAMESRC__(hyper1f1)(__a - __b + 1,
						2 - __b, __z) +
				__FNAMESRC__(exp)(__FNAMESRC__(lgamma)(1 - __b) -
						__FNAMESRC__(lgamma)(__a - __b + 1)) *
				__FNAMESRC__(hyper1f1)(__a, __b, __z);
	}
#	else
	if(__FNAMESRC__(fmod)(__b, 1) == 0) {
		__TYPENAME__ sum = 0,
				step = ((__TYPENAME__) 2 * __LIMIT_DIST__ * __FNAMESRC__(epsilon)(__b)) / (__LIMIT_POINTS__ - 1);

		for(int i = 0; i < __LIMIT_POINTS__; i++) {
			__TYPENAME__ boff = __b + step * (((__TYPENAME__) 2 * i) / (__LIMIT_POINTS__ - 1) - 1);
			sum += __FNAMESRC__(creal)(__FNAMESRC__(cexp)(__FNAMESRC__(clgamma)(boff - 1) -
					__FNAMESRC__(clgamma)(__a))) *
					__FNAMESRC__(pow)(__z, 1 - boff) *
					__FNAMESRC__(hyper1f1)(__a - boff + 1,
							2 - boff, __z) +
					__FNAMESRC__(creal)(__FNAMESRC__(cexp)(__FNAMESRC__(clgamma)(1 - boff) -
							__FNAMESRC__(clgamma)(__a - boff + 1))) *
					__FNAMESRC__(hyper1f1)(__a, boff, __z);
		}
		return sum / __LIMIT_POINTS__;
	} else {
		return __FNAMESRC__(creal)(__FNAMESRC__(cexp)(__FNAMESRC__(clgamma)(__b - 1) -
				__FNAMESRC__(clgamma)(__a))) *
				__FNAMESRC__(pow)(__z, 1 - __b) *
				__FNAMESRC__(hyper1f1)(__a - __b + 1,
						2 - __b, __z) +
				__FNAMESRC__(creal)(__FNAMESRC__(cexp)(__FNAMESRC__(clgamma)(1 - __b) -
						__FNAMESRC__(clgamma)(__a - __b + 1))) *
				__FNAMESRC__(hyper1f1)(__a, __b, __z);
	}
#	endif
}

EXTRAMATH_FUNDEF(whittakerm, (__TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ __z)) {
	return __FNAMESRC__(pow)(__z, __b + 0.5) * __FNAMESRC__(exp)(-__z / 2) *
			__FNAMESRC__(hyper1f1)(__b - __a + 0.5, 2 * __b + 1, __z);
}

EXTRAMATH_FUNDEF(whittakerw, (__TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ __z)) {
	return __FNAMESRC__(pow)(__z, __b + 0.5) * __FNAMESRC__(exp)(-__z / 2) *
			__FNAMESRC__(hyperu)(__b - __a + 0.5, 2 * __b + 1, __z);
}

EXTRAMATH_FUNDEF(hypergeo2f1, (__TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ __c, __TYPENAME__ __z)) {
	if(__FNAMESRC_PREF__(abs)(__z) > 1 && !__FNAMESRC__(isinteger)(__a - __b)) {
		__TYPENAME__ sum1 = 0, sum2 = 0, coef = 1, pz = 1;

		for(int i = 0; !__FNAMESRC__(absconv)(sum1, coef * pz); i++) {
			sum1 += coef * pz;
			pz /= __z;
			coef *= (__a + i) * (__a - __c + i + 1);
			coef /= (i + 1) * (__a - __b + 1 + i);
		}
		pz = 1;
		coef = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum2, coef * pz); i++) {
			sum2 += coef * pz;
			pz /= __z;
			coef *= (__b + i) * (__b - __c + i + 1);
			coef /= (i + 1) * (__b - __a + 1 + i);
		}

#ifdef __IS_COMPLEX__
		return __FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__b - __a) + __FNAMESRC__(lgamma)(__c) -
				__FNAMESRC__(lgamma)(__b) - __FNAMESRC__(lgamma)(__c - __a)) * __FNAMESRC__(pow)(-__z, -__a) * sum1 +
				__FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__a - __b) + __FNAMESRC__(lgamma)(__c) -
						__FNAMESRC__(lgamma)(__a) - __FNAMESRC__(lgamma)(__c - __b)) * __FNAMESRC__(pow)(-__z, -__b) * sum2;
#else
		__TYPENAME__ _Complex out = __FNAMESRC__(cexp)(__FNAMESRC__(clgamma)(__b - __a) + __FNAMESRC__(clgamma)(__c) -
				__FNAMESRC__(clgamma)(__b) - __FNAMESRC__(clgamma)(__c - __a)) * __FNAMESRC__(cpow)(-__z, -__a) * sum1 +
				__FNAMESRC__(cexp)(__FNAMESRC__(clgamma)(__a - __b) + __FNAMESRC__(clgamma)(__c) -
						__FNAMESRC__(clgamma)(__a) - __FNAMESRC__(clgamma)(__c - __b)) *
						__FNAMESRC__(cpow)(-__z, -__b) * sum2;
		if(__FNAMESRC__(cimag)(out) != 0) {
			return NAN;
		} else {
			return __FNAMESRC__(creal)(out);
		}
#endif
	} else if(__FNAMESRC_PREF__(abs)(__z) > 1) {
#ifdef __IS_COMPLEX__
		__TYPENAME__ sum = 0, root = __FNAMESRC_SCAL__(cos)(2 * M_PI / __LIMIT_POINTS__ + __PHASE__) +
				I * __FNAMESRC_SCAL__(sin)(2 * M_PI / __LIMIT_POINTS__ + __PHASE__),
				direction = __FNAMESRC__(epsilon)(__b) * __LIMIT_DIST__;
		for(int i = 0; i < __LIMIT_POINTS__; i++) {
			sum += __FNAMESRC__(hypergeo2f1)(__a, __b + direction, __c, __z);
			direction *= root;
		}
		return sum / __LIMIT_POINTS__;
#else
		__TYPENAME__ sum = 0,
				step = ((__TYPENAME__) 2 * __LIMIT_DIST__ * __FNAMESRC__(epsilon)(__b)) / (__LIMIT_POINTS__ - 1);

		for(int i = 0; i < __LIMIT_POINTS__; i++) {
			__TYPENAME__ boff = __b + step * (((__TYPENAME__) 2 * i) / (__LIMIT_POINTS__ - 1) - 1);
			sum += __FNAMESRC__(hypergeo2f1)(__a, boff, __c, __z);
		}
		return sum / __LIMIT_POINTS__;
#endif
	} else if(__z == 1) {
		__TYPENAME__ sum = 0;
		__SCALARTYPE__ step = ((__SCALARTYPE__) 2 * __LIMIT_DIST__ * __FNAMESRC_SCAL__(epsilon)(1)) / (__LIMIT_POINTS__ - 1);

		for(int i = 0; i < __LIMIT_POINTS__; i++) {
			__SCALARTYPE__ oneoff = 1 + step * (((__TYPENAME__) 2 * i) / (__LIMIT_POINTS__ - 1) - 1);
			sum += __FNAMESRC__(hypergeo2f1)(__a, __b, __c, __z * oneoff);
		}
		return sum / __LIMIT_POINTS__;
	} else if(__FNAMESRC__(isinteger)(__a) && __FNAMESRC__(real)(__a) <= 0 &&
			__FNAMESRC__(isinteger)(__c) && __FNAMESRC__(real)(__c) <= 0 &&
			__FNAMESRC__(real)(__a) <= __FNAMESRC__(real)(__c)) {
		__TYPENAME__ sum = 0, coef = 1, pz = 1;

		for(int i = 0; i <= __FNAMESRC__(real)(__a); i++) {
			sum += coef * pz;
			coef *= (-__a + i) * (__b + i);
			coef /= (-__c + i) * (i + 1);
		}
		return sum;
	} else if(__FNAMESRC__(isinteger)(__b) && __FNAMESRC__(real)(__b) <= 0 &&
			__FNAMESRC__(isinteger)(__c) && __FNAMESRC__(real)(__c) <= 0 &&
			__FNAMESRC__(real)(__b) <= __FNAMESRC__(real)(__c)) {
		__TYPENAME__ sum = 0, coef = 1, pz = 1;

		for(int i = 0; i <= __FNAMESRC__(real)(__b); i++) {
			sum += coef * pz;
			coef *= (-__b + i) * (__a + i);
			coef /= (-__c + i) * (i + 1);
		}
		return sum;
	} else if(__FNAMESRC__(isinteger)(__a) && __FNAMESRC__(real)(__a) <= 0 &&
			__FNAMESRC__(isinteger)(__c) && __FNAMESRC__(real)(__c) <= 0 &&
			__FNAMESRC__(real)(__a) > __FNAMESRC__(real)(__c)) {
		return INFINITY;
	} else if(__FNAMESRC__(isinteger)(__b) && __FNAMESRC__(real)(__b) <= 0 &&
			__FNAMESRC__(isinteger)(__c) && __FNAMESRC__(real)(__c) <= 0 &&
			__FNAMESRC__(real)(__b) > __FNAMESRC__(real)(__c)) {
		return INFINITY;
	} else {
		__TYPENAME__ sum = 0, coef = 1, pz = 1;

		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
			sum += coef * pz;
			coef *= (__a + i) * (__b * i);
			coef /= (__c + i) * (i + 1);
			pz *= __z;
		}
		return sum;
	}
}

static int condition1(unsigned int __p, unsigned int __q, const __TYPENAME__ *__a,
		const __TYPENAME__ *__b, __TYPENAME__ __z) {
	__TYPENAME__ sum = -__a[__p - 1];

	for(int i = 0; i < __p - 1; i++) {
		sum += __b[i] - __a[i];
	}
#ifdef __IS_COMPLEX__
	return (__FNAMESRC__(real)(sum) > 0 && __FNAMESRC_PREF__(abs)(__z) == 1 && __q == __p - 1);
#else
	return (sum > 0 && __FNAMESRC_PREF__(abs)(__z) == 1 && __q == __p - 1);
#endif
}

static int condition2(unsigned int __p, unsigned int __q, const __TYPENAME__ *__a,
		const __TYPENAME__ *__b, __TYPENAME__ __z) {
	long mina = 1, minb = 1;
	for(int i = 0; i < __p; i++) {
		if(__FNAMESRC__(isinteger)(__a[i]) && (long) __a[i] <= 0) {
			if(mina == 1) {
				mina = (long) __a[i];
			} else if((long) __a[i] > mina) {
				mina = (long) __a[i];
			}
		}
	}

	for(int i = 0; i < __q; i++) {
		if(__FNAMESRC__(isinteger)(__b[i]) && (long) __b[i] <= 0) {
			if(minb == 1) {
				minb = (long) __b[i];
			} else if((long) __b[i] > minb) {
				minb = (long) __b[i];
			}
		}
	}

	if(minb == 1 && mina != 1) {
		return 1;
	} else if(minb != 1 && mina != 1) {
		return (minb < mina);
	} else {
		return 0;
	}
}

EXTRAMATH_FUNDEF(hypergeopfq,(unsigned int __p, unsigned int __q, const __TYPENAME__ *__a,
		const __TYPENAME__ *__b, __TYPENAME__ __z)) {

	if((__q >= __p) || (__q == __p - 1 && __FNAMESRC_PREF__(abs)(__z) <= 1) ||
			condition1(__p, __q, __a, __b, __z)) {
		__TYPENAME__ sum = 0, coef = 1, pz = 1;
		for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz); i++) {
			sum += coef * pz;
			pz *= __z;
			for(int j = 0; j < __q; j++) {
				if(j < __p) {
					coef *= __a[j] + i;
				}
				coef /= __b[j] + i;
			}
			for(int j = __q; j < __p; j++) {
				coef *= __a[j] + i;
			}
			coef /= i + 1;
		}
		return sum;
	} else if(condition2(__p, __q, __a, __b, __z)) {
		__TYPENAME__ sum = 0, coef = 1, pz = 1;

		long mina = 1;
		for(int i = 0; i < __p; i++) {
			if(__FNAMESRC__(isinteger)(__a[i]) && (long) __a[i] <= 0) {
				if(mina == 1) {
					mina = (long) __a[i];
				} else if((long) __a[i] > mina) {
					mina = (long) __a[i];
				}
			}
		}

		for(int i = 0; i < -mina; i++) {
			sum += coef * pz;
			pz *= __z;
			for(int j = 0; j < __q; j++) {
				if(j < __p) {
					coef *= __a[j] + i;
				}
				coef /= __b[j] + i;
			}
			for(int j = __q; j < __p; j++) {
				coef *= __a[j] + i;
			}
			coef /= i + 1;
		}
		return sum;
	} else {
		return NAN;
	}
}


EXTRAMATH_FUNDEF(meijerg, (unsigned int __p, unsigned int __q, unsigned int __m, unsigned int __n,
		const __TYPENAME__ *__a, const __TYPENAME__ *__b, __TYPENAME__ __z)) {

	int cond;
	for(int i = 0; i < __m; i++) {
		for(int j = 0; j < i; j++) {
			if(__FNAMESRC__(isinteger)(__b[j] - __b[i])) {
				return NAN;
			}
		}
	}

	if(__p > __q) {
		return NAN;
	}
	if(__p == __q && __FNAMESRC_PREF__(abs)(__z) >= 1) {
		return NAN;
	}

	__TYPENAME__ sum = 0, innersum = 0, coef = 1, innercoef = 1, pz, innerpz;
	int alter = ((-__m - __n + __p) % 2)? -1: 1;

	for(int k = 0; k < __m; k++) {
		coef = 1;
		pz = __FNAMESRC__(pow)(__z, __b[k]);
		__SCALARTYPE__ _Complex cofsum = 0;
		for(int j = 0; j < __m; j++) {
			if(j == k) {
				continue;
			}
			cofsum += __FNAMESRC_SCAL__(clgamma)(__b[j] - __b[k]);
		}
		for(int j = __n; j < __p; j++) {
			cofsum -= __FNAMESRC_SCAL__(clgamma)(__a[j] - __b[k]);
		}
#		ifdef __IS_COMPLEX__
		coef = __FNAMESRC_SCAL__(cexp)(cofsum);
#		else
		coef = __FNAMESRC__(creal)(__FNAMESRC__(cexp)(cofsum));
#		endif

		cofsum = 0;
		for(int j = 0; j < __n; j++) {
			cofsum += __FNAMESRC_SCAL__(clgamma)(1 - __a[j] + __b[k]);
		}
		for(int j = __m; j < __q; j++) {
			cofsum -= __FNAMESRC_SCAL__(clgamma)(1 - __b[j] + __b[k]);
		}
#		ifdef __IS_COMPLEX__
		innercoef = __FNAMESRC_SCAL__(cexp)(cofsum);
#		else
		innercoef = __FNAMESRC__(creal)(__FNAMESRC__(cexp)(cofsum));
#		endif

		innerpz = 1;
		innersum = 0;

		for(int i = 0; !__FNAMESRC__(absconv)(innersum, innerpz * innercoef); i++) {
			innersum += innerpz * innercoef;
			innerpz *= alter * __z;

			for(int j = 0; j < __n; j++) {
				innercoef *= 1 - __a[j] + __b[k] + i;
			}
			for(int j = __m; j < __q; j++) {
				innercoef /= 1 - __b[j] + __b[k] + i;
			}
			for(int j = __n; j < __p; j++) {
				innercoef *= 1 - __a[j] + __b[k] + i;
			}
			for(int j = 0; j < __n; j++) {
				innercoef /= 1 - __b[j] + __b[k] + i;
			}
		}
		sum += innersum * coef * pz;
	}
	return sum;

}

EXTRAMATH_FUNDEF(appellf1, (__TYPENAME__ __a, __TYPENAME__ __b1, __TYPENAME__ __b2, __TYPENAME__ __c,
		__TYPENAME__ __z1, __TYPENAME__ __z2)) {
	if(__FNAMESRC_PREF__(abs)(__z1) >= 1 || __FNAMESRC_PREF__(abs)(__z2) >= 1) {
		return NAN;
	}

	__TYPENAME__ sum = 0, innersum = 0, coef = 1, innercoef = 1, pz1 = 1, pz2 = 1;

	for(int i = 0; !__FNAMESRC__(absconv)(sum, coef * pz1); i++) {
		innercoef = coef * __b2;
		pz2 = 1;
		innersum = innercoef * pz2;
		for(int j = 1; !__FNAMESRC__(absconv)(innersum, innercoef * pz2); j++) {
			innersum += innercoef * pz2;
			pz2 *= __z2;
			innercoef *= (__a + i + j) * (__b2 + j);
			innercoef /= (__c + i + j) * (j + 1);
		}
		sum += innersum * pz1;
		pz1 *= __z1;
		coef *= (__a + i);
		coef /= (__c + i) * (i + 1);
	}
	return sum;
}
