/*
 * polynomials.c
 *
 *  Created on: Jan 31, 2023
 *      Author: connor
 */

#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>

EXTRAMATH_FUNDEF(fbinom, (int __n, int __k)) {
	return __FNAMESRC_SCAL__(exp)(__FNAMESRC_SCAL__(lgamma)(__n + 1) -
			__FNAMESRC_SCAL__(lgamma)(__k + 1) - __FNAMESRC_SCAL__(lgamma)(__n - __k + 1));
}

EXTRAMATH_FUNDEF(polynomial, (__TYPENAME__ __x, int __length, const __TYPENAME__ *coefs)) {
	__TYPENAME__ sum = 0, px = 1;

	for(int i = 0; i < __length; i++) {
		sum += px * coefs[i];
		px *= __x;
	}
	return sum;
}

EXTRAMATH_FUNDEF(hermite, (int __n, __TYPENAME__ __x)) {
  if(__n < 0) {
    return NAN;
  }
  __TYPENAME__ sum = 0, px = __FNAMESRC__(pow)(2 * __x, __n), coef = 1;
  
  for(int i = 0; i <= __n / 2; i++) {
    sum += coef * px;
    px /= -__x * __x * 4;
    coef *= (__n - 2 * i) * (__n - 2 * i - 1);
    coef /= i + 1;
  }
  return sum;
}

EXTRAMATH_FUNDEF(hermitederiv, (int __n, __TYPENAME__ __x, int __order)) {
  __TYPENAME__ sum = 0, px = __FNAMESRC__(pow)(2 * __x, __n - __order), coef = 1;
  if(__n < 0 || __order < 0) {
    return NAN;
  }
  for(int i = 0; i < __order; i++) {
    coef *= __n - i;
  }
  
  for(int i = 0; i <= __n / 2; i++) {
    if(__n - 2 * i - __order <= 0) {
      break;
    }
    sum += coef * px;
    px /= -__x * __x * 4;
    coef *= (__n - 2 * i) * (__n - 2 * i - 1) * (__n - 2 * i - __order + 1) * (__n - 2 * i - __order);
    coef /= (i + 1) * (__n - 2 * i) * (__n - 2 * i - 1);
  }
  return sum;
}

#ifndef __IS_COMPLEX__
EXTRAMATH_ARRFUNDEF(hermitecof, (int __n, __TYPENAME__ *__out)) {
  __TYPENAME__ coef = 1;
  if(__n < 0) {
    return 0;
  }
  for(int i = 0; i <= __n / 2; i++) {
    coef *= (__n - 2 * i) * (__n - 2 * i - 1);
    coef /= -4 * (i + 1);
    __out[__n - 2 * i] = coef;
    if(__n - 2 * i - 1 >= 0) {
      __out[__n - 2 * i - 1] = 0;
    }
  }
  return __n + 1;
}
#endif

EXTRAMATH_FUNDEF(laguerre, (int __n, __TYPENAME__ __x)) {
  __TYPENAME__ sum = 0, px = 1, coef = 1;

  if(__n < 0) {
    return NAN;
  }
  
  for(int i = 0; i <= __n; i++) {
    sum += coef * px;
    px *= __x;
    coef /= (i + 1) * (i + 1);
    coef *= (-__n + i);
  }
  return sum;
}

EXTRAMATH_FUNDEF(laguerrederiv, (int __n, __TYPENAME__ __x, int __order)) {
	__TYPENAME__ sum = 0, px = 1, coef = 1;
	if(__n < 0 || __order < 0) {
	  return NAN;
	}

	for(int i = 1; i <= __order; i++) {
		coef *= i;
	}

	for(int i = __order; i <= __n; i++) {
		sum += coef * px;
		px *= __x;
		coef /= (i + 1) * (i - __order + 1);
		coef *= (-__n + i);
	}
	return sum;
}

#ifndef __IS_COMPLEX__
EXTRAMATH_ARRFUNDEF(laguerrecof, (int __n, __TYPENAME__ *__out)) {
  __TYPENAME__ coef = 1;
  if(__n < 0) {
    return 0;
  }
  
  for(int i = 0; i <= __n; i++) {
    __out[i] = coef;
    coef /= (i + 1) * (i + 1);
    coef *= (-__n + i);
  }
  return __n + 1;
}
#endif

EXTRAMATH_FUNDEF(assoclaguerre, (int __n, __TYPENAME__ __alpha, __TYPENAME__ __x)) {
  if(__n < 0) {
    return NAN;
  }
  __TYPENAME__ sum = 0, px = 1, coef = __FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__alpha + __n + 1) -
							 __FNAMESRC_SCAL__(lgamma)(__n + 1) - __FNAMESRC__(lgamma)(__alpha + 1));
  
  for(int i = 0; i <= __n; i++) {
    sum += coef * px;
    px *= __x;
    coef /= (i + 1) * (__alpha + i + 1);
    coef *= (-__n + i);
  }
  return sum;
}

EXTRAMATH_FUNDEF(assoclaguerrederiv, (int __n, __TYPENAME__ __alpha, __TYPENAME__ __x, int __order)) {
  __TYPENAME__ sum = 0, px = 1, coef = __FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__alpha + __n + 1) -
							 __FNAMESRC_SCAL__(lgamma)(__n + 1) - __FNAMESRC__(lgamma)(__alpha + 1));
  if(__n < 0 || __order < 0) {
    return NAN;
  }
  
  for(int i = 1; i <= __order; i++) {
    coef *= i;
  }
  
  for(int i = __order; i <= __n; i++) {
    sum += coef * px;
    px *= __x;
    coef /= (__alpha + i + 1) * (i - __order + 1);
    coef *= (-__n + i);
  }
  return sum;
}

EXTRAMATH_ARRFUNDEF(assoclaguerrecof, (int __n, __TYPENAME__ __alpha, __TYPENAME__ *__out)) {
  __TYPENAME__ coef = __FNAMESRC__(exp)(__FNAMESRC__(lgamma)(__alpha + __n + 1) -
					__FNAMESRC_SCAL__(lgamma)(__n + 1) - __FNAMESRC__(lgamma)(__alpha + 1));
  if(__n < 0) {
    return 0;
  }
  
  for(int i = 0; i <= __n; i++) {
    __out[i] = coef;
    coef /= -(i + 1) * (__alpha + i + 1);
    coef *= (-__n + i);
  }
  return __n + 1;
}

EXTRAMATH_FUNDEF(legendre, (int __n, __TYPENAME__ __x)) {
  __TYPENAME__ sum = 0, px = __FNAMESRC__(pow)(__x, __n), coef = __n * __FNAMESRC__(fbinom)(2 * __n, __n);

  if(__n < 0) {
    return NAN;
  }
  for(int i = 0; i <= __n / 2; i++) {
    sum += coef * px;
    px /= -__x * __x;
    coef *= (__n - i) * (__n - 2 * i) * (__n - 2 * i - 1);
    coef /= (i + 1) * (2 * __n - 2 * i) * (2 * __n - 2 * i - 1);
  }
  return sum;
}

EXTRAMATH_FUNDEF(legendrederiv, (int __n, __TYPENAME__ __x, int __order)) {
  __TYPENAME__ sum = 0, px = __FNAMESRC__(pow)(__x, __n - __order), coef = __n * __FNAMESRC__(fbinom)(2 * __n, __n);

  if(__n < 0 || __order < 0) {
    return NAN;
  }
  for(int i = 0; i < __order; i++) {
    coef *= __n - i;
  }
  
  for(int i = 0; i <= __n / 2; i++) {
    sum += coef * px;
    if(__n - 2 * i - __order <= 0) {
      break;
    }
    px /= -__x * __x;
    coef *= (__n - i) * (__n - 2 * i) * (__n - 2 * i - 1) * (__n - __order - 2 * i) * (__n - __order - 2 * i + 1);
    coef /= (i + 1) * (2 * __n - 2 * i) * (2 * __n - 2 * i - 1) * (__n - 2 * i) * (__n - 2 * i - 1);
  }
  return sum;
}

#ifndef __IS_COMPLEX__
EXTRAMATH_ARRFUNDEF(legendrecof, (int __n, __TYPENAME__ *__out)) {
  __TYPENAME__ coef = __n * __FNAMESRC__(fbinom)(2 * __n, __n);

  if(__n < 0) {
    return 0;
  }
  
  for(int i = 0; i <= __n / 2; i++) {
    __out[__n - 2 * i] = coef;
    if(__n - 2 * i - 1 >= 0) {
      __out[__n - 2 * i - 1] = 0;
    }
    coef *= -(__n - i) * (__n - 2 * i) * (__n - 2 * i - 1);
    coef /= (i + 1) * (2 * __n - 2 * i) * (2 * __n - 2 * i - 1);
  }
  return __n + 1;
}
#endif

EXTRAMATH_FUNDEF(assoclegendre, (int __n, int __m, __TYPENAME__ __x)) {
#	ifndef __IS_COMPLEX__
  // This would be a square root of a negative.
  if(__m % 2 == 1 && (__x >= 1 || __x < -1)) {
    return NAN;
  }
#	endif
  if(__n < 0) {
    return NAN;
  }
  __TYPENAME__ sum = 0, px = 1, coef = 1;
  if(__m < 1) {
    coef /= __FNAMESRC_SCAL__(tgamma)(1 - __m);
  }
  
  for(int i = 0; i <= __n; i++) {
    if(i - __m + 1 > 0) {
      sum += coef * px;
    }
    px *= (1 - __x) / 2;
    if(i - __m + 1 > 0) {
      coef /= i - __m + 1;
    }
    coef *= (__n + 1 + i) * (-__n + i);
    coef /= (i + 1);
  }
  return sum * __FNAMESRC__(pow)((1 + __x) / (1 - __x), __m / 2.0);
}

#ifdef __IS_COMPLEX__
EXTRAMATH_FUNDEF(spherharm,(int __l, int __m, __TYPENAME__ __theta, __TYPENAME__ __phi)) {
  if(__l < 0 || __m > __l || -__m > __l) {
    return NAN;
  }
  return __FNAMESRC_SCAL__(exp)(0.5 * (__FNAMESRC_SCAL__(lgamma)(__l - __m + 1) - __FNAMESRC_SCAL__(lgamma)(__l + __m + 1))) *
    __FNAMESRC_SCAL__(sqrt)((2 * __l + 1) / (4 * M_PI)) * __FNAMESRC__(exp)(I * __m * __phi) *
    __FNAMESRC__(assoclegendre)(__l, __m, __FNAMESRC__(cos)(__theta));
}
#endif

EXTRAMATH_FUNDEF(realspherharm,(int __l, int __m, __TYPENAME__ __theta, __TYPENAME__ __phi)) {
  if(__l < 0 || __m > __l || -__m > __l) {
    return NAN;
  }
  if(__m == 0) {
    return __FNAMESRC_SCAL__(sqrt)((2 * __l + 1) / (4 * M_PI)) *
      __FNAMESRC__(assoclegendre)(__l, __m, __FNAMESRC__(cos)(__theta));
  } else if(__m < 0) {
    if(__m % 2 == 0) {
      return __FNAMESRC_SCAL__(sqrt)((2 * __l + 1) / (2 * M_PI)) *
	__FNAMESRC_SCAL__(exp)(0.5 * (__FNAMESRC_SCAL__(lgamma)(__l - abs(__m)) -
				      __FNAMESRC_SCAL__(lgamma)(__l + abs(__m)))) *
	__FNAMESRC__(assoclegendre)(__l, __m, __FNAMESRC__(cos)(__theta)) *
	__FNAMESRC__(sin)(abs(__m) * __phi);
    } else {
      return -__FNAMESRC_SCAL__(sqrt)((2 * __l + 1) / (2 * M_PI)) *
	__FNAMESRC_SCAL__(exp)(0.5 * (__FNAMESRC_SCAL__(lgamma)(__l - abs(__m)) -
				      __FNAMESRC_SCAL__(lgamma)(__l + abs(__m)))) *
	__FNAMESRC__(assoclegendre)(__l, __m, __FNAMESRC__(cos)(__theta)) *
	__FNAMESRC__(sin)(abs(__m) * __phi);
    }
  } else {
    if(__m % 2 == 0) {
      return __FNAMESRC_SCAL__(sqrt)((2 * __l + 1) / (2 * M_PI)) *
	__FNAMESRC_SCAL__(exp)(0.5 * (__FNAMESRC_SCAL__(lgamma)(__l - abs(__m)) -
				      __FNAMESRC_SCAL__(lgamma)(__l + abs(__m)))) *
	__FNAMESRC__(assoclegendre)(__l, __m, __FNAMESRC__(cos)(__theta)) *
	__FNAMESRC__(sin)(abs(__m) * __phi);
    } else {
      return -__FNAMESRC_SCAL__(sqrt)((2 * __l + 1) / (2 * M_PI)) *
	__FNAMESRC_SCAL__(exp)(0.5 * (__FNAMESRC_SCAL__(lgamma)(__l - abs(__m)) -
				      __FNAMESRC_SCAL__(lgamma)(__l + abs(__m)))) *
	__FNAMESRC__(assoclegendre)(__l, __m, __FNAMESRC__(cos)(__theta)) *
	__FNAMESRC__(cos)(abs(__m) * __phi);
    }
  }
}

EXTRAMATH_FUNDEF(chebyshevt, (int __n, __TYPENAME__ __x)) {
  if(__n < 0) {
    return NAN;
  }
  // These two cases will break the algorithm.
  if(__n == 0) {
    return 1;
  } else if(__n == 1) {
		return __x;
  }
  __TYPENAME__ sum = __FNAMESRC__(pow)(2 * __x, __n) / 2, coef = -__n / 2.0,
    px = __FNAMESRC__(pow)(2 * __x, __n - 2);
  
  for(int i = 1; i <= __n / 2; i++) {
    sum += coef * px;
    px /= 4 * __x * __x;
    coef *= -(__n - 2 * i) * (__n - 2 * i - 1);
    coef /= (i + 1) * (__n - i - 1);
  }
  return sum;
}

EXTRAMATH_FUNDEF(chebyshevtderiv, (int __n, __TYPENAME__ __x, int __order)) {
  if(__n < 0 || __order < 0) {
    return NAN;
  }
  // These cases will break the algorithm.
  if(__n == 0 && __order > 0) {
    return 0;
  } else if(__n == 0) {
    return 1;
  } else if(__n == 1 && __order == 0) {
		return __x;
  } else if(__n == 1 && __order == 1) {
    return 1;
  } else if(__n == 1) {
    return 0;
  }
  __TYPENAME__ sum, coef = -__n / 2.0,
    px = __FNAMESRC__(pow)(2 * __x, __n - 2 - __order);
  
  sum = __FNAMESRC__(pow)(2, __n - 1) * __FNAMESRC__(pow)(__x, __n - __order);
  for(int i = 0; i < __order; i++) {
    sum *= __n - i;
    coef *= i + 1;
  }
  
  for(int i = 1; i <= __n / 2; i++) {
    sum += coef * px;
    if(coef == 0) {
      break;
    }
    px /= __x * __x;
    coef *= -(__n - 2 * i) * (__n - 2 * i - 1) * (__n - __order - 2 * i) * (__n - __order - 2 * i + 1);
    coef /= (i + 1) * (__n - i - 1) * (__n - 2 * i) * (__n - 2 * i - 1) * 4;
  }
  return sum;
}

#ifndef __IS_COMPLEX__
EXTRAMATH_ARRFUNDEF(chebyshevtcof, (int __n, __TYPENAME__ *__out)) {
  if(__n < 0) {
    return 0;
  }
  // These two cases will break the algorithm.
  if(__n == 0) {
    __out[0] = 1;
    return 1;
  } else if(__n == 1) {
    __out[0] = 0;
    __out[1] = 1;
    return 2;
  }
  __TYPENAME__ coef = -__n / 2.0 * __FNAMESRC__(pow)(2, __n - 2);
  __out[0] = 0;
  __out[__n] = __FNAMESRC__(pow)(2, __n - 1);
  
  for(int i = 1; i <= __n / 2; i++) {
    __out[__n - 2 * i] = coef;
    if(__n - 2 * i - 1 >= 0) {
      __out[__n - 2 * i - 1] = 0;
    }
    coef *= -(__n - 2 * i) * (__n - 2 * i - 1);
    coef /= (i + 1) * (__n - i - 1) * 4;
  }
  return __n + 1;
}
#endif

EXTRAMATH_FUNDEF(chebyshevu, (int __n, __TYPENAME__ __x)) {
  if(__n < 0) {
    return NAN;
  }
  // These two cases will break the algorithm.
  if(__n == 0) {
    return 1;
  }
  __TYPENAME__ sum = 0, coef = 1,
    px = __FNAMESRC__(pow)(2 * __x, __n);
  
  for(int i = 0; i <= __n / 2; i++) {
    sum += coef * px;
    px /= 4 * __x * __x;
    coef *= -(__n - 2 * i) * (__n - 2 * i - 1);
    coef /= (i + 1) * (__n - i);
  }
  return sum;
}

EXTRAMATH_FUNDEF(chebyshevuderiv, (int __n, __TYPENAME__ __x, int __order)) {
  if(__n < 0 || __order < 0) {
    return NAN;
  }
  // These two cases will break the algorithm.
  if(__n == 0 && __order == 0) {
    return 1;
  } else if(__n == 0) {
    return 0;
  }
  
  __TYPENAME__ sum = 0, coef = __FNAMESRC__(pow)(2, __n),
    px = __FNAMESRC__(pow)(__x, __n - __order);
  
  for(int i = 0; i < __order; i++) {
    coef *= __n - i;
  }
  
  for(int i = 0; i <= __n / 2; i++) {
    sum += coef * px;
    px /= __x * __x;
    coef *= -(__n - 2 * i) * (__n - 2 * i - 1) * (__n - 2 * i) * (__n - 2 * i + 1);
    coef /= (i + 1) * (__n - i) * (__n - 2 * i - __order + 1) * (__n - 2 * i - __order) * 4;
  }
  return sum;
}

#ifndef __IS_COMPLEX__
EXTRAMATH_ARRFUNDEF(chebyshevucof, (int __n, __TYPENAME__ *__out)) {
  if(__n < 0) {
    return 0;
  }
  // These two cases will break the algorithm.
  if(__n == 0) {
    __out[0] = 1;
    return 1;
  } else if(__n == 1) {
    __out[0] = 0;
    __out[1] = 2;
    return 2;
  }
  __TYPENAME__ coef = 1;
  
  for(int i = 1; i <= __n / 2; i++) {
    __out[__n - 2 * i] = coef;
    if(__n - 2 * i - 1 >= 0) {
      __out[__n - 2 * i - 1] = 0;
    }
    coef *= -(__n - 2 * i) * (__n - 2 * i - 1);
    coef /= (i + 1) * (__n - i) * 4;
  }
  return __n + 1;
}
#endif

EXTRAMATH_FUNDEF(jacobi,(int __n, __TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ __x)) {
  if(__n < 0) {
    return NAN;
  }
  __TYPENAME__ sum = 0, coef = 1 / __FNAMESRC_SCAL__(tgamma)(__n + 1), px = 1;
  
  for(int i = 0; i <= __n; i++) {
    sum += coef * px;
    px *= (1 - __x) / 2;
    coef *= -(__a + __b + __n + 1 + i) * (__n + i) * (__a + 1 + __n);
    coef /= (i + 1) * (__a + i) * (__a + i + 1) * (__a + __n);
	}
  return sum;
}

EXTRAMATH_FUNDEF(jacobideriv,(int __n, __TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ __x, int __order)) {
  if(__n < 0 || __order < 0) {
    return 0;
  }
  __TYPENAME__ sum = 0, coef = 1 / __FNAMESRC_SCAL__(tgamma)(__n + 1), px = 1;
  
  for(int i = 1; i <= __order; i++) {
    coef *= -i / 2.0;
  }
  
  for(int i = __order; i <= __n; i++) {
    sum += coef * px;
    px *= (1 - __x) / 2;
    coef *= -(__a + __b + __n + 1 + i) * (__n + i) * (__a + 1 + __n);
    coef /= (__a + i) * (__a + i + 1) * (__a + __n) * (i - __order + 1);
  }
  return sum;
}

EXTRAMATH_ARRFUNDEF(altjacobicof,(int __n, __TYPENAME__ __a, __TYPENAME__ __b, __TYPENAME__ *__coefs)) {
  if(__n < 0) {
    return 0;
  }
  __TYPENAME__ coef = 1 / __FNAMESRC_SCAL__(tgamma)(__n + 1);
  
  for(int i = 0; i <= __n; i++) {
    __coefs[i] = coef;
    coef *= -(__a + __b + __n + 1 + i) * (__n + i) * (__a + 1 + __n);
    coef /= (i + 1) * (__a + i) * (__a + i + 1) * (__a + __n);
  }
  return __n + 1;
}

