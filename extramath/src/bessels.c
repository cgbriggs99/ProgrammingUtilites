/*
 * bessels.c
 *
 *  Created on: Jan 31, 2023
 *      Author: connor
 */


#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>
#include <stdio.h>

#ifdef __IS_COMPLEX__

EXTRAMATH_FUNDEF(jn, (int __n, __TYPENAME__ __x)) {
  __TYPENAME__ sum = 0, coef = 1, pz = 1;
  
  if(__x == 0 && __n != 0) {
    return 0;
  } else if(__x == 0 && __n == 0) {
    return 1;
  } else {
    int k = 0;
    while(!__FNAMESRC__(absconv)(sum, coef)) {
      sum += coef * pz;
      coef *= -1.0 / (4.0 * (__n + 1.0 + k) * (k + 1.0));
      pz *= __x * __x;
      k++;
    }
    return __FNAMESRC__(exp)(__n * (__FNAMESRC__(log)(__x) - M_LN2) - __FNAMESRC_SCAL__(lgamma)(__n + 1)) * sum;
  }
}

EXTRAMATH_FUNDEF(j0, (__TYPENAME__ __x)) {
  if(__x == 0) {
    return 1;
  } else {
    __TYPENAME__ sum = 0, coef = 1, pz = 1;
    int k = 0;
    while(!__FNAMESRC__(absconv)(sum, coef)) {
      sum += coef * pz;
      coef *= 1.0 / -(4 * (k + 1) * (k + 1));
      pz *= __x * __x;
      k++;
    }
    return sum;
  }
    
}

EXTRAMATH_FUNDEF(j1, (__TYPENAME__ __x)) {
	return __FNAMESRC__(jn)(1, __x);
}


#define LIMIT_POINTS 4
#define WIDTH 0.01
EXTRAMATH_FUNDEF(yn, (int __n, __TYPENAME__ __x)) {
  if(__x == 0) {
    return NAN;
  }

  // Chebychev nodes.
  __TYPENAME__ nodes[LIMIT_POINTS] = {0.9238795325112867 * WIDTH,
    0.38268343236508984 * WIDTH,
    -0.3826834323650897 * WIDTH,
    -0.9238795325112867 * WIDTH};
  __TYPENAME__ funcs[LIMIT_POINTS];

  // Calculate the interpolant.
  __TYPENAME__ sum = 0;
  for(int i = 0; i < LIMIT_POINTS; i++) {
    __TYPENAME__ prod = __FNAMESRC__(ynu)(__n + nodes[i], __x);
    for(int j = 0; j < LIMIT_POINTS; j++) {
      if(i != j) {
	prod *= (__n - nodes[j]) / (nodes[i] - nodes[j]);
      }
    }
    sum += prod;
  }
  return sum;
}

EXTRAMATH_FUNDEF(y1, (__TYPENAME__ __x)) {
	return __FNAMESRC__(yn)(1, __x);
}

EXTRAMATH_FUNDEF(y0, (__TYPENAME__ __x)) {
	return __FNAMESRC__(yn)(0, __x);
}
#endif

EXTRAMATH_FUNDEF(jnu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
  __TYPENAME__ sum = 0, coef = 1, pz = 1;

#	ifndef __IS_COMPLEX__
  if(__x == 0 && (__nu > 0 || __FNAMESRC__(fmod)(__nu, 1) == 0)) {
    return 0;
  } else if(__x == 0 && __nu < 0 && __FNAMESRC__(fmod)(__nu, 1) != 0) {
    return INFINITY;
  } else {
    // Skips the complex stuff.
#	else
    if(__x == 0 && (__FNAMESRC__(real)(__nu) > 0 || (__FNAMESRC__(imag)(__nu) == 0 &&
						     __FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) == 0))) {
      return 0;
    } else if(__x == 0 && __FNAMESRC__(real)(__nu) < 0 && (__FNAMESRC__(imag)(__nu) != 0 ||
							   __FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) != 0)) {
      return INFINITY;
    } else if(__x == 0 && __FNAMESRC__(real)(__nu) == 0 && __nu != 0) {
      return NAN;
    } else {
#	endif
      int k = 0;
      while(!__FNAMESRC__(absconv)(sum, coef)) {
	sum += coef * pz;
	coef *= -1.0 / (4 * (__nu + 1 + k) * (k + 1));
	pz *= __x * __x;
	k++;
      }
      return __FNAMESRC__(exp)(__nu * (__FNAMESRC__(log)(__x) - M_LN2) - __FNAMESRC__(lgamma)(__nu + 1)) * sum;
    }
}

EXTRAMATH_FUNDEF(ynu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
  __TYPENAME__ sum = 0, coef = 1, pz = 1;
  
#	ifndef __IS_COMPLEX__
  // Avoid calculating twice when __nu is a half-integer.
  if(__FNAMESRC__(fmod)(__FNAMESRC_PREF__(abs)(__nu), 1) == 0.5) {
    if(__FNAMESRC__(fmod)(2 * __nu, 4) == 1 ||
       __FNAMESRC__(fmod)(2 * __nu, 4) == -3) {
      return -__FNAMESRC__(jnu)(-__nu, __x);
    } else {
      return __FNAMESRC__(jnu)(-__nu, __x);
    }
  } else if(__FNAMESRC__(fmod)(__nu, 1) != 0) {
    return (__FNAMESRC__(jnu)(__nu, __x) / __FNAMESRC__(tan)(M_PI * __nu) -
	    __FNAMESRC__(jnu)(-__nu, __x) / __FNAMESRC__(sin)(M_PI * __nu));
  } else {
    return __FNAMESRC__(yn)((int) __nu, __x);
  }
#	else
  // Avoid calculating twice when __nu is a half-integer.
  if(__FNAMESRC__(imag)(__nu) == 0 &&
     __FNAMESRC_SCAL__(fmod)(__FNAMESRC_SCAL_PREF__(abs)(__FNAMESRC__(real)(__nu)), 1) == 0.5) {
    return -__FNAMESRC__(jnu)(-__nu, __x) * __FNAMESRC__(sin)(__nu * M_PI);
  } else if(__FNAMESRC__(imag)(__nu) != 0 || __FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) != 0) {
    return 1 /__FNAMESRC__(sin)(M_PI * __nu) * (__FNAMESRC__(cos)(M_PI * __nu) * __FNAMESRC__(jnu)(__nu, __x) -
						__FNAMESRC__(jnu)(-__nu, __x));
  } else {
    return __FNAMESRC__(yn)((int) __nu, __x);
  }
#	endif
}

EXTRAMATH_FUNDEF(inu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

#	ifndef __IS_COMPLEX__
	if(__x == 0 && (__nu > 0 || __FNAMESRC__(fmod)(__nu, 1) == 0)) {
		return 0;
	} else if(__x == 0 && __nu < 0 && __FNAMESRC__(fmod)(__nu, 1) != 0) {
		return INFINITY;
	} else {
#	else
	if(__x == 0 && (__FNAMESRC__(real)(__nu) > 0 || (__FNAMESRC__(imag)(__nu) == 0 &&
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) == 0))) {
		return 0;
	} else if(__x == 0 && __FNAMESRC__(real)(__nu) < 0 && (__FNAMESRC__(imag)(__nu) != 0 ||
			__FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) != 0)) {
		return INFINITY;
	} else if(__x == 0 && __FNAMESRC__(real)(__nu) == 0 && __nu != 0) {
		return NAN;
	} else {
#	endif
		int k = 0;
		while(!__FNAMESRC__(absconv)(sum, coef)) {
			sum += coef * pz;
			coef *= 1.0 / (4 * (__nu + 1 + k) * (k + 1));
			pz *= __x * __x;
			k++;
		}
		return __FNAMESRC__(exp)(__nu * (__FNAMESRC__(log)(__x) - M_LN2) - __FNAMESRC__(lgamma)(__nu + 1)) * sum;
	}
}

EXTRAMATH_FUNDEF(i0, (__TYPENAME__ __x)) {
  if(__x == 0) {
    return 0;
  }
  __TYPENAME__ sum = 0, coef = 1, pz = 1;
  int k = 0;
  while(!__FNAMESRC__(absconv)(sum, coef)) {
    sum += coef * pz;
    coef *= 1.0 / (4 * (k + 1) * (k + 1));
    pz *= __x * __x;
    k++;
  }
  return sum;
}

EXTRAMATH_FUNDEF(i1, (__TYPENAME__ __x)) {
  if(__x == 0) {
    return 0;
  }
  __TYPENAME__ sum = 0, coef = 1, pz = 1;
  int k = 0;
  while(!__FNAMESRC__(absconv)(sum, coef)) {
    sum += coef * pz;
    coef *= 1.0 / (4 * (k + 2) * (k + 1));
    pz *= __x * __x;
    k++;
  }
  return __x * sum / 2;
}

EXTRAMATH_FUNDEF(in, (int __n, __TYPENAME__ __x)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

	if(__x == 0) {
		return 0;
	} else {
		int k = 0;
		while(!__FNAMESRC__(absconv)(sum, coef)) {
			sum += coef * pz;
			coef *= 1 / (4 * (abs(__n) + 1 + k) * (k + 1));
			pz *= __x * __x;
			k++;
		}
		return __FNAMESRC__(exp)(__n * (__FNAMESRC__(log)(__x) - M_LN2) - __FNAMESRC__(lgamma)(__n + 1)) * sum;
	}
}

EXTRAMATH_FUNDEF(kn, (int __n, __TYPENAME__ __x)) {
#	ifndef __IS_COMPLEX__
	__TYPENAME__ _Complex ret, is;

	switch(__n % 4) {
	case 0:
		is = 1;
		break;
	case 1:
		is = I;
		break;
	case 2:
		is = -1;
		break;
	case 3:
		is = -I;
		break;
	}
	ret = -M_PI_2 * is * __FNAMESRC__(cyn)(__n, I * __x) + ((__n % 2)? 1: -1) *
			(__FNAMESRC__(clog)(I * __x) - __FNAMESRC__(clog)(__x)) * __FNAMESRC__(cin)(__n, __x);
	if(__FNAMESRC__(cimag)(ret) == 0) {
		return __FNAMESRC__(creal)(ret);
	} else {
		return NAN;
	}
#	else
	__TYPENAME__ ret, is;

		switch(__n % 4) {
		case 0:
			is = 1;
			break;
		case 1:
			is = I;
			break;
		case 2:
			is = -1;
			break;
		case 3:
			is = -I;
			break;
		}
		ret = -M_PI_2 * is * __FNAMESRC__(yn)(__n, I * __x) + ((__n % 2)? 1: -1) *
				(__FNAMESRC__(log)(I * __x) - __FNAMESRC__(log)(__x)) * __FNAMESRC__(in)(__n, __x);
		return ret;
#	endif
}

EXTRAMATH_FUNDEF(knu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
	__TYPENAME__ sum = 0, coef = 1, pz = 1;

#	ifndef __IS_COMPLEX__
	if(__FNAMESRC__(fmod)(__nu, 1) != 0) {
		return M_PI_2 /__FNAMESRC__(sin)(M_PI * __nu) * (__FNAMESRC__(inu)(-__nu, __x) -
				__FNAMESRC__(inu)(__nu, __x));
	} else {
		return __FNAMESRC__(kn)((int) __nu, __x);
	}
#	else
	if(__FNAMESRC__(imag)(__nu) != 0 || __FNAMESRC_SCAL__(fmod)(__FNAMESRC__(real)(__nu), 1) != 0) {
		return M_PI_2 /__FNAMESRC__(sin)(M_PI * __nu) * (__FNAMESRC__(inu)(-__nu, __x) -
						__FNAMESRC__(inu)(__nu, __x));
	} else {
		return __FNAMESRC__(kn)((int) __nu, __x);
	}
#	endif
}

EXTRAMATH_FUNDEF(k0, (__TYPENAME__ __x)) {
	return __FNAMESRC__(kn)(0, __x);
}

EXTRAMATH_FUNDEF(k1, (__TYPENAME__ __x)) {
	return __FNAMESRC__(kn)(1, __x);
}

EXTRAMATH_FUNDEF(ai,(__TYPENAME__ __x)) {
	__TYPENAME__ sum1 = 0, sum2 = 0, px = 1, coef1 = 1, coef2 = 1;
	int k = 0;
	while(!__FNAMESRC__(absconv)(sum1, coef1) &&
			!__FNAMESRC__(absconv)(sum2, coef2)) {
		sum1 += px * coef1;
		sum2 += px * coef2;
		px *= __x * __x * __x / 9;
		coef1 /= ((__TYPENAME__) 2 / 3 + k) * (k + 1);
		coef2 /= ((__TYPENAME__) 4 / 3 + k) * (k + 1);
		k += 1;
	}

	return sum1 / (__FNAMESRC_SCAL__(cbrt)(9) * __FNAMESRC_SCAL__(tgamma)((__SCALARTYPE__) 2 / 3)) -
			__x * sum2 / (__FNAMESRC_SCAL__(cbrt)(3) * __FNAMESRC_SCAL__(tgamma)((__SCALARTYPE__) 1 / 3));
}

EXTRAMATH_FUNDEF(bi,(__TYPENAME__ __x)) {
	__TYPENAME__ sum1 = 0, sum2 = 0, px = 1, coef1 = 1, coef2 = 1;
	int k = 0;
	while(!__FNAMESRC__(absconv)(sum1, coef1) &&
			!__FNAMESRC__(absconv)(sum2, coef2)) {
		sum1 += px * coef1;
		sum2 += px * coef2;
		px *= __x * __x * __x / 9;
		coef1 /= ((__TYPENAME__) 2 / 3 + k) * (k + 1);
		coef2 /= ((__TYPENAME__) 4 / 3 + k) * (k + 1);
		k += 1;
	}

	return sum1 / (__FNAMESRC__(pow)(3, 1.0/6) * __FNAMESRC__(tgamma)((__TYPENAME__) 2 / 3)) +
			__FNAMESRC__(pow)(3, 1.0/6) * __x * sum2 / (__FNAMESRC__(tgamma)((__TYPENAME__) 1 / 3));
}

EXTRAMATH_FUNDEF(aip,(__TYPENAME__ __x)) {
	__TYPENAME__ sum1 = 0, sum2 = 0, px = 1, coef1 = 1, coef2 = 1;
	int k = 0;
	while(!__FNAMESRC__(absconv)(sum1, coef1) &&
			!__FNAMESRC__(absconv)(sum2, coef2)) {
		sum1 += px * coef1;
		sum2 += px * coef2;
		px *= __x * __x * __x / 9;
		coef1 /= ((__TYPENAME__) 5 / 3 + k) * (k + 1);
		coef2 /= ((__TYPENAME__) 1 / 3 + k) * (k + 1);
		k += 1;
	}

	return __x * __x * sum1 / (2 * __FNAMESRC_SCAL__(cbrt)(9) * __FNAMESRC_SCAL__(tgamma)((__TYPENAME__) 2 / 3)) -
			sum2 / (__FNAMESRC_SCAL__(cbrt)(3) * __FNAMESRC_SCAL__(tgamma)((__TYPENAME__) 1 / 3));
}

EXTRAMATH_FUNDEF(bip,(__TYPENAME__ __x)) {
	__TYPENAME__ sum1 = 0, sum2 = 0, px = 1, coef1 = 1, coef2 = 1;
	int k = 0;
	while(!__FNAMESRC__(absconv)(sum1, coef1) &&
			!__FNAMESRC__(absconv)(sum2, coef2)) {
		sum1 += px * coef1;
		sum2 += px * coef2;
		px *= __x * __x * __x / 9;
		coef1 /= ((__TYPENAME__) 5 / 3 + k) * (k + 1);
		coef2 /= ((__TYPENAME__) 1 / 3 + k) * (k + 1);
		k += 1;
	}

	return __x * __x * sum1 / (2 * __FNAMESRC__(pow)(3, 1.0/6) * __FNAMESRC__(tgamma)((__TYPENAME__) 2 / 3)) +
			__FNAMESRC__(pow)(3, 1.0/6) * sum2 / (__FNAMESRC__(tgamma)((__TYPENAME__) 1 / 3));
}

EXTRAMATH_FUNDEF(struvehnu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
	__TYPENAME__ sum = 0, px = 1, coef = 1;
	int k = 0;

	while(!__FNAMESRC__(absconv)(sum, coef)) {
		sum += px * coef;
		px *= -__x * __x;
		coef /= ((__TYPENAME__) 3 / 2 + k) * ((__TYPENAME__) 3 / 2 + k + __nu) * 4;
		k++;
	}
	return __FNAMESRC__(pow)(__x / 2, __nu + 1) * M_2_SQRTPI / __FNAMESRC__(tgamma)((__TYPENAME__) 3 / 2 + __nu) * sum;
}

EXTRAMATH_FUNDEF(struvehn, (int __n, __TYPENAME__ __x)) {
	__TYPENAME__ sum = 0, px = 1, coef = 1;
	int k = 0;

	while(!__FNAMESRC__(absconv)(sum, coef)) {
		sum += px * coef;
		px *= -__x * __x;
		coef /= ((__TYPENAME__) 3 / 2 + k) * ((__TYPENAME__) 3 / 2 + k + __n) * 4;
		k++;
	}
	return __FNAMESRC__(pow)(__x / 2, __n + 1) * M_2_SQRTPI / __FNAMESRC__(tgamma)((__TYPENAME__) 3 / 2 + __n) * sum;
}

EXTRAMATH_FUNDEF(struveh0, (__TYPENAME__ __x)) {
	return __FNAMESRC__(struvehn)(0, __x);
}

EXTRAMATH_FUNDEF(struveh1, (__TYPENAME__ __x)) {
	return __FNAMESRC__(struvehn)(1, __x);
}

EXTRAMATH_FUNDEF(struvelnu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
	__TYPENAME__ sum = 0, px = 1, coef = 1;
	int k = 0;

	while(!__FNAMESRC__(absconv)(sum, coef)) {
		sum += px * coef;
		px *= __x * __x;
		coef /= ((__TYPENAME__) 3 / 2 + k) * ((__TYPENAME__) 3 / 2 + k + __nu) * 4;
		k++;
	}
	return __FNAMESRC__(pow)(__x / 2, __nu + 1) * M_2_SQRTPI / __FNAMESRC__(tgamma)((__TYPENAME__) 3 / 2 + __nu) * sum;
}

EXTRAMATH_FUNDEF(struveln, (int __n, __TYPENAME__ __x)) {
	__TYPENAME__ sum = 0, px = 1, coef = 1;
	int k = 0;

	while(!__FNAMESRC__(absconv)(sum, coef)) {
		sum += px * coef;
		px *= __x * __x;
		coef /= ((__TYPENAME__) 3 / 2 + k) * ((__TYPENAME__) 3 / 2 + k + __n) * 4;
		k++;
	}
	return __FNAMESRC__(pow)(__x / 2, __n + 1) * M_2_SQRTPI / __FNAMESRC__(tgamma)((__TYPENAME__) 3 / 2 + __n) * sum;
}

EXTRAMATH_FUNDEF(struvel0, (__TYPENAME__ __x)) {
	return __FNAMESRC__(struveln)(0, __x);
}

EXTRAMATH_FUNDEF(struvel1, (__TYPENAME__ __x)) {
	return __FNAMESRC__(struveln)(1, __x);
}

EXTRAMATH_FUNDEF(spj0, (__TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(jnu)(0.5, __x);
}

EXTRAMATH_FUNDEF(spj1, (__TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(jnu)(1.5, __x);
}

EXTRAMATH_FUNDEF(spjn, (int __n, __TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(jnu)(__n + 0.5, __x);
}

EXTRAMATH_FUNDEF(spjnu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(jnu)(__nu + 0.5, __x);
}

EXTRAMATH_FUNDEF(spy0, (__TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(ynu)(0.5, __x);
}

EXTRAMATH_FUNDEF(spy1, (__TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(ynu)(1.5, __x);
}

EXTRAMATH_FUNDEF(spyn, (int __n, __TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(ynu)(__n + 0.5, __x);
}

EXTRAMATH_FUNDEF(spynu, (__TYPENAME__ __nu, __TYPENAME__ __x)) {
	return __FNAMESRC__(sqrt)(M_PI / 2 / __x) * __FNAMESRC__(ynu)(__nu + 0.5, __x);
}

