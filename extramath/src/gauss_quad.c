/*
 * gauss_quad.c
 *
 *  Created on: Feb 3, 2023
 *      Author: connor
 */


#include "../include/extramath.h"
#include "../include/extramath_srcdefs.h"
#include <math.h>

static inline __SCALARTYPE__ square(__SCALARTYPE__ x) {
	return x * x;
}

EXTRAMATH_FUNDEF(gausslegendreint,(__FNAMESRC__(kernel_) __func,
		__SCALARTYPE__ __start, __SCALARTYPE__ __end, unsigned int __points, const void *__extra)) {

	__TYPENAME__ sum = 0;
	__SCALARTYPE__ scale = __end - __start;

	__SCALARTYPE__ *polycoefs = calloc(__points + 1, sizeof(__SCALARTYPE__)),
			*roots = calloc(__points, sizeof(__SCALARTYPE__)),
			weight;

	__FNAMESRC_SCAL__(legendrecof)(__points, polycoefs);
	__FNAMESRC_SCAL__(polyroots)(polycoefs, __points, roots);

	for(int i = 0; i < __points; i++) {
		weight = 2 / ((1 - roots[i] * roots[i]) * square(__FNAMESRC_SCAL__(legendrederiv)(__points, roots[i], 1)));
		sum += __func((roots[i] + 1) / 2 * scale + __start, __extra) * weight;
	}

	free(polycoefs);
	free(roots);

	return sum;

}

EXTRAMATH_FUNDEF(gaussjacobiint,(__FNAMESRC__(kernel_) __func,
		__SCALARTYPE__ __start, __SCALARTYPE__ __end, __SCALARTYPE__ __alpha, __SCALARTYPE__ __beta,
		unsigned int __points, const void *__extra)) {
	if(__alpha <= -1 || __beta <= -1) {
		return NAN;
	}

	__TYPENAME__ sum = 0;
	__SCALARTYPE__ scale = __end - __start;

	__SCALARTYPE__ *polycoefs = calloc(__points + 1, sizeof(__SCALARTYPE__)),
			*roots = calloc(__points, sizeof(__SCALARTYPE__)),
			weight;

	// The terms here are (1 - x) / 2
	__FNAMESRC_SCAL__(altjacobicof)(__points, __alpha, __beta, polycoefs);
	__FNAMESRC_SCAL__(polyroots)(polycoefs, __points, roots);

	for(int i = 0; i < __points; i++) {
		roots[i] = 1 - 2 * roots[i];	// Bring into canonical form.
		weight = -(2 * __points + __alpha + __beta + 2) / (__points + __alpha + __beta + 1) *
				__FNAMESRC_SCAL__(exp)(__FNAMESRC_SCAL__(lgamma)(__points + __alpha + 1) +
						__FNAMESRC_SCAL__(lgamma)(__points + __beta + 1) -
						__FNAMESRC_SCAL__(lgamma)(__points + __alpha + __beta + 1) -
						__FNAMESRC_SCAL__(lgamma)(__points + 2)) *
						__FNAMESRC_SCAL__(pow)(2, __alpha + __beta)  /
						__FNAMESRC_SCAL__(jacobi)(__points + 1, __alpha, __beta, roots[i]) /
						__FNAMESRC_SCAL__(jacobideriv)(__points, __alpha, __beta, roots[i], 1);
		sum += __func((roots[i] + 1) / 2 * scale + __start, __extra) * weight;
	}

	free(polycoefs);
	free(roots);

	return sum;

}

EXTRAMATH_FUNDEF(gausschebychevint1,(__FNAMESRC__(kernel_) __func,
		__SCALARTYPE__ __start, __SCALARTYPE__ __end, unsigned int __points, const void *__extra)) {

	__TYPENAME__ sum = 0;
	__SCALARTYPE__ scale = __end - __start;

	__SCALARTYPE__ weight, root;

	for(int i = 0; i < __points; i++) {
		root = __FNAMESRC_SCAL__(cos)((2 * i + 1) * M_PI / (2 * __points));
		weight = M_PI / __points;
		sum += __func((root + 1) / 2 * scale + __start, __extra) * weight;
	}

	return sum;

}

EXTRAMATH_FUNDEF(gausschebychevint2,(__FNAMESRC__(kernel_) __func,
		__SCALARTYPE__ __start, __SCALARTYPE__ __end, unsigned int __points, const void *__extra)) {

	__TYPENAME__ sum = 0;
	__SCALARTYPE__ scale = __end - __start;

	__SCALARTYPE__ weight, root;

	for(int i = 0; i < __points; i++) {
		root = __FNAMESRC_SCAL__(cos)((i + 1) * M_PI / (__points + 1));
		weight = M_PI / (__points + 1) * square(__FNAMESRC_SCAL__(sin)((i + 1) * M_PI / (__points + 1)));
		sum += __func((root + 1) / 2 * scale + __start, __extra) * weight;
	}

	return sum;

}

EXTRAMATH_FUNDEF(gausslaguerreint,(__FNAMESRC__(kernel_) __func,
				   __SCALARTYPE__ __start, unsigned int __points, const void *__extra)) {
  
  __TYPENAME__ sum = 0;

  __SCALARTYPE__ *polycoefs = calloc(__points + 1, sizeof(__SCALARTYPE__)),
    *roots = calloc(__points, sizeof(__SCALARTYPE__)),
    weight;
  
  __FNAMESRC_SCAL__(laguerrecof)(__points, polycoefs);
  __FNAMESRC_SCAL__(polyroots)(polycoefs, __points, roots);
  
  for(int i = 0; i < __points; i++) {
    weight = roots[i] / square(__FNAMESRC_SCAL__(laguerre)(__points + 1, roots[i]) * (__points + 1));
    sum += __func(roots[i] + __start, __extra) * weight;
  }
  
  free(polycoefs);
  free(roots);
  
  return sum;

}

EXTRAMATH_FUNDEF(gaussgenlaguerreint,(__FNAMESRC__(kernel_) __func,
		__SCALARTYPE__ __start, __SCALARTYPE__ __alpha, unsigned int __points, const void *__extra)) {
	if(__alpha <= -1) {
		return NAN;
	}

	__TYPENAME__ sum = 0;

	__SCALARTYPE__ *polycoefs = calloc(__points + 1, sizeof(__SCALARTYPE__)),
			*roots = calloc(__points, sizeof(__SCALARTYPE__)),
			weight;

	__FNAMESRC_SCAL__(assoclaguerrecof)(__points, __alpha, polycoefs);
	__FNAMESRC_SCAL__(polyroots)(polycoefs, __points, roots);

	for(int i = 0; i < __points; i++) {
		weight = roots[i] / square(__FNAMESRC_SCAL__(laguerre)(__points + 1, roots[i]) * (__points + 1)) *
				__FNAMESRC_SCAL__(exp)(__FNAMESRC_SCAL__(lgamma)(__points + __alpha + 1) -
						__FNAMESRC_SCAL__(lgamma)(__points + 1));
		sum += __func(roots[i] + __start, __extra) * weight;
	}

	free(polycoefs);
	free(roots);

	return sum;

}

EXTRAMATH_FUNDEF(gausshermiteint,(__FNAMESRC__(kernel_) __func,
		unsigned int __points, const void *__extra)) {

	__TYPENAME__ sum = 0;

	__SCALARTYPE__ *polycoefs = calloc(__points + 1, sizeof(__SCALARTYPE__)),
			*roots = calloc(__points, sizeof(__SCALARTYPE__)),
			weight;

	__FNAMESRC_SCAL__(hermitecof)(__points, polycoefs);
	__FNAMESRC_SCAL__(polyroots)(polycoefs, __points, roots);

	for(int i = 0; i < __points; i++) {
		weight = __FNAMESRC_SCAL__(exp)(__FNAMESRC_SCAL__(lgamma)(__points + 1) + (__points) * M_LN2)
				/ (M_2_SQRTPI * square(__points * hermite(__points - 1, roots[i])));
		sum += __func(roots[i], __extra) * weight;
	}

	free(polycoefs);
	free(roots);

	return sum;

}

