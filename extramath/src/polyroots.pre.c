/*
 * polyroots.c
 *
 *  Created on: Jan 31, 2023
 *      Author: connor
 */

#include "${CMAKE_SOURCE_DIR}/extramath/include/extramath.h"
#include "${CMAKE_SOURCE_DIR}/extramath/include/extramath_srcdefs.h"
#include <math.h>
#include <${LAPACKE_HEADER}>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef __IS_COMPLEX__
static __TYPENAME__ evaluate(__TYPENAME__ x, const __TYPENAME__ *coefs, unsigned int __degree) {
	__TYPENAME__ sum = 0, xp = 1;
#else
static __TYPENAME__ _Complex evaluate(__TYPENAME__ _Complex x, const __TYPENAME__ *coefs, unsigned int __degree) {
	__TYPENAME__ _Complex sum = 0, xp = 1;
#endif

	for(int i = 0; i <= __degree; i++) {
		sum += xp * coefs[i];
		xp *= x;
	}
	return sum;
}

#ifdef __IS_COMPLEX__
static __TYPENAME__ derivative(__TYPENAME__ x, const __TYPENAME__ *coefs, unsigned int __degree,
		unsigned int order) {
	__TYPENAME__ sum = 0, xp = 1, prod = 1;
#else
static __TYPENAME__ _Complex derivative(__TYPENAME__ _Complex x, const __TYPENAME__ *coefs, unsigned int __degree,
		unsigned int order) {
	__TYPENAME__ sum = 0, xp = 1;
	__TYPENAME__ prod = 1;
#endif
	// Get the coefficients due to the powers.
	for(int i = 1; i <= order; i++) {
		prod *= i;
	}

	for(int i = order; i <= __degree; i++) {
		sum += prod * xp * coefs[i];
		xp *= x;
		prod /= i - order + 1;
		prod *= i + 1;
	}
	return sum;
}

#define __CYCLES__ 10
#define __SMALL_VAL__ 1e-6
#ifdef __IS_COMPLEX__
static void newtonsteps(const __TYPENAME__ *coefs, unsigned int __order, __TYPENAME__ *roots) {
	__TYPENAME__ *x0 = calloc(__order, sizeof(__TYPENAME__)),
			*x1 = calloc(__order, sizeof(__TYPENAME__)),
			*x2 = calloc(__order, sizeof(__TYPENAME__)),
			yval, deriv;
	size_t elemsize = sizeof(__TYPENAME__);
#else
static void newtonsteps(const __TYPENAME__ *coefs, unsigned int __order, __TYPENAME__ _Complex *roots) {
	__TYPENAME__ _Complex *x0 = calloc(__order, sizeof(__TYPENAME__ _Complex)),
			*x1 = calloc(__order, sizeof(__TYPENAME__ _Complex)),
			*x2 = calloc(__order, sizeof(__TYPENAME__ _Complex)),
			yval, deriv;
	size_t elemsize = sizeof(__TYPENAME__ _Complex);
#endif
	memcpy(x2, roots, (__order) * elemsize);

	for(int cycle = 0; cycle < __CYCLES__; cycle++) {
		memcpy(x0, x1, (__order) * elemsize);
		memcpy(x1, x2, (__order) * elemsize);
		for(int i = 0; i < __order; i++) {
			yval = evaluate(x1[i], coefs, __order);
			if(yval == 0) {	// Already found a root.
				continue;
			}
			deriv = derivative(x1[i], coefs, __order, 1);
			// Zero derivative, try to resolve.
			if(deriv == 0) {
				for(int j = 2; j < __order && (deriv == 0 || yval == 0); j++) {
					yval = deriv;
					deriv = derivative(x1[i], coefs, __order, j);
				}
				// Could not resolve. Taking the same step as last time to try to escape.
				if(deriv == 0 || yval == 0) {
					if(cycle < 2) {
						// This should alternate so that complex numbers go in different directions.
						x2[i] += __SMALL_VAL__ * x2[i] * ((i % 2)? -1: 1);
					} else {
						x2[i] += x1[i] - x0[i];
					}
				} else {	// Could resolve. Take the step.
					x2[i] = x1[i] - yval / deriv;
				}
			} else {	// Non-zero derivative. Step.
				x2[i] = x1[i] - yval / deriv;
			}
		}
	}

	memcpy(roots, x2, (__order) * elemsize);
	free(x0);
	free(x1);
	free(x2);
}

EXTRAMATH_ARRFUNDEF(polyroots,(const __TYPENAME__ *__coefs, unsigned int __order, __TYPENAME__ *__roots)) {
  __TYPENAME__ *companion = calloc((__order) * (__order), sizeof(__TYPENAME__)),
    *real = calloc(__order, sizeof(__TYPENAME__)),
    *imag = calloc(__order, sizeof(__TYPENAME__)),
    *beta = calloc(__order, sizeof(__TYPENAME__)),
    *eye = calloc((__order) * (__order), sizeof(__TYPENAME__)),
    *lv = calloc((__order), sizeof(__TYPENAME__)),
    *rv = calloc((__order), sizeof(__TYPENAME__));
  
	// Set up the companion matrix.
	for(int i = 0; i < __order; i++) {
		for(int j = 0; j < __order; j++) {
			if(i == __order - 1) {
				companion[i + j * (__order)] = -__coefs[j] / __coefs[__order];
			} else if(j == i + 1) {
				companion[i + j * (__order)] = 1;
			} else {
				companion[i + j * (__order)] = 0;
			}
			if(i == j) {
				eye[i + j * (__order)] = 1;
			} else {
				eye[i + j * (__order)] = 0;
			}
		}
	}
#	if __TYPEVAL__ == __TYPEVAL_DOUBLE__
	LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'N', __order, companion, __order, real, imag,  lv, 1, rv, 1);
#	elif __TYPEVAL__ == __TYPEVAL_FLOAT__
	LAPACKE_sgeev(LAPACK_COL_MAJOR, 'N', 'N', __order, companion, __order, real, imag, lv, 1, rv, 1);
#	elif __TYPEVAL__ == __TYPEVAL_COMPLEXDOUBLE__
	LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'N', __order, companion, __order, real, lv, 1, rv, 1);
#	elif __TYPEVAL__ == __TYPEVAL_COMPLEXFLOAT__
	LAPACKE_cgeev(LAPACK_COL_MAJOR, 'N', 'N', __order, companion, __order, real, lv, 1, rv, 1);
#	elif defined(__IS_COMPLEX__)
	__FNAMESRC__(geneigen)('N', 'N', __order, companion, eye, real, beta, lv, rv);
#	else
	__FNAMESRC__(geneigen)('N', 'N', __order, companion, eye, real, imag, beta, lv, rv);
#	endif
	int ret_val = 0;	
#	ifdef __IS_COMPLEX__
#       if !defined(NO_NEWTON) || NO_NEWTON == 0
	newtonsteps(__coefs, __order, real);
#       endif
	memcpy(__roots, real, (__order) * sizeof(__TYPENAME__));
	ret_val = __order;
#	else
	__TYPENAME__ _Complex *eigs = calloc(__order, sizeof(__TYPENAME__ _Complex));
	for(int i = 0; i < __order; i++) {
		eigs[i] = __builtin_complex((__TYPENAME__) real[i], (__TYPENAME__) imag[i]);
	}
#       if !defined(NO_NEWTON) || NO_NEWTON == 0
	newtonsteps(__coefs, __order, eigs);
#       endif
	for(int i = 0; i < __order; i++) {
	  if(__FNAMESRC_PREF__(abs)(__FNAMESRC__(cimag)(eigs[i])) <= 1e-8) {
			__roots[ret_val] = __FNAMESRC__(creal)(eigs[i]);
			ret_val++;
		}
	}
	free(eigs);
#	endif
	free(companion);
	free(real);
	free(imag);
	free(beta);
	free(eye);
	free(lv);
	free(rv);

	return ret_val;

}

