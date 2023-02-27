/*
 * longlinalg.h
 *
 *  Created on: Jan 31, 2023
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

/*
 * Basic arithmetic
 */
EXTRAMATH_ARRDECL(matmult,(const __TYPE_NAME__ *left, const __TYPE_NAME__ *right, __TYPE_NAME__ *out,
		int lrows, int lc_rr, int rcols));
EXTRAMATH_ARRDECL(matsum,(__TYPE_NAME__ *left_inout, const __TYPE_NAME__ *right, int rows, int cols));
EXTRAMATH_ARRDECL(matdiff,(__TYPE_NAME__ *left_inout, const __TYPE_NAME__ *right, int rows, int cols));
EXTRAMATH_ARRDECL(scalmult,(__TYPE_NAME__ *inout, __TYPE_NAME__ __scal, int rows, int cols));
EXTRAMATH_DECL(dotprod,(const __TYPE_NAME__ *left, const __TYPE_NAME__ *right, int size));
EXTRAMATH_DECL(norm,(const __TYPE_NAME__ *arr, int size));
EXTRAMATH_ARRDECL(transpose, (__TYPE_NAME__ *inout, int rows, int cols));
EXTRAMATH_ARRDECL(conjtranspose, (__TYPE_NAME__ *inout, int rows, int cols));

EXTRAMATH_ARRDECL(matsumtimes,(__TYPE_NAME__ *left_inout, const __TYPE_NAME__ *right, __TYPE_NAME__ scal, int rows, int cols));
EXTRAMATH_ARRDECL(matdifftimes,(__TYPE_NAME__ *left_inout, const __TYPE_NAME__ *right, __TYPE_NAME__ scal, int rows, int cols));

/*
 * Norms and conditions
 */

/*
 * More complicated linear algebra stuff.
 */
EXTRAMATH_DECL(det,(const __TYPE_NAME__ *in, int dim));
EXTRAMATH_DECL(trace,(const __TYPE_NAME__ *in, int dim));
EXTRAMATH_ARRDECL(linearsolve,(__TYPE_NAME__ *coef_inout, __TYPE_NAME__ *right_inout, int eqs, int vars, int rightcols));
EXTRAMATH_ARRDECL(matinv,(__TYPE_NAME__ *inout, int dim));

/*
 * Decomposition schemes.
 */
EXTRAMATH_ARRDECL(qrdecompose,(__TYPE_NAME__ *inout,  __TYPE_NAME__ *r_out, int rows, int cols));
EXTRAMATH_ARRDECL(ludecompose,(__TYPE_NAME__ *l_inout, __TYPE_NAME__ *u_out, int dim));




/*
 * Computes similar to *ggev.
 */
#ifdef __IS_COMPLEX_HEADER__
EXTRAMATH_ARRDECL(geneigen,(unsigned char __do_left, unsigned char __do_right,
		unsigned int __order, __TYPE_NAME__ *__A, __TYPE_NAME__ *__B,
		__TYPE_NAME__ *__eigs, __TYPE_NAME__ *__beta,
		__TYPE_NAME__ *__left, __TYPE_NAME__ *__right));
#else
EXTRAMATH_ARRDECL(geneigen,(unsigned char __do_left, unsigned char __do_right,
		unsigned int __order, __TYPE_NAME__ *__A, __TYPE_NAME__ *__B,
		__TYPE_NAME__ *__real_eigs, __TYPE_NAME__ *__imag_eigs, __TYPE_NAME__ *__beta,
		__TYPE_NAME__ *__left, __TYPE_NAME__ *__right));
#endif

