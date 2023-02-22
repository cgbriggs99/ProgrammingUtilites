/*
 * hypergeometric.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

/*
 * Parabolic cyllinder function.
 */
EXTRAMATH_DECL(parabolicd, (__TYPE_NAME__ __nu, __TYPE_NAME__ __z));

/*
 * Hypergeometric functions.
 */
EXTRAMATH_DECL(hyper1f1, (__TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __z));
EXTRAMATH_DECL(hyperu, (__TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __z));
EXTRAMATH_DECL(whittakerm, (__TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __z));
EXTRAMATH_DECL(whittakerw, (__TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __z));
EXTRAMATH_DECL(hyper2f1, (__TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __c, __TYPE_NAME__ __z));
EXTRAMATH_DECL(hyperpfq, (unsigned int __p, unsigned int __q, const __TYPE_NAME__ *__a,
		const __TYPE_NAME__ *__b, __TYPE_NAME__ __z));
EXTRAMATH_DECL(meijerg, (unsigned int __p, unsigned int __q, unsigned int __m, unsigned int __n,
		const __TYPE_NAME__ *__a, const __TYPE_NAME__ *__b, __TYPE_NAME__ __z));

EXTRAMATH_DECL(appellf1, (__TYPE_NAME__ __a, __TYPE_NAME__ __b1, __TYPE_NAME__ __b2, __TYPE_NAME__ __c,
		__TYPE_NAME__ __z1, __TYPE_NAME__ __z2));

EXTRAMATH_DECL(hyper1f1reg, (__TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __z));
EXTRAMATH_DECL(hyper2f1reg, (__TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __c, __TYPE_NAME__ __z));
EXTRAMATH_DECL(hyperpfqreg, (unsigned int __p, unsigned int __q, const __TYPE_NAME__ *__a,
		const __TYPE_NAME__ *__b, __TYPE_NAME__ __z));
