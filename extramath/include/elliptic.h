/*
 * elliptic.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif
/*
 * Complete and incomplete elliptic integrals.
 */
EXTRAMATH_DECL(elliptice, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(elliptick, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(ellipticpi, (__TYPE_NAME__ __n, __TYPE_NAME__ __m));
EXTRAMATH_DECL(incelliptice, (__TYPE_NAME__ __x, __TYPE_NAME__ __m));
EXTRAMATH_DECL(incellipticf, (__TYPE_NAME__ __x, __TYPE_NAME__ __m));
EXTRAMATH_DECL(incellipticpi, (__TYPE_NAME__ __n, __TYPE_NAME__ __x, __TYPE_NAME__ __m));
EXTRAMATH_DECL(jacobiz, (__TYPE_NAME__ __x, __TYPE_NAME__ __m));

