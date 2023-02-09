/*
 * zeta.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

/*
 * Zeta functions.
 */
EXTRAMATH_DECL(riemannzeta, (__TYPE_NAME__ __s));
EXTRAMATH_DECL(riemannzeta2, (__TYPE_NAME__ __s, __TYPE_NAME__ __a));

/*
 * Lerch transcendent.
 */
EXTRAMATH_DECL(lerchtrans, (__TYPE_NAME__ __z, __TYPE_NAME__ __s, __TYPE_NAME__ __a));

/*
 * Polylogarithms.
 */
EXTRAMATH_DECL(dilog, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(polylogn, (__TYPE_NAME__ __n, __TYPE_NAME__ __x));

