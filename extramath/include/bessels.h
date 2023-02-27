/*
 * bessels.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

/*
 * Extensions of jn and yn to include non-integer subscripts.
 */
#if defined(__IS_COMPLEX_HEADER__)
EXTRAMATH_DECL(j0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(j1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(y0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(y1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(jn, (int __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(yn, (int __nu, __TYPE_NAME__ __x));
#endif
EXTRAMATH_DECL(jnu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(ynu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));

/*
 * Modified bessel functions.
 */
EXTRAMATH_DECL(i0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(i1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(in, (int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(k0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(k1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(kn, (int __n, __TYPE_NAME__ __x));

EXTRAMATH_DECL(inu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(knu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));


/*
 * Airy functions.
 */

EXTRAMATH_DECL(ai, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(bi, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(aip, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(bip, (__TYPE_NAME__ __x));

/*
 * Struve functions.
 */
EXTRAMATH_DECL(struveh0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(struveh1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(struvehn, (int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(struvehnu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(struvel0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(struvel1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(struveln, (int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(struvelnu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));

/*
 * Spherical bessel functions.
 */
EXTRAMATH_DECL(spj0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(spj1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(spjn, (int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(spy0, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(spy1, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(spyn, (int __n, __TYPE_NAME__ __x));

EXTRAMATH_DECL(spjnu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(spynu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));
