/*
 * polynomials.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

/*
 * Evaluate a polynomial given coefficients.
 */
EXTRAMATH_DECL(polynomial, (__TYPE_NAME__ __x, int __len, const __TYPE_NAME__ *coefs));

/*
 * Hermite polynomials. There is one more coefficient than the order, so __out should have length __n + 1.
 * Uses the physicist's definition.
 */
EXTRAMATH_DECL(hermite, (int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(hermitederiv, (int __n, __TYPE_NAME__ __x, int __order));
#ifndef __IS_COMPLEX_HEADER__
EXTRAMATH_ARRDECL(hermitecof, (int __n, __TYPE_NAME__ *__out));
#endif

/*
 * Laguerre polynomials and associated laguerre polynomials.
 */
EXTRAMATH_DECL(laguerre, (int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(laguerrederiv, (int __n, __TYPE_NAME__ __x, int __order));
#ifndef __IS_COMPLEX_HEADER__
EXTRAMATH_ARRDECL(laguerrecof, (int __n, __TYPE_NAME__ *__out));
#endif
EXTRAMATH_DECL(assoclaguerre, (int __n, __TYPE_NAME__ __alpha, __TYPE_NAME__ __x));
EXTRAMATH_DECL(assoclaguerrederiv, (int __n, __TYPE_NAME__ __alpha, __TYPE_NAME__ __x, int __order));
EXTRAMATH_ARRDECL(assoclaguerrecof, (int __n, __TYPE_NAME__ __alpha, __TYPE_NAME__ *__out));

/*
 * Legendre polynomials, associated legendre polynomials, and spherical harmonics. Spherical harmonic coefficients will
 * be a three dimensional array. Each dimension should be l + 1 long for a total size of (l + 1)^3.
 */
EXTRAMATH_DECL(legendre, (int __l, __TYPE_NAME__ __x));
EXTRAMATH_DECL(legendrederiv, (int __l, __TYPE_NAME__ __x, int __order));
#ifndef __IS_COMPLEX_HEADER__
EXTRAMATH_ARRDECL(legendrecof, (int __l, __TYPE_NAME__ *__out));
#endif
EXTRAMATH_DECL(assoclegendre, (int __l, int __m, __TYPE_NAME__ __x));
//EXTRAMATH_DECL(assoclegendrederiv, (unsigned int __l, signed int __m, __TYPE_NAME__ __x, unsigned int __order));
#ifdef __IS_COMPLEX_HEADER__
EXTRAMATH_DECL(spherharm, (int __l, int __m, __TYPE_NAME__ __theta, __TYPE_NAME__ __phi));
#endif
EXTRAMATH_DECL(realspherharm, (int __l, int __m, __TYPE_NAME__ __theta, __TYPE_NAME__ __phi));

/*
 * Chebyshev polynomials.
 */
EXTRAMATH_DECL(chebyshevt, (int __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(chebyshevu, (int __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(chebyshevtderiv, (int __nu, __TYPE_NAME__ __x, int order));
EXTRAMATH_DECL(chebyshevuderiv, (int __nu, __TYPE_NAME__ __x, int order));
#ifndef __IS_COMPLEX_HEADER__
EXTRAMATH_ARRDECL(chebyshevtcof, (int __n, __TYPE_NAME__ *__out));
EXTRAMATH_ARRDECL(chebyshevucof, (int __n, __TYPE_NAME__ *__out));
#endif

EXTRAMATH_DECL(jacobi, (int __n, __TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __x));
EXTRAMATH_DECL(jacobideriv, (int __n, __TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __x, int __order));
EXTRAMATH_ARRDECL(altjacobicof, (int __n, __TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ *__out));

/*
 * Binomial coefficient, since it will be useful.
 */
EXTRAMATH_DECL(fbinom, (int __n, int __k));
