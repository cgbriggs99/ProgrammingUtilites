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
EXTRAMATH_DECL(polynomial, (__TYPE_NAME__ __x, unsigned int __len, const __TYPE_NAME__ *coefs));

/*
 * Hermite polynomials. There is one more coefficient than the order, so __out should have length __n + 1.
 * Uses the physicist's definition.
 */
EXTRAMATH_DECL(hermite, (unsigned int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(hermitederiv, (unsigned int __n, __TYPE_NAME__ __x, unsigned int __order));
EXTRAMATH_ARRDECL(hermitecof, (unsigned int __n, __TYPE_NAME__ *__out));

/*
 * Laguerre polynomials and associated laguerre polynomials.
 */
EXTRAMATH_DECL(laguerre, (unsigned int __n, __TYPE_NAME__ __x));
EXTRAMATH_DECL(laguerrederiv, (unsigned int __n, __TYPE_NAME__ __x, unsigned int __order));
EXTRAMATH_ARRDECL(laguerrecof, (unsigned int __n, __TYPE_NAME__ *__out));
EXTRAMATH_DECL(assoclaguerre, (unsigned int __n, __TYPE_NAME__ __alpha, __TYPE_NAME__ __x));
EXTRAMATH_DECL(assoclaguerrederiv, (unsigned int __n, __TYPE_NAME__ __alpha, __TYPE_NAME__ __x, unsigned int __order));
EXTRAMATH_ARRDECL(assoclaguerrecof, (unsigned int __n, __TYPE_NAME__ __alpha, __TYPE_NAME__ *__out));

/*
 * Legendre polynomials, associated legendre polynomials, and spherical harmonics. Spherical harmonic coefficients will
 * be a three dimensional array. Each dimension should be l + 1 long for a total size of (l + 1)^3.
 */
EXTRAMATH_DECL(legendre, (unsigned int __l, __TYPE_NAME__ __x));
EXTRAMATH_DECL(legendrederiv, (unsigned int __l, __TYPE_NAME__ __x, unsigned int __order));
EXTRAMATH_ARRDECL(legendrecof, (unsigned int __l, __TYPE_NAME__ *__out));
EXTRAMATH_DECL(assoclegendre, (unsigned int __l, signed int __m, __TYPE_NAME__ __x));
//EXTRAMATH_DECL(assoclegendrederiv, (unsigned int __l, signed int __m, __TYPE_NAME__ __x, unsigned int __order));
#ifdef __IS_COMPLEX_HEADER__
EXTRAMATH_DECL(spherharm, (unsigned int __l, signed int __m, __TYPE_NAME__ __theta, __TYPE_NAME__ __phi));
#endif
EXTRAMATH_DECL(realspherharm, (unsigned int __l, signed int __m, __TYPE_NAME__ __theta, __TYPE_NAME__ __phi));

/*
 * Chebyshev polynomials.
 */
EXTRAMATH_DECL(chebyshevt, (unsigned int __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(chebyshevu, (unsigned int __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(chebyshevtderiv, (unsigned int __nu, __TYPE_NAME__ __x, unsigned int order));
EXTRAMATH_DECL(chebyshevuderiv, (unsigned int __nu, __TYPE_NAME__ __x, unsigned int order));
EXTRAMATH_ARRDECL(chebyshevtcof, (unsigned int __n, __TYPE_NAME__ *__out));
EXTRAMATH_ARRDECL(chebyshevucof, (unsigned int __n, __TYPE_NAME__ *__out));

EXTRAMATH_DECL(jacobi, (unsigned int __n, __TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __x));
EXTRAMATH_DECL(jacobideriv, (unsigned int __n, __TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ __x, unsigned int __order));
EXTRAMATH_ARRDECL(altjacobicof, (unsigned int __n, __TYPE_NAME__ __a, __TYPE_NAME__ __b, __TYPE_NAME__ *__out));

/*
 * Binomial coefficient, since it will be useful.
 */
EXTRAMATH_DECL(fbinom, (int __n, int __k));
