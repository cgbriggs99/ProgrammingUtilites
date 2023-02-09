/*
 * integrated.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

#if defined(__IS_COMPLEX_HEADER__)
/*
 * C does not provide complex gamma functions.
 */
EXTRAMATH_DECL(lgamma, (__TYPE_NAME__ __z));
EXTRAMATH_DECL(tgamma, (__TYPE_NAME__ __z));
#endif

/*
 * Incomplete gamma functions.
 */

/*
 * inGamma(a, z)
 *
 * The incomplete Gamma function, which is the integral from z to infinity of t^(a-1)e^-t.
 */
EXTRAMATH_DECL(inGamma, (__TYPE_NAME__ __a, __TYPE_NAME__ __z));

/*
 * ingamma(a, z)
 *
 * The incomplete gamma function, which is the integral from 0 to z of t^(a-1)e^-t.
 */
EXTRAMATH_DECL(ingamma, (__TYPE_NAME__ __a, __TYPE_NAME__ __z));

/*
 * ingamma2(a, z1, z2)
 *
 * The incomplete gamma function between two values.
 */
EXTRAMATH_DECL(ingamma2, (__TYPE_NAME__ __a, __TYPE_NAME__ __z1, __TYPE_NAME__ __z2));
EXTRAMATH_DECL(inGammareg, (__TYPE_NAME__ __a, __TYPE_NAME__ __z));
EXTRAMATH_DECL(ingammareg, (__TYPE_NAME__ __a, __TYPE_NAME__ __z));
EXTRAMATH_DECL(ingamma2reg, (__TYPE_NAME__ __a, __TYPE_NAME__ __z1, __TYPE_NAME__ __z2));

/*
 * Subfactorial.
 */
EXTRAMATH_DECL(subfact, (__TYPE_NAME__ __z));

/*
 * Beta functions.
 */
EXTRAMATH_DECL(beta, (__TYPE_NAME__ __a, __TYPE_NAME__ __b));
EXTRAMATH_DECL(betainc, (__TYPE_NAME__ __z, __TYPE_NAME__ __a, __TYPE_NAME__ __b));
EXTRAMATH_DECL(betareg, (__TYPE_NAME__ __z, __TYPE_NAME__  __a, __TYPE_NAME__ __b));

/*
 * Polygammas.
 */
EXTRAMATH_DECL(polygamma, (__TYPE_NAME__ __z));
EXTRAMATH_DECL(polygamman, (int __n, __TYPE_NAME__ __z));
EXTRAMATH_DECL(polygammanu, (__TYPE_NAME__ __nu, __TYPE_NAME__ __z));

#if defined(__IS_COMPLEX_HEADER__)
EXTRAMATH_DECL(erf, (__TYPE_NAME__ __z));
#endif

/*
 * Inverse error function.
 */
EXTRAMATH_DECL(ierf, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(ierfc, (__TYPE_NAME__ __x));

/*
 * Other integrals.
 */
EXTRAMATH_DECL(fresnels, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(fresnelc, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(integrale, (__TYPE_NAME__ __nu, __TYPE_NAME__ __x));
EXTRAMATH_DECL(integralei, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(logint, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(sinint, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(cosint, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(sinhint, (__TYPE_NAME__ __x));
EXTRAMATH_DECL(coshint, (__TYPE_NAME__ __x));


