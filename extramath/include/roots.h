/*
 * roots.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

/*
 * Find the roots of the polynomial represented by __coefs__. The order in this array is the constant term first,
 * then the linear, and so on.
 */
EXTRAMATH_ARRDECL(polyroots,(const __TYPE_NAME__ *__coefs__, unsigned int __length__, __TYPE_NAME__ *roots));

