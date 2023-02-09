/*
 * utils.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_H_
#error "Do not import this file on its own."
#endif

extern __SCALAR_TYPE__ __FNAME__(epsilon)(__TYPE_NAME__ __x);

extern int __FNAME__(isinteger)(__TYPE_NAME__ __z);

extern int __FNAME__(isnatural)(__TYPE_NAME__ __z);

extern int __FNAME__(absconv)(__TYPE_NAME__ __center, __TYPE_NAME__ __diff);

#if !defined(__IS_COMPLEX__) && !defined(__IS_COMPLEX_HEADER__)
EXTRAMATH_DECL(real, (__TYPE_NAME__ __z));
EXTRAMATH_DECL(imag, (__TYPE_NAME__ __z));
#endif
//extern int __FORTRAN__(__FNAME__(absconv))(__TYPE_NAME__ __center, __TYPE_NAME__ __diff);

