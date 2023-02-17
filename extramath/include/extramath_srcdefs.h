/*
 * extramath_srcdefs.h
 *
 *  Created on: Oct 13, 2021
 *      Author: connor
 */

#ifndef EXTRAMATH_SRCDEFS_H_
#define EXTRAMATH_SRCDEFS_H_

#define EXTRAMATH_FUNDEF(fname, args) EXTRAMATH_FUNDEF1(fname, __TYPESUFF__, args)
#define EXTRAMATH_FUNDEF1(fname, suff, args) __attribute__((nothrow)) __TYPENAME__ __FNAMESRC_SUFF__(fname,suff) args
#define EXTRAMATH_ARRFUNDEF(fname, args) EXTRAMATH_ARRFUNDEF1(fname, __TYPESUFF__, args)
#define EXTRAMATH_ARRFUNDEF1(fname, suff, args)  __attribute__((nothrow)) int __FNAMESRC_SUFF__(fname,suff) args

#define EXTRAMATH_FUNDEFFOR(fname, args) EXTRAMATH_FUNDEFFOR1(fname, __TYPESUFF__, args)
#define EXTRAMATH_FUNDEFFOR1(fname, suff, args) __TYPENAME__ __FNAMEFOR_SUFF__(fname,suff) args
#define EXTRAMATH_ARRFUNDEFFOR(fname, args) EXTRAMATH_FUNDEF1(fname, __TYPESUFF__, args)
#define EXTRAMATH_ARRFUNDEFFOR1(fname, suff, args) int __FNAMEFOR_SUFF__(fname,suff) args
#define __CONCAT1(x, y) __CONCAT(x,y)

#define MAX_ITERS 1000

#define __TYPEVAL_FLOAT__ 0
#define __TYPEVAL_DOUBLE__ 1
#define __TYPEVAL_LONGDOUBLE__ 2
#define __TYPEVAL_COMPLEXFLOAT__ 3
#define __TYPEVAL_COMPLEXDOUBLE__ 4
#define __TYPEVAL_COMPLEXLONGDOUBLE__ 5
#define __TYPEVAL_INT__ 6
#define __TYPEVAL_LONG__ 7
#define __TYPEVAL_LONGLONG__ 8

/*
 * Define macros for manipulating function names and calls.
 */
#if !defined(__IS_COMPLEX__) && !defined(__IS_COMPLEX_HEADER__)
/*
 * These give the various function call macros.
 */

/*
 * Gives the name of a C function. Has no prefix for real, c prefix for complex.
 * Has no suffix for double, f for float, and l for long double.
 * Works with most math.h functions, with the notable exception of fabs.
 */
#define __FNAMESRC__(fname) __FNAMESRC_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_SUFF__(fname,suff) __FNAMESRC_SUFF1__(fname,suff)
#define __FNAMESRC_SUFF1__(fname,suff) fname##suff

/*
 * Gives the name of the C function wrapping a FORTRAN procedure. Has no prefix for real, c prefix for complex.
 * No first suffix for double, f for float, l for long double. A suffix of "_fortran" follows the type suffix.
 */
#define __FNAMEFOR__(fname) __FNAMEFOR_SUFF__(fname,__TYPESUFF__)
#define __FNAMEFOR_SUFF__(fname,suff) __FNAMEFOR_SUFF1__(fname,suff)
#define __FNAMEFOR_SUFF1__(fname,suff) fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure using the standard naming convention. Has no prefix for real, c prefix for complex.
 * No suffix for double (real*8), f for float (real*4), l for long double (real*16). Always ends in an underscore per
 * FORTRAN binding conventions.
 */
#define __FNAMEBIND__(fname) __FNAMEBIND_SUFF__(fname,__TYPESUFF__)
#define __FNAMEBIND_SUFF__(fname,suff) __FNAMEBIND_SUFF1__(fname,suff)
#define __FNAMEBIND_SUFF1__(fname,suff) fname##suff##_

/*
 * Gives the name of a C function. Only for the scalar type, that is if the type is _Complex float, then this will
 * give the function for float.
 */
#define __FNAMESRC_SCAL__(fname) __FNAMESRC_SCAL_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_SCAL_SUFF__(fname,suff) __FNAMESRC_SCAL_SUFF1__(fname,suff)
#define __FNAMESRC_SCAL_SUFF1__(fname,suff) fname##suff

/*
 * Gives the name of a C wrapper function for a FORTRAN procedure. Only for the scalar type, that is if the type
 * is _Complex float, then this will give the function for float.
 */
#define __FNAMEFOR_SCAL__(fname) __FNAMEFOR_SCAL_SUFF__(fname, __TYPESUFF__)
#define __FNAMEFOR_SCAL_SUFF__(fname,suff) __FNAMEFOR_SCAL_SUFF1__(fname,suff)
#define __FNAMEFOR_SCAL_SUFF1__(fname,suff) fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure. Only for the scalar type, that is if the type is _Complex float, then
 * this will give the procedure for float.
 */
#define __FNAMEBIND_SCAL__(fname) __FNAMEBIND_SCAL_SUFF__(fname,__TYPESUFF__)
#define __FNAMEBIND_SCAL_SUFF__(fname,suff) __FNAMEBIND_SCAL_SUFF1__(fname,suff)
#define __FNAMEBIND_SCAL_SUFF1__(fname,suff) fname##suff##_

/*
 * Gives the name of a C function. Has f prefix for real, c prefix for complex.
 * Has no suffix for double, f for float, and l for long double.
 * Works with fabs and cabs.
 */
#define __FNAMESRC_PREF__(fname) __FNAMESRC_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_PREF_SUFF__(fname,suff) __FNAMESRC_PREF_SUFF1__(fname,suff)
#define __FNAMESRC_PREF_SUFF1__(fname,suff) f##fname##suff

/*
 * Gives the name of the C function wrapping a FORTRAN procedure. Has f prefix for real, c prefix for complex.
 * No first suffix for double, f for float, l for long double. A suffix of "_fortran" follows the type suffix.
 */
#define __FNAMEFOR_PREF__(fname) __FNAMEFOR_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEFOR_PREF_SUFF__(fname,suff) __FNAMEFOR_PREF_SUFF1__(fname,suff)
#define __FNAMEFOR_PREF_SUFF1__(fname,suff) f##fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure using the standard naming convention. Has f prefix for real, c prefix for complex.
 * No suffix for double (real*8), f for float (real*4), l for long double (real*16). Always ends in an underscore per
 * FORTRAN binding conventions.
 */
#define __FNAMEBIND_PREF__(fname) __FNAMEBIND_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEBIND_PREF_SUFF__(fname,suff) __FNAMEBIND_PREF_SUFF1__(fname,suff)
#define __FNAMEBIND_PREF_SUFF1__(fname,suff) f##fname##suff##_

/*
 * Gives the name of a C function. Has f prefix. Only for the scalar. That is, if the type is _Complex float,
 * then this will give the function for float. Works with fabs and cabs.
 */
#define __FNAMESRC_SCAL_PREF__(fname) __FNAMESRC_SCAL_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_SCAL_PREF_SUFF__(fname,suff) __FNAMESRC_SCAL_PREF_SUFF1__(fname,suff)
#define __FNAMESRC_SCAL_PREF_SUFF1__(fname,suff) f##fname##suff

/*
 * Gives the name of a C wrapper function for a FORTRAN procedure. Only for the scalar type, that is if the type
 * is _Complex float, then this will give the function for float. Has an f prefix.
 */
#define __FNAMEFOR_SCAL_PREF__(fname) __FNAMESRC_SCAL_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEFOR_SCAL_PREF_SUFF__(fname,suff) __FNAMESRC_SCAL_PREF_SUFF1__(fname,suff)
#define __FNAMEFOR_SCAL_PREF_SUFF1__(fname,suff) f##fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure. Only for the scalar type, that is if the type is _Complex float, then
 * this will give the procedure for float. Has f prefix.
 */
#define __FNAMEBIND_SCAL_PREF__(fname) __FNAMEBIND_SCAL_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEBIND_SCAL_PREF_SUFF__(fname,suff) __FNAMEBIND_SCAL_PREF_SUFF1__(fname,suff)
#define __FNAMEBIND_SCAL_PREF_SUFF1__(fname,suff) f##fname##suff##_
#else
/*
 * These give the various function call macros.
 */

/*
 * Gives the name of a C function. Has no prefix for real, c prefix for complex.
 * Has no suffix for double, f for float, and l for long double.
 * Works with most math.h functions, with the notable exception of fabs.
 */
#define __FNAMESRC__(fname) __FNAMESRC_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_SUFF__(fname,suff) __FNAMESRC_SUFF1__(fname,suff)
#define __FNAMESRC_SUFF1__(fname,suff) c##fname##suff

/*
 * Gives the name of the C function wrapping a FORTRAN procedure. Has no prefix for real, c prefix for complex.
 * No first suffix for double, f for float, l for long double. A suffix of "_fortran" follows the type suffix.
 */
#define __FNAMEFOR__(fname) __FNAMEFOR_SUFF__(fname,__TYPESUFF__)
#define __FNAMEFOR_SUFF__(fname,suff) __FNAMEFOR_SUFF1__(fname,suff)
#define __FNAMEFOR_SUFF1__(fname,suff) c##fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure using the standard naming convention. Has no prefix for real, c prefix for complex.
 * No suffix for double (real*8), f for float (real*4), l for long double (real*16). Always ends in an underscore per
 * FORTRAN binding conventions.
 */
#define __FNAMEBIND__(fname) __FNAMEBIND_SUFF__(fname,__TYPESUFF__)
#define __FNAMEBIND_SUFF__(fname,suff) __FNAMEBIND_SUFF1__(fname,suff)
#define __FNAMEBIND_SUFF1__(fname,suff) c##fname##suff##_

/*
 * Gives the name of a C function. Only for the scalar type, that is if the type is _Complex float, then this will
 * give the function for float.
 */
#define __FNAMESRC_SCAL__(fname) __FNAMESRC_SCAL_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_SCAL_SUFF__(fname,suff) __FNAMESRC_SCAL_SUFF1__(fname,suff)
#define __FNAMESRC_SCAL_SUFF1__(fname,suff) fname##suff

/*
 * Gives the name of a C wrapper function for a FORTRAN procedure. Only for the scalar type, that is if the type
 * is _Complex float, then this will give the function for float.
 */
#define __FNAMEFOR_SCAL__(fname) __FNAMEFOR_SCAL_SUFF__(fname, __TYPESUFF__)
#define __FNAMEFOR_SCAL_SUFF__(fname,suff) __FNAMEFOR_SCAL_SUFF1__(fname,suff)
#define __FNAMEFOR_SCAL_SUFF1__(fname,suff) fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure. Only for the scalar type, that is if the type is _Complex float, then
 * this will give the procedure for float.
 */
#define __FNAMEBIND_SCAL__(fname,suff) __FNAMEBIND_SCAL_SUFF__(fname,suff)
#define __FNAMEBIND_SCAL_SUFF__(fname,suff) __FNAMEBIND_SCAL_SUFF1__(fname,suff)
#define __FNAMEBIND_SCAL_SUFF1__(fname,suff) fname##suff##_

/*
 * Gives the name of a C function. Has f prefix for real, c prefix for complex.
 * Has no suffix for double, f for float, and l for long double.
 * Works with fabs and cabs.
 */
#define __FNAMESRC_PREF__(fname) __FNAMESRC_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_PREF_SUFF__(fname,suff) __FNAMESRC_PREF_SUFF1__(fname,suff)
#define __FNAMESRC_PREF_SUFF1__(fname,suff) c##fname##suff

/*
 * Gives the name of the C function wrapping a FORTRAN procedure. Has f prefix for real, c prefix for complex.
 * No first suffix for double, f for float, l for long double. A suffix of "_fortran" follows the type suffix.
 */
#define __FNAMEFOR_PREF__(fname) __FNAMEFOR_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEFOR_PREF_SUFF__(fname,suff) __FNAMEFOR_PREF_SUFF1__(fname,suff)
#define __FNAMEFOR_PREF_SUFF1__(fname,suff) c##fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure using the standard naming convention. Has f prefix for real, c prefix for complex.
 * No suffix for double (real*8), f for float (real*4), l for long double (real*16). Always ends in an underscore per
 * FORTRAN binding conventions.
 */
#define __FNAMEBIND_PREF__(fname) __FNAMEBIND_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEBIND_PREF_SUFF__(fname,suff) __FNAMEBIND_PREF_SUFF1__(fname,suff)
#define __FNAMEBIND_PREF_SUFF1__(fname,suff) c##fname##suff##_

/*
 * Gives the name of a C function. Has f prefix. Only for the scalar. That is, if the type is _Complex float,
 * then this will give the function for float. Works with fabs and cabs.
 */
#define __FNAMESRC_SCAL_PREF__(fname) __FNAMESRC_SCAL_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMESRC_SCAL_PREF_SUFF__(fname,suff) __FNAMESRC_SCAL_PREF_SUFF1__(fname,suff)
#define __FNAMESRC_SCAL_PREF_SUFF1__(fname,suff) f##fname##suff

/*
 * Gives the name of a C wrapper function for a FORTRAN procedure. Only for the scalar type, that is if the type
 * is _Complex float, then this will give the function for float. Has an f prefix.
 */
#define __FNAMEFOR_SCAL_PREF__(fname) __FNAMEFOR_SCAL_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEFOR_SCAL_PREF_SUFF__(fname,suff) __FNAMEFOR_SCAL_PREF_SUFF1__(fname,suff)
#define __FNAMEFOR_SCAL_PREF_SUFF1__(fname,suff) f##fname##suff##_fortran

/*
 * Gives the name of a FORTRAN procedure. Only for the scalar type, that is if the type is _Complex float, then
 * this will give the procedure for float. Has f prefix.
 */
#define __FNAMEBIND_SCAL_PREF__(fname) __FNAMEBIND_SCAL_PREF_SUFF__(fname,__TYPESUFF__)
#define __FNAMEBIND_SCAL_PREF_SUFF__(fname,suff) __FNAMEBIND_SCAL_PREF_SUFF1__(fname,suff)
#define __FNAMEBIND_SCAL_PREF_SUFF1__(fname,suff) f##fname##suff##_
#endif


#endif /* EXTRAMATH_SRCDEFS_H_ */
